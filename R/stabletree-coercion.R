### -- as.stabletree generic and methods --------------------------------------

### S3
as.stabletree <- function(x, ...) {
  UseMethod("as.stabletree")
}

### set S4 (for party's RandomForest), default is S3
setGeneric("as.stabletree")

### function for extracting split information from trees (partykit::partysplit like list)
extract_breaks <- function(x, x_names, x_classes, x_levels, x_nlevels, extract_split_fun, start_0 = TRUE) {
  sp <- extract_split_fun(x, x_names = x_names, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels)
  if (!is.null(sp)) {
    vi <- sapply(sp, "[[", "varid")
    if (start_0 == TRUE) {
      vi <- vi - 1L
    }
  } else vi <- NULL
  br <- lapply(sp, "[[", "breaks")
  id <- lapply(sp, "[[", "index")
  names(id) <- names(br) <- x_names[vi]
  br <- lapply(x_names, function(n) {
    brs <- br[names(br) == n]
    ids <- id[names(id) == n]
    if (length(brs) > 0L || length(ids) > 0L) {
      if (is.null(brs[[n]])) {
        ans <- do.call("rbind", ids)
        if (!is.null(ans)) {
          rownames(ans) <- NULL
          ## sometimes the following fails for weird reasons, see also below
          tmp <- try(colnames(ans) <- x_levels[[n]], silent = TRUE)
          if (class(tmp) == "try-error") class(ans) <- "try-error"
        }
      } else {
        ans <- unlist(brs)
        names(ans) <- NULL
      }
    } else ans <- NULL
    return(ans)
  })
  names(br) <- x_names
  return(br)
}



### function to add levels to list with breakpoints
add_levels <- function(x, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels) {
  nm <- names(x)
  ans <- lapply(nm, function(n) {
    br <- x[[n]]
    if (!is.null(br)) {
      if (x_classes[n] == "ordered") {
        br <- ordered(br, levels = 1L:x_nlevels[[n]], labels = x_levels[[n]])
      }
      br
    } else NULL
  })
  names(ans) <- nm
  return(ans)
}

### as.stabletree.randomForest (randomForest)
as.stabletree.randomForest <- function(x, applyfun = NULL, cores = NULL, ...) {
  call <- try(getCall(x), silent = TRUE)
  sampler <- list(method = "randomForest::randomForest", sampler = "randomForest::randomForest")
  B <- x$forest$ntree

  ## get terms
  tr <- terms(x)

  ## get envoronment of x
  env <- try(environment(tr), silent = TRUE)
  if (inherits(env, "try-error")) env <- NULL

  ## extract information from call
  sfit <- call$subset
  wfit <- call$weights
  dfit <- call$data

  ## get data
  data <- x$data
  if (is.null(data)) {
    if (is.null(dfit)) {
      ## there is no data object
      data <- NULL
    } else
      ## get local copy of data object from where x was generated
      data <- eval(dfit, envir = env, enclos = parent.frame())
      if (!is.null(sfit)) {
        sfit <- eval(sfit, envir = env, enclos = parent.frame())
        data <- subset(data, subset = sfit)
      }
  }

  ## facilitate parallelization
  if (is.null(applyfun)) {
    applyfun <- if (is.null(cores)) {
      lapply
    } else {
      function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
    }
  }

  ## get trees of the forest
  xx <- applyfun(seq_len(B), FUN = function(b) randomForest::getTree(x, b, TRUE))

  ## extract names of all variables and omit response
  mf <- model.frame(tr, data = data)
  yi <- attr(tr, "response")
  x_classes <- sapply(mf[, - yi, drop = FALSE], function(x) class(x)[1])
  x_levels <- sapply(mf[, - yi, drop = FALSE], levels, simplify = FALSE)
  x_nlevels <- sapply(mf[, - yi, drop = FALSE], nlevels)
  x_names <- names(mf[- yi])

  ## there is no "original" tree so we set vs0 to 0 and br0 to NULL
  nvar <- length(x_names)
  vs0 <- rep(0, nvar)
  names(vs0) <- x_names
  br0 <- vector("list", nvar)
  names(br0) <- x_names

  ## function for computing the binary expansion of an integer
  ## important: see ?randomForest::getTree for binary expansion of categoricals
  binary_expansion <- function(x, n) {
    remainder <- numeric(n)
    for (i in seq_len(n)) {
      remainder[i] <- x %% 2L
      x <- floor(x / 2L)
    }
    stopifnot(identical(x, 0)) ## safety check
    return(remainder)
  }

  ## function for extracting splits from trees (generate partykit::partysplit like list)
  ## 1s get send to the left daughter node, 0s to the right; mapped to 1 and 2 (for consistency)
  ## info contains a level of the splitting variable if this split sends this level alone to the
  ## left or right daughter node (and subsequent daughter nodes cannot use this level for splitting)
  ## is_factor additionally states whether the splitting variable is a factor
  extract_split <- function(x, x_names, x_classes, x_levels, x_nlevels) {
    ids <- as.numeric(rownames(x[x$status != -1, ]))
    ids_c <- as.character(ids)

    splits <- lapply(ids, function(id) {
      varid <- match(x$"split var"[id], x_names)
      is_factor <- x_classes[varid] == "factor"
      is_ordered <- x_classes[varid] == "ordered" ## note that ordered factors are then treated as numerics
      index <- if (is_factor) {
        as.integer(!as.logical(binary_expansion(x$"split point"[id], n = x_nlevels[varid]))) + 1L
        #as.integer(!as.logical(binaryLogic::as.binary(x$"split point"[id], littleEndian = TRUE, n = x_nlevels[varid]))) + 1L
      } else c(1L, 2L)
      list(varid = varid,
           breaks = if (is_factor) NULL else if (is_ordered) floor(x$"split point"[id]) else x$"split point"[id],
           index = index,
           right = NULL, prob = NULL,
           info = if (is_factor) which(index == which(table(index) == 1L)) else NULL,
           is_factor = is_factor)
    })
    names(splits) <- ids_c

    ## FIXME: this does not handle the scenario of a split not being possible
    ## because no individuals were left fulfilling the splitting criterion in
    ## the subsample and there is no easy fix, at least it is documented
    ## properly
    for (i in seq_len(length(ids))) { ## set relevant indices to NA if levels could not be used
      id <- ids[i]
      id_c <- ids_c[i]
      if (splits[[id_c]]$is_factor && length(splits[[id_c]]$info > 0L)) {
        daughter_ids <- c(x$"left daughter"[id], x$"right daughter"[id])
        varids <- match(x$"split var"[daughter_ids], x_names)
        daughter_ids <- daughter_ids[which(varids == splits[[id_c]]$varid)]
        if (length(daughter_ids > 0L)) {
          daughter_ids_c <- as.character(daughter_ids)
          for (daughter_id_c in daughter_ids_c) {
            splits[[daughter_id_c]]$index[splits[[id_c]]$info] <- NA
          }
        }
      }
    }
    return(splits)
  }

  ## function for extracting variable id from trees
  extract_varid <- function(x, x_names) {
    vi <- sort(unique(as.numeric(factor(x$"split var", levels = x_names))))
    vi <- as.numeric(x_names %in% x_names[vi])
    names(vi) <- x_names
    return(vi)
  }

  ## selection proportions
  vi <- applyfun(xx, FUN = extract_varid, x_names = x_names)
  vi_mat <- do.call("rbind", vi)

  ## breakpoints
  br <- applyfun(xx, FUN = extract_breaks, x_names = x_names, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels, extract_split_fun = extract_split, start_0 = FALSE)
  ## weird internal error handling
  ## if some errors occured within extract_split, extract_varid_, or extract_breaks
  ## drop these resamples from vs and br
  tmp <- which(sapply(br, function(x) any(sapply(x, class) == "try-error")))
  if(length(tmp)) {
    vi_mat <- vi_mat[-tmp, ]
    br[tmp] <- NULL
    warning("Due to internal coercion errors, only the results of ", B - length(tmp), "resamples are returned.")
  }
  br <- lapply(x_names, function(n) {
    if (x_classes[n] == "factor") {
      do.call("rbind", lapply(br, "[[", n))
    } else {
      unlist(lapply(br, "[[", n))
    }
  })
  names(br) <- x_names

  ## build stabletree object
  rval <- list(
    call = call,
    B = B - length(tmp),
    sampler = sampler,
    vs0 = vs0,
    br0 = br0,
    vs = vi_mat,
    br = add_levels(br, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels),
    classes = x_classes,
    trees = NULL # not really usable, so we set trees to NULL
  )
  class(rval) <- "stabletree"
  return(rval)
}

### as.stabletree.RandomForest (party)
as.stabletree.RandomForest <- function(x, applyfun = NULL, cores = NULL, ...) {
  mf <- x@data
  call <- mf
  sampler <- list(method = "party::cforest", sampler = "party::cforest")

  B <- length(x@ensemble)

  ## get terms not needed

  ## get envoronment of x
  env <- try(mf@env, silent = TRUE)
  if (inherits(env, "try-error")) env <- NULL

  ## extract information from call not needed

  ## get data (should always work)
  data <- mf@get("response")

  ## facilitate parallelization
  if (is.null(applyfun)) {
    applyfun <- if (is.null(cores)) {
      lapply
    } else {
      function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
    }
  }

  ## get trees of the forest
  xx <- x@ensemble

  ## extract names of all variables and omit response
  x_classes <- sapply(mf@get("input"), function(x) class(x)[1L])
  x_levels <- sapply(mf@get("input"), levels, simplify = FALSE)
  x_nlevels <- sapply(mf@get("input"), nlevels)
  x_names <- names(mf@get("input"))

  ## there is no "original" tree so we set vs0 to 0 and br0 to NULL
  nvar <- length(x_names)
  vs0 <- rep(0, nvar)
  names(vs0) <- x_names
  br0 <- vector("list", nvar)
  names(br0) <- x_names

  ## function to get nodeids and nodes of a tree
  extract_nodeids <- function(x) {
    root_id <- x[[1L]]
    root_term <- x[[4L]]
    left <- if (!is.null(x[[8L]])) {
      extract_nodeids(x[[8L]])
    } else {
      NULL
    }
    left_id <- left$nodeids
    left_term <- left$terminal
    right <- if (!is.null(x[[9L]])) {
      extract_nodeids(x[[9L]])
    } else {
      NULL
    }
    right_id <- right$nodeids
    right_term <- right$terminal
    return(data.frame(nodeids = c(root_id, left_id, right_id),
      terminal = c(root_term, left_term, right_term)))
  }

  ## function to switch the index if "split"$toleft ist FALSE
  switch_right <- function(index, toleft, table) {
    ## switch the index if "split"$toleft ist FALSE
    if (!toleft) {
      index[index == 1L] <- 0L
      index[index == 2L] <- 1L
      index[index == 0L] <- 2L
    }
    ## check which levels are not available for the daughter nodes
    if (length(table) > 0L) {
      index[which(table == 0L)] <- NA
    }
    return(index)
  }

  ## function to make a "BinaryTree" with no slots but the tree slots out of a
  ## tree; we need this so we can later extract the nodes via party::nodes
  as_BinaryTree <- function(tree) {
    new("BinaryTree", tree = tree)
  }

  ## function for extracting splits from trees (generate partykit::partysplit like list)
  ## 1s get send to the left daughter node, 0s to the right; mapped to 1 and 2 (for consistency)
  ## If "split"$toleft ist FALSE, we switch the "index"; moreover, for factors we need to check
  ## if the level is actually available for the daugther nodes (see "split"$table)
  extract_split <- function(x, x_names, x_classes, x_levels, x_nlevels) {
    ids <- extract_nodeids(x)
    ids <- ids$nodeids[ids$terminal == FALSE]
    ids_c <- as.character(ids)
    nodes <- lapply(ids, function(i) party::nodes(as_BinaryTree(x), i)[[1]])
    #nodes <- lapply(ids, function(i) .Call("R_get_nodebynum", x, i, PACKAGE = "party"))

    splits <- lapply(nodes, function(node) {
      split <- if (is.null(node[[5L]])) node[[6L]] else node[[5L]] ## psplit/split
      split[[5L]] <- if (is.null(split[[5L]])) FALSE else if (split[[5L]] == 1L) TRUE ## terminal
      varid <- split[[1L]]
      is_factor <- !split[[2L]] ## i.e., not ordered
      index <- if (is_factor) {
        (!split[[3L]]) + 1L ## map the index to 1 and 2 (see comment above)
      } else c(1L, 2L)
      list(varid = varid,
           breaks = if (is_factor) NULL else split[[3L]],
           index = switch_right(index, split[[5L]], if (is_factor) split[[6L]] else NULL),
           right = NULL, prob = NULL, info = NULL)
    })
    names(splits) <- ids_c
    return(splits)
  }

  ## function for extracting variable id from trees
  extract_varid <- function(x, x_names, x_classes, x_levels, x_nlevels) {
    sp <- extract_split(x, x_names = x_names, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels)
    if (!is.null(sp)) {
      vi <- sapply(sp, "[[", "varid") ## no -1L see extract_split above
      vi <- sort(unique(vi))
    } else vi <- NULL
    vi <- as.numeric(x_names %in% x_names[vi])
    names(vi) <- x_names
    return(vi)
  }

  ## selection proportions
  vi <- applyfun(xx, FUN = extract_varid, x_names = x_names, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels)
  vi_mat <- do.call("rbind", vi)

  ## breakpoints
  br <- applyfun(xx, FUN = extract_breaks, x_names = x_names, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels, extract_split_fun = extract_split, start_0 = FALSE)
  ## weird internal error handling
  ## if some errors occured within extract_split, extract_varid_, or extract_breaks
  ## drop these resamples from vs and br
  tmp <- which(sapply(br, function(x) any(sapply(x, class) == "try-error")))
  if(length(tmp)) {
    vi_mat <- vi_mat[-tmp, ]
    br[tmp] <- NULL
    warning("Due to internal coercion errors, only the results of ", B - length(tmp), "resamples are returned.")
  }
  br <- lapply(x_names, function(n) {
    if (x_classes[n] == "factor") {
      do.call("rbind", lapply(br, "[[", n))
    } else {
      unlist(lapply(br, "[[", n))
    }
  })
  names(br) <- x_names

  ## build stabletree object
  rval <- list(
    call = call,
    B = B - length(tmp),
    sampler = sampler,
    vs0 = vs0,
    br0 = br0,
    vs = vi_mat,
    br = add_levels(br, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels),
    classes = x_classes,
    trees = NULL # not really usable, so set trees to NULL
  )
  class(rval) <- "stabletree"
  return(rval)
}

setMethod("as.stabletree", "RandomForest", as.stabletree.RandomForest)

### as.stabletree.cforest (partykit)
as.stabletree.cforest <- function(x, applyfun = NULL, cores = NULL, savetrees = FALSE, ...) {
  call <- try(getCall(x), silent = TRUE)
  sampler <- list(method = "partykit::cforest", sampler = "partykit::cforest")
  B <- length(x$nodes)

  ## get terms
  tr <- terms(x)

  ## get envoronment of x
  env <- try(environment(tr), silent = TRUE)
  if (inherits(env, "try-error")) env <- NULL

  ## extract information from call
  sfit <- call$subset
  wfit <- call$weights
  dfit <- call$data

  ## get data
  data <- x$data
  if (is.null(data)) {
    if (is.null(dfit)) {
      ## there is no data object
      data <- NULL
    } else
      ## get local copy of data object from where x was generated
      data <- eval(dfit, envir = env, enclos = parent.frame())
      if (!is.null(sfit)) {
        sfit <- eval(sfit, envir = env, enclos = parent.frame())
        data <- subset(data, subset = sfit)
      }
  }

  ## facilitate parallelization
  if (is.null(applyfun)) {
    applyfun <- if (is.null(cores)) {
      lapply
    } else {
      function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
    }
  }

  ## get trees of the forest
  xx <- applyfun(seq_len(B), FUN = function(b) partykit::gettree(x, b))

  ## extract names of all variables and omit response
  mf <- model.frame(tr, data = data)
  yi <- attr(tr, "response")
  x_classes <- sapply(mf[, - yi, drop = FALSE], function(x) class(x)[1])
  x_levels <- sapply(mf[, - yi, drop = FALSE], levels, simplify = FALSE)
  x_nlevels <- sapply(mf[, - yi, drop = FALSE], nlevels)
  x_names <- names(mf[- yi])

  ## there is no "original" tree so we set vs0 to 0 and br0 to NULL
  nvar <- length(x_names)
  vs0 <- rep(0, nvar)
  names(vs0) <- x_names
  br0 <- vector("list", nvar)
  names(br0) <- x_names

  ## function for extracting splits from trees
  extract_split <- function(x, x_names, x_classes, x_levels, x_nlevels) {
    ids <- nodeids(x)
    ids <- ids[- nodeids(x, terminal = TRUE)]
    nodeapply(x, ids = ids, FUN = split_node)
  }

  ## function for extracting variable id from trees
  extract_varid <- function(x, x_names, x_classes, x_levels, x_nlevels) {
    sp <- extract_split(x)
    if (!is.null(sp)) {
      vi <- sapply(sp, "[[", "varid") - 1L
      vi <- sort(unique(vi))
    } else vi <- NULL
    vi <- as.numeric(x_names %in% x_names[vi])
    names(vi) <- x_names
    return(vi)
  }

  ## selection proportions
  vi <- applyfun(xx, FUN = extract_varid, x_names = x_names, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels)
  vi_mat <- do.call("rbind", vi)

  ## breakpoints
  br <- applyfun(xx, FUN = extract_breaks, x_names = x_names, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels, extract_split_fun = extract_split, start_0 = TRUE)
  ## weird internal error handling
  ## if some errors occured within extract_split, extract_varid_, or extract_breaks
  ## drop these resamples from vs and br
  tmp <- which(sapply(br, function(x) any(sapply(x, class) == "try-error")))
  if(length(tmp)) {
    vi_mat <- vi_mat[-tmp, ]
    br[tmp] <- NULL
    warning("Due to internal coercion errors, only the results of ", B - length(tmp), "resamples are returned.")
  }
  br <- lapply(x_names, function(n) {
    if (x_classes[n] == "factor") {
      do.call("rbind", lapply(br, "[[", n))
    } else {
      unlist(lapply(br, "[[", n))
    }
  })
  names(br) <- x_names

  ## build stabletree object
  rval <- list(
    call = call,
    B = B - length(tmp),
    sampler = sampler,
    vs0 = vs0,
    br0 = br0,
    vs = vi_mat,
    br = add_levels(br, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels),
    classes = x_classes,
    trees = if (savetrees) {
      xx
    } else NULL
  )
  class(rval) <- "stabletree"
  return(rval)
}

### -- as.stabletree.ranger (ranger) ----------------------------------------
as.stabletree.ranger <- function(x, applyfun = NULL, cores = NULL, ...) {
  call <- try(getCall(x), silent = TRUE)
  formula <- if (is.null(call$formula)) as.formula(call[[2L]]) else as.formula(call$formula)

  ## if no ranger.forest object was saved this does not work
  if (!is.null(call$write.forest) && call$write.forest == FALSE) {
    stop("Refit the forest using `write.forest = TRUE`.")
  }

  sampler <- list(method = "ranger::ranger", sampler = "ranger::ranger")
  B <- x$num.trees

  ## get envoronment of x
  env <- try(environment(formula), silent = TRUE)
  if (inherits(env, "try-error")) env <- NULL

  ## extract information from call
  #sfit <- call$subset
  #wfit <- call$weights
  dfit <- call$data

  ## get data
  #data <- x$data
  #if (is.null(data)) {
  #  if (is.null(dfit)) {
  #    ## there is no data object
  #    data <- NULL
  #  } else
  #    ## get local copy of data object from where x was generated
  #    data <- eval(dfit, envir = env, enclos = parent.frame())
  #    if (!is.null(sfit)) {
  #      sfit <- eval(sfit, envir = env, enclos = parent.frame())
  #      data <- subset(data, subset = sfit)
  #    }
  #}
  data <- eval(dfit, envir = env, enclos = parent.frame())

  # check for the terms here because due to no terms slot and "~" we need to
  # supply the data
  tr <- terms.formula(formula, data = data)

  # interactions are not supported (yet)
  if (length(grep(":", x = attr(tr, which = "term.labels")))) {
    stop("Interaction terms of variables are not supported (yet).")
  }

  # ranger behaves differently with respect to unordered factors:
  # options "ignore, "order", "partition", see "respect.unordered.factors"
  # "extratrees" splitrule default: "partition"
  # for all other splitrules default: "ignore"
  respect_unordered_factors <- if (is.null(call$respect.unordered.factors)) {
    if (!is.null(call$splitrule) && call$splitrule == "extratrees") {
     "partition"
    } else {
      "ignore"
    }
  } else {
    call$respect.unordered.factors
  }

  ## facilitate parallelization
  if (is.null(applyfun)) {
    applyfun <- if (is.null(cores)) {
      lapply
    } else {
      function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
    }
  }

  ## get trees of the forest
  xx <- applyfun(seq_len(B), FUN = function(b) ranger::treeInfo(x, b))

  ## extract names of all variables and omit response
  mf <- model.frame(tr, data = data)
  yi <- attr(tr, "response")
  x_classes <- sapply(mf[, - yi, drop = FALSE], function(x) class(x)[1])
  x_levels <- sapply(mf[, - yi, drop = FALSE], levels, simplify = FALSE)
  x_nlevels <- sapply(mf[, - yi, drop = FALSE], nlevels)
  x_names <- if (!is.null(x$forest$covariate.levels)) {
    x$forest$covariate.levels
  } else {
    names(mf[- yi])
  }

  ## there is no "original" tree so we set vs0 to 0 and br0 to NULL
  nvar <- length(x_names)
  vs0 <- rep(0, nvar)
  names(vs0) <- x_names
  br0 <- vector("list", nvar)
  names(br0) <- x_names

  ## function for extracting splits from trees
  extract_split <- function(x, x_names, x_classes, x_levels, x_nlevels) {
    ids <- x$nodeID[!x$terminal] + 1L
    ids_c <- as.character(ids)

    splits <- lapply(ids, function(id) {
      varid <- match(x$splitvarNam[id], x_names)
      is_factor <- x_classes[varid] == "factor"
      is_ordered <- x_classes[varid] == "ordered" ## note that ordered factors are then treated as numeric
      index <- if (is_factor) {
        if (respect_unordered_factors == "partition") {
          left <- rep.int(1L, x_nlevels[varid])
          left[as.integer(strsplit(x$splitval[id], split = ",")[[1L]])] <- 2
          left
        } else {
          left <- rep(1L, times = floor(x$splitval[id]))
          c(left, rep.int(2L, times = x_nlevels[varid] - length(left)))
        }
      } else c(1L, 2L)
      list(varid = varid,
           breaks = if (is_factor) NULL else if (is_ordered) floor(as.numeric(x$splitval[id])) else as.numeric(x$splitval[id]),
           index = index,
           right = NULL, prob = NULL,
           info = if (is_factor) which(index == which(table(index) == 1L)) else NULL,
           is_factor = is_factor)
    })

    split(x, x$nodeID)
    names(splits) <- ids_c

    ## FIXME: this does not handle the scenario of a split not being possible
    ## because no individuals were left fulfilling the splitting criterion in
    ## the subsample and there is no easy fix, at least it is documented
    ## properly
    for (i in seq_len(length(ids))) { ## set relevant indices to NA if levels could not be used
      id <- ids[i]
      id_c <- ids_c[i]
      if (splits[[id_c]]$is_factor && length(splits[[id_c]]$info > 0L)) {
        daughter_ids <- c(x$leftChild[id], x$rightChild[id]) + 1L
        varids <- match(x$splitvarName[daughter_ids], x_names)
        daughter_ids <- daughter_ids[which(varids == splits[[id_c]]$varid)]
        if (length(daughter_ids > 0L)) {
          daughter_ids_c <- as.character(daughter_ids)
          for (daughter_id_c in daughter_ids_c) {
            splits[[daughter_id_c]]$index[splits[[id_c]]$info] <- NA
          }
        }
      }
    }
    return(splits)

  }

  ## function for extracting variable id from trees
  extract_varid <- function(x, x_names, x_classes, x_levels, x_nlevels) {
    sp <- extract_split(x, x_names, x_classes, x_levels, x_nlevels)
    if (!is.null(sp)) {
      vi <- sapply(sp, "[[", "varid") ## no -1L see extract_split above
      vi <- sort(unique(vi))
    } else vi <- NULL
    vi <- as.numeric(x_names %in% x_names[vi])
    names(vi) <- x_names
    return(vi)
  }

  ## selection proportions
  vi <- applyfun(xx, FUN = extract_varid, x_names = x_names, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels)
  vi_mat <- do.call("rbind", vi)

  ## breakpoints
  br <- applyfun(xx, FUN = extract_breaks, x_names = x_names, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels, extract_split_fun = extract_split, start_0 = FALSE)
  ## weird internal error handling
  ## if some errors occured within extract_split, extract_varid_, or extract_breaks
  ## drop these resamples from vs and br
  tmp <- which(sapply(br, function(x) any(sapply(x, class) == "try-error")))
  if(length(tmp)) {
    vi_mat <- vi_mat[-tmp, ]
    br[tmp] <- NULL
    warning("Due to internal coercion errors, only the results of ", B - length(tmp), "resamples are returned.")
  }
  br <- lapply(x_names, function(n) {
    if (x_classes[n] == "factor") {
      do.call("rbind", lapply(br, "[[", n))
    } else {
      unlist(lapply(br, "[[", n))
    }
  })
  names(br) <- x_names

  ## build stabletree object
  rval <- list(
    call = call,
    B = B - length(tmp),
    sampler = sampler,
    vs0 = vs0,
    br0 = br0,
    vs = vi_mat,
    br = add_levels(br, x_classes = x_classes, x_levels = x_levels, x_nlevels = x_nlevels),
    classes = x_classes,
    trees = NULL # not really usable, so set trees to NULL
  )
  class(rval) <- "stabletree"
  return(rval)
}

### ----------------------------------------------------------------------------
