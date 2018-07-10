## -----------------------------------------------------------------------------
## Stability
## -----------------------------------------------------------------------------

stab_control <- function(B = 500, measure = list(tvdist, ccc),
                         sampler = "bootstrap", evaluate = "OOB", 
                         holdout = 0.25, seed = NULL, na.action = na.exclude, 
                         savepred = TRUE, silent = TRUE, ...) {
  
  ## process evaluaton method
  if(!is.null(evaluate)) {
    evaluate <- match.arg(toupper(evaluate), 
                          choices = c("OOB", "OOS", "ALL"), 
                          several.ok = TRUE)
  }
  
  ## process sampler
  if (is.character(sampler))
    sampler <- get(sampler, mode = "function", envir = parent.frame())
  if (is.function(sampler)) {
    samp <- try(sampler(B = B, ...), silent = TRUE)
    if (!inherits(samp, "try-error") && is.list(samp) && c("method", "sampler") %in% 
        names(samp)) {
      sampler <- samp
    } else {
      sampler <- list(method = "User-defined", sampler = sampler)
    }
  }
  
  if(is.null(measure)) stop("Please choose an appropriate measure!")
  
  ## process measure
  if(is.function(measure)) measure <- list(measure)
  measure <- lapply(measure, function(mfun) {
    if (is.character(mfun)) 
      mfun <- get(mfun, mode = "function", envir = parent.frame())
    if (is.function(mfun)) {
      samp <- try(mfun(), silent = TRUE)
      if (!inherits(samp, "try-error") && is.list(samp) && c("name", "measure") %in% 
          names(samp)) {
        mfun <- samp
      } else {
        mfun <- list(name = "User-defined", measure = mfun, classes = NULL, 
                     range = NULL, reverse = NULL)
      }
    }
    mfun
  })
  
  rval <- list(
    B = B,
    sampler = sampler,
    evaluate = evaluate,
    measure = measure,
    na.action = na.action,
    holdout = holdout,
    seed = seed,
    savepred = savepred,
    silent = silent
  )
  
  return(rval)
  
}

stability <- function(x, ..., data = NULL, control = stab_control(), 
                      weights = NULL, applyfun = NULL, cores = NULL, 
                      names = NULL) {
  
  if (!inherits(x, "list")) x <- list(x, ...)
  
  ## no sampling when data argument is a function
  if(is.function(data))
    control$sampler <- control$evaluate <- NULL
  
  ## determine B when weights are submitted
  if(length(dim(weights))==3L) {
    control$B <- dim(weights)[2]
    control$sampler <- control$evaluate <- NULL
  }
  
  ## get learner infos
  nlearner <- length(x)
  m <- lapply(x, stablelearner::getLearner)
  
  rng <- RNGkind()
  
  rval <- lapply(1L:nlearner, function(k) {
    ## process random number generator
    if(!is.null(control$seed))
      set.seed(control$seed, kind = "L'Ecuyer-CMRG")
    ## process stability computation
    stability_internal(x = x[[k]], learner = m[[k]], data = data, 
                       weights = weights, control = control, 
                       applyfun = applyfun, cores = cores)
  })
  
  ## reset random number generator
  set.seed(NULL, kind = rng[1])
  
  if(is.null(names)) 
    names(rval) <- 1L:nlearner 
  else names(rval) <- names
  
  class(rval) <- "stablelearnerList"
  
  if(nlearner==1L)
    return(rval[[1]]) 
  else return(rval)
  
}

process_sampler <- function(n, sampler, evaluate, holdout, weights) {
  if(is.null(weights)|isTRUE(weights)) {
    allobs <- 1L:n
    if(sampler$method=="User-defined") {
      rout <- sampler$sampler(n = n)
      ls1 <- rout[[1L]]
      ls2 <- rout[[2L]]
      es <- rout[[3L]]
      # stop("Currently not implemented. Please use weights argument!")
    } else {
      if(evaluate=="OOS") {
        test <- sample(allobs, floor(holdout*n))
        train  <- setdiff(allobs, test)
      } else {
        train <- allobs
      }
      r <- length(train)
      ## sampling
      ls1 <- sampler$sampler(r)
      B <- ncol(ls1)
      ls1 <- matrix(train[ls1], ncol = B)
      if(sampler$method=="Split-half sampling") {
        ls2 <- sapply(1L:B, function(b) setdiff(train, ls1[,b]))
      } else {
        ls2 <- matrix(train[sampler$sampler(r)], ncol = B)
      }
      ## evaluate
      es <- switch(evaluate,
                   "OOB" = lapply(1L:B, function(b) setdiff(train, union(ls1[,b], ls2[,b]))),
                   "OOS" = lapply(1L:B, function(b) test),
                   "ALL" = lapply(1L:B, function(b) train))
      if(any(sapply(es, length) == 0)) warning("Repetitions with empty evaluation samples.")
    }
    if(isTRUE(weights)) {
      ls1 <- apply(ls1, 2L, tabulate, nbins = n)
      ls2 <- apply(ls2, 2L, tabulate, nbins = n)
      es  <- sapply(es, tabulate, nbins = n)
    }
    rval <- if(isTRUE(weights)) {
      array(c(ls1, ls2, es), dim = c(dim(ls1)[1], B, 3))
    } else {
      list(ls = array(c(ls1, ls2), dim = c(dim(ls1)[1], B, 2)), es = es)
    }
  } else {
    ## FIXME: currently only case weights implemented!
    rval <- round(weights)
  }
  return(rval)
}

stability_internal <- function(x, learner, data, weights, control, 
                               applyfun, cores, ...) {
  
  B        <- control$B
  measure  <- control$measure
  sampler  <- control$sampler
  evaluate <- control$evaluate
  holdout  <- control$holdout
  na.fun   <- control$na.action
  savepred <- control$savepred
  silent   <- control$silent
  
  prfun <- if(!is.null(learner$predict)) learner$predict else predict
  upfun <- if(!is.null(learner$update)) learner$update else update
  trfun <- if(!is.null(learner$terms)) learner$terms else terms
  
  ## extract information from call
  call <- try(getCall(x), silent = TRUE)
  if(inherits(call, "try-error")) {
    # stop("Learner objects that do not provide the model call are not supported.")
    sfit <- wfit <- dfit <- NULL
  } else {
    sfit <- call$subset
    wfit <- call$weights
    dfit <- call$data
  }
  
  # ## process call to receive response class
  # ## FIXME: probably not working for all learner in general
  # call <- getCall(x)
  # f <- formula(x)
  # resp <- f[[2]]
  # if(!is.name(resp)) stop("response is not a name")
  # yname <- as.character(resp)
  # mf <- call
  # mf <- mf[c(1L, match(c("formula", "data"), names(mf), 0L))]
  # mf[[1L]] <- as.name("model.frame")
  # mf[[2L]] <- f
  # mf <- eval(mf, parent.frame())
  # yclass <- class(mf[,yname])
  
  ## check if weights are supported
  # wsup <- "weights" %in% formalArgs(as.character(call[[1L]]))
  # if(!is.null(weights)&!wsup) stop("Weights not supported by the learner! Please use sampler argument instead.")
  
  # call <- getCall(x)
  # if (is.null(data)) data <- eval(call$data, parent.frame())  #FIXME# more elegant default?
  # n <- nrow(data)
  
  ## get envoronment of x
  env <- try(environment(trfun(x)), silent = TRUE)
  if(inherits(env, "try-error")) env <- NULL
  
  ## data extraction and resampling
  if(!is.function(data)) {
    ## get data
    if (is.null(data)) {
      if (is.null(dfit)) {
        ## there is no data object
        data <- NULL
      } else
        ## get local copy of data object from where x was generated
        data <- eval(dfit, envir = env, enclos = parent.frame())
      if(!is.null(sfit)) {
        sfit <- eval(sfit, envir = env, enclos = parent.frame())
        data <- subset(data, subset = sfit)
      }
    }
    
    if(is.null(data)) {
      stop("Could not extract the data from the model object. Please provide 
           the data via the 'data' argument.")
    }
    
    ## apply na.action
    data <- na.fun(data)
    
    ## get sample size
    n <- NROW(data)
    
    samp <- process_sampler(n, sampler, evaluate, holdout, weights)
    
    ## are weights supported?
    if(!inherits(call, "try-error")) {
      fitfun <- get(as.character(call[[1L]]), mode = "function", 
                    envir = asNamespace(learner$package))
      wsup <- "weights" %in% formalArgs(fitfun)
    } else wsup <- FALSE
    
  }
  
  ## extract names and class of response variable via terms
  yclass <- yname <- NULL
  try({
    tr <- trfun(x)
    yi <- attr(tr, "response")
    yclass <- attr(tr, "dataClasses")[yi]
    yname <- names(yclass)
  }, silent = TRUE)
  
  ## if the class could not be extracted from terms
  if(inherits(yclass, "try-error") | is.null(yclass) | is.null(yname)) {
    fmla <- eval(call$formula, envir = env, enclos = parent.frame())
    yname <- all.vars(fmla)[1L]
    yval <- if(is.function(data)) {
      do.call("data", list(...))[,yname] 
    } else { 
      data[,yname] 
    }
    yclass <- class(yval)
  }
  
  ## extract names and class of response variable
  # fmla <- eval(call$formula, envir = env, enclos = parent.frame())
  # mf <- model.frame(fmla, data = data)
  # mf <- try(model.frame(tr, data = data))
  # if(inherits(mf, "try-error")) stop("Could not extract terms from fitted model object.")
  # yi     <- attr(tr, "response")
  # yclass <- class(mf[,yi])
  # yname  <- colnames(mf)[yi]
  
  ## exclude measures that are not valid for the response class of the model
  measure <- lapply(measure, function(x) {
    if(yclass %in% x$classes) x else NULL
  })
  measure <- measure[!sapply(measure, is.null)]
  control$measure <- measure
  if(length(measure)==0L) stop("None of the measures cann be applied for class ", yclass, "!")
  
  ## facilitate parallelization
  if (is.null(applyfun)) {
    applyfun <- if (is.null(cores)) {
      lapply
    } else {
      function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
    }
  }
  
  ## generate function for repetitions
  repeatfun <- function(b) {
    
    rval <- NULL
    
    try({
      
      if(!is.function(data)) {
        if(is.null(weights)) {
          ## generate subsamples
          ls1 <- data[samp$ls[,b,1L], , drop = FALSE]
          ls2 <- data[samp$ls[,b,2L], , drop = FALSE]
          es  <- data[samp$es[[b]], , drop = FALSE]
          ew  <- rep(1L, nrow(es))
        } else {
          ## extract weight
          lw1 <- samp[,b,1L]
          lw2 <- samp[,b,2L]
          ew  <- samp[,b,3L]
          es  <- data[ew>0L,]
          ew  <- ew[ew>0L]
        }
      } else {
        ls1 <- do.call("data", list(...))
        ls2 <- do.call("data", list(...))
        es  <- do.call("data", list(...))
        ew  <- rep(1L, nrow(es))
      }
      
      if(nrow(es)>0L) {
        
        ## update results
        if(is.null(weights)) {
          x1 <- upfun(x, data = ls1)
          x2 <- upfun(x, data = ls2)
        } else {
          x1 <- upfun(x, weights = lw1)
          x2 <- upfun(x, weights = lw2)
        }
        
        ## predict from both results
        p1 <- prfun(x1, newdata = na.fun(es), yclass = yclass)
        p2 <- prfun(x2, newdata = na.fun(es), yclass = yclass)
        
        ## make matrix
        if(yclass %in% c("ordered", "factor") & NCOL(p1)<2) {
          p1 <- cbind(1-p1, p1)
          p2 <- cbind(1-p2, p2)
          colnames(p1) <- levels(data[,yname])
          colnames(p2) <- levels(data[,yname])
        } else {
          p1 <- as.matrix(p1)
          p2 <- as.matrix(p2)
        }
        
        ## correct order of factors in matrix
        if(yclass %in% c("ordered", "factor")) {
          c1 <- colnames(p1)
          c2 <- colnames(p2)
          if(!identical(c1, c2)) p2 <- p2[,c1]
        }
        
        ## replicate according to weights
        p1 <- p1[rep(1L:nrow(p1), times = ew), , drop = FALSE]
        p2 <- p2[rep(1L:nrow(p2), times = ew), , drop = FALSE]
        
        sval <- sapply(measure, function(mfun) mfun$measure(p1, p2))
        
        rval <- list(stab = sval, 
                     pred = if (savepred) {
                       resp <- es[,yname]
                       list(data.frame(y = resp, p1), data.frame(y = resp, p2))
                     } else NA)
        
        ## clear memory
        rm(x1, x2, p1, p2, ls1, ls2, es, ew)
        invisible(gc(FALSE))

      } else {
        stop("Nothing to predict.")
      }

    }, silent = silent)

    if(is.null(rval))
      return(list(stab = rep(NA, length(measure)), pred = NA))
    else
      return(rval)
    
  }
  
  ## process resampling and compute similarities
  xx   <- applyfun(1L:B, repeatfun)
  pval <- lapply(xx, function(x) if(is.null(x)) NA else x$pred)
  sval <- sapply(xx, function(x) if(is.null(x)) NA else x$stab)
  sval <- if(is.vector(sval)) as.matrix(sval) else t(sval)
  colnames(sval) <- sapply(measure, function(x) x$name)
  
  if(!is.function(data)) {
    ## learning size
    ## evaluation size
    ls <- if(is.null(weights)) {
      dim(samp$ls)[1]
    } else {
      c(colSums(samp[,,1L]), colSums(samp[,,2L]))
    }
    
    ## learning overlap
    lo <- sapply(1L:B, function(b) {
      if(is.null(weights)) {
        ls1 <- samp$ls[,b,1L]
        ls2 <- samp$ls[,b,2L]
      } else {
        ls1 <- samp[,b,1L]
        ls2 <- samp[,b,2L]
        ls1 <- which(ls1 > 0L)
        ls2 <- which(ls2 > 0L)
      }
      length(intersect(ls1, ls2))
    })
    
    ## evaluation size
    es <- sapply(1L:B, function(b) {
      if(is.null(weights)) {
        length(samp$es[[b]])
      } else {
        sum(samp[,b,3L])
      }
    })
    
    ## evaluation overlap
    eo <- sapply(1L:B, function(b) {
      if(is.null(weights)) {
        ls1 <- samp$ls[,b,1L]
        ls2 <- samp$ls[,b,2L]
      } else {
        ls1 <- samp[,b,1L]
        ls2 <- samp[,b,2L]
        ls1 <- which(ls1 > 0L)
        ls2 <- which(ls2 > 0L)
      }
      es <- if(is.null(weights)) samp$es[[b]] else which(samp[,b,3L]>0L)
      length(intersect(es, union(ls1, ls2)))
    })
    sampstat <- list(ls = ls, lo = lo, es = es, eo = eo)
  } else {
    sampstat <- NULL
  }
  
  res <- list(
    call = if(inherits(call, "try-error")) NULL else call,
    learner = learner,
    B = B,
    sval = sval,
    pred = pval,
    sampstat = sampstat,
    control = control
  )
  class(res) <- "stablelearner"
  
  return(res)
  
}

accuracy <- function(x, measure = "kappa", na.action = na.exclude, 
                     applyfun = NULL, cores = NULL) {
  
  if(class(x) == "stablelearner") x <- list(x)
  
  ## facilitate parallelization
  if (is.null(applyfun)) {
    applyfun <- if (is.null(cores)) {
      lapply
    } else {
      function(X, FUN) parallel::mclapply(X, FUN, mc.cores = cores)
    }
  }
  
  rval <- applyfun(x, function(m) {
    rval <- sapply(m$pred, function(rep.outer) {
      sapply(rep.outer, function(rep.inner) {
        rep.inner <- na.action(rep.inner)
        if(nrow(rep.inner)<1L) return(rep(NA, 4))
        lev <- levels(rep.inner[,1])
        cnm <- colnames(rep.inner)[-1L]
        prd <- factor(cnm[apply(rep.inner[,-1], 1, which.max)], levels = lev)
        tab <- table(rep.inner[,1], prd)
        unlist(e1071::classAgreement(tab))
      })
    }, simplify = "array")
    rval[measure,,]
  })
  
  attr(rval, "measure") <- measure
  
  return(rval)
  
}

tuner <- function(method, tunerange, ...) {
  
  m <- match.call(expand.dots = FALSE)
  
  range_df <- expand.grid(tunerange)
  
  rval <- lapply(1L:nrow(range_df), function(i) {
    return(do.call(method, c(m$..., as.list(range_df[i, , drop = FALSE]))))
  })
  
  attr(rval, "range") <- range_df
  
  return(rval)
  
}

## -----------------------------------------------------------------------------
## Plotting and printing
## -----------------------------------------------------------------------------

similarity_values <- function(x, reverse = TRUE)  {
  if(class(x) == "stablelearner") {
    x <- list(x)
    class(x) <- "stablelearnerList"
  }
  if(reverse) x <- reverse(x)
  rval <- sapply(x, function(xx) xx$sval, simplify = "array")
  dimnames(rval) <- list(NULL, colnames(x[[1]]$sval), names(x))
  return(aperm(rval, perm = c(1L,3L,2L)))
}

reverse <- function(x, ...) {
  if(class(x) == "stablelearner") x <- list(x)
  rval <- lapply(x, function(xx) {
    sval <- xx$sval
    rvalues <- sapply(seq(ncol(sval)), function(k) {
      rfun <- xx$control$measure[[k]]$reverse
      if(!is.null(rfun)) 
        rfun(sval[,k])
      else sval[,k]
    })
    cn <- colnames(sval)
    rev <- sapply(xx$control$measure, function(z) !is.null(z$reverse))
    cn[rev] <- paste0(cn[rev], " (reversed)")
    colnames(rvalues) <- cn
    xx$sval <- rvalues
    return(xx)
  })
  return(rval)
}

summary.stablelearner <- function(object, ...) {
  obj_list <- list(object)
  names(obj_list) <- object$learner$class
  class(obj_list) <- "stablelearnerList"
  summary.stablelearnerList(obj_list, ...)
}

summary.stablelearnerList <- function(object, ..., reverse = TRUE, 
                                      probs = c(0.05, 0.25, 0.5, 0.75, 0.95),
                                      digits = 3, names = NULL)
{
  if(!is.null(names)) names(object) <- names
  nlearner <- length(object)
  ans <- list()
  ans$digits <- digits
  ans$B <- object[[1]]$B
  ans$sampmethod <- object[[1]]$control$sampler$method
  ans$evalmethod <- object[[1]]$control$evaluate
  ans$sampstat <- NULL
  if(!is.null(object[[1]]$sampstat)) {
    tmp <- sapply(object, function(x) x$sampstat$ls)
    ls <- c(mean(tmp), min(tmp), max(tmp))
    tmp <- sapply(object, function(x) x$sampstat$lo)
    lo <- c(mean(tmp), min(tmp), max(tmp))
    tmp <- sapply(object, function(x) x$sampstat$es)
    es <- c(mean(tmp), min(tmp), max(tmp))
    tmp <- sapply(object, function(x) x$sampstat$eo)
    eo <- c(mean(tmp), min(tmp), max(tmp))
    ans$sampstat <- rbind(ls, lo, es, eo)
    rownames(ans$sampstat) <- c(
      "Learning sample size",
      "Learning overlap",
      "Evaluation sample size",
      "Evaluation overlap"
    )
    colnames(ans$sampstat) <- c("avg", "min", "max")
  }
  ans$data <- object[[1]]$data
  ans$call <- lapply(object, function(x) x$call)
  simv <- similarity_values(object, reverse = reverse)
  qtab <- apply(simv, c(2,3), quantile, probs = probs, na.rm = TRUE)
  qtab <- array(qtab, dim = c(length(probs), dim(simv)[-1]))
  dimnames(qtab) <- c(list(sprintf("%i%%", probs*100)), dimnames(simv)[2:3])
  ans$qTable <- aperm(qtab, perm = c(2,1,3))
  ans$missings <- sapply(object, function(x) max(colSums(is.na(x$sval))), simplify = "array")
  names(ans$missings) <- names
  class(ans) <- "summary.stablelearnerList"
  ans
}

print.summary.stablelearnerList <- function(x, ...)
{
  cat("\nCall:\n")
  sapply(x$call, function(y) print(y))
  if(is.null(x$sampstat)) {
    cat("\nData-generating process: ")
    cat(as.character(x$data), "\n")
  } else {
    cat("\nSampler:\n")
    cat("B =", x$B, "\n")
    cat("Resampling method =", x$sampmethod, "\n")
    cat("Evaluation method =", x$evalmethod, "\n")
    cat("\nSampling summary:\n")
    print(x$sampstat)
  }
  # cat("   = ")
  # cat(sprintf("%03.2f obs [%i, %i]\n", x$ls[1], x$ls[2], x$ls[3]))
  # cat("     = ")
  # cat(sprintf("%03.2f obs [%i, %i]\n", x$lo[1], x$lo[2], x$lo[3]))
  # cat(" = ")
  # cat(sprintf("%03.2f obs [%i, %i]\n", x$es[1], x$es[2], x$es[3]))
  # cat("        = ")
  # cat(sprintf("%03.2f obs [%i, %i]\n", x$eo[1], x$eo[2], x$eo[3]))
  cat("\nSimilarity summary:\n\n")
  print(format(x$qTable, digits = x$digits), quote = FALSE)
  if(any(x$missings>0L)) {
    cat("Missing values (caused by errors in repetitions):\n\n")
    print(x$missings)
  }
}

print.stablelearner <- function(x, ...) {
  obj <- list(x)
  names(obj) <- x$learner$class
  print.stablelearnerList(obj, ...)
}

print.stablelearnerList <- function(x, ...) {
  ans <- summary.stablelearnerList(x, ...)
  cat("\nCall:\n")
  sapply(ans$call, function(y) print(y))
  if(!is.null(ans$sampstat)){
    cat("\nSampler:\n")
    cat("B =", ans$B, "\n")
    cat("Resampling method =", ans$sampmethod, "\n")
    cat("Evaluation method =", ans$evalmethod, "\n")
  } else {
    cat("\nData-generating process: ")
    cat(as.character(ans$data), "\n")
  }
}

boxplot.stablelearner <- function(x, ...) {
  obj <- list(x)
  names(obj) <- x$learner$class
  boxplot.stablelearnerList(obj, ...)
}

boxplot.stablelearnerList <- 
  function(x, ..., main = NULL, xlab = NULL,  ylab = NULL, 
           reverse = TRUE)
  {
    
    sval <- similarity_values(x, reverse = reverse)
    nplt <- dim(sval)[3]
    
    if(is.null(ylab)) ylab <- dimnames(sval)[[3]]
    
    if(nplt>1L)
      par(mfrow = n2mfrow(nplt))
    
    sapply(1L:nplt, function(x) {
      boxplot(sval[,,x], main = main, xlab = xlab, ylab = ylab[x], ...)
      return(invisible(NULL))
    })
    
    if(nplt>1L)
      par(mfrow = c(1, 1))
    
  }

# plot.stablelearner <- function(s, ...) {
#   obj <- list(s)
#   names(obj) <- s$learner$class
#   plot.stablelearnerList(obj, ...)
# }
# 
# plot.stablelearnerList <- function(obj, reverse.values = TRUE, 
#                                    col = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) {
#   require("ggplot2", quietly = TRUE)
#   require("reshape", quietly = TRUE)
#   obj <- if(reverse.values) reverse.stablelearnerList(obj) else obj
#   plotdata <- lapply(1L:length(obj), function(k)
#     data.frame(names(obj)[k], melt(obj[[k]]$sval)))
#   plotdata <- do.call("rbind", plotdata)
#   colnames(plotdata) <- c("learner", "rep", "measure", "sval")
#   plotdata$measure <- factor(plotdata$measure, colnames(obj[[1]]$sval))
#   ggplot(plotdata, aes(learner, sval)) + 
#     geom_violin(aes(fill = learner)) + geom_boxplot(width = 0.1) + 
#     scale_fill_manual(name = "Learner", values = col[1L:length(obj)], labels = names(obj)) + 
#     facet_wrap(~ measure, ncol = 3) + 
#     xlab(label = "") + 
#     ylab(label = 
#            if(reverse.values)
#              "Similarity"
#          else 
#            "Interpretation corresponding to original definition of the measure(s)")
# }

## -----------------------------------------------------------------------------
## DGPs
## -----------------------------------------------------------------------------

dgp_twoclass <- function(n = 100, p = 4, noise = 16, rho = 0, b0 = 0, 
                         b = rep(1, p), fx = identity) {
  d <- p + noise
  S <- diag(d)
  S[1L:p,1L:p] <- rho^(toeplitz(1L:p)-1)
  X <- MASS::mvrnorm(n, mu = rep(0, d), Sigma = S)
  colnames(X) <- paste0("x", 1L:d)
  xx <- fx(X[, 1L:p, drop = FALSE])
  pi <- plogis(b0 + xx %*% b)
  resp <- sapply(pi, rbinom, n = 1, size = 1)
  ans <- data.frame(class = factor(resp, labels = c("A", "B")), X)
  attr(ans, "signal") <- 1L + c(1L:p)
  attr(ans, "prob") <- pi
  class(ans) <- c("data.frame", "dgp")
  return(ans)
}
