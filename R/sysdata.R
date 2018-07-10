LearnerList <- list(
  "constparty" = list(
    class = "constparty",
    package = "partykit",
    predict = function(x, newdata, yclass = NULL) {
      type <- ifelse(match(yclass, c("ordered", "factor"), nomatch = FALSE), "prob", "response")
      predict(x, newdata = newdata, type = type)
    }
  ),
  "cforest" = list(
    class = "cforest",
    package = "partykit",
    predict = function(x, newdata, yclass = NULL) {
      type <- ifelse(match(yclass, c("ordered", "factor"), nomatch = FALSE), "prob", "response")
      predict(x, newdata = newdata, type = type)
    }
  ),
  "rpart" = list(
    class = "rpart",
    package = "rpart",
    predict = function(x, newdata, yclass = NULL) {
      type <- ifelse(match(yclass, c("ordered", "factor"), nomatch = FALSE), "prob", "vector")
      predict(x, newdata = newdata, type = type)
    }
  ),
  "J48" = list(
    class = "J48",
    package = "RWeka",
    predict = function(x, newdata, yclass = NULL) {
      type <- ifelse(match(yclass, c("ordered", "factor"), nomatch = FALSE), "prob", "class")
      predict(x, newdata = newdata, type = type)
    }
  ),
  "C5.0" = list(
    class = "C5.0",
    package = "C50",
    predict = function(x, newdata, yclass = NULL) {
      type <- ifelse(match(yclass, c("ordered", "factor"), nomatch = FALSE), "prob", "class")
      predict(x, newdata = newdata, type = type)
    }
  ),
  "tree" = list(
    class = "tree",
    package = "tree",
    predict = function(x, newdata, yclass = NULL) {
      predict(x, newdata = newdata, type = "vector")
    }
  ),
  "randomForest" = list(
    class = "randomForest",
    package = "randomForest",
    predict = function(x, newdata, yclass = NULL) {
      type <- ifelse(match(yclass, c("ordered", "factor"), nomatch = FALSE), "prob", "response")
      predict(x, newdata = newdata, type = type)
    }
  ),
  "svm" = list(
    class = "svm",
    package = "e1071",
    predict = function(x, newdata, yclass = NULL) {
      if(match(yclass, c("ordered", "factor"))) {
        attr(predict(x, newdata = newdata, probability = TRUE), "probabilities")
      } else {
        predict(x, newdata = newdata)
      }
    }
  ),
  "ksvm" = list(
    class = "ksvm",
    package = "kernlab",
    predict = function(x, newdata, yclass = NULL) {
      if(match(yclass, c("ordered", "factor"))) {
        predict(x, newdata = newdata, type = "probabilities")
      } else {
        predict(x, newdata = newdata)
      }
    },
    update = function(x, data = NULL, weights = NULL) {
      if(!is.null(weights)) stop("Weights are not supported for this class. Please use data argument.")
      call <- as.list(x@kcall)[-1]
      call$x <- formula(x@terms)
      call$data <- data
      do.call("ksvm", call)
    }
  ),
  "mboost" = list(
    class = "mboost",
    package = "mboost",
    predict = function(x, newdata, yclass = NULL) {
      p <- predict(x, newdata = newdata, type = "response")
      return(p)
    },
    update = function(x, data = NULL, weights = NULL) {
      call <- as.list(getCall(x))[-1]
      call$data <- data
      call$weights <- weights
      do.call("mboost", call)
    }
  ),
  "boosting" = list(
    class = "boosting",
    package = "adabag",
    predict = function(x, newdata, yclass = NULL) {
      predict(x, newdata = newdata)$prob
    }
  ),
  "adaboost" = list(
    class = "adaboost",
    package = "fastAdaboost",
    predict = function(x, newdata, yclass = NULL) {
      predict(x, newdata = newdata)$prob
    }
  ),
  "gbm" = list(
    class = "gbm",
    package = "gbm",
    predict = function(x, newdata, yclass = NULL) {
      predict(x, newdata = newdata, n.trees = x$n.trees, type = "response")[,,1]
    },
    terms = function(x) x$Terms
  ),
  "nnet" = list(
    class = "nnet",
    package = "nnet",
    predict = function(x, newdata, yclass = NULL) {
      predict(x, newdata = newdata, type = "raw")
    }
  ),
  "multinom" = list(
    class = "multinom",
    package = "nnet",
    predict = function(x, newdata, yclass = NULL) {
      predict(x, newdata = newdata, type = "probs")
    }
  ),
  "lda" = list(
    class = "lda",
    package = "MASS",
    predict = function(x, newdata, yclass = NULL)  {
      predict(x, newdata = newdata)$posterior
    }
  ),
  "lm" = list(
    class = "lm",
    package = "stats",
    predict = function(x, newdata, yclass = NULL) {
      predict(x, newdata = newdata)
    }
  ),
  "glm" = list(
    class = "glm",
    package = "stats",
    predict = function(x, newdata, yclass = NULL) {
      predict(x, newdata = newdata, type = "response")
    }
  ),
  "glmrob" = list(
    class = "glmrob",
    package = "robustbase",
    predict = function(x, newdata, yclass = NULL) {
      predict(x, newdata = newdata, type = "response")
    }
  ),
  "glmnet.formula" = list(
    class = "glmnet.formula",
    package = "glmnetUtils",
    predict = function(x, newdata, yclass = NULL) {
      p <- predict(x, newdata = newdata, type = "response")
      if(getCall(x)$family == "multinomial") p <- p[,,1]
      return(p)
    },
    update = function(x, data = NULL, weights = NULL) {
      call <- as.list(getCall(x))[-1]
      call$data <- data
      call$weights <- weights
      do.call("glmnet.formula", call, envir = asNamespace("glmnetUtils"))
    }
  )
)

## getCal function for constparties
getCall.cforest <- function(x, ...) x$info$call

## getCall function for class ksvm
# getCall.ksvm <- function(object) {
#   call <- as.list(x@kcall)
#   call[[1]] <- as.symbol("ksvm")
#   call$x <- formula(x@terms)
#   call$data = NULL
#   return(as.call(call))
# }

# save("LearnerList", file = "sysdata.rda")