## -----------------------------------------------------------------------------
## Learner
## -----------------------------------------------------------------------------

## The function stability uses update() which uses getCall() to get x$call 
## to refit the model with new data sets. Therefore, only learners that 
## store the call in their resulting object will work. Otherwise, a function 
## getCall must be provided for the that class. For details, see ?update.
# 
getLearner <- function(x) {
  LearnerList <- get("LearnerList", envir = .GlobalEnv)
  if(inherits(x, names(LearnerList))) {
    pos <- inherits(x, names(LearnerList), which = TRUE)
    pos[pos==0L] <- NA
    LearnerList[[which.min(pos[pos>0L])]]
  } else {
    stop("Learner of class: ", class(x), " not in LearnerList. See ?LearnerList for help.")
  }
}

addLearner <- function(x) {
  eval.parent(substitute(LearnerList <- c(list(x), LearnerList)))
  eval.parent(substitute(names(LearnerList)[1] <- x$class))
  #eval.parent(LearnerList <- (names(LearnerList)[1] <- x$class))
}
