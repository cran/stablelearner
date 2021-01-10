bootstrap <- function(B = 500, v = 1) {
  sampfun <- function(n) replicate(B, sample(1L:n, floor(v * n), replace = TRUE))
  list(method = sprintf("Bootstrap sampling with %1.1f%% data", 100 * v), 
       sampler = sampfun)
}

# bootstrap <- function(B = 500) {
#   sampfun <- function(n) replicate(B, sample(1L:n, n, replace = TRUE))
#   list(method = "Bootstrap sampling", sampler = sampfun)
# }

subsampling <- function(B = 500, v = 0.632) {
  sampfun <- function(n) replicate(B, sample(1L:n, floor(v * n), replace = FALSE))
  list(method = sprintf("Subsampling with %1.1f%% data", 100 * v), 
    sampler = sampfun)
}

samplesplitting <- function(k = 5) {
  sampfun <- function(n) {
    ret <- matrix(NA, ncol = ceiling(n/k), nrow = k)
    ret[1L:n] <- sample(1L:n)
    t(ret)
  }
  list(method = sprintf("%s-fold sample splitting", k), sampler = sampfun)
}

jackknife <- function(d = 1, maxrep = 5000) {
  sampfun <- function(n) {
    if (choose(n, d) > maxrep) 
      stop("Maximum number of repetitions allowed exceeded. Reduce d!")
    apply(utils::combn(1L:n, d), 2, function(x) (1L:n)[-x])
  }
  list(method = sprintf("Leave-%s-out jackknife", d), sampler = sampfun)
} 

splithalf <- function (B = 500)
{
  sampfun <- function(n) replicate(B, sample(1L:n, floor(1/2 * n), replace = FALSE))
  list(method = "Split-half sampling", sampler = sampfun)
}
