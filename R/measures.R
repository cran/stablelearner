## Classification
clagree <- function()
{
  measure <- function(p, q) {
    p <- apply(p, 1, which.max)
    q <- apply(q, 1, which.max)
    mean(p==q)
  } 
  list(name    = "Class agreement", 
       measure = measure, 
       classes = c("ordered", "factor", "logical"),
       range   = list(lower = 0, upper = 1),
       reverse = NULL)
}

ckappa <- function()
{
  measure <- function(p, q) {
    p <- apply(p, 1, which.max)
    q <- apply(q, 1, which.max)
    H <- table(p, q)
    N <- sum(H)
    p0 <- sum(diag(H))/N
    pc <- sum(colSums(H) * rowSums(H)) / N^2
    rval <- (p0-pc)/(1-pc)
  } 
  list(name    = "Cohens kappa", 
       measure = measure, 
       classes = c("ordered", "factor", "logical"),
       range   = list(lower = -1, upper = 1),
       reverse = NULL)
}

bdist <- function()
{
  measure <- function(p, q) mean(1 - rowSums(sqrt(p*q)))
  list(name    = "Bhattacharyya distance",
       measure = measure, 
       classes = c("ordered", "factor", "logical"),
       range   = list(lower = 0, upper = 1),
       reverse = function(x) 1 - x)
} 

tvdist <- function()
{
  measure <- function(p, q) mean(rowSums(abs(p - q))/2)
  list(name    = "Total variation distance", 
       measure = measure, 
       classes = c("ordered", "factor", "logical"),
       range   = list(lower = 0, upper = 1),
       reverse = function(x) 1 - x)
} 

hdist <- function() 
{
  measure <- function(p, q) mean(sqrt(rowSums((sqrt(p) - sqrt(q))^2))/sqrt(2))
  list(name    = "Hellinger distance", 
       measure = measure, 
       classes = c("ordered", "factor", "logical"),
       range   = list(lower = 0, upper = 1),
       reverse = function(x) 1 - x)
} 

jsdiv <- function(base = 2)
{
  measure <- function(p, q) {
    m <- (p + q) / 2
    d1 <- p*(log(p/m, base = base))
    d1[p==0] <- 0
    d2 <- p*(log(p/m, base = base))
    d2[p==0] <- 0
    return(mean(rowSums(d1+d2)/2))
  }
  list(name    = "Jensen-Shannon divergence", 
       measure = measure, 
       classes = c("ordered", "factor", "logical"),
       range   = list(lower = 0, upper = 1),
       reverse = function(x) 1 - x)
} 

## Regression
edist <- function() {
  measure <- function(x, y) 
    return(sqrt(sum((x-y)^2)))
  list(name    = "Euclidean distance",
       measure = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = 0, upper = Inf),
       reverse = function(x) -x)
}

msdist <- function()
{
  measure <- function(x, y) mean((x - y)^2)
  list(name    = "Mean squared distance", 
       measure = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = 0, upper = Inf),
       reverse = function(x) -x)
} 

rmsdist <- function()
{
  measure <- function(x, y) sqrt(mean((x - y)^2))
  list(name    = "Root mean squared distance", 
       measure = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = 0, upper = Inf),
       reverse = function(x) -x)
} 

madist <- function()
{
  measure <- function(x, y) mean(abs(x - y))
  list(name    = "Mean absolute distance", 
       measure = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = 0, upper = Inf),
       reverse = function(x) -x)
} 

qadist <- function(p = 0.95)
{
  measure <- function(x, y) quantile(abs(x - y), probs = p)
  list(name    = "Quantile of absolute deviation", 
       measure = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = 0, upper = Inf),
       reverse = function(x) -x)
} 

cprob <- function(kappa = 0.1)
{
  measure <- function(x, y) mean(abs(x - y) < kappa)
  list(name    = "Coverage probability of absolute deviation", 
       measure  = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = 0, upper = 1),
       reverse = function(x) x)
} 

ccc <- function() {
  measure <- function(x, y) {
    n <- length(x)
    s2.x <- var(x) * (n-1)/n
    s2.y <- var(y) * (n-1)/n
    s.xy <- var(x, y) * (n-1)/n
    mu.x <- mean(x)
    mu.y <- mean(y)
    return(2 * s.xy / (s2.x + s2.y + (mu.x - mu.y)^2))
  }
  list(name    = "Concordance correlation coefficient",
       measure = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = 0, upper = 1),
       reverse = NULL)
}

pcc <- function() {
  measure <- function(x, y) return(cor(x, y, method = "pearson"))
  list(name    = "Pearson correlation coefficient",
       measure = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = 0, upper = 1),
       reverse = NULL)
}

cosine <- function() {
  measure <- function(x, y) 
    return(sum(x*y) / (sqrt(sum(x^2))*sqrt(sum(y^2))))
  list(name    = "Cosine similarity",
       measure = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = -1, upper = 1),
       reverse = NULL)
}

rbfkernel <- function() {
  measure <- function(x, y, sigma = 1)
    return(exp(-sum((x-y)^2) / (2 * sigma^2)))
  list(name    = "Radial basis function kernel",
       measure = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = 0, upper = 1),
       reverse = NULL)
}

tanimoto <- function() {
  measure <- function(x, y) 
    return(sum(x*y) / (sum(x^2) + sum(y^2) - sum(x*y)))
  list(name    = "Tanimoto coefficient",
       measure = measure,
       classes = c("numeric", "integer"),
       range   = list(lower = -1/3, upper = 1),
       reverse = NULL)
}