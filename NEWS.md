# stablelearner 0.1-3

* Small improvements for CRAN checks.


# stablelearner 0.1-2

* Changed default sampling method in `stabletree()` from `bootstrap()` to
  `subsampling()` with default fraction of `v = 0.632`.

* `as.stabletree()` coercion generic added which allows to coerce a
  `randomForest` (_randomForest_ package), `RandomForest` (_party_ package), 
  `cforest` (_partykit_ package) or `ranger` (_ranger_ package) to a
  `stabletree` object.

* Added a vignette on the variable and cutpoint selection analysis of
  random forests.


# stablelearner 0.1-1

* Project `stablelearner` has been launched and a stable version of the 
  package has been uploaded to CRAN.

* `stability()` is available to estimate the stability of the results
  from a given supervised statistical learning method.

* `stabletree()` is available to estimate the stability of the results
  from recursive partitioning.
