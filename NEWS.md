# stablelearner 0.1-7

* Fix small bug in `splithalf()` sampling for `stability()`. In case of
  an uneven number of observations in the learning sample some small overlap
  in the two splitted halves would occur. This is avoided now by dropping
  one random observation now and restricting both halves to be always of
  the same size. Reported by Constantin Wiegand along with two further
  small improvements: Warnings about lack of variance in predictions on
  learning samples and avoiding manual triggering of the garbage collector.


# stablelearner 0.1-6

* Improve non-anchored links in manual pages (prompted by CRAN).


# stablelearner 0.1-5

* The `as.stabletree()` method for `RandomForest` objects (_party_ package)
  is now registered as an S3 (rather than S4) method.

* Environment `.stabEnv` is only used internally by `stabletree()` and
  hence not exported anymore in `NAMESPACE`.


# stablelearner 0.1-4

* Z takes over maintenance from MP.

* Fix checks for `sampler` specifications where a logical of length 2 is
  now correctly aggregated with `all()`.


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
