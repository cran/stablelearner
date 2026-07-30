# Stability Assessment of Statistical Learning Methods in R

## Overview

The R package [stablelearner](https://zeileis.codeberg.page/stablelearner/) provides:

* Stability assessment for tree learners: `stabletree` and accompanying methods,
  including coercion functions for various random forest objects to `stabletree` objects.

* Stability assessment for supervised statistical learners in general: `stability` and
  accompanying methods, including a broad range of similarity measures for both
  classification and regression problems.


## References

Philipp M, Zeileis A, Strobl C (2016).
  "A Toolkit for Stability Assessment of Tree-Based Learners."
  In Colubi A, Blanco A, Gatu C (eds.),
  _Proceedings of COMPSTAT 2016 - 22nd International Conference on Computational Statistics_, 315-325.
  ISBN 978-90-73592-36-0.
  Preprint available at <https://EconPapers.RePEc.org/RePEc:inn:wpaper:2016-11>

Philipp M, Rusch T, Hornik K, Strobl C (2018).
  "Measuring the Stability of Results from Supervised Statistical Learning."
  _Journal of Computational and Graphical Statistics_, **27**(4), 685-700.
  [doi:10.1080/10618600.2018.1473779](https://doi.org/10.1080/10618600.2018.1473779)


## Installation

The stable version of `stablelearner` is available from
[CRAN](https://CRAN.R-project.org/package=stablelearner):

``` r
install.packages("stablelearner")
```

The latest development version can be installed from
[R-universe](https://zeileis.R-universe.dev/stablelearner):

``` r
install.packages("stablelearner", repos = "https://zeileis.R-universe.dev")
```


## License

The package is available under the
[General Public License version 3](https://www.gnu.org/licenses/gpl-3.0.html)
or [version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
