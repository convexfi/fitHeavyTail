## Changes in fitHeavyTail version 0.1.3 (2022-4-14)

* New method for skewed t distributions: fit_mvst()

* fit_mvt() and fit_mvst(): Now the bounds for nu estimation can be set as a global option, e.g.: options(nu_min = 4.2).

* Fixed description regarding covariance matrix for Cauchy distribution.

* fit_mvt(): It accepts weights as argument to weight differently the samples (as opposed to uniform weights).

* fit_mvt(): Many more methods to estimate nu iteratively (via argument nu_iterative_method).

* fit_mvt(): New argument scale_minMSE to include a correction factor in the covariance matrix for minimum MSE (still in development).


## Changes in fitHeavyTail version 0.1.2 (2020-1-7)

* Vignette revised: detailed descriptions of the algorithms included.

* Comparison with additional existing benchmark sn::selm() included in the vignette.

* Now the three fitting functions also return the number of iterations and elapsed cpu_time.

* Significant revision of the fitting function fit_mvt(); in particular:

  - the nu_target for the estimation of nu has been removed since it was not effective;
  - several new options for the initial value of nu or fixed value of nu have been included; and
  - improved and more robust estimation of nu at each EM iteration.

* Function fit_mvt() now allows the choice (via the argument na_rm) to drop the observations with NAs 
  or impute them.


## Changes in fitHeavyTail version 0.1.1 (2019-11-22)

* Initial release is on CRAN.

* It includes three functions for heavy tails fitting: fit_mvt(), fit_Tyler(), and fit_Cauchy().

* Vignette illustrates its use and comparison with existing packages.

* Tests are included.

* fit_mvt() can deal with NAs and a factor model structure on the covariance matrix.


## Changes in fitHeavyTail version 0.1.0 (2017-11-04)

* Initial bare-bone implementation (not publicly released).
