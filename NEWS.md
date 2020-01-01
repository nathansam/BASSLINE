## BASSLINE 0.0.0.9004 (YYYY-MM-DD)

### Patch notes 

* Changed N = 1000 to N = 100 (and thin to 2 and burn to 20) in example for
`LML_LEP` to reduce running times. 
* Added additonal tests

## BASSLINE 0.0.0.9003 (2019-12-28)

The main focus of this update has been to improve documentation and ease of use. 

### Patch notes 

* Functions have been made considerably easier to use (such as sampling initial 
values if they have not been provided by the user) 
* A proper vignette has been produced to show the user how to use the package.
* All functions now have examples in their documentation
* The `cancer` dataset has been wrangled into a format which can be more 
directly used with the BASSLINE functions.
* Parameters for the MCMC functions are now checked before the main algorithm 
runs.

## BASSLINE 0.0.0.9002 (2019-11-28)

### Patch notes

* Renamed package from SMLN to BASSLINE
* Began using Rcpp & RcppArmadillo for some basic functions

## SMLN 0.0.0.9001 (2019-11-19)

`R CMD check` runs with No errors, warnings, or notes. All functions are
currently implemented 100% in R.

## SMLN 0.0.0.9000 (2019-11-15)

Beginning of project to convert supplementary code from Vallejos and Steel to a
a fully functioning R library with robust documentation with Rcpp implementation
for improved efficiency. 
