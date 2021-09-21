# stagedtrees

<!-- [![stagedtrees](https://www.r-pkg.org/badges/version/stagedtrees)](https://cran.r-project.org/package=stagedtrees) -->
<!-- [![Build Status](https://travis-ci.com/gherardovarando/stagedtrees.svg?branch=main)](https://travis-ci.com/gherardovarando/stagedtrees) -->
<!--   [![R-CMD-check](https://github.com/gherardovarando/stagedtrees/workflows/R-CMD-check/badge.svg)](https://github.com/gherardovarando/stagedtrees/actions) -->
<!--  [![Coverage status](https://codecov.io/gh/gherardovarando/stagedtrees/branch/master/graph/badge.svg)](https://codecov.io/github/gherardovarando/stagedtrees?branch=main) -->
<!--  [![](https://cranlogs.r-pkg.org/badges/stagedtrees)](https://cran.r-project.org/package=stagedtrees) -->


### Overview 

This repository contains a modification of the original `stagedtrees` package 
that implements staged event trees, a probability model for categorical random variables. 
The original code can be found [here](https://github.com/gherardovarando/stagedtrees).

### Modifications
This package extends the original `stagedtrees` package with the following functions:

* `stages_ordered_bhc()` - performs Backward Hill-Climbing search over ordinal 
variables, restricting the merges to only "adjacent" categories.

* `exhaustive_ordered_search()` -  performs a full search over an ordinal 
variable and returns the model with best partition of the variable into n stages.

* `join_multiple_stages()` - joins multiple stages corresponding to one variable

* `stages_bhc_plot()`- simultaniously runs the backwatd Hill-Climbing algorithm 
and plots the stagedtree together with the probability barplot for every single merge of two stages.
