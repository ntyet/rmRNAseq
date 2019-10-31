
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rmRNAseq

<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/ntyet/rmRNAseq.svg?branch=master)](https://travis-ci.org/ntyet/rmRNAseq)
<!-- badges: end -->

The goal of rmRNAseq is to conduct differential expression analysis
using RNA-seq data from repeated-measures designs, where mRNA samples
are obtained repeatedly at different times from same experimental units.
Our method is developed based on a general linear model framework with
continuous autoregressive correlation structure of order one,
accompanied by a parametric bootstrap inference strategy to conduct
general hypothesis testings.

## Installation

You can install the released version of rmRNAseq from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("rmRNAseq")
```

(The package has just been submitted to
[CRAN](https://CRAN.R-project.org) on June 26, 2019; it usually takes
about 10 days to receive their feedback.)

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ntyet/rmRNAseq")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(rmRNAseq)
data(dat)
data(design)
data(covset)
Subject <- covset$ear # identity of experimental units
Time <- covset$time # times at which mRNA samples are taken
Nboot <- 2  # for real data analysis, use Nboot at least 100
ncores <- 1 # for real data analysis and if the computer allows, increase ncores to save time
print.progress <- FALSE
saveboot <- FALSE
counts <- dat[1:3,]
C.matrix <- list()
# test for Line main effect
C.matrix[[1]] <- limma::makeContrasts(line2, levels = design)
# test for Time main effect
C.matrix[[2]] <- limma::makeContrasts(time2, time6, time24, levels = design)
names(C.matrix) <- c("line2", "time")
TCout <- rmRNAseq:::TC_CAR1(counts, design, Subject, Time, C.matrix,
Nboot, ncores, print.progress, saveboot)
#> running bootstrap sample nrep =  1 
#> ---------------------------------------------
#> running bootstrap sample nrep =  2 
#> ---------------------------------------------
names(TCout)
#> [1] "NewTime"    "TimeMinOut" "ori.res"    "pqvalue"
TCout$NewTime[1:4]
#> [1] 0.1492027 0.8028973 1.0000000 0.0000000
TCout$pqvalue$pv
#>       line2      time
#> 1 0.2857143 0.1428571
#> 2 0.1428571 0.1428571
#> 3 1.0000000 0.1428571
TCout$pqvalue$qv
#>       line2      time
#> 1 0.4285714 0.1428571
#> 2 0.4285714 0.1428571
#> 3 1.0000000 0.1428571
```
