
[![Travis build status](https://travis-ci.org/USCbiostats/partition.svg?branch=master)](https://travis-ci.org/USCbiostats/partition) [![Build status](https://ci.appveyor.com/api/projects/status/9rnot7opu66xpyns/branch/master?svg=true)](https://ci.appveyor.com/project/gvegayon/partition-nqr95/branch/master) [![Coverage status](https://codecov.io/gh/USCbiostats/partition/branch/master/graph/badge.svg)](https://codecov.io/github/USCbiostats/partition?branch=master)

<!-- README.md is generated from README.Rmd. Please edit that file -->
partition
=========

The goal of partition is to ...

Installation
------------

You can install partition from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("USCbiostats/partition")
```

Example
-------

This is a basic example which shows you how to solve a common problem:

``` r
library(partition)

blk.vec = 2:20
c.lb = .2
c.ub = .4
n = 200

dat = sim_blk_diag_mvn( blk.vec, c.lb, c.ub, n  )

rslts = plot_dr( dat, method="PC1")
#> observed 2018-02-20 08:49:33 pct.var = 0.10
#> Percent variance explained: 0.10
#> observed 2018-02-20 08:49:34 pct.var = 0.14
#> Percent variance explained: 0.14
#> observed 2018-02-20 08:49:36 pct.var = 0.19
#> Percent variance explained: 0.19
#> observed 2018-02-20 08:49:38 pct.var = 0.23
#> Percent variance explained: 0.23
#> observed 2018-02-20 08:49:40 pct.var = 0.28
#> Percent variance explained: 0.28
#> observed 2018-02-20 08:49:41 pct.var = 0.32
#> Percent variance explained: 0.32
#> observed 2018-02-20 08:49:43 pct.var = 0.37
#> Percent variance explained: 0.37
#> observed 2018-02-20 08:49:45 pct.var = 0.41
#> Percent variance explained: 0.41
#> observed 2018-02-20 08:49:47 pct.var = 0.46
#> Percent variance explained: 0.46
#> observed 2018-02-20 08:49:50 pct.var = 0.50
#> Percent variance explained: 0.50
#> permuted 2018-02-20 08:49:53 pct.var=0.10
#> Percent variance explained: 0.10
#> permuted 2018-02-20 08:49:55 pct.var=0.14
#> Percent variance explained: 0.14
#> permuted 2018-02-20 08:49:57 pct.var=0.19
#> Percent variance explained: 0.19
#> permuted 2018-02-20 08:49:59 pct.var=0.23
#> Percent variance explained: 0.23
#> permuted 2018-02-20 08:50:00 pct.var=0.28
#> Percent variance explained: 0.28
#> permuted 2018-02-20 08:50:02 pct.var=0.32
#> Percent variance explained: 0.32
#> permuted 2018-02-20 08:50:04 pct.var=0.37
#> Percent variance explained: 0.37
#> permuted 2018-02-20 08:50:06 pct.var=0.41
#> Percent variance explained: 0.41
#> permuted 2018-02-20 08:50:08 pct.var=0.46
#> Percent variance explained: 0.46
#> permuted 2018-02-20 08:50:11 pct.var=0.50
#> Percent variance explained: 0.50
```

<img src="man/figures/README-example-1.png" width="100%" />
