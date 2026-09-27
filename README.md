# xtrec

**Panel Unit Root Test Based on Recursive Detrending**

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/xtrec)](https://CRAN.R-project.org/package=xtrec)
<!-- badges: end -->

## Overview

`xtrec` implements the recursively detrended panel unit root tests proposed by
Westerlund (2015). Two variants are provided:

- **t-REC**: Basic test assuming iid errors
- **t-RREC**: Robust test accounting for serial correlation, cross-sectional
  dependence, and heteroskedasticity

Both tests have a standard normal limiting distribution, requiring no mean or
variance correction factors.

## Installation

```r
# CRAN
install.packages("xtrec")

# Development version
# remotes::install_github("muhammedalkhalaf/xtrec")
```

## Usage

```r
library(xtrec)

# Load example data
dat <- grunfeld_data()

# t-REC test (constant only)
res <- xtrec(dat, var = "invest", panel_id = "firm",
             time_id = "year", trend = 0L, robust = FALSE)
print(res)

# t-RREC test (with linear trend, robust)
res2 <- xtrec(dat, var = "invest", panel_id = "firm",
              time_id = "year", trend = 1L, robust = TRUE)
summary(res2)
```

## Reference

Westerlund, J. (2015). The effect of recursive detrending on panel unit root
tests. *Journal of Econometrics*, 185(2), 453–467.
<https://doi.org/10.1016/j.jeconom.2014.06.015>

## Author

Muhammad Alkhalaf <muhammedalkhalaf@gmail.com>
