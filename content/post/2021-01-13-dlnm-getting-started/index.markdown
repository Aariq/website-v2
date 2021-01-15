---
title: 'DLNMs: building and visualizing crossbasis functions'
author: Eric R. Scott
date: '2021-01-13'
slug: dlnm-getting-started
categories: []
tags: [DLNMs, GAMs, R]
subtitle: ''
summary: ''
authors: []
lastmod: '2021-01-13T14:56:14-05:00'
featured: no
image:
  caption: ''
  focal_point: ''
  preview_only: no
projects: [heliconia]
draft: no
---

This is the second post in a series about distributed lag non-linear models. Please read the [first post](/post/dlnm) for an introduction and a disclaimer.

## The dlnm package

The `dlnm` package offers two ways of fitting crossbasis functions: an "internal" and an "external" method. The "external" method involves fitting the crossbasis function outside of a model, using some functions in the `dlnm` package, then including the results as a predictor in a model such as a generalized linear model (GLM). I'm going to focus entirely on the "internal" method that fits the crossbasis function in the context of a generalied additive model (GAM) to take advantage of the penalization and other stuff the `mgcv` package offers.

## The data

Throughout this series, I'm going to use a subset of data from my postdoc project on *Heliconia acuminata*. In this subset, 100 plants were tracked over a decade. Every year in February their size was recorded as height and number of shoots, and it was recorded whether or not they flowered. Any dead plants were marked as such. The goal is to determine how drought impacted growth, survival, and flowering probability with a potentially delayed and/or non-linear relationship. To that end, I've calculated SPEI, a measure of drought, where more negative numbers represent more severe drought. SPEI is monthly while the demography data is yearly. For every observation of a plant, there is an entire history of SPEI for the past 36 months from that observation.


```r
head(ha)
```

```
## # A tibble: 6 x 11
##   plot  ha_id_number  year  size size_next log_size log_size_next  flwr  surv
##   <fct> <fct>        <dbl> <dbl>     <dbl>    <dbl>         <dbl> <dbl> <dbl>
## 1 2107  107           1998    66       112     4.19          4.72     0     1
## 2 2107  107           1999   112       102     4.72          4.62     0     1
## 3 2107  107           2000   102        68     4.62          4.22     0     1
## 4 2107  107           2001    68        96     4.22          4.56     0     1
## 5 2107  107           2002    96       164     4.56          5.10     0     1
## 6 2107  107           2003   164       114     5.10          4.74     0     1
## # … with 74 more variables: spei_history[,1] <dbl>, [,2] <dbl>, [,3] <dbl>,
## #   [,4] <dbl>, [,5] <dbl>, [,6] <dbl>, [,7] <dbl>, [,8] <dbl>, [,9] <dbl>,
## #   [,10] <dbl>, [,11] <dbl>, [,12] <dbl>, [,13] <dbl>, [,14] <dbl>,
## #   [,15] <dbl>, [,16] <dbl>, [,17] <dbl>, [,18] <dbl>, [,19] <dbl>,
## #   [,20] <dbl>, [,21] <dbl>, [,22] <dbl>, [,23] <dbl>, [,24] <dbl>,
## #   [,25] <dbl>, [,26] <dbl>, [,27] <dbl>, [,28] <dbl>, [,29] <dbl>,
## #   [,30] <dbl>, [,31] <dbl>, [,32] <dbl>, [,33] <dbl>, [,34] <dbl>,
## #   [,35] <dbl>, [,36] <dbl>, [,37] <dbl>, L[,1] <int>, [,2] <int>, [,3] <int>,
## #   [,4] <int>, [,5] <int>, [,6] <int>, [,7] <int>, [,8] <int>, [,9] <int>,
## #   [,10] <int>, [,11] <int>, [,12] <int>, [,13] <int>, [,14] <int>,
## #   [,15] <int>, [,16] <int>, [,17] <int>, [,18] <int>, [,19] <int>,
## #   [,20] <int>, [,21] <int>, [,22] <int>, [,23] <int>, [,24] <int>,
## #   [,25] <int>, [,26] <int>, [,27] <int>, [,28] <int>, [,29] <int>,
## #   [,30] <int>, [,31] <int>, [,32] <int>, [,33] <int>, [,34] <int>,
## #   [,35] <int>, [,36] <int>, [,37] <int>
```

- `plot` (factor): A plot ID
- `ha_id_number` (factor): A unique plant ID
- `year` (numeric): year of census
- `size` (numeric): number of shoots x height in cm
- `size_next` (numeric): size in the next year
- `log_size` (numeric): log transformed size
- `log_size_next` (numeric): log transformed size next year
- `flwr` (numeric): Did a plant flower? 1 = yes, 0 = no
- `surv` (numeric): Did a plant survive? 1 = yes, 0 = no
- `spei_history` (c("matrix", "array")): A matrix column of the drought history starting in the current month (`spei_history[,1]` = February) and going back 24 months (`spei_history[,25]` = February 2 years ago)
- `L` (c("matrix", "array")): A matrix column describing the lag structure of `spei_history`.  Literally  just `0:24` for every row.



## Fit a DLNM


```r
library(mgcv) #for gam()
library(dlnm) #for the "cb" basis
```



```r
growth <-
  gam(log_size_next ~ 
        log_size +
        s(spei_history, L, # <- the two dimensions
          bs = "cb", # <- fit as crossbasis
          k = c(3, 24), # <- knots for each dimension
          xt = list(bs = "cr")), # <- what basis to use for each dimension
      family = gaussian(link = "identity"),
      method = "REML",
      data = ha)
```

Above is a simple DLNM with survival modeled as a function of number of shoots and the crossbasis function of SPEI over the past 36 months. `shts` is a fixed effect (i.e. not a smooth, but to be fit as a straight line), and the crossbasis is defined in `s(spei_history, L, …)`. `spei_history` and `L` are the two dimensions of the crossbasis function, `bs = "cb"` tells `gam()` that this is a crossbasis function from the `dlnm` package (calls `dlnm::smooth.construct.cb.smooth.spec` behind the scenes). `xt = list(bs = "cr")` tells it to use a cubic regression spline as the basis for both dimensions of the crossbasis function (but you can also mix and match marginal basis functions by providing a length 2 vector here).

## Problem 1: visualizing the results

Unfortunately `plot.gam()` does not work with these crossbasis functions. 


```r
plot.gam(growth)
```

```
## Error in plot.gam(growth): No terms to plot - nothing for plot.gam() to do.
```

The `dlnm` package provides some functions for visualizing the results of a DLNM, though I don't like them much.

First you use `crosspred()` to get predicted values for the DLNM.


```r
pred_dat <- crosspred("spei_history", growth)
```

```
## centering value unspecified. Automatically set to 0
```

Then you plot those with `plot.crosspred()`. The default is a 3D plot.


```r
plot(pred_dat)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-6-1.png" width="672" />

I prefer a heatmap, although the one produced here has some issues.


```r
plot(pred_dat, ptype = "contour", xlab = "SPEI", ylab = "lag(months)")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" width="672" />

First obvious problem is the colors.  The range is the same for red and blue, despite different number of breaks.  Second, the units are not what I'd expect.  For a marginal effects plot these should be the size of an average plant in year t+1, all else being equal.  This is plotting the size relative to the size at an average value of SPEI, which is a weird thing to think about.  That's because the package was built with epidemiology and relative risk in mind.  Here is the plot relative to SPEI = 1.5


```r
pred_dat <- crosspred("spei_history", growth, cen = 1.5)
plot(pred_dat, ptype = "contour", xlab = "SPEI", ylab = "lag(months)")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-8-1.png" width="672" />
## Solution (?)

So, I spent a lot of time writing a complicated function, `cb_margeff()`, to create data for a marginal effects plot. It creates a `newdata` data frame to be passed to `predict()` and loops across different matrixes with all columns of `spei_history` set to average except for one, representing a range of possible SPEI values.





```r
plotdata <- cb_margeff(growth, spei_history, L)
ggplot(plotdata, aes(x = x, y = lag, fill = fitted)) +
  geom_raster() +
  scale_fill_viridis_c("size in year t+1", option = "A") +
  scale_x_continuous("SPEI", expand = c(0,0)) +
  scale_y_continuous("lag (months)", expand = c(0,0))
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-10-1.png" width="672" />
Yeah, this is looking better.

The interpretation of this type of plot (which I would describe as a marginal effects plot, but correct me if I'm wrong) makes more sense to me.  If there was drought (low SPEI) about 8 months prior to the census, that's bad for growth.  Drought 20 months prior is good for growth though.

**BUT WAIT**

I poked around in `plot.gam` with `debug()` and it turns out the reason the plotting doesn't work is only because the author of `dlnm`, Gasparrini, didn't want it to work.

I can change a simple flag inside the `growth` model, and then it produces something very similar (identical?) to what I have above: 


```r
growth$smooth[[1]]$plot.me
```

```
## [1] FALSE
```

```r
growth$smooth[[1]]$plot.me <- TRUE
plot.gam(growth, scheme = 2)
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-11-1.png" width="672" />
Why is this default plot not available? It's literally **exactly** what I wanted, and I'm pretty sure there's nothing incorrect about it, but it worries me that the author of `dlnm` didn't want me to make it.


