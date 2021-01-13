---
title: 'DLNMs: Getting started'
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
draft: yes
---

This is the second post in a series about distributed lag non-linear models. Please read the [first post](/post/dlnm) for an introduction and a disclaimer.

# The dlnm package

The `dlnm` package offers two ways of fitting crossbasis functions: an "internal" and an "external" method. The "external" method involves fitting the crossbasis function outside of a model, using some functions in the `dlnm` package, then including the results as a predictor in a model such as a generalized linear model (GLM). I'm going to focus entirely on the "internal" method that fits the crossbasis function in the context of a generalied additive model (GAM) to take advantage of the penalization and other stuff the `mgcv` package offers.

# The data

Throughout this series, I'm going to use a subset of data from my postdoc project on *Heliconia acuminata*. In this subset, 100 plants were tracked over \_\_\_ years. Every year in February their size was recorded as height and number of shoots, and it was recorded whether or not they flowered. Any dead plants were marked as such. In this subset, all the plants were growing in plots in continuous, intact rainforest in the Biological Dynamics of Forest Fragments Project (BDFFP) near Manaus, Brazil. The goal is to determine how drought impacted growth, survival, and flowering probability with a potentially delayed and/or non-linear relationship. To that end, I've calculated SPEI, a measure of drought, where more negative numbers represent more severe drought. SPEI is monthly while the demography data is yearly. For every observation of a plant, there is an entire history of SPEI for the past 24 months from that observation (24 is arbitrary, but what I've chosen as the maximum lag time to model).


```r
head(ha)
```

```
## # A tibble: 6 x 11
##   plot  ha_id_number  year    ht ht_next  shts shts_next  flwr  surv
##   <fct> <fct>        <dbl> <dbl>   <dbl> <dbl>     <dbl> <dbl> <dbl>
## 1 5750  1020          1998    71      25     1         2     0     1
## 2 5750  1020          1999    25      57     2         2     0     1
## 3 5750  1020          2000    57      31     2         2     0     1
## 4 5750  1020          2001    31      45     2         4     0     1
## 5 5750  1020          2002    45      43     4         3     0     1
## 6 5750  1020          2003    43      51     3         4     0     1
## # … with 50 more variables: spei_history[,1] <dbl>, [,2] <dbl>, [,3] <dbl>,
## #   [,4] <dbl>, [,5] <dbl>, [,6] <dbl>, [,7] <dbl>, [,8] <dbl>, [,9] <dbl>,
## #   [,10] <dbl>, [,11] <dbl>, [,12] <dbl>, [,13] <dbl>, [,14] <dbl>,
## #   [,15] <dbl>, [,16] <dbl>, [,17] <dbl>, [,18] <dbl>, [,19] <dbl>,
## #   [,20] <dbl>, [,21] <dbl>, [,22] <dbl>, [,23] <dbl>, [,24] <dbl>,
## #   [,25] <dbl>, L[,1] <int>, [,2] <int>, [,3] <int>, [,4] <int>, [,5] <int>,
## #   [,6] <int>, [,7] <int>, [,8] <int>, [,9] <int>, [,10] <int>, [,11] <int>,
## #   [,12] <int>, [,13] <int>, [,14] <int>, [,15] <int>, [,16] <int>,
## #   [,17] <int>, [,18] <int>, [,19] <int>, [,20] <int>, [,21] <int>,
## #   [,22] <int>, [,23] <int>, [,24] <int>, [,25] <int>
```

- `plot` (factor): A plot ID
- `ha_id_number` (factor): A unique plant ID
- `year` (numeric): year of census
- `ht` (numeric): height in cm
- `ht_next` (numeric): height in the next year in cm
- `shts` (numeric): number of shoots
- `shts_next` (numeric): number of shoots in the next year
- `flwr` (numeric): Did a plant flower? 1 = yes, 0 = no
- `surv` (numeric): Did a plant survive? 1 = yes, 0 = no
- `spei_history` (c("matrix", "array")): A matrix column of the drought history starting in the current month (`spei_history[,1]` = February) and going back 24 months (`spei_history[,25]` = February 2 years ago)
- `L` (c("matrix", "array")): A matrix column describing the lag structure of `spei_history`.  Literally  just `0:24` for every row.



# Fit a DLNM


```r
library(mgcv) #for gam()
library(dlnm) #for the "cb" basis
```



```r
survival <-
  gam(surv ~ 
        shts +
        s(spei_history, L,
          bs = "cb",
          k = c(3, 24),
          xt = list(bs = "cr")),
      family = binomial(link = "logit"),
      method = "REML",
      data = ha)
```

Above is a simple DLNM with survival modeled as a function of number of shoots and the crossbasis function of SPEI over the past 24 months. `shts` is a fixed effect (i.e. not a smooth, but to be fit as a straight line), and the crossbasis is defined in `s(spei_history, L, …)`. `spei_history` and `L` are the two dimensions of the crossbasis function, `bs = "cb"` tells `gam()` that this is a crossbasis function from the `dlnm` package (calls `dlnm::smooth.construct.cb.smooth.spec` behind the scenes). `xt = list(bs = "cr")` tells it to use a cubic regression spline as the basis for both dimensions of the crossbasis function (but you can also mix and match marginal basis functions by providing a length 2 vector here).

# Problem 1: visualizing the results

Unfortunately `plot.gam()` does not work with these crossbasis functions. 


```r
plot.gam(survival)
```

```
## Error in plot.gam(survival): No terms to plot - nothing for plot.gam() to do.
```

The `dlnm` package provides some functions for visualizing the results of a DLNM, though I don't like them much.

First you use `crosspred()` to get predicted values for the DLNM.


```r
pred_dat <- crosspred("spei_history", survival)
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
plot(pred_dat, ptype = "contour")
```

<img src="{{< blogdown/postref >}}index_files/figure-html/unnamed-chunk-7-1.png" width="672" />
So, with this **very cursory** model, we might say that high SPEI about 15 months ago results in the highest survival.  This is probably **NOT** correct, I'm just trying to give you an idea of how to interpret this plot.

# So... PROBLEMS:

These are some of the problems I'm going to need to deal with.

- The plots
  - The color scale on the heatmap is bad and misleading (same range for red and blue despite different number of breaks)
  - The color scale is relative risk (I think?), not probability.
  - Plot shows response *relative* to some arbitrary SPEI (0, the mean, by default).  I.e. this is not really a marginal effects plot showing the response to drought, all else being average.
  - Some areas may not be well supported by the data (because of missing combinations of lag and SPEI), but that's not easy to see.
- There are parameters to be tweaked in the `s()` call
  - number of knots, `k`
  - marginal basis functions, `xt = list(bs = )`
- There are parameters to be tweaked in the `gam()` call
  - `gamma` adjusts smoothness, but how to choose it?
  - `select = TRUE` adds a shrinkage penalty, but what does it mean?
- Interpretation
  - How to get a reliable unbiased p-value? (help files for `gam()` indicate that p-values from `anova()` are probably wrong!)
  - Or should I not use p-values and instead use AIC-based model selection?
  - Can I somehow say something about in *which months* drought has a statistically significant effect on growth/survival/flowering?
  
  
Some of these are already solved, and I'll try to be brief about it in the next couple of posts, but a few of these I'm still very confused about.
