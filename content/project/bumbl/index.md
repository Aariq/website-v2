---
title: bumbl
date: '2017-12-13'
summary: An R package for modeling bumblebee colony growth
slug: bumbl
image:
  caption: Example output of bumbl package
  focal_point: Smart
tags:
- R
- modeling
- demography
math: false
---

The [bumbl package](https://github.com/Aariq/bumbl) implements a model for bumblebee colony growth described in Crone and Williams 2019. It models colony growth as having a change point at some time, tau, where the colony switches from growth and worker production to gyne production. The `bumbl()` function applies this model to data from multiple colonies, allowing for each colony to have it’s own tau and returns the original data augmented with coefficients from the changepoint model.

In addition, I'm working on implementing an integral projection model for bumblebee colony growth with variation in worker size distribution and resource return based on the (currently unpublished) work of [Natlie Kerr](https://nataliezoekerr.com/) and others.
