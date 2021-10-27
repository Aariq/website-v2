---
title: Assessing the reliability of an R package
author: Eric R. Scott
date: '2021-10-27'
slug: []
categories:
  - Blog
tags:
  - R
subtitle: ''
summary: ''
authors: []
lastmod: '2021-10-27T16:59:54-04:00'
featured: no
image:
  caption: ''
  focal_point: ''
  preview_only: no
projects: []
---

In the most recent [rOpenSci community call](https://ropensci.org/commcalls/oct2021-statsreview02/), [Juliane Manitz](http://www.manitz.org/) presented on standards for statistical software in the pharmaceutical industry. One of her slides in particular was about how the pharmaceutical industry assesses the reliability of a package.  This inspired me to write a little bit about how I assess the reliability of an R package that's new to me, especially if I have the choice of two or more packages that do similar things.

I think the most common ways people in my field (ecology) choose an R package or function to learn to solve a particular problem are:

- It's what they were exposed to in a course or training
- It's widely used in their field (e.g. cited often)
- It's written or maintained by a "big-shot" in the field

After gaining experience in R package development, I've begun to place less importance on these factors and more importance on other factors that I'll discuss briefly below.

## Active Development

It's really important to me that a package I'm using is being actively developed. There are a number of reasons to choose a package with active development over one that is not updated often. One obvious reason is that if you encounter bugs, you can be more confident that they'll get fixed if you report them. Even mature, well established packages need active maintenance to ensure they remain functional as their dependencies get updated. I also like choosing packages with active development because I may have an opportunity help improve the package through my feedback and suggestions. 

So, how to assess if a package is being actively developed?

1) **Check CRAN release date.**  How recent was the latest stable version published to CRAN?
2) **Check bug reports.**  How many issues are open and how many have been closed?  How many really old bug reports are there, if any?  Has the package author at least responded to bug reports made over a month ago?
3) **Check GitHub development activity.**  Have there been somewhat recent commits or pull requests?  Is there a NEWS file documenting changes in the development version?


## Development Lifecycle

Where an R package or function is in its [lifecycle](https://lifecycle.r-lib.org/index.html) can help me make a decision about whether to use it over an alternative.  You might have noticed lifecycle badges in tidyverse package help files (for example, the "superseded" badge on the help file for `gather()` from the `tidyr` package) or on README files on GitHub (e.g. the "maturing" badge on the [`ipmr` package](https://github.com/levisc8/ipmr#readme) I've been learning lately).

<figure>
<img src="https://lifecycle.r-lib.org/reference/figures/lifecycle-superseded.svg" alt="a badge that reads 'lifecycle: superseded'">
<img src="https://lifecycle.r-lib.org/articles/figures/lifecycle-maturing.svg" alt = "a badge that read 'lifecycle: maturing'">
<img src="https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg" alt = "a badge that read 'lifecycle: maturing'">
<figcaption align = "center"><b>Examples of lifecycle badges </b></figcaption>
</figure>


I generally want to avoid learning anything that has been superseded (no longer being developed) and I definitely don't want to learn anything that's been deprecated. Instead I focus my efforts on learning functions or packages that are in the maturing or stable stages. If a package is still in the experimental stage, I may choose to learn it, but probably only if I'm interested in contributing to the package development and it's really the best option for what I'm trying to do.

These specific lifecycle stages are not always used in all packages, and sometimes the stage must be sussed out through some exploration.  For example, if you check the [website for the `raster` package](https://rspatial.org/), you'll see that it has been superseded by the `terra` package although the lifecycle badges are not used by these packages. 

## Testing

Unit tests are code that is written to check for the correct behavior of functions under a range of possible inputs from users. Tests can help catch bugs, check that users get informative error messages, and check for correctness. Most ecologists probably don't know to look for tests when choosing to use an R package. I think there is the assumption that if it's on CRAN, and has been cited before, then it's legit.  That simply is not true.  For example, the `spi` package was on CRAN (it's now been archived) and had been cited in papers, but didn't calculate anything *even close* to what it claimed to calculate. If the `spi` package had unit tests, perhaps this error would have been caught before it went to CRAN.
  The easiest way to check if a package has tests is to find the GitHub repository (often linked to from the CRAN page). If there is a "tests" folder, and it has some code in it, then that's a start.  You can also look for a badge in the README that describes the test coverage (how much of the package code is tested), such as the codecov badge [here](https://github.com/levisc8/ipmr#readme).  If test coverage is > 75% you're in really good shape. 
  
  <figure>
<img src="https://codecov.io/gh/levisc8/ipmr/branch/master/graph/badge.svg">
<figcaption align = "center"><b>Example of a Codecov badge </b></figcaption>
</figure>
  
  As [Noam Ross](https://www.noamross.net/) pointed out in the collaborative notes for the call, having a big user base *can* catch bugs without formal tests, but 
"tests are the 'ratchet', though---they make sure you don't go backwards, introducing old bugs again when you fix new ones".
