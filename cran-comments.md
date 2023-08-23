---
title: "cran-comments"
author: "Simon Dedman"
date: "23 August 2023"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

***

## Test environments
* local R installation, R 4.2.2, xubuntu 23.04
* win-builder (devel and release)

***

## R CMD check results

0 errors | 2 warnings | 2 notes

Warnings:

* "Rd files with duplicated alias 'movegroup': ‘movegroup-package.Rd’ ‘movegroup.Rd’". Have found no way to remove this, caused by lifecycle and usethis autogeneration of movegroup-package.R and Rd.

* ‘qpdf’ is needed for checks on size reduction of PDFs. Not yet installed. No PDFs.

Notes:

* Possible code problems: no visible binding for global variables: named variables are column names in a csv exported by another function.

* Non-standard files/directories found at top level: README cran-comments

***

## Downstream dependencies

There are currently no downstream dependencies for this package
