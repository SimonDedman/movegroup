
<!-- README.md is generated from README.Rmd. Please edit that file -->

# movegroup

<!-- badges: start -->

[![R-CMD-check](https://github.com/SimonDedman/movegroup/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/SimonDedman/movegroup/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->
<!-- badgeplacer(location = ".", status = "active", githubaccount = SimonDedman, githubrepo = gbm.auto, branch = master, name = "README.Rmd") -->

Especially on Linux systems it is recommended to type, in terminal:

``` r
sudo apt install libgeos-dev
sudo apt install libproj-dev
sudo apt install libgdal-dev
```

then manually install rgeos and rgdal in R/RStudio.

Also see each script’s Details section in the manual pages, as these
frequently contain tips or common bugfixes.

I strongly recommend that you download papers:

Kranstauber, B., Kays, R., LaPoint, S. D., Wikelski, M. and Safi, K.
(2012), A dynamic Brownian bridge movement model to estimate utilization
distributions for heterogeneous animal movement. Journal of Animal
Ecology. doi: 10.1111/j.1365-2656.2012.01955.x

Kranstauber, B., M. Smolla & A. K. Scharf. 2019. Move: visualizing and
analyzing animal track data. R package version 4.2.4 (at 2023-08-15).
<https://CRAN.R-project.org/package=move>.

Also it’s imperative you read the R help files for each function before
you use them. In RStudio: Packages tab, scroll to movegroup, click its
name, the click the function to see its man (manual) page. Read the
whole thing. Function man pages can also be accessed from the console by
typing

``` r
?function
```

------------------------------------------------------------------------

### movegroup

Visualizing and Quantifying Space Use Data for groups of animals

Automates dynamic Brownian bridge movement model calculation for
utilization distribution (UD) estimation for multiple individuals
simultaneously, using functions in the ‘move’ package. The authors are
indebted to the move package authors Bart Kraunstauber, Marco Smolla,
and Anne K Scharf, and to Sarah Becker for seed code which inspired the
development of this function.

------------------------------------------------------------------------

### scaleraster

Scales individual utilization distribution rasters and volume area
estimates

Scales individual-level utilization distribution (UD) rasters from 0 to
1 to facilitate interpretation as relative intensity of utilization (as
opposed to absolute), making comparisons across individuals and
interpretations at the group level more straightforward. Subsequently,
scaled individual-level rasters are aggregated to create a single
group-level UD raster. See www.GitHub.com/SimonDedman/movegroup for
issues, feedback, and development suggestions. There is an option to
account for bias in acoustic receiver array spatial representation (see
Details).

------------------------------------------------------------------------

### alignraster

Combines region-specific group-level UD rasters into a single raster.

Extends the spatial extent of each area-specific group-level raster to
the spatial extent shared by all rasters. This will only be required if
you have multiple individuals (e.g. different sharks) divided amongst a
few discrete areas (e.g. around different islands) and the effort
(e.g. receiver coverage) is different among islands. Not required for
multiple individuals all within the same region or sampling regime.

------------------------------------------------------------------------

### plotraster

Plots a group-level utilization distribution

Plots 50 and 95pct contours of a group-level utilization distribution
raster on a spatial map background. Contains functionality to also
visualize geographic locations of individual listening stations (e.g.,
acoustic receivers) as well as the entire surface UD.

------------------------------------------------------------------------

### moveLocErrorCalc

moveLocError calculator for ARGOS or state space models resulting in
95percent latlon confidence intervals

Builds a dataframe of original locations plus rowmeans of mean distance
of location extremities lon975, lat; lon025, lat; lon, lat975; lon,
lat025 from the centre point lon, lat.

------------------------------------------------------------------------

## Installation

You can install the released version of movegroup from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("movegroup")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
remotes::install_github("SimonDedman/movegroup")
```

------------------------------------------------------------------------

## Example

(See each function’s help file for specific examples, and the documents
listed above)

------------------------------------------------------------------------

## ToDo List

See GitHub issues section
<https://github.com/SimonDedman/movegroup/issues> Feel free to
contribute to this!

<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/master/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
