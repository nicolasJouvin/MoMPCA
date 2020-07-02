
<!-- badges: start -->

[![Travis build
status](https://travis-ci.org/nicolasJouvin/MMPCA.svg?branch=master)](https://travis-ci.org/nicolasJouvin/MMPCA)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/MoMPCA)](https://cran.r-project.org/package=MoMPCA)
<!-- badges: end -->

> Cluster any count data matrix with a fixed number of variables.
> Implements the branch & bound Classification-Variational
> Expectation-Maximisation of this
> [paper](https://arxiv.org/abs/1909.00721) (to appear in Computational
> Statistics).

## Installation

**MoMPCA** is available on
[CRAN](https://cran.r-project.org/web/packages/MoMPCA/index.html) and
the development version available on
[Github](https://github.com/nicolasJouvin/MoMPCA).

### R Package installation

#### CRAN dependencies

**MoMPCA** needs the following CRAN R packages, so check that they are
are installed on your computer.

``` r
required_CRAN <- c("methods",
    "topicmodels",
    "tm",
    "Matrix",
    "slam",
    "magrittr",
    "dplyr",
    "stats",
    "doParallel",
    "foreach")
not_installed_CRAN <- setdiff(required_CRAN, rownames(installed.packages()))
if (length(not_installed_CRAN) > 0) install.packages(not_installed_CRAN)
```

#### Installing MoMPCA

  - For the last stable version, use the CRAN version

<!-- end list -->

``` r
install.packages("MoMPCA")
```

  - For the development version, use the github install

<!-- end list -->

``` r
remotes::install_github("nicolasJouvin/MoMPCA")
```

<!-- - For a specific tagged release, use -->

<!-- ```{r package tag, eval = FALSE} -->

<!-- remotes::install_github("nicolasJouvin/MoMPCA@tag_number") -->

<!-- ``` -->

## Usage and main fitting functions

The package comes with the BBCmsg data set and a `simulate_BBC()`
function wich allows to reproduce the simulation of the
[paper](https://arxiv.org/abs/1909.00721).

``` r
library(MoMPCA)
simu <- simulate_BBC(N = 400, L = 200, epsilon = 0, lambda = 1)
dtm <- simu$dtm.full
Ytruth <- simu$Ytruth # true clustering
```

The `dtm` is a `tm::DocumentTermMatrix()` object. The main fitting
function is `mmpca_clust()`, which allow for a parralel backend via its
argument `mc.cores`. There is a simple wrapper around this function
called `mmpca_clust_modelselect()` which allows for model selection of
\((Q, K)\) with ICL criterion. Please be aware that the greedy nature of
the algorithm may induce quite intensive computations.

### Greedy clustering with a mixture of multinomial PCA

``` r
res <- mmpca_clust(simu$dtm.full, Q = 6, K = 4,
                          Yinit = 'random',
                          method = 'BBCVEM',
                          max.epochs = 7,
                          keep = 1,
                          verbose = 2,
                          nruns = 2,
                          mc.cores = 1)
```

### Model selection

``` r
res <- mmpca_clust_modelselect(simu$dtm.full, Qs = 5:7, Ks = 3:5,
                          Yinit = 'kmeans_lda',
                          init.beta = 'lda',
                          method = 'BBCVEM',
                          max.epochs = 7,
                          nruns = 3,
                          verbose = 1)
               

best_model = res$models
```

## References

Please cite our work using the following reference:

  - N. Jouvin, P. Latouche, C. Bouveyron, A. Livartowski, G. Bataillon,
    Greedy clustering of count data through a mixture of multinomial PCA
    (To appear in [Computational
    Statistics](https://www.springer.com/journal/180))

and consider citing this package

``` r
citation('MoMPCA')
## Warning in citation("MoMPCA"): no date field in DESCRIPTION file of package
## 'MoMPCA'
## 
## To cite package 'MoMPCA' in publications use:
## 
##   Nicolas Jouvin (2020). MoMPCA: Inference and Clustering for
##   Mixture of Multinomial Principal Component Analysis. R package
##   version 1.0.0.
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {MoMPCA: Inference and Clustering for Mixture of Multinomial Principal
## Component Analysis},
##     author = {Nicolas Jouvin},
##     year = {2020},
##     note = {R package version 1.0.0},
##   }
## 
## ATTENTION: This citation information has been auto-generated from
## the package DESCRIPTION file and may need manual editing, see
## 'help("citation")'.
```
