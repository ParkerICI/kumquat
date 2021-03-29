This package is a stripped-down version of the [Citrus](https://github.com/nolanlab/citrus) package, which removes all the clustering and GUI code, and focuses exclusively on model building

It can be used to build models that associate features from single-cell data (e.g. the abundance of a specific cluster or gated cell population) to a clinical endpoint of interest (e.g. response to therapy)

The original reference for Citrus is

```
Robert V Bruggner, Bernd Bodenmiller, David L Dill, Robert J Tibshirani, Garry P Nolan
Automated identification of stratifying signatures in cellular subpopulations
Proc Natl Acad Sci U S A. 2014 Jul 1;111(26):E2770-7. doi: 10.1073/pnas.1408792111. Epub 2014 Jun 16.
```



# Installation

To install `kumquat` first intall the Bioconductor `impute` package

```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("impute")
```

and then install `kumquat` using the `devtools` package as follows

```R
devtools::install_github("ParkerICI/kumquat")
```





