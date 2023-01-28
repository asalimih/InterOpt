
<!-- README.md is generated from README.Rmd. Please edit that file -->

# InterOpt

<!-- badges: start -->
<!-- badges: end -->

The goal of InterOpt is to improve qPCR data normalization by providing
optimal weights for weighted mean of multiple internal controls. The
common method when dealing with multiple internal controls is taking a
geometric mean of them.

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("asalimih/InterOpt")
```

## Use Cases

1)  You have multiple internal controls (usually 2 or 3) and their
    corresponding raw CT values. Here is how you can aggregate them into
    one new internal control using weighted geometric mean:

``` r
library(InterOpt)

controls = c('RNU48','hsa-miR-16-5p')
x = data_GSE78870[controls,]
# x is a matrix of raw CT values, each row is an internal control and columns are samples
w = calcWeight(x, ctVal = TRUE, weight_method = 'geom_sd')
new_control = colSums(w*x)
# new_control is the weighted mean of the internal controls which can be used like a new internal control
```

2)  You have a dataset containing lots of genes. Here is how you can
    check all possible pair combinations (k=2) and calculate their
    corresponding aggregation weights (optimal weights to use in
    weighted mean) along with the stability measures (SD and CV) of each
    resulting aggregated combination.

``` r
library(InterOpt)

# only check miRNAs which have CT values less than 37 in all samples
mirs_for_combinations = rownames(data_GSE50013)[rowSums(data_GSE50013>37)==0]

result = run_experiment(data_source = data_GSE50013,
                      gr_source = groups_GSE50013,
                      ctVal_source = T,
                      sub_names = mirs_for_combinations,
                      norm_method = 'high_exp',
                      norm_method_exp_thr = 35,
                      k = 2,
                      weight_methods=c('geom', 'geom_sd'),
                      mc.cores=4,
                      verbose=F)
# Note: the data is automatically normalized using the norm_method and the aggregation weights are calculated based on the normalized data by default. moreover the stability measures are also calculated based on the normalized data.
head(result$res_source$geom_sd)
#>          Gene1           Gene2          w1        w2        CV       SD
#> 1 has-miR-1305     has-miR-155  0.13664283 0.8633572 0.4980741 1.223680
#> 2 has-miR-1305 hsa-miR-106b-5p  0.01469799 0.9853020 0.9600015 1.622915
#> 3 has-miR-1305  hsa-miR-126-3p  0.16942566 0.8305743 1.5454219 1.800675
#> 4 has-miR-1305   hsa-miR-1274A  0.07901806 0.9209819 1.2866663 2.375019
#> 5 has-miR-1305   hsa-miR-1274B  0.11235760 0.8876424 1.0162513 1.424279
#> 6 has-miR-1305    hsa-miR-1290 -0.04124614 1.0412461 1.3745626 2.228133
```

## Datasets

The following pre-processed datasets and their sample groups are
available in the package. for more information please refer to the
manual (e.g. ?data_GSE78870):  
- data_GSE78870, groups_GSE78870  
- data_GSE50013, groups_GSE50013  
- data_GSE67075, groups_GSE67075  
- data_GSE90828, groups_GSE90828  
- data_TCGA_BRCA, groups_TCGA_BRCA

## Note

The documentation for the package is in progress …
