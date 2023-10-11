
<!-- README.md is generated from README.Rmd. Please edit that file -->

# InterOpt

<!-- badges: start -->
<!-- badges: end -->

<img src="https://github.com/asalimih/InterOpt/blob/8eaf86e459fc14f54594429d44097b8a4e252ef4/man/figures/TheGraphics_v2.jpg" width="512" height="512" />

The goal of InterOpt is to improve qPCR data normalization by providing
optimal weights for weighted mean of multiple internal controls
(reference genes). It can be easily utilized in combination with usual
ΔΔCT method. Here instead of taking average of multiple internal
controls (which is the common method) you use weighted mean.

## Online web application

A Shiny App is provided at [interopt.ir](https://interopt.ir/). Its step
by step tutorial is described [here](README_WebApp.md)

## Installation

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("asalimih/InterOpt")
```

## Use Cases

### 1) qPCR with multiple reference genes

You have multiple internal controls (usually 2 or 3) and their
corresponding raw CT values. Here is how you can aggregate them into one
new internal control using weighted geometric mean:

``` r
library(InterOpt)

controls = c('RNU48','hsa-miR-16-5p')
x = data_GSE78870[controls,]
# x is a matrix of raw CT values, each row is an internal control and columns are samples
w = calcWeight(x, ctVal = TRUE, weight_method = 'geom_sd_plus')
new_control = colSums(w*x)
# new_control is the weighted mean of the internal controls which can be used like a new internal control
```

### 2) Aggregated reference gene selection

You have a dataset containing lots of genes. Here is how you can check
all possible pair combinations (k=2) and calculate their corresponding
aggregation weights (optimal weights to use in weighted mean) along with
the stability measures (SD and CV) of each resulting aggregated
combination.

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
                      weight_methods=c('geom', 'geom_sd_hybrid'),
                      output_agg_refs = T,
                      mc.cores=4,
                      verbose=F)
# Note: the data is automatically normalized using the norm_method and the aggregation weights
#       are calculated based on the normalized data by default. Moreover, the stability measures
#       are also calculated based on the normalized data.

# The calculated weights:
head(result$res_source$geom_sd_hybrid)
#>          Gene1           Gene2          w1        w2        CV       SD
#> 1 has-miR-1305     has-miR-155  0.12348869 0.8765113 0.4845656 1.225702
#> 2 has-miR-1305 hsa-miR-106b-5p  0.03695359 0.9630464 0.9646868 1.626087
#> 3 has-miR-1305  hsa-miR-126-3p  0.15735982 0.8426402 1.5436056 1.801823
#> 4 has-miR-1305   hsa-miR-1274A  0.09881178 0.9011882 1.2826722 2.376691
#> 5 has-miR-1305   hsa-miR-1274B  0.09709854 0.9029015 1.0153004 1.426435
#> 6 has-miR-1305    hsa-miR-1290 -0.01028321 1.0102832 1.3578631 2.231676

# The corresponding weighted mean of the miRNAs:
head(result$aggregated_refs$geom_sd_hybrid[,1:7])
#>          Gene1           Gene2 SAMPLE.1 SAMPLE.2 SAMPLE.3 SAMPLE.4 SAMPLE.5
#> 1 has-miR-1305     has-miR-155 30.01693 29.46874 27.23921 31.16995 25.41960
#> 2 has-miR-1305 hsa-miR-106b-5p 30.88982 30.83069 27.57111 30.28982 27.07828
#> 3 has-miR-1305  hsa-miR-126-3p 24.66599 24.07716 22.03350 24.06599 24.58630
#> 4 has-miR-1305   hsa-miR-1274A 29.12176 31.21664 27.02058 31.13521 23.37273
#> 5 has-miR-1305   hsa-miR-1274B 27.15234 28.53192 24.86886 28.80959 22.65435
#> 6 has-miR-1305    hsa-miR-1290 26.76684 30.72340 28.49666 33.13780 26.65450
```

In order to have NormFinder and Genorm stability measures in the output,
you need to first add `InterOptCuda` to your PATH (Please check
[InterOptCuda](https://github.com/asalimih/InterOptCuda)). Then just add
`algors=c('SD','CV','Genorm','NormFinder')` argument in `run_experiment`

## Datasets

The following pre-processed datasets and their sample groups are
available in the package. for more information please refer to the
manual (e.g. ?data_GSE78870):  
- `data_GSE78870`, `groups_GSE78870`  
- `data_GSE50013`, `groups_GSE50013`  
- `data_GSE57661`, `groups_GSE57661`  
- `data_GSE59520`, `groups_GSE59520`  
- `data_GSE67075`, `groups_GSE67075`  
- `data_GSE90828`, `groups_GSE90828`  
- `data_TCGA_BRCA`, `groups_TCGA_BRCA`

## Note

The documentation for the package is in progress …

## Citation

If you used this method in your studies please cite the following
paper:  
[Salimi, Adel, Saeid Rahmani, and Ali Sharifi-Zarchi. “InterOpt:
Improved gene expression quantification in qPCR experiments using
weighted aggregation of reference genes.” iScience 26.10
(2023)](https://doi.org/10.1016/j.isci.2023.107945)
