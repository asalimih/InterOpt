#' GSE78870 Primary breast cancer qPCR array (Tissue)
#'
#' A subset of GSE78870
#' miRs with ct values more than 39 in at least one sample and
#' miRs with ct values more than 35 in all samples are removed.
#' miRs with at most 3 missing values were imputed using nonedetects package
#'
#' @format ## `data_GSE78870`
#' A matrix with 364 rows and 106 columns:
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78870>
"data_GSE78870"

#' GSE50013 Liver Cancer qPCR array (Plasma)
#'
#' A subset of GSE50013
#' miRs with ct values more than 39 in at least one sample and
#' miRs with ct values more than 35 in all samples are removed.
#' miRs with at most 3 missing values were imputed using nonedetects package
#'
#' @format ## `data_GSE50013`
#' A matrix with 241 rows and 40 columns:
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50013>
"data_GSE50013"

#' GSE67075 Colon Cancer qPCR array (Plasma)
#'
#' A subset of GSE67075
#' miRs with ct values more than 39 in at least one sample and
#' miRs with ct values more than 35 in all samples are removed.
#' miRs with at most 3 missing values were imputed using nonedetects package
#'
#' @format ## `data_GSE67075`
#' A matrix with 554 rows and 48 columns:
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE67075>
"data_GSE67075"

#' GSE57661 Early-stage breast cancer qPCR array (Plasma)
#'
#' A subset of GSE57661
#' miRs with ct values more than 39 in at least one sample and
#' miRs with ct values more than 35 in all samples are removed.
#' miRs with at most 3 missing values were imputed using nonedetects package
#'
#' @format ## `data_GSE57661`
#' A matrix with 176 rows and 48 columns:
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57661>
"data_GSE57661"

#' GSE59520 Embryonic tumors of testis qPCR array (Plasma)
#'
#' A subset of GSE59520
#' miRs with ct values more than 39 in at least one sample and
#' miRs with ct values more than 35 in all samples are removed.
#' miRs with at most 3 missing values were imputed using nonedetects package
#'
#' @format ## `data_GSE59520`
#' A matrix with 597 rows and 36 columns:
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59520>
"data_GSE59520"

#' GSE90828 Colon Cancer qPCR array (Plasma)
#'
#' A subset of GSE90828
#' miRs with ct values more than 39 in at least one sample and
#' miRs with ct values more than 35 in all samples are removed.
#' miRs with at most 3 missing values were imputed using nonedetects package
#'
#' @format ## `data_GSE90828`
#' A matrix with 336 rows and 53 columns:
#' @source <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE90828>
"data_GSE90828"

#' TCGA Breast Cancer miRNA-Seq data (Tissue)
#'
#' the data is normalized in counts per million (cpm)
#'
#' @format ## `data_TCGA_BRCA`
#' A matrix with 185 rows and 1201 columns:
"data_TCGA_BRCA"

#' sample groups of data_TCGA_BRCA data
#'
#' each row represents the tumor subtype of each sample
#'
#' @format ## `groups_TCGA_BRCA`
#' A matrix with 1201 rows and 7 columns:
"groups_TCGA_BRCA"

#' sample groups of data_GSE78870 data
#'
#' each element represents the group of the sample
#'
#' @format ## `groups_GSE78870`
#' A factor with 106 elements
"groups_GSE78870"

#' sample groups of data_GSE50013 data
#'
#' each element represents the group of the sample
#'
#' @format ## `groups_GSE50013`
#' A factor with 40 elements
"groups_GSE50013"


#' sample groups of data_GSE67075 data
#'
#' each element represents the group of the sample
#'
#' @format ## `groups_GSE67075`
#' A factor with 48 elements
"groups_GSE67075"

#' sample groups of data_GSE90828 data
#'
#' each element represents the group of the sample
#'
#' @format ## `groups_GSE90828`
#' A factor with 53 elements
"groups_GSE90828"
