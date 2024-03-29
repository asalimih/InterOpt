% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{data_GSE59520}
\alias{data_GSE59520}
\title{GSE59520 Embryonic tumors of testis qPCR array (Plasma)}
\format{
## `data_GSE59520`
A matrix with 597 rows and 36 columns:
}
\source{
<https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59520>
}
\usage{
data_GSE59520
}
\description{
A subset of GSE59520
miRs with ct values more than 39 in at least one sample and
miRs with ct values more than 35 in all samples are removed.
miRs with at most 3 missing values were imputed using nonedetects package
}
\keyword{datasets}
