% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_features.R
\name{extract_features}
\alias{extract_features}
\title{Extract features used for identifying contamination-derived variants.}
\usage{
extract_features(mut.inf)
}
\arguments{
\item{mut.inf}{The matrix of mtDNA vatiant informatrion which is generated from mtDNA NGS
data containing information of sample name, variant position, reference allele, variant allele
and variant allele frequency.}
}
\value{
The function outputs a of matrix about the features used in random forest model.
}
\description{
This function extract features including variant-level features, sample-level
features and haplotype-matching features.
}
\examples{
extract_features(mut.inf = variants)
}
