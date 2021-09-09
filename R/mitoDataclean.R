#' @title mitoDataclean pipeline.
#' @description A machine learning approach to retrieve genuine mtDNA mutations from sample
#' contamination-derived chaos
#' @details Calls the functions to run each of the three steps of the pipeline
#' (extract features, classification, visualization) with the specified parameters.
#' See the individual function documentation for more details and required arguments.
#' Required steps: \code{\link{extract_features}}, \code{\link{classification}}.
#' \code{\link{visualization}}.
#' @param variants The variant information used as input for feature extraction;
#' this argument is required
#' @param visualization Should the plot of probability of contamination origin
#' for all variants in a sample be output?
#' this argument is optional
#' @return The sample-level and variant-level contamination information.
#' @export
#' @examples
#' data(variants)
#' result<-mitoDataclean(variants=variants)
#' variant.level<-result$variant.level
#' sample.level<-result$sample.level
mitoDataclean <- function(variants,visualization = FALSE){
  dat <- extract_features(mut.inf = variants)
  result <- classification(features = dat)
  visualization(variant.level.matrix = result$variant.level, verbose = visualization)
  return(result)
}
if (getRversion()>="2.15.1") {
  utils::globalVariables(c("dbsnp.info","hap.info",'mitomap.info','mitomap.var.freq','rf_model'))
}
