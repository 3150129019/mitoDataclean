#' Detect sample cross-contamination, evaluate contamination level, and identify
#' contamination-derived variants
#' @description This function identifies the contamination-derived mtDNA variants by
#' random forest and detects sample contamination status and sample contamination level.
#' @param features The features extracted by \code{\link{extract_features}};
#' @return The function outputs a list: the first element is a matrix of sample-level
#' contamination results. The second list element is a matrix of variant-level contamination results.
#' @importFrom grDevices dev.off png
#' @importFrom graphics abline box legend par
#' @importFrom stats median predict sd
#' @importFrom randomForest randomForest
#' @export
#' @examples
#' features <- data.frame("Sample"=c("SampleA","SampleB"),
#'                        "Variant_allele_frequency"=c(0.9835,0.1036),
#'                        "Homoplasmic_status" = c(0,0),
#'                        "dbSNP_status" = c(1, 0),
#'                        "Number_of_contamination_derived_variant" = c(4, 0),
#'                        "Variant_frequency_in_MITOMAP" = c(0.439734004,0.213280188),
#'                        "Haplogroup_defining_variants" = c(1, 0),
#'                        "Heteroplasmic_variant_number" = c(6,3),
#'                        "CV_of_heteroplasmic_variant_allele_frequency" = c(0.373981559, 0.330185292),
#'                        "Sample_contamination_status" = c(1, 0),
#'                        "Predicted_contamination_status_for_major_heteroplasmic_variants" = c(0, 0),
#'                        "Predicted_contamination_status_for_minor_heteroplasmic_variants" = c(0, 1))
#' result <- classification(features = features)
classification <- function(features){
  vaf <- 1
  data<-as.matrix(features[which(as.numeric(features[,'Variant_allele_frequency'])>=vaf*0.01),])
  dat<-as.matrix(features[which(as.numeric(features[,'Variant_allele_frequency'])>=vaf*0.01),
                          -which(colnames(features)%in%'Sample')])
  mode(dat)<-"numeric"
  dat<-as.data.frame(dat)
  set.seed(1)
  Predict.contamination.status<-as.character(predict(rf_model,newdata=dat))
  Probability.of.contamination.origin<-predict(rf_model,dat,type = "prob")[,'1']
  Variant.matrix<-cbind(data,Predict.contamination.status,
                        Probability.of.contamination.origin)
  result<-list()
  result$variant.level<-Variant.matrix
  CL<-calculate_contamination_level(variant.level.matrix = Variant.matrix)
  result$sample.level<-CL
  return(result)
}
