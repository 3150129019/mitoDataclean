#' Plot the classification results against variant allele frequency and probability of contamination origin
#' @description Creates the scatter plots by variant allele frequency and 
#' probability of contamination origin for each sample.
#' @param variant.level.matrix matrix containing variant allele frequency and probability of contamination origin, 
#' as calculated by \code{\link{classification}}
#' @param verbose Should the plot be output? default is FALSE.
#' @return plots for each sample(returned silently)
#' @export
#' @examples
#' variant.level.matrix <- data.frame("Sample" = c("SampleA","SampleB"),
#'                                    "Variant_allele_frequency" = c(0.9850, 0.1036),
#'                                    "Predict.contamination.status" = c(0, 0),
#'                                    "Probability.of.contamination.origin" = c(0.1025, 0.3655))
#' visualization(variant.level.matrix, verbose=FALSE)
visualization <- function(variant.level.matrix, verbose){
  if(verbose){
    lwd<-2
    par(lwd = lwd,cex = 1.2)
    mat<-variant.level.matrix
    sample.names<-unique(mat[,1])
    dir.create('Plot')
    for(i in 1:length(sample.names))
    {
      dat<-mat[which(mat[,1]==sample.names[i]),]
      class.freq<-dat[,'Probability.of.contamination.origin']
      freq<-as.numeric(dat[,'Variant_allele_frequency'])
      cols<-rep(nrow(dat),0)
      cols[which(class.freq < 0.5)]<-'red'
      cols[which(class.freq >= 0.5)]<-'DarkTurquoise'
      png(filename =paste0('plot\\',sample.names[i],".png"),
          width = 400, height = 400)
      plot(freq,class.freq,xlim = c(0,1),
           ylim = c(0,1),col = cols,
           xlab = 'Variant allele frequency',
           ylab = 'Probability of contamination origin',
           cex.lab = 1.5,cex.axis = 1.2,cex = 1.5,
           lwd = 2,font.lab = 2,font.axis = 2)
      abline(v = 0.5,lwd = 1,lty = 2)
      box(lwd = 2)
      legend('topright',legend = c('Genuine','Contamination-derived'),
             bty = 'n',col = c('red','DarkTurquoise'),
             pch = 1,cex = 1.2,text.font = 2)
      dev.off()
    }
  }
}