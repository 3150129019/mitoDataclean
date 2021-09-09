#' Extract features used for identifying contamination-derived variants.
#' @description This function extract features including variant-level features, sample-level
#' features and haplotype-matching features.
#' @param mut.inf The matrix of mtDNA vatiant informatrion which is generated from mtDNA NGS
#' data containing information of sample name, variant position, reference allele, variant allele
#' and variant allele frequency.
#' @return The function outputs a of matrix about the features used in random forest model.
#' @export
#' @examples
#' extract_features(mut.inf = variants)
extract_features<-function(mut.inf){
  paste.snp.inf<-paste(mitomap.var.freq[,'pos'],mitomap.var.freq[,'ref'],mitomap.var.freq[,'alt'])
  mitomap.sample<-unique(mitomap.info[,'sample'])
  mitomap.var.num<-table(mitomap.info[,'sample'])[mitomap.sample]
  mitomap.sample.num<-table(mitomap.info[,'sample'])
  dbsnp.info.p<-paste(dbsnp.info[,'pos'],dbsnp.info[,'ref'],dbsnp.info[,'alt'])

  res<-list()
  vaf<-1
  mut.inf<-as.matrix(mut.inf)
  colnames(mut.inf)<-c('name','pos','ref','alt','freq')
  sampleID<-unique(mut.inf[,1])

  for(i in 1:length(sampleID))
  {
    dat<-mut.inf[which(sampleID[i]==mut.inf[,1] & mut.inf[,'pos']!=''),,drop=F]
    mode(dat[,'pos'])<-'numeric'

    total<-length(which(as.numeric(dat[,'freq'])>=vaf*0.01))
    homo.num<-nrow(dat[which(as.numeric(dat[,'freq'])>=(1-vaf*0.01)),,drop=F])
    hetro.num<-nrow(dat[which(as.numeric(dat[,'freq'])<(1-vaf*0.01)
                              & as.numeric(dat[,'freq'])>=(vaf*0.01) ),,drop=F])
    hetro.1.freq<-as.numeric(dat[which(as.numeric(dat[,'freq'])<(1-vaf*0.01)
                                       & as.numeric(dat[,'freq'])>0.5 ),'freq',drop=F])
    hetro.2.freq<-as.numeric(dat[which(as.numeric(dat[,'freq'])<=0.5
                                       & as.numeric(dat[,'freq'])>=(vaf*0.01) ),'freq',drop=F])
    hetro.1.2.freq<-c(1-hetro.1.freq,hetro.2.freq)
    if(length(hetro.1.2.freq)<=1){major.minor.cv<-0}else{major.minor.cv<-sd(hetro.1.2.freq)/mean(hetro.1.2.freq)}

    pos<-as.numeric(dat[,'pos'])
    freq<-as.numeric(dat[,'freq'] )
    phylotree.state<-snp.status<-snp.freqs<-predict_contamination_status<-homo.status<-rep(0,length(pos))
    phylotree.state[which(!is.na(match(paste0(dat[,'ref'],dat[,'pos'],dat[,'alt']),hap.info[,1])))]<-1
    homo.status[which(as.numeric(dat[,'freq'])>=(1-vaf*0.01))]<-1

    homo.dat<-dat[which(as.numeric(dat[,'freq'])>=(1-vaf*0.01)),,drop=F]
    major.dat<-dat[which(as.numeric(dat[,'freq'])< (1-vaf*0.01) & as.numeric(dat[,'freq'])>0.5 ),,drop=F]
    minor.dat<-dat[which(as.numeric(dat[,'freq'])<0.5),,drop=F]
    homo.minor.dat<-rbind(homo.dat,minor.dat)
    homo.major.dat<-rbind(homo.dat,major.dat)

    homo.major.var<-paste0(homo.major.dat[,'ref'],homo.major.dat[,'pos'],homo.major.dat[,'alt'])
    homo.minor.var<-paste0(homo.minor.dat[,'ref'],homo.minor.dat[,'pos'],homo.minor.dat[,'alt'])
    minor.var<-paste0(minor.dat[,'ref'],minor.dat[,'pos'],minor.dat[,'alt'])
    major.var<-paste0(major.dat[,'ref'],major.dat[,'pos'],major.dat[,'alt'])

    major.inter.sample<-mitomap.info[mitomap.info[,'variant'] %in% major.var,'sample']
    major.inter.sample.num<-table(major.inter.sample)
    major.inter.sample.num<-major.inter.sample.num[mitomap.sample]
    major.inter.sample.num[which(is.na(major.inter.sample.num))]<-0
    names(major.inter.sample.num)<-mitomap.sample

    minor.inter.sample<-mitomap.info[mitomap.info[,'variant'] %in% minor.var,'sample']
    minor.inter.sample.num<-table(minor.inter.sample)
    minor.inter.sample.num<-minor.inter.sample.num[mitomap.sample]
    minor.inter.sample.num[which(is.na(minor.inter.sample.num))]<-0
    names(minor.inter.sample.num)<-mitomap.sample

    homo.major.inter.sample<-mitomap.info[mitomap.info[,'variant'] %in% homo.major.var,'sample']
    homo.major.inter.sample.num<-table(homo.major.inter.sample)[mitomap.sample]
    homo.major.inter.sample.num[which(is.na(homo.major.inter.sample.num))]<-0
    names(homo.major.inter.sample.num)<-mitomap.sample
    homo.major.score<-major.inter.sample.num/length(homo.major.var)+(homo.major.inter.sample.num*(length(homo.major.var)+mitomap.var.num))/(2*length(homo.major.var)*mitomap.var.num)

    homo.minor.inter.sample<-mitomap.info[mitomap.info[,'variant'] %in% homo.minor.var,'sample']
    homo.minor.inter.sample.num<-table(homo.minor.inter.sample)[mitomap.sample]
    homo.minor.inter.sample.num[which(is.na(homo.minor.inter.sample.num))]<-0
    names(homo.minor.inter.sample.num)<-mitomap.sample
    homo.minor.score<-minor.inter.sample.num/length(homo.minor.var)+(homo.minor.inter.sample.num*(length(homo.minor.var)+mitomap.var.num))/(2*length(homo.minor.var)*mitomap.var.num)

    major.sample<-names(which(homo.major.score==max(homo.major.score,na.rm=T)))
    predict_major_contamination_status<-rep(0,nrow(dat))
    pre.conted.major.var<-unique(gsub('A|T|C|G','',
                                      mitomap.info[which(mitomap.info[,'sample']%in%major.sample
                                                         & mitomap.info[,'variant']%in%major.var),'variant']))
    predict_major_contamination_status[match(as.numeric(pre.conted.major.var),as.numeric(dat[,'pos']))]<-1

    minor.sample<-names(which(homo.minor.score==max(homo.minor.score,na.rm=T)))
    predict_minor_contamination_status<-rep(0,nrow(dat))
    pre.conted.minor.var<-unique(gsub('A|T|C|G','',
                                      mitomap.info[which(mitomap.info[,'sample']%in%minor.sample
                                                         & mitomap.info[,'variant']%in%minor.var),'variant']))
    predict_minor_contamination_status[match(as.numeric(pre.conted.minor.var),as.numeric(dat[,'pos']))]<-1

    major_inter<-length(pre.conted.major.var)
    minor_inter<-length(pre.conted.minor.var)
    major.value<-homo.major.score[major.sample[1]]
    minor.value<-homo.minor.score[minor.sample[1]]
    major.minor.inter<-major_inter+minor_inter

    p<-paste(dat[,'pos'],dat[,'ref'],dat[,'alt'])
    snp.status[match(dbsnp.info.p,p,nomatch=0)]<-1
    snp.freqs<-mitomap.var.freq[match(p,paste.snp.inf),'freq']
    snp.freqs[which(is.na(snp.freqs))]<-0

    predict_contamination_status<-predict_minor_contamination_status+predict_major_contamination_status
    if(all(predict_contamination_status==0))
    {
      predict_sample_contamination<-0
    }else{
      predict_sample_contamination<-1
    }

    res[[i]]<-cbind(sampleID[i],pos, freq , total,
                    hetro.num , homo.num,homo.status,
                    phylotree.state,
                    predict_major_contamination_status,
                    predict_minor_contamination_status,
                    major_inter,
                    minor_inter,
                    major.minor.inter,
                    major.value,minor.value,
                    snp.status ,snp.freqs,
                    major.minor.cv,
                    predict_contamination_status)
  }
  re<-do.call(rbind,res)
  colnames(re)<-c('Sample','Position','Variant_allele_frequency','Total_variant_number',
                  'Heteroplasmic_variant_number','Homoplasmic_variant_number','Homoplasmic_status',
                  'Haplogroup_defining_variants',
                  'Predicted_contamination_status_for_major_heteroplasmic_variants',
                  'Predicted_contamination_status_for_minor_heteroplasmic_variants',
                  'Number_of_contamination_derived_variant_for_major',
                  'Number_of_contamination_derived_variant_for_minor',
                  'Number_of_contamination_derived_variant',
                  'Matching_score(major)','Matching_score(minor)',
                  'dbSNP_status' ,'Variant_frequency_in_MITOMAP' ,
                  'CV_of_heteroplasmic_variant_allele_frequency',
                  'Sample_contamination_status'
  )
  return(re)
}
