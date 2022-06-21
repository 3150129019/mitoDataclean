# mitoDataclean
  mitoDataclean, a random forest-based machine learning package for accurate identification of cross-contamination, evaluation of contamination level, and detection of contamination-derived variants in mtDNA NGS data.<br>
The main function of mitoDataclean are:<br>

  * Identification of sample cross-contamination.<br>
  * Evaluation of sample contamination level.<br>
  * Detection of contamination-derived variants.<br>



# Installation 

The mitoDatalclean R package is installed from github using devtools. Installation requires the ability to compile R packages. This means that the R tool-chain must be installed, which requires the [Xcode command-line tools](http://railsapps.github.io/xcode-command-line-tools.html) on Mac and [Rtools](https://cran.r-project.org/bin/windows/Rtools/) on Windows. To install the most recent version of the mitoDataclean package, start R (version "4.0") and enter:

```
library(devtools)
install_github("3150129019/mitoDataclean")
```



# Usage

## Example

```
library(mitoDataclean)
data(variants)
result<-mitoDataclean(variants=variants)
variant.level<-result$variant.level
sample.level<-result$sample.level
```

### Input

  variants is a tab-formatted mtDNA variant table, which was generated from mtDNA NGS data containing information of sample name, variant position, 
reference allele, variant allele and variant allele frequency(VAF), was applied as input file into mitoDataclean package. 

| Name    | Position | Reference allele | Alteration allele | variant allele frequency |
| ------- | -------- | ---------------- | ----------------- | ------------------------ |
| SampleA | 10398    | A                | G                 | 1                        |
| SampleA | 10400    | C                | T                 | 1                        |
| SampleB | 12091    | T                | C                 | 0.1014                   |
| SampleB | 12308    | A                | G                 | 0.9167                   |
| SampleC | 11719    | G                | A                 | 1                        |
| SampleC | 13327    | A                | G                 | 0.0849                   |

### Output 

  mitoDataclean outputs variant-level table and sample-level table including informations about contamination.

#### 1. variant-level table 

- **Contamination-derived status：** If a variant is predicted as contamination-derived variants by mitoDataclean, this value is set to 1.
- **Probability of contamination origin**: According to majority voting, the probability of contamination origin was calculated by random forest model. When the probability was no less than 0.5, a variant was predicted as contamination-derived.
- **Features :** Features used in random forest model was described below.

| Features                                                  | Type of value | Number of      distinct values | Description                                                  |
| --------------------------------------------------------- | ------------- | ------------------------------ | ------------------------------------------------------------ |
| Homoplasmic status                                        | Boolean       | 2                              | The value of  homoplasmic variants status is set to 1.       |
| dbSNP status                                              | Boolean       | 2                              | The value  of dbSNP status is set to 1, if a variant was collection in dbSNP database. |
| Variant frequency in MitoMAP                              | Double        | Numeric                        | It is  calculated as the ratio of number of one variant over the total number of all  variants in 54234 full-length mtDNA in Mitomap database, while this value is assigned to 0 if a  variant is not present in the full-length mtDNA. |
| Haplogroup-defining variants                              | Boolean       | 2                              | Haplogroup-defining  variants are stable mtDNA polymorphisms that identify individuals to  mitochondrial haplogroups. if  a variant is a haplogroup-defining variant, this value is set to 1. |
| Heteroplasmic variant number                              | Integer       | Numeric                        | Number of  Heteroplasmic variants ( variant allele frequency < 99% ) in a sample. |
| CV of heteroplasmic variant allele frequency              | Double        | Numeric                        | Calculating the coefficient of  variation of heteroplasmic variant allele frequency in  a sample. |
| Predict contamination status for major                    | Boolean       | 2                              | It is a Boolean  feature based on whether the major heteroplasmic variants is presented in the  most-matched sequence in Mitomap or not. |
| Predict contamination status for minor                    | Boolean       | 2                              | It is a Boolean  feature based on whether the minor heteroplasmic variants is presented in the  most-matched sequence in Mitomap or not. |
| Sample contamination status                               | Boolean       | 2                              | It is a Boolean  feature based on the presence of predicted heteroplasmic variants. |
| Number of predicted contamination-derived major and minor | Integer       | Numeric                        | Number of  predicted contamination-derived major and minor variants in a sample. |

#### 2. sample-level table 

- **Contamination status：** Depending on number of predicted contamination-derived variants(≥2) in a sample, mitoDataclean assigns a contamination status to each sample. This value can either be **YES** or **NO**.
- **Contamination level:** The detectable contamination level for each sample. Contamination level was estimated as the median heteroplasmy level of all detected contamination-derived variants obtained from mitoDataclean.
- **Total variant number:** Total number of variants in a sample.
- **Homoplasmy number:** Total  number of homoplasmic variants in a sample.
- **Heteroplasmy number:** Total number of  heteroplasmic variants in a sample.
- **Contaminated-derived variants number(major):**  Total number of  contaminated-derived major Heteroplasmic variants number in a sample.
- **Contaminated-derived variants number(minor):**  Total number of  contaminated-derived minor Heteroplasmic variants number in a sample.

# Citation

Su L, Guo S, Guo W, Ji X, Liu Y, Zhang H, Huang Q, Zhou K, Guo X, Gu X, Xing J. mitoDataclean: A machine learning approach for the accurate identification of cross-contamination-derived tumor mitochondrial DNA mutations. Int J Cancer. 2022 May 15;150(10):1677-1689. doi: 10.1002/ijc.33927. Epub 2022 Jan 31. PMID: 35001369.

# Contact 

Bugs and difficulties in using mitoDataclean are welcome on [the issue tracker](https://github.com/3150129019/mitoDataclean/issues). 

Or you can address comments/questions/suggestions regarding this R package to: 929025191@qq.com<br>

