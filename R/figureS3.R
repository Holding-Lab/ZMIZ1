#' Figure S3_METABRIC
#'
#' This function allows you generate figure s3_METABRIC.
#'
#' @keywords RNAseq ZMIZ1
#' @examples figureS3_METABRIC()
#' @import survival
#' @import survminer
#' @export

figureS3_METABRIC <- function() {


    filename<-system.file("extdata",
                                  "PATIENT_DATA_oncoprint.tsv",
                                  package = "ZMIZ1")

    survival<-data.frame(t(read.csv(filename,sep="\t")),stringsAsFactors = FALSE)
    survival<-survival[-1:-2,]

    survival[,5]<-as.numeric(survival[,5])
    colnames(survival)[5]<-'time'
    colnames(survival)[6]<-'status'
    colnames(survival)[4]<-'ER'
    colnames(survival)[12]<-'Expression'

    survivalFiltered<-survival[survival$ER=="Positive",]
    survivalFiltered<-survivalFiltered[!is.na(survivalFiltered$status==''),]


    survivalFiltered[survivalFiltered$status=='LIVING',]$status<-'0'
    survivalFiltered[survivalFiltered$status=='DECEASED',]$status<-'1'

    survivalFiltered$status<- as.numeric(survivalFiltered$status)


    survivalFiltered$Expression<-as.numeric(survivalFiltered$Expression)

    #<-median(survivalFiltered$Expression)-1>survivalFiltered$Expression


    res.cut<-surv_cutpoint(survivalFiltered, time = "time", event = "status", "Expression",
                           minprop = 0.1, progressbar = TRUE)
    summary(res.cut)

    plot(res.cut,"Expression",palette="npg")

    res.cat <- surv_categorize(res.cut)
    head(res.cat)


    fit <- survfit(Surv(time, status) ~Expression, data = res.cat)

    p<-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE)
    p
    summary(fit)

    surv_pvalue(fit )
    #variable         pval   method   pval.txt
    #1 Expression 6.819809e-07 Log-rank p < 0.0001
    return(p);
}

#' Figure S3_TCGA
#'
#' This function allows you generate figure s3_TCGA.
#'
#' @keywords RNAseq ZMIZ1
#' @examples figureS3_TCGA()
#' @import survival
#' @import survminer
#' @export

figureS3_TCGA <- function() {


    filename<-system.file("extdata",
                          "PATIENT_DATA_oncoprint_tcga.tsv",
                          package = "ZMIZ1")

    survival<-data.frame(t(read.csv(filename,sep="\t")),stringsAsFactors = FALSE)
    survival<-survival[-1:-2,]


    survival[,4]<-as.numeric(survival[,4])
    colnames(survival)[4]<-'time'
    colnames(survival)[5]<-'status'
    colnames(survival)[3]<-'ER'
    colnames(survival)[11]<-'Expression'

    survivalFiltered<-survival[survival$ER=="Positive",]
    survivalFiltered<-survivalFiltered[!is.na(survivalFiltered$status==''),]


    survivalFiltered[survivalFiltered$status=='LIVING',]$status<-'0'
    survivalFiltered[survivalFiltered$status=='DECEASED',]$status<-'1'

    survivalFiltered$status<- as.numeric(survivalFiltered$status)



    survivalFiltered$Expression<-as.numeric(survivalFiltered$Expression)

    #<-median(survivalFiltered$Expression)-1>survivalFiltered$Expression



    res.cut<-surv_cutpoint(survivalFiltered, time = "time", event = "status", "Expression",
                           minprop = 0.1, progressbar = TRUE)
    summary(res.cut)

    plot(res.cut,"Expression",palette="npg")

    res.cat <- surv_categorize(res.cut)
    head(res.cat)

    library("survival")
    fit <- survfit(Surv(time, status) ~Expression, data = res.cat)
    p<-ggsurvplot(fit, data = res.cat, risk.table = TRUE, conf.int = TRUE)
    p
    summary(fit)

    surv_pvalue(fit )
    #    variable        pval   method   pval.txt
    #1 Expression 0.001781372 Log-rank p = 0.0018
    return(p)
}

