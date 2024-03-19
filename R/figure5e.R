#' Figure 5e
#
#' This function allows you generate figure4e - ZMIZ1 and KI67
#'
#' @keywords VIPER ER ZMIZ1
#' @examples figure4e()
#' @import ggpubr
#' @import reshape
#' @export


figure5e_ki67 <- function() {

    proteinData

    sum(is.na(proteinData[,'P03372'])) #0 missing ESR
    sum(is.na(proteinData[,'Q9ULJ6'])) #55 missing values in ZMIZ1
    sum(is.na(proteinData[,'P46013'])) #0 missing MKI67


    plot(log(proteinData[,'P03372']),log(proteinData[,'P46013']), xlab="ERa", ylab="KI67")
    plot(log(proteinData[,'P03372']),log(proteinData[,'Q9ULJ6']), xlab="ERa", ylab="ZMIZ1")
    plot(log(proteinData[,'Q9ULJ6']),log(proteinData[,'P46013']), xlab="ZMIZ1", ylab="KI67")


    library(ggpubr)
    library(reshape)
    df<-data.frame(id=(1:nrow(proteinData)),log(proteinData))
    #dfLong<-reshape::melt(df, id="id")

    dfLong<-data.frame(
        rbind(
            cbind("ERa-KI67", df$P03372, df$P46013),
            cbind("ERa-ZMIZ1", df$P03372, df$Q9ULJ6),
            cbind("KI67-ZMIZ1", df$P46013, df$Q9ULJ6)
        )
    )
    dfLong[,2:3]<-sapply(dfLong[,2:3], as.numeric)
    dfLong[,1]<-sapply(dfLong[,1], as.factor)
    sapply(dfLong[2:3,], class) #check numeric/factor

    colnames(dfLong)<-c('Comparison','Protein1','Protein2')


    #ggscatter(df,x='P03372',y='P46013', conf.int=TRUE) +
    #stat_cor(label.x = 3)

    p<-ggscatter(dfLong, x='Protein1',
                 y='Protein2',
                 add="reg.line",
                 conf.int=TRUE,
                 color='Comparison',
                 palette="jco",
                 shape='Comparison'
    ) +
        stat_cor(aes(color = Comparison), label.x = 1.5)  +
        xlab('Log(intensity Protein 1)') +
        ylab('Log(intensity Protein 2)') +
        ggtitle('ER, ZMIZ1 and KI67 Protein Abundance \n Correlates in TCGA (PXD024322)')



    #pdf(file="KI76vZMIZ1_Protein.pdf", width=5, height=6)
    p


}
