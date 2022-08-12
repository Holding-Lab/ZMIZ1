#' Figure 3a
#
#' This function allows you generate figure1e
#'
#' @keywords VIPER ER ZMIZ1
#' @examples figure3a()
#' @importFrom stats setNames cor.test lm
#' @importFrom graphics plot grid legend abline mtext
#' @export


figure3a <- function() {

    ZMIZ1<-any2entrez("ZMIZ1")
    ESR1<-any2entrez("ESR1")

    tumAcros<-c("blca","brca","coad","gbm","hnsc","kirc","kirp","laml","lgg","lihc","luad","lusc","ov","prad","read","sarc","skcm","stad","thca","ucec")
    niceAcros<-c(
        "Bladder Carcinoma",
        "Breast Carcinoma",
        "Colon Adenocarcinoma",
        "Glioblastoma",
        "Head and Neck Squamous Carcinoma",
        "Kidney Renal Clear Cell Carcinoma",
        "Kidney Papillary Carcinoma",
        "Acute Myeloid Leukemia",
        "Low Grade Glioma",
        "Liver Hepatocellular Carcinoma",
        "Lung Adenocarcinoma",
        "Lung Squamous Carcinoma",
        "Ovarian Carcinoma",
        "Prostate Carcinoma",
        "Rectal Adenocarcinoma",
        "Sarcoma",
        "Skin Carcinoma",
        "Stomach Adenocarcinoma",
        "Thyroid Carcinoma",
        "Utherine Corpus Endometroid Carcinoma"
    )
    names(niceAcros)<-tumAcros


    ### Plot function
    plotfun<-function(tums){
        for(tum in tums){
            message(tum)
            breast<-FALSE
            if(tum%in%c("brca","metabric")){breast<-TRUE}

            expmat_brca<-rbind(expmat_brca_a,expmat_brca_b)

            expmat<-get(paste0("expmat_",tum))
            vipermat<-get(paste0("vipermat_",tum))
            subtypes<-get(paste0("subtypes_",tum))

            if(breast){
                common<-intersect(colnames(vipermat),names(subtypes))
                vipermat<-vipermat[,common]
                subtypes<-subtypes[common]
            }
            expmat<-expmat[,colnames(vipermat)]

            # COlor annotation
            if(breast){
                colors<-setNames(rep("black",length(subtypes)),subtypes)
                colors[subtypes=="Basal"]<-"#FF0000AA"
                colors[subtypes=="LumA"]<-"#00FFFFAA"
                colors[subtypes=="LumB"]<-"#0000FFAA"
                colors[subtypes=="Her2"]<-"#FFFF00AA"
                colors[subtypes=="Tumor, Normal-Like"]<-"#00FF00AA"
            } else {
                colors<-"black"
            }

            vZ<-vipermat[ZMIZ1,]
            vE<-vipermat[ESR1,]
            eZ<-expmat[ZMIZ1,]
            eE<-expmat[ESR1,]




            # Activity
            plot(vZ,vE,pch=20,col=colors,xlab="ZMIZ1 activity",ylab="ESR1 activity",cex.lab=2)
            grid()
            if(breast){
                legend("topleft",col=c("red","cyan","blue","yellow","green"),legend=c("Basal","LumA","LumB","Her2","Normal-Like"),pch=16,bg="white")
            }
            abline(lm(vE~vZ)$coef)
            scc<-cor.test(vE,vZ)
            mtext(paste0("SCC=",signif(scc$estimate,3)," (p=",signif(scc$p.value,3),")"),cex=0.8,font=2)


            title<-niceAcros[tum]
            if(tum=="brca"){title<-"TCGA"}
            if(tum=="metabric"){title<-"METABRIC"}

            title(title, outer=TRUE, line=-1.5,cex.main=2)


        }
    }



    ### Plot setup
    tums<-c("brca","metabric")
    plotfun(tums)



}
