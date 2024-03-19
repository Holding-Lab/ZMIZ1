#' Figure S15
#'
#' This function allows you generate figure S15 metabric - luminal vs basal ZMIZ1 expression
#'
#' @keywords RNAseq ZMIZ1
#' @examples figureS15_METABRIC()
#' @import ggpubr
#' @export

figureS15_METABRIC <- function() {
    ZMIZ1 <- any2entrez("ZMIZ1")
    ESR1 <- any2entrez("ESR1")
    tumAcros <- c("blca", "brca", "coad", "gbm", "hnsc", "kirc",
                  "kirp", "laml", "lgg", "lihc", "luad", "lusc", "ov",
                  "prad", "read", "sarc", "skcm", "stad", "thca", "ucec")
    niceAcros <- c("Bladder Carcinoma", "Breast Carcinoma", "Colon Adenocarcinoma",
                   "Glioblastoma", "Head and Neck Squamous Carcinoma", "Kidney Renal Clear Cell Carcinoma",
                   "Kidney Papillary Carcinoma", "Acute Myeloid Leukemia",
                   "Low Grade Glioma", "Liver Hepatocellular Carcinoma",
                   "Lung Adenocarcinoma", "Lung Squamous Carcinoma", "Ovarian Carcinoma",
                   "Prostate Carcinoma", "Rectal Adenocarcinoma", "Sarcoma",
                   "Skin Carcinoma", "Stomach Adenocarcinoma", "Thyroid Carcinoma",
                   "Utherine Corpus Endometroid Carcinoma")
    names(niceAcros) <- tumAcros
    plotfun <- function(tums) {
        for (tum in tums) {
            message(tum)
            breast <- FALSE
            if (tum %in% c("brca", "metabric")) {
                breast <- TRUE
            }
            expmat_brca <- rbind(expmat_brca_a, expmat_brca_b)
            expmat <- get(paste0("expmat_", tum))
            vipermat <- get(paste0("vipermat_", tum))
            subtypes <- get(paste0("subtypes_", tum))
            if (breast) {
                common <- intersect(colnames(vipermat), names(subtypes))
                vipermat <- vipermat[, common]
                subtypes <- subtypes[common]
            }
            expmat <- expmat[, colnames(vipermat)]
            if (breast) {
                colors <- setNames(rep("black", length(subtypes)),
                                   subtypes)
                colors[subtypes == "Basal"] <- "#FF0000AA"
                    colors[subtypes == "LumA"] <- "#00FFFFAA"
                        colors[subtypes == "LumB"] <- "#0000FFAA"
                            colors[subtypes == "Her2"] <- "#FFFF00AA"
                                colors[subtypes == "Tumor, Normal-Like"] <- "#00FF00AA"
            }
            else {
                colors <- "black"
            }

            if (tum == "brca") {
                title <- "TCGA"
                my_comparisons <- list( c("Basal", "Her2"), c("Basal", "LumA"), c("Basal", "LumB"),  c("Basal", "Normal") )
                labHeight=16.5
            }
            if (tum == "metabric") {
                title <- "METABRIC"
                my_comparisons <- list( c("Basal", "Her2"), c("Basal", "LumA"), c("Basal", "LumB"),  c("Basal", "Normal") )
                labHeight=14.5
            }

            vZ <- vipermat[ZMIZ1, ]
            vE <- vipermat[ESR1, ]
            eZ <- expmat[ZMIZ1, ]
            eE <- expmat[ESR1, ]


            longBRCA<-as.data.frame(cbind(eZ,subtypes))
            longBRCA$subtypes[longBRCA$subtypes == "Tumor, Normal-Like"] <-'Normal'
            longBRCA$eZ<-as.numeric(longBRCA$eZ)
            longBRCA$subtypes<-as.factor(longBRCA$subtypes)



            p <- ggboxplot(longBRCA, x = "subtypes", y = "eZ",
                           fill = "subtypes", palette =c("#FF0000AA","#FFFF00AA", "#00FFFFAA", "#0000FFAA","#00FF00AA"),
                           add = "jitter")
            q<- p + stat_compare_means(label.y=labHeight)+ # Add pairwise comparisons p-value
                stat_compare_means(method="wilcox.test",comparisons = my_comparisons) +ggtitle(title) +xlab("Subtype")+ylab("ZMIZ1 Expression")+ theme(legend.position = "none")


            return(q)

        }
    }





    table(get(paste0("subtypes_", 'brca')))
    # Basal               Her2               LumA               LumB             Normal
    # 173                 60                572                223                 18
    # Normal, Tumor-Like Tumor, Normal-Like
    # 92                 16

    table(get(paste0("subtypes_", 'metabric')))
    #
    # Basal   Her2   LumA   LumB Normal
    # 168     11    590    203     24


    return(plotfun('metabric'))
}

#' Figure S15 TCGA
#'
#' This function allows you generate figure S15 TCGA - luminal vs basal ZMIZ1 expression
#'
#' @keywords RNAseq ZMIZ1
#' @examples figureS15_TCGA()
#' @import ggpubr
#' @export

figureS15_TCGA <- function() {
    ZMIZ1 <- any2entrez("ZMIZ1")
    ESR1 <- any2entrez("ESR1")
    tumAcros <- c("blca", "brca", "coad", "gbm", "hnsc", "kirc",
                  "kirp", "laml", "lgg", "lihc", "luad", "lusc", "ov",
                  "prad", "read", "sarc", "skcm", "stad", "thca", "ucec")
    niceAcros <- c("Bladder Carcinoma", "Breast Carcinoma", "Colon Adenocarcinoma",
                   "Glioblastoma", "Head and Neck Squamous Carcinoma", "Kidney Renal Clear Cell Carcinoma",
                   "Kidney Papillary Carcinoma", "Acute Myeloid Leukemia",
                   "Low Grade Glioma", "Liver Hepatocellular Carcinoma",
                   "Lung Adenocarcinoma", "Lung Squamous Carcinoma", "Ovarian Carcinoma",
                   "Prostate Carcinoma", "Rectal Adenocarcinoma", "Sarcoma",
                   "Skin Carcinoma", "Stomach Adenocarcinoma", "Thyroid Carcinoma",
                   "Utherine Corpus Endometroid Carcinoma")
    names(niceAcros) <- tumAcros
    plotfun <- function(tums) {
        for (tum in tums) {
            message(tum)
            breast <- FALSE
            if (tum %in% c("brca", "metabric")) {
                breast <- TRUE
            }
            expmat_brca <- rbind(expmat_brca_a, expmat_brca_b)
            expmat <- get(paste0("expmat_", tum))
            vipermat <- get(paste0("vipermat_", tum))
            subtypes <- get(paste0("subtypes_", tum))
            if (breast) {
                common <- intersect(colnames(vipermat), names(subtypes))
                vipermat <- vipermat[, common]
                subtypes <- subtypes[common]
            }
            expmat <- expmat[, colnames(vipermat)]
            if (breast) {
                colors <- setNames(rep("black", length(subtypes)),
                                   subtypes)
                colors[subtypes == "Basal"] <- "#FF0000AA"
                    colors[subtypes == "LumA"] <- "#00FFFFAA"
                        colors[subtypes == "LumB"] <- "#0000FFAA"
                            colors[subtypes == "Her2"] <- "#FFFF00AA"
                                colors[subtypes == "Tumor, Normal-Like"] <- "#00FF00AA"
            }
            else {
                colors <- "black"
            }

            if (tum == "brca") {
                title <- "TCGA"
                my_comparisons <- list( c("Basal", "Her2"), c("Basal", "LumA"), c("Basal", "LumB"),  c("Basal", "Normal") )
                labHeight=16.5
            }
            if (tum == "metabric") {
                title <- "METABRIC"
                my_comparisons <- list( c("Basal", "Her2"), c("Basal", "LumA"), c("Basal", "LumB"),  c("Basal", "Normal") )
                labHeight=14.5
            }

            vZ <- vipermat[ZMIZ1, ]
            vE <- vipermat[ESR1, ]
            eZ <- expmat[ZMIZ1, ]
            eE <- expmat[ESR1, ]


            longBRCA<-as.data.frame(cbind(eZ,subtypes))
            longBRCA$subtypes[longBRCA$subtypes == "Tumor, Normal-Like"] <-'Normal'
            longBRCA$eZ<-as.numeric(longBRCA$eZ)
            longBRCA$subtypes<-as.factor(longBRCA$subtypes)



            p <- ggboxplot(longBRCA, x = "subtypes", y = "eZ",
                           fill = "subtypes", palette =c("#FF0000AA","#FFFF00AA", "#00FFFFAA", "#0000FFAA","#00FF00AA"),
                           add = "jitter")
            q<- p + stat_compare_means(label.y=labHeight)+ # Add pairwise comparisons p-value
                stat_compare_means(method="wilcox.test",comparisons = my_comparisons) +ggtitle(title) +xlab("Subtype")+ylab("ZMIZ1 Expression")+ theme(legend.position = "none")


            return(q)

        }
    }





    table(get(paste0("subtypes_", 'brca')))
    # Basal               Her2               LumA               LumB             Normal
    # 173                 60                572                223                 18
    # Normal, Tumor-Like Tumor, Normal-Like
    # 92                 16

    table(get(paste0("subtypes_", 'metabric')))
    #
    # Basal   Her2   LumA   LumB Normal
    # 168     11    590    203     24


    return(plotfun('brca'))
}



