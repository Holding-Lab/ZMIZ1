#' Figure S9
#'
#' This function allows you generate figure S9.
#'
#' @keywords RNAseq ZMIZ1
#' @examples figureS9_siZMIZ1_ESR1_ZMIZ1_reads()
#' @import DESeq2
#' @import vulcan
#' @import ggpubr
#' @export

figureS9_siZMIZ1_ESR1_ZMIZ1_reads <- function() {

    result_list <- list()
    for (cells in c("MCF7")) {
        for (time in unique(annotationTable$Treatment_Duration)) {
            resname <- paste0(cells, "_", time)
            message(resname)
            subsamples <- annotationTable$Sample[annotationTable$Cells ==
                                                     cells & annotationTable$Treatment_Duration ==
                                                     time & annotationTable$Treatment == "Oestrogen"]
            subraw <- rawcounts[, subsamples]
            subannot <- annotationTable[annotationTable$Sample %in%
                                            subsamples, c("Cells", "Condition", "Treatment")]
            rownames(subannot) <- annotationTable$Sample[annotationTable$Sample %in%
                                                             subsamples]
            subannot <- subannot[subsamples, ]
            dds <- DESeqDataSetFromMatrix(countData = subraw,
                                          colData = subannot, design = ~Condition)
            dds <- dds[rowSums(counts(dds)) > 1, ]
            dds$Condition <- relevel(dds$Condition, ref = "siCTRL")
            dea <- DESeq(dds, parallel = TRUE)
            res <- results(dea, contrast = c("Condition", "siZMIZ1",
                                             "siCTRL"))
            resannot <- cbind(rownames(res), eg2sym(rownames(res)))
            annotations <- annotategene(rownames(res))
            resannot <- cbind(as.matrix(resannot), as.matrix(annotations),
                              as.matrix(res))
            colnames(resannot)[1:3] <- c("ENTREZID", "SYMBOL",
                                         "NAME")
            resannot <- as.data.frame(resannot)
            resannot$log2FoldChange <- as.numeric(as.character(resannot$log2FoldChange))
            resannot$stat <- as.numeric(as.character(resannot$stat))
            resannot$pvalue <- as.numeric(as.character(resannot$pvalue))
            resannot$padj <- as.numeric(as.character(resannot$padj))
            resannot <- resannot[order(resannot$pvalue), ]
            result_list[[resname]] <- resannot
            rm(dea, resannot, res)
        }
    }


    #ZMIZ1 is knocked down.
    result_list[["MCF7_3h"]]["57178",]
    result_list[["MCF7_6h"]]["57178",]
    result_list[["MCF7_12h"]]["57178",]
    result_list[["MCF7_24h"]]["57178",]


    #ESR1 not signifantly different
    result_list[["MCF7_3h"]]["2099",]
    result_list[["MCF7_6h"]]["2099",]
    result_list[["MCF7_12h"]]["2099",]
    result_list[["MCF7_24h"]]["2099",]


    for (cells in c("MCF7")) {
        ZMIZ1expression<-matrix(ncol=5,nrow=0)
        colnames(ZMIZ1expression)<-c("Cells"  ,   "Condition","Treatment" ,"ZMIZ1"  ,   "Time"  )
        for (time in unique(annotationTable$Treatment_Duration)) {
            resname <- paste0(cells, "_", time)
            message(resname)
            subsamples <- annotationTable$Sample[annotationTable$Cells ==
                                                     cells & annotationTable$Treatment_Duration ==
                                                     time & annotationTable$Treatment == "Oestrogen"]
            subraw <- rawcounts[, subsamples]
            subrawZMIZ1<-subraw["57178",]
            subannot <- annotationTable[annotationTable$Sample %in%
                                            subsamples, c("Cells", "Condition", "Treatment")]
            rownames(subannot) <- annotationTable$Sample[annotationTable$Sample %in%
                                                             subsamples]
            subannot <- subannot[subsamples, ]
            ZMIZ1expression<-rbind(ZMIZ1expression,cbind(subannot,ZMIZ1=subrawZMIZ1[rownames(subannot)],Time=time))
        }
    }

    ZMIZ1expression<-ZMIZ1expression[ZMIZ1expression$Condition %in% c("siCTRL" , "siZMIZ1" ),]
    ZMIZ1expression<-cbind(ZMIZ1expression,"LogZMIZ1"=log(ZMIZ1expression$ZMIZ1 ))


    for (cells in c("MCF7")) {
        ESR1expression<-matrix(ncol=5,nrow=0)
        colnames(ESR1expression)<-c("Cells"  ,   "Condition","Treatment" ,"ESR1"  ,   "Time"  )
        for (time in unique(annotationTable$Treatment_Duration)) {
            resname <- paste0(cells, "_", time)
            message(resname)
            subsamples <- annotationTable$Sample[annotationTable$Cells ==
                                                     cells & annotationTable$Treatment_Duration ==
                                                     time & annotationTable$Treatment == "Oestrogen"]
            subraw <- rawcounts[, subsamples]
            subrawESR1<-subraw["2099",]
            subannot <- annotationTable[annotationTable$Sample %in%
                                            subsamples, c("Cells", "Condition", "Treatment")]
            rownames(subannot) <- annotationTable$Sample[annotationTable$Sample %in%
                                                             subsamples]
            subannot <- subannot[subsamples, ]
            ESR1expression<-rbind(ESR1expression,cbind(subannot,ESR1=subrawESR1[rownames(subannot)],Time=time))
        }
    }

    ESR1expression<-ESR1expression[ESR1expression$Condition %in% c("siCTRL" , "siZMIZ1" ),]
    ESR1expression<-cbind(ESR1expression,"LogESR1"=log(ESR1expression$ESR1 ))

    colnames(ESR1expression)[4]<-"RawCounts"
    colnames(ZMIZ1expression)[4]<-"RawCounts"
    colnames(ESR1expression)[6]<-"LogCounts"
    colnames(ZMIZ1expression)[6]<-"LogCounts"

    expression<-rbind(cbind(ESR1expression,Gene="ESR1"),cbind(ZMIZ1expression,Gene="ZMIZ1"))

    expression$Gene<-as.factor(expression$Gene)

    p <- ggboxplot(expression, x = "Time", y = "RawCounts", palette = c("#00AFBB",
                                            "#E7B800", "#FC4E07"), outlier.shape=NA, facet.by = "Gene", fill = "Condition")

    p <- p      + geom_point(aes(fill=factor(expression$Condition)),
                             size=1,
                             position = position_jitterdodge(dodge.width=0.8))   +ylab('Raw counts')

    p
}
