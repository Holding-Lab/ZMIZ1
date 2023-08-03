#' Figure 3e
#'
#' This function allows you generate figure3e - E2F2 expression by timepoint
#'
#' @keywords siZMIZ1 E2F2 RNA timecourse
#' @examples figure3e()
#' @import DESeq2
#' @import ggpubr
#' @export

figure3e <- function() {


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
#E2F2
result_list[["MCF7_3h"]]["1870",]
result_list[["MCF7_6h"]]["1870",]
result_list[["MCF7_12h"]]["1870",]
result_list[["MCF7_24h"]]["1870",]




for (cells in c("MCF7")) {
  E2F2expression<-matrix(ncol=5,nrow=0)
  colnames(E2F2expression)<-c("Cells"  ,   "Condition","Treatment" ,"E2F2"  ,   "Time"  )
  for (time in unique(annotationTable$Treatment_Duration)) {
    resname <- paste0(cells, "_", time)
    message(resname)
    subsamples <- annotationTable$Sample[annotationTable$Cells ==
                                           cells & annotationTable$Treatment_Duration ==
                                           time & annotationTable$Treatment == "Oestrogen"]
    subraw <- rawcounts[, subsamples]
    subrawE2F2<-subraw["1870",]
    subannot <- annotationTable[annotationTable$Sample %in%
                                  subsamples, c("Cells", "Condition", "Treatment")]
    rownames(subannot) <- annotationTable$Sample[annotationTable$Sample %in%
                                                   subsamples]
    subannot <- subannot[subsamples, ]
    E2F2expression<-rbind(E2F2expression,cbind(subannot,E2F2=subrawE2F2[rownames(subannot)],Time=time))
  }
}

E2F2expression<-E2F2expression[E2F2expression$Condition %in% c("siCTRL" , "siZMIZ1" ),]
E2F2expression<-cbind(E2F2expression,"LogE2F2"=log(E2F2expression$E2F2 ))
colnames(E2F2expression)[4]<-"RawCounts"
colnames(E2F2expression)[6]<-"LogCounts"


expression<-rbind(cbind(E2F2expression,Gene="E2F2"))

expression$Gene<-as.factor(expression$Gene)

my_comparisons <- list( c(".siCTRL", "siZMIZ1"))


p <- ggline(expression, x = "Time", y = "RawCounts",
            palette =c("#00AFBB", "#E7B800", "#FC4E07"),
            color='Condition', add =c("mean_se","jitter"), ylab="E2F2 Counts")

p

}

