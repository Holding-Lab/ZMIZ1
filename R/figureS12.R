#' Figure S12
#'
#' This function allows you generate figure S12 - E2 response of gene sets.
#'
#' @keywords RNAseq ZMIZ1 GSEA
#' @examples figureS12_E2_gene_set_volcano()
#' @import DESeq2
#' @import vulcan
#' @import EnhancedVolcano
#' @export

figureS12_E2_gene_set_volcano <- function() {


    result_list <- list()
    for (cells in c("MCF7")) {
        for (time in c("24h")) {
            resname <- paste0(cells, "_", time)
            message(resname)
            subsamples <- annotationTable$Sample[annotationTable$Cells ==
                                                     cells & annotationTable$Treatment_Duration ==
                                                     time & annotationTable$Condition == "siCTRL"]
            subraw <- rawcounts[, subsamples]
            subannot <- annotationTable[annotationTable$Sample %in%
                                            subsamples, c("Cells", "Condition", "Treatment")]
            rownames(subannot) <- annotationTable$Sample[annotationTable$Sample %in%
                                                             subsamples]
            subannot <- subannot[subsamples, ]
            dds <- DESeqDataSetFromMatrix(countData = subraw,
                                          colData = subannot, design = ~Treatment)
            dds <- dds[rowSums(counts(dds)) > 1, ]
            dds$Condition <- relevel(dds$Treatment, ref = "EtOH")
            dea <- DESeq(dds, parallel = TRUE)
            res <- results(dea, contrast = c("Treatment", "Oestrogen",
                                             "EtOH"))
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


    williams <- msigdb[["c2_cgp;_;WILLIAMS_ESR1_TARGETS_UP"]]
    stein <- msigdb[["c2_cgp;_;STEIN_ESR1_TARGETS"]]
    bhat <- msigdb[["c2_cgp;_;BHAT_ESR1_TARGETS_NOT_VIA_AKT1_UP"]]

    go_cc <- msigdb[["c5_bp;_;GO_CELL_CYCLE"]]
    kegg_cell_cycle <- msigdb[["c2_cp;_;KEGG_CELL_CYCLE"]]
    react_cc <- msigdb[["c2_cpreactome;_;REACTOME_CELL_CYCLE"]]



    overlapped_reactEG<-react_cc[react_cc %in% c(bhat, williams, stein)]
    resMCF7reactEG<-result_list$MCF7_24h[overlapped_reactEG,]

    overlapped_goEG<-go_cc[go_cc %in% c(bhat, williams, stein)]
    resMCF7goEG<-result_list$MCF7_24h[overlapped_goEG,]

    overlapped_keggEG<-kegg_cell_cycle[kegg_cell_cycle %in% c(bhat, williams, stein)]
    resMCF7keggEG<-result_list$MCF7_24h[overlapped_keggEG,]




    EnhancedVolcano(resMCF7reactEG,
                    lab = resMCF7reactEG$SYMBOL,
                    x = 'log2FoldChange',
                    y = 'padj',
                    title = 'Reactome / ER Cell Cycle Gene set',
                    drawConnectors = TRUE,
                    pCutoff = 0.05,
                    subtitle = 'MCF7 cells, 24 hours Estrogen treatment vs EtOH control',
                    selectLab = c('E2F2','FEN1','MCM2','MCM4','MCM6', 'PCNA', 'POLE2', 'MCM10', 'RFC4',
                                  'SKP2','CCNA2','CCNE2', 'CDC20','CDC25A','CDC6'))


    EnhancedVolcano(resMCF7goEG,
                    lab = resMCF7goEG$SYMBOL,
                    x = 'log2FoldChange',
                    y = 'padj',
                    pCutoff = 0.05,
                    title = 'GO / ER Cell Cycle Gene set',
                    subtitle = 'MCF7 cells, 24 hours Estrogen treatment vs EtOH control',
                    drawConnectors = TRUE,
                    selectLab = c('E2F2','FEN1','MCM2','MCM6', 'PCNA', 'POLE2', 'MCM10', 'RFC4',
                                  'SKP2','RRM1', 'CCNA2','CCNE2', 'CDC20','CDC25A','CDC6', 'E2F6', 'CENPU'))


    EnhancedVolcano(resMCF7keggEG,
                    lab = resMCF7keggEG$SYMBOL,
                    x = 'log2FoldChange',
                    y = 'padj',
                    pCutoff = 0.05,
                    drawConnectors = TRUE,
                    title = 'Kegg / ER Cell Cycle Gene set',
                    subtitle = 'MCF7 cells, 24 hours Estrogen treatment vs EtOH control')

}

