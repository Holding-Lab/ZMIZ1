#' Figure S13mcf7
#'
#' This function allows you generate figure S13. GSEA results Part 2. P-value will vary as seed is
#' not set.
#' @keywords RNAseq ZMIZ1
#'
#' @export
#' @import DESeq2
#' @import vulcan
#' @importFrom graphics barplot
#' @importFrom stats relevel
#' @examples figureS13mcf7()


figureS13mcf7 <- function() {




    result_list <- list()
    for (cells in c("MCF7")) {
        for (time in c("6h")) {
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



    geneList <- result_list[["MCF7_6h"]]$log2FoldChange
    names(geneList) <- rownames(result_list[["MCF7_6h"]])

    williams <- msigdb[["c2_cgp;_;WILLIAMS_ESR1_TARGETS_UP"]]
    stein <- msigdb[["c2_cgp;_;STEIN_ESR1_TARGETS"]]
    bhat <- msigdb[["c2_cgp;_;BHAT_ESR1_TARGETS_NOT_VIA_AKT1_UP"]]

    go_cc <- msigdb[["c5_bp;_;GO_CELL_CYCLE"]]
    kegg_cell_cycle <- msigdb[["c2_cp;_;KEGG_CELL_CYCLE"]]
    react_cc <- msigdb[["c2_cpreactome;_;REACTOME_CELL_CYCLE"]]

    overlapped <- react_cc[react_cc %in% c(bhat, williams, stein)]

    overlapped_react <- ZMIZ1:::eg2sym(react_cc[react_cc %in% c(bhat, williams, stein)])
    overlapped_go <- ZMIZ1:::eg2sym(go_cc[go_cc  %in% c(bhat, williams, stein)])
    overlapped_kegg <- ZMIZ1:::eg2sym(kegg_cell_cycle[kegg_cell_cycle  %in% c(bhat, williams, stein)])

    obj <- gsea(sort(geneList, decreasing = TRUE), set = names(overlapped_react),
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h", title = "Reactome Cell Cycle/ER Response overlap")

    obj <- gsea(sort(geneList, decreasing = TRUE), set = names(overlapped_go),
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h", title = "Go Cell Cycle/ER Response overlap")

    obj <- gsea(sort(geneList, decreasing = TRUE), set = names(overlapped_kegg),
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h", title = "Kegg Cell Cycle/ER Response overlap")

}

#' Figure S13t47d
#'
#' This function allows you generate figure S3. GSEA results Part 2. P-value will vary as seed is
#' not set.
#' @keywords RNAseq ZMIZ1
#'
#' @export
#' @import DESeq2
#' @import vulcan
#' @importFrom graphics barplot
#' @importFrom stats relevel
#' @examples figureS13mcf7()


figureS13t47d <- function() {




    result_list <- list()
    for (cells in c("T47D")) {
        for (time in c("6h")) {
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



    geneList <- result_list[["T47D_6h"]]$log2FoldChange
    names(geneList) <- rownames(result_list[["T47D_6h"]])

    williams <- msigdb[["c2_cgp;_;WILLIAMS_ESR1_TARGETS_UP"]]
    stein <- msigdb[["c2_cgp;_;STEIN_ESR1_TARGETS"]]
    bhat <- msigdb[["c2_cgp;_;BHAT_ESR1_TARGETS_NOT_VIA_AKT1_UP"]]

    go_cc <- msigdb[["c5_bp;_;GO_CELL_CYCLE"]]
    kegg_cell_cycle <- msigdb[["c2_cp;_;KEGG_CELL_CYCLE"]]
    react_cc <- msigdb[["c2_cpreactome;_;REACTOME_CELL_CYCLE"]]

    overlapped <- react_cc[react_cc %in% c(bhat, williams, stein)]

    overlapped_react <- ZMIZ1:::eg2sym(react_cc[react_cc %in% c(bhat, williams, stein)])
    overlapped_go <- ZMIZ1:::eg2sym(go_cc[go_cc  %in% c(bhat, williams, stein)])
    overlapped_kegg <- ZMIZ1:::eg2sym(kegg_cell_cycle[kegg_cell_cycle  %in% c(bhat, williams, stein)])

    obj <- gsea(sort(geneList, decreasing = TRUE), set = names(overlapped_react),
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h", title = "Reactome Cell Cycle/ER Response overlap")

    obj <- gsea(sort(geneList, decreasing = TRUE), set = names(overlapped_go),
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h", title = "Go Cell Cycle/ER Response overlap")

    obj <- gsea(sort(geneList, decreasing = TRUE), set = names(overlapped_kegg),
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h", title = "Kegg Cell Cycle/ER Response overlap")

}

