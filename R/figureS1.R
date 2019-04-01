#' Figure S1
#'
#' This function allows you generate figure s1. P-value will vary as seed is
#' set.
#' @keywords RNAseq ZMIZ1
#'
#' @export
#' @import DESeq2
#' @import vulcan
#' @examples figures1()


figures1 <- function() {
    require(DESeq2)
    require(vulcan)


    result_list <- list()


    for (cells in c("MCF7")) {
        #unique(annotationTable$Cells)
        for (time in unique(annotationTable$Treatment_Duration)) {
            resname <- paste0(cells, "_", time)
            message(resname)

            # Subset Cells and Time  # Keep only oestrogen (ethanol is control)
            subsamples <-
                annotationTable$Sample[annotationTable$Cells == cells &
                            annotationTable$Treatment_Duration == time &
                               annotationTable$Treatment == "Oestrogen"]
            subraw <- rawcounts[, subsamples]
            subannot <-
                annotationTable[annotationTable$Sample %in% subsamples,
                                c("Cells", "Condition", "Treatment")]
            rownames(subannot) <-
                annotationTable$Sample[annotationTable$Sample %in% subsamples]
            subannot <- subannot[subsamples, ]

            # DESeq2 analysis
            dds <-
                DESeqDataSetFromMatrix(
                    countData = subraw,
                    colData = subannot,
                    design =  ~ Condition
                )
            dds <- dds[rowSums(counts(dds)) > 1, ]
            dds$Condition <- relevel(dds$Condition, ref = "siCTRL")
            dea <- DESeq(dds, parallel = TRUE)
            res <-
                results(dea, contrast = c("Condition", "siZMIZ1", "siCTRL"))

            # Annotate results
            resannot <- cbind(rownames(res), eg2sym(rownames(res)))
            annotations <- annotategene(rownames(res))
            resannot <-
                cbind(as.matrix(resannot),
                      as.matrix(annotations),
                      as.matrix(res))
            colnames(resannot)[1:3] <- c("ENTREZID", "SYMBOL", "NAME")
            resannot <- as.data.frame(resannot)
            resannot$log2FoldChange <-
                as.numeric(as.character(resannot$log2FoldChange))
            resannot$stat <- as.numeric(as.character(resannot$stat))
            resannot$pvalue <-
                as.numeric(as.character(resannot$pvalue))
            resannot$padj <- as.numeric(as.character(resannot$padj))
            resannot <- resannot[order(resannot$pvalue), ]

            result_list[[resname]] <- resannot
            rm(dea, resannot, res)
        }
    }



    ###Graph to find maxium effect.
    sig_genes <- matrix()
    sig_genes[1] <- sum(result_list[["MCF7_3h"]]$padj < 0.05, na.rm = T)
    sig_genes[2] <- sum(result_list[["MCF7_6h"]]$padj < 0.05, na.rm = T)
    sig_genes[3] <- sum(result_list[["MCF7_12h"]]$padj < 0.05, na.rm = T)
    sig_genes[4] <- sum(result_list[["MCF7_24h"]]$padj < 0.05, na.rm = T)

    names(sig_genes) <- c("3h", "6h", "12h", "24h")

    barplot(sig_genes,
            ylab = "Number of genes with signficantly different expresion",
            xlab = "Time after addition of E2",
            main = "Gene expression in MCF7 +/-siZMIZ1")


    williams <- msigdb[["c2_cgp;_;WILLIAMS_ESR1_TARGETS_UP"]]
    geneList <- result_list[["MCF7_6h"]]$log2FoldChange
    names(geneList) <- rownames(result_list[["MCF7_6h"]])
    obj <- gsea(sort(geneList, dec = TRUE),
                set = williams,
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h",
              title = "Williams ESR1 Targets")

    stein <- msigdb[["c2_cgp;_;STEIN_ESR1_TARGETS"]]
    obj <- gsea(sort(geneList, dec = TRUE),
                set = stein,
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h",
              title = "Stein ESR1 Targets")

    bhat <- msigdb[["c2_cgp;_;BHAT_ESR1_TARGETS_NOT_VIA_AKT1_UP"]]
    obj <- gsea(sort(geneList, dec = TRUE),
                set = bhat,
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h",
              title = "Bhat ESR1 Targets Up")


    ##Cell Cycle
    #GO
    go_cc <- msigdb[["c5_bp;_;GO_CELL_CYCLE"]]
    obj <- gsea(sort(geneList, dec = TRUE),
                set = go_cc,
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h",
              title = "GO: Cell Cycle")

    #Kegg
    kegg_cell_cycle <- msigdb[["c2_cp;_;KEGG_CELL_CYCLE"]]
    obj <-
        gsea(sort(geneList, dec = TRUE),
             set = kegg_cell_cycle,
             method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h",
              title = "Kegg Cell Cycle")

    #Reactome
    react_cc <- msigdb[["c2_cpreactome;_;REACTOME_CELL_CYCLE"]]
    obj <- gsea(sort(geneList, dec = TRUE),
                set = react_cc,
                method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h",
              title = "Reactome Cell Cycle")

    ##Overlap
    overlapped <- react_cc[react_cc %in% c(bhat, williams, stein)]
    obj <-
        gsea(sort(geneList, dec = TRUE),
             set = overlapped,
             method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h",
              title = "Reactome Cell Cycle/ER Response overlap")


    #Repead for T47D at just 6hours
    for (cells in c("T47D")) {
        #unique(annotation$Cells)
        for (time in c('6h')) {
            resname <- paste0(cells, "_", time)
            message(resname)

            # Subset Cells and Time  # Keep only oestrogen (ethanol is control)
            subsamples <-
                annotation$Sample[annotation$Cells == cells &
                                      annotation$Treatment_Duration == time &
                                      annotation$Treatment == "Oestrogen"]
            subraw <- rawcounts[, subsamples]
            subannot <-
                annotation[annotation$Sample %in% subsamples,
                           c("Cells", "Condition", "Treatment")]
            rownames(subannot) <-
                annotation$Sample[annotation$Sample %in% subsamples]
            subannot <- subannot[subsamples, ]

            # DESeq2 analysis
            dds <-
                DESeqDataSetFromMatrix(
                    countData = subraw,
                    colData = subannot,
                    design =  ~ Condition
                )
            dds <- dds[rowSums(counts(dds)) > 1, ]
            dds$Condition <- relevel(dds$Condition, ref = "siCTRL")
            dea <- DESeq(dds, parallel = TRUE)

            res <-
                results(dea, contrast = c("Condition", "siZMIZ1", "siCTRL"))

            # Annotate results
            resannot <- cbind(rownames(res), eg2sym(rownames(res)))
            annotations <- annotategene(rownames(res))
            resannot <-
                cbind(as.matrix(resannot),
                      as.matrix(annotations),
                      as.matrix(res))
            colnames(resannot)[1:3] <- c("ENTREZID", "SYMBOL", "NAME")
            resannot <- as.data.frame(resannot)
            resannot$log2FoldChange <-
                as.numeric(as.character(resannot$log2FoldChange))
            resannot$stat <- as.numeric(as.character(resannot$stat))
            resannot$pvalue <-
                as.numeric(as.character(resannot$pvalue))
            resannot$padj <- as.numeric(as.character(resannot$padj))
            resannot <- resannot[order(resannot$pvalue), ]

            result_list[[resname]] <- resannot
            rm(dea, resannot, res)
        }
    }


    ###Validate siZMIZ1/siCTRL are E2 genes at 6h in T47D

    geneList <- result_list[["T47D_6h"]]$log2FoldChange
    names(geneList) <- rownames(result_list[["T47D_6h"]])
    overlapped <- react_cc[react_cc %in% c(bhat, williams, stein)]
    obj <-
        gsea(sort(geneList, dec = TRUE),
             set = overlapped,
             method = "pareto")
    plot_gsea(obj, bottomYtitle = "siZMIZ/siCTRL at 6h",
              title = "Reactome Cell Cycle/ER Response overlap")




}
