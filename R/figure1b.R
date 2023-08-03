#' Figure 1b
#'
#' This function allows you generate figure1b - qPLEX-RIME results.
#' @keywords qPlexRIME ER estrogen
#' @export
#' @import EnhancedVolcano
#' @importFrom utils read.table
#' @examples figure1b()


figure1b <- function() {



    qplexrimeData<-system.file("extdata", "MCF7_qPLEX-RIME.txt", package = "ZMIZ1")

    volcano<-read.table(qplexrimeData,sep="\t",header = T)

    df<-data.frame(volcano)

    df$GeneSymbol
    df$GeneSymbol<-as.character(df$GeneSymbol)

    EnhancedVolcano(df,
                    title = "MCF7, ER interactors +/- E2",
                    subtitle = "",
                    lab = df$GeneSymbol,
                    selectLab = c("HSP90AA1","HSP90AB1","FKBP5","FOXA1","ZMIZ1","RARA","ESR1","PIAS3","SUMO1","SUMO2","SUMO3","GATA3","NCOA3","EP300","MED16","MED8","GRHL2","MED24","RXRA","MED15"),
                    x = 'log2FC',
                    y = 'adj.P.Val',
                    xlim = c(-4, 4),
                    ylim=c(0,6),
                    FCcutoff=1,
                    pCutoff = 0.05,
                    pointSize = 1,
                    labSize = 5.0,
                    drawConnectors = TRUE)



}
