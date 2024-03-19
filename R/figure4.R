#' Figure 4
#'
#' This function allows you generate figure2d - siZMIZ1 cell growth
#'
#' @keywords siZMIZ1 cell growth
#' @examples figure4_cellFrowth()
#' @import ggpubr
#' @import ggplot2
#' @import data.table
#' @export

figure4_cellGrowth <- function() {


    #//////////////////////////////////////////////////
    # MDA
    #/////////////////////////////////////////////////

    MDA.growth.1 <- read.table(system.file("extdata",
                                           "MDA_MB_231_growth.txt",
                                           package = "ZMIZ1"),
                               skip = 1,
                               header = TRUE)
    MDA.growth.rawdata <-
        MDA.growth.1[, c(1:3, seq(6, 15, by = 3), seq(5, 15, by = 3), seq(4, 15, by =
                                                                              3))]
    columnNames <-     c(
        "Date",
        "Time",
        "Elapsed",
        "siZMIZ1 1",
        "siZMIZ1 2",
        "siZMIZ1 3",
        "siZMIZ1 4",
        "siCTRL 1",
        "siCTRL 2",
        "siCTRL 3",
        "siCTRL 4",
        "RNAiMAX 1",
        "RNAiMAX 2",
        "RNAiMAX 3",
        "RNAiMAX 4"
    )


    colnames(MDA.growth.rawdata) <- columnNames


    #//////////////////////////////////////////////////
    # T47D
    #/////////////////////////////////////////////////
    T47D.growth.1 <- read.table(system.file("extdata",
                                            "T47D_growth.txt",
                                            package = "ZMIZ1"),
                                skip = 1,
                                header = TRUE)
    T47D.growth.rawdata <-
        T47D.growth.1[, c(1:3, seq(6, 15, by = 3), seq(5, 15, by = 3), seq(4, 15, by =
                                                                               3))] #, T47D.growth.2[,4:9])
    colnames(T47D.growth.rawdata) <- columnNames

    #//////////////////////////////////////////////////
    # MCF7
    #/////////////////////////////////////////////////
    MCF7.growth.1 <- read.table(system.file("extdata",
                                            "MCF7_growth.txt",
                                            package = "ZMIZ1"),
                                skip = 1,
                                header = TRUE)
    MCF7.growth.rawdata <- cbind(MCF7.growth.1)#, MCF7.growth.2[,4:9])
    MCF7.growth.rawdata <-
        MCF7.growth.1[, c(1:3, seq(6, 15, by = 3), seq(5, 15, by = 3), seq(4, 15, by =
                                                                               3))]

    colnames(MCF7.growth.rawdata) <- columnNames



    #MCF7 in ggplot format
    df <- melt(MCF7.growth.rawdata[, -1:-2], id = "Elapsed")
    df_split <-
        cbind(df[, c(1, 3)], do.call(rbind, strsplit(as.character(df$variable), split =
                                                         " ")))

    #T47D in ggplot format
    colnames(df_split) <- c("Elapsed", "Value", "Condition", "Rep")
    df_complete_MCF7 <- cbind(df_split, "MCF7")
    colnames(df_complete_MCF7)[5] <- "Cell line"

    df <- melt(T47D.growth.rawdata[, -1:-2], id = "Elapsed")
    df_split <-
        cbind(df[, c(1, 3)], do.call(rbind, strsplit(as.character(df$variable), split =
                                                         " ")))
    colnames(df_split) <- c("Elapsed", "Value", "Condition", "Rep")
    df_complete_T47D <- cbind(df_split, "T47D")
    colnames(df_complete_T47D)[5] <- "Cell line"

    #231s in ggplot format
    df <- melt(MDA.growth.rawdata[, -1:-2], id = "Elapsed")
    df_split <-
        cbind(df[, c(1, 3)], do.call(rbind, strsplit(as.character(df$variable), split =
                                                         " ")))
    colnames(df_split) <- c("Elapsed", "Value", "Condition", "Rep")
    df_complete_MDA <- cbind(df_split, "MDA")
    colnames(df_complete_MDA)[5] <- "Cell line"


    df_complete <-
        rbind(df_complete_T47D, df_complete_MCF7, df_complete_MDA)
    df_complete$Elapsed[df_complete$Elapsed == 1e-100] <- 0

    #df_complete$Elapsed <- as.numeric(as.character(df_complete$Elapsed))

    colnames(df_complete)[5] <- "Cellline"
    p <-
        ggline(
            df_complete,
            x = "Elapsed",
            xlab = "Elapsed/hours",
            ylab = "Confluence/percent",
            add = c("mean_se"),
            xlim = c(5, 110),
            y = "Value",
            color = "Condition",
            palette = "jco",
            facet.by = "Cellline"
        )
    p + scale_x_discrete(breaks = c(seq(21, 330, by = 30)), labels = c(seq(21, 330, by =
                                                                               30) - 21))



}
