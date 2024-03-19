#' Figure 2c
#'
#' This function allows you generate the old figure2c - luciferase activity assay
#'
#' @keywords luciferase siZMIZ1 ER estrogen
#' @examples figure_removed_2c
#' @export
#' @import ggpubr
#' @importFrom utils read.table read.csv


figure_removed_2c <- function() {

    luciferaseData<-system.file("extdata",
                                "luciferase.csv",
                                package = "ZMIZ1")
    df_rep_cell<-read.csv(file=luciferaseData)

    p <- ggboxplot(df_rep_cell,
                   x="condition",
                   y="activity",
                   color="condition",
                   add="jitter",
                   outlier.shape=NA,
                   shape="condition",
                   palette =c("#00AFBB", "#E7B800", "#FC4E07","#5F7EF7"),
                   ylab="Log(Relative Luciferase Activity)",
                   xlab="Condition",
                   legend="none",
                   facet.by="cell"
    )
    p<-p + stat_compare_means(comparisons=list(c("siCTRL","siZMIZ1")),
                              geom="label",
                              method="t.test",
                              paired=T,
                              method.args=list(alternative="greater")
    )
    p




}
