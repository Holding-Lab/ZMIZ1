#' Figure 1d
#'
#' This function allows you generate figure1d
#'
#' @keywords luciferase siZMIZ1 ER estrogen
#' @examples figure1d()
#' @export
#' @import ggpubr


figure1d <- function() {

    library(ggpubr)

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
