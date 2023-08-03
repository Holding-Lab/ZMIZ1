#' Figure 3c
#'
#' This function allows you generate figure3c - ZMIZ1 ChIP-qpcr at E2F2 promoter
#'
#' @keywords siZMIZ1 ChIP
#' @examples figure3c()
#' @import ggpubr
#' @export

figure3c <- function() {

  df<-data.frame(ct=c(0,-0.4514596571154,0,-0.8532381171768,0,-1.1463684138967),
                 fc=c(1,1.36742305725018,	1,1.80655116741969,	1,2.21355989985898),
                 rep=c(1,1,2,2,3,3),
                 Factor=c("IgG","ZMIZ1","IgG","ZMIZ1","IgG","ZMIZ1"))



  p <- ggbarplot(df, x = "Factor", y = "fc",
                 fill = "Factor", palette = "jco",
                 line.color = "gray", line.size = 0.4,
                 xlab="Factor",
                 ylab="Fold Change",
                 ylim=c(0,3),
                 add = c("mean_se", "jitter"))


  my_comparisons <- list( c("ZMIZ1", "IgG") )

  p+ stat_compare_means(method = "t.test",
                        paired=TRUE,
                        comparisons = my_comparisons,
                        method.args = list(alternative = "greater")
  ) + theme(legend.position = "none")
}

