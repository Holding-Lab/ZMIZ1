#' Figure 2a
#'
#' This function allows you generate figure2a - PLA results.
#'
#' @keywords PLA Analysis
#' @examples figure2a_PLA()
#' @export
#' @import ggpubr
#' @import rstatix
#' @import dplyr
#' @importFrom utils read.table read.csv


figure2a_PLA <- function() {

        zmiz1PLAfilename<-system.file("extdata",
                          "ZMIZ1plaData.csv",
                          package = "ZMIZ1")

        zmiz1PLAoutlier<-read.csv(zmiz1PLAfilename)

        #Remove point 15 as it's big outlier.
        zmiz1PLA<-zmiz1PLAoutlier[zmiz1PLAoutlier$Dots.Nuclea != zmiz1PLAoutlier$Dots.Nuclea[15],]


        orderStats <- c("MCF7", "T47D", "MDA-MB-231")

        zmiz1PLA$Cell.line<-factor(zmiz1PLA$Cell.line,levels = orderStats)

        # Add p-values onto the box plots
        stat.test  <- zmiz1PLA %>%
                group_by(Cell.line, Target.2) %>%
                t_test(Dots.Nuclea ~ Treatment, alternative = "t" ) %>%
                adjust_pvalue(method = "bonferroni") %>%
                add_significance("p.adj")
        stat.test <- stat.test %>%
                add_xy_position(fun = "mean_sd", x = "Cell.line")


        #stat.test$Cell.line <- factor(stat.test$Cell.line, levels = orderStats)
        stat.test<-arrange(stat.test, Target.2, Cell.line)


        p<-ggboxplot(zmiz1PLA, x = "Cell.line", y = "Dots.Nuclea",
                     color = "Black", fill = "Treatment", palette =c("#00AFBB", "#E7B800", "#FC4E07"),
                     facet.by="Target.2",outlier.shape=NA ) + geom_point(aes(fill=factor(Treatment)),
                                                                         size=0.6,
                                                                         position = position_jitterdodge(dodge.width=0.8)
                     )

        p

        Q<- p +   stat_pvalue_manual(
                stat.test, tip.length = 0.02,
                step.increase = 0.05,bracket.nudge.y = -3,
                label = "{p.adj.signif}",
                hide.ns = TRUE, size=3,
        ) +
                xlab("Cell Line") +ylab("Dots/Nuclea") +
                theme(axis.text.x = element_text(size=7),strip.text.x = element_text(size = 8),axis.text.y = element_text(size=8), axis.title.x=element_text(size=8.1),axis.title.y=element_text(size=8.1),
                      legend.title=element_text(size=8), legend.text=element_text(size=8)) +
                theme(legend.margin=margin(t = 0, unit='cm')) +
                theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"))  +
                scale_y_continuous(limits = c(0,42))

        a<-Q

        a


}
