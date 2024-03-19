#' Figure S2
#'
#' This function allows you generate figureS2 - PLA controls.
#'
#' @keywords PLA Analysis
#' @examples figureS2()
#' @export
#' @import ggpubr
#' @import rstatix
#' @import dplyr
#' @import tidyr
#' @importFrom utils read.table read.csv


figureS2 <- function() {

        #MCF7 (n=4)
        ER_ZMIZ1_NT <-c(171/13, 200/16, 127/18, 170/18)
        ER_ZMIZ1_siZMIZ1 <-c(0,0, 21/13, 19/13)
        ER_IgG_NT <- c(43/12, 26/12, 62/11, 19/15)
        ER_IgG_siZMIZ1 <-c(3/12,7/11,6/20,3/5)

        MCF7_PLA <- data.frame(ER_ZMIZ1_NT,ER_ZMIZ1_siZMIZ1,ER_IgG_NT,ER_IgG_siZMIZ1)
        MCF7_PLA_long <- pivot_longer(MCF7_PLA, cols = everything(), names_to = "Condition", values_to = "Value")

        #T47D (n=3)
        T4ER_ZMIZ1_NT <-c(108/12, 104/12, 150/14, NA)
        T4ER_ZMIZ1_siZMIZ1 <-c(140/23,46/19, 29/18, 71/25)
        T4ER_IgG_NT <- c(29/28, 61/34, 13/43, NA)
        T4ER_IgG_siZMIZ1 <-c(56/23,32/17,27/18, NA)

        T47D_PLA <- data.frame(T4ER_ZMIZ1_NT,T4ER_ZMIZ1_siZMIZ1,T4ER_IgG_NT,T4ER_IgG_siZMIZ1)
        T47D_PLA_long <- pivot_longer(T47D_PLA, cols = everything(), names_to = "Condition", values_to = "Value")
        T47D_PLA_long<-T47D_PLA_long[!is.na(T47D_PLA_long$Value),]

        T47D_PLA_long<-cbind(T47D_PLA_long,Cell='T47D')
        MCF7_PLA_long<-cbind(MCF7_PLA_long,Cell='MCF7')

        PLA_long<-rbind(MCF7_PLA_long,T47D_PLA_long)

        conditions<-do.call("rbind",strsplit(PLA_long$Condition,"_"))
        colnames(conditions)<-c("Factor1","Factor2","siRNA")
        conditions[,1]<-"ER"

        PLA_long<-cbind(PLA_long,conditions)

        orderStats <- c("MCF7", "T47D")
        PLA_long$Cell<-factor(PLA_long$Cell,levels = orderStats)

        #palette =c("#00AFBB", "#E7B800", "#FC4E07"),
        p<-ggboxplot(PLA_long, x = "Cell", y = "Value",
                     color = "Black", fill = "siRNA",
                     palette =c("#00AFBB", "#FC4E07"),
                     facet.by="Factor2",outlier.shape=NA ) +
                geom_point(aes(fill=factor(siRNA)),
                           size=0.6,
                           position = position_jitterdodge(dodge.width=0.8)
                )

        # Add p-values onto the box plots
        stat.test  <- PLA_long %>%
                group_by(Cell, Factor2) %>%
                t_test(Value ~ siRNA, alternative = "t" ) %>%
                adjust_pvalue(method = "bonferroni") %>%
                add_significance("p.adj")

        stat.test<-arrange(stat.test, Factor2, Cell)

        stat.test <- stat.test %>%
                add_xy_position(fun = "max", x = "Cell")

        stat.test$y.position<-c(0,0,25,15)

        Q<- p +   stat_pvalue_manual(
                stat.test, tip.length = 0.02,
                step.increase = 0.05,bracket.nudge.y = -3,
                label = "{p}",
                hide.ns = TRUE, size=3,
        )
        a<-Q+
                xlab("Cell Line") +ylab("Dots/Nuclea") +
                theme(axis.text.x = element_text(size=7),strip.text.x = element_text(size = 8),axis.text.y = element_text(size=8), axis.title.x=element_text(size=8.1),axis.title.y=element_text(size=8.1),
                      legend.title=element_text(size=8), legend.text=element_text(size=8)) +
                theme(legend.margin=margin(t = 0, unit='cm')) +
                theme(plot.margin = margin(0.25,0.25,0.25,0.25, "cm"))

        a

}
