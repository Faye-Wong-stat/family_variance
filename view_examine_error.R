setwd("~/family_variance/")
library(abind)
library(ggplot2)
library(cowplot)



ind_names <- readRDS("create_marker_list/indiv_names.rds")
chrom_names <- readRDS("create_marker_list/chrom_names.rds")
phasing_error <- c(0, 0.005, 0.01, 0.025, 0.05, 0.1, 0.15)

score <- readRDS("examine_error/score.rds")

score_summ <- data.frame(error_level = rep(phasing_error, each=5), 
                         cM = rep(c(">100cM", ">10cM & <100cM", ">1cM & <10cM", 
                                    ">0.1cM & 1cM", ">0.01cM"), 7), 
                         num_compr = rep(NA, 7*5), 
                         num_corr = rep(NA, 7*5))

for (h in 1:length(phasing_error)){
  score_summ[,"num_compr"][((h-1)*5 + 1) : ((h-1)*5 + 5)] = 
    apply(score[[h]][3:4, , , ,], c(3), FUN=function(x){
      sum(!is.na(x)) 
    }) #+
    # apply(score[[2]][[h]][3:4, , , ,], c(3), FUN=function(x){
    #   sum(!is.na(x))
    # })
  
  score_summ[,"num_corr"][((h-1)*5 + 1) : ((h-1)*5 + 5)] = 
    apply(score[[h]][3:4, , , ,], c(3), FUN=function(x){
      sum(x, na.rm=T)
    }) #+ 
    # apply(score[[2]][[h]][3:4, , , ,], c(3), FUN=function(x){
    #   sum(x, na.rm=T)
    # })

  score[[h]] = NA 
}

# saveRDS(score_summ, "view_examine_error/score_summ.rds")
# score_summ <- readRDS("view_examine_error3/score_summ.rds")

score_summ$error_rate <- 1 - score_summ$num_corr/score_summ$num_compr
score_summ$cM <- factor(score_summ$cM, levels=c(">0.01cM", ">0.1cM & 1cM", ">1cM & <10cM", 
                                                ">10cM & <100cM", ">100cM"))
score_summ$error_level <- as.character(score_summ$error_level)



score_real <- readRDS("get_phase_error/score.rds")
# score_not_consc <- readRDS("codes/get_phase_error6_2_score.rds")

comp_real <- apply(score_real[3:4, , , , ], c(3), function(x){sum(!is.na(x))})
# comp_not_consc <- apply(score_not_consc[3:4, , , , ], c(3), function(x){sum(!is.na(x))})

error_real <- apply(score_real[3:4, , , , ], c(3), function(x){sum(x, na.rm=T)})
# error_not_consc <- apply(score_not_consc[3:4, , , , ], c(3), function(x){sum(x, na.rm=T)})

# real_esti_comp <- comp_consc + comp_not_consc
# real_esti_erro <- error_consc + error_not_consc

errors <- data.frame(error_level=rep("real_trios", 5), 
                     cM=c(">100cM", ">10cM & <100cM", ">1cM & <10cM", 
                          ">0.1cM & 1cM", ">0.01cM"), 
                     num_compr=comp_real, 
                     num_corr=error_real, 
                     error_rate=(1 - error_real/comp_real))
errors$cM <- factor(c(">100cM", ">10cM & <100cM", ">1cM & <10cM", 
                      ">0.1cM & 1cM", ">0.01cM"), 
                    levels=c(">0.01cM", ">0.1cM & 1cM", ">1cM & <10cM", 
                             ">10cM & <100cM", ">100cM"))

score_summ <- rbind(score_summ, errors)
score_summ$error_level <- as.factor(score_summ$error_level)

saveRDS(score_summ, "view_examine_error/score_summ.rds")
# score_summ <- readRDS("view_examine_error/score_summ.rds")



gg_color_hue <- function(n){
  hues = seq(15, 375, length=n + 1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols <- gg_color_hue(length(levels(score_summ$error_level)))
cols[levels(score_summ$error_level)=="real_trios"] <- "black"

p1 <- ggplot(data=score_summ, aes(x=cM, y=error_rate, group=error_level, color=error_level)) + 
  geom_line() + 
  xlab("pairwise distance") + ylab("fraction of inconsistent phasing") + 
  ylim(0, 0.5) + 
  scale_color_manual(values=cols) + 
  scale_x_discrete(guide=guide_axis(angle=45)) + 
  labs(colour="phase error rate \nper 1 cM") + 
  theme_minimal_grid(font_size=10)
  # guides(fill=guide_legend(title="phase error rate per 100cM"))
# p2 <- ggplot(data=errors, aes(x=cM, y=error_rate, group=1)) + 
#   geom_line() + 
#   facet_wrap(~consc) + 
#   xlab("log10 cM") + ylab("inconsistency") + 
#   ylim(0, 0.5) + 
#   theme_minimal_grid(font_size=10)
# prow <- plot_grid(p2, p1 + theme(legend.position="none"), labels="auto", align="h", axis="tb")
# legend <- get_legend(
#   p1
# )

# 6.5/(1+0.175)
# # [1] 5.531915
save_plot("view_examine_error/plots/introduced_error.pdf", 
          plot_grid(p1),
          base_width=6.5, base_height=4, base_asp=NULL)
# save_plot("view_examine_error/plots/testplot2.pdf", 
#           plot_grid(p1 + theme(legend.position="none"), legend, nrow=1, rel_widths=c(1, 0.175)),
#           base_width=6.5, base_height=5.53, base_asp=NULL)







