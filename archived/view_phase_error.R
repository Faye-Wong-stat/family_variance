setwd("~/family_variance/")
library(ggplot2)



score <- readRDS("get_phase_error/score.rds")
comp <- apply(score[3:4, , , , ], c(3), function(x){sum(!is.na(x))})
error <- apply(score[3:4, , , , ], c(3), function(x){sum(x, na.rm=T)})

# score_consc <- readRDS("codes/get_phase_error6_1_score.rds")
# score_not_consc <- readRDS("codes/get_phase_error6_2_score.rds")

# comp_consc <- apply(score_consc[3:4, , , , ], c(3), function(x){sum(!is.na(x))})
# comp_not_consc <- apply(score_not_consc[3:4, , , , ], c(3), function(x){sum(!is.na(x))})

# error_consc <- apply(score_consc[3:4, , , , ], c(3), function(x){sum(x, na.rm=T)})
# error_not_consc <- apply(score_not_consc[3:4, , , , ], c(3), function(x){sum(x, na.rm=T)})



errors <- data.frame(cM=c("2", "1", "0", "-1", "-2"), 
                     # consc=c(rep("consecutive", 5), rep("not consecutive", 5)), 
                     no.comp=comp, 
                     error=(1 - error/comp))
errors$cM <- factor(errors$cM, 
                    levels=c("-2", "-1", "0", "1", "2"))
pdf("view_phase_error/plots/error against cM.pdf")
ggplot(data=errors, aes(x=cM, y=error, group=1)) + 
  geom_line() + 
  # facet_wrap(~consc) + 
  xlab("log10 cM") + ylab("phasing error rate") + 
  ggtitle("phasing error increases as marker distance increases")
dev.off()
pdf("view_phase_error/plots/no.comparisons against cM.pdf")
ggplot(data=errors, aes(x=cM, y=no.comp)) + 
  geom_bar(stat="identity") + 
  # facet_wrap(~consc) + 
  xlab("log10 cM") + ylab("number of comparisons") + 
  ggtitle("number of comparisons agains cM distance")
dev.off()


