setwd("~/family_variance/")
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 512, 1024)
h2s <- c(0.9, 0.7, 0.5, 0.3, 0.1)
sf <- c(0.8, 0.9, 0.95, 0.98)
b <- qnorm(sf, 0, 1)
si <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)



fam_mean_sd <- readRDS("view_usefulness/fam_mean_sd.rds")
fam_mean_sd <- fam_mean_sd[, 1:43]
fam_mean_sd <- fam_mean_sd[fam_mean_sd$h2s==0.9, ][, 2:43]
fam_mean_sd$effective_marker_sizes <- as.numeric(as.character(fam_mean_sd$effective_marker_sizes))
fam_mean_sd$family <- as.character(fam_mean_sd$family)

fam_mean_sd_all <- data.frame(h2s=NA, number_of_QTL=NA, family=NA, replic=NA, mean=NA, sd=NA)
for (i in 1:nrow(fam_mean_sd)){
  fam_mean_sd_all = rbind(fam_mean_sd_all, 
                          data.frame(h2s=rep(NA, 20),
                                     number_of_QTL=rep(fam_mean_sd$effective_marker_sizes[i], 20), 
                                     family=rep(fam_mean_sd$family[i], 20), 
                                     replic=1:20, 
                                     mean=as.numeric(fam_mean_sd[i, 3:22][1, ]), 
                                     sd=as.numeric(fam_mean_sd[i, 23:42][1, ])))
}
fam_mean_sd_all <- fam_mean_sd_all[-1, -1]
fam_mean_sd_all[, c("parent1", "parent2")] <- 
  t(sapply(fam_mean_sd_all$family[530:535], FUN=function(x){
  strsplit(x, "_")[[1]]
}))
fam_mean_sd_all$parents_type <- "random"
fam_mean_sd_all <- fam_mean_sd_all[, c(1, 3, 6, 7, 8, 4, 5)]



fam_mean_sd_use <- fam_mean_sd_all[rep(1:nrow(fam_mean_sd_all), 4), ]
fam_mean_sd_use$si <- rep(si, each=nrow(fam_mean_sd_all))
fam_mean_sd_use$si_sd <- fam_mean_sd_use$si * fam_mean_sd_use$sd
fam_mean_sd_use$use <- fam_mean_sd_use$mean + fam_mean_sd_use$si * fam_mean_sd_use$sd



fam_cor <- data.frame(number_of_QTL=rep(as.character(effective_marker_sizes), each=4*20), 
                      si=rep(rep(paste("si = ", round(si, 2), sep=""), 6), each=20),
                      replic=rep(1:20, 4*6),
                      cor_use_mean=rep(NA, 4*6*20))

for (i in 1:length(effective_marker_sizes)){
  for (j in 1:length(si)){
    for (k in 1:20){
      temp = fam_mean_sd_use[(fam_mean_sd_use$number_of_QTL == effective_marker_sizes[i]) & 
                               (fam_mean_sd_use$si == si[j]) & 
                               (fam_mean_sd_use$replic == k), ]
      fam_cor[(i-1)*4*20 + (j-1)*20 + k, "cor_use_mean"] = cor(temp$mean, temp$use)
    }
  }
}
fam_cor$number_of_QTL <- factor(fam_cor$number_of_QTL, levels=as.factor(effective_marker_sizes))
# fam_cor$si <- factor(round(fam_cor$si, 2), levels=round(si, 2))
fam_cor$replic <- factor(fam_cor$replic, levels=1:20)

# pdf("examine_usefulness_real_data/plots/cor_use_mean_BV.pdf")
# ggplot(fam_cor, aes(x=as.numeric(number_of_QTL), y=cor_use_mean)) + 
#   geom_point() + 
#   geom_smooth() + 
#   facet_wrap(~si) + 
#   xlab("number of causal loci") + 
#   ylab("correlation between BV mean and usefulness") + 
#   scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) + 
#   theme_minimal_grid(font_size=10)
# dev.off()



fam_var_mean_var_sd <- aggregate(fam_mean_sd_use$mean, list(fam_mean_sd_use$number_of_QTL, 
                                                            fam_mean_sd_use$replic, 
                                                            fam_mean_sd_use$si), var)
colnames(fam_var_mean_var_sd) <- c("number_of_QTL", "replic", "si", "var_mean")
fam_var_mean_var_sd$var_si_sd <- aggregate(fam_mean_sd_use$si_sd, list(fam_mean_sd_use$number_of_QTL, 
                                                            fam_mean_sd_use$replic, 
                                                            fam_mean_sd_use$si), var)[, 4]
fam_var_mean_var_sd$si <- round(fam_var_mean_var_sd$si, 2)
fam_var_mean_var_sd$si <- paste("si = ", (fam_var_mean_var_sd$si), sep="")

# fam_var_mean_var_sd_4 <- fam_var_mean_var_sd[fam_var_mean_var_sd$number_of_QTL==4, ]
# summary(fam_var_mean_var_sd_4$var_mean)
# seqnc <- seq(0.1, 1.8, length.out=5)
# fam_var_mean_var_sd_4$group <- .bincode(fam_var_mean_var_sd_4$var_mean, seqnc)
# xlab <- cut(seq(0.1, 1.8, length.out=5), seqnc)
# xlab <- as.character(xlab)
# fam_var_mean_var_sd_4$xlab <- NA
# for (i in 1:nrow(fam_var_mean_var_sd_4)){
#   fam_var_mean_var_sd_4$xlab[i] = xlab[fam_var_mean_var_sd_4$group[i]+1]
# }
# fam_var_mean_var_sd_4$xlab <- factor(fam_var_mean_var_sd_4$xlab, levels=xlab)

fam_var_mean_var_sd_1024 <- fam_var_mean_var_sd[fam_var_mean_var_sd$number_of_QTL==1024, ]
summary(fam_var_mean_var_sd_1024$var_mean)
seqnc <- seq(100, 255, length.out=5)
fam_var_mean_var_sd_1024$group <- .bincode(fam_var_mean_var_sd_1024$var_mean, seqnc)
xlab <- cut(seq(100, 255, length.out=5), seqnc)
xlab <- as.character(xlab)
fam_var_mean_var_sd_1024$xlab <- NA
for (i in 1:nrow(fam_var_mean_var_sd_1024)){
  fam_var_mean_var_sd_1024$xlab[i] = xlab[fam_var_mean_var_sd_1024$group[i]+1]
}
fam_var_mean_var_sd_1024$xlab <- factor(fam_var_mean_var_sd_1024$xlab, levels=xlab)

# pdf("examine_usefulness_real_data/plots/cor_var_mean_var_si_sd_BV.pdf")
# for (i in 1:length(effective_marker_sizes)){
#   print(
#     ggplot(fam_var_mean_var_sd[fam_var_mean_var_sd$number_of_QTL == effective_marker_sizes[i], ], 
#            aes(x=var_mean, y=var_si_sd)) + 
#       geom_point() + 
#       # geom_smooth() + 
#       facet_wrap(~si) + 
#       xlab("var_mean") + 
#       ylab("var_si_sd") + 
#       ggtitle(paste("plot of variance of family SD against variance of family mean, \nnumber of effective QTL: ", effective_marker_sizes[i], sep="") )
#   )
# }
# dev.off()

# pdf("examine_usefulness_real_data/plots/cor_var_mean_var_si_sd_BV.pdf")
# ggplot(fam_var_mean_var_sd_1024, aes(xlab, var_si_sd)) + 
#   geom_boxplot() + 
#   facet_wrap(~si) + 
#   xlab("variance of mean") + 
#   ylab("variance of si*sd")
# dev.off()



p1 <- ggplot(fam_cor, aes(x=as.numeric(number_of_QTL), y=cor_use_mean)) + 
  geom_point() + 
  geom_smooth() + 
  facet_wrap(~si) + 
  xlab("number of causal loci") + 
  ylab("correlation between BV mean and usefulness") + 
  scale_x_continuous(breaks=1:6, labels=as.character(effective_marker_sizes)) +
  theme_minimal_grid(font_size=10)
p2 <- ggplot(fam_var_mean_var_sd_1024, aes(var_mean, var_si_sd)) + 
  geom_point() + 
  facet_wrap(~si) + 
  xlab("variance of mean") + 
  ylab("variance of si*sd") +
  theme_minimal_grid(font_size=10)

save_plot("examine_usefulness_real_data/plots/newplot.pdf", 
          plot_grid(p1, p2, nrow=1, labels="auto"),
          base_width=6.5)










