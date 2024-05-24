setwd("~/family_variance/")
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s <- c(0.8, 0.5, 0.2)
sf <- c(0.8, 0.99)
b <- qnorm(sf, 0, 1)
si <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)



results_cor_use <- readRDS("view_usefulness/results_cor_use.rds")
results_cor_use_best <- readRDS("view_usefulness_best_parents/results_cor_use.rds")

results_cor_use_mean <- matrix(NA, nrow=30, ncol=5)
results_cor_use_mean <- as.data.frame(results_cor_use_mean)
colnames(results_cor_use_mean) <- c("si", "h2s", "effective_marker_sizes", "mean", "se")
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      A = results_cor_use[results_cor_use$si == si[h] & 
                            results_cor_use$h2s == h2s[i] & 
                            results_cor_use$effective_marker_sizes == effective_marker_sizes[j], ]
      results_cor_use_mean[(h-1)*3*5 + (i-1)*5 + j, ] = 
        c(si[h], h2s[i], effective_marker_sizes[j], 
          mean(A$realuse_realmean_cor), sd(A$realuse_realmean_cor)/sqrt(20))
    }
  }
}

results_cor_use_mean_best <- matrix(NA, nrow=30, ncol=5)
results_cor_use_mean_best <- as.data.frame(results_cor_use_mean_best)
colnames(results_cor_use_mean_best) <- c("si", "h2s", "effective_marker_sizes", "mean", "se")
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      A = results_cor_use_best[results_cor_use_best$si == si[h] & 
                                 results_cor_use_best$h2s == h2s[i] & 
                                 results_cor_use_best$effective_marker_sizes == effective_marker_sizes[j], ]
      results_cor_use_mean_best[(h-1)*3*5 + (i-1)*5 + j, ] = 
        c(si[h], h2s[i], effective_marker_sizes[j], 
          mean(A$realuse_realmean_cor), sd(A$realuse_realmean_cor)/sqrt(20))
    }
  }
}

results_cor_use$`parent type` <- "random parents"
results_cor_use_best$`parent type` <- "best parents"
results_cor_use_mean$`parent type` <- "random parents"
results_cor_use_mean_best$`parent type` <- "best parents"
results_cor_use <- rbind(results_cor_use, results_cor_use_best)
results_cor_use_mean <- rbind(results_cor_use_mean, results_cor_use_mean_best)

results_cor_use$i <- paste("i == ", round(results_cor_use$si, 2), sep="")
results_cor_use$h2 <- paste("h^2 == ", results_cor_use$h2s, sep="")
results_cor_use$effective_marker_sizes <- factor(results_cor_use$effective_marker_sizes, 
                                                 levels=as.factor(effective_marker_sizes))
results_cor_use_mean$i <- paste("i == ", round(results_cor_use_mean$si, 2), sep="")
results_cor_use_mean$h2 <- paste("h^2 == ", results_cor_use_mean$h2s, sep="")
results_cor_use_mean$effective_marker_sizes <- factor(results_cor_use_mean$effective_marker_sizes, 
                                                 levels=as.factor(effective_marker_sizes))



p1 <- ggplot(results_cor_use_mean[results_cor_use_mean$h2s == 0.5, ], 
             aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=mean, colour=`parent type`)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=`parent type`), width=0.2) + 
  geom_line(aes(y=mean, colour=`parent type`), linewidth=0.5) + 
  facet_grid(~i, labeller = label_parsed) +
  xlab("number of causal loci") + 
  ylab("correlation between BV mean \nand usefulness") + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  scale_colour_manual(values=c("blue", "gold2")) + 
  theme_minimal_grid(font_size=8) + 
  theme(legend.position="bottom") 

save_plot("plot_BV_family_mean_sd/trans_cor_mean_use.pdf", 
          # prow,
          plot_grid(p1),
          base_width=6.5, base_height=3.3)

p2 <- ggplot(results_cor_use_mean, aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=mean, colour=`parent type`)) + 
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, colour=`parent type`), width=0.2) + 
  geom_line(aes(y=mean, colour=`parent type`), linewidth=0.5) + 
  facet_grid(i~h2, labeller = label_parsed) +
  xlab("number of causal loci") + 
  ylab("correlation between BV mean \nand usefulness") + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  scale_colour_manual(values=c("blue", "gold2")) + 
  theme_minimal_grid(font_size=8) + 
  theme(legend.position="bottom") 

save_plot("plot_BV_family_mean_sd/trans_cor_mean_use2.pdf", 
          # prow,
          plot_grid(p2),
          base_width=6.5, base_height=4.33)







fam_mean_sd <- readRDS("view_usefulness/fam_mean_sd.rds")
fam_mean_sd <- fam_mean_sd[, 1:43]
# any(fam_mean_sd[fam_mean_sd$h2s==0.9, ][, 2:43] != fam_mean_sd[fam_mean_sd$h2s==0.7, ][, 2:43], 
#     fam_mean_sd[fam_mean_sd$h2s==0.9, ][, 2:43] != fam_mean_sd[fam_mean_sd$h2s==0.5, ][, 2:43], 
#     fam_mean_sd[fam_mean_sd$h2s==0.9, ][, 2:43] != fam_mean_sd[fam_mean_sd$h2s==0.3, ][, 2:43], 
#     fam_mean_sd[fam_mean_sd$h2s==0.9, ][, 2:43] != fam_mean_sd[fam_mean_sd$h2s==0.1, ][, 2:43])
fam_mean_sd <- fam_mean_sd[fam_mean_sd$h2s==0.9, ][, 2:43]
fam_mean_sd$effective_marker_sizes <- as.numeric(as.character(fam_mean_sd$effective_marker_sizes))
fam_mean_sd$family <- as.character(fam_mean_sd$family)

fam_mean_sd_all <- data.frame(number_of_QTL=NA, family=NA, replic=NA, BV_fammean=NA, BV_famsd=NA)
for (i in 1:nrow(fam_mean_sd)){
  fam_mean_sd_all[(nrow(fam_mean_sd_all) + 1) : (nrow(fam_mean_sd_all) + 20), ] = 
    data.frame(number_of_QTL=rep(fam_mean_sd$effective_marker_sizes[i], 20), 
               family=rep(fam_mean_sd$family[i], 20), 
               replic=1:20, 
               BV_fammean=as.numeric(fam_mean_sd[i, 3:22][1, ]), 
               BV_famsd=as.numeric(fam_mean_sd[i, 23:42][1, ]))
} 
fam_mean_sd_all <- fam_mean_sd_all[-1, ]
fam_mean_sd_all$parent1 <- strsplit(fam_mean_sd_all$family, "_")[[1]][1]
fam_mean_sd_all$parent2 <- strsplit(fam_mean_sd_all$family, "_")[[1]][2]

set.seed(1)
random_190_par_names <- sample(unique(fam_mean_sd_all$family), 190)
random_190_par_names[1:20]
# [1] "12C014P001_16C087P017"  "NA.1_16C096P025"        "87C112P010_16C036P022" 
# [4] "16C054P020_16C091P001"  "16C045P017_16C014P006"  "12C178P005_16C098P025" 
# [7] "16C058P010_16C005P007"  "16C046P037_65C272P104"  "16C090P015_09C061P602" 
# [10] "16C565P028_16C033P021"  "97C192P001_16C517P034"  "16C067P036_57C004P001" 
# [13] "16C042P016_01C206P005"  "16C067P027_PI270471"    "11C032P0021_16C089P012"
# [16] "16C017P024_12C112P004"  "08C043P001_16C565P028"  "11C113P0061_16C546P011"
# [19] "16C018P011_16C104P029"  "16C092P005_12C072P001" 

fam_mean_sd_190 <- fam_mean_sd_all[fam_mean_sd_all$family %in% random_190_par_names, ]

fam_mean_sd_all <- fam_mean_sd_all[, c(1, 3, 6, 7, 4, 5)]
fam_mean_sd_all$parents_type <- "random"

fam_mean_sd_190 <- fam_mean_sd_190[, c(1, 3, 6, 7, 4, 5)]
fam_mean_sd_190$parents_type <- "random 190"

fam_mean_sd_all <- rbind(fam_mean_sd_all, fam_mean_sd_190)


# info_df <- readRDS("view_crosses_good_parents/info_df.rds")
# info_df <- info_df[, c(1:5, 7, 8)]
# # any(info_df[info_df$h2s==0.9, ][, 2:7] != info_df[info_df$h2s==0.7, ][, 2:7],
# #     info_df[info_df$h2s==0.9, ][, 2:7] != info_df[info_df$h2s==0.5, ][, 2:7],
# #     info_df[info_df$h2s==0.9, ][, 2:7] != info_df[info_df$h2s==0.3, ][, 2:7],
# #     info_df[info_df$h2s==0.9, ][, 2:7] != info_df[info_df$h2s==0.1, ][, 2:7])
# info_df$parents_type <- "good parents predicted from models"
# 
# fam_mean_sd_all <- rbind(fam_mean_sd_all, info_df)



info_df <- readRDS("view_crosses_best_parents/info_df.rds")
info_df <- info_df[, 1:6]
info_df$parents_type <- "best parents"

fam_mean_sd_all <- rbind(fam_mean_sd_all, info_df)



# info_df <- readRDS("view_crosses_medium_parents/info_df.rds")
# info_df <- info_df[, 1:6]
# info_df$parents_type <- "medium parents"
# 
# fam_mean_sd_all <- rbind(fam_mean_sd_all, info_df)


unique(fam_mean_sd_all$parents_type)
fam_mean_sd_all$parents_type <- factor(fam_mean_sd_all$parents_type, 
                                       levels=unique(fam_mean_sd_all$parents_type))
fam_mean_sd_all$number_of_QTL <- paste("number of QTL = ", fam_mean_sd_all$number_of_QTL, sep="")
fam_mean_sd_all$number_of_QTL <- factor(fam_mean_sd_all$number_of_QTL, 
                                        levels=paste("number of QTL = ", effective_marker_sizes, sep=""))

saveRDS(fam_mean_sd_all, "plot_BV_family_mean_sd/fam_mean_sd_all.rds")
# fam_mean_sd_all <- readRDS("plot_BV_family_mean_sd/fam_mean_sd_all.rds")
dim(fam_mean_sd_all)
# [1] 108000      7
6*20*190*2 + 6*20*520
# [1] 108000


fam_mean_sd_all_si <- do.call("rbind", replicate(4, fam_mean_sd_all, simplify=F))
fam_mean_sd_all_si$si <- rep(si, each=dim(fam_mean_sd_all)[1])
fam_mean_sd_all_si$use <- 
  fam_mean_sd_all_si$BV_fammean + fam_mean_sd_all_si$si*fam_mean_sd_all_si$BV_famsd

fam_mean_sd_all_summ_cor <- sapply(split(fam_mean_sd_all_si, 
                                         list(fam_mean_sd_all_si$number_of_QTL, 
                                              fam_mean_sd_all_si$replic, 
                                              fam_mean_sd_all_si$parents_type, 
                                              fam_mean_sd_all_si$si)), 
                                   FUN=function(x){
                                     c(var(x$BV_fammean), 
                                       var(x$BV_famsd), 
                                       cov(x$BV_fammean, x$BV_famsd), 
                                       cor(x$BV_fammean, x$us))
                                   })
fam_mean_sd_all_summ_cor <- t(fam_mean_sd_all_summ_cor)
colnames(fam_mean_sd_all_summ_cor) <- c("var_mean", "var_sd", "cov_mean_sd", "cor_mean_use")
fam_mean_sd_all_summ_cor <- as.data.frame(fam_mean_sd_all_summ_cor)
fam_mean_sd_all_summ_cor[, c("number_of_QTL", "replic", "parents_type", "si")] <- 
  as.data.frame(t(sapply(rownames(fam_mean_sd_all_summ_cor), 
                         FUN=function(x){
                           c(strsplit(x, "\\.")[[1]][1:3], 
                             paste(strsplit(x, "\\.")[[1]][4:5], collapse="."))
                         })))
fam_mean_sd_all_summ_cor$cor_fishertrans <- 
  0.5 * log((1+fam_mean_sd_all_summ_cor$cor_mean_use) / (1-fam_mean_sd_all_summ_cor$cor_mean_use))

fam_mean_sd_all_summ_cor_agg <- aggregate(cor_fishertrans ~ 
                                            number_of_QTL + parents_type + si, 
                                          data=fam_mean_sd_all_summ_cor, 
                                          FUN=function(x){
                                            c(mean=mean(x), se=sd(x)/sqrt(length(x)))
                                          })
fam_mean_sd_all_summ_cor_agg[, 4:5] <- data.frame(fam_mean_sd_all_summ_cor_agg[, 4])
fam_mean_sd_all_summ_cor_agg$number_of_QTL <- 
  factor(fam_mean_sd_all_summ_cor_agg$number_of_QTL, 
         levels=paste("number of QTL = ", effective_marker_sizes, sep=""))
fam_mean_sd_all_summ_cor_agg$parents_type <- 
  factor(fam_mean_sd_all_summ_cor_agg$parents_type, levels=parents_types)
fam_mean_sd_all_summ_cor_agg$si <- as.numeric(fam_mean_sd_all_summ_cor_agg$si)
fam_mean_sd_all_summ_cor_agg$upper_trans <- 
  fam_mean_sd_all_summ_cor_agg$cor_fishertrans + 2*fam_mean_sd_all_summ_cor_agg$se
fam_mean_sd_all_summ_cor_agg$lower_trans <- 
  fam_mean_sd_all_summ_cor_agg$cor_fishertrans - 2*fam_mean_sd_all_summ_cor_agg$se
fam_mean_sd_all_summ_cor_agg$mean <- 
  (exp(2*fam_mean_sd_all_summ_cor_agg$cor_fishertrans) - 1) / 
  (exp(2*fam_mean_sd_all_summ_cor_agg$cor_fishertrans) + 1)
fam_mean_sd_all_summ_cor_agg$upper <- 
  (exp(2*fam_mean_sd_all_summ_cor_agg$upper_trans) - 1) / 
  (exp(2*fam_mean_sd_all_summ_cor_agg$upper_trans) + 1)
fam_mean_sd_all_summ_cor_agg$lower <- 
  (exp(2*fam_mean_sd_all_summ_cor_agg$lower_trans) - 1) / 
  (exp(2*fam_mean_sd_all_summ_cor_agg$lower_trans) + 1)

A <- fam_mean_sd_all_summ_cor_agg[fam_mean_sd_all_summ_cor_agg$number_of_QTL == "number of QTL = 4" & 
                                    fam_mean_sd_all_summ_cor_agg$parents_type == "best parents", ]
A
# number_of_QTL parents_type       si cor_fishertrans        se
# 4  number of QTL = 4 best parents 1.399810       0.5387682 0.2048167
# 22 number of QTL = 4 best parents 1.754983       0.1448880 0.2178758
# 40 number of QTL = 4 best parents 2.062713      -0.1005408 0.2237328
# 58 number of QTL = 4 best parents 2.420907      -0.3009145 0.2291969
# upper_trans lower_trans       mean     upper      lower
# 4    0.9484015   0.1291348  0.4920550 0.7390585  0.1284218
# 22   0.5806396  -0.2908637  0.1438825 0.5231301 -0.2829295
# 40   0.3469247  -0.5480063 -0.1002034 0.3336454 -0.4990245
# 58   0.1574793  -0.7593083 -0.2921493 0.1561903 -0.6406694



pdf(paste("plot_BV_family_mean_sd/", "trans_cor_mean_use.pdf", sep=""), width=10)
ggplot(fam_mean_sd_all_summ_cor_agg, aes(x=si)) + 
  geom_line(aes(y=mean, colour=parents_type)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, 
                  colour=parents_type, fill=parents_type), alpha=0.3, linetype=0) + 
  facet_wrap(~number_of_QTL) + 
  xlab("selection intensity") + 
  ylab("correlation between BV family mean and BV family usefulness") + 
  ylim(0, 1)
dev.off()

gg_color_hue <- function(n){
  hues = seq(15, 375, length=n + 1)
  hcl(h=hues, l=65, c=100)[1:n]
}
cols <- gg_color_hue(length(levels(fam_mean_sd_all_summ_cor_agg$parents_type)))

p1 <- ggplot(fam_mean_sd_all_summ_cor_agg, aes(x=si)) + 
  geom_line(aes(y=mean, colour=parents_type)) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, 
                  colour=parents_type, fill=parents_type), alpha=0.3, linetype=0) + 
  facet_wrap(~number_of_QTL) + 
  xlab("selection intensity") + 
  ylab("correlation between family BV mean and family BV UC") + 
  coord_cartesian(ylim = c(0, 1)) + 
  theme_minimal_grid(font_size=10)

save_plot("plot_BV_family_mean_sd/trans_cor_mean_use.pdf", 
          p1,
          base_width=6.5, base_height = 4.33)



# don't run the rest of the script
Zalphas <- readRDS("simulate_phenotypes_real_data/Zalphas.rds")

pdf(paste("plot_BV_family_mean_sd/", "hist_BV.pdf", sep=""))
hist(Zalphas[[4]][[2]], breaks=100, 
     xlab="breeding values", ylab="counts", 
     main="histogram of BV of potential parents")
dev.off()












fam_mean_sd_all_summ_si <- aggregate(cbind(BV_fammean, use) ~ 
                                       number_of_QTL + replic + parents_type + si, 
                                     data=fam_mean_sd_all_si, FUN=var)
colnames(fam_mean_sd_all_summ_si)[5:6] <- c("var_of_mean", "var_of_use")
fam_mean_sd_all_summ_si$ratio <- fam_mean_sd_all_summ_si$var_of_mean / fam_mean_sd_all_summ_si$var_of_use

fam_mean_sd_all_summ_si <- 
  fam_mean_sd_all_summ_si[(!is.na(fam_mean_sd_all_summ_si$ratio)) & 
                            (fam_mean_sd_all_summ_si$ratio<1), ]

fam_mean_sd_all_summ_si$ratio_logit <- 
  log(fam_mean_sd_all_summ_si$ratio / (1-fam_mean_sd_all_summ_si$ratio))

sum(is.na(fam_mean_sd_all_summ_si$ratio))
# [1] 396
sum(fam_mean_sd_all_summ_si$ratio>1, na.rm=T)
# [1] 17740
# too many ratios > 1

fam_mean_sd_all_summ_si_agg <- 
  aggregate(fam_mean_sd_all_summ_si$ratio_logit, 
            list(fam_mean_sd_all_summ_si$number_of_QTL, 
                 fam_mean_sd_all_summ_si$parents_type, 
                 fam_mean_sd_all_summ_si$si), 
            FUN=function(x){
              c(mean(x), sd(x)/sqrt(length(x)))
            })




# fam_mean_sd_all_summ_2 <- data.frame(number_of_QTL=NA, replic=NA, parents_type=NA, 
#                                      var_of_mean=NA, var_of_sd=NA, cov_of_mean_sd=NA)
# for (i in 1:length(effective_marker_sizes)){
#   for (j in 1:20){
#     for (k in 1:length(parents_types)){
#       temp = fam_mean_sd_all[fam_mean_sd_all$number_of_QTL==effective_marker_sizes[i] & 
#                                fam_mean_sd_all$replic==j & 
#                                fam_mean_sd_all$parents_type==parents_types[k], ]
#       fam_mean_sd_all_summ_2 = rbind(fam_mean_sd_all_summ_2, 
#                                      data.frame(number_of_QTL=effective_marker_sizes[i], 
#                                                 replic=j, 
#                                                 parents_type=parents_types[k], 
#                                                 var_of_mean=var(temp$BV_fammean), 
#                                                 var_of_sd=var(temp$BV_famsd), 
#                                                 cov_of_mean_sd=cov(temp$BV_fammean, temp$BV_famsd)))
#     }
#   }
# }
# fam_mean_sd_all_summ_2 <- fam_mean_sd_all_summ_2[-1, ]

# dim(split(fam_mean_sd_all, list(fam_mean_sd_all$number_of_QTL, 
#                                    fam_mean_sd_all$replic, 
#                                    fam_mean_sd_all$parents_type))[[1]])
# split(fam_mean_sd_all, list(fam_mean_sd_all$number_of_QTL, 
#                             fam_mean_sd_all$replic, 
#                             fam_mean_sd_all$parents_type))[[1]][1:10, ]


fam_mean_sd_all_summ <- sapply(split(fam_mean_sd_all, list(fam_mean_sd_all$number_of_QTL, 
                                                           fam_mean_sd_all$replic, 
                                                           fam_mean_sd_all$parents_type)), 
                               FUN=function(x){
                                 c(var(x$BV_fammean), 
                                   var(x$BV_famsd), 
                                   cov(x$BV_fammean, x$BV_famsd), 
                                   cor(x$BV_fammean, x$BV_famsd))
                               })
fam_mean_sd_all_summ <- t(fam_mean_sd_all_summ)
colnames(fam_mean_sd_all_summ) <- c("var_of_mean", "var_of_sd", "cov_of_mean_sd", "cor_of_mean_sd")
fam_mean_sd_all_summ <- as.data.frame(fam_mean_sd_all_summ)
fam_mean_sd_all_summ[c("number_of_QTL", "replic", "parents_type")] <- 
  as.data.frame(t(sapply(rownames(fam_mean_sd_all_summ), FUN=function(x){
    strsplit(x, "\\.")[[1]]
  })))

# fam_mean_sd_all_summ_funct_2 <- do.call("rbind", replicate(
#   99, fam_mean_sd_all_summ_2, simplify=F
# ))

fam_mean_sd_all_summ_funct <- do.call("rbind", replicate(
  4, fam_mean_sd_all_summ, simplify=F
))
# for (i in 2:99){
#   print(c(i, ((i-1)*480+1), 
#           sum((fam_mean_sd_all_summ_funct[1:480, ]) != 
#                 (fam_mean_sd_all_summ_funct[((i-1)*480+1) : ((i-1)*480+480), ]))))
# }
fam_mean_sd_all_summ_funct$si <- rep(si, each=(6*20*4))
fam_mean_sd_all_summ_funct$f <- 
  ((fam_mean_sd_all_summ_funct$si^2) * (fam_mean_sd_all_summ_funct$var_of_sd)) / 
  ((fam_mean_sd_all_summ_funct$si^2) * (fam_mean_sd_all_summ_funct$var_of_sd) + 
     (fam_mean_sd_all_summ_funct$var_of_mean) + 
     (2 * fam_mean_sd_all_summ_funct$si * fam_mean_sd_all_summ_funct$cov_of_mean_sd))
# fam_mean_sd_all_summ_funct$use <- 
  
# fam_mean_sd_all_summ_funct <- fam_mean_sd_all_summ_funct[!is.na(fam_mean_sd_all_summ_funct$f), ]

# fam_mean_sd_all_summ_funct_2$si <- rep(si, each=(6*20*4))
# fam_mean_sd_all_summ_funct_2$f <- 
#   ((fam_mean_sd_all_summ_funct_2$si^2) * (fam_mean_sd_all_summ_funct_2$var_of_sd)) / 
#   ((fam_mean_sd_all_summ_funct_2$si^2) * (fam_mean_sd_all_summ_funct_2$var_of_sd) + 
#      (fam_mean_sd_all_summ_funct_2$var_of_mean) + 
#      (2 * fam_mean_sd_all_summ_funct_2$si * fam_mean_sd_all_summ_funct_2$cov_of_mean_sd))
# fam_mean_sd_all_summ_funct_2 <- fam_mean_sd_all_summ_funct_2[!is.na(fam_mean_sd_all_summ_funct_2$f), ]
# 
# fam_mean_sd_all_summ_funct_2_ratio_big <- 
#   fam_mean_sd_all_summ_funct_2[which(fam_mean_sd_all_summ_funct_2$f>1), ]

# dim(fam_mean_sd_all_summ_funct[fam_mean_sd_all_summ_funct$cov_of_mean_sd==0, ])
# dim(fam_mean_sd_all_summ_funct[fam_mean_sd_all_summ_funct$f>1, ])

fam_mean_sd_all_summ_funct_ratio_big <- 
  fam_mean_sd_all_summ_funct[which(fam_mean_sd_all_summ_funct$f>1), ]

# position_1 <- 
#   fam_mean_sd_all_summ_funct_ratio_big[1, 
#                                        c("number_of_QTL", "replic", "parents_type", "si")]
# fam_mean_sd_all_position_1 <- 
#   fam_mean_sd_all[fam_mean_sd_all$number_of_QTL=="4" & fam_mean_sd_all$replic==5 & 
#                     fam_mean_sd_all$parents_type=="medium parents", ]
# cov(fam_mean_sd_all_position_1$BV_fammean, fam_mean_sd_all_position_1$BV_famsd)

# pdf(paste("plot_BV_family_mean_sd/", "hist_cov_of_mean_sd.pdf", sep=""))
# hist(fam_mean_sd_all_summ_funct$cov_of_mean_sd, breaks=200)
# dev.off()

# fam_mean_sd_all_summ_funct$f_logit

fam_mean_sd_all_summ_funct_plot <- aggregate(fam_mean_sd_all_summ_funct$f, 
                                             list(fam_mean_sd_all_summ_funct$number_of_QTL, 
                                                  fam_mean_sd_all_summ_funct$parents_type, 
                                                  fam_mean_sd_all_summ_funct$si), 
                                             mean)
colnames(fam_mean_sd_all_summ_funct_plot) <- c("number_of_QTL", "parents_type", "si", "mean")
fam_mean_sd_all_summ_funct_plot$sd <- aggregate(fam_mean_sd_all_summ_funct$f, 
                                               list(fam_mean_sd_all_summ_funct$number_of_QTL, 
                                                    fam_mean_sd_all_summ_funct$parents_type, 
                                                    fam_mean_sd_all_summ_funct$si), 
                                               sd)[, 4]
fam_mean_sd_all_summ_funct_plot$number_of_QTL <- 
  factor(fam_mean_sd_all_summ_funct_plot$number_of_QTL, levels=effective_marker_sizes)

pdf(paste("plot_BV_family_mean_sd/", "ratio_var_i*sd_var_use_cont.pdf", sep=""))
ggplot(fam_mean_sd_all_summ_funct_plot, 
       aes(x=si, colour=parents_type, fill=parents_type)) + 
  geom_line(aes(y=mean)) + 
  geom_ribbon(aes(ymin=mean-2*sd, ymax=mean+2*sd), alpha=0.3, linetype=0) +
  facet_wrap(~number_of_QTL)
  
dev.off()






















# fam_mean_sd_all$`0.25_i*sd` <- 0.25*fam_mean_sd_all$BV_famsd
# fam_mean_sd_all$`0.50_i*sd` <- 0.5*fam_mean_sd_all$BV_famsd
# fam_mean_sd_all$`0.75_i*sd` <- 0.75*fam_mean_sd_all$BV_famsd
# 
# fam_mean_sd_all$`0.25_use` <- fam_mean_sd_all$BV_fammean + 0.25*fam_mean_sd_all$BV_famsd
# fam_mean_sd_all$`0.50_use` <- fam_mean_sd_all$BV_fammean + 0.5*fam_mean_sd_all$BV_famsd
# fam_mean_sd_all$`0.75_use` <- fam_mean_sd_all$BV_fammean + 0.75*fam_mean_sd_all$BV_famsd



# fam_mean_sd_all_var <- aggregate(fam_mean_sd_all[, 8], list(fam_mean_sd_all$number_of_QTL, 
#                                                                 fam_mean_sd_all$replic, 
#                                                                 fam_mean_sd_all$parents_type), 
#                                  var)
# for (i in 9:13){
#   fam_mean_sd_all_var = cbind(fam_mean_sd_all_var, 
#                               aggregate(fam_mean_sd_all[, i], list(fam_mean_sd_all$number_of_QTL, 
#                                                                        fam_mean_sd_all$replic, 
#                                                                        fam_mean_sd_all$parents_type), 
#                                         var)[, 4])
# }
# 
# colnames(fam_mean_sd_all_var) <- c("number_of_QTL", "replic", "parents_type", 
#                                    "0.25_i*sd", "0.50_i*sd", "0.75_i*sd", 
#                                    "0.25_use", "0.50_use", "0.75_use")
# 
# fam_mean_sd_all_var_long <- data.frame(number_of_QTL=NA, replic=NA, parents_type=NA, si=NA, 
#                                       `i*sd_var`=NA, `use_var`=NA)
# for (i in 1:nrow(fam_mean_sd_all_var)){
#   fam_mean_sd_all_var_long = rbind(fam_mean_sd_all_var_long, 
#                                   data.frame(number_of_QTL=rep(fam_mean_sd_all_var[i, 1], 3), 
#                                              replic=rep(fam_mean_sd_all_var[i, 2], 3), 
#                                              parents_type=rep(fam_mean_sd_all_var[i, 3], 3), 
#                                              si=si, 
#                                              `i*sd_var`=as.numeric(fam_mean_sd_all_var[i, 4:6]), 
#                                              use_var=as.numeric(fam_mean_sd_all_var[i, 7:9])))
# }
# fam_mean_sd_all_var_long <- fam_mean_sd_all_var_long[-1, ]
# fam_mean_sd_all_var_long$ratio <- fam_mean_sd_all_var_long$`i.sd_var` / fam_mean_sd_all_var_long$use_var
# fam_mean_sd_all_var_long$number_of_QTL <- factor(fam_mean_sd_all_var_long$number_of_QTL, 
#                                                  levels=effective_marker_sizes)
# fam_mean_sd_all_var_long$replic <- factor(fam_mean_sd_all_var_long$replic, levels=1:20)
# fam_mean_sd_all_var_long$si <- factor(fam_mean_sd_all_var_long$si, levels=si)
# sum(is.na(fam_mean_sd_all_var_long$ratio))
# # [1] 12
# fam_mean_sd_all_var_long$replic <- as.character(fam_mean_sd_all_var_long$replic)
# 
# 
# 
# pdf(paste("plot_BV_family_mean_sd/", "ratio_var_i*sd_var_use.pdf", sep=""))
# for (i in 1:20){
#   print(
#     ggplot(fam_mean_sd_all_var_long[fam_mean_sd_all_var_long$replic==as.character(i), ]) + 
#       geom_point(aes(si, ratio, color=parents_type)) + 
#       facet_wrap(~number_of_QTL) + 
#       ggtitle(paste("replicate ", i, sep=""))
#   )
# }
# dev.off()



# pdf(paste("plot_BV_family_mean_sd/", "ratio_BV_family_sd_mean.pdf", sep=""))
# ggplot(fam_mean_sd_all, aes(as.numeric(parents_type), ratio)) + 
#   geom_point() + 
#   geom_smooth(method="loess") + 
#   facet_wrap(~number_of_QTL) + 
#   xlab("parents type") +
#   ylab("ratio of BV family sd and mean") +
#   scale_x_continuous(breaks=1:4, labels=as.character(unique(fam_mean_sd_all$parents_type))) + 
#   theme(axis.text.x = element_text(angle = 90)) + 
#   ggtitle("ratio between BV family sd and mean")
# dev.off()
# 
# pdf(paste("plot_BV_family_mean_sd/", "BV_family_mean_sd.pdf", sep=""))
# for (i in 1:length(effective_marker_sizes)){
#   print(ggplot(fam_mean_sd_all[fam_mean_sd_all$number_of_QTL==effective_marker_sizes[i], ], 
#                aes(BV_fammean, BV_famsd)) + 
#           geom_point() + 
#           # geom_smooth(method="loess") + 
#           facet_wrap(~parents_type) +
#           xlab("BV family mean") +
#           ylab("BV family sd") +
#           # scale_x_continuous(breaks=1:4, labels=as.character(unique(fam_mean_sd_all$parents_type))) + 
#           # theme(axis.text.x = element_text(angle = 90)) + 
#           ggtitle(paste("BV family sd against mean, # of causal loci = ", effective_marker_sizes[i], 
#                         sep="")))
# }
# dev.off()

# saveRDS(fam_mean_sd_all, "plot_BV_family_mean_sd/fam_mean_sd_all.rds")

# fam_mean_sd_all[fam_mean_sd_all$ratio < (-400), ]

# fam_mean_sd_all <- readRDS("plot_BV_family_mean_sd/fam_mean_sd_all.rds")
# fam_mean_sd_df <- aggregate(fam_mean_sd_all$BV_fammean, list(fam_mean_sd_all$parents_type, 
#                                                              fam_mean_sd_all$number_of_QTL), var)
# colnames(fam_mean_sd_df) <- c("parents_type", "number_of_QTL", "var_of_mean")
# fam_mean_sd_df$var_of_sd <- 
#   aggregate(fam_mean_sd_all$BV_famsd, list(fam_mean_sd_all$parents_type, 
#                                            fam_mean_sd_all$number_of_QTL), var)[, 3]
# fam_mean_sd_df$ratio <- fam_mean_sd_df$var_of_sd / fam_mean_sd_df$var_of_mean
# 
# pdf(paste("plot_BV_family_mean_sd/", "ratio_var_sd_var_mean.pdf", sep=""))
# ggplot(fam_mean_sd_df, aes(parents_type, ratio)) + 
#   geom_point() + 
#   facet_wrap(~number_of_QTL) + 
#   xlab("parents type") +
#   ylab("ratio of var sd and var mean") + 
#   theme(axis.text.x = element_text(angle = 90)) + 
#   ggtitle("ratio between var sd and var mean")
# dev.off()






