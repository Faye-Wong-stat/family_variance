# view_usefulness_best_parents_gametes.R

setwd("~/family_variance/")
library(ggplot2)
library(cowplot)



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s                    <- c(0.8, 0.5, 0.2)
sf                     <- c(0.8, 0.99)
b                      <- qnorm(sf, 0, 1)
si                     <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)
types <- c("best predicted mean from all cross", "best predicted usefulness from all cross", 
           "best predicted mean from elite parents", "best predicted usefulness from elite parents")

best_pred_mean  <- readRDS("select_best_parents_gametes/best_pred_mean.rds")
best_pred_use   <- readRDS("select_best_parents_gametes/best_pred_use.rds")
best_pred_mean2 <- readRDS("select_best_parents_gametes/best_pred_mean2.rds")
best_pred_use2  <- readRDS("select_best_parents_gametes/best_pred_use2.rds")

best_pred_mean$type  <- "best predicted mean from all cross"
best_pred_use$type   <- "best predicted usefulness from all cross"
best_pred_mean2$type <- "best predicted mean from elite parents"
best_pred_use2$type  <- "best predicted usefulness from elite parents"

best_pred_long <- rbind(best_pred_mean, best_pred_use, best_pred_mean2, best_pred_use2)
best_pred_long$effective_marker_sizes <- factor(best_pred_long$effective_marker_sizes, 
                                                levels=as.factor(effective_marker_sizes))
best_pred_long$h2 <- paste("h^2 == ", best_pred_long$h2s, sep="")
best_pred_long$i  <- paste("i == ", round(best_pred_long$si, 2), sep="")
best_pred_long[, 6:9] <- sapply(best_pred_long[, 6:9], as.numeric)
# sapply(best_pred_long, function(x){sum(is.na(x))})

# 2*3*5*20*4 = 2400
usefulness_mean <- data.frame(si=rep(si, each=3*5*20*4), 
                   h2s=rep(rep(h2s, each=5*20*4), 2), 
                   effective_marker_sizes=rep(rep(effective_marker_sizes, each=20*4), 2*3), 
                   trait_number=rep(rep(1:20, each=4), 2*3*5), 
                   type=rep(types, 2*3*5*20),
                   mean_BV_use=NA)
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (k in 1:20){
        for (l in 1:length(types)){
          A = best_pred_long[best_pred_long$si==si[h] & 
                               best_pred_long$h2s==h2s[i] & 
                               best_pred_long$effective_marker_sizes==effective_marker_sizes[j] & 
                               best_pred_long$trait_number==k & 
                               best_pred_long$type==types[l], ]
          usefulness_mean[(h-1)*3*5*20*4 + (i-1)*5*20*4 + (j-1)*20*4 + (k-1)*4 + l, "mean_BV_use"] = 
            mean(A$BV_use)
        }
      }
    }
  }
}

# 2*3*5*4 = 120
usefulness_mean_se <- data.frame(si=rep(si, each=3*5*4), 
                              h2s=rep(rep(h2s, each=5*4), 2), 
                              effective_marker_sizes=rep(rep(effective_marker_sizes, each=4), 2*3), 
                              type=rep(types, 2*3*5),
                              mean_BV_use_mean=NA, 
                              mean_BV_use_se=NA)
for (h in 1:length(si)){
  for (i in 1:length(h2s)){
    for (j in 1:length(effective_marker_sizes)){
      for (l in 1:length(types)){
        A = usefulness_mean[usefulness_mean$si==si[h] & 
                              usefulness_mean$h2s==h2s[i] & 
                              usefulness_mean$effective_marker_sizes==effective_marker_sizes[j] & 
                              usefulness_mean$type==types[l], ]
        usefulness_mean_se[(h-1)*3*5*4 + (i-1)*5*4 + (j-1)*4  + l, 
                           c("mean_BV_use_mean", "mean_BV_use_se")] = 
          c(mean(A$mean_BV_use), sd(A$mean_BV_use)/sqrt(20))
      }
    }
  }
}
usefulness_mean_se$effective_marker_sizes <- factor(usefulness_mean_se$effective_marker_sizes, 
                                                    levels=as.factor(effective_marker_sizes))
usefulness_mean_se$h2 <- paste("h^2 == ", usefulness_mean_se$h2s, sep="")
usefulness_mean_se$i  <- paste("i == ", round(usefulness_mean_se$si, 2), sep="")



p1 <- ggplot(usefulness_mean_se, aes(as.numeric(effective_marker_sizes))) + 
  geom_point(aes(y=mean_BV_use_mean, color=type)) + 
  geom_errorbar(aes(ymin=mean_BV_use_mean-mean_BV_use_se, 
                    ymax=mean_BV_use_mean+mean_BV_use_se, 
                    color=type), width=0.2) + 
  geom_line(aes(y=mean_BV_use_mean, color=type), linewidth=0.5, alpha = 0.5) +
  facet_grid(i~h2, labeller = label_parsed) +
  xlab("number of causal loci") + 
  ylab("true usefulness") + 
  # coord_cartesian(ylim=c(30, 50)) + 
  scale_x_continuous(breaks=1:5, labels=as.character(effective_marker_sizes)) + 
  guides(color=guide_legend(title="selection method", ncol=1)) + 
  theme_minimal_grid(font_size=10) +
  theme(legend.position="bottom") 
save_plot(paste("view_usefulness_best_parents_gametes/plots/", "true_use_under_4selection_type.pdf", sep=""),
          plot_grid(p1),
          base_width=6.5, base_height=4.33)



for (i in 1:length(types)){
  print(types[i])
  print(summary(usefulness_mean_se[usefulness_mean_se$type==types[i], ]$mean_BV_use_mean))
  print(summary(usefulness_mean_se[usefulness_mean_se$type==types[i], ]$mean_BV_use_se))
}

for (i in 1:length(h2s)){
  for (j in 1:length(effective_marker_sizes)){
    print(
      usefulness_mean_se[usefulness_mean_se$h2s==h2s[i] & 
                           usefulness_mean_se$effective_marker_sizes==effective_marker_sizes[[j]], ]
    )
  }
}







