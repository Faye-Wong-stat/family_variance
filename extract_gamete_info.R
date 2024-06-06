# extract_gamete_info.R
setwd("~/family_variance/")



effective_marker_sizes <- c(4, 16, 64, 256, 1024)
h2s                    <- c(0.8, 0.5, 0.2)
sf                     <- c(0.8, 0.99)
b                      <- qnorm(sf, 0, 1)
si                     <- dnorm(b, 0, 1) / pnorm(b, 0, 1, lower.tail=F)
indiv_names            <- readRDS("create_marker_list/indiv_names.rds")

Z                        <- readRDS("simulate_phenotypes/Z.rds")
effective_marker_indices <- readRDS("simulate_phenotypes/effective_marker_indices.rds")
results                  <- readRDS("view_correlation_gametes/results.rds")



# get the pairwise mean and sd
results_BV <- results[results$h2s == h2s[1], ]

# this should be arranged in the order of indiv_names
BV_mean_family <- vector(mode="list", length=length(effective_marker_sizes))
BV_mean_family <- lapply(BV_mean_family, function(x){
  vector(mode="list", length=20)
})
BV_sd_family <- vector(mode="list", length=length(effective_marker_sizes))
BV_sd_family <- lapply(BV_sd_family, function(x){
  vector(mode="list", length=20)
})
for (g in 1:length(effective_marker_sizes)){
  for (h in 1:20){
    BV_mean_family[[g]][[h]]           = matrix(NA, nrow=1007, ncol=1007)
    colnames(BV_mean_family[[g]][[h]]) = indiv_names
    rownames(BV_mean_family[[g]][[h]]) = indiv_names
    
    BV_sd_family[[g]][[h]]             = matrix(NA, nrow=1007, ncol=1007)
    colnames(BV_sd_family[[g]][[h]])   = indiv_names
    rownames(BV_sd_family[[g]][[h]])   = indiv_names
    
    print(c(g, h))
    
    A = results_BV[results_BV$effective_marker_sizes==effective_marker_sizes[g] & 
                     results_BV$trait_number==h , ]
    
    BV_mean_family[[g]][[h]] = outer(A$BV_mean, A$BV_mean, "+")/2
    BV_sd_family[[g]][[h]]   = sqrt(outer((A$BV_sd)^2, (A$BV_sd)^2, "+"))
    
    BV_mean_family[[g]][[h]][lower.tri(BV_mean_family[[g]][[h]], diag = T)] = NA
    BV_sd_family[[g]][[h]][lower.tri(BV_sd_family[[g]][[h]], diag = T)] = NA
  }
}

BV_use_family <- vector(mode="list", length=length(si))
BV_use_family <- lapply(BV_use_family, function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(BV_use_family)){
  BV_use_family[[i]] = lapply(BV_use_family[[i]], function(x){
    vector(mode="list", length=20)
  })
}
for (e in 1:length(si)){
  for (g in 1:length(effective_marker_sizes)){
    for (h in 1:20){
      BV_use_family[[e]][[g]][[h]]           = matrix(NA, nrow=1007, ncol=1007)
      colnames(BV_use_family[[e]][[g]][[h]]) = indiv_names
      rownames(BV_use_family[[e]][[g]][[h]]) = indiv_names
      
      BV_use_family[[e]][[g]][[h]] = BV_mean_family[[g]][[h]] + si[e] * BV_sd_family[[g]][[h]]
    }
  }
}

saveRDS(BV_mean_family, "extract_gamete_info/BV_mean_family.R")
saveRDS(BV_sd_family, "extract_gamete_info/BV_sd_family.R")
saveRDS(BV_use_family, "extract_gamete_info/BV_use_family.R")





predY_RR_mean_family <- vector(mode="list", length=length(h2s))
predY_RR_mean_family <- lapply(predY_RR_mean_family, function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(predY_RR_mean_family)){
  predY_RR_mean_family[[i]] = lapply(predY_RR_mean_family[[i]], function(x){
    vector(mode="list", length=20)
  })
}
predY_RR_sd_family <- vector(mode="list", length=length(h2s))
predY_RR_sd_family <- lapply(predY_RR_sd_family, function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(predY_RR_sd_family)){
  predY_RR_sd_family[[i]] = lapply(predY_RR_sd_family[[i]], function(x){
    vector(mode="list", length=20)
  })
}
for (f in 1:length(h2s)){
  for (g in 1:length(effective_marker_sizes)){
    for (h in 1:20){
      predY_RR_mean_family[[f]][[g]][[h]]           = matrix(NA, nrow=1007, ncol=1007)
      colnames(predY_RR_mean_family[[f]][[g]][[h]]) = indiv_names
      rownames(predY_RR_mean_family[[f]][[g]][[h]]) = indiv_names
      
      predY_RR_sd_family[[f]][[g]][[h]]             = matrix(NA, nrow=1007, ncol=1007)
      colnames(predY_RR_sd_family[[f]][[g]][[h]])   = indiv_names
      rownames(predY_RR_sd_family[[f]][[g]][[h]])   = indiv_names
      
      A = results[results$h2s==h2s[f] & 
                     results$effective_marker_sizes==effective_marker_sizes[g] & 
                     results$trait_number==h , ]
      
      predY_RR_mean_family[[f]][[g]][[h]] = outer(A$predY_RR_mean, A$predY_RR_mean, "+")/2
      predY_RR_sd_family[[f]][[g]][[h]]   = sqrt(outer((A$predY_RR_sd)^2, (A$predY_RR_sd)^2, "+"))
      
      predY_RR_mean_family[[f]][[g]][[h]][lower.tri(predY_RR_mean_family[[f]][[g]][[h]], diag = T)] = NA
      predY_RR_sd_family[[f]][[g]][[h]][lower.tri(predY_RR_sd_family[[f]][[g]][[h]], diag = T)] = NA
    }
  }
  
}

predY_RR_use_family <- vector(mode="list", length=length(si))
predY_RR_use_family <- lapply(predY_RR_use_family, function(x){
  vector(mode="list", length=length(h2s))
})
for (i in 1:length(predY_RR_use_family)){
  predY_RR_use_family[[i]] = lapply(predY_RR_use_family[[i]], function(x){
    vector(mode="list", length=length(effective_marker_sizes))
  })
  for (j in 1:length(predY_RR_use_family[[i]])){
    predY_RR_use_family[[i]][[j]] = lapply(predY_RR_use_family[[i]][[j]], function(x){
      vector(mode="list", length=20)
    })
  }
}
for (e in 1:length(si)){
  for (f in 1:length(h2s)){
    for (g in 1:length(effective_marker_sizes)){
      for (h in 1:20){
        predY_RR_use_family[[e]][[f]][[g]][[h]]           = matrix(NA, nrow=1007, ncol=1007)
        colnames(predY_RR_use_family[[e]][[f]][[g]][[h]]) = indiv_names
        rownames(predY_RR_use_family[[e]][[f]][[g]][[h]]) = indiv_names
        
        predY_RR_use_family[[e]][[f]][[g]][[h]]           = 
          predY_RR_mean_family[[f]][[g]][[h]] + si[e] * predY_RR_sd_family[[f]][[g]][[h]]
      }
    }
  }
}

saveRDS(predY_RR_mean_family, "extract_gamete_info/predY_RR_mean_family.R")
saveRDS(predY_RR_sd_family, "extract_gamete_info/predY_RR_sd_family.R")
saveRDS(predY_RR_use_family, "extract_gamete_info/predY_RR_use_family.R")



predY_mean_family <- vector(mode="list", length=length(h2s))
predY_mean_family <- lapply(predY_mean_family, function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(predY_mean_family)){
  predY_mean_family[[i]] = lapply(predY_mean_family[[i]], function(x){
    vector(mode="list", length=20)
  })
}
predY_sd_family <- vector(mode="list", length=length(h2s))
predY_sd_family <- lapply(predY_sd_family, function(x){
  vector(mode="list", length=length(effective_marker_sizes))
})
for (i in 1:length(predY_sd_family)){
  predY_sd_family[[i]] = lapply(predY_sd_family[[i]], function(x){
    vector(mode="list", length=20)
  })
}
for (f in 1:length(h2s)){
  for (g in 1:length(effective_marker_sizes)){
    for (h in 1:20){
      predY_mean_family[[f]][[g]][[h]]           = matrix(NA, nrow=1007, ncol=1007)
      colnames(predY_mean_family[[f]][[g]][[h]]) = indiv_names
      rownames(predY_mean_family[[f]][[g]][[h]]) = indiv_names
      
      predY_sd_family[[f]][[g]][[h]]             = matrix(NA, nrow=1007, ncol=1007)
      colnames(predY_sd_family[[f]][[g]][[h]])   = indiv_names
      rownames(predY_sd_family[[f]][[g]][[h]])   = indiv_names
      
      A = results[results$h2s==h2s[f] & 
                    results$effective_marker_sizes==effective_marker_sizes[g] & 
                    results$trait_number==h , ]
      
      predY_mean_family[[f]][[g]][[h]] = outer(A$predY_mean, A$predY_mean, "+")/2
      predY_sd_family[[f]][[g]][[h]]   = sqrt(outer((A$predY_sd)^2, (A$predY_sd)^2, "+"))
      
      predY_mean_family[[f]][[g]][[h]][lower.tri(predY_mean_family[[f]][[g]][[h]], diag = T)] = NA
      predY_sd_family[[f]][[g]][[h]][lower.tri(predY_sd_family[[f]][[g]][[h]], diag = T)] = NA
    }
  }
  
}

predY_use_family <- vector(mode="list", length=length(si))
predY_use_family <- lapply(predY_use_family, function(x){
  vector(mode="list", length=length(h2s))
})
for (i in 1:length(predY_use_family)){
  predY_use_family[[i]] = lapply(predY_use_family[[i]], function(x){
    vector(mode="list", length=length(effective_marker_sizes))
  })
  for (j in 1:length(predY_use_family[[i]])){
    predY_use_family[[i]][[j]] = lapply(predY_use_family[[i]][[j]], function(x){
      vector(mode="list", length=20)
    })
  }
}
for (e in 1:length(si)){
  for (f in 1:length(h2s)){
    for (g in 1:length(effective_marker_sizes)){
      for (h in 1:20){
        predY_use_family[[e]][[f]][[g]][[h]]           = matrix(NA, nrow=1007, ncol=1007)
        colnames(predY_use_family[[e]][[f]][[g]][[h]]) = indiv_names
        rownames(predY_use_family[[e]][[f]][[g]][[h]]) = indiv_names
        
        predY_use_family[[e]][[f]][[g]][[h]]           = 
          predY_mean_family[[f]][[g]][[h]] + si[e] * predY_sd_family[[f]][[g]][[h]]
      }
    }
  }
}

saveRDS(predY_mean_family, "extract_gamete_info/predY_mean_family.R")
saveRDS(predY_sd_family, "extract_gamete_info/predY_sd_family.R")
saveRDS(predY_use_family, "extract_gamete_info/predY_use_family.R")
