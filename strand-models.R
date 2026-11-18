# ##############################################################################################################################
# #
# #   MULTIPLEX Binomial Analyses
# #
# ##############################################################################################################################
#  library(devtools)
#  install_github('ctross/STRAND@ctross-patch-1') # Need to install a patch that supports numpyro multiplex models with masking

##############################################################################################################################
# Load libraries
 library(STRAND)
 library(stringr)
 library(ggplot2)
 library(psych)

##############################################################################################################################

## Impute NAs. This is just a robustness check
impute_sample = function(x) {
  na_idx = is.na(x)

  if(any(na_idx)){
    x[na_idx] = sample(
      x[!na_idx],
      size = sum(na_idx),
      replace = TRUE
    )
  }

  return(x)
}

## This rescales beta variable off of end points, but we dont actually use it here
ScaleBeta = function(X){
 a = 0
 b = 100
 y = X
 Samp = 50

 y2 = (y-a)/(b-a)
 y3 = (y2*(Samp - 1) + 0.5)/Samp
 return(y3)
}

## This will build a STRAND data object for the f_cd model
Make_STRAND_Data_From_RDS = function(d){
 dist = d$allocations$dist
 dist[is.infinite(dist)] = max(dist[which(!is.infinite(dist))],na.rm=TRUE) + 1

 model_dat = list(
    N = nrow(d$allocations),

    N_ID = max(as.numeric(factor(d$allocations$ID))),
    N_AID = max(as.numeric(factor(d$allocations$AID))),
    N_RID = max(as.numeric(factor(d$allocations$RID))),
    N_ID_SUID = max(as.numeric(factor(d$allocations$ID_SUID))),
    N_AID_SUID = max(as.numeric(factor(d$allocations$AID_SUID))),

    percGiven = ScaleBeta(d$allocations$percGiven),
    # OutcomeBin = round(d$allocations$percGiven/10, 0),
    OutcomeBin = round(d$allocations$percGiven/0.1, 0),
    Condition = as.numeric(factor(d$allocations$Condition)), # 1 anon, 2 revealed

    ones = rep(1, nrow(d$allocations)),
    anon_first = ifelse(d$allocations$anon_first, 1,0),
    dist = standardize_strand(dist),
    ego_cent = standardize_strand(d$allocations$ego_cent),
    alter_cent = standardize_strand(d$allocations$alter_cent),

    ID = as.numeric(factor(d$allocations$ID)),
    AID = as.numeric(factor(d$allocations$AID)),
    RID = as.numeric(factor(d$allocations$RID)),
    ID_SUID = as.numeric(factor(d$allocations$ID_SUID)),
    AID_SUID = as.numeric(factor(d$allocations$AID_SUID))
  )

 for(q in 7:14){
  model_dat[[q]] = impute_sample(model_dat[[q]])  # Randomly fills NAs via re-sampling. Not the best thing to do, but just to draft a model.
 }

 # Now reshape to STRAND style
 IDs = unique(c(d$allocations$ID, d$allocations$AID))
 NN = length(IDs)
 AF = Dist = Revealed = Anon = matrix(-99, nrow = NN, ncol = NN)
 Exposure = matrix(10, nrow = NN, ncol = NN)

 colnames(AF) = colnames(Dist) = colnames(Revealed) = colnames(Anon) = IDs
 rownames(AF) = rownames(Dist) = rownames(Revealed) = rownames(Anon) = IDs

 for(n in 1:model_dat$N){
    if(model_dat$Condition[n] == 2){
      Revealed[which(IDs == d$allocations$ID[n]), which(IDs == d$allocations$AID[n])] = model_dat$OutcomeBin[n]
    }
    if(model_dat$Condition[n] == 1){
      Anon[which(IDs == d$allocations$ID[n]), which(IDs == d$allocations$AID[n])] = model_dat$OutcomeBin[n]
    }

    AF[which(IDs == d$allocations$ID[n]), which(IDs == d$allocations$AID[n])] = model_dat$anon_first[n]
    Dist[which(IDs == d$allocations$ID[n]), which(IDs == d$allocations$AID[n])] = model_dat$dist[n]
  }

 ind = data.frame(ID=IDs, Cent=NA)
 rownames(ind) = IDs
 for(n in 1:model_dat$N){
   ind$Cent[which(ind$ID == d$allocations$ID[n])] = model_dat$ego_cent[n]
   ind$Cent[which(ind$ID == d$allocations$AID[n])] = model_dat$alter_cent[n]
 }

 # Outcomes stored as a labeled list
 outcome = list(
  Revealed = Revealed,
  Anon = Anon
 )

 # Exposure stored as a labeled list
 exposure = list(
  Revealed = Exposure,
  Anon = Exposure
 )
 colnames(exposure$Revealed) = colnames(exposure$Anon) = IDs
 rownames(exposure$Revealed) = rownames(exposure$Anon) = IDs

 # The mask lets us cut out NAs, now recorded as -99, from the adjacency matrices
 mask = list(
  Revealed = ifelse((Revealed == -99), 1, 0),
  Anon = ifelse((Anon == -99), 1, 0)
 )

 colnames(mask$Revealed) = colnames(mask$Anon) = IDs
 rownames(mask$Revealed) = rownames(mask$Anon) = IDs

 dyadic = list(
  AF = AF,
  Dist = Dist
 )

 # Merge data
 dat = make_strand_data(
  outcome = outcome,
  block_covariates = NULL,
  individual_covariates = ind,
  dyadic_covariates = dyadic,
  exposure = exposure,
  mask = mask,
  outcome_mode="binomial",
  link_mode="logit",
  multiplex = TRUE,
  check_standardization = FALSE
 )

 return(dat)
}

# Load data

allocs <- readRDS("allocs.rds")

N_sites = length(unique(allocs$site))
d_list = res = fit = vector("list", N_sites)
site_names = unique(allocs$site)

d_list <- lapply(site_names, \(s) list(allocations = allocs[allocs$site == s,]))

# And we can fit the models
for(s in 1:N_sites){
 fit[[s]] = fit_multiplex_model(
 data = Make_STRAND_Data_From_RDS(d_list[[s]]),
 block_regression = ~ 1,
 focal_regression = ~ Cent,
 target_regression = ~ Cent,
 dyad_regression = ~ AF + Dist,
 mode="numpyro",
 mcmc_parameters = list(
   chains = 1,
   parallel_chains = 1,
   refresh = 1,
   iter_warmup = 1000,
   iter_sampling = 1000,
   max_treedepth = 11,
   adapt_delta = 0.96)
)
}

for(s in 1:N_sites){
  temp = summarize_strand_results(fit[[s]])

  temp2 = rbind(
    sum_stats("Revelation",temp$sample$srm_model_samples$block_parameters$Revealed[[1]][,1,1] - temp$sample$srm_model_samples$block_parameters$Anon[[1]][,1,1],0.9),
    sum_stats("Revelation_Focal_Cent",temp$sample$srm_model_samples$focal_coeffs[,1,1] - temp$sample$srm_model_samples$focal_coeffs[,2,1],0.9),
    sum_stats("Revelation_Target_Cent",temp$sample$srm_model_samples$target_coeffs[,1,1] - temp$sample$srm_model_samples$target_coeffs[,2,1],0.9),
    sum_stats("Revelation_Dyad_AF",temp$sample$srm_model_samples$dyadic_coeffs[,1,1] - temp$sample$srm_model_samples$dyadic_coeffs[,2,1],0.9),
    sum_stats("Revelation_Dyad_Dist",temp$sample$srm_model_samples$dyadic_coeffs[,1,2] - temp$sample$srm_model_samples$dyadic_coeffs[,2,2],0.9)
  )

  temp3 = rbind(
    sum_stats("Base",temp$sample$srm_model_samples$block_parameters$Anon[[1]][,1,1],0.9),
    sum_stats("Base_Focal_Cent",temp$sample$srm_model_samples$focal_coeffs[,2,1],0.9),
    sum_stats("Base_Target_Cent",temp$sample$srm_model_samples$target_coeffs[,2,1],0.9),
    sum_stats("Base_Dyad_AF",temp$sample$srm_model_samples$dyadic_coeffs[,2,1],0.9),
    sum_stats("Base_Dyad_Dist",temp$sample$srm_model_samples$dyadic_coeffs[,2,2],0.9)
  )
  colnames(temp2) = colnames(temp3) = colnames(temp$summary)
  res[[s]] = rbind(temp$summary, temp2, temp3)
  res[[s]]$Site = site_names[s]
}

# Plot results

res_all = do.call(rbind, res)

 colors = c("#0D205C", "#507BAF", "#C5570E", "#DF913E", "#13667B", "#35ABC0",
"#E39913", "#F3BB57", "#721712", "#AF5249", "#0B231A", "#4E7A65",
"#431008", "#A56137")

in_set = c(
"Generalized reciprocity - Revealed (out) - Anon (out)",
"Generalized reciprocity - Revealed (in) - Anon (in)",
"Dyadic reciprocity - Revealed (i to j) - Anon (i to j)",
"Dyadic reciprocity - Revealed (i to j) - Revealed (j to i)",
"Dyadic reciprocity - Anon (i to j) - Anon (j to i)",
"Revelation",
"Revelation_Focal_Cent",
"Revelation_Target_Cent",
"Revelation_Dyad_Dist",
"Base_Focal_Cent",
"Base_Target_Cent",
"Base_Dyad_Dist")

d_sub = d[which(d$Variable %in% in_set),]

d_sub$Variable2 = d_sub$Variable

d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Generalized reciprocity - Revealed (out) - Anon (out)"), "Gen. Recip. (Rev. Out - Anon. Out)")
d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Generalized reciprocity - Revealed (in) - Anon (in)"), "Gen. Recip. (Rev. In - Anon. In)")

d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Dyadic reciprocity - Revealed (i to j) - Anon (i to j)"), "Dyad. Recip. (Rev. i-to-j - Anon. i-to-j)")
d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Dyadic reciprocity - Revealed (i to j) - Revealed (j to i)"), "Dyad. Recip. (Rev. i-to-j - Rev. j-to-i)")
d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Dyadic reciprocity - Anon (i to j) - Anon (j to i)"), "Dyad. Recip. (Anon. i-to-j - Anon. j-to-i)")

d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Base_Dyad_Dist"), "Distance")
d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Base_Focal_Cent"), "Donor centrality")
d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Base_Target_Cent"), "Recipient centrality")

d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Revelation_Dyad_Dist"), "Revealed X Distance")
d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Revelation_Focal_Cent"), "Revealed X Donor centrality")
d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Revelation_Target_Cent"), "Revealed X Recipient centrality")

d_sub$Variable2 = str_replace_all(d_sub$Variable2, fixed("Revelation"), "Revealed")

d_sub$L = as.numeric(d_sub$"HPDI:0.05")
d_sub$H = as.numeric(d_sub$"HPDI:0.95")
d_sub$M = as.numeric(d_sub$Median)


d_sub$Layer = ifelse(str_detect(d_sub$Variable2,"Dyad."), "Dyadic Recip.",
                 ifelse(str_detect(d_sub$Variable2,"Gen."), "Generalized Recip.",
                 "Base Coeffs."
                  ))

d_sub$Layer = factor(d_sub$Layer)
d_sub$Layer = factor(d_sub$Layer, levels = c("Base Coeffs.", "Generalized Recip.", "Dyadic Recip."))

d_sub$Variable3 = factor(d_sub$Variable2)
d_sub$Variable3 = factor(d_sub$Variable3, levels = rev(c(
        "Revealed", "Distance", "Donor centrality", "Recipient centrality",
        "Revealed X Distance", "Revealed X Donor centrality", "Revealed X Recipient centrality",
        "Dyad. Recip. (Rev. i-to-j - Anon. i-to-j)",
        "Dyad. Recip. (Anon. i-to-j - Anon. j-to-i)",
        "Dyad. Recip. (Rev. i-to-j - Rev. j-to-i)",
        "Gen. Recip. (Rev. In - Anon. In)",
        "Gen. Recip. (Rev. Out - Anon. Out)"
          )))

p = ggplot(d_sub, aes(x = Variable3, y = M, group = Site, color=Site,
        ymin = L, ymax = H)) + geom_linerange(size = 1,, position = position_dodge(width = 0.6)) +
        geom_point(size = 2,, position = position_dodge(width = 0.6)) + facet_grid(Layer ~ . , scales = "free", space="free_y") +
         geom_hline(aes(yintercept = 0),
        color = "black", linetype = "dashed") + labs(y = "Regression parameters",
        x = "") + theme(strip.text.x = element_text(size = 12,
        face = "bold"), strip.text.y = element_text(size = 12,
        face = "bold"), axis.text = element_text(size = 12),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.title.x = element_blank()) +
        guides(color = guide_legend(nrow = 1)) +
        theme(strip.text.y = element_text(angle = 360)) +
        coord_flip() + theme(panel.spacing = grid::unit(1,
        "lines")) + theme(legend.position="bottom") + scale_color_manual(values=c("az" = colors[8],
                                                                                  "np" = colors[4],
                                                                                  "te" = colors[3],
                                                                                  "tm" = colors[2],
                                                                                  "tz" = colors[1],
                                                                                  "up" = colors[6],
                                                                                  "wh" = colors[7],
                                                                                  "wl" = colors[5]
                                                                                  ))


ggsave("STRAND_Multiplex_Robustness_Check.pdf", p, height=8, width=18)
