## ---- setup, include=FALSE------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
library(readr)
library(stringr)
library('coda')
library(here)
source(here::here('State_space_growth_models/Code/ssm_config.R'))

## ---- inputs, echo=FALSE--------------------------------------------------------------------------------


## Uncomment to run for a specified site. 
## Otherwise, these inputs will be passed in from state_space_model_template.
# site <- 'Blodgett'
# pft <- 'fir'
# model <- 'ddbh_model'
# chain <- 'c1'


# Load input data for this site and PFT
source(file.path(here::here('State_space_growth_models/Code', 'ssm_functions.R')))
source(file.path(here::here('State_space_growth_models/Code', 'process_inputs.R')))


# Set burn-in
if(exists('burn_in')){
  burn_in <- burn_in
}else{
  burn_in <- 1
}

###____________________________________________________________________________
## Generate a set of random trees for which to visualize the SSM modeled sizes
###____________________________________________________________________________

## Empty list; can replace with list of specific tree tags if desired.
rndm_tags <- c() 

## Define max size thresholds (cm) for small and medium trees
sm_size_max <- 30 # smallest tree at UPFU is 15 cm
med_size_max <- 60

set.seed(1)

## List of trees and sizes
treeList <- obs_pft

## Randomly select three tree tags in each size class
if(is.null(rndm_tags)){
  # Make size class labels
  sm_class <- paste0('<',sm_size_max,' cm')
  md_class <- paste0(sm_size_max,'-',med_size_max,' cm')
  lg_class <- paste0('>=',med_size_max,' cm')
  
  # Add a column defining the size class associated with each tree at its first observation
  treeList$Size_class = NA
  # Small
  treeList[treeList[,1] <sm_size_max, 'Size_class'] = sm_class
  # Medium
  treeList[(treeList[,1]>=sm_size_max) & (treeList[,1]<med_size_max), 'Size_class'] = md_class
  # Large
  treeList[treeList[,1]>=med_size_max, 'Size_class'] = lg_class
  
  for(i in c(sm_class,md_class,lg_class)){
    rndm_tags <- append(rndm_tags, sample_trees(treeList,i))
  }
  
}

print(rndm_tags)

## Make a list of column names of all latent states (DBHs) associated with each random tree
end_yr <- max(treeData$end_year)
col_nm <- list()

for (tree_nm in rndm_tags) {
  tree_no <- which(rownames(treeData)==tree_nm)
  tree_str <- gsub('x', as.character(tree_no), 'size[x, y]')

  for (i in 1:end_yr) {col_nm <- append(col_nm, gsub('y', as.character(i), tree_str))
  }
}

## Add these columns to the reported parameters
to_report_rndmTrees <- c(to_report, unlist(col_nm))

## ---- load_mcmc, echo=FALSE, message=FALSE, warning=FALSE----------------------------------------------

###____________________________________________________________________________
## Load MCMC samples from all runs
###____________________________________________________________________________


## Load and combine into one data table all MCMC iterations associated with this chain
all_samples <- load_mcmc(site = site,
                         pft = pft, 
                         model = model, 
                         burn_in = burn_in,
                         chain = chain, 
                         mcmc_results_dir = mcmc_results_dir, 
                         to_report = to_report_rndmTrees)


## ---- summary, echo=FALSE-------------------------------------------------------------------------------

###____________________________________________________________________________
## Generate summary information about MCMC chain
###____________________________________________________________________________

## Read in the list of prior runs in this chain.
if(chain == 'c1'){
  ssm_model <- readRDS(here::here('State_space_growth_models/Outputs', runs_c1_fn))
  }else if(chain == 'c2'){
    ssm_model <- readRDS(here::here('State_space_growth_models/Outputs', runs_c2_fn))
  }

## List, in chronological order, of the output sub-directories that contain the 
## MCMC iterations for this site and PFT.
outdirs <- ssm_model[[model]][[site]][[pft]] 

print(paste0("Dimensions of all_samples: ", dim(all_samples)))
print(paste0('Site: ', site))
print(paste0('PFT: ', pft))
print(paste0('Model: ', model))
print(paste0('No. Trees: ', nrow(treeData)))
print(paste0('No. Plots:', nrow(plotData)))
print(paste0('Treatments: ', unique(plotData$treatment)))
print(paste0('Start date: ',   unlist(strsplit(outdirs[1], '_'))[1]))
print(paste0('End date: ', unlist(strsplit(outdirs[length(outdirs)], '_'))[1]))
print(paste0('No. segments: ', length(outdirs)))
print(paste0('Burn-in: ', burn_in))

## For each run, make a list of the time elapsed, thin interval, and number of iterations
elapsed <- vector(mode="numeric", length=length(outdirs))
thin_itv <- vector(mode="numeric", length=length(outdirs))
n_iters <- vector(mode="numeric", length=length(outdirs))

for(i in 1:length(outdirs)){
  substrs <- unlist(strsplit(outdirs[i], '_'))
  fn <- paste0(pft, '_', site, '_', model, '_', substrs[length(substrs)], '.Rout')
  summary_file <- readLines(here::here(mcmc_results_dir, outdirs[i], fn))
  print(paste0('Segment ', i))

  elapsed[i] <- as.numeric(unlist(strsplit(summary_file[grep("sec elapsed", summary_file)], ' '))[1])

  thn_ln <- unlist(strsplit(summary_file[grep("\"thin_interval: [[:digit:]]", summary_file)], ':'))[2]
  thin_itv[i] <- as.numeric(unlist(strsplit(thn_ln, "\\\"")))

  iter_ln <- unlist(strsplit(summary_file[grep("\"MCMC iterations: [[:digit:]]", summary_file)], ': '))[2]
  n_iters[i] <- as.numeric(unlist(strsplit(iter_ln, "\\\"")))
}

## Report metrics across all runs
print(paste0('Total number of MCMC iterations: ', sum(n_iters)))
print(paste0('Thin interval: ', thin_itv))
print(paste0('Total time elapsed: ', round(sum(elapsed)/3600, 2), ' hours'))

# Print WAIC from final segment
waic_ln <- grep("Field \"WAIC\"", summary_file)
print(summary_file[waic_ln:(waic_ln + 5)])


## ---- traceplots, echo=FALSE----------------------------------------------------------------------------

###____________________________________________________________________________
## Trace plots
###____________________________________________________________________________

# Remove random trees from parameters to report (but keep what we added in load_mcmc)
to_report <- colnames(all_samples)[!(colnames(all_samples) %in% unlist(col_nm))]

for(v in to_report){
  if(is.na(all_samples[1, v])){
    print(paste("Null value in ", v))
  }else{
    plot(all_samples[, v])
    title(v)
  }
}


## ---- paramest, echo=FALSE------------------------------------------------------------------------------

###____________________________________________________________________________
## Create summary statistics for posterior distributions
###____________________________________________________________________________


## Loop through and summarize all parameters of interest
means <- c()
lower_ci <- c()
upper_ci <- c()
sds <- c()
ses <- c()
tsse <- c()
p<- c()

for(i in 1:(length(to_report))){
  if(is.na(all_samples[1, to_report[i]])){
    print(paste('Null value in', to_report[i]))
  }else{
    p <- append(p, to_report[i])
    var <- summarize_mcmc(all_samples, to_report[i])
    means <- append(means,as.numeric(var['Mean']))
    lower_ci <- append(lower_ci, as.numeric(var['2.5%']))
    upper_ci <- append(upper_ci, as.numeric(var['97.5%']))
    sds <- append(sds, as.numeric(var['SD']))
    ses <- append(ses, as.numeric(var['Naive SE']))
    tsse <- append(tsse, as.numeric(var['Time-series SE']))
  }
}


## Format outputs as a dataframe
out_df <- data.frame(variable = p,mean = means, low_ci = lower_ci, high_ci = upper_ci, sd = sds,
                     naive_se = ses, timeseries_se = tsse)


out_df

## Save outputs dataframe 
write.csv(out_df, 
          here::here(mcmc_results_dir, 
                     ### CHANGE
                     outdirs[length(outdirs)], # 'Parameter_summaries',
                           stri_replace_all_regex('params_Z_Y_X.csv',
                                                  pattern = c('Z','Y','X'), 
                                                  replacement = c(as.character(site), 
                                                                  as.character(pft),
                                                                  as.character(model)),
                                                  vectorize=FALSE)), row.names = FALSE)

# ---- ts, echo=FALSE------------------------------------------------------------------------------------

###____________________________________________________________________________
## Plot observed and modeled DBH vs. time for random trees
###____________________________________________________________________________

## Identify size and deltaDBH columns associated with example trees
ex_trees_size <- c()
ex_trees_ddbh <- c()
for (i in rndm_tags) {
 tree_no <- which(rownames(treeData)==i)
 yr_rng <- (treeData[tree_no, 'start_year']+1):(treeData[tree_no, 'end_year'])
 smp_yr <- sample(yr_rng,1)
 ex_trees_size<- append(ex_trees_size, stri_replace_all_regex('size[x, y]',pattern = c('x','y'),
                                                    replacement = c(as.character(tree_no),
                                                                    as.character(smp_yr)),
                                                    vectorize=FALSE))
 ex_trees_ddbh<- append(ex_trees_ddbh, stri_replace_all_regex('ddbh[x, y]',pattern = c('x','y'),
                                                              replacement = c(as.character(tree_no),
                                                                              as.character(smp_yr)),
                                                              vectorize=FALSE))
}

# Graph individual tree timeseries
end_yr <- max(treeData$end_year)

for (tree_nm in rndm_tags) {
  tree_no <- which(rownames(treeData)==tree_nm)
  tree_str <- gsub('x', as.character(tree_no), 'size[x, y]')
  col_nm <- list()
  for (i in 1:end_yr) {col_nm <- append(col_nm, gsub('y', as.character(i), tree_str))
  }

  annual_size_mean <- c()
  lower_ci <- c()
  upper_ci <- c()

  for (i in 1:end_yr) {
    if(is.na(all_samples[1, as.character(col_nm[i])])){
      print(paste('Null value in', col_nm[i]))
      all_samples[1, as.character(col_nm[i])] <- 0
    }
    out_yr <- summary(all_samples[,as.character(col_nm[i])])
    annual_size_mean[i] = out_yr$statistics['Mean']
    lower_ci[i] = out_yr$quantiles['2.5%']
    upper_ci[i] = out_yr$quantiles['97.5%']
  }
  plot_df <- data.frame(mean = annual_size_mean, low_ci = lower_ci, high_ci = upper_ci)
  plot_df[plot_df==0]=NA
  plot_df <- plot_df 

  ## Plot observations
  plot(as.numeric(year_obs[tree_no,]), as.numeric(obs_pft[tree_no,]), 
       ylim = c(min((min(plot_df$low_ci, na.rm = TRUE)), 
                    min(as.numeric(obs_pft[tree_no,]), na.rm = TRUE)),
                max(max(plot_df$high_ci, na.rm = TRUE), 
                    max(as.numeric(obs_pft[tree_no,]), na.rm = TRUE))),
       xlab = "Observation year", ylab = "Diameter (cm)")


  ## Plot modeled DBH (mean and 95% credible interval)
  lines(seq(1, end_yr, by=1), plot_df$mean, col='grey', lwd=1)
  lines(seq(1, end_yr, by=1), plot_df$low_ci, col='grey', lwd=1, lty="dashed")
  lines(seq(1, end_yr, by=1), plot_df$high_ci, col='grey', lwd=1, lty="dashed")
  title(tree_no)

  legend("topleft", legend=c("Observation", "Mean", "95% CI"),
       lty=c(NA, "solid", "dashed"), pch = c('o', NA, NA), cex=0.8)
}



## ---- ess, echo=FALSE-----------------------------------------------------------------------------------

###____________________________________________________________________________
## Calculate effective sample size (ESS) for each parameter
###____________________________________________________________________________

## Calculate ESS using coda functions
ESS <- effectiveSize(all_samples[, to_report])

## Order from lowest to highest ESS
ESS[order(ESS)]


## ------------------------------------------------------------------------------------------------------
# # Option to move output file
# outdir <- file.path('/global','scratch','users','jessicakatz','forest_growth_synthesis','Outputs','SS_outputs', outdirs[length(outdirs)])
# file.copy(from = 'show_outputs.html', to = file.path(outdir, paste0(site, '_', pft, '_', model, '.html')), overwrite = TRUE)


