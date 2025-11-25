### State-space model for tree growth (Equations 1-4 in main text)
### Modified to set an informative prior for the standard deviation on the 
### unit (plot) random effect

library(nimble)

# nimble code
stateSpaceCode <- nimbleCode({
  for(treenum in 1:n_trees) {
    for(yr in (start_year[treenum] + 1):end_year[treenum]){
      
      # Mean log DBH increment (deterministic)
      log_ddbh_predicted[treenum,yr] <- intercept[treenum, yr] + 
        beta1[treenum, yr] * (log(size[treenum, yr-1]) - obs_ln_mean)/obs_ln_stdv + # center and scale log(dbh)
        beta2[treenum, yr] * ((size[treenum, yr-1])^2 - obs_sq_mean)/obs_sq_sd + # center and scale log(dbh)
        betaCWD * cwd[treenum, yr -1] + 
        treat_mean[treatment[plotnum[treenum]]] +
        treat_cwd[treatment[plotnum[treenum]]] * cwd[treenum, yr -1] +
        treat_size[treatment[plotnum[treenum]]] * (log(size[treenum, yr-1]) - obs_ln_mean)/obs_ln_stdv
      
      # True DBH increment follows a lognormal distribution
      ddbh[treenum, yr] ~ dlnorm(log_ddbh_predicted[treenum, yr], sdlog = res_sd)
      
      # True size in yr is size in previous year plus true log DBH increment
      size[treenum,yr] <- size[treenum, yr-1] + ddbh[treenum,yr]
    }
    
    # Observed log DBH increment follows a normal distribution
    for(obs in 1:nobs[treenum]){
      observation[treenum, obs] ~ dnorm(size[treenum, year[treenum, obs]], 
                                        sd = obs_sd)
    }
    
    # Prior on initial size
    size[treenum, start_year[treenum]] ~  dnorm(obs_mean, 1e-6) 
  }
  
  for(treenum in 1:n_trees){
    # Tree random effect
    int_tree[treenum] ~ dnorm(0, sd = int_tree_sd)
    
    for(yr in 1:n_years){
      # Overall intercept
      intercept[treenum, yr] <- int_overall + int_tree[treenum] + 
        int_plot[plotnum[treenum]] 
      
      # beta1 and beta2 are constant over time and tree
      beta1[treenum, yr] <- log_size
      beta2[treenum, yr] <- d2_size
    }
  }
  
  # Unit (plot) random effect
  for (j in 1:n_plots) {
    int_plot[j] ~ dnorm(0, sd = int_plot_sd)
  }
  
  # Priors
  int_overall ~ dnorm(0, 1.0E-6)
  log_size ~ dnorm(0, 1.0E-6)
  d2_size ~ dnorm(0, 1.0E-6)
  betaCWD ~ dnorm(0, 1.0E-6)
  int_tree_sd ~ dunif(0, 100)
  int_plot_sd ~ dnorm(0.19, (1/0.03)^2) # informative prior based on cedar and fir at other sites
  res_sd ~ dunif(0, 100)
  obs_sd ~ dunif(0.0289, 100)
  
  treat_mean[1] <-0 # set to zero, because this is reference
  treat_cwd[1] <- 0 # set to zero, because this is reference
  treat_size[1] <- 0 # set to zero, because this is reference
  for(ind in 2:n_treats){
    treat_mean[ind] ~ dnorm(0, 1.0E-6)
    treat_cwd[ind] ~ dnorm(0, 1.0E-6)
    treat_size[ind] ~ dnorm(0, 1.0E-6)
  }
}

)

# Observations
data = list(observation = obs_pft)

# Constants
constants = list(n_trees = nrow(treeData),
                 n_plots = nrow(plotData),
                 n_years = max(treeData$end_year),
                 n_treats = length(unique(plotData$fTreatment)),
                 start_year = treeData$start_year,
                 end_year = treeData$end_year,
                 plotnum = as.numeric(treeData$fPlot),
                 nobs = treeData$n_obs,
                 year = year_obs,
                 obs_mean = obs_mean,
                 obs_sq_mean = obs_sq_mean,
                 obs_sq_sd = obs_sq_sd,
                 obs_ln_mean = obs_ln_mean,
                 obs_ln_stdv = obs_ln_stdv,
                 cwd = cwd_scaled,
                 treatment = as.numeric(plotData$fTreatment))

# Initial values
inits = list(int_overall = rnorm(1, 0, sd=5), 
             res_sd = runif(1, 0.00001, 1), 
             obs_sd= runif(1, 0.00001, 2), 
             int_tree_sd = runif(1, 0.00001, 1), 
             int_plot_sd = rnorm(1, 0.19, 0.03), # Corresponds to informative prior
             int_tree = rnorm(nrow(treeData),0, sd=1), 
             int_plot = rnorm(nrow(plotData),0, sd=1), 
             log_size = rnorm(1, 0, sd=1), 
             d2_size = rnorm(1, 0, sd=1), 
             betaCWD = rnorm(1, 0, sd=1),
             treat_mean = c(NA, rnorm(3, 0, 1)), 
             treat_cwd = c(NA, rnorm(3, 0, 1)),
             treat_size = c(NA, rnorm(3, 0, 1)),
             size = as.matrix(matching_size_init),
             ddbh = as.matrix(matching_ddbh_init)
)

# Record the tree and unit random effect for each MCMC iteration
monitors2 = c("int_tree", "int_plot") 
