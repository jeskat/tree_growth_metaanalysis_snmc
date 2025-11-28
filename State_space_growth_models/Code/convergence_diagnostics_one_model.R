## Gets Gelman-Rubin diagnostics for two chains of the same model.

library('coda')
library(readr)


## Uncomment and specify variables below to run for only one model. 
# site <- 'Blodgett'
# pft <- 'white_pine'
# model <- 'ddbh_model'
# burn_in <- 250

# Get correct name of site to pull input data

source(here::here('State_space_growth_models/Code/process_inputs.R'))
source(here::here('State_space_growth_models/Code/ssm_functions.R'))

chain1 <- load_mcmc(site = site, 
                    pft = pft, 
                    model = model, 
                    burn_in = burn_in, 
                    chain = 'c1',
                    mcmc_results_dir = mcmc_results_dir,
                    to_report = to_report
                    )
chain2 <- load_mcmc(site = site, 
                    pft = pft, 
                    model = model, 
                    burn_in = burn_in, 
                    chain = 'c2',
                    mcmc_results_dir = mcmc_results_dir,
                    to_report = to_report
                    )

# Make chains the same length by trimming the longer one
if(dim(chain1)[1]>dim(chain2)[1]){
  chain1 <- as.mcmc(chain1[(dim(chain1)[1]-dim(chain2)[1]+1):dim(chain1)[1],])
}else if(dim(chain2)[1]>dim(chain1)[1]){
  chain2 <- as.mcmc(chain2[(dim(chain2)[1]-dim(chain1)[1]+1):dim(chain2)[1],])
}

# Remove columns that only contain zeros
to_keep <- colnames(chain1[, colSums(chain1 != 0, na.rm=TRUE) > 0])

## Get Gelman-Rubin statistics
combinedChains <- mcmc.list(chain1[, to_keep], chain2[, to_keep])
g <- gelman.diag(combinedChains, autoburnin = TRUE)
