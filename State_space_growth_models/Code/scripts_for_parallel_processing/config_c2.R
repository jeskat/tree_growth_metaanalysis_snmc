###_____________________________________________________________________________
### CONFIG FILE: CHAIN 2
###_____________________________________________________________________________

### Use this file to provide thinning intervals and number of MCMC iterations 
### for each site- and PFT-specific state-space model run. Update for each
### run in a chain.
### If it's the first run in the chain, specify first_run=TRUE; otherwise
### first_run = FALSE. 

### Blodgett
Blodgett <- list(
  cedar = list(thin_interval = 500,
               first_run = FALSE,
               n_iter = 2000000), 
  
  yellow_pine = list(thin_interval = 500,
                     first_run = FALSE,
                     n_iter = 10000000), 
  
  fir = list(thin_interval = 500, 
             first_run = FALSE,
             n_iter = 3000000), 
  
  white_pine = list(thin_interval = 500,
                    first_run = FALSE,
                    n_iter = 10000000) 
)

### LaTour
LaTour <- list(
  fir = list(thin_interval = 500,
             first_run = TRUE,
             n_iter = 5000000),
  
  white_pine = list(thin_interval = 500,
                    first_run = TRUE,
                    n_iter = 10000000)
)

### Sequoia
Sequoia <- list(
  cedar = list(thin_interval = 500,
               first_run = TRUE,
               n_iter = 5500000),
  
  fir = list(thin_interval = 500,
             first_run = TRUE,
             n_iter = 2500000),
  
  white_pine =list(thin_interval = 500,
                   first_run = TRUE,
                   n_iter = 10000000)
  
)

### STEF
STEF <- list(
  cedar = list(thin_interval = 500,
               first_run = TRUE,
               n_iter = 5000000),
  
  yellow_pine = list(thin_interval = 500,
                     first_run = TRUE,
                     n_iter = 27000000),
  
  fir = list(thin_interval = 500, 
             first_run = TRUE,
             n_iter = 6000000),
  
  white_pine = list(thin_interval = 500,
                    first_run = TRUE,
                    n_iter = 10000000)
  
)


### Teakettle
Teakettle <- list(
  cedar = list(thin_interval = 100,
               first_run = TRUE,
               n_iter = 600000), 
  
  yellow_pine = list(thin_interval = 100,
                     first_run = TRUE,
                     n_iter = 3000000), 
  
  fir = list(thin_interval = 100,
             first_run = TRUE,
             n_iter = 200000),
  
  white_pine = list(thin_interval = 500,
                    first_run = TRUE,
                    n_iter = 1000000)
)


### Tharp's Creek
TharpsCreek <- list(
  fir = list(thin_interval = 500,
             first_run = TRUE,
             n_iter = 700000),
  
  white_pine = list(thin_interval = 500,
                    first_run = TRUE,
                    n_iter = 10000000)
)


### West Lake Tahoe
WLakeTahoe <- list(
  fir = list(thin_interval = 500,
             first_run = TRUE,
             n_iter = 2500000)
)


configs <- list(
  Blodgett = Blodgett,
  LaTour = LaTour,
  Sequoia = Sequoia,
  STEF = STEF,
  Teakettle = Teakettle,
  TharpsCreek = TharpsCreek,
  WLakeTahoe = WLakeTahoe
)
