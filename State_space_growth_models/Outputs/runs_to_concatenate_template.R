#### Template for starting a runs_to_concatenate_c*.RData file from scratch
#### Instantiates lists of lists; the format is model[[site]][[pft]]

###____________________________________________________________________________
## Sites using the base model (ddbh_model)
###____________________________________________________________________________

Blodgett = list(
  cedar = c(),
  fir = c(),
  yellow_pine = c(),
  white_pine = c()
)

LaTour = list(
  fir = c(),
  white_pine = c()
)

STEF = list(
  cedar = c(),
  fir = c(),
  yellow_pine = c(),
  white_pine = c()
)

Teakettle = list(
  cedar = c(),
  fir = c(),
  yellow_pine = c(),
  white_pine = c()
)

### Name of list needs to match the model
ddbh_model = list(
  Blodgett = Blodgett,
  LaTour = LaTour,
  STEF = STEF,
  Teakettle = Teakettle
)

###____________________________________________________________________________
## Sites using an informative unit prior (ddbh_model_UnitPrior)
###____________________________________________________________________________

Sequoia = list(
  cedar = c(),
  fir = c(),
  white_pine = c()
)

TharpsCreek = list(
  fir = c(),
  white_pine = c()
)

ddbh_model_UnitPrior = list(
  Sequoia = Sequoia,
  TharpsCreek = TharpsCreek
)

###____________________________________________________________________________
## Sites using an informative tree prior (ddbh_model_TreePrior)
###____________________________________________________________________________

WLakeTahoe = list(
  fir = c()
)

ddbh_model_TreePrior = list(
  WLakeTahoe = WLakeTahoe
)


###____________________________________________________________________________
## Concatenate into a single list of model[[site]][[pft]]
###____________________________________________________________________________

ssm_model = list(
  ddbh_model = ddbh_model,
  ddbh_model_UnitPrior = ddbh_model_UnitPrior,
  ddbh_model_TreePrior = ddbh_model_TreePrior
)

# ### Save
# chain <- 'c1'
# saveRDS(ssm_model, file.path('State_space_growth_models/Outputs', paste0('runs_to_concatenate_', chain, '.RData')))
