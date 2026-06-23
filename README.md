# Meta-analysis of tree growth responses to fuel reduction treatments across Sierra Nevada mixed-conifer forests

This repository contains code used to generate the results, figures, and tables in the article "Tree growth varies with fuel treatment type, tree size, and climate in Sierra Nevada mixed conifer forests" by Jessica R. Katz, Perry de Valpine, Alice Cima, John J. Battles, Adrian J. Das, Eric E. Knapp, Kathryn E. Low, Patricia N. Manley, Malcolm P. North, Polly C. Thornton, Angela M. White, and Lara M. Kueppers.

For any questions, please contact Jessica Katz (jessica_katz@berkeley.edu).

## Citation
[Placeholder]

## Overview
In this study, we compiled and analyzed data from seven field experiments and continuous forest inventories to synthesize tree growth responses to mechanical thinning and prescribed burning across the geographically- and climatically-diverse Sierra Nevada mixed-conifer (SNMC) forests of California. The full analysis consists of multiple components, each represented by a folder in this repository: 1) processing forest inventory data; 2) processing covariate and meta-regressor data; 3) estimating 19 state-space growth models (SSGMs), each model specific to a unique combination of plant functional type [PFT] and site; 4) using mixed-effects meta-analysis and meta-regression to estimate pooled effects across sites; and 5) producing figures and tables for the publication. 

The seven sites included in this analysis include:

* “Blodgett:” The Fire and Fire Surrogate Study located at Blodgett Forest Research Station (Stephens and Moghaddas 2005);
* "LaTour:" The Continuous Forest Inventory program at LaTour Demonstration State Forest (Ritchie and Withrow, Jr. 1964);
* "Sequoia:" The Fire and Fire Surrogate Study located at Sequoia National Park (Knapp et al. 2005);
* "STEF:" The Variable Density Thinning Experiment at Stanislaus-Tuolumne Experimental Forest (Knapp et al. 2017); 
* "Teakettle:" The Teakettle Ecosystem Experiment (North et al. 2002); 
* "Tharp's Creek:" A prescribed fire study located in the Tharp’s Creek and Log’s Creek watersheds of Sequoia National Park (Mutch and Parsons 1998); and 
* “W. Lake Tahoe:” A thinning study located in West Lake Tahoe (Low et al. 2021).

The figure below illustrates the main components of this repository and the general workflow. For a detailed description of the study and methods, see Katz et al. ([placeholder]). Due to the large size of the Markov chain Monte Carlo outputs of the SSGMs, we do not include these outputs in the repository. However, they can be reproduced by running the code within the "State_space_growth_models" folder.   


![Overview of analysis workflow. Folders are represented by boxes; those containing data files (no code) have a tab in the lower right corner. The purple dashed box denotes the contents of this repository; items outside of the dashed box are external to this repository.](code_workflow.png)

### Forest inventory data availability
All site-level forest inventory data files needed to run the SSGMs are saved in the `State_space_growth_models/Input_data` directory. We processed these data files from raw forest inventory data for each field site. Raw inventory data are publicly available for Blodgett at [Stephens et al. 2025](https://doi.org/10.6073/pasta/06a3a2556f5289b89b9f8683becfe08e) and for Sequoia and Tharp's Creek at [Das 2026](https://doi.org/10.5066/P1MHZKJY). For the four other sites, the raw data are not available publicly but can be requested from the authors. As an example of our data pre-processing procedure, we provide the code needed to convert the Blodgett, Sequoia, and Tharp's Creek inventory data to SSGM inputs in the `Inventory_data_processing/` directory of this repository. Data pre-processing for the other sites followed a very similar procedure, with minor changes to reflect the unique structure of each raw dataset. 

## Contents
### Inventory_data_processing
Contains code for processing raw inventory data to create the inputs needed to run the SSGMs. Key processes include: data cleaning (removing duplicated, null, and filled records; harmonizing tree species and life status across each tree’s remeasurements; and correcting typos); removing trees with less than two post-treatment diameter at breast height (DBH) measurements; removing trees with DBH < 15cm; removing remeasurements on dead trees; removing remeasurements associated with growth outliers; and formatting data tables for SSGM ingestion. 

We provide processing scripts only for sites with permanently archived inventory data (i.e., Blodgett, Sequoia, and Tharp's Creek).


* `Code/`
  + `config.py`: provides values and list-style mappings to enable consistent data processing (e.g., consistent treatment naming conventions and minimum diameter thresholds) across all sites. 
  + `inventory_processing_functions.py`: provides all common functions needed to process inventory data for each site. These functions are called in `process_inventory_data_[site].ipynb` and `get_standStructure_[site].ipynb` notebooks.
  + `process_inventory_data_Blodgett.ipynb`: pipeline for processing raw inventory data [(Stephens et al. 2025)](https://doi.org/10.6073/pasta/06a3a2556f5289b89b9f8683becfe08e) and producing SSGM input data for Blodgett. Saves processed data to `State_space_growth_models/Input_data/Blodgett/`.
  + `process_inventory_data_Sequoia_Tharp.ipynb`: pipeline for processing raw inventory data [Das 2026](https://doi.org/10.5066/P1MHZKJY) and producing SSGM input data for Sequoia and Tharp's Creek. Notebook must be run separately for each site, specifying "Sequoia" or "TharpsCreek" in cell 2. Saves processed data to `State_space_growth_models/Input_data/Sequoia/` or `State_space_growth_models/Input_data/TharpsCreek/`. 
  
### Covariate_and_metaregressor_processing
Contains code for processing the climatic water deficit (CWD) covariate in the SSGMs, along with the mean climate, soil, and stand structure variables that are used in the meta-regressions. 

* `Code/`
  + `covariate_functions.py`: contains functions for associating the geospatial coordinates of treatment units with raw covariate and/or meta-regressor data. 
  + `process_terraclimate_data.R`: takes in raw data from Terraclimate (which must be downloaded separately); shifts calendar year to water year; and aggregates each climate variable (winter precipitation, CWD, Palmer Drought Severity Index [PDSI]) from monthly to annual. Outputs the `*_by_waterYr.nc` netCDF files saved in the `Covariate_and_metaregressor_processing/Processed_data/Climate/` subfolder.
  + `get_CWD_timeseries_for_SSM.ipynb`: takes in the output of `process_terraclimate_data.R` and extracts the CWD timeseries associated with each tree within each site. Produces the `cwd_tree_year.csv` file (one of the inputs to the SSGM). Works for Blodgett, Sequoia, and Tharp's Creek; other sites have this file saved in the permanent archive. 
  + `create_climate_metaregressors.ipynb`: takes in the output of `process_terraclimate_data.R` and calculates the winter precipitation, CWD, and PDSI metaregressor for each treatment unit within each site. Outputs `site_mean_clim.csv`.
  + `process_ssurgo_variables.R`: adapts the workflow of [Chasen 2020](https://rpubs.com/emchasen/SSURGOcleaning) to import and clean soil survey (SSURGO) data for each survey area in which a study site falls. Outputs `all_unit_soils.csv`.
  + `create_soil_metaregressors.ipynb`: takes in the output of `process_ssurgo_variables.R` and other SSURGO data (which must be downloaded and saved) and calculates the available water capacity (AWC) and soil depth meta-regressors for each treatment unit within each site. Outputs `unit_mean_soil.csv`. 
  + `get_standStructure_Blodgett.ipynb`: creates relative stand density index (rSDI) and treatment intensity metaregressors along with pre-treatment species composition for each treatment unit at Blodgett, using raw inventory data [(Stephens et al. 2025)](https://doi.org/10.6073/pasta/06a3a2556f5289b89b9f8683becfe08e) as input. Creates `Blodgett_intensity.csv`, `Blodgett_rsdi.csv`, and `Blodgett_pretrt_comp.csv`. Other sites use an analogous process to produce similar outputs.
  + `get_standStructure_Sequoia_Tharp.ipynb`: creates relative stand density index (rSDI) and treatment intensity metaregressors along with pre-treatment species composition for each treatment unit at Sequoia or Tharp's Creek (user must specify in cell 2), using raw inventory data [(Das 2026)](https://doi.org/10.5066/P1MHZKJY) as input. Creates `Sequoia_intensity.csv`, `Sequoia_rsdi.csv`, and `Sequoia_pretrt_comp.csv` or `TharpsCreek_intensity.csv`, `TharpsCreek_rsdi.csv`, and `TharpsCreek_pretrt_comp.csv`.
  + `compile_all_metaregressors.ipynb`: takes in the outputs of `create_climate_metaregressors.ipynb`, `create_soil_metaregressors.ipynb`, and site-specific stand structure processing scripts (e.g., `get_standStructure_Blodget.ipynb`) and compiles all meta-regressors for all sites into one table. Outputs `Covariate_and_metaregressor_processing/Processed_data/all_regressors_by_unit.csv`. 
* `Processed_data/`:
  + `all_regressors_by_unit.csv`: Tabular output of `compile_all_metaregressors.ipynb`, containing unit-specific climate, soil, and stand structure meta-regressors for all sites. 
  + `Climate/`
    - `*_by_waterYr.nc`: NetCDF outputs of `process_terraclimate_data.R`, providing gridded values of winter precipitation ("ppt"), CWD, and PDSI for each water year.
    - `site_mean_clim.csv`: tabular output of `create_climate_metaregressors.ipynb`, containing estimates of each climate variable for all sites.
  + `Soil/`:
    - `all_unit_soils.csv`: tabular output of `process_ssurgo_variables.R`, containing estimates of soil depth and AWC for each SSURGO map unit.
    - `unit_mean_soil.csv`: tabular output of `create_soil_metaregressors.ipynb`, containing treatment unit estimates of soil depth and AWC for all sites. 
  + `Stand_structure/`
    - `*_intensity.csv`: percent reduction in basal area between the first post-treatment and last pre-treatment measurement for each treatment unit within the site (the site is identified in the filename).
    - `*_rdsi.csv`: estimated SDI, site-specific maximum SDI, and rSDI for each treatment unit within the site (the site is identified in the filename). 
    - `*_pretrt_comp.csv`: species fraction of pre-treatment basal area (trees >=15cm) in each treatment unit within the site (the site is identified in the filename). 
  
## State_space_growth_models
Contains the code needed to run site- and PFT-specific SSGMs (in parallel if needed) and run diagnostics for the resulting Markov chain Monte Carlo (MCMC) outputs. 

* `Code/`
  + `scripts_for_parallel_processing`: contains the files needed to run multiple SSGMs in parallel using a slurm workload manager.
    - `parallel_ss_cores.sh`: shell script that dispatches the parallel SSGM model runs. When called by `sbatch`, it will run a SSGM for each combination of site, PFT, and model given in the three `.lst` files within this sub-folder.
    - `args_ssm.sh`: shell script that calls `State_space_growth_models/Code/run_ssm.R` with the arguments passed from `parallel_ss_cores.sh` and renames the log file. New users should replace the parameters (job name, account, partition, etc.) with those appropriate for their slurm computing environment. 
    - `pft.lst`: list of PFTs (e.g., cedar, fir, yellow_pine, white_pine) for which to run SSGMs. 
    - `site.lst`: list of sites (e.g., Blodgett, Latour, STEF) for which to run SSGMs. Site names must be among the sub-folder names in `State_space_growth_models/Input_data`. 
    - `model.list`: list of models (e.g., ddbh_model, corresponding to one or more file names within `State_space_growth_models/Code/nimble_models`) for which to run SSGMs. 
    - `config_c1.R` and `config_c2.R`: R scripts with lists that configure each SSGM for the first ("c1") and second ("c2") MCMC chains. Users may manually edit the MCMC thinning interval and number of iterations for each site- and PFT-specific SSGM. 
  + `nimble_models`: contains different SSGM formulations using the nimble package.
    - `ddbh_model.R`: default state-space model for tree growth (corresponds to Equations 1-4 in Katz et al. (placeholder)).
    - `ddbh_model_TreePrior.R`: modification of `ddbh_model.R` to set an informative prior for the standard deviation on the tree random effect. Used only for the W. Lake Tahoe SSGMs.
    - `ddbh_model_UnitPrior.R`: modification of `ddbh_model.R` to set an informative prior for the standard deviation on the unit (plot) random effect. Used only for the Sequoia and Tharp's Creek SSGMs.
  + `ssm_config.R`: defines file paths for inputs and outputs of SSGMs. New users should specify the directory in which the MCMC outputs of the SSGMs should be saved.
  + `process_inputs.R`: takes in data from the `State_space_growth_models/Input_data` subfolders and creates the data and constant variables needed to initialize the nimble model for a given site and PFT. 
  + `ssm_restart_functions.R`: functions written by [Daniel Turek](https://danielturek.github.io/public/saveMCMCstate/saveMCMCstate.html) that allow a nimble model's internal states to be saved. These saved states can be reloaded into the nimble model to restart it where it left off (which needs to be done for very long chains during which one or more R sessions are terminated before the MCMC has converged).
  + `run_ssm.R`: takes in arguments passed by the shell scripts in `scripts_for_parallel_processing/` and either starts or resumes the production of an MCMC for an SSGM for a given site and PFT. The MCMCs are configured using ``scripts_for_parallel_processing/config_*.R`. Calls `create_ssm_diagnostics.R` after the run has completed to produce a diagnostic file. 
  + `ssm_functions.R`: provides a set of functions needed to process the MCMC outputs of an SSGM.
  + `create_ssm_diagnostics.R`: creates an html document with figures and tables useful for observing and diagnosing the progress of the MCMC outputs produced by an SSGM. Summarizes details across multiple restarts of the same chain; creates traceplots; summarizes parameter estimates; plots observed and modeled growth for a set of randomly selected trees; and calculates the effective sample size for each parameter. 
  + `batch_ssm_diagnostics.R`: loops through all SSGMs in the list of completed models (e.g., `State_space_growth_model_Outputs/authors_complete_models.RData`), applies the specified burn-in for each SSGM, and creates and saves final diagnostic figures and parameter estimates for each SSGM in  `State_space_growth_models/Outputs/Final_SSM_diagnostics/`.
  + `convergence_diagnostics_one_model.R`: gets Gelman-Rubin diagnostics for two chains of the same SSGM.
  + `batch_convergence_diagnostics.R`: loops through all SSGMs in the list of completed models (e.g., `State_space_growth_model_Outputs/authors_complete_models.RData`) and compiles Gelman-Rubin statistics for all parameters. Output saved in `State_space_growth_models/Outputs/Final_SSM_diagnostics/all_gr_diagnostics.csv`. 
  + `make_convergence_diagnostic_table.R`: reads the csv produced by `batch_convergence_diagnostics.R` and creates a summary table of Gelman-Rubin statistics (Table S5). 
  + `compile_ssm_params.R`: loops through all  SSGMs in the list of completed models (e.g., `State_space_growth_model_Outputs/authors_complete_models.RData`) and compiles the final parameter estimates into one dataframe, saved as `State_space_growth_model_Outputs/allParams_[date].csv`.
  + `summarize_parameter_estimates`: reads the csv produced by `compile_ssm_params.R` and creates Tables S6-S7 and Figure S8, summarizing parameter estimates from each SSGM.
* `Outputs/`: 
  + `allParams_2025-12-02.csv`: tabular output of `compile_ssm_params.R`, summarizing parameter estimates for each site- and PFT-specific SSGM.
  + `authors_complete_models.RData`: R data object containing the final list of models (specified by site, PFT, and nimble model name), along with burn-in values. Users should create their own version of this object for new analyses.
  + `runs_to_concatenate_c*.RData`: R data object providing a structure in which to store the lists of MCMC filenames associated with each SSGM (this is important if the MCMC needs to be restarted one or more times). The filename indicates the chain (c1=chain 1; c2 = chain 2) associated with the SSGM. These files are automatically updated when SSGMs are run via `run_ssm.R`. 
  + `runs_to_concatenate_template.R`: R script providing a template to create a new, blank version of `runs_to_concatenate_c*.RData` at the beginning of a new project.
  + `Final_SSM_diagnostics/`:
    - `all_gr_diagnostics.csv`: tabular output of `batch_convergence_diagnostics.R`, containing the Gelman-Rubin statistics for each parameter of each site- and PFT-specific SSGM. 
    - `[site]_[pft]_[chain].html`: each html file (produced by `create_ssm_diagnostics.R`) contains MCMC diagnostic figures and tables for a SSGM for a given combination of site and PFT (specified in the filename). The chain specified in the filename corresponds to the chain from which results in Katz et al. (placeholder) were reported. 
  
## Metaanalysis_models
Contains the code needed to post-process outputs of the SSGMs and estimate all meta-analysis and meta-regression models.

* `Code/`
  + `metaanalysis_config.R`: defines study-wide significance thresholds, list of meta-regressors, and other inputs to meta-analyses and meta-regressions.
  + `metaanalysis_functions.R`: suite of functions that are used to run meta-analyses and meta-regressions.
  + `postprocess_fixed_effects.R`: loops through all  SSGMs in the list of completed models (e.g., `State_space_growth_model_Outputs/authors_complete_models.RData`) and reparameterizes the posterior distributions for DBH = 30 cm and CWD = 0. Outputs unit-level estimates of growth (`Metaanalysis_models/Reparameterized_SSM_outputs/unit_growth_effects_by_trt_*.csv`); site-level parameter estimates for the log(DBH) and CWD effects under each treatment (`Metaanalysis_models/Reparameterized_SSM_outputs/fixed_effects_by_trt_*.csv`); and the variance-covariance matrices for each of these outputs (*.RData files). 
  + `postprocess_unit_growth_CWD_DBH.R`: loops through all  SSGMs in the list of completed models (e.g., `State_space_growth_model_Outputs/authors_complete_models.RData`), a list of DBH values, and a list of CWD values, and reparameterizes the SSGM posterior distributions for each combination of DBH and CWD. Outputs `Metaanalysis_models/Reparameterized_SSM_outputs/growthOutcomesByCWDandDBH*`, including growth outcomes by treatment (.csv file) and variance-covariance matrices (.RData file).
  + `get_pooled_effects.R`: takes in the outputs of `postprocess_fixed_effects.R`. Runs meta-analyses and makes Figure 2 and Supplementary Tables S8-S10 (pooled effects of treatment, tree size, measurement period CWD, and their interactions on tree growth). 
  + `calculate_trt_effects_by_CWD_DBH.R`: takes in the outputs of `postprocess_unit_growth_CWD_DBH.R` and runs meta-analyses of treatment effects on growth for different combinations of PFT, DBH, and CWD. Makes Figures 3, 4, and 5.
  + `run_metaregressions.R`: runs meta-regressions for mean climate, soil, and stand structure variables. Makes Table 2, Figures S9-S15, and Tables S11-21. 
* `Reparameterized_SSM_outputs/`:
  + `unit_growth_effects_by_trt_2025-12-09.csv`: tabular output of `postprocess_fixed_effects.R`, containing estimates of growth in each treatment unit based on SSGM posterior distributions re-parameterized for DBH = 30 and CWD = 0.
  + `cov_unit_growth_by_trt_2025-12-09.RData`: R data object output of `postprocess_fixed_effects.R`, containing variance-covariance matrices for growth in each treatment unit based on SSGM posterior distributions re-parameterized for DBH = 30 and CWD = 0.
  + `fixed_effects_by_trt_2025_12_09.csv`: tabular output of `postprocess_fixed_effects.R`, containing estimates of growth and the log(DBH) and CWD effects by treatment for each site.
  + `cov_cwd_by_trt_2025-12-09.RData`: R data object output of `postprocess_fixed_effects.R`, containing variance-covariance matrices for the CWD effects by treatment for each site.
  + `cov_size_by_trt_2025-12-09.RData`:  R data object output of `postprocess_fixed_effects.R`, containing variance-covariance matrices for the tree size (log DBH) effects by treatment for each site.
  + `growthOutcomesByCWDAndDBH.csv`: tabular output of `postprocess_unit_growth_CWD_DBH.R`, containing estimates of growth in each treatment unit based on SSGM posterior distributions re-parameterized for each combination of DBH and CWD.
  + `growthOutcomesByCWDAndDBH.RData`: R data object output of `postprocess_unit_growth_CWD_DBH.R`, containing variance-covariance matrices for growth in each treatment unit based on SSGM posterior distributions re-parameterized for each combination of DBH and CWD.

## Figures_and_Tables
Contains figures and tables included in the manuscript, and code to make supplementary figures and tables not produced elsewhere in this repository.

* `Code/`
  + `plotting_config.R`: specifies values used in creating tables and graphics, including colors, styles, ordering, and naming conventions.
  + `plotting_functions.R`: suite of functions used specifically for making figures and tables.
  + `make_site_map.R`: creates Figure 1 (map of study sites).
  + `make_metaregressor_boxplots.R`: creates Figures S1-S3 (scatterplots and boxplots of meta-regressors).
  + `make_composition_barchart.R`: creates Figure S4 (pre-treatment species composition by site).
  + `characterize_PFTs_across_sites.R`: creates Table S1 (number of trees and observations per PFT and site) and Figure S5 (distribution of post-treatment PFT sizes at each site).
  + `make_metaregressor_corrplot.R`: creates Figure S7 (correlation plot of meta-regressors).
* `Figures/`: contains jpeg objects, with filenames corresponding to Figure numbers in Katz et al. (placeholder).
* `Tables/`: contains csv files, with filenames corresponding to Table numbers in Katz et al. (placeholder).
  

## References
Das, A.J. 2026. Sequoia National Park tree growth data for thinning and burning study: U.S. Geological Survey data release. https://doi.org/10.5066/P1MHZKJY.

Knapp, E. E., J. E. Keeley, E. A. Ballenger, and T. J. Brennan. 2005. “Fuel Reduction and Coarse Woody Debris Dynamics with Early Season and Late Season Prescribed Fire in a Sierra Nevada Mixed Conifer Forest.” Forest Ecology and Management 208 (1): 383–97. https://doi.org/10.1016/j.foreco.2005.01.016.

Knapp, E. E., J. M. Lydersen, M. P. North, and B. M. Collins. 2017. “Efficacy of Variable Density Thinning and Prescribed Fire for Restoring Forest Heterogeneity to Mixed-Conifer Forest in the Central Sierra Nevada, CA.” Forest Ecology and Management 406: 228–41. https://doi.org/10.1016/j.foreco.2017.08.028.

Low, K. E., B. M. Collins, A. A. Bernal, et al. 2021. “Longer-Term Impacts of Fuel Reduction Treatments on Forest Structure, Fuels, and Drought Resistance in the Lake Tahoe Basin.” Forest Ecology and Management 479: 118609. https://doi.org/10.1016/j.foreco.2020.118609.

Mutch, L. S., and D. J. Parsons. 1998. “Mixed Conifer Forest Mortality and Establishment Before and After Prescribed Fire in Sequoia National Park, California.” Forest Science 44 (3): 341–55. https://doi.org/10.1093/forestscience/44.3.341.

North, M. P., B. Oakley, J. Chen, et al. 2002. Vegetation and Ecological Characteristics of Mixed-Conifer and Red Fir Forests at the Teakettle Experimental Forest. U.S. Department of Agriculture, Forest Service, Pacific Southwest Research Station. https://www.fs.fed.us/psw/publications/documents/psw_gtr186/.

Ritchie, R. W., and R. N. Withrow, Jr. 1964. Latour State Forest Continuous Forest Inventory Plan: Field Instructions. State of California, The Resources Agency, Department of Conservation, Division of Forestry.

Stephens, S. L., and J. J. Moghaddas. 2005. “Experimental Fuel Treatment Impacts on Forest Structure, Potential Fire Behavior, and Predicted Tree Mortality in a California Mixed Conifer Forest.” Forest Ecology and Management 215 (1): 21–36. https://doi.org/10.1016/j.foreco.2005.03.070.

Stephens, S. L., R. A. York, A. Roughton, J. J. Battles, and B. M. Collins. 2025. “The Fire and Fire Surrogate Study: Berkeley Forests 2001-2020.” Edi. 2104.1. Version 1. Environmental Data Initiative, August 12. https://doi.org/10.6073/pasta/06a3a2556f5289b89b9f8683becfe08e.
