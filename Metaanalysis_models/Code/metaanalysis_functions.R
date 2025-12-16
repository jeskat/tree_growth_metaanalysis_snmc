
###-----------------------------------------------------------------------------
### Functions for formatting inputs to the meta-analyses
###-----------------------------------------------------------------------------

### Function takes in dataframe df, treatment trt (Burn or Thin), and the name of 
### the column (var_col) that designates treatment. Returns the same dataframe 
### with a "Burn" or "Thin" column in which the values are Thin/no_Thin and 
### Burn/no_Burn, respectively. 
addBurnThinCols <- function(df, trt, var_col){
  df[[trt]] <- paste0('no_', trt)
  df[grepl(trt, df[[var_col]], ignore.case = TRUE), trt] = trt
  df[[trt]] <- factor(df[[trt]])
  df[[trt]] <- relevel(df[[trt]], ref = paste0('no_', trt))
  return(df)
}

### Function that takes a dataframe df (with `site` and `pft` columns) and 
### returns same dataframe with a column with site and PFT aliases. 
addPrettyNames <- function(df){
  # Add pretty names for plotting
  df$common_nm <- NaN
  df$PFT <- NaN
  
  for(i in 1:length(df$site)){
    df$common_nm[i] <- common_nms[df$site[i]]
    df$PFT[i] <- pft_dict[[df$pft[i]]]
  }
  
  # Make Common name and PFT factors
  df$common_nm <- factor(df$common_nm, site_order)
  df$PFT <- factor(df$PFT, level = c('Cedar','Fir','Yellow Pine', 'Sugar Pine'))
  df$PFT <- relevel(df$PFT, ref = 'Fir')
  
  return(df)
}


### Function that takes in a dataframe df, adds "Burn" and "Thin" columns, and
### adds site and PFT aliases. Dataframes need to be formatted in this way for
### subsequent analyses in metafor.
formatForMetafor <- function(df){
  df <- addBurnThinCols(df = df, trt = 'Burn', var_col = 'variable')
  df <- addBurnThinCols(df = df, trt = 'Thin', var_col = 'variable')
  
  # Treatment is a factor. None is reference level.
  ref_level <- unique(df$variable)[grepl('None', unique(df$variable))]
  df$variable <- factor(df$variable)
  df$variable <- relevel(df$variable, ref = ref_level)
  
  df <- addPrettyNames(df)
  
  # Site is also a factor
  df$site <- factor(df$site)
  
  return(df)
}

### Function that takes in a variance-covariance object, turns it into a matrix,
### and makes the rows and column names match. Variance-covariance matrices
### need to be in this format for subsequent analyses in metafor.
get_cv_w_names <- function(m){
  new_mat <- as.matrix(m[, 2:length(m)])
  rownames(new_mat) = colnames(new_mat) = m$X
  return(new_mat)
}

###-----------------------------------------------------------------------------
### Meta-analysis functions
###-----------------------------------------------------------------------------

## Function to run a metafor model for one pft
mf_meta_one_pft <- function(full_data, # dataframe formatted according to `formatForMetafor`
                            pft_to_plot, # modeled pft (raw, not alias for plotting)
                            cv_list, # list of var-cov matrices associated with full_data
                            formula, # meta-analysis or meta-regression formula
                            by_unit, # Meta-analysis by unit? T/F. 
                            mtrgrsn, # Meta-regression? T/F
                            regressor, # If meta-regression, specify metaregressor column
                            rgr_df, # If meta-regression, provide metaregressor dataframe
                            optimizer, # metafor optimizer
                            p_dict = pft_dict, # Mapping of raw pft names to aliases
                            clubSandwich = TRUE # Apply cluster-robust inference methods (robust variance estimation)
){
  
  # Filter PFT of interest
  to_plot <- full_data[full_data$pft==pft_to_plot,]
  cv_pft <- cv_list[grepl(pft_to_plot, names(cv_list))]
  
  if(by_unit){
    ## Prepare data for metaregression
    if(mtrgrsn){
      ## Add metaregressor data to the dateframe
      to_plot <- left_join(to_plot, 
                           rgr_df, 
                           join_by(site == Site, 
                                   unit == UnitID,
                                   Treatment == Treatment))
      
      ## Center regressor
      to_plot[[regressor]] = scale(to_plot[[regressor]], center = TRUE, scale = FALSE)
    }
    
    ## Metafor setup: 
    mfmod <- rma.mv(yi = mean, V = cv_pft, mods = formula, 
                    # random = list(~variable|site, ~unit|site),
                    # struct = c('CS', 'ID'),
                    
                    ## Alternative formulation
                    random = list(~1|site, ~variable|site, ~unit|site),
                    struct = c('ID', 'ID', 'ID'),
                    
                    data = to_plot, 
                    test = 't', 
                    dfs = 'contain',
                    control=list(optimizer=optimizer)) # control=list(iter.max=10000, rel.tol=1e-9), verbose = TRUE) 
  }else{
    if(mtrgrsn){
      site_rgrs <- rgr_df[, !(names(rgr_df) %in% c('Treatment', 'UnitID', 'UniqueID'))] %>% 
        group_by(Site) %>% 
        summarise(across(everything(), mean))
      
      to_plot <- left_join(to_plot, 
                           site_rgrs, 
                           join_by(site == Site))
      
      ## Center regressor
      to_plot[[regressor]] = scale(to_plot[[regressor]], center = TRUE, scale = FALSE)
    }
    
    
    ## Metafor setup: 
    mfmod <- rma.mv(yi = mean, V = cv_pft, mods = formula, 
                    # random = ~variable|site, 
                    # struct = 'CS',
                    
                    # Alternative formulation
                    random = list(~1|site, ~variable|site),
                    struct = c('ID', 'ID'),
                    
                    data = to_plot, test = 't', dfs = 'contain',
                    control=list(optimizer=optimizer)) #iter.max=1000, rel.tol=1e-5, , verbose = TRUE)
  }
  
  ## Apply cluster-robust inference methods (robust variance estimation)
  ## Note: use the improved methods from the clubSandwich package
  
  r_mfmod <- robust(mfmod, cluster = site, clubSandwich = clubSandwich)
  
  return(r_mfmod)
}

### Function that takes a single metafor model object and a specified pft  
### and returns summary statistics, including mean parameter values, and standard 
### deviation, confidence intervals, and pValues. 
summarize_mfmod <- function(mf_mod, pft, p_dict = pft_dict){
  mvsum <- summary(mf_mod)
  df <- data.frame(Treatment = rownames(mvsum$beta), 
                   Estimate = mvsum$beta[,1],
                   ci.lb = mvsum$ci.lb,
                   ci.ub = mvsum$ci.ub,
                   pval = mvsum$pval,
                   PFT = p_dict[[pft]])
  rownames(df) <- NULL
  return(df)
}

### Function that runs meta-analyses and summarizes results for all PFTs in the 
### input dataframe (full data).
make_meta_df <- function(full_data, 
                         cv_list, # list of var-cov matrices associated with full_data
                         formula, # Meta-analysis or meta-regression formula
                         mtrgrsn = FALSE, # Meta-regression? T/F
                         regressor = NULL, # If meta-regression, specify metaregressor column
                         rgr_df = NULL, # If meta-regression, provide metaregressor dataframe
                         p_dict=pft_dict, # Mapping of raw pft names to aliases
                         by_unit=FALSE, # Meta-analysis by unit? T/F,
                         optimizer = 'Nelder-Mead'
){
  ## Instantiate lists to hold all meta-analysis models
  allOutcomes <- vector('list', length = length(unique(full_data$pft)))
  allModels <- list()
  
  ## Loop through PFTs
  for(p in unique(full_data$pft)){
    ## Run meta-analysis for one PFT
    one_pft <- mf_meta_one_pft(full_data = full_data, 
                               pft_to_plot = p, 
                               formula = formula, 
                               cv_list = cv_list, 
                               p_dict = p_dict, 
                               by_unit = by_unit,
                               mtrgrsn = mtrgrsn,
                               regressor = regressor, 
                               rgr_df = rgr_df,
                               optimizer = optimizer
    )
    allModels[[p]] <- one_pft
    
    # Make a neat dataframe
    mvsum <- summary(one_pft)
    pft_df <- data.frame(Treatment = rownames(mvsum$beta), 
                         Estimate = mvsum$beta[,1],
                         ci.lb = mvsum$ci.lb,
                         ci.ub = mvsum$ci.ub,
                         pval = mvsum$pval,
                         PFT = p_dict[[p]])
    rownames(pft_df) <- NULL
    allOutcomes[[p]] <- pft_df
  }
  
  ## Combine results for all PFTs
  all_pfts <- do.call(rbind, allOutcomes)
  
  # Clean treatment names
  all_pfts$Effect = NaN
  all_pfts[all_pfts$Treatment=='intrcpt', 'Effect'] = 'Control'
  all_pfts[all_pfts$Treatment=='BurnBurn', 'Effect'] = 'Burn'
  all_pfts[all_pfts$Treatment=='ThinThin', 'Effect'] = 'Thin'
  all_pfts[all_pfts$Treatment=='BurnBurn:ThinThin', 'Effect'] = 'Burn:Thin Interaction'
  
  out <- list(all_pfts, allModels)
  
  return(out)
}

### Function that estimates metaanalysis models and uses the models to predict
### growth for each treatment.
get_growth_predictions <- function(full_data, 
                                   cv_list, # list of var-cov matrices associated with full_data
                                   formula, # Meta-analysis or meta-regression formula
                                   by_unit = TRUE, # Meta-analysis by unit? T/F
                                   pfts = c('cedar', 'fir', 'white_pine', 'yellow_pine')) # PFTs to return
  {
  
  ## Run meta-analysis for all PFTs
  meta_results <- make_meta_df(full_data = full_data,
                               cv_list = cv_list,
                               formula = formula,
                               by_unit = by_unit)
  
  ## Create a dataframe of treatments (independent variables) to input to growth predictions
  to_predict <- data.frame(BurnBurn=c(0,1,0,1), 
                           ThinThin=c(0,0,1,1), 
                           BurnBurnThinThin=c(0,0,0,1))
  
  ## Loop over PFTs and predict growth for each treatment
  out <- vector(mode = 'list', length = 4)
  
  for(p in 1:length(pfts)){
    new_pred <- predict(meta_results[[2]][[pfts[p]]], newmods = as.matrix(to_predict))
    new_pred$Treatment <- c('Control', 'Burn', 'Thin', 'Burn+Thin')
    new_pred$PFT <- pft_dict[[pfts[p]]]
    out[[p]] <- as.data.frame(new_pred)
  }
  
  ## Bind PFT-specific outputs into one table
  out_df <- do.call(rbind, out)
  
  return(out_df)
}

### Wrapper function that estimates meta-analyses and returns only the summary 
### table, not the individual metafor models.
get_trt_effects_tbl <- function(full_data, 
                         cv_list, # list of var-cov matrices associated with full_data
                         formula, # Meta-analysis or meta-regression formula
                         by_unit = TRUE) # Meta-analysis by unit? T/F)
  {
  meta_results <- make_meta_df(full_data = full_data,
                               cv_list = cv_list,
                               formula = formula,
                               by_unit = TRUE)
  
  return(meta_results[[1]])
}

