
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
                           join_by(common_nm == common_nm, 
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
  all_pfts[all_pfts$Treatment=='BurnBurn:ThinThin', 'Effect'] = 'Thin:Burn Interaction'
  
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
    new_pred$Treatment <- c('Control', 'Burn', 'Thin', 'Thin+Burn')
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


#' Creates a dataframe defining significant X-ranges (e.g., DBH or CWD) for plotting.
#'
#' @param df The main dataframe containing all data.
#' @param x_var The unquoted column name for the X-variable (e.g., DBH).
#' @param alpha The significance threshold.
#' @return A dataframe with columns: treatment, xstart, xend.
identify_significant_ranges <- function(df, x_var, alpha = sgnf) {
  data <- df[df$Effect!='Control',]
  
  # 1. Flag significance
  data_sig <- data %>%
    mutate(
      is_significant = signif(pval, digits = 2) <= alpha
    )
  
  # 2. Generate RLE ID for continuous runs of TRUE/FALSE
  significance_df <- data_sig %>%
    group_by(PFT, Effect) %>%
    mutate(
      rle_id = {
        # Use the flagged column for RLE
        rle_result <- rle(is_significant)
        rep(seq_along(rle_result$lengths), rle_result$lengths)
      }
    ) %>%
    # 3. Filter for significant runs and summarise to find segment boundaries
    filter(is_significant == TRUE) %>%
    group_by(PFT, Effect, rle_id) %>%
    summarise(
      # The Embrace operator {{}} is used here to dynamically select the x_var column
      xstart = min({{ x_var }}),
      xend = max({{ x_var }}),
      .groups = 'drop'
    ) %>%
    # 4. Select and return
    select(PFT, Effect, xstart, xend) %>%
    rename(Treatment = Effect)
  
  return(significance_df)
} 

###-----------------------------------------------------------------------------
### Meta-regression functions
###-----------------------------------------------------------------------------

### Function to extract meta-regression weights associated with each study
get_weights <- function(model, ## Metafor model object
                        rgr ## name of regressor
){
  mod_dat <- model$data
  
  # Make sure that order of points is the same
  if(identical(mod_dat$mean, as.vector(model$yi))==FALSE){
    return(print('Order of weights and yi are not the same.'))
  }else{
    mod_dat$weights <- weights(model)
    
    ## Add regressor center to dataframe for easy transforming later
    mod_dat$rgr_cntr <- attr(mod_dat[[rgr]], 'scaled:center')
    
    return(mod_dat)
  }
}

## Function loops through all PFT-specific model and produces a dataframe of 
## weights associated with each data point in the meta-analysis. 
make_weights_df <- function(
    allPFTmodels, ## Second output of `make_meta_df`
    rgr ## name of regressor
){
  
  allPFTweights <- vector('list', length = length(names(allPFTmodels)))
  weight_ranges <- vector('numeric', length = length(names(allPFTmodels)))
  
  ## Loop through all PFTs and get model weights
  for(p in names(allPFTmodels)){
    pft_mod <- allPFTmodels[[p]]
    pft_df <- get_weights(pft_mod, rgr)
    
    allPFTweights[[p]] <- pft_df
    weight_ranges[[p]] <- max(pft_df$weights) - min(pft_df$weights)
  }
  
  ## Compine weights from all PFT-specific models into one table
  out <- do.call(rbind, allPFTweights)
  
  ## For point sizes, scale weights to the max weight for the PFT with the largest weight range
  pft_to_ref <- names(which.max(unlist(weight_ranges)))
  ref_max <- max(allPFTweights[[pft_to_ref]]$weights)
  
  out$psize <- out$weights
  
  for(p in names(allPFTmodels)){
    if(p!=pft_to_ref){
      pft_max <- max(out[out$pft == p, 'weights'])
      adj <- ref_max / pft_max
      out[out$pft == p, 'psize'] = out[out$pft ==p, 'weights'] * adj
    }
  }
  
  return(out)
}

### Function takes a PFT-specific meta-analysis model and produces a dataframe with 
### the predicted meta-analysis growth model outputs for different combinations of
### regressor and treatment. Only produces lines if the control or a given 
### treatments has a significant interaction with the meta-regressor.
regLine_onePFT <- function(
    allPFToutputs, ## full set of outputs from `make_meta_df`
    rgr, ## name of regressor
    p, ## pft
    sgnf ## significance level
){
  
  ## Subset PFT model of interest
  pmod <- allPFToutputs[[2]][[p]]
  
  ## Get regressor min and max value in control units and create a regressor vector x
  temp_rgr <- pmod$data[pmod$data$Treatment == 'None',]
  x <- seq(min(temp_rgr[[rgr]]), max(temp_rgr[[rgr]]), length.out = 100)
  
  ## Subset p-values for PFT of interest
  out_pft <- allPFToutputs[[1]][allPFToutputs[[1]]$PFT==pft_dict[[p]],]
  
  ## Instantiate an predictor dataframe consisting of a columns for Burn (1 or 0), 
  ## Thin (1 or 0), and regressor (x)
  ## For the control, set Burn = 0 and Thin = 0.
  df <- data.frame(cbind(BurnBurn = 0, ThinThin = 0, var = x))
  colnames(df) <- c('BurnBurn', 'ThinThin', rgr)
  
  ## Loop through pooled effects for the Burn and Thin treatment
  for(trt_lbl in c('BurnBurn', 'ThinThin')){
    
    ## Subset data from Burned or Thinned units
    if(trt_lbl=='BurnBurn'){
      rgr_trt <- pmod$data[pmod$data$Burn == 'Burn',] 
    }else if(trt_lbl=='ThinThin'){
      rgr_trt <- pmod$data[pmod$data$Thin == 'Thin',]
    }
    ## Reconstruct x based on the min and max value of the regressor in the 
    ## treated units
    x <- seq(min(rgr_trt[[rgr]]), max(rgr_trt[[rgr]]), length.out = 100)
    
    ## If the model includes an interaction between the treatment and the metaregressor
    if(paste0(trt_lbl, ":", rgr) %in% out_pft$Treatment){
      
      ## If the interaction term is significant 
      if(signif(out_pft[out_pft$Treatment == paste0(trt_lbl, ":", rgr), 'pval'], 2) <= sgnf){
        
        ## Add rows to the predictor dataframe for the treatment
        df_trt <- data.frame(cbind(BurnBurn = 0, ThinThin = 0, var = x))
        colnames(df_trt) <- c('BurnBurn', 'ThinThin', rgr)
        ## Reset the treatment of interest to 1
        df_trt[[trt_lbl]] = 1
        df <- rbind(df, df_trt)
      }
      ## Otherwise, if a treatment:regressor interaction is not in the model
    }else{ 
      ## Check whether the treatment effect is in the model
      if(trt_lbl %in% out_pft$Treatment){
        ## Check whether the treatment is significant
        if(signif(out_pft[out_pft$Treatment == paste0(trt_lbl), 'pval'], 2) <= sgnf){
          ## If significant, add rows to the predictor dataframe for the treatment
          df_trt <- data.frame(cbind(BurnBurn = 0, ThinThin = 0, var = x))
          colnames(df_trt) <- c('BurnBurn', 'ThinThin', rgr)
          df_trt[[trt_lbl]] = 1
          df <- rbind(df, df_trt)
        }
      }
    }
  }
  
  ## If the model contains a Burn:Thin:regressor interaction
  if(paste0("BurnBurn:ThinThin:", rgr) %in% out_pft$Treatment){
    ## Subset input data from Burn+Thin units
    rgr_trt <- pmod$data[pmod$data$Treatment == 'Burn+Thin',] 
    
    ## Reconstruct x based on the min and max value of the regressor in the 
    ## Burn+Thin units
    x <- seq(min(rgr_trt[[rgr]]), max(rgr_trt[[rgr]]), length.out = 100)
    
    ## If the Burn:Thin:regressor interaction is significant
    if(signif(out_pft[out_pft$Treatment==paste0("BurnBurn:ThinThin:", rgr), 
                      'pval'], 2) <= sgnf){
      
      ## Add rows to the predictor dataframe for the Burn+Thin treatment
      df_trt <- data.frame(cbind(BurnBurn = 0, ThinThin = 0, var = x))
      colnames(df_trt) <- c('BurnBurn', 'ThinThin', rgr)
      df_trt$BurnBurn = 1
      df_trt$ThinThin = 1
      df <- rbind(df, df_trt)
    }
  }
  
  ## If none of the regressor effects are significant, return an empty dataframe
  if(nrow(df) == length(x)){
    if(signif(out_pft[out_pft$Treatment == rgr, 'pval'], 2) > sgnf){
      return(data.frame())
    }
  }
  
  ## If there is a Burn:Thin interaction in the model, add a Burn:Thin column 
  ## to the predictor dataframe by multiplying the Burn and Thin (1 or 0) values
  if('BurnBurn:ThinThin' %in% out_pft$Treatment){
    df[['BurnBurn:ThinThin']] = df$BurnBurn * df$ThinThin
  }
  
  ## If there was a treatment:regressor interaction in the model, add interaction
  ## columns to the predictor dataframe by multiplying the Burn or Thin indicator
  ## by the regressor
  for(trt_lbl in c('BurnBurn', 'ThinThin')){
    if(paste0(trt_lbl, ":", rgr) %in% out_pft$Treatment){
      df[[paste0('BurnBurn:', rgr)]] = df$BurnBurn * df[[rgr]]
      df[[paste0('ThinThin:', rgr)]] = df$ThinThin * df[[rgr]]
    }
  }
  
  ## If there was a Burn:Thin:regressor interaction, add this column to the 
  ## predictor dataframe
  if(paste0('BurnBurn:ThinThin:',rgr) %in% out_pft$Treatment){
    df[[paste0('BurnBurn:ThinThin:',rgr)]] = df$BurnBurn * df[[rgr]] * df$ThinThin
  }
  
  ## Add a Treatment column
  df$Treatment <- "Control"
  df[df$BurnBurn==1 & df$ThinThin == 0, 'Treatment'] = 'Burn'
  df[df$BurnBurn==1 & df$ThinThin == 1, 'Treatment'] = 'Thin+Burn'
  df[df$BurnBurn==0 & df$ThinThin == 1, 'Treatment'] = 'Thin'
  
  ## Use the predictor dataframe to predict growth 
  pred <- predict(allPFToutputs[[2]][[p]], 
                  newmods = as.matrix(df[1:(length(df)-1)])) # Drop treatment column for the prediction
  pred$var <- df[[rgr]]
  pred$Treatment <- df$Treatment ## Add treatment column back to the predicted values
  pred$PFT <- pft_dict[[p]] ## Use PFT aliases
  
  ## Add the regressor center to the dataframe
  pred$rgr_cntr <- attr(pmod$data[[rgr]], 'scaled:center')
  
  ## Return a dataframe
  out <- as.data.frame(pred)
  colnames(out)[colnames(out) == 'var'] = rgr
  return(out)
}


### Function to loop through PFT models and compile a dataframe containing the
### x and y values needed to plot a regression line for each model. 
make_regLine_allPFTs <- function(
    allPFToutputs, ## full set of outputs from `make_meta_df`
    rgr, ## name of regressor
    sgnf ## significance level
){
  
  out_df <- data.frame()
  
  for(p in names(allPFToutputs[[2]])){
    ## Get the significant regression line dataframe for each PFT
    out_df <- rbind(out_df, regLine_onePFT(allPFToutputs = allPFToutputs, 
                                           rgr = rgr, p = p, sgnf = sgnf))
  }
  
  ## Return if the dataframe is empty
  if(nrow(out_df)==0){
    return(data.frame())
  }
  
  ## Beautify treatment names
  out_df[out_df$Treatment =='BurnBurn', 'Treatment'] = 'Burn'
  out_df[out_df$Treatment == 'ThinThin', 'Treatment'] = 'Thin'
  out_df[out_df$Treatment == 'BurnBurn:ThinThin', 'Treatment'] = 'Thin+Burn'
  
  return(out_df)
}
  

