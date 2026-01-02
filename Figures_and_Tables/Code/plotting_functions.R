## Function that takes in the a set of model outputs and PFT and returns a 
## PFT-specific column of for a summary table of meta-analysis estimates and uncertainty
make_pft_results_column <- function(
    all_pft_models, # This is the first output of `make_meta_full`
    p, # pft
    by_unit # Meta-analysis was at the unit-level? T/F
){
  
  ## Subset PFT model of interest
  mod <- all_pft_models[[p]]
  
  ## Make a vector containing asterisk notation of significance levels for each parameter
  signf <- vector('character', length(mod$pval))
  for(i in 1:length(mod$pval)){
    if(signif(mod$pval[i], 2) <= 0.05){
      signf[i] = '**'
    }else if(signif(mod$pval[i], 2) <= 0.1){
      signf[i] = '*'
    }else{
      signf[i] = ''
    }
  }
  
  ## Each cell of the table summarizes the mean parameter estimate, significance
  ## level, confidence interval, and p-value
  pft_out <- paste0(signif(mod$beta, 2), 
                    signf, '\n[', 
                    signif(mod$ci.lb, 2), ', ', signif(mod$ci.ub, 2), 
                    '] \np-value=', signif(mod$pval, 2),  sep = '')
  
  ## Add a row for tau^2 estimates for random effects variance-covariance matrices  
  pft_out[length(pft_out) + 1] <- signif(mod$sigma2, 2)
  pft_out[length(pft_out) + 1] <- signif(mod$tau2, 2)
  
  if(by_unit){
    pft_out[length(pft_out) + 1] <- signif(mod$gamma2, 2)
  }
  
  return(pft_out)
}



### Function to loop through all PFT models and construct a summary table
make_results_table <- function(
    all_pft_models, # This is the first output of `make_meta_full`
    by_unit = FALSE # Meta-analysis was at the unit-level? T/F
){
  
  ## Instantiate a list of columns
  col_list <- list()
  
  ## Create a column for each PFT
  for(p in names(all_pft_models)){
    col_list[[pft_dict[[p]]]] <- make_pft_results_column(all_pft_models = all_pft_models, p = p,
                                                         by_unit=by_unit)
  }
  
  # Get row names
  if(by_unit){
    rnames <- c(rownames(all_pft_models[[1]]$beta), 'sigma2', 'tau2', 'gamma2')
    rnames[rnames == 'sigma2'] = "$\\tau^2_{site}$"
    rnames[rnames == 'tau2'] = "$\\tau^2_{treatment}$"
    rnames[rnames == 'gamma2'] = "$\\tau^2_{unit}$"
  }else{
    rnames <- c(rownames(all_pft_models[[1]]$beta), 'sigma2', 'tau2')
    rnames[rnames == 'sigma2'] = "$\\tau^2_{site}$"
    rnames[rnames == 'tau2'] = "$\\tau^2_{treatment}$"
  }
  
  rnames[rnames == 'intrcpt'] = 'Intercept'
  rnames <- str_replace_all(rnames, 'BurnBurn', 'Burn')
  rnames <- str_replace_all(rnames, 'ThinThin', 'Thin')
  
  ## Compile as a dataframe
  out <- data.frame(Parameter = rnames, col_list)
  
  ## Clean up
  colnames(out)[colnames(out) == 'Yellow.Pine'] = 'Yellow Pine'
  colnames(out)[colnames(out) == 'Sugar.Pine'] = 'Sugar Pine'
  
  ## Subset columns of intrest
  out <- out[c('Parameter', 'Cedar', 'Fir', 'Sugar Pine', 'Yellow Pine')]
  
  ## Flip naming convention
  out[out$Parameter == 'Burn:Thin', 'Parameter'] = 'Thin:Burn'
  
  return(out)
}



### Function to plot fixed effects meta-analysis (each panel of Figure 2)
plot_meta_analysis <- function(effects_summary, # This is the first output of `make_meta_full`
                               axis_label, # x-axis label
                               control=TRUE, # Plot results in the control (T) or in the treatments (F)
                               zero_line = TRUE # Include a vertical line to delineate 0? T/F
){ 
  ## Subset data to either plot the control or to plot all treatments
  if(control){
    df_to_plot <- effects_summary[effects_summary$Effect == 'Control',]
  }else{
    df_to_plot <- effects_summary[effects_summary$Effect != 'Control',]
  }
  
  p <- ggplot(df_to_plot,
              aes(x= '', y=Estimate, 
                  ymin=ci.lb, ymax=ci.ub, 
                  col=PFT, fill=PFT)) +
    
    geom_linerange(linewidth=5, position=position_dodge(width = 0.5)) +
    geom_point(size=3, shape=21, colour="white", stroke = 0.5, 
               position=position_dodge(width = 0.5)) +
    
    scale_fill_manual(values=pft_cols, name='PFT')+
    scale_color_manual(values=pft_cols, name='PFT') +
    scale_x_discrete(name="") +
    ylab(axis_label) +
    
    # Styling
    coord_flip() + theme_minimal(base_size = 14) +
    theme(axis.line = element_line(color='black'), 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = 'bottom')
  
  if(zero_line){
    p <- p + geom_hline(yintercept=0, lty=2, col = 'black')
  }
  
  if(control){
    p <- p + 
      labs(subtitle = 'Control') +
      theme(plot.subtitle = element_text(hjust = 0.5))
  }else{
    p <- p + facet_wrap(~factor(Effect, 
                                c('Burn', 'Thin', 'Thin:Burn Interaction')))
  }
  
  return(p)
}

## Function to make a forest-like plot the estimated SSM parameters by site, 
## Option to include meta-analysis
plotSSMparam <- function(df, ## Dataframe that includes parameter estimates
                         effect_to_plot, ## Parameter to plot
                         label_dict = pretty_nms, ## Mapping of parameter names in df to plotting names
                         plot_metanalysis=FALSE, ## Plot meta-analysis? T/F
                         meta_df = NULL, ## If plot_metanalysis is TRUE, provide meta-analysis results DF
                         meta_effect = NULL ## If plot_metanalysis is TRUE, provide name of effect to plot
){
  
  ## Subset parameter to plot
  to_plot <- df[df$variable == effect_to_plot,]
  
  p1 <- ggplot(to_plot,
               aes(x=common_nm, y=mean, ymin=low_ci, ymax=high_ci,col=PFT,fill=PFT)) + 
    
    # black line is 0
    geom_hline(yintercept=0, lty=2, col = 'black') #+
  
  ## If plotting meta-analysis, include lines for meta-analysis means and CIs
  if(plot_metanalysis){
    for(p in unique(df$PFT)){
      p1 <- p1 +     
        geom_hline(yintercept = meta_df[meta_df$PFT==p & meta_df$Effect==meta_effect, 'Estimate'], 
                   lty = 1, col = pft_cols[p]) +
        geom_hline(yintercept = meta_df[meta_df$PFT==p & meta_df$Effect==meta_effect, 'ci.lb'], 
                   lty = 3, col = pft_cols[p]) +
        geom_hline(yintercept = meta_df[meta_df$PFT==p & meta_df$Effect==meta_effect, 'ci.ub'], 
                   lty = 3, col = pft_cols[p])
    }
  }
  
  ## Subset meta-analysis DF to effect of interest
  meta_bounds <- meta_df[meta_df$Effect == meta_effect,]
  
  ## Plot dots (means) and bars (CIs) for parameter estimate at each site
  p1 <- p1 + 
    geom_linerange(linewidth=5,position=position_dodge(width = 0.75)) +
    geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.75)) +
    scale_fill_manual(values=pft_cols, name='PFT')+
    scale_color_manual(values=pft_cols, name='PFT') +
    scale_x_discrete(name="Site") +
    scale_y_continuous(TeX(label_dict[[effect_to_plot]]), limits = c(min(to_plot$low_ci, meta_bounds$ci.lb), max(to_plot$high_ci, meta_bounds$ci.ub))) +
    
    # Styling
    coord_flip() + theme_minimal(base_size = 14) +
    theme(axis.line = element_line(color='black'), 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank())
  
  return(p1)
}

#' Plots the X vs growth relationship with significance bars at the bottom.
#'
#' @param sig_df The significance dataframe generated by create_significance_df.
#' @param all_predictions All growth predictions to plot
#' @param x_var The unquoted column name for the X-variable (e.g., DBH or CWD).
#' @return A ggplot object.
plot_growth_by_xVar <- function(sig_df,
                                all_predictions,
                                x_var){
  
  ## Make markers where only one X value is significant visible
  sig_df[sig_df$xstart == sig_df$xend, "xend"] = 
    sig_df[sig_df$xstart == sig_df$xend, "xend"] - 2
  
  ## Add y axis locations for plotting the significance bars
  sig_df$yloc <- -2.5
  sig_df$yloc[sig_df$Treatment=='Burn'] = -2.4
  sig_df$yloc[sig_df$Treatment=='Thin:Burn Interaction'] = -2.3
  
  ## Plot X vs log(growth)
  p1 <- ggplot() + 
    geom_point(data = all_predictions,
               aes(x =  {{ x_var }}, y = pred, colour = Treatment), 
               alpha = 0.6) +
    geom_smooth(data = all_predictions, 
                aes(x = {{ x_var }}, y = pred, color = Treatment), se = FALSE) +
    
    ## Segments indicate ranges over which treatment effects are significant
    geom_segment(data = sig_df,
                 aes(x = xstart, 
                     xend = xend, y=yloc, yend = yloc, 
                     color = Treatment),
                 inherit.aes = FALSE, linewidth = 3, linetype = 1, alpha = 0.3,
                 show.legend = FALSE) +
    scale_color_manual(values=unlist(trt_cols), name='Treatment', 
                       breaks = c('Control', 'Burn', 'Thin', 'Thin+Burn')) +
    
    ylab("Log DBH increment") +
    
    # Styling
    theme_minimal(base_size = 14) +
    theme(axis.line = element_line(color='black'), 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = 'bottom')
  
  return(p1 + facet_wrap(~PFT))
  
}

## Function to plot meta-regression lines, faceted by PFT, on top of data points
## Option to map point size to weight of data point in the meta-analysis
plot_weightedRegression_allTreats <- function(df_to_plot, ## dataframe containing input data and weights to the meta-regression model
                                              rgr, ## name of regressor
                                              xlab, ## x-axis label
                                              plim=c(0.5, 2), ## min and max marker size
                                              regLine_df, ## dataframe with regression line data
                                              by_unit, ## Meta-analysis by unit? T/F
                                              uncenter = TRUE, ## Uncenter metaregressor? T/F
                                              oneRegLine = FALSE ## Only TRUE if treatment is not in the meta-regression model
){
  
  ## Scale point sizes to the plim
  wi <- df_to_plot$psize
  rng <- max(wi) - min(wi)
  psize <- (wi - min(wi)) / rng
  psize <- (psize * (plim[2] - plim[1])) + plim[1]
  df_to_plot$psizeScaled <- psize
  
  ## Change name of control plot
  df_to_plot$Treatment[df_to_plot$Treatment == 'None'] = 'Control'
  df_to_plot$Treatment[df_to_plot$Treatment == 'Burn+Thin'] = 'Thin+Burn'
  
  ## Uncenter metaregressor
  if(uncenter){
    df_to_plot$rgr_uncntred <- df_to_plot[[rgr]] + df_to_plot$rgr_cntr
    var_nm <- 'rgr_uncntred'
  }else{
    var_nm <- rgr
  }
  
  ## Plot scatter plot
  if(by_unit){
    p <- ggplot(df_to_plot,
                aes(x=get(var_nm), y=mean, 
                    color=Treatment, 
                    shape = common_nm)
    ) +
      # geom_hline(yintercept=0, lty=2) + ## Option to include a horizontal line at y = 0
      geom_point(position=position_dodge(width = sd(df_to_plot[[var_nm]])/100), 
                 alpha = 0.4, 
                 size = 3) +
      scale_color_manual(values = unlist(trt_cols), 
                         breaks = c('Control', 'Burn', 'Thin', 'Thin+Burn') 
      ) + 
      ## Option to map marker shape to treatment 
      # scale_shape_manual(values=trt_shapes, breaks = c('Control', 'Burn', 'Thin', 'Burn+Thin'))
      
      ## Option to map marker shape to site
      scale_shape_manual(values=site_shapes, name = 'Site', 
                         breaks = c('Blodgett', 'LaTour', 'Sequoia', 'STEF', 'Teakettle',
                                    "Tharp's Creek", "W. Lake Tahoe"))
    
    
  }else{
    p <- ggplot(df_to_plot,
                aes(x=get(var_nm), y=mean, 
                    # size=psizeScaled, 
                    color=Treatment,
                    shape = common_nm)) +
      # geom_hline(yintercept=0, lty=2) + ## Option to include a horizontal line at y = 0
      geom_point(position=position_dodge(width = sd(df_to_plot[[var_nm]])/100), 
                 alpha = 0.4, 
                 size = 3) +
      
      scale_color_manual(values = unlist(trt_cols), 
                         breaks = c('Control', 'Burn', 'Thin', 'Thin+Burn') 
      ) + 
      
      scale_shape_manual(values=site_shapes, name = 'Site', 
                         breaks = c('Blodgett', 'LaTour', 'Sequoia', 'STEF', 'Teakettle',
                                    "Tharp's Creek", "W. Lake Tahoe"))  
  }
  
  ## Plot regression lines
  if(nrow(regLine_df!=0)){
    if(uncenter){
      regLine_df[[var_nm]] <- regLine_df[[rgr]] + regLine_df$rgr_cntr
    }
    
    ## If treatment is not a variable in the meta-regression model
    if(oneRegLine){ 
      p <- p + geom_line(data = regLine_df, inherit.aes = FALSE,
                         aes(x=get(var_nm), y=pred), 
                         linewidth = 1.0, color = 'black') 
      ## Option to include a ribbon
      # geom_ribbon(data = regLine_df, aes(x =get(var_nm), y=pred, ymin = ci.lb, ymax = ci.ub, fill=PFT), alpha = 0.1, color=NA, show.legend = FALSE, size = 0.5) +
      
      ## Option to color by PFT
      # scale_color_manual(values=unlist(pft_cols), guide = FALSE) #+
      # scale_fill_manual(values=unlist(pft_cols), guide = FALSE)
    }else{
      p <- p + geom_line(data = regLine_df, inherit.aes = FALSE,
                         aes(x=get(var_nm), y=pred, color=Treatment), # Alternative: linetype = Treatment
                         linewidth = 1.0) + 
        scale_color_manual(values = unlist(trt_cols), 
                           breaks = c('Control', 'Burn', 'Thin', 'Thin+Burn') 
        ) #+
      ## Option to style lines by treatment
      # scale_linetype_manual(values = unlist(trt_styles), breaks = c('Control', 'Burn', 'Thin', 'Burn+Thin'))
    }
    ## If there are not significant regression lines, just return the scatterplot
  }else{
    p <- p + scale_color_manual(values = unlist(trt_cols), 
                                breaks = c('Control', 'Burn', 'Thin', 'Thin+Burn') 
    )
  }
  
  ## Stylying
  p <- p + xlab(xlab) + ylab(paste0("Log DBH increment")) +
    theme_minimal(base_size = 14) +
    theme(axis.line = element_line(color='black'),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.position = 'bottom', 
          legend.direction = 'horizontal',
          legend.box = 'vertical')
  
  return(p) #+ rremove('ylab'))}
}


## Function takes in a set of meta-regression model outputs and plots the meta-regression results
plotMetaRegression <- function(
    allPFToutputs, ## full set of outputs from `make_meta_df`
    rgr, ## name of regressor
    sgnf, ## significance level
    uncenter = TRUE, ## Uncenter metaregressor? T/F
    by_unit=FALSE ## Metaregression by unit? T/F
){
  
  ## Get weighted input data
  df_w_weights <- make_weights_df(allPFTmodels = allPFToutputs[[2]], 
                                  rgr = rgr)
  
  ## Construct regression line data
  reg_lines <- make_regLine_allPFTs(
    allPFToutputs = allPFToutputs, 
    rgr = rgr, 
    sgnf = sgnf)
  
  ## Plot meta-regression; facet by PFT
  p <- plot_weightedRegression_allTreats(
    df_to_plot = df_w_weights, 
    rgr = rgr , 
    xlab = pretty_mrs[[rgr]], 
    regLine_df = reg_lines, 
    uncenter = uncenter, 
    by_unit = by_unit)
  
  return(p+facet_wrap(~factor(PFT, 
                              levels = c('Cedar', 'Fir', 'Sugar Pine', 'Yellow Pine')), 
                      ncol=2))
}

## Function to clean up regression labels in summary tables
beautify_table_labels <- function(
    mr_tbl, ## metaregression results table
    rgr ## name of regressor
){
  mr_tbl$Parameter <- str_replace_all(mr_tbl$Parameter, rgr, 
                                      corrplot_shorthand[[rgr]])
  return(mr_tbl)
}

