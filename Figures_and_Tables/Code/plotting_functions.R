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


