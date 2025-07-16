#load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

#load data
dat <- read.csv("../data/MDR_and_pLA_Data.csv", stringsAsFactors = T)

selected_TPC <- c("quadratic_2008","briere1simplified_1999","thomas_2017", "deutsch_2008", "atkin_2005") #this worked without # as a list
selected_species <- c("Aalb","Cpip","Ctar","Cqui")
selected_trait <- "pLA"


# Define a function
shiny_TPC = function(trait, species){
  
  #create df of selected trait/species - from interface, will change in shiny
  filt_data <- dat |>
  filter(trait.name == selected_trait)|>
  filter(host.code %in% selected_species) 
  
  #create predictions df and results table
  preds <- data.frame(T.C = seq(min(filt_data$T.C), max(filt_data$T.C), length.out = 100))  
  
  col_names <- c("model","AIC", "ct_min", "ct_max", "T_opt", "breadth", "r_max", "thermal_tolerance", "safetey_margin", "parameters")
  r1 <- rep(c(NA), times = length(col_names))
  output_table <- rbind.data.frame(r1)
  colnames(output_table) <- col_names
  
  #plot data point only
  plot_bare <- ggplot(filt_data, aes(T.C, trait)) +
          geom_point() +
          theme_bw(base_size = 12) +
          labs(x = 'Temperature (ºC)',
               y = 'probability of survival from larvae to adult',
               title = '')

  
  #loop through each model
  for (tpc in selected_TPC){
    if(tpc == "atkin_2005"){
      #fitting model
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = 'atkin_2005')
      
      # fit model
      mod <- nls_multstart(trait~atkin_2005(temp = T.C, r0, a, b),
                           data = filt_data,
                           iter = 200,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = 'atkin_2005'),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = 'atkin_2005'),
                           supp_errors = 'Y',
                           convergence_count = FALSE)

      

      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, atkin_2005 = .fitted)
    }
    else if(tpc == "quadratic_2008") {
      #fitting model
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = 'quadratic_2008')
      
      # fit model
      mod2 <- nls_multstart(trait~quadratic_2008(temp = T.C, a, b, c),
                           data = filt_data,
                           iter = 200,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = 'quadratic_2008'),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = 'quadratic_2008'),
                           supp_errors = 'Y',
                           convergence_count = FALSE)
     
      #make a get params function
      
      get_parameters <- function(model){
        mod_params <- c(tpc, AIC(model), get_ctmin(model), get_ctmax(model), get_topt(model), get_breadth(model), get_rmax(model), get_thermaltolerance(model), get_thermalsafetymargin(model), "param list")
    return(mod_params)    
      }
      
     save_mod_params <- get_parameters(mod2)
     output_table <- rbind(output_table, save_mod_params)
      
      # get predictions
      preds <- broom::augment(mod2, newdata = preds)
      preds <- rename(preds, quadratic_2008 = .fitted)
    }
        # else {print("not atkins")
  }
    
 
      # plot
      curve_plot <- ggplot(preds, aes(x = T.C, y = trait)) +
        geom_point(aes(T.C, trait), filt_data) +
        geom_line(aes(T.C, atkin_2005), col = "blue") +
        geom_line(aes(T.C, quadratic_2008), col = "red")+
        geom_hline(aes(yintercept = 0), linetype = 2) +
        theme_bw(base_size = 12)+
        labs(x = 'Temperature (ºC)',
             y = selected_trait,
             title = '')   
  
      output_table2 <- output_table[-1,]
  return(curve_plot)
  #return(output_table2)
}

saved_output <- shiny_TPC(selected_trait,selected_species)
saved_output



# are the loops correct 


