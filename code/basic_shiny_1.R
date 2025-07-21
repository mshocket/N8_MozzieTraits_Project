#load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

library(paletteer)

#data in####

#load data
dat <- read.csv("../data/MDR_and_pLA_Data.csv", stringsAsFactors = T)

selected_TPC <- c("quadratic_2008","briere1simplified_1999","thomas_2017", "deutsch_2008", "atkin_2005") #this worked without # as a list
selected_species <- c("Aalb","Cpip","Ctar","Cqui")
selected_trait <- "pLA"

#define small functions ####

#make get pararmeters function
get_parameters <- function(model){
  mod_params <- c(AIC(model), get_ctmin(model), get_ctmax(model), get_topt(model), get_breadth(model), get_rmax(model), get_thermaltolerance(model), get_thermalsafetymargin(model), "param list")
  return(mod_params)    
}

#define function to extract string of model coeficients
coef_string <- function(model){
vec <- character(length(coef(model)))

for(i in 1:length(coef(model))){
  vec[i] <- paste(names(coef(model))[i], signif(coef(model)[i], digits = 3), sep = ": ")
}
coefs_out <- paste(vec, collapse = ", ")
return(coefs_out)
}


# Define the main function ####

shiny_TPC = function(trait, species){
 
  # set up ----
   
  #create df of selected trait/species - from interface, will change in shiny
  filt_data <- dat |>
  filter(trait.name == selected_trait)|>
  filter(host.code %in% selected_species) 
  
  
  #create first column of predictions df 
  preds <- data.frame(T.C = seq(min(filt_data$T.C), max(filt_data$T.C), length.out = 100))  
  
  
  #create empty output table,specify data type to get empty rows
  output_table <- data.frame(Model = character(),
                             AIC = numeric(),
                             ct_min = numeric(),
                             ct_max = numeric(),
                             T_opt = numeric(),
                             breadth = numeric(),
                             r_max = numeric(),
                             thermal_tolerance = numeric(),
                             safetey_margin = numeric(),
                             coeficients = character())


  
  # loop through each model----
  
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
   

      output_table <- output_table |>
        add_row(Model = tpc,
                AIC = AIC(mod),
                ct_min = get_ctmin(mod),
                ct_max = get_ctmax(mod),
                T_opt = get_topt(mod),
                breadth = get_breadth(mod),
                r_max = get_rmax(mod),
                thermal_tolerance = get_thermaltolerance(mod),
                safetey_margin = get_thermalsafetymargin(mod),
                coeficients = coef_string(mod))
      
      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, atkin_2005 = .fitted)
    }
    else if(tpc == "quadratic_2008") {
      
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = "quadratic_2008")
      
      # fit model
      mod <- nls_multstart(trait~quadratic_2008(temp = T.C, a, b, c),
                           data = filt_data,
                           iter = 200,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = 'quadratic_2008'),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = 'quadratic_2008'),
                           supp_errors = 'Y',
                           convergence_count = FALSE)
      
      
       # save_mod_params <- c("tpc",get_parameters(mod))
       # output_table |> add_row(save_mod_params)
      
      output_table <- output_table |>
        add_row(Model = tpc,
                AIC = AIC(mod),
                ct_min = get_ctmin(mod),
                ct_max = get_ctmax(mod),
                T_opt = get_topt(mod),
                breadth = get_breadth(mod),
                r_max = get_rmax(mod),
                thermal_tolerance = get_thermaltolerance(mod),
                safetey_margin = get_thermalsafetymargin(mod),
                coeficients = coef_string(mod))
      
      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, quadratic_2008 = .fitted)
    }
    else if(tpc == "quadratic_2008") {
      
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = "quadratic_2008")
      
      # fit model
      mod <- nls_multstart(trait~quadratic_2008(temp = T.C, a, b, c),
                           data = filt_data,
                           iter = 200,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = 'quadratic_2008'),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = 'quadratic_2008'),
                           supp_errors = 'Y',
                           convergence_count = FALSE)
      
      
      # save_mod_params <- c("tpc",get_parameters(mod))
      # output_table |> add_row(save_mod_params)
      
      output_table <- output_table |>
        add_row(Model = tpc,
                AIC = AIC(mod),
                ct_min = get_ctmin(mod),
                ct_max = get_ctmax(mod),
                T_opt = get_topt(mod),
                breadth = get_breadth(mod),
                r_max = get_rmax(mod),
                thermal_tolerance = get_thermaltolerance(mod),
                safetey_margin = get_thermalsafetymargin(mod),
                coeficients = coef_string(mod))
      
      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, quadratic_2008 = .fitted)
    }
  }
    
  # create outputs ------
 
  #plot data point only
  plot_bare <- ggplot(filt_data, aes(T.C, trait)) +
    geom_point() +
    theme_bw(base_size = 12) +
    labs(x = 'Temperature (ºC)',
         y = 'probability of survival from larvae to adult',
         title = '')
  
  # plot with curves
  curve_plot <- ggplot(preds, aes(x = T.C, y = trait)) +
    geom_point(aes(T.C, trait), data = filt_data) +
    geom_line(aes(T.C, atkin_2005), col = "red") +
    geom_line(aes(T.C, quadratic_2008), col = "blue")+
    geom_hline(aes(yintercept = 0), linetype = 2) +
    theme_bw(base_size = 12) +
    labs(x = 'Temperature (ºC)',
         y = selected_trait,
         title = '')   
  
  return(list(plot_bare, curve_plot, output_table))
}

saved_output <- shiny_TPC(selected_trait,selected_species)
plot(saved_output[[1]])
plot(saved_output[[2]])
print(saved_output[[3]])






