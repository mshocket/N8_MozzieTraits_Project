#load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(viridis) 
library(RColorBrewer)
library(ggrepel)


#to add a new models - add to model options
                  #  - add to mod rename function
                  #  - add fitting model to list of ifs

#data in####

#load data
dat <- read.csv("../data/MDR_and_pLA_Data.csv", stringsAsFactors = T)



show_lowest_aic <- "yes"

selected_TPC <- c("deutsch_2008", "atkin_2005", "briere1simplified_1999") #this worked without # as a list
selected_species <- c("Aalb","Cpip","Ctar","Cqui")
selected_trait <- "MDR"

model_options <- c("atkin_2005",  "quadratic_2008","thomas_2017", "briere1simplified_1999", "deutsch_2008", "flextpc_2024")

#"sharpeschoolhigh_1981", "sharpeschoollow_1981", #"sharpeschoolfull_1981")

species_options <- levels(dat$host.code)
trait_options <- levels(dat$trait.name)

#define small functions ####


add_row_tpc = function(table_to_add_to, TPC, model) {
  output_table_added <- table_to_add_to |>
    add_row(Model = TPC,
            AIC = AIC(model),
            ct_min = get_ctmin(model),
            ct_max = get_ctmax(model),
            T_opt = get_topt(model),
            breadth = get_breadth(model, level = 0.1),
            r_max = get_rmax(model),
            thermal_tolerance = get_thermaltolerance(model),
            safety_margin = get_thermalsafetymargin(model),
            coeficients = coef_string(model))
  return(output_table_added)
}

#make get pararmeters function
get_parameters <- function(model){
  mod_params <- c(AIC(model), get_ctmin(model), get_ctmax(model), get_topt(model), get_breadth(model, val = 0.1), get_rmax(model), get_thermaltolerance(model), get_thermalsafetymargin(model), "param list")
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

#function to create a vector of tidy names
model_rename <- function(model_name){
  if(model_name == "atkin_2005") {"Atkin"}
  else if(model_name == "quadratic_2008") {"Quadratic"}
  else if(model_name == "thomas_2017") {"Thomas"}
  else if(model_name == "briere1simplified_1999") {"Briere"}
  else if(model_name == "deutsch_2008") {"Deutsch"}
  else if(model_name == "flextpc_2024") {"Flex TPC"}
  # else if(model_name == "sharpeschoolhigh_1981") {"Sharpe-School high"}
  # else if(model_name == "sharpeschoollow_1981") {"Sharpe-School low"}
  # else if(model_name == "sharpeschoolfull_1981") {"Sharpe-School full"}
  
}

trait_rename <- function(trait_name) {
  if(trait_name == "pLA") {"Proportion survived from larvae to adulthood"}
  else if(trait_name == "MDR") {"Mosquito development rate (days^1)"}
}
expression(paste(" ",R^2 ,"= 0.647"))
# Define the main function ####

shiny_TPC = function(TRAIT, SPECIES, MODEL){
 
  # set up ----
   
  #create df of selected trait/species - from interface, will change in shiny
  filt_data <- dat |>
  filter(trait.name == TRAIT)|>
  filter(host.code %in% SPECIES)
  

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
                             safety_margin = numeric(),
                             coeficients = character())

  n_iterations <- 200
  
  # loop through each model----
  
  for (tpc in MODEL){
    
    #atkin still has mod name instead of tpc for checking 
    if(tpc == "atkin_2005"){

      #fitting model
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = "atkin_2005")

      # fit model
      mod <- nls_multstart(trait~atkin_2005(temp = T.C, r0, a, b),
                           data = filt_data,
                           iter = n_iterations,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = "atkin_2005"),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = "atkin_2005"),
                           supp_errors = 'Y',
                           convergence_count = FALSE)
      
      # add parameter row to table
      output_table <- add_row_tpc(output_table, "atkin_2005", mod)

      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, atkin_2005 = .fitted)
    }
    else if(tpc == "quadratic_2008") {
      
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = tpc)
      
      # fit model
      mod <- nls_multstart(trait~quadratic_2008(temp = T.C, a, b, c),
                           data = filt_data,
                           iter = n_iterations,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           supp_errors = 'Y',
                           convergence_count = FALSE)
      
      
      # add parameter row to table
      output_table <- add_row_tpc(output_table, tpc, mod)
      
      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, quadratic_2008 = .fitted)
    }
    else if(tpc == "briere1simplified_1999"){

      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = tpc)

      # fit model
      mod <- nls_multstart(trait~briere1simplified_1999(temp = T.C, tmin, tmax, a),
                           data = filt_data,
                           iter = n_iterations,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           supp_errors = 'Y',
                           convergence_count = FALSE)


      # add parameter row to table
      output_table <- add_row_tpc(output_table, tpc, mod)

      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, briere1simplified_1999 = .fitted)
    }
    else if (tpc == "thomas_2017") {

    start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = tpc)

    # fit model
    mod <- nls_multstart(trait~thomas_2017(temp = T.C, a, b, c, d, e),
                         data = filt_data,
                         iter = n_iterations,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                         upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                         supp_errors = 'Y',
                         convergence_count = FALSE)

    # add parameter row to table
    output_table <- add_row_tpc(output_table, tpc, mod)


    # get predictions
    preds <- broom::augment(mod, newdata = preds)
    preds <- rename(preds, thomas_2017 = .fitted)
    }
    else if (tpc == "deutsch_2008") {
      
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = tpc)
      
      # fit model
      mod <- nls_multstart(trait~deutsch_2008(temp = T.C, rmax, topt, ctmax, a),
                           data = filt_data,
                           iter = n_iterations,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           supp_errors = 'Y',
                           convergence_count = FALSE)
      
      # add parameter row to table
      output_table <- add_row_tpc(output_table, tpc, mod)
      
      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, deutsch_2008 = .fitted)
      }
    else if (tpc == "flextpc_2024") {      
      
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = tpc)
    
    # fit model
    mod <- nls_multstart(trait~flextpc_2024(temp = T.C, tmin, tmax, rmax, alpha, beta),
                         data = filt_data,
                         iter = n_iterations,
                         start_lower = start_vals - 10,
                         start_upper = start_vals + 10,
                         lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                         upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                         supp_errors = 'Y',
                         convergence_count = FALSE)
    
    # add parameter row to table
    output_table <- add_row_tpc(output_table, tpc, mod)
    
    # get predictions
    preds <- broom::augment(mod, newdata = preds)
    preds <- rename(preds, flextpc_2024 = .fitted)}
    else if (tpc == "sharpeschoolhigh_1981") {      
      
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = tpc)
      
      # fit model
      mod <- nls_multstart(trait~sharpeschoolhigh_1981(temp = T.C, r_tref, e, eh, th, tref),
                           data = filt_data,
                           iter = n_iterations,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           supp_errors = 'Y',
                           convergence_count = FALSE)
      
      # add parameter row to table
      output_table <- add_row_tpc(output_table, tpc, mod)
      
      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, sharpeschoolhigh_1981 = .fitted)}
    else if (tpc == "sharpeschoollow_1981") {      
      
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = tpc)
      
      # fit model
      mod <- nls_multstart(trait~sharpeschoollow_1981(temp = T.C, r_tref, e, el, tl, tref =20),
                           data = filt_data,
                           iter = n_iterations,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           supp_errors = 'Y',
                           convergence_count = FALSE)
      
      # add parameter row to table
      output_table <- add_row_tpc(output_table, tpc, mod)
      
      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, sharpeschoollow_1981 = .fitted)}
    else if (tpc == "sharpeschoolfull_1981") {      
      
      start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = tpc)
      
      # fit model
      mod <- nls_multstart(trait~sharpeschoolfull_1981(temp = T.C, r_tref, e, el, tl, eh, th, tref = 20),
                           data = filt_data,
                           iter = n_iterations,
                           start_lower = start_vals - 10,
                           start_upper = start_vals + 10,
                           lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = tpc),
                           supp_errors = 'Y',
                           convergence_count = FALSE)
      
      # add parameter row to table
      output_table <- add_row_tpc(output_table, tpc, mod)
      
      # get predictions
      preds <- broom::augment(mod, newdata = preds)
      preds <- rename(preds, sharpeschoolfull_1981 = .fitted)}
    
  }
    

  # create outputs ------
  
  #create long df with col of mod names and col of preds (temps repeated)
  long_preds <- preds |>
    pivot_longer(cols = - T.C,
                 names_to = "model_name",
                 values_to = "pred_val")
  
  #create new column with tidy version of model name, uses model rename function defined above
  long_preds <- long_preds |>
    mutate(tidy_name = map_chr(long_preds$model_name, model_rename)
    )
  
  filt_data <- filt_data |>
    mutate(y_axis_name = map_chr(filt_data$trait.name, trait_rename)
    )
  
  
  #plot 
  
  base_plot <-  ggplot(filt_data, aes(T.C, trait)) +
    geom_point(alpha = 1/5, shape = 16) +
    geom_hline(aes(yintercept = 0), linetype = 2, color = "grey30") +
    theme_bw(base_size = 16) +
    labs(x = 'Temperature (ºC)',
         y = filt_data$y_axis_name) 
    
  
  curve_plot <- base_plot +
    geom_line(data=long_preds, aes(x = T.C, y = pred_val, colour= tidy_name), lwd = 0.8) +
    labs(title = '',
         color = "Model") +
    scale_color_brewer(type = 'qual', palette = 2)
    # scale_color_viridis(discrete= T, option = "viridis")

  #not important but trying to figure out, maybe add a col to long preds? is it usefull to still be able to tell the other models appart
  best_model = filter(output_table, AIC == min(AIC)) |> 
    pull(Model)
  
  for (i in 1:nrow(long_preds)) {
    if(long_preds$model_name[i] == best_model) {"best"}
    else {"not best"}
  }
  
  best_mod_plot <- base_plot +
    geom_line(data = filter(long_preds, model_name == best_model), aes(x = T.C, y = pred_val), lwd = 2, col = "yellow") +
    geom_line(data=long_preds, aes(x = T.C, y = pred_val, col = tidy_name), lwd = 0.8, linetype = 2 ) +
    labs(title = 'best mod',
         color = "Model") +
    scale_color_brewer(type = 'qual', palette = 2)
    #scale_color_viridis(discrete= T, option = "viridis")
  
    
  return(list(curve_plot,
              output_table, 
              best_mod_plot
              ))
  }


saved_output <- shiny_TPC(selected_trait, selected_species, selected_TPC)
print(saved_output[[1]])
print(saved_output[[2]])


print(saved_output)




