#load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(viridis) 
library(RColorBrewer)

#data in####

#load data
dat <- read.csv("../data/MDR_and_pLA_Data.csv", stringsAsFactors = T)
# Remove cars with less than 6 cylinders from the mtcars dataset
dat <- dat[dat$trait <= 100, ]

selected_TPC <- c("atkin_2005",  "quadratic_2008","thomas_2017") #this worked without # as a list
selected_species <- c("Aalb","Cpip","Ctar","Cqui")
selected_trait <- "pLA"

model_options <- c("atkin_2005",  "quadratic_2008","thomas_2017", "briere1simplified_1999", "deutsch_2008", "flextpc_2024", "sharpeschoolhigh_1981", "sharpeschoollow_1981")
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
            safetey_margin = get_thermalsafetymargin(model),
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

#function for adding lines to plot - this might need both model and plot passing to it, does that defeat the point 
add_plot_line = function(model) { curve_plot <- curve_plot + geom_line(aes(T.C, model), col = pal[[i]])
i <- i+1
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

  n_iterations <- 200
  
  # loop through each model----
  
  for (tpc in selected_TPC){
    
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
      mod <- nls_multstart(trait~sharpeschoollow_1981(temp = T.C, r_tref, e, el, tl, tref),
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
  }
    

  # create outputs ------
  
  #pal <- viridis(n = length(selected_TPC), option = "plasma")
  pal <- brewer.pal(n = 5,"Set1")
  
  #plot data point only
  plot_bare <- ggplot(filt_data, aes(T.C, trait)) +
    geom_point() +
    theme_bw(base_size = 12) +
    labs(x = 'Temperature (ºC)',
         y = selected_trait,
         title = '')
  
  # plot with curves
  curve_plot <- ggplot(preds, aes(x = T.C, y = trait)) +
    geom_point(aes(T.C, trait), data = filt_data) +
    geom_hline(aes(yintercept = 0), linetype = 2) +
    theme_bw(base_size = 12) +
    labs(x = 'Temperature (ºC)',
         y = selected_trait,
         title = '') 
   
  
  #i was trying to do this as a for loop, can it be done and is there a problem with this way using while 
  #use while for colours 
  
  
  i <- 1
  while(i <= length(selected_TPC)) {
    
    #could this section be a function???
    
    if("atkin_2005" %in% selected_TPC){
      curve_plot <- curve_plot + 
        geom_line(aes(T.C, atkin_2005), col = pal[[i]]) +
        geom_label_repel(aes(T.C, atkin_2005,
          label="Look at this!"
          ),
          # nudge_y = 0.8,
          # nudge_x = 12,
          segment.size = 0.2, 
          segment.colour = 'red')
      i <- i+1
      }
    if("quadratic_2008" %in% selected_TPC) {
      curve_plot <- curve_plot + geom_line(aes(T.C, quadratic_2008), col = pal[[i]])
    
      i <- i+1
      }
    if("thomas_2017" %in% selected_TPC) {
      curve_plot <- curve_plot + geom_line(aes(T.C, thomas_2017), col = pal[[i]])
    
      i <- i+1
      }
    if("briere1simplified_1999" %in% selected_TPC) {
      curve_plot <- curve_plot + geom_line(aes(T.C, briere1simplified_1999), col = pal[[i]])
      i <- i+1
      }
    if("deutsch_2008" %in% selected_TPC) {
      curve_plot <- curve_plot + geom_line(aes(T.C, deutsch_2008), col = pal[[i]])
      i <- i+1
    }
    if("flextpc_2024" %in% selected_TPC) {
      curve_plot <- curve_plot + geom_line(aes(T.C, flextpc_2024), col = pal[[i]])
      i <- i+1
    }
    if("sharpeschoolhigh_1981" %in% selected_TPC) {
      curve_plot <- curve_plot + geom_line(aes(T.C, sharpeschoolhigh_1981), col = pal[[i]])
      i <- i+1
    }
    if("sharpeschoollow_1981" %in% selected_TPC) {
      curve_plot <- curve_plot + geom_line(aes(T.C, sharpeschoollow_1981), col = pal[[i]])
      i <- i+1
    }
  
    
    
  return(list(plot_bare, curve_plot, output_table))
  }
}

saved_output <- shiny_TPC(selected_trait,selected_species)
plot(saved_output[[1]])
plot(saved_output[[2]])
print(saved_output[[3]])






