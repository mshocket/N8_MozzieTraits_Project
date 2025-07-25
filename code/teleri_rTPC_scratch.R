
# install.packages("nls.multstart")
# install.packages("remotes")
# remotes::install_github("padpadpadpad/rTPC")

#set working directory
setwd("~/Documents/GitHub/N8_MozzieTraits_Project/code")
dat <- read.csv("../data/MDR_and_pLA_Data.csv", stringsAsFactors = T)

selected_TPC <- c("quadratic_2008","briere1simplified_1999","thomas_2017", "deutsch_2008", "atkin_2005") #this worked without # as a list
selected_species <- c("Aalb","Cpip","Ctar","Cqui")
selected_trait <- "1/MDR"

filt_data <- dat |>
  filter(trait.name == selected_trait)|>
  filter(host.code %in% selected_species) 

preds <- data.frame(T.C = seq(min(filt_data$T.C), max(filt_data$T.C), length.out = 100))  

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

add_row_tpc = function(table_to_add_to, TPC, model) {
  output_table_added <- table_to_add_to |>
    add_row(Model = TPC,
            AIC = AIC(model),
            ct_min = get_ctmin(model),
            ct_max = get_ctmax(model),
            T_opt = get_topt(model),
            breadth = get_breadth(model),
            r_max = get_rmax(model),
            thermal_tolerance = get_thermaltolerance(model),
            safetey_margin = get_thermalsafetymargin(model),
            coeficients = "coef_string(model)")
  return(output_table_added)
}

output_table <- add_row_tpc(output_table, "quadratic_2008", mod)


# get predictions
preds <- broom::augment(mod, newdata = preds)
preds <- rename(preds, quadratic_2008 = .fitted)



start_vals <- get_start_vals(filt_data$T.C, filt_data$trait, model_name = "thomas_2017")

# fit model
mod <- nls_multstart(trait~thomas_2017(temp = T.C, a, b, c, d, e),
                     data = filt_data,
                     iter = 200,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = get_lower_lims(filt_data$T.C, filt_data$trait, model_name = "thomas_2017"),
                     upper = get_upper_lims(filt_data$T.C, filt_data$trait, model_name = "thomas_2017"),
                     supp_errors = 'Y',
                     convergence_count = FALSE)

# add parameter row to table
output_table <- add_row_tpc(output_table, "thomas_2017", mod)


# get predictions
preds <- broom::augment(mod, newdata = preds)
preds <- rename(preds, thomas_2017 = .fitted)


#pivot_longer(preds, names_to = 'model_name', values_to = 'fit')
# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

# MDR ----
# MDR_Data <- read.csv("data/MDR_Data.csv", stringsAsFactors = T)
# names(MDR_Data)
# summary(MDR_Data)

# #take reciprical so all in days?
# for(i in 1:nrow(MDR_Data)){
#   
#   if(MDR_Data$trait.name[i] == "MDR"){
#     MDR_Data$trait[i] <- 1/MDR_Data$trait[i]
#   }
# }

# For Teleri - example of how to structure taking reciprical of some values
# data_test <- d |>
#   mutate(rate = if_else(temp > 28, rate, rate*100))


# chlorella----
# load in data
data("chlorella_tpc")

# keep just a single curve
d <- filter(chlorella_tpc, curve_id == 1)

# show the data
ggplot(d, aes(temp, rate)) +
  geom_point() +
  theme_bw(base_size = 12) +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

get_model_names()

# get start values and fit model
start_vals <- get_start_vals(d$temp, d$rate, model_name = 'atkin_2005')

# fit model
mod <- nls_multstart(rate~atkin_2005(temp = temp, r0, a, b),
                     data = d,
                     iter = 200,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = get_lower_lims(d$temp, d$rate, model_name = 'atkin_2005'),
                     upper = get_upper_lims(d$temp, d$rate, model_name = 'atkin_2005'),
                     supp_errors = 'Y',
                     convergence_count = FALSE)
# look at model fit
summary(mod)

# get predictions
preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
preds <- broom::augment(mod, newdata = preds)

# plot
ggplot(preds) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  theme_bw(base_size = 12)+
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

# multiple models ----
#from https://cran.r-project.org/web/packages/rTPC/vignettes/model_averaging_selection.html
# fit five chosen model formulations in rTPC

# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)
library(ggrepel)

d_fits <- nest(d, data = c(temp, rate)) %>%
  mutate(boatman = map(data, ~nls_multstart(rate~boatman_2017(temp = temp, rmax, tmin, tmax, a,b),
                                            data = .x,
                                            iter = c(4,4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'boatman_2017') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'boatman_2017'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         gaussian = map(data, ~nls_multstart(rate~gaussian_1987(temp = temp, rmax, topt, a),
                                             data = .x,
                                             iter = c(4,4,4),
                                             start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') - 10,
                                             start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'gaussian_1987') + 10,
                                             lower = get_lower_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             upper = get_upper_lims(.x$temp, .x$rate, model_name = 'gaussian_1987'),
                                             supp_errors = 'Y',
                                             convergence_count = FALSE)),
         oneill = map(data, ~nls_multstart(rate~oneill_1972(temp = temp, rmax, ctmax, topt, q10),
                                           data = .x,
                                           iter = c(4,4,4,4),
                                           start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') - 10,
                                           start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'oneill_1972') + 10,
                                           lower = get_lower_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                           upper = get_upper_lims(.x$temp, .x$rate, model_name = 'oneill_1972'),
                                           supp_errors = 'Y',
                                           convergence_count = FALSE)),
         quadratic = map(data, ~nls_multstart(rate~quadratic_2008(temp = temp, a, b, c),
                                              data = .x,
                                              iter = c(4,4,4),
                                              start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') - 0.5,
                                              start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'quadratic_2008') + 0.5,
                                              lower = get_lower_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                              upper = get_upper_lims(.x$temp, .x$rate, model_name = 'quadratic_2008'),
                                              supp_errors = 'Y',
                                              convergence_count = FALSE)),
         rezende = map(data, ~nls_multstart(rate~rezende_2019(temp = temp, q10, a,b,c),
                                            data = .x,
                                            iter = c(4,4,4,4),
                                            start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') - 10,
                                            start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'rezende_2019') + 10,
                                            lower = get_lower_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                            upper = get_upper_lims(.x$temp, .x$rate, model_name = 'rezende_2019'),
                                            supp_errors = 'Y',
                                            convergence_count = FALSE)),
         sharpeschoolhigh = map(data, ~nls_multstart(rate~sharpeschoolhigh_1981(temp = temp, r_tref,e,eh,th, tref = 15),
                                                     data = .x,
                                                     iter = c(4,4,4,4),
                                                     start_lower = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') - 10,
                                                     start_upper = get_start_vals(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981') + 10,
                                                     lower = get_lower_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     upper = get_upper_lims(.x$temp, .x$rate, model_name = 'sharpeschoolhigh_1981'),
                                                     supp_errors = 'Y',
                                                     convergence_count = FALSE)))

# stack models
d_stack <- select(d_fits, -data) %>%
  pivot_longer(., names_to = 'model_name', values_to = 'fit', boatman:sharpeschoolhigh)

# get predictions using augment
newdata <- tibble(temp = seq(min(d$temp), max(d$temp), length.out = 100))
d_preds <- d_stack %>%
  mutate(., preds = map(fit, augment, newdata = newdata)) %>%
  select(-fit) %>%
  unnest(preds)

# take a random point from each model for labelling
d_labs <- filter(d_preds, temp < 30) %>%
  group_by(., model_name) %>%
  sample_n(., 1) %>%
  ungroup()

# plot
ggplot(d_preds, aes(temp, .fitted)) +
  geom_line(aes(col = model_name)) +
  geom_label_repel(aes(temp, .fitted, label = model_name, col = model_name), fill = 'white', nudge_y = 0.8, segment.size = 0.2, segment.colour = 'grey50', d_labs) +
  geom_point(aes(temp, rate), d) +
  theme_bw(base_size = 12) +
  theme(legend.position = 'none') +
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures') +
  geom_hline(aes(yintercept = 0), linetype = 2) +
  scale_color_brewer(type = 'qual', palette = 2)

AIC() # can't work out how to do this for the multiple models version

d_ic <- d_stack %>%
  mutate(., info = map(fit, glance), #this line gets aic and bic from glance, what is.?, don't fully understand what map does 
         AICc =  map_dbl(fit, MuMIn::AICc)) %>%
  select(-fit) %>%
  unnest(info) %>%
  select(model_name, AIC, AICc, BIC)

d_ic
?mutate
# trying different models -----------
get_model_names()
# get start values and fit model
start_vals <- get_start_vals(d$temp, d$rate, model_name = 'quadratic_2008')

# fit model
mod2 <- nls_multstart(rate~quadratic_2008(temp = temp, r0, a, b),
                     data = d,
                     iter = 200,
                     start_lower = start_vals - 10,
                     start_upper = start_vals + 10,
                     lower = get_lower_lims(d$temp, d$rate, model_name = 'quadratic_2008'),
                     upper = get_upper_lims(d$temp, d$rate, model_name = 'quadratic_2008'),
                     supp_errors = 'Y',
                     convergence_count = FALSE)
# look at model fit
summary(mod2)

# get predictions
preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
preds <- broom::augment(mod2, newdata = preds)

# plot
ggplot(preds) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  theme_bw(base_size = 12)+
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')


# try a different model
get_model_names()
# get start values and fit model
start_vals <- get_start_vals(d$temp, d$rate, model_name = 'boatman_2017')

# fit model
mod3 <- nls_multstart(rate~boatman_2017(temp = temp,  rmax, tmin, tmax, a, b),
                      data = d,
                      iter = 200,
                      start_lower = start_vals - 10,
                      start_upper = start_vals + 10,
                      lower = get_lower_lims(d$temp, d$rate, model_name = 'boatman_2017'),
                      upper = get_upper_lims(d$temp, d$rate, model_name = 'boatman_2017'),
                      supp_errors = 'Y',
                      convergence_count = FALSE)
# look at model fit
summary(mod3)

# get predictions
preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
preds <- broom::augment(mod3, newdata = preds)

# plot
ggplot(preds) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  theme_bw(base_size = 12)+
  labs(x = 'Temperature (ºC)',
       y = 'Metabolic rate',
       title = 'Respiration across temperatures')

get_model_names()

?rezende_2019

dat <- read.csv("data/Fecundity_Data.csv", stringsAsFactors = T)
levels(dat$Trait.Name)
summary(dat)

#for loops
vec <- c(3:4,6)
numbers <- c()

counts <- function(vec){
for(val in vec) {
  if(val == 1){ 
    numbers <- paste(numbers, "one")
    }
  else if(val == 2){
    numbers <- paste(numbers, "two")
  }
  else if(val == 3){
    numbers <- paste(numbers, "three")
  }
  else if(val == 4){
    numbers <- paste(numbers, "four")
  }
  else if(val == 5){
    numbers <- paste(numbers, "five")
  }
  else if(val == 6){
    numbers <- paste(numbers, "six")
  }
  else if(val == 7){
    numbers <- paste(numbers, "seven")
  }
}
  return(numbers)
}

counts(vec)

#get params function

get_parameters <- function(model){
  mod_params <- c(AIC(model), get_ctmin(model), get_ctmax(model), get_topt(model), get_breadth(model), get_rmax(model), get_thermaltolerance(model), get_thermalsafetymargin(model), "param list")
}



output_table <- rbind(output_table, mod_params)


# add_row_function <- function(results_table, model, tpc){
#   results_table <- 
#   add_row(results_table, 
#           Model = tpc,
#           AIC = AIC(model),
#           ct_min = get_ctmin(model),
#           ct_max = get_ctmax(model),
#           T_opt = get_topt(model),
#           breadth = get_breadth(model),
#           r_max = get_rmax(model),
#           thermal_tolerance = get_thermaltolerance(model),
#           safetey_margin = get_thermalsafetymargin(model),
#           coeficients = coef_string(model))
# return(output_table)
# }





# Turn this into function - takes coef(mod) and turn it into a single string

coef_string <- function(model){
  vec <- character(length(coef(mod)))
  
  for(i in 1:length(coef(mod))){
    vec[i] <- paste(names(coef(mod))[i], signif(coef(mod)[i], digits = 3), sep = ": ")
  }
  coefs_out <- paste(vec, collapse = ", ")
  return(coefs_out)
}

coef_string(mod)

coefs_out <- paste(vec, collapse = ", ")
coefs_out





#all data for papers ----

a_Data <- read.csv("data/a_Data.csv", stringsAsFactors = T)
AdultSurvival_Data <- read.csv("data/AdultSurvival_Data.csv", stringsAsFactors = T)
Fecundity_Data <- read.csv("data/Fecundity_Data.csv", stringsAsFactors = T)
MDR_Data <- read.csv("data/MDR_Data.csv", stringsAsFactors = T)
PDR_Data <- read.csv("data/PDR_Data.csv", stringsAsFactors = T)
pLA_Data <- read.csv("data/pLA_Data.csv", stringsAsFactors = T)

levels(a_Data$Citation)
levels(AdultSurvival_Data$Citation)
levels(Fecundity_Data$Citation)
levels(MDR_Data$Citation)
levels(PDR_Data$Citation)
levels(pLA_Data$Citation)

cits <- as.vector(c(levels(a_Data$Citation),
                    levels(AdultSurvival_Data$Citation),
                    levels(Fecundity_Data$Citation),
                    levels(MDR_Data$Citation),
                    levels(PDR_Data$Citation),
                    levels(pLA_Data$Citation)))
sum(duplicated(cits))

citations <- unique(cits)


#previous plot code #####
i <- 1
while(i <= length(selected_TPC)) {
  
  #could this section be a function???
  
  if("atkin_2005" %in% selected_TPC){
    curve_plot <- curve_plot + 
      geom_line(aes(T.C, atkin_2005), col = pal[[i]]) 
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
  
}
