# install.packages("rTPC")
# install.packages("nls.multstart")
# install.packages("remotes")
# remotes::install_github("padpadpadpad/rTPC")

# load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)

# write function to label ggplot2 panels
label_facets_num <- function(string){
  len <- length(string)
  string = paste('(', 1:len, ') ', string, sep = '')
  return(string)
}

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

#Get TPC characteristics
get_topt(mod)
get_ctmax(mod)
get_ctmin(mod)
get_breadth(mod, level = 0.00009)
get_rmax(mod)
get_thermaltolerance(mod)

#Get parameters and names
coef(mod)
coef(mod)[1]
coef(mod)[1]+coef(mod)[2]
as.numeric(coef(mod)[1])
names(coef(mod))[1]

# concatenate text - cat only prints to console, does not save as a variable
printitem <- paste("ant", "bee", "cow", sep = "")
paste(paste(names(coef(mod))[1], signif(coef(mod)[1], digits = 3), sep = ": "),
      paste(names(coef(mod))[2], signif(coef(mod)[2], digits = 3), sep = ": "),
      paste(names(coef(mod))[3], signif(coef(mod)[3], digits = 3), sep = ": "), sep = ", ")
?paste

# Turn this into function - takes coef(mod) and turn it into a single string
vec <- character(length(coef(mod)))

for(i in 1:length(coef(mod))){
 vec[i] <- paste(names(coef(mod))[i], signif(coef(mod)[i], digits = 3), sep = ": ")
}
  
coefs_out <- paste(vec, collapse = ", ")


str(coefs_out)

# Get AIC score
AIC(mod)
BIC(mod)
glance(mod)

# get predictions
preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
preds <- broom::augment(mod, newdata = preds)

# plot
ggplot(preds) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  theme_bw()

get_model_names()

1 - as.numeric(c(TRUE, FALSE, FALSE, TRUE))

for(i in 1:nrow(d)){
  
  if(d$temp[i] == 16){
    d$rate[i] <- 1/d$rate[i]
  }
}

# For Teleri - example of how to structure taking reciprical of some values
data_test <- d |>
  mutate(rate = if_else(temp > 28, rate, rate*100))

fake_data <- data.frame(name = c("ant", "ant", "ant", "bee", "bee", "cow", "cow", "cow", "cow"),
                        value = c(seq(1,9,1)))

filt_data <- fake_data |>
  filter(name == "ant" | name == "cow")

filt_data <- fake_data |>
  filter(name != c("bee"))

filt_data <- fake_data |>
  filter(name %in% c("ant", "bee"))


##### Writing your own functions

# Define a function
VectorTPCshiny = function(species, trait, TPC){
  
  VTS_out <- paste(species, trait, TPC)
  
  return(VTS_out)
  
}

# Give inputs for a function
species <- "Aaeg"
trait <- "1/MDR"
TPC <- "atkins2005"

# Run the function
saved_output <- VectorTPCshiny(species, trait, TPC)
