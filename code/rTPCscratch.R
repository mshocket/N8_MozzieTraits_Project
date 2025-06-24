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

# get predictions
preds <- data.frame(temp = seq(min(d$temp), max(d$temp), length.out = 100))
preds <- broom::augment(mod, newdata = preds)

# plot
ggplot(preds) +
  geom_point(aes(temp, rate), d) +
  geom_line(aes(temp, .fitted), col = 'blue') +
  theme_bw()

get_model_names()
