#load packages
library(rTPC)
library(nls.multstart)
library(broom)
library(tidyverse)


models <- c("quadratic_2008","briere1simplified_1999","thomas_2017", "deutsch_2008", "atkin_2005") #this worked without # as a list
species <- c("Aalb","Cpip","Ctar","Cqui")
trait <- "1/MDR"

#load data
dat <- read.csv("../data/MDR_and_pLA_Data.csv", stringsAsFactors = T)

