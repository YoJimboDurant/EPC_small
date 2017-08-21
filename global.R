## loading source/packages needed application server.R
# package check
packs <- c("shiny", "lattice", "plyr", "EnvStats", "Metrics", "reshape2", "NADA")
packs <- packs[!packs %in% rownames(installed.packages())]
if(length(packs) > 0 ) sapply(packs, install.packages)

# load needed libraries:
library(shiny)
library(lattice)
library(plyr)
library(EnvStats)
library(Metrics)
library(reshape2)
library(NADA)
source("./utilityFunctions.R")
source("./EH_Bayes.R")
