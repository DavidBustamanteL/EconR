############################### Testing FE and RE Models EconR ###############################
##############################################################################################

# let us begin as usual by setting our wd
rm(list = ls())
gc()

# English Format and Avoiding scientific notation
Sys.setlocale("LC_TIME", "English")
options(scipen = 999)

# Packages
library(tidyverse)
library(magrittr)
library(estimatr)
library(sandwich)
library(clubSandwich)
library(plm)
library(fixest)

# Installing test pack if required
if (!requireNamespace("EconR", quietly = TRUE)) {
    install.packages(
        "G:/My Drive/5. CURSOS ONLINE/R Practice FULL/Econ Own Package",
        repos = NULL,
        type = "source"
    )
}

library(EconR)




































#### Dev Tool Testing ####
devtools::load_all()
devtools::test(filter = "panel")

devtools::document()

devtools::check()