# Load libraries ----------------------------------------------------------
library("tidyverse")
library("patchwork")
library("scales")
library("nnet")
library("pROC")
library("caret")
library("broom")
library("knitr")
library("kableExtra")
library("cutoff")
library("cowplot")

# Load functions ----------------------------------------------------------
source(file = "PancRISK_project/R/99_proj_functions.R")
# Run all scripts ---------------------------------------------------------
source(file = "PancRISK_project/R/01_load.R")
source(file = "PancRISK_project/R/02_clean.R")
source(file = "PancRISK_project/R/03_data_visualization.R")
source(file = "PancRISK_project/R/04_augment.R")
source(file = "PancRISK_project/R/05_model_i.R")
source(file = "PancRISK_project/R/06_model_ii.R")
source(file = "PancRISK_project/R/07_model_iii.R")
rmarkdown::render(input = "PancRISK_project/Documents/trial_presentation_slide.Rmd", 
     output_file = "PancRISK_project/Documents/trial_presentation_slide.html",
     output_dir = "PancRISK_project/Documents/")