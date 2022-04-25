# Run all scripts ---------------------------------------------------------
source(file = "PancRISK_project/R/01_load.R")
source(file = "PancRISK_project/R/02_clean.R")
source(file = "PancRISK_project/R/03_augment.R")
source(file = "PancRISK_project/R/04_model_i.R")
source(file = "PancRISK_project/R/05_model_ii.R")
source(file = "PancRISK_project/R/05_model_iii.R")
rmarkdown::render(input = "PancRISK_project/Documents/trial_presentation_slide.Rmd", 
     output_file = "PancRISK_project/Documents/trial_presentation_slide.html",
     output_dir = "PancRISK_project/Documents/")