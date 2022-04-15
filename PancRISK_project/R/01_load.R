# Load libraries ----------------------------------------------------------
library("tidyverse")


# Load data ---------------------------------------------------------------
my_data_raw <- read_csv(file = "PancRISK_project/data/_raw/PancRISK_raw.csv", 
                        na = c("", "NA"))
# Optional to drop NAs of plasma_CA19_9, since it is important to compare with urinary values.


# Write data --------------------------------------------------------------
write_tsv(x = my_data_raw,
          file = "PancRISK_project/data/01_my_data.tsv")