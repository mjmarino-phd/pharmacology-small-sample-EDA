#-------------------------------------------------------------------------------
# TITLE: EffectSIzeCalculator.R
# AUTHOR: Michael J. Marino, PhD
# DATE: 2026-04-20
# MANUSCRIPT: "A Robust Workflow for Exploratory Data Analysis and Outlier 
#   Management in Small-Sample Pharmacological Datasets" (Submitted)
#
# REPOSITORY: https://github.com/mjmarino-phd/pharmacology-small-sample-EDA
# 
# LICENSE: MIT License (c) 2026 Michael J. Marino
#   Permission is granted to use, modify, and distribute this code for any 
#   purpose, provided that the manuscript above is cited in any resulting 
#   publications or derivative works.
#
# DESCRIPTION: This script imports a CSV file formatted like the provided sample 
# data set CaffeineData.csv and calculates effect size estimates including Cohen's d, 
# Hedges' g, and Cliff's Delta.  
# The script is annotated to provide a beginner friendly introduction to the R 
# programming language.
#-------------------------------------------------------------------------------



# Clear environment: Prevents variables from previous projects 
# from interfering with your current analysis.
rm(list = ls())

# This block checks for the 'effsize' package. If it is missing, it installs it.
if (!require("effsize", quietly = TRUE)) {
  install.packages("effsize", repos = "https://cloud.r-project.org/")
  library(effsize)
}

# install.packages
library(effsize)

# Data Import
# We define the filename as a variable to make the script easily adaptable 
# for different experimental runs.
filename <- "CaffeineData.csv"  # file must be in the current working directory

# Read the raw data. 
data <- read.csv(filename)

# Log-Transformation
log_data <- log10(data)

#The step below rounds the log transformed data to 2 decimal places.  
log_data<-round(log_data,2)

sample1 <- log_data[, 1]
sample2 <- log_data[, 2]


#----------Cohen's d----------
d_res <- cohen.d(sample1, sample2)

#----------Hedges' g----------
g_res <- cohen.d(sample1, sample2, hedges.correction = TRUE)

#----------Cliff's Delta----------
cliff_res <- cliff.delta(sample1, sample2)


# 4. Print results for your table
print('Effect Size Estimates')
print(d_res)
print(g_res)
print(cliff_res)

# Create the Human-Friendly Summary Table ---
# We extract the estimates, format the CIs, and grab the magnitude.
effect_summary <- data.frame(
  Metric    = c("Cohen's d", "Hedges' g", "Cliff's Delta"),
  
  Estimate  = round(c(d_res$estimate, 
                      g_res$estimate, 
                      cliff_res$estimate), 3),
  
  '95% CI'  = c(
    paste0("[", round(d_res$conf.int[1], 3), ", ", round(d_res$conf.int[2], 3), "]"),
    paste0("[", round(g_res$conf.int[1], 3), ", ", round(g_res$conf.int[2], 3), "]"),
    paste0("[", round(cliff_res$conf.int[1], 3), ", ", round(cliff_res$conf.int[2], 3), "]")
  ),
  
  Magnitude = c(as.character(d_res$magnitude), 
                as.character(g_res$magnitude), 
                as.character(cliff_res$magnitude))
)

# --- 4. Final Output ---
cat("\n--- Statistical Effect Size Summary ---\n")
print(effect_summary, row.names = FALSE)