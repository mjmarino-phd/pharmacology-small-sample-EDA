#-------------------------------------------------------------------------------
# TITLE: SimpleDixonQ.R
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
# DESCRIPTION: This script provides a template for running Dixon's Q-test on 
# pharmacological data. It is designed for researchers who need to identify 
# statistical outliers in small datasets (n=3 to n=30).  The script is 
# annotated to provide a beginner friendly introduction to the R programming language.
#-------------------------------------------------------------------------------


# Clear environment: Prevents variables from previous projects 
# from interfering with your current analysis.
rm(list = ls())

# This block checks for the 'outliers' package. If it is missing, it installs it.
if (!require("outliers", quietly = TRUE)) {
  install.packages("outliers", repos = "https://cloud.r-project.org/")
  library(outliers)
}

# Load necessary libraries
library(outliers)

#-------Using your own data---------------------------
#This script is a simple example of running Dixon's Q
#to run on your own CSV file you would need to change the filename variable 
#below and choose the appropriate group name in the function call to dixon.test()
# (e.g. dixon.test(log_data$MyGroupName,type=0,opposite = FALSE,two.sided = TRUE)


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


# The R Dixon's Q test function takes the form:
# dixon.test(x, type = 0, opposite = FALSE, two.sided = TRUE)
# Where:
## x is the sample data set to be evaluated (e.g. caffeine data or vehicle data)

## type is an integer value that sets the Ratio value (R10, R14, etc) 
###(default =0)

## opposite is a flag that determines if the test will evaluate the most extreme 
### value relative to the mean (opposite=FALSE) or the value at the other end of 
### the sample distribution (opposite=TRUE)
###(default=FALSE)

## two.sided determines one vs two sided p-value stringency
### (default = TRUE)


#This will run the Dixon's Q and output the results in the console. 
#  (output explained below)

dixon.test(log_data$Caffeine,type=0,opposite = FALSE,two.sided = TRUE)

# log_data$Caffeine chooses the Caffeine vector from the dataframe

# type 10 is appropriate for small (n=3 - n=7) datasets

# type = 0 (chosen here as the smart default)
## This automatically selects the most powerful Dixon ratio for your sample size.
## For n < 30, it ensures the test uses the statistically appropriate formula 
## (e.g., it may switch from type 10 to type 11 or 21 as your 'n' increases).

# two.sided is TRUE for statistical rigor (no reason to expect high vs low outliers)

#Here is the expected output:

# Dixon test for outliers
# 
# data:  log_data$Caffeine
# Q = 0.6, p-value = 0.03341
# alternative hypothesis: highest value 3.1 is an outlier

# Using the provided CaffeineData.csv sample data file, The function provides the 
# value for Qexp (0.6), and the exact p-value (0.03341) if p<0.05, we reject the 
# null hypothesis and can exclude the value identified in the alternative hypothesis


# Beginners tip
# in RStudio if you highlight any function or keyword and press F1 
# you will be taken to the corresponding help screen




