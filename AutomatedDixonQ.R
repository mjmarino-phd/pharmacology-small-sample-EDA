#-------------------------------------------------------------------------------
# TITLE: AutomatedDixonQ.R
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
# DESCRIPTION: This script provides an automated statistical screen that mimics the 
# workflow outlined in the manuscript. It first looks at the data to identify 
# outliers that would be flagged in a Tukey boxplot, then runs the appropriate 
# Dixon's Q-test and outputs a table of identified values and associated p-values
# It is designed for researchers who need to identify 
# statistical outliers in small datasets (n=3 to n=30).
# The script is annotated to provide a beginner friendly introduction to the R 
# programming language, however this is a more advanced script
#
# This uses the example file CaffeineData.csv
# See section titled "-Using your own data-" to apply to any CSV file
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

#--------------------FUNCTIONS--------------------

#-------------------------------------------------------------------------------
# FUNCTION: generate_extreme_outliers
# Purpose: Scans a dataset to identify the most extreme high and low values 
# for each group, based on the standard Tukey Boxplot criteria.
#
# This is useful for identifying the specific "candidates" that might 
# require formal outlier testing (like Dixon's Q).
#
# Input:  data - A data frame or list where each column is a group.
# Output: A summary data frame listing the group name and its most extreme 
#         upper and lower outliers (if any).
#-------------------------------------------------------------------------------

generate_extreme_outliers <- function(data) {
  
  # Extract Boxplot Statistics
  # Setting 'plot = FALSE' allows us to use the math of a boxplot (fences, 
  # whiskers, and outliers) without actually drawing a figure.
  boxplot_stats <- boxplot(data, plot = FALSE)
  
  # Prepare the Results Table
  # Initialize an empty data frame to store our findings as we loop.
  extreme_outliers <- data.frame(
    group_name = character(),
    lower_outlier = numeric(),
    upper_outlier = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Analyze Each Group
  # 'boxplot_stats$names' contains the list of groups (e.g., Vehicle, Caffeine).
  for (i in 1:length(boxplot_stats$names)) {
    group_name <- boxplot_stats$names[i]
    
    # R stores all outliers in one long vector ($out).
    # We filter that vector to find only the points belonging to the current group index (i).
    group_outliers <- boxplot_stats$out[boxplot_stats$group == i]
    
    # We use the median ($stats[3,i]) as the "anchor" to determine if a point 
    # is a high outlier or a low outlier.
    group_median <- boxplot_stats$stats[3, i]
    
    lower_outlier <- NA
    upper_outlier <- NA
    
    # Filter for Extremes
    if (length(group_outliers) > 0) {
      # Lower outliers: Values located below the median. 
      # We take the 'min' because it is the point furthest away on the low end.
      lower_outlier_values <- group_outliers[group_outliers < group_median]
      if (length(lower_outlier_values) > 0) {
        lower_outlier <- min(lower_outlier_values)
      }
      
      # Upper outliers: Values located above the median.
      # We take the 'max' because it is the point furthest away on the high end.
      upper_outlier_values <- group_outliers[group_outliers > group_median]
      if (length(upper_outlier_values) > 0) {
        upper_outlier <- max(upper_outlier_values)
      }
    }
    
    # Compile Results
    # We 'rbind' (row-bind) the group's extremes to our master table.
    extreme_outliers <- rbind(extreme_outliers, data.frame(
      group_name = group_name,
      lower_outlier = lower_outlier,
      upper_outlier = upper_outlier
    ))
  }
  
  return(extreme_outliers)
}

#-------------------------------------------------------------------------------
# FUNCTION: generate_dixon_truth
# Purpose: Creates a "Truth Table" to guide the use of the 'opposite' parameter 
# in Dixon's Q-test. 
#
# Because Dixon's test automatically targets the value furthest from the mean,
# this function calculates whether your specific outlier (high or low) requires
# the 'opposite = TRUE' flag to be correctly evaluated.
#
# Inputs: 
#   data             : The full data frame (log-transformed).
#   extreme_outliers : The output table from the 'generate_extreme_outliers' function.
# Output: A summary data frame listing the 'opposite'(TRUE/FALSE) flags for identified potential 
# outliers.
#-------------------------------------------------------------------------------

generate_dixon_truth <- function(data, extreme_outliers) {
  
  # Initialize the Output Table
  # Create a table that will state if 'opposite' should 
  # be TRUE or FALSE for each potential outlier.
  dixon_truth <- data.frame(
    group_name = extreme_outliers$group_name,
    test_low   = NA,  # Recommended 'opposite' setting for the lower outlier
    test_high  = NA,  # Recommended 'opposite' setting for the upper outlier
    stringsAsFactors = FALSE
  )
  
  # Compare Distances from the Mean
  for (i in 1:nrow(extreme_outliers)) {
    group <- extreme_outliers$group_name[i]
    current_data <- na.omit(data[[group]])
    grp_mean <- mean(current_data)
    
    # --- Logic for Lower Outlier ---
    if (!is.na(extreme_outliers$lower_outlier[i])) {
      val <- extreme_outliers$lower_outlier[i]
      
      # Is the Max value further from the mean than our Low outlier?
      # If TRUE, then the Low outlier is NOT the most extreme value,
      # and we must set 'opposite = TRUE' to force Dixon to look at it.
      dixon_truth$test_low[i] <- abs(grp_mean - max(current_data)) > abs(grp_mean - val)
    }
    
    # --- Logic for Upper Outlier ---
    if (!is.na(extreme_outliers$upper_outlier[i])) {
      val <- extreme_outliers$upper_outlier[i]
      
      # Is the Min value further from the mean than our High outlier?
      # If TRUE, then the High outlier is NOT the most extreme value,
      # and we must set 'opposite = TRUE' to force Dixon to look at it.
      dixon_truth$test_high[i] <- abs(grp_mean - min(current_data)) > abs(grp_mean - val)
    }
  }
  
  return(dixon_truth)
}

#-------------------------------------------------------------------------------
# FUNCTION: run_dixon
# Purpose: Automates the execution of Dixon's Q-test for every group in a dataset.
# It uses the Truth Table generated by the generate_dixon_truth function 
# to ensure each test is targeted at the correct value using the proper 'opposite' flag.
#
# Inputs: 
#   extreme_outliers : Table of candidate outlier values (from generate_extreme_outliers).
#   dixon_truth      : Logic map for the 'opposite' flag (from generate_dixon_truth).
#   dixon_type       : The Dixon ratio type (use 0 for automatic selection).
#   data             : The full dataset being analyzed.
#
# Output: A list containing two data frames:
#         1. p_values: The statistical probability for each candidate.
#         2. q_statistics: The calculated Qexp value for each candidate.
#-------------------------------------------------------------------------------

run_dixon <- function(extreme_outliers, dixon_truth, dixon_type, data) {
  
  # 1. Initialize the p_value output table as a blank slate.
  p_extreme_outliers <- data.frame(
    group_name = extreme_outliers$group_name, # Pre-fill names to keep tables identical
    lower_outlier = as.numeric(rep(NA, nrow(extreme_outliers))),
    upper_outlier = as.numeric(rep(NA, nrow(extreme_outliers))),
    stringsAsFactors = FALSE
  )
  
  # 2. Initialize the Q-Statistic (Qexp) Results Table by copying the P table structure
  q_results <- p_extreme_outliers 
  
  for (i in 1:nrow(extreme_outliers)) {
    
    group_name <- extreme_outliers[i, "group_name"]
    
    # --- 1. Test for a LOWER outlier ---
    if (!is.na(extreme_outliers[i, "lower_outlier"])) {
      is_opp_low <- dixon_truth[i, "test_low"]
      
      test_result_lower <- dixon.test(data[[group_name]], 
                                      type = dixon_type, 
                                      opposite = is_opp_low, 
                                      two.sided = TRUE)
      
      # Assign both P and Q values while the test object is "active"
      p_extreme_outliers[i, "lower_outlier"] <- test_result_lower$p.value
      q_results[i, "lower_outlier"]          <- test_result_lower$statistic 
    }
    
    # --- 2. Test for an UPPER outlier ---
    if (!is.na(extreme_outliers[i, "upper_outlier"])) {
      is_opp_high <- dixon_truth[i, "test_high"]
      
      test_result_upper <- dixon.test(data[[group_name]], 
                                      type = dixon_type, 
                                      opposite = is_opp_high, 
                                      two.sided = TRUE)
      
      # Assign both P and Q values
      p_extreme_outliers[i, "upper_outlier"] <- test_result_upper$p.value
      q_results[i, "upper_outlier"]          <- test_result_upper$statistic
    }
  }
  
  # Return both tables as a named list for easy access
  return(list(p_values = p_extreme_outliers, q_statistics = q_results))
}

#--------------------END FUNCTIONS--------------------



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

#Run the test
#For a full explanation of the 
#basic dixon.test() function see SimpleDixonQ.R

dixon_type=0
extreme_outliers<-generate_extreme_outliers(log_data)
dixon_truth<-generate_dixon_truth(log_data,extreme_outliers)
results<-run_dixon(extreme_outliers,dixon_truth,dixon_type,log_data)

#print outputs to terminal
print(extreme_outliers) #prints log-transformed outlier
print(dixon_truth) # Prints 'opposite=' flag
print(results$q_statistics) # prints Qexp
print(results$p_values) # prints p-values





