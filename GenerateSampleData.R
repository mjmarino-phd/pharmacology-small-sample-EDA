#-------------------------------------------------------------------------------
# TITLE: GenerateSampleData.R
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
# DESCRIPTION: This is a more advanced script provided for those who might be interested 
# This script was designed to randomly draw from the two populations to find a
# pedagogical data set that met the following criteria:
## 1) Samples not significantly different by t-test
## 2) Log10-transformed samples not significantly different by t-test
## 3) Contained at least 1 statistical outlier assessed by Dixon's Q
## 4) When outlier is removed, the Log10-transformed samples are significantly 
#   different by t-test
#-------------------------------------------------------------------------------



# Clear environment: Prevents variables from previous projects 
# from interfering with your current analysis.
rm(list = ls())

# Note: A specific seed was not used for the original manuscript dataset. 
# As written, each execution will draw fresh samples from the populations. 
# Uncomment the line below to ensure repeatable results for your own testing.
# set.seed(42)

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
  

  #for Sample generation script, just rerun p-value data frame
  return(p_extreme_outliers)
}

#--------------------END FUNCTIONS--------------------



# Set up populations
# Define the Vehicle and Caffeine Populations
# In rlnorm, meanlog and sdlog are the mean and SD of the data ON THE LOG SCALE

# Historical targets
m <- 800  # Vehicle population mean
s <- 200  # Vehicle population standard deviation
caf_m <-600 # Caffeine population mean
caf_s <-200 # Caffein population standard deviation

# Math to convert Arithmetic Mean/SD to Log-Normal meanlog/sdlog
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))

caf_location <- log(caf_m^2 / sqrt(caf_s^2 + caf_m^2))
caf_shape <- sqrt(log(1 + (caf_s^2 / caf_m^2)))

# Generate the populations
pop_veh <- rlnorm(n = 10000, meanlog = location, sdlog = shape)
pop_caf <- rlnorm(n = 10000, meanlog = caf_location, sdlog = caf_shape)





#set up a loop to run a maximum of loopMax times but halt if criteria are met
data_loop<-0
loopMax<-1000
sample_size <- 7  # Number of sampels to draw from each population
ceiling_limit <- 2000  # Upper practical limit on tumor volume
while (data_loop<loopMax){
  counter<-0  #This keeps track of the criteria
  
  # Sample Vehicle population
  veh_sample <- c()
  while(length(veh_sample) < sample_size) {
    # Pull a fresh candidate from the population
    candidate <- sample(pop_veh, 1)
    # Only retain it if it's within the limit
    if(candidate <= ceiling_limit) {
      veh_sample <- c(veh_sample, candidate)
    }
  }
  
  # Sample Caffeine population
  caf_sample <- c()
  while(length(caf_sample) < sample_size) {
    candidate <- sample(pop_caf, 1)
    if(candidate <= ceiling_limit) {
      caf_sample <- c(caf_sample, candidate)
    }
  }
  
  #Combine samples in a data frame
  data <- data.frame(
    Vehicle = veh_sample,
    Caffeine = caf_sample
  )
  
  #check for linear significant  difference
  
  #t-test
  t_results <- t.test(veh_sample, caf_sample, alternative = "two.sided", var.equal = TRUE)
  p_val<-t_results$p.value
  
  if (!is.na(p_val) && p_val > 0.05) {
    counter <- counter + 1}  # Increment counter if not significant
  
  
  #check for log significant  difference
  
  #Log Transform
  log_data<-log10(data)
  
  #t-test
  t_log_results <- t.test(log_data$Vehicle, log_data$Caffeine, alternative = "two.sided", var.equal = TRUE)
  p_log_val<-t_log_results$p.value
  
  if (!is.na(p_log_val) && p_log_val > 0.05) {
    counter <- counter + 1}  # Increment counter if not significant
  
  #Check for outliers
  
  dixon_type=0
  extreme_outliers<-generate_extreme_outliers(log_data)
  dixon_truth<-generate_dixon_truth(log_data,extreme_outliers)
  p_extreme_outliers<-run_dixon(extreme_outliers,dixon_truth,dixon_type,log_data)
  
  
  #Based on the dixon's p-values remove any potential outliers
  outlier_check<-0
  
  if (any(p_extreme_outliers[, c("lower_outlier", "upper_outlier")] < 0.05, na.rm = TRUE)) {
    
    counter <- counter + 1
    outlier_check<-1}
  
  if (outlier_check==1){
    # Create a clean copy so we don't overwrite raw data
    log_data_cleaned <- log_data
    
    # Iterate through the rows of  p-value table
    for (i in 1:nrow(p_extreme_outliers)) {
      
      # Get the group name (e.g., "Vehicle")
      grp <- p_extreme_outliers$group_name[i]
      
      # Check LOWER Outlier: Is the p-value significant?
      if (!is.na(p_extreme_outliers$lower_outlier[i]) && p_extreme_outliers$lower_outlier[i] < 0.05) {
        # Find the smallest value in the actual data for this group
        val_to_remove <- min(log_data_cleaned[[grp]], na.rm = TRUE)
        # Replace it with NA
        log_data_cleaned[[grp]][log_data_cleaned[[grp]] == val_to_remove] <- NA
        message("Significant Lower Outlier removed from: ", grp)
      }
      
      # Check UPPER Outlier: Is the p-value significant?
      if (!is.na(p_extreme_outliers$upper_outlier[i]) && p_extreme_outliers$upper_outlier[i] < 0.05) {
        # Find the largest value in the actual data for this group
        val_to_remove <- max(log_data_cleaned[[grp]], na.rm = TRUE)
        # Replace it with NA
        log_data_cleaned[[grp]][log_data_cleaned[[grp]] == val_to_remove] <- NA
        message("Significant Upper Outlier removed from: ", grp)
      }
    }
    
    
    #t-test
    t_log_cleaned_results <- t.test(log_data_cleaned$Vehicle, log_data_cleaned$Caffeine, alternative = "two.sided", var.equal = TRUE)
    p_log_cleaned_val<-t_log_cleaned_results$p.value
    
    if (!is.na(p_log_cleaned_val) && p_log_cleaned_val < 0.05) {
      counter <- counter + 1}
  }
  # If all 4 criteria are met, we are done!
  if (counter == 4) {
    # Store the successful data so we can see it after the loop ends
    final_data <- data 
    break # This exits the while loop immediately
  }
  
  # Increment the loop counter
  data_loop <- data_loop + 1
}

# Print the results
if (counter == 4) {
  message("Dataset found on iteration: ", data_loop)
  print(final_data)
} else {
  message("Searched ", loopMax, " iterations but did not find a dataset meeting all 4 criteria.")
}

