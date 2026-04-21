#-------------------------------------------------------------------------------
# TITLE: Cum_ProbabilityObservingOutlier.R
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
# the script calculates and visualizes the risk of outliers using the Binomial Theorem.
# This script creates a 3-panel figure showing the probability of obtaining 
# 1, 2, or 3+ outliers as sample size (n) increases.
#-------------------------------------------------------------------------------


# Clear environment
rm(list = ls())



#--------------------FUNCTIONS--------------------

#-------------------------------------------------------------------------------
# FUNCTION: cum_probability_vector
# Purpose: Use the binomial theorem to calculate the probability that a sample 
#    will contain k or more outliers as a function of sample size
# Input:  k: number of extreme values, n_vector: n values to calculate across
# Output: vector p_k_or_more_extremes_vector containing p-values for each value 
#   in n_vector
#-------------------------------------------------------------------------------

cum_probability_vector <- function(k,n_vector) {

  # We want to calculate the probability that a sample will contain k or more outliers 
  # as a function of sample size.  100*prob_k_extremes can be interpreted as
  # how many samples out of 100 would be expected to contain exactly k extreme values.  
  # Extreme values are >2 population SD outside of the population mean
  
  # P(exactly k extreme values) = [n! / (k! * (n-k)!)] * p^k * q^(n-k)
  
  
  q=0.95 #probability of not being >2 population SD outside population mean
  p=1-q #probability of 1 value being extreme

  p_k_or_more_extremes_vector <- numeric(length(n_vector))

  for (i in 1:length(n_vector)) {
  current_n <- n_vector[i]
  cumulative_prob <- 0  # Initialize the cumulative probability for this n
  
  for (j in k:current_n) {  # Loop from k to n (inclusive)
    prob_j_extremes <- choose(current_n, j) * p^j * q^(current_n - j)
    cumulative_prob <- cumulative_prob + prob_j_extremes
  }
  
  p_k_or_more_extremes_vector[i] <- cumulative_prob
  }
  return (p_k_or_more_extremes_vector)
}

#--------------------END FUNCTIONS--------------------


# Setup Parameters & Canvas 
n_vector <- seq(3, 100, by = 1)
par(mfrow = c(1, 3), mar = c(5, 4, 4, 1))

#  Panel 1: k=1 
k <- 1
p1 <- cum_probability_vector(k, n_vector)

# Construct the title
plot_title <- paste("Probability of k=", k," or More Extreme Values") 

plot(n_vector, p1, type = "l", col = "blue", lwd = 2, ylim = c(0, 1),
     xlab = "Sample Size (n)", ylab = "Probability",
     main=(plot_title))

abline(v = 7, col = "red", lty = 2)
# Index 5 is n=7
text(7, 1, paste("n=7:", "p=",round(p1[5], 3)), pos = 4, col = "red", cex = 1.2)
abline(h = 0.5, col = "darkgreen", lty = 3)

# Panel 2: k=2
k <- 2
p2 <- cum_probability_vector(k, n_vector)

# Construct the title 
plot_title <- paste("Probability of k=", k," or More Extreme Values") # Dynamic title

plot(n_vector, p2, type = "l", col = "blue", lwd = 2, ylim = c(0, 1),
     xlab = "Sample Size (n)", ylab = "", 
     main=(plot_title))

abline(v = 7, col = "red", lty = 2)
text(7, 1, paste("n=7:", "p=",round(p2[5], 3)), pos = 4, col = "red", cex = 1.2)
abline(h = 0.5, col = "darkgreen", lty = 3)

# 5. Panel 3: k=3
k <- 3
p3 <- cum_probability_vector(k, n_vector)

# Construct the title
plot_title <- paste("Probability of k=", k," or More Extreme Values") # Dynamic title

plot(n_vector, p3, type = "l", col = "blue", lwd = 2, ylim = c(0, 1),
     xlab = "Sample Size (n)", ylab = "", 
     main=(plot_title))

abline(v = 7, col = "red", lty = 2)
text(7, 1, paste("n=7:", "p=",round(p3[5], 3)), pos = 4, col = "red", cex = 1.2)
abline(h = 0.5, col = "darkgreen", lty = 3)

# Reset layout
par(mfrow = c(1, 1))