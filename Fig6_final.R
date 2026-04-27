#-------------------------------------------------------------------------------
# TITLE: Fig6_final.R
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
# DESCRIPTION: This script generates figure 6 from the manuscript and is 
#   annotated to provide a beginner friendly introduction to the R programming language.
#-------------------------------------------------------------------------------


# Clear environment: Prevents variables from previous projects 
# from interfering with your current analysis.
rm(list = ls())

# Tip for Publication: Journals usually require high-resolution TIFFs.
# The code below is a "file device." When active, R redirects the plot 
# from the screen directly to a high-quality file in your working directory.
# Remove the "#" in the line below and change the Figure file name as needed
# You will also need to remove the '#' from the last line of code in this file

# tiff("Figure6_Final.tif", width = 8, height = 4.5, units = "in", res = 300, compression = "lzw")

#--------------------Figure 6:  Visualizing Confidence Limits--------------------------



#--------------------FUNCTIONS--------------------

#-------------------------------------------------------------------------------
# FUNCTION: plot_caterpillar
# Purpose: Simulates multiple experiments to visualize Confidence Interval (CI) 
# coverage. It demonstrates that a CI is a property of the sampling procedure, 
# not a single experiment.
#
# Inputs: 
#   CI         : The confidence level (e.g., 0.95 for 95%).
#   n          : Sample size for each individual simulation.
#   mu_true    : The known population mean (the "truth" we are trying to capture).
#   sigma      : The population standard deviation.
#   n_sims     : The number of simulated experiments to stack.
#   main_title : Title for the plot.
#-------------------------------------------------------------------------------

plot_caterpillar <- function(CI,n,mu_true,sigma,n_sims,main_title) {
  set.seed(25)# Ensures the "random" simulation is reproducible for the paper.
  
  # Statistical Setup
  # We calculate the upper quantile (e.g., 0.975 for a 95% CI) to find the 
  # critical t-value. This accounts for the uncertainty in small sample sizes (df = n-1).
  quant <- (((1 - CI) / 2) + CI)
  t_crit <- qt(quant, df = n - 1)
  
  # The Simulation Loop
  # Pre-allocate vectors to store results (memory efficient R)
  means <- numeric(n_sims)
  lowers <- numeric(n_sims)
  uppers <- numeric(n_sims)

  for(i in 1:n_sims) {
    # Generate random sample from the theoretical population
    samp <- rnorm(n, mean = mu_true, sd = sigma)
    m <- mean(samp)
    se <- sd(samp) / sqrt(n)
    
    # Define the bounds of the interval using the t-distribution
    means[i] <- m
    lowers[i] <- m - t_crit * se
    uppers[i] <- m + t_crit * se
  }

  # Determine "Capture" status (Color Coding)
  # A 'capture' occurs only if the true mean (mu_true) lies between the lower 
  # and upper bounds of the sample's CI.
  # Blue (4) for capture, Orange-ish Red (2) for miss
  colors <- ifelse(lowers <= mu_true & uppers >= mu_true, "steelblue", "orangered")

  #Create the Visualization
  # Initialize an empty plot. The Y-axis represents each unique 'Simulation ID'.
  plot(NULL, xlim = c(200, 1400), ylim = c(1, n_sims), 
       xlab = expression("Tumor Volume (mm"^3*")"), ylab = "Simulation ID",
       main = main_title)

  # Add a vertical reference line for the ground truth population mean
  abline(v = mu_true, col = "black", lty = 2, lwd = 2)

  # Use segments() to draw the horizontal CI lines for every simulation
  # x0/x1 define the CI width; y0/y1 stack them vertically from 1 to n_sims.
  segments(x0 = lowers, y0 = 1:n_sims, x1 = uppers, y1 = 1:n_sims, col = colors)

  # Overlay the calculated mean for each experiment
  points(means, 1:n_sims, pch = 20, cex = 0.6, col = colors)

  # Add a legend
  legend("topright", 
         legend = c("Captured", "Missed"), 
         text.col = c("steelblue", "orangered"), 
         bty = "n")
}

#--------------------End FUNCTIONS--------------------


# Set up for a 2-panel figure to compare different CIs
par(mfrow = c(1, 2))

#--------------------Figure 6A--------------------------------------------------
#Define Population and Simulation Parameters for 95% CIs

n       <- 7    # Small sample size, common in pilot studies.
mu_true <- 800  # The "Ground Truth" we want our CIs to capture.
sigma   <- 200  # The biological variability (Standard Deviation).
n_sims  <- 100  # We repeat the entire experiment 100 times.
CI      <- 0.95 # We want 95% confidence.

main_title <- 'A) 95% CI Simulation (n=7)'

# Execute the Simulation and Plot
# This calls our custom function which:
#  a) Generates random data
#  b) Calculates t-based CIs
#  c) Plots the 'Caterpillar' stack
plot_caterpillar(CI, n, mu_true, sigma, n_sims, main_title)


#--------------------Figure 6b--------------------------------------------------
#Define Population and Simulation Parameters for 68% (SEM) CIs

#Same logic as fig 6A
#Setup Parameters for 68% CI (SEM)
n <- 7
mu_true <- 800
sigma <- 200
n_sims <- 100
CI<-0.68
main_title<-'B) 68% CI (SEM) Simulation (n=7)'

plot_caterpillar(CI,n,mu_true,sigma,n_sims,main_title)


# Final step for high-quality export
# If you opened the tiff() device at the start of the script, 
# you MUST run dev.off() by removing the '#' below to finish writing the file and close it.
# dev.off()