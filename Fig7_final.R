#-------------------------------------------------------------------------------
# TITLE: Fig7_final.R
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
# DESCRIPTION: This script generates figure 7 from the manuscript and is 
#   annotated to provide a beginner friendly introduction to the R programming language.
#-------------------------------------------------------------------------------


# Clear environment: Prevents variables from previous projects 
# from interfering with your current analysis.
# rm(list = ls())

# Tip for Publication: Journals usually require high-resolution TIFFs.
# The code below is a "file device." When active, R redirects the plot 
# from the screen directly to a high-quality file in your working directory.
# Remove the "#" in the line below and change the Figure file name as needed
# You will also need to remove the '#' from the last line of code in this file

 tiff("Figure7v2_Final.tif", width = 8, height = 4.5, units = "in", res = 300, compression = "lzw")

#--------------------Figure 7: Confidence Limits for our simulated study--------


#-------------------------------------------------------------------------------
# FUNCTION: plot_caterpillar_compare
# Purpose: Simulates a two-arm study (Vehicle vs. Caffeine) to visualize 
# how Confidence Intervals overlap across multiple experimental iterations.
# This demonstrates the variability of 'significance' in small-n studies.
# See figure 6_Final.R for more detail on the caterpillar plot
#
# Inputs: 
#   CI         : Confidence level (e.g., 0.95).
#   n          : Sample size per group.
#   mu_veh     : True population mean for the Vehicle group.
#   mu_caf     : True population mean for the Caffeine group.
#   sigma      : Population standard deviation (assumed equal for both).
#   n_sims     : Number of simulated "head-to-head" experiments.
#   main_title : Title for the plot.
#-------------------------------------------------------------------------------

plot_caterpillar_compare <- function(CI,n,mu_veh,mu_caf,sigma,n_sims,main_title) {
  set.seed(123) # New seed for the comparative simulation
  
  # Calculate t-distribution critical value
  quant <- ((1 - CI) / 2) + CI
  t_crit <- qt(quant, df = n - 1)

  # Initialize plot
  # xlim expanded to 1400 to accommodate both distributions.
  # las = 1: Rotates y-axis labels to be horizontal (easier to read).
  plot(NULL, xlim = c(0, 1400), ylim = c(1, n_sims), 
       xlab = expression("Tumor Volume (mm"^3*")"), ylab = "Simulation ID",
       main = main_title,
       las = 1)
  
  # Add Population "Ground Truth" Lines
  # These represent the 'True' biological effect we are trying to detect.

  # Add Population Means as reference lines
  abline(v = mu_veh, col = "steelblue", lty = 2, lwd = 1.5)
  abline(v = mu_caf, col = "darkgreen", lty = 2, lwd = 1.5)
  
  # The Dual-Simulation Loop
  # For each 'experiment', we pull samples for both groups and plot their CIs.
  for(i in 1:n_sims) {
    # --- Vehicle Group ---
    s_v <- rnorm(n, mean = mu_veh, sd = sigma)
    m_v <- mean(s_v)
    se_v <- sd(s_v) / sqrt(n)
    segments(m_v - t_crit * se_v, i, m_v + t_crit * se_v, i, 
             col = rgb(70/255, 130/255, 180/255, alpha = 0.5)) # Transparent SteelBlue
    
    # --- Caffeine Group ---
    s_c <- rnorm(n, mean = mu_caf, sd = sigma)
    m_c <- mean(s_c)
    se_c <- sd(s_c) / sqrt(n)
    segments(m_c - t_crit * se_c, i, m_c + t_crit * se_c, i, 
             col = rgb(0/255, 100/255, 0/255, alpha = 0.5))    # Transparent DarkGreen
  }
  

  # Add a legend
  legend("topright", 
         legend = c("Vehicle", "Caffeine"), 
         text.col = c("steelblue", "darkgreen"), 
         bty = "n")
}




#Set up for 2 panel figure
par(mfrow = c(1, 2))


#--------------------Figure 7A--------------------------------------------------
# This simulation demonstrates the 'Reliability' of detecting a difference.
# We know the groups differ by 200 units, but with n=7 and high variance (200),
# we can visualize how often the 95% CIs overlap (failing to show 'significance').

#Define Population and Simulation Parameters for 95% CIs



n <- 7          # Small sample size, common in pilot studies.
mu_veh <- 800   #Mean of Vehicle population
mu_caf <- 600   #Mean of Caffeine population
sigma <- 200    # SD of relevant samples
n_sims <- 100   # We repeat the entire experiment 100 times.
CI<-0.95
main_title<-'A) Vehicle vs. Caffeine\n95% CIs (n=7)'

plot_caterpillar_compare(CI,n,mu_veh,mu_caf,sigma,n_sims,main_title)

#--------------------Figure 7B--------------------------------------------------
# Here we use the same data and groups as Panel A, but we report them using 
# 68% CIs (representing +/- 1 SEM). 
# This illustrates why SEM can be visually deceptive.
# The logic is the same as described for figure 7A above.
#Setup Parameters for 68% CI (SEM)
n <- 7
mu_veh <- 800
mu_caf <- 600
sigma <- 200
n_sims <- 100
CI<-0.68      # 68% Confidence is the equivalent of reporting +/- 1 SEM
main_title<-'B) Vehicle vs. Caffeine\n68% CIs (n=7)'

plot_caterpillar_compare(CI,n,mu_veh,mu_caf,sigma,n_sims,main_title)


# Final step for high-quality export
# If you opened the tiff() device at the start of the script, 
# you MUST run dev.off() by removing the '#' below to finish writing the file and close it.
# dev.off()