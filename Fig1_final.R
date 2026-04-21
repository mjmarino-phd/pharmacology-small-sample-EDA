#-------------------------------------------------------------------------------
# TITLE: Fig1_final.R
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
# DESCRIPTION: This script generates figure 1 from the manuscript and is 
#   annotated to provide a beginner friendly introduction to the R programming language.
#-------------------------------------------------------------------------------


# Clear environment: Prevents variables from previous projects 
# from interfering with your current analysis.
rm(list = ls())

# Set a random seed for reproducibility. 
# The number (123) is arbitrary, but using a seed ensures that the 
# 'random' sampling draws the exact same values every time the script is run.
# This is essential for scientific transparency and creating consistent figures.
set.seed(123)


# Tip for Publication: Journals usually require high-resolution TIFFs.
# The code below is a "file device." When active, R redirects the plot 
# from the screen directly to a high-quality file in your working directory.
# Remove the "#" in the line below and change the Figure file name as needed
# You will also need to remove the '#' from the last line of code in this file

# tiff("Figure1_Final.tif", width = 8, height = 4.5, units = "in", res = 300, compression = "lzw")


#--------------------Figure 1: Simulation of 2 normal distributions----------------

# Set population parameters (the 'Ground Truth' we are simulating)

mu <- 800    # Control mean (e.g., untreated tumor volume in mm^3)
popSD <- 200 # Control variability (Standard Deviation)

exp_mu <- 600    # Experimental mean (e.g., 25% reduction after caffeine)
exp_popSD <- 200 # Experimental SD (assumed equal for this simulation)


# Generate the 'Bell Curve' shapes

# seq creates 1000 evenly spaced numbers (x-values) across the width of the distribution (+/- 4SD)
x_pop <- seq(mu - 4 * popSD, mu + 4 * popSD, length.out = 1000) 

# dnorm calculates the 'Probability Density' (the height of the curve or y-values) at every point in x_pop
y_pop <- dnorm(x_pop, mean = mu, sd = popSD)


# Experimental population gets a similar treatment
# seq creates 1000 evenly spaced numbers (x-values) across the width of the distribution (+/- 4SD)
exp_x_pop <- seq(exp_mu - 4 * exp_popSD, exp_mu + 4 * exp_popSD, length.out = 1000) #generates x-axis +/-4SD

# dnorm calculates the 'Probability Density' (the height of the curve or y-values) at every point in exp_x_pop
exp_y_pop <- dnorm(exp_x_pop, mean = exp_mu, sd = exp_popSD) #generates y-axis frequency values


#Set up for 3 panel figure
# par stands for 'parameters.' mfrow = c(1, 3) tells R to arrange plots 
# in 1 row and 3 columns.
par(mfrow = c(1, 3))


#--------------------Figure 1A--------------------------------------------------

#Plot theoretical populations

# The plot() function creates a NEW coordinate system and window.
# We use x_pop and y_pop to define the blue 'Vehicle' curve.
plot(x_pop, y_pop, type = "l", col = "blue", lwd = 2, #x-y blue line plot for vehicle
     main = "A) Theoretical Population Distributions", xlab = "Tumor Volume", ylab = "Frequency", #set up labels
     yaxt = "n", # yaxt = "n" hides the y-axis (frequency isn't critical here)
     xaxt = "n") # xaxt = "n" hides the default x-axis so we can customize it below

# Customizing the X-axis for better readability
# side = 1 refers to the bottom; 'at' defines where the tick marks go
axis(side = 1, at = seq(0, 1600, by = 200))  #change x-axis tics

# The lines() function is different from plot(); it layers data onto the 
# existing window without clearing it. Here we add the 'Caffeine' curve.
lines(exp_x_pop, exp_y_pop, col = "red", lwd = 2) #adds an x-y red line plot for treatment


# abline adds straight lines. 'v' stands for vertical.
# These dashed lines (lty=2) highlight the true 'center' of each population.
abline(v = mu, col = "skyblue2", lwd = 2,lty=2)
abline(v = exp_mu, col = "lightcoral", lwd = 2,lty=2)

# Legend configuration
legend("topright", #placement of legend
       legend = c("Vehicle (Control)", "Caffeine (Treated)"), #vector of text labels for legend
       col = c("blue", "red"), #legend colors
       lty = c(1, 1),  # Tells the legend to show solid lines
       bty = "n", # 'bty' stands for 'box type'; "n" means no border box
       cex = 0.8) # 'cex' scales the character size (font size)

#--------------------Figure 1B--------------------------------------------------
#Generate stacked sample plots

#Set up parameters for the simulation
small_offset <- 0.00007 # Tiny vertical step to stack each 'experiment' so they don't overlap
n_samples_to_plot <- 25 # We will simulate 25 separate 'pilot studies'
sample_size <- 7 # Each study will have n=7

# Set up the base plot (same structure as 1A above)
plot(x_pop, y_pop, type = "l", col = "skyblue2", lwd = 1,
     main = "B) Theoretical n=7 Samples", xlab = "Tumor Volume", ylab = "Frequency",
     yaxt = "n",
     xaxt = "n")
axis(side = 1, at = seq(0, 1600, by = 200))

# Shading the 'Outlier Regions' (+/-2 SD)
## adjustcolor adds transparency (alpha) so we can see the stacked samples through the shading.
col_blue_trans <- adjustcolor("skyblue2", alpha.f = 0.3)

# Drawing the left shaded tail (< 2 SD from mean)
x_left <- seq(min(x_pop), mu - 2 * popSD, length.out = 500)  #Generates a sequence of values between the minimum population value and the mean -2SD
y_left <- dnorm(x_left, mean = mu, sd = popSD) #samples the distribution 

# polygon() works by connecting coordinates in order. 
# rev(x_left) and rep(0...) tell R to trace the curve, then 'walk back' 
# along the x-axis to close the shape.
polygon(c(x_left, rev(x_left)), c(y_left, rep(0, length(y_left))), col = col_blue_trans, border = NA)

# Drawing the right shaded tail (> 2 SD from mean) same logic as teh left tail above
x_right <- seq(mu + 2 * popSD, max(x_pop), length.out = 500)
y_right <- dnorm(x_right, mean = mu, sd = popSD)
polygon(c(x_right, rev(x_right)), c(y_right, rep(0, length(y_left))), col = col_blue_trans, border = NA)


#add caffeine treated population line (as above)
lines(exp_x_pop, exp_y_pop, col = "lightcoral", lwd = 1)

# Shading the 'Outlier Regions' (+/-2 SD) (as above)
col_red_trans <- adjustcolor("lightcoral", alpha.f = 0.3)

x_left <- seq(min(exp_x_pop), exp_mu - 2 * exp_popSD, length.out = 500)
y_left <- dnorm(x_left, mean = exp_mu, sd = exp_popSD)

polygon(c(x_left, rev(x_left)), c(y_left, rep(0, length(y_left))), col = col_red_trans, border = NA)

x_right <- seq(exp_mu + 2 * exp_popSD, max(exp_x_pop), length.out = 500)
y_right <- dnorm(x_right, mean = exp_mu, sd = exp_popSD)

polygon(c(x_right, rev(x_right)), c(y_right, rep(0, length(y_left))), col = col_red_trans, border = NA)


#Simulation Loop
#This loop simulates 25 runs of the the same n=7 experiment.
for (i in 1:n_samples_to_plot) {
  
  # rnorm() generates 'observed data' by drawing random numbers from our populations
  sample_data <- rnorm(sample_size, mean = mu, sd = popSD)
  exp_sample_data <- rnorm(sample_size, mean = exp_mu, sd = popSD)
  
  # Calculate the vertical 'shelf' for this specific row
  y_height <- small_offset * i
  
  # Draw a horizontal segment representing the RANGE of the 7 samples
  # x0/y0 is the start point, x1/y1 is the end point
  # Vehicle population visualization
  segments(x0 = min(sample_data), y0 = y_height, 
           x1 = max(sample_data), y1 = y_height, 
           col = "skyblue4", lwd = 1)
  
  # Add a vertical tick (pch="|") at the sample MEAN
  points(mean(sample_data), y_height, col = "blue", pch = "|", cex = 1.2)
  
  #Repeat for Treatment  Population
  segments(x0 = min(exp_sample_data), y0 = y_height, 
           x1 = max(exp_sample_data), y1 = y_height, 
           col = "firebrick", lwd = 1)
  points(mean(exp_sample_data), y_height, col = "red", pch = "|", cex = 1.2)
}

## Legend configuration (see Fig1A above)
legend("topright", 
       legend = c("Vehicle (Control)", "Caffeine (Treated)"),
       col = c("skyblue2", "lightcoral"), 
       lty = c(1, 1, NA), 
       bty = "n", 
       cex = 0.8)




#--------------------Figure 1C--------------------------------------------------
#Fig 1C is generated by increasing sample_size to 17 below.
# Comparing Fig1B and Fig1C allows the researcher to see how 'n' influences 
# the stability of the range and the mean estimate.

small_offset <- 0.00007
n_samples_to_plot <- 25
sample_size <- 17 # Increased from n=7 to n=17

# Set up the base plot
plot(x_pop, y_pop, type = "l", col = "skyblue2", lwd = 1,
     main = "C) Theoretical n=17 Samples", xlab = "Tumor Volume", ylab = "Frequency",
     yaxt = "n",
     xaxt = "n") 
axis(side = 1, at = seq(0, 1600, by = 200)) 

# Shading the 'Outlier Regions' (+/-2 SD)
# [Shading logic remains identical to 1B to ensure visual consistency]
col_blue_trans <- adjustcolor("skyblue2", alpha.f = 0.3)

x_left <- seq(min(x_pop), mu - 2 * popSD, length.out = 500)
y_left <- dnorm(x_left, mean = mu, sd = popSD)

polygon(c(x_left, rev(x_left)), c(y_left, rep(0, length(y_left))), col = col_blue_trans, border = NA)

x_right <- seq(mu + 2 * popSD, max(x_pop), length.out = 500)
y_right <- dnorm(x_right, mean = mu, sd = popSD)

polygon(c(x_right, rev(x_right)), c(y_right, rep(0, length(y_left))), col = col_blue_trans, border = NA)


#add caffeine treated population line 
lines(exp_x_pop, exp_y_pop, col = "lightcoral", lwd = 1)

# Shading the 'Outlier Regions' (+/-2 SD)
col_red_trans <- adjustcolor("lightcoral", alpha.f = 0.3) #define transparent color

x_left <- seq(min(exp_x_pop), exp_mu - 2 * exp_popSD, length.out = 500)
y_left <- dnorm(x_left, mean = exp_mu, sd = exp_popSD)

polygon(c(x_left, rev(x_left)), c(y_left, rep(0, length(y_left))), col = col_red_trans, border = NA)

x_right <- seq(exp_mu + 2 * exp_popSD, max(exp_x_pop), length.out = 500)
y_right <- dnorm(x_right, mean = exp_mu, sd = exp_popSD)

polygon(c(x_right, rev(x_right)), c(y_right, rep(0, length(y_left))), col = col_red_trans, border = NA)

#Simulation Loop
#This loop simulates 25 runs of the the same n=17 experiment.
for (i in 1:n_samples_to_plot) {
  
  sample_data <- rnorm(sample_size, mean = mu, sd = popSD)
  exp_sample_data <- rnorm(sample_size, mean = exp_mu, sd = popSD)
  
  y_height <- small_offset * i
  
  # Vehicle population visualization
  segments(x0 = min(sample_data), y0 = y_height, 
           x1 = max(sample_data), y1 = y_height, 
           col = "skyblue4", lwd = 1)
  
  points(mean(sample_data), y_height, col = "blue", pch = "|", cex = 1.2)
  
  #Repeat for Treatment  Population
  segments(x0 = min(exp_sample_data), y0 = y_height, 
           x1 = max(exp_sample_data), y1 = y_height, 
           col = "firebrick", lwd = 1)
  
  points(mean(exp_sample_data), y_height, col = "red", pch = "|", cex = 1.2)
}

## Legend configuration
legend("topright", 
       legend = c("Vehicle (Control)", "Caffeine (Treated)"),
       col = c("skyblue2", "lightcoral"), 
       lty = c(1, 1, NA), 
       #pch = c(NA, NA, 15), # Solid square for the shaded area
       bty = "n", 
       cex = 0.8)


# Final step for high-quality export
# If you opened the tiff() device at the start of the script, 
# you MUST run dev.off() by removing the '#' below to finish writing the file and close it.
# dev.off()