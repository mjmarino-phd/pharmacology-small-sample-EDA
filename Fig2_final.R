#-------------------------------------------------------------------------------
# TITLE: Fig2_final.R
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
# DESCRIPTION: This script generates figure 2 from the manuscript and is 
#   annotated to provide a beginner friendly introduction to the R programming language.
#-------------------------------------------------------------------------------


# Clear environment: Prevents variables from previous projects 
# from interfering with your current analysis.
rm(list = ls())

# Set a random seed for reproducibility. 
# The number (42) is arbitrary, but using a seed ensures that the 
# 'random' sampling draws the exact same values every time the script is run.
# This is essential for scientific transparency and creating consistent figures.
set.seed(42)


# Tip for Publication: Journals usually require high-resolution TIFFs.
# The code below is a "file device." When active, R redirects the plot 
# from the screen directly to a high-quality file in your working directory.
# Remove the "#" in the line below and change the Figure file name as needed
# You will also need to remove the '#' from the last line of code in this file

# tiff("Figure2_Final.tif", width = 8, height = 4.5, units = "in", res = 300, compression = "lzw")


#--------------------Figure 2: CLT and Parametric Pivot----------------


#--------------------FUNCTIONS--------------------

#-------------------------------------------------------------------------------
# FUNCTION: get_means
# Purpose: Simulates the 'Sampling Distribution of the Mean.' 
# This is the engine for demonstrating the Central Limit Theorem.
#
# Inputs: 
#   n_size     : The number of individuals in each sample.
#   iterations : How many times to repeat the experiment (Default = 1000).
#
# Output: A vector of 1000 sample means.
#-------------------------------------------------------------------------------

get_means <- function(n_size, iterations = 1000) {
  # replicate() runs the expression 1000 times and collects the results into a vector.
  replicate(iterations, mean(sample(pop, n_size, replace = TRUE)))
}


#-------------------------------------------------------------------------------
# FUNCTION: plot_clt_overlay
# Purpose: Generates a density histogram of sample means and overlays a 
# theoretical Normal distribution to visualize the Central Limit Theorem.
#
# Inputs: 
#   data    : Vector of sample means (calculated via get_means).
#   n_label : The sample size (n) to be displayed in the title.
#   color   : The fill color for the histogram bars.
#
# Dependencies: Requires global variables 'my_breaks' and 'x_range'.
#-------------------------------------------------------------------------------
plot_clt_overlay <- function(data, n_label, color) {
  
  # Calculate statistics for the overlay
  mu <- mean(data)      # The average of our sample means
  sem_val <- sd(data)   # The Standard Error of the Mean (SD of the sampling distribution)
  
  # Draw the Density Histogram
  # prob = TRUE ensures the y-axis is 'Density' so the Normal curve fits the scale.
  hist(data, prob = TRUE, breaks = my_breaks, xlim = x_range, 
       ylim = c(0, 0.004),
       main = "", # We leave this blank to use the more flexible title() command below
       xlab = "Value", 
       col = color, 
       border = "white")
  
  # Add a custom title with the sample size
  title(main = paste("Sample Means (n=", n_label, ")", sep = ""), line = 0.8)
  
  # Generate coordinates for a theoretical Normal Curve based on our sample stats
  x_seq <- seq(x_range[1], x_range[2], length = 500)
  y_seq <- dnorm(x_seq, mean = mu, sd = sem_val)
  
  # Draw the Normal Curve in dark red to contrast with the histogram bars
  lines(x_seq, y_seq, col = "darkred", lwd = 2)
  
  # Add a legend box to show the numeric results for comparison
  legend("topright", 
         legend = c(paste("Mean:", round(mu, 1)), 
                    paste("SEM:", round(sem_val, 1))),
         bty = "n", 
         cex = 0.9, 
         text.col = "darkred")
}

#--------------------END FUNCTIONS--------------------



# Creating a non-normal 'weird' population
# Real-world data is rarely perfect. Here we build a "weird" distribution 
# using samples and repetitions to create multiple peaks (multimodal).
a <- sample(0:200, 20)     # A small cluster at the low end
b <- sample(200:500, 200)  # A larger middle cluster
c <- rep(1200, 100)        # A sharp peak at 1200
d <- rep(2000, 100)        # A sharp peak at 2000
pop <- c(a, b, c, d)       # Combine them into one population
mean_pop <- mean(pop)      # The 'True' mean we want to estimate

# Set the visual stage
# To compare plots fairly, they must all share the same X-axis and 'bins'.
x_range <- range(pop)
bin_width <- 50 
my_breaks <- seq(floor(x_range[1]), ceiling(x_range[2]) + bin_width, by = bin_width)


# Setup complex layout

# The layout() function is more flexible than mfrow. 
# We define a matrix to tell R exactly which plot goes in which 'cell'.
# Here, we create a 2x2 grid:
# -------------------------------------------------------------------------
# The numbers in this matrix represent the 'Order of Appearance' in your code.
# [1, 4]  <- Row 1: The 1st plot called (2A) and the 4th plot called (2C).
# [2, 3]  <- Row 2: The 2nd plot called (2B-1) and the 3rd plot called (2B-2).
#
# This 'out-of-order' numbering is a powerful trick. It allows us to 
# write the code in a logical narrative (Population -> Samples -> Theory) 
# while placing the 'Theory' plot (2C) in the top-right corner.
# -------------------------------------------------------------------------
layout(matrix(c(1, 4, 2, 3), nrow = 2, byrow = TRUE))

# par(mar = ...) adjusts the 'Margins' around each individual plot.
# The order is: c(bottom, left, top, right).
# We give a bit more space to the bottom (4.5) and left (4.5) 
# to ensure the 'Value' and 'Density' labels aren't cut off.
par(mar = c(4.5, 4.5, 3, 2))




#--------------------Figure 2A--------------------------------------------------

# Plotting The 'Weird' Population 
# This plot shows the distribution of every individual item in our population.
# Notice how it has multiple peaks and doesn't look like a 'Bell Curve' at all.
hist(pop, breaks = my_breaks, xlim = x_range,
     main = "A) Weird Population Distribution", 
     xlab = "Measured Value", col = "orange", border = "white")

# Adding the 'True Mean' of the population for reference.
# No matter how many samples we take, this is the value we are trying to estimate.
legend("topleft", 
       legend = paste("Mean:", round(mean_pop, 1)),
       bty = "n", 
       cex = 0.9, 
       text.col = "darkorange")


#--------------------Figure 2B--------------------------------------------------

# Draw samples and calculate sample statistics
# We use our 'get_means' function to simulate 1000 samples.
# Each sample (n=5) is drawn and the average is recorded.
means_5 <- get_means(5)
means_5_mean <- mean(means_5) # The 'Mean of Means' for n=5
means_5_sd <- sd(means_5)     # This SD of means is actually the SEM 

# Now we simulate a larger sample (n=30).
means_30 <- get_means(30)
means_30_mean <- mean(means_30) # The 'Mean of Means' for n=30
means_30_sd <- sd(means_30)     # SEM fro n=30

# Plotting the n=5 Results (Fig 2B-1)

plot_clt_overlay(means_5, "5", "lightblue")

# Master Label: This adds the 'B) Sampling' header above the individual panel.
# line = 2 moves the text higher so it doesn't crowd the plot title.
title(main = "B) Sampling", adj = 0, line = 2) 

# Plotting the n=30 Results (Fig 2B-2)
# This demonstrates the Central Limit Theorem in action. The distribution 
# becomes Gaussian, centering precisely on the population mean.
plot_clt_overlay(means_30, "30", "salmon")



#--------------------Figure 2C--------------------------------------------------

# This curve demonstrates the relationship between sample size and SEM and 
# defines the 'Parametric Pivot.' at n=30

# Define the relationship
n <- 1:100            # We look at sample sizes from 1 up to 100
se_relative <- 1/sqrt(n) # The Relative Standard Error 

# Draw the Curve
plot(n, se_relative, type="l", lwd=3, col="navy",
     main="C) The Relationship: Size vs. Precision",
     xlab="Sample Size (n)", ylab="Relative Std Error (1/√n)",
     frame.plot=FALSE, # Removes the top and right box borders for a cleaner look
     xaxs="i", yaxs="i") # 'i' ensures the axes touch the (0,0) point exactly

# Add Vertical Reference Lines
# These connect this mathematical theory back to our actual simulations in Panel B.
abline(v=5, col="red", lty=2)      # The 'Unstable' zone (from Fig 2B-1)
abline(v=30, col="darkgreen", lty=2) # The 'Pivot' zone (from Fig 2B-2)

# Text Annotations
# pos=4 places the text to the right of the coordinates
text(5, 0.9, "n=5", pos=4, cex=0.8, col="red")
text(30, 0.4, "n=30: Pivot", pos=4, cex=0.8, col="darkgreen")


# Final step for high-quality export
# If you opened the tiff() device at the start of the script, 
# you MUST run dev.off() by removing the '#' below to finish writing the file and close it.
# dev.off()