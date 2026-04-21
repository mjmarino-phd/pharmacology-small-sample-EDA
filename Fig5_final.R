#-------------------------------------------------------------------------------
# TITLE: Fig5_final.R
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
# DESCRIPTION: This script generates figure 5 from the manuscript and is 
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

# tiff("Figure5_Final.tif", width = 8, height = 4.5, units = "in", res = 300, compression = "lzw")

#--------------------Figure 5: Boxplots of simulated caffeine study----------------

# Data Import
# We define the filename as a variable to make the script easily adaptable 
# for different experimental runs.
filename <- "CaffeineData.csv"

# Read the raw data. 
data <- read.csv(filename)

# Log-Transformation
log_data <- log10(data)

# Setup Visual Layout
# mfrow = c(2, 2): We are building a 'Quadrant View' for this figure.
# This allows us to compare raw vs. log-transformed data with and without scatter
# mar = c(4, 5, 3, 2): We tighten the bottom margin slightly (4) while 
# maintaining ample space on the left (5) for .
par(mfrow = c(2, 2), mar = c(4, 5, 3, 2))


#--------------------Figure 5A--------------------------------------------------
# Boxplot of raw data
# When 'data' is a data frame with multiple columns (Vehicle, Caffeine), 
# boxplot() automatically treats each column as a separate group on the X-axis.
boxplot(data, main = "A) Raw Data", ylab=expression("Tumor Volume (mm"^3*")"))

#--------------------Figure 5B--------------------------------------------------
# Boxplot of log-transformed data
# log10() is vectorized; when applied to a data frame, it transforms every 
# value in every column simultaneously while preserving the column names.
boxplot(log_data, main = "B) Log Transformed", 
        ylab=expression(Log[10]*"[Tumor Volume (mm"^3*")]"))

#--------------------Figure 5C--------------------------------------------------
# Boxplot raw data plus scatter

# The 'pars' List
# This argument allows us to pass specific graphical parameters to the 
# underlying 'bxp' function. 
# outcol: sets the color of the formal Tukey outliers.
# outcex: scales the size of the outlier symbol.
boxplot(data, main = "C) Raw Data + Scatter", ylab=expression("Tumor Volume (mm"^3*")"),
        pars = list(pch = 18, col = "blue", outcol = "blue",outcex=1.5))

# The Iterative Overlay Loop
# We use ncol(data) to ensure the loop runs once for every group in the data frame.
for (i in 1:ncol(data)) {
  
  # rep(i, nrow(data)): 
  # Creates a vector of 'x' coordinates matching the group index (1, 2, etc.).
  # jitter(..., amount = 0.1): 
  # Adds a small amount of random 'noise' to the x-coordinates.
  # This prevents the points from forming a single vertical line, which
  # would hide points with the same or similar Y-values.
  x_coords <- jitter(rep(i, nrow(data)), amount = 0.1)
  
  
  # data[,i]: 
  # Pulls the actual measurements for the current group being processed.
  # rgb(1, 0, 0, 0.5): 
  # Red color with 50% transparency. This 'alpha' channel is a standard 
  # coding technique to visualize 'overplotting' (high-density areas).
  points(x_coords, data[,i], pch = 16, col = rgb(1, 0, 0, 0.5))
}

#--------------------Figure 5D--------------------------------------------------
# Boxplot of log-transformed data plus scatter
#Uses same logic as Fig 5C
boxplot(log_data, main = "D) Log Data + Scatter", ylab=expression(Log[10]*"[Tumor Volume (mm"^3*")]"),
        pars = list(pch = 18, col = "blue", outcol = "blue",outcex=1.5))
for (i in 1:ncol(log_data)) {
  x_coords <- jitter(rep(i, nrow(log_data)), amount = 0.1)
  points(x_coords, log_data[,i], pch = 16, col = rgb(1, 0, 0, 0.5))
}


# Final step for high-quality export
# If you opened the tiff() device at the start of the script, 
# you MUST run dev.off() by removing the '#' below to finish writing the file and close it.
# dev.off()