#-------------------------------------------------------------------------------
# TITLE: Fig3_final.R
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
# DESCRIPTION: This script generates figure 3 from the manuscript and is 
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

# tiff("Figure3_Final.tif", width = 8, height = 4.5, units = "in", res = 300, compression = "lzw")


#--------------------Figure 3: Simulation log-normal population 
                          #and effect of log transformation----------------


# Modeling Log-Normal Tumor Data

# Defining our Population Statistics
# We defined the  mean and SD in 
# arithmetic units (e.g., 800 mm3).
m <- 800 # Target Arithmetic Mean
s <- 200 # Target Arithmetic Standard Deviation

# Converting to Log-Scale Parameters
# The rlnorm() function does not take '800' and '200' directly.
# It requires the mean (location) and SD (shape) of the data 
# AFTER it has been log-transformed. 
location <- log(m^2 / sqrt(s^2 + m^2))
shape <- sqrt(log(1 + (s^2 / m^2)))

# Generating the Population
# We create 10,000 virtual tumors. Note that the log-normal 
# distribution is naturally bounded at zero,consistent with biological 
# volumes that cannot be negative.
pop_tumors <- rlnorm(n = 10000, meanlog = location, sdlog = shape)


# Setup Visual Layout
# mfrow = c(1, 2) creates a side-by-side comparison (1 row, 2 columns).
# We use generous margins (5, 5...) to ensure the superscript 
# units (mm3) aren't cut off by the edge of the window.
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))




#--------------------Figure 3A--------------------------------------------------


# The Linear Scale (Testing for Normality)

# Draw the Histogram
# We use 'prob = TRUE' so the Y-axis represents Density rather than raw counts.
# This is required to overlay the mathematical Normal curve on top.
hist(pop_tumors, breaks = 50, prob = TRUE, col = "grey80", border = "white",
     main = "A) Population of Tumor\nDensity Linear Scale",
     xlab = expression("Tumor Volume (mm"^3*")"), ylab = "Density", 
     xlim = c(0, 2000), ylim = c(0, 0.0025))

# Fit a Theoretical Normal Curve
# We generate a sequence of 500 points along the X-axis to draw a smooth line.
x_seq <- seq(0, 2000, length = 500)

# We calculate the 'dnorm' (Density of Normal) using the data's own mean and SD.
# This asks: "If this data were Normal, what would the curve look like?"
y_normal <- dnorm(x_seq, mean = mean(pop_tumors), sd = sd(pop_tumors))

# Draw the Fit Line
lines(x_seq, y_normal, col = "darkred", lwd = 3)

#  Add the 'Rug' Plot
# A rug plot places a tiny vertical line for every single data point at the bottom.
# This helps the viewer see exactly where the individual tumors are clustering,
# highlighting the density in the 'tail' that the red curve ignores.
rug(pop_tumors, lwd = 0.2, col = "red")



#--------------------Figure 3B--------------------------------------------------
# Same log-normal distribution based on tumor volume example plotted on a log scale

#Transform the Data
# We apply a base-10 logarithm to our population.
# This 'pulls in' the long right tail and spreads out the crowded 
# points at the low end, making the distribution symmetrical.
log_pop <- log10(pop_tumors)

# Draw the Histogram on the New Scale
# Notice the X-axis now represents powers of 10 (e.g., 2 = 100, 3 = 1000)
hist(log_pop, breaks = 50, prob = TRUE, col = "skyblue", border = "white",
     main = "B) Population of Tumor
     Density Log10 Scale",
     xlab = expression(Log[10]*"[Tumor Volume (mm"^3*")]"), ylab = "Density")

# Fit a New Normal Curve
# We generate a sequence based on our log-transformed values.
x_seq_log <- seq(min(log_pop), max(log_pop), length = 500)

# We calculate the Normal density based on the mean and SD of the LOGS.
y_log_normal <- dnorm(x_seq_log, mean = mean(log_pop), sd = sd(log_pop))
# Draw the Fit Line
# Observe how the dark blue line now fits the data perfectly. 
# This confirms that our tumor data is 'Log-Normal.'
lines(x_seq_log, y_log_normal, col = "darkblue", lwd = 3)

# Final step for high-quality export
# If you opened the tiff() device at the start of the script, 
# you MUST run dev.off() by removing the '#' below to finish writing the file and close it.
# dev.off()