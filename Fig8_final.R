#-------------------------------------------------------------------------------
# TITLE: Fig8_final.R
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
# DESCRIPTION: This script generates figure 8 from the manuscript and is 
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

# tiff("Figure8_Final.tif", width = 8, height = 4.5, units = "in", res = 300, compression = "lzw")


#--------------------Figure 8:  Introducing Confidence Funnels------------------


#--------------------FUNCTIONS--------------------


#-------------------------------------------------------------------------------
# FUNCTION: plot_conf_funnel_with_points
# Purpose: Creates a "Confidence Funnel"
# It shows how the Confidence Interval expands as you demand higher levels 
# of certainty (from 50% up to 99.9%). And allows a visual assessment of effect size 
# and confidence in effect
#
# Inputs: 
#   data        : A data frame or list where each column is an experimental group.
#   main_title  : The text displayed at the top of the plot.
#   x_title     : The label for the horizontal axis (e.g., "Tumor Volume").
#   show_legend : Logical (TRUE/FALSE) to toggle the group legend.
#-------------------------------------------------------------------------------

plot_conf_funnel_with_points <- function(data,main_title,x_title,show_legend) {
  # Set up groups and manage colors
  # Use 'hcl.colors' for professional, distinct palettes.
  groups <- names(data)
  n_groups <- length(groups)
  line_colors  <- hcl.colors(n_groups, palette = "Set 2")
  
  # 'adjustcolor' is a powerful tool for novices. 
  # Adjust or modify a vector of colors by “turning knobs” on one or more 
  #coordinates in  (r,g,b,alpha) space
  # 'alpha.f' sets transparency:
  # 0.15 is very faint (good for large areas), 0.6 is bolder (good for points).
  fill_colors  <- adjustcolor(line_colors, alpha.f = 0.15)
  point_colors <- adjustcolor(line_colors, alpha.f = 0.6)
  
  # Setup Plot Window
  # We calculate the min/max of the entire dataset to ensure every group fits.
  x_min <- min(unlist(data), na.rm = TRUE)
  x_max <- max(unlist(data), na.rm = TRUE)
  pad   <- (x_max - x_min) * 0.1 # Adds 10% 'breathing room' on the sides.
  
  plot(NULL, xlim = c(x_min - pad, x_max + pad), ylim = c(50, 100),
       xlab = x_title, ylab = "Confidence Level (%)",
       main = main_title, frame.plot = FALSE)
  
  # Add a subtle grid to help the eye track horizontal percentages.
  grid(nx = NULL, ny = NULL, col = "gray90", lty = "solid")
  
  # Loop through each group
  for (i in 1:n_groups) {
    vals <- na.omit(data[[groups[i]]]) # Remove 'NA' (missing) values
    n <- length(vals)
    if(n < 3) next # Statistics require at least 3 samples to calculate variance.
    
    # Calculate key statistics
    m <- mean(vals)
    se <- sd(vals) / sqrt(n) # the standard error
    
    # Generate a sequence of confidence levels from 50% to 99.9%
    conf_levels <- seq(0.50, 0.999, length.out = 100)
    
    # Calculate the widening bounds using the Student's t-distribution.
    # As 'conf_levels' increases, the 't-critical' value gets much larger.
    lower_bound <- m - qt(conf_levels + (1 - conf_levels) / 2, df = n - 1) * se
    upper_bound <- m + qt(conf_levels + (1 - conf_levels) / 2, df = n - 1) * se
    
    # Drawing the Funnel (The Polygon)
    # A polygon needs a 'path' of points that goes out along the bottom 
    # and comes back along the top to close the loop.
    px <- c(lower_bound, rev(upper_bound))
    py <- c(conf_levels * 100, rev(conf_levels * 100))
    polygon(px, py, col = fill_colors[i], border = NA)
    
    # Overlaying Data Points
    # We place the raw data at the 50% mark (the 'floor' of our plot).
    # This shows the reader exactly what raw data generated the funnel.
    points(x = vals, y = rep(50, n),
           pch = 19, # Solid circles
           col = point_colors[i],
           cex = 1.3) # Make them slightly larger for readability
   
    # Draw a vertical dashed line at the mean for reference.
    abline(v = m, col = line_colors[i], lty = 2, lwd = 2)
    
    # Highlight the standard 95% threshold with a solid horizontal bar.
    error_95 <- qt(0.975, df = n - 1) * se
    segments(x0 = m - error_95, y0 = 95, 
             x1 = m + error_95, y1 = 95, 
             col = line_colors[i], lwd = 2)
    
  }
  
  #Show legend
  if(show_legend) {
    legend("bottomright", legend = groups,
           fill = fill_colors, border = line_colors, 
           bty = "n", cex = 0.9, y.intersp = 1.2,
           inset = c(0.05, 0.05)) # Nudges it 5% away from the right and bottom edges
  }
}

#--------------------END FUNCTIONS--------------------



# Data Import
# We define the filename as a variable to make the script easily adaptable 
# for different experimental runs.
filename <- "CaffeineData.csv"

# Read the raw data. 
data <- read.csv(filename)

# Log-Transformation
log_data <- log10(data)

#Based on the boxplots (Figure5) and Dixon's Q (see publication text), we know 
#that the highest value in the caffeine group is a statistical outlier.
#Therefore we will remove the value
cleaned_log_data<-log_data
cleaned_log_data$Caffeine[which.max(cleaned_log_data$Caffeine)] <- NA

# Setup Visual Layout
# mfrow = c(1, 3) creates a side-by-side comparison for 3 plots (1 row, 3 columns).
par(mfrow = c(1, 3))

#--------------------Figure 8A--------------------------------------------------
#plot Confidence Funnel from raw data
# We begin by plotting the un-transformed data. This provides a baseline 
# to show impact of refinement (Log-transformation, Outlier exclusion)
main_title<-'A) Raw Data'
x_title<-expression("Tumor Volume (mm"^3*")")
# Note: We set show_legend to FALSE here because we will add a 
# single, unified legend to the final panel (8C) to save space.
plot_conf_funnel_with_points (data, main_title,x_title,show_legend=FALSE)

#--------------------Figure 8B--------------------------------------------------
#plot Confidence Funnel from Log10 transformed data
# By applying a Log10 scale, we transform 'multiplicative' biological noise 
# into 'additive' normal noise. This makes the funnels more symmetrical, 
# which is a requirement for the math behind parametric statistical testing
main_title<-'B) Log Data'
x_title<-expression(Log[10]*"[Tumor Volume (mm"^3*")]")
# We plot 'log_data' here, which still includes the suspected outlier.
plot_conf_funnel_with_points (log_data,main_title,x_title,show_legend=FALSE)

#--------------------Figure 8c--------------------------------------------------
#plot Confidence Funnel from outlier removed Log10 transformed data
# This represents our most statistically sound and robust analysis.
main_title<-'C) Outlier Removed'
x_title<-expression(Log[10]*"[Tumor Volume (mm"^3*")]")
# We use show_legend = TRUE here to provide a 'key' for the entire 3-panel figure
plot_conf_funnel_with_points (cleaned_log_data,main_title,x_title,show_legend=TRUE)


# Final step for high-quality export
# If you opened the tiff() device at the start of the script, 
# you MUST run dev.off() by removing the '#' below to finish writing the file and close it.
# dev.off()