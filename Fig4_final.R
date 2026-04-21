#-------------------------------------------------------------------------------
# TITLE: Fig4_final.R
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
# DESCRIPTION: This script generates figure 4 from the manuscript and is 
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

# tiff("Figure4_Final.tif", width = 8, height = 4.5, units = "in", res = 300, compression = "lzw")


# Resetting the Canvas:
# R 'remembers' the layout from the previous figure (e.g., the 4-panel Fig 3).
# par(mfrow = c(1, 1)) forces the plotting window back to a single, 
# full-page layout. This ensures Fig 4 doesn't get 'squashed' into a corner.
par(mfrow = c(1, 1), mar = c(3, 1, 3, 1)) # Tighten margins for text


#--------------------Figure 4:  Anatomy of a boxplot--------------------------


#--------------------FUNCTIONS--------------------

#-------------------------------------------------------------------------------
# FUNCTION: label_feature
# Purpose: Draws an arrow and a text label pointing to specific boxplot features.
# This ensures that our 'Anatomy' labels don't overlap the plot itself.
#-------------------------------------------------------------------------------

label_feature <- function(y_val, label, side = 1) {
  # The 'side' argument determines if the label sits on the right or left.
  # side = 1: Right-hand side (defaults)
  # side = -1: Left-hand side
  
  if (side == 1) {
    # RIGHT SIDE: Arrow points from 1.4 back toward the box (1.15)
    x_start <- 1.4
    x_end <- 1.15 
    position <- 4   # pos=4: Places text to the RIGHT of the x_start coordinate
  
  } else {
    # LEFT SIDE: Arrow points from 0.675 forward toward the box (0.9)
    # These coordinates are shifted to provide a clean, symmetrical look.
    x_start <- 0.675 
    x_end <- 0.9
    position <- 2   # pos=2: Places text to the LEFT of the x_start coordinate
  }
  
  # Draw the pointer arrow
  # length = 0.08 defines the size of the arrowhead.
  arrows(x_start, y_val, x_end, y_val, length = 0.08, col = "darkblue")
  
  # 2. Add the text label
  # cex = 0.8 slightly shrinks the text
  text(x_start, y_val, label, pos = position, cex = 0.8)
}

#--------------------End FUNCTIONS--------------------




# Generate a generic Tukey boxplot
# Defining a Typical 'Small-n' Pharmacological Dataset
# This vector represents a standard experimental group (n=11). 
# Note the final value (95.0), which appears significantly higher than the rest.
data_skewed <- c(63.7, 34.4, 53.6, 56.3, 54.0, 48.9, 65.1, 49.1, 70.2, 49.4, 95.0)

# boxplot.stats() returns 5 values: 
# [Lower Whisker, 1st Quartile, Median, 3rd Quartile, Upper Whisker]
# Note: Whiskers are NOT the min/max if outliers are present!
stats <- boxplot.stats(data_skewed)$stats

# Defining the Interquartile Range (IQR)
# The IQR (the height of the box) represents the middle 50% of your data.
# stats[4] is the 3rd Quartile (75th percentile)
# stats[2] is the 1st Quartile (25th percentile)
# Using R's 1-indexing: Index 4 is the top of the box, Index 2 is the bottom.
iqr <- stats[4] - stats[2]

# Calculating the Tukey 'Fences'
# John Tukey defined outliers as points falling more than 1.5x the IQR 
# away from the edge of the box. 
# These fences determine where the 'whiskers' end and 'outliers' begin.

upper_fence <- stats[4] + 1.5 * iqr
lower_fence <- stats[2] - 1.5 * iqr

# Setup the Blank Canvas
# We use plot(NULL) to define our axis limits without drawing any data yet.
# xlim = c(0, 2): Since the boxplot sits at x=1, this gives us equal space 
# on the left and right for our descriptive labels.
# ylim = c(20, 100): We set the scale to comfortably fit our outlier (95.0) 
# and the lower fence.
plot(NULL, xlim = c(0, 2), 
     ylim = c(20, 100), 
     axes = FALSE, 
     xlab = "", 
     ylab = "", 
     main = "Anatomy of a Tukey Boxplot")


# Draw the actual Boxplot
# at = 1: Places the plot in the center of our 0-2 canvas.
# add = TRUE: This is critical! It tells R to draw ON TOP of our blank plot.
# col = "#FFD70044": An 'alpha-blended' gold color. The '44' at the end 
# makes it transparent, a professional look that prevents the box from 
# obscuring any background gridlines or labels.
# boxwex = 0.3: We slim the box down to 30% width to maximize the space 
# for our anatomical pointers.
boxplot(data_skewed, 
        at = 1, 
        add = TRUE, 
        col = "#FFD70044", 
        boxwex = 0.3, 
        axes = FALSE, 
        frame.plot = FALSE)


# Adding the Annotations

# Left Side Labels (side = -1)
# These identify the "Box" components. Note that we use the 
# stats[2] and stats[4] values we extracted earlier.
label_feature(stats[2], "1st Quartile (Q1)", side = -1)
label_feature(stats[4], "3rd Quartile (Q3)", side = -1)
label_feature(stats[3], "Median", side = -1)

# Right Side Labels (side = 1)
# These identify the "Whiskers" and extreme points.
# stats[1] and stats[5] are the boundaries of 'typical' data.
label_feature(stats[1], "Lower Whisker (Min)", side = 1)
label_feature(stats[5], "Upper Whisker", side = 1)

# Here we manually point to our known high value (95) 
# to label it as a 'Potential Outlier.'
label_feature(95, "Potential Outlier", side = 1)

# Visualizing the "Tukey Fences"
# These red dotted lines (lty = 3) represent the invisible 'thresholds' 
# that the boxplot() function uses to decide what gets a whisker 
# and what gets a dot.
abline(h = upper_fence, col = "red", lty = 3)
text(0.6, upper_fence + 5, "Tukey's Upper Fence (Q3 + 1.5*IQR)", 
     col = "red", cex = 0.7)

abline(h = lower_fence, col = "red", lty = 3)
text(0.6, lower_fence - 5, "Tukey's Lower Fence (Q1 - 1.5*IQR)", 
     col = "red", cex = 0.7)


# IQR Bracket (Manual drawing for clarity)
# This 'manual' bracket emphasizes that the box height is a 
# statistical unit (the Interquartile Range).
brackets_x <- 1.2 # Position the bracket just to the right of the box

# Draw the main vertical line of the bracket
segments(brackets_x, stats[2], brackets_x, stats[4], lwd = 1.5)

# Draw the tiny horizontal 'ticks' to cap the bracket at Q1 and Q3
segments(brackets_x - 0.02, stats[2], brackets_x, stats[2])
segments(brackets_x - 0.02, stats[4], brackets_x, stats[4])

# Label the bracket
# We use the average of Q1 and Q3 to center the text 'IQR' vertically.
text(brackets_x + 0.05, (stats[2] + stats[4]) / 2, "IQR", pos = 4, font = 2)

# 2. Restoring the Y-Axis
# Since we suppressed the axes earlier for a clean look, we add 
# a custom axis now. 
# pos = 0.3: This shifts the axis slightly inward from the edge of 
# the plot for a modern, 'floating' aesthetic
axis(2, at = seq(0, 140, 20), pos = 0.3)

# Final step for high-quality export
# If you opened the tiff() device at the start of the script, 
# you MUST run dev.off() by removing the '#' below to finish writing the file and close it.
# dev.off()