#-------------------------------------------------------------------------------
# TITLE: ConfFunnel_with_Points.R
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
# DESCRIPTION: This script provides a tool that will plot confidence funnels 
# from any data set set in a CSv file With the the individual group data 
# arranged in columns with column titles are placed in The top row.  
# See sample data set in Caffeine.csv
# The script is annotated to provide a beginner friendly introduction to the R 
# programming language, however this is a more advanced script
#-----------------------------------------------------------------------------



# Clear environment
rm(list = ls())


# ----------------change file name and data source---------------------#
filename<-"CaffeineData.csv"
x_title<- "Log10 [Tumor Volume (mm3)]"
y_title<- "Confidence Level (%)"
main_title<-"Log Data Confidence Funnel"

#--------now ctrl-A to select all and ctrl-enter to run--------#

#Read data
data <- read.csv(filename)
log_data<-log10(data)

# Setup Layout
par(mfrow=c(1,1))

#Function to plot confidence funnel with sample points
plot_conf_funnel_with_points <- function(data,main_title,x_title,show_legend) {
  # 1. Identify groups and setup dynamic colors
  groups <- names(data)
  n_groups <- length(groups)
  
  # Set distinct line colors, but keep them bold
  line_colors <- hcl.colors(n_groups, palette = "Set 2")
  # Fill colors get high transparency (alpha=0.15) for the shaded areas
  fill_colors <- adjustcolor(line_colors, alpha.f = 0.15)
  # NEW: Point colors with moderate transparency to handle overlap
  point_colors <- adjustcolor(line_colors, alpha.f = 0.6)
  
  # 2. Setup the plot area with global limits
  x_min <- min(unlist(data), na.rm = TRUE)
  x_max <- max(unlist(data), na.rm = TRUE)
  
  # A bit of horizontal padding (10%) makes the plot look better
  pad <- (x_max - x_min) * 0.1
  
  plot(NULL, xlim = c(x_min - pad, x_max + pad), ylim = c(50, 100),
       xlab = x_title, ylab = y_title,
       main = main_title, frame.plot = FALSE)
  
  # 3. Add background grid for easier reading
  grid(nx = NULL, ny = NULL, col = "gray90", lty = "solid")
  
  # 4. Loop through the variable number of groups
  for (i in 1:n_groups) {
    # Isolate and clean the data for this group
    vals <- na.omit(data[[groups[i]]])
    n <- length(vals)
    if(n < 3){
      warning(paste("Group '", groups[i], "' has only", n, 
                            "points. Skipping plot (minimum n=3 required for t-distribution)."))
      next} # Dixon and t-dist need at least n=3
    
    # Calculate key statistics (using t-distribution for small n)
    m <- mean(vals)
    se <- sd(vals) / sqrt(n)
    
    # Generate the smooth probability gradient
    conf_levels <- seq(0.50, 0.999, length.out = 100)
    # qt gives the t-critical value, using n-1 degrees of freedom
    lower_bound <- m - qt(conf_levels + (1 - conf_levels) / 2, df = n - 1) * se
    upper_bound <- m + qt(conf_levels + (1 - conf_levels) / 2, df = n - 1) * se
    
    # Draw the shaded funnel (polygon)
    # The high transparency lets groups overlap without hiding information
    px <- c(lower_bound, rev(upper_bound))
    py <- c(conf_levels * 100, rev(conf_levels * 100))
    polygon(px, py, col = fill_colors[i], border = NA)
    
    # ---  Overlay Individual Data Points ---
    # We plot the 'vals' on the X-axis. For the Y-axis, we place all
    # points on the 50% line. This visually positions them where
    # the probability is highest (near the median).
    points(x = vals, y = rep(50, n),
           pch = 19, # Solid circles
           col = point_colors[i],
           cex = 1.3) # Make them slightly larger for readability
    # --------------------------------------------
    
    # Draw the vertical mean line (dashed)
    abline(v = m, col = line_colors[i], lty = 2, lwd = 2)
    
    # Draw the 95% reference segment (bold horizontal line)
    # A t-statistic of qt(0.975) corresponds to a two-tailed 95% CI
    error_95 <- qt(0.975, df = n - 1) * se
    segments(x0 = m - error_95, y0 = 95, 
             x1 = m + error_95, y1 = 95, 
             col = line_colors[i], lwd = 2)
    
    # Add a small label "95%" near the segment
    text(x = m + error_95, y = 96, labels = "95%", col = line_colors[i], cex=0.7, adj=0)
  }
  
  # Plot Legend
  # Uses point and shaded color markers
  if(show_legend) {
    legend("bottomright", legend = groups,
           fill = fill_colors, border = line_colors, 
           bty = "n", cex = 0.9, y.intersp = 1.2,
           inset = c(0.05, 0.05)) # Nudges it 5% away from the right and bottom edges
  }
}


# Execute

plot_conf_funnel_with_points (log_data, main_title,x_title,show_legend=TRUE)
