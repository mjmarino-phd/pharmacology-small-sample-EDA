# pharmacology-small-sample-EDA
R scripts and workflows for the Biochemical Pharmacology manuscript:  "A Robust Workflow for Exploratory Data Analysis and Outlier Management in Small-Sample Pharmacological Datasets" (Submitted to Biochemical Pharmacology)

This repository contains the complete suite of R-scripts associated with the manuscript:  
**"A Robust Workflow for Exploratory Data Analysis and Outlier Management in Small-Sample Pharmacological Datasets"** (Marino, 2026).

##  Overview
In pharmacological research, experimental constraints often limit sample sizes to $n=7$ to $30$. These scripts provide a practical, beginner-friendly toolkit for:

* **Outlier Detection:** Applying Dixon’s Q-test and Grubbs' test with modern precision.
* **Effect Size Estimation:** Calculating Cohen’s $d$, Hedges’ $g$, and Cliff’s Delta.
* **Data Visualization:** Generating high-quality, reproducible figures for publication.
* In addition, the R-scripts employed to generate figures 1 - 9 of the manuscript are included
with substantial annotation sufficient for the novice to use as an entry to R programming

##  Repository Contents
The repository includes 16 annotated R-scripts. If you are new to R, we recommend starting with:

1.  `Fig1_final.R`: A simulation showing how small-sample stochasticity affects mean estimates
2.  `Fig2_final.R`: A demonstration of the Central Limit Theorom in action and the basis for the n=30 'parametric pivot'
3.  `Fig3_final.R`: A simulation of a log-normal population distribution and log transformation
4.  `Fig4_final.R`: An overview of the anatomy of a Tukey Boxplot
5.  `Fig5_final.R`: Application of Tukey boxplots and overlaid scatter plots to simulated experimental data
6.  `Fig6_final.R`: Use of 'caterpillar' plots to demonstrate the concept of confidence intervals
7.  `Fig7_final.R`: Application of caterpillar plots to simulated experimental populations
8.  `Fig8_final.R`: Introduction of 'Confidence Funnel' plots as a novel data visualization tool
9.  `Fig9_final.R`: A final summary figure suggested as an ideal method to report small data set results
10.  `GenerateSampleData.R`: Script used to search for the type of data set used as a case study in the manuscript.
11.  `SimpleDixonQ.R`: A highly annoted example of how the `dixon.test` function works in R
12.  `AutomatedDixonQ.R`: A fully automated implementation of the outlier screen logic presented in the manuscript provided as a tool for analyzing data in CSV files formatted like `CaffeineData.csv` containing any number of groups
13.  `Cum_ProbabilityObservingOutlier.R`: Application of the binomial theorem to estimate the probability of observing 1, 2, or 3 outliers based on sample size
14.  `EffectSIzeCalculator.R`: Calculates  Cohen's d, Hedges' g, and Cliff's Delta all with 95% CIs for data in CSV files formatted like `CaffeineData.csv` containing two groups
15.  `ConfFunnel_with_Points.R`: Generates Confidence Funnel plots for data in CSV files formatted like `CaffeineData.csv` containing any number of groups
16.  `CaffeineData.csv`: Data for case study in the manuscript.  File is required for generating Fig5,8-9 and serves as an example   template analysis scripts


##  Requirements
These scripts were developed using **R (version 4.3.3 (2024-02-29))** and utilize the following libraries:
* `effsize` (for Hedges' $g$ and Cliff's Delta)
* `outliers` (for Dixon and Grubbs tests)

##  Citation & License
**License:** MIT License.  
This code is open-source. You are free to use, modify, and distribute these scripts. 

**Citation:** If you use this code in your research, please cite the original manuscript:
> Marino, M. J. (2026). A Robust Workflow for Exploratory Data Analysis and Outlier Management in Small-Sample Pharmacological Datasets. *Biochemical Pharmacology*. [Submitted DOI TBD]
