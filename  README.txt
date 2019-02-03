This directory has the R-code used for the BACI and trend power analysis for Alberta Environment.

Key code and directories are


Functions to actually compute power for various designs

baci.power.r                <- power analysis function for BACI designs.
before-after-power-stroup.r <- power analysis function for before/after designs
slr-power-stroup.r          <- power analysis function for trend analyses
baci-power-pool-noise.r     <- power analysis function for BACI design pooling watersheds


Code to analyze the current data and create figures/graphs for report

RawData directory - contains the raw data
read.data.R  <- read in the data from the RawData directory

baci.anal.r      <- reads data, calls baci-power function, and creates plots and figures
                    for separate BACI design for each watershed
baci.anal-pool.r <- reads data, calls baci-power-pool-noise.r function and creates
                    plots and figures for combined BACI analysis over ALL watersheds
trend.anal.r     <- reads data, calls the slr-power-stroup and creates
                    plots and figures for trend analysis


Figures/Charts in report
Report  directory contains the figures/charts and various versions of the report


AnalysisReport
- A Rmarkdown document that shows how to analyze the data once it is collected
  based on simulated data.


