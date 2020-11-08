# PM2.5-Nonattainment
This repository contains files for constructing and analyzing data sets used in the paper "An empirical evaluation of the causal impact of NAAQS nonattainment designations on particulate pollution and health" by Zigler, Choirat, and Dominici

Supporting data files can be obtained at: https://dataverse.harvard.edu/dataverse/pm97naaqs
The scripts contained herein also require the following R packages:
1) AREPA, https://github.com/czigler/arepa
2) HEIfunctions, https://github.com/czigler/HEIfunctions

# -- IMPORTANT NOTE -- #
The actual analysis in the paper relies on Medicare data that cannot be openly shared.  The data available on the Dataverse contain fake, simulated health data that mimics the format of the Medicare data used for the analysis.  The R scripts in this repository are the same as those used for the analysis, but analysis results cannot be exactly reproduced due to reliance on Medicare data. This means that effect estimates will differ, but it also means that other data summaries (even sample sizes) will differ because some data restrictions were employed with Medicare data.


# -- CONSTRUCTING THE DATA -- #
Raw data files used to construct the data are available at https://dataverse.harvard.edu/dataverse/pm97naaqs and should be obtained prior to running the scripts in this repository.

The folder 'Data Construction' contains a series of R scripts to construct the a data set of the form used in the analysis (using fake Medicare data).  These scripts should be run in numerical order according to the file name, and will produce a series of intermediate data sets and finally a data set called 'pmanalysis.RData',  which can read into the R script to implement the analysis (described below).

Note that the missing baseline pollution imputations used in the analysis are provided here in the form of `20151210pm2002_2004preds.Rda'.  These can be reconstructed with the file '03 Impute Missing Baseline Pollution.R'


# -- ANALYSIS OF THE DATA -- #
The 'Analysis Files' portion of this repository contains:
- 'Paper Analysis.Rmd' (and .pdf).  This file will start with the 'pmanalysis.RData' file produced with the Data Construction files.  This file generates basic data summaries, estimates the propensity score, and outputs a data file called 'pm_withps.Rda'.  'Paper Analysis.Rmd' also summarizes the pollution- and health-outcomes analyses, implementation and output of which are described below.  

- 'fake_med'.  This folder contains the R script to run the principal stratification analysis of health outcomes (PS_med.R), the batch script to submit the job to the computing cluster (PSbatch.txt), the R output files for the analysis of each health outcome (PS_medX.Rout), and the statistical output used to produce analysis results for each health outcome (pstratamod_XX.Rda).  The latter files are read into the 'Paper Analysis.Rmd' script to summarize results.  Note that for reproducibility (and because this is an analysis using fake Medicare data anyway), the outputs provided here contain only 51 MCMC iterations.

- 'fake_poll_med'.  This folder contains the R script to run the analysis of the pollution outcome (Poll_med.R), the batch script to submit the job to the computing cluster (Pollbatch.txt), the R output files for the analysis (Poll_med.Rout), and the statistical output used to produce analysis results for the pollution outcome (pollutionmodel.RData).  The latter file is read into the 'Paper Analysis.Rmd' script to summarize results.  Note that for reproducibility (and because this is an analysis using fake Medicare data anyway), the outputs provided here contain only 51 MCMC iterations.

- 'fake_med_noprune'.  Analogous to 'fake_med', but without propensity score pruning.

-'fake_poll_med_noprune',  Analogous to 'fake_poll_med', but without propensity score pruning. 

NOTE: a data file called 'pm_withps_nomed.Rda' is available on the Dataverse.  This file can also be read into the 'Paper Analysis.Rmd' script.  It contains the actual data set used for the analysis in the paper, but with each Medicare variable replaced with '0'.  It can be used to reproduce some basic data summaries, but it cannot be used to re-implement the entire analysis.



