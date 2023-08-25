Species Ditribution Models for Alps Data <a name="TOP"></a>
===================
In this repository you will find some code revolving around producing statistical models aimed at understanding what are the most important factors to predict animal distribution; the code is free to download but beware that, due to the large amount of data used, it is not advised to reproduce without proper instrumentation.

# Contents:  #

## code       ##

the *code* folder contains the scripts, written in R and Python languages at version 4.2 and 3.10.
Here, the *species_distribution_model_Script_V1.Rmd* file represents a first run of SDM applied to Capra Ibex; 
the script is divided in sections, with the following steps being mandatory:
  - Loading the __environmental variables__, joined in a **rasterStack** item and the __species presence records__, in **SpatialPointsDataFrame format**. In this case, we used two DataFrames of presence-only records of Capra Ibex and several species of Bats within the Alpine territory;
  - __BIOMOD_FormatingData__ --> used to format the loaded data to make it compatible with **biomod2** and, if necessary, generate pseudo-absence data to be used as training set. Four different selection strategies: disk, sre, random, user-defined; random was used for this script;
  - **BIOMOD_Modeling** --> define all the algorithms to implement and their associated parameters; here, RF, GLM, GBM.
  - **get_evaluations** --> used to perform the statistical analysis, applying the specified models to our compatible data.

The script sdm_b_c.Rmd applies RF, GBM, GLM and MaxEnt to Capra Ibex and Bats data.
The future conditions were taken from WorldClim, choosing scenario's SSP126 and SSP370 of the algorithm EC-Earth3-Veg, the best performing for simulation on European territory.

From here, many information can be obtained to perform statistical analysis and distribution assessments: 
- **Variable Importance** ranking with corresponding contribution for each factor.
- **Model Evaluations**, based on statistics such as TSS, ROC curves or Cohen's KAPPA.
- **Response Curves** corresponding to each of our factors.
- **Model projections**, seeing how close are the model's results with empirical distribution.
- **Ensemble Modeling**, joining all the results of the different algorithms to produce a more accurate prediction of distributions.
- **Future Estimations**, if given the adequate settings concerning how the environment is expected to change.

## pictures      ##

This folder contains images which where produced using the models (outcome, references and much more).
They were all produced by me, using R and Python's plot functions.



## resources     ##

Notes, links and anything useful for these kind of scripts.
  ### sdm references ###

Here I uploaded all graphs, tables and diagrams found while documenting myself on SDMs, always quoting the source in the name in the format of (Surname, Year).
