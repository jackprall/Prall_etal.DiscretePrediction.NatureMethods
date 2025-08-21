#!/bin/bash
################################################################################
## Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
## June 2025

################################################################################
## Running the script below will test multiple phylogenetic and non-phylogenetic prediction models. 
## This should be placed in a directory with the BayesTraits executable and a
## a directory 'Scripts' which houses the scripts below.

## Running these scripts will simulate trees, data, and directories to store them.
## It will test multiple phylogenetic and non-phylogenetic prediction models.
## It will also create other files that are necessary to run the tests.

## Outputs will include any files made by BayesTraits and results tables.
## Results are stored within a 'Results' folder found within the folder for each matrix.

################################################################################

## Simulation settings
## Simulation version
sim_version="V15"

## Set the number of iterations
num_iterations=1000

## Set the number of taxa that are used to build the tree and in the prediction trials
##    "Pop size" is the number of taxa in the dataset
pop_size=500

## Set the MCMC algorithms being used.
##    Use 'MCMC', 'RJMCMC', or 'BOTH'
RJmodel=BOTH

## Define the types and trials arrays (it may be helpful to split these up for longer simulations)
types=("ER.L" "ER.M" "ER.H" "DR.LM" "DR.LH" "DR.MH" "DEP1.L" "DEP1.M" "DEP1.H" "DEP2.L" "DEP2.M" "DEP2.H")

## Set the rates of simulated evolution
low_rate=0.1
medium_rate=0.25
high_rate=0.5

## Determine whether we are testing variable rates models of evolution.
## Set a vector of scale factors which will be used to create rate-heterogeneity when simulating data
variable_rates=true        # Must be "true" or "false"
tree_scales=(0.5 1 1 1 1 2 3)

## Set the degree of dependency (multiplier on rates) and adjustment so that the sum of the rates remains constant
depend_scale1=0.5
depend_adj1=1.16667

depend_scale2=0.1
depend_adj2=1.3

## Set the settings for the prediction target
unknown_size=10            # Only used when either of the next two are "true"
multiple_prediction=true   # Must be "true" or "false"
clade_prediction=true      # Must be "true" or "false"

## Set whether or not to add an additional BayesTraits test that only uses one trait
multistate_prediction=true # Must be "true" or "false"



## Call the scripts that will perform the simulation
## Start the simulation by setting up the directories
Rscript Scripts/SetupDirectories.${sim_version}.R

## Second, we need to generate the phylogenetic trees
Rscript Scripts/TreeGeneration.${sim_version}.R

## Next, generate the full datasets that we will be using for this simulation
Rscript Scripts/FullDataGeneration.${sim_version}.R

## Run the scripts responsible for selecting the unknown taxa
Rscript Scripts/SampleSingleTaxa.${sim_version}.R

if [[ "$multiple_prediction" == "true" || "$clade_prediction" == "true" ]]; then 
  Rscript Scripts/SampleMultipleTaxa.${sim_version}.R
fi

## Run the next script, which will build the final (empty) results tables
Rscript Scripts/ResultsMatrixGeneration.${sim_version}.R

## Gather the information from the nearest sister taxon
## This will get stored in the results tables
Rscript Scripts/GatherSisterInfo.${sim_version}.R

## Run this script, which will write all necessary files for BayesTraits
Rscript Scripts/FilesForBayesTraits.${sim_version}.R

## This script will make predictions using Beta Binomial and Naive Bayes
## The results will get stored in the results tables
Rscript Scripts/BBandNBPrediction.${sim_version}.R

## This script will run any BayesTraits runs that you have in the settings
Rscript Scripts/RunBayesTraits.${sim_version}.R

## This script will read the output files from BayesTraits and compile them into the results tables
Rscript Scripts/CompileBayesTraits.${sim_version}.R

## This script will summarize all of the results files into more readable summaries
## These will be found in the "Results" folder in the main directory
Rscript Scripts/SummarizeResults.${sim_version}.R

## This script will summarize ONLY the results of trials where Trait_A is 1
## These will be found in the "Results" folder in the main directory
Rscript Scripts/SummarizeA1Results.${sim_version}.R

## This final script removes all of the files created at intermediate steps
## These files are either superfluous or unreadable to humans, thus deemed unnecessary
## You will be left will all trees, full data, prediction data, and BayesTraits prediction log files.
Rscript Scripts/RemoveSuperfluousFiles.${sim_version}.R
