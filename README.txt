# Discrete Trait Prediction

This repository contains the scripts for "Predicting Discrete Traits in Evolving Systems".

## Contents

- Discrete_Simulation.V#.sh
    - This shell script calls the others to perform the whole simulation.
    - You may run this to recreate our study as we performed it.
    - This holds the primary settings for this simulation.
- DiscreteFunctions.V#.R
    - This R script details all sub-functions used in the following scripts.
- SetupDirectories.V#.R
    - This R script creates the infrastructure for the study
- TreeGeneration.V#.R
    - This R script generates the phylogenetic trees for the study
- FullDataGeneration.V#.R
    - This R script simulates the data used for this study
- SampleSingleTaxa.V#.R
    - This R script samples a taxon for single tip prediction
- SampleMultipleTaxa.V#.R
    - This R script samples taxa for multiple prediction
- ResultsMatrixGeneration.V#.R
    - This R script generates empty results tables
- GatherSisterInfo.V#.R
    - This R script gathers the trait information for the sister taxa of the taxon being predicted
- FilesForBayesTraits.V#.R
    - This R script writes extra files necessary for BayesTraits
- BBandNBPrediction.V#.R
    - This file makes predicitons using the non-phylogenetic Beta Binomial method and the Naive Bayes classifier
- RunBayesTraits.V#.R
    - This R script performs all BayesTraits runs
- CompileBayesTraits.V#.R
    - This R script compiles the predictions made by BayesTraits from its output files
- SummarizeResults.V#.R
    - This R script summarizes the results tables for easier reading
- SummarizeA1Results.V#.R
    - This R script summarizes the results for trials where the unknown taxon has a 1 for Trait A
- RemoveSuperfluousFiles.V#.R
    - This R script deletes all extra files made by/for BayesTraits

## System and Program Requirements

- R 4.4.3
    - This can be found at the following URL: https://www.r-project.org/
- BayesTraits 5.0.2 (or newer)
    - This can be found at the following URL: https://www.evolution.reading.ac.uk/BayesTraitsV5.0.2/BayesTraitsV5.0.2.html
- 5 Gigabytes of extra storage space on the device you run the simulation on.

## Instructions

1. To repeat this study, download alll scripts into a folder on the computer which is to run the simulation.
2. Create a folder, name it 'Scripts', and move all R scripts in to this folder, leaving the shell script in the main directory.
3. Move the BayesTraits executable into the main folder with the shell script and the 'Scripts' folder.
4. Open the shell script, and ensure the settings are your preferred settings.
5a. From here you have two options, the first is to navigate to your powershell, make the folder with the shell script the working directory, and run the shell script.
5b. The second option is to run the R scripts individually. Open RStudio, your shell script, or your preferred method of running R code, and set the workding directory to the folder with the shell script. Then, you can run each R script individually from this main folder to see what each script is doing or to rerun specific parts of the study.
Ensure that however you run these scripts, they are written to be run from the main directory, not from the 'Scripts' directory where they should be stored. This allows additional file tidiness and streamlines the simulation study.

## Expected Outputs

When you return to the main directory, you will find several new folders. Those folders containing the simulation trees, data, and individual results which will be named based on the matrix with which they were simulated. 
Many files will be named with this pattern: "__.__", where the first letters represent the model of evolution and the second represents the rates.

- Models
    - "ER"      - This denotes an equal rates model
    - "DR"      - This denotes a direction rates model
    - "DEP1"     - This denotes a dependent rates model with lower dependency
    - "DEP2"     - This denotes a dependent rates model with higher dependency
    - "Random"  - This denotes a model where data is simulated randomly
- Rates
    - "L"   - This denotes that the base rate is LOW for an equal or dependent model.
    - "M"   - This denotes that the base rate is MEDIUM for an equal or dependent model.
    - "H"   - This denotes that the base rate is HIGH for an equal or dependent model.
    - "LM"  - This denotes that the beta rate is LOW and the alpha rate is MEDIUM for a directional model.
    - "LH"  - This denotes that the beta rate is LOW and the alpha rate is HIGH for a directional model.
    - "MH"  - This denotes that the beta rate is MEDIUM and the alpha rate is HIGH for a directional model.

There will also be a Results folder where you can find the full results table for each test.
The four columns are split into four categories: Trial information, predictive probabilities for each method, 
    accuracy rates for each method, the Log-Loss scores for each method, and information about the other taxa in the tree.

- "Trial_#" - This is the trial number for that Matrix
- "Taxon_#" - This is the tip sampled from that trial's tree
- "Trait_A" - This is the unknown taxon's character state for the FIRST trait
- "Trait_B" - This is the unknown taxon's character state for the SECOND trait
- "Terminal_Branch_Length" - The length of the unknown taxon's terminal branch

- "Beta_Bin_Prob" - This is the predictive probability that the state of Trait B is a 1, using the Beta Binomial method.
- "Naive_Prob" -    This is the predictive probability that the state of Trait B is a 1, using the Naive Bayes method.
- "MS_MCMC_Prob" -  This is the predictive probability that the state of Trait B is a 1, using the Multistate model.
- "Ind_MCMC_Prob" - This is the predictive probability that the state of Trait B is a 1, using Pagel's Independent model.
- "Dep_MCMC_Prob" - This is the predictive probability that the state of Trait B is a 1, using Pagel's Dependent model.
- "MS_RJ_Prob" -    This is the predictive probability that the state of Trait B is a 1, using the Multistate model and a Reversible-Jump algorithm.
- "Ind_RJ_Prob" -   This is the predictive probability that the state of Trait B is a 1, using Pagel's Independent model and a Reversible-Jump algorithm.
- "Dep_RJ_Prob" -   This is the predictive probability that the state of Trait B is a 1, using Pagel's Dependent model and a Reversible-Jump algorithm.

- "Beta_Bin_Acc" - This reports whether or not a researcher would make the correct prediction, using the Beta Binomial method.
- "Naive_Acc" -    This reports whether or not a researcher would make the correct prediction, using the Naive Bayes method.
- "MS_MCMC_Acc" -  This reports whether or not a researcher would make the correct prediction, using the Multistate model.
- "Ind_MCMC_Acc" - This reports whether or not a researcher would make the correct prediction, using Pagel's Independent model.
- "Dep_MCMC_Acc" - This reports whether or not a researcher would make the correct prediction, using Pagel's Dependent model.
- "MS_RJ_Acc" -    This reports whether or not a researcher would make the correct prediction, using the Multistate model and a Reversible-Jump algorithm.
- "Ind_RJ_Acc" -   This reports whether or not a researcher would make the correct prediction, using Pagel's Independent model and a Reversible-Jump algorithm.
- "Dep_RJ_Acc" -   This reports whether or not a researcher would make the correct prediction, using Pagel's Dependent model and a Reversible-Jump algorithm.

- "Beta_Bin_LL" - This reports the Log-Loss score for the Beta Binomial method.
- "Naive_LL" -    This reports the Log-Loss score for the Naive Bayes method.
- "MS_MCMC_LL" -  This reports the Log-Loss score for the Multistate model.
- "Ind_MCMC_LL" - This reports the Log-Loss score for Pagel's Independent model.
- "Dep_MCMC_LL" - This reports the Log-Loss score for the Pagel's Dependent model.
- "MS_RJ_LL" -    This reports the Log-Loss score for the Multistate model and a Reversible-Jump algorithm.
- "Ind_RJ_LL" -   This reports the Log-Loss score for Pagel's Independent model and a Reversible-Jump algorithm.
- "Dep_RJ_LL" -   This reports the Log-Loss score for Pagel's Dependent model and a Reversible-Jump algorithm.

- "00" - The number of "known" taxa that have a character state of 0,0
- "01" - The number of "known" taxa that have a character state of 0,1
- "10" - The number of "known" taxa that have a character state of 1,0
- "11" - The number of "known" taxa that have a character state of 1,1
- "Avg_Sister_A" - The average value of Trait A for the sister taxon/taxa of the unknown taxon
- "Avg_Sister_B" - The average value of Trait B for the sister taxon/taxa of the unknown taxon
