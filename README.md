# Discrete Trait Prediction

This repository contains the scripts for "Predicting Discrete Traits in Evolving Systems".

## Contents

- Discrete_Simulation.sh
    - This is a shell script that calls the others to perform the whole simluation.
    - This script also runs the phylogenetic analysis through the software BayesTraits.
    - You may run this to recreate our study as we performed it.
    - This holds the primary settings for this simulation.
- DiscreteFunctions.R
    - This is an R script that details all sub-functions used in the next two scripts.
- TestModels.R
    - This is an R script that defines a function that will test two non-phylogenetic prediction methods as well produce the input files necessary to run BayesTraits
    - The function defined in this script can be used to test other evolutionary models or phylogenetic tree types.
- SimInstructions.R
    - This is an R script that calls 'DiscreteFunctions' and 'TestModels' to test different prediction models over thirteen different evolutionary models.
    - The matrices being tested are all found in this script.
- CompileResults.R
    - This is an R script that calls 'DiscreteFunctions' to compile the outputs of 'Discrete_Simulation' and 'TestModels' into Results files
