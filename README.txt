# Discrete Trait Prediction

This repository contains the scripts for "Predicting Discrete Traits in Evolving Systems".

## Contents

- Discrete_Simulation.V#.sh
    - This shell script calls the others to perform the whole simulation.
    - This script also runs the phylogenetic analysis through the software BayesTraits.
    - You may run this to recreate our study as we performed it.
    - This holds the primary settings for this simulation.
- DiscreteFunctions.V#.R
    - This R script details all sub-functions used in the next two scripts.
- TestModels.V#.R
    - This is an R script that defines a function that will test two non-phylogenetic prediction methods as well as produce the input files necessary to run BayesTraits
    - The function defined in this script can be used to test other evolutionary models or phylogenetic tree types.
- SimInstructions.V#.R
    - This R script calls 'DiscreteFunctions' and 'TestModels' to test different prediction models over thirteen different evolutionary models.
    - The matrices being tested are all found in this script.
- CompileResults.V#.R
    - This R script calls 'DiscreteFunctions' to compile the outputs of 'Discrete_Simulation' and 'TestModels' into Results files
 
## System and Program Requirements

- R 4.4.3
    - This can be found at the following URL: https://www.r-project.org/
- BayesTraits 4.1.3
    - This can be found at the following URL: https://www.evolution.reading.ac.uk/BayesTraitsV4.1.3/BayesTraitsV4.1.3.html
- 5 Gigabytes of extra storage space on the device you run the simulation on.

These scripts were run using the HPC research cluster at Montana State University known as Tempest. Using the settings found in 'Discrete_Simulation.sh', the code will download any required R packages, run the full simulation, and summarize the results in a specified file. 
More information on the Tempest research cluster can be found at https://www.montana.edu/uit/rci/tempest/

## Instructions

To repeat the study, simply download all five scripts into a folder on the computer which is to run the simulation. Next, open the command line of your computer and navigate to the folder containing these scripts. Then run 'Discrete_Simulation.sh' using a line that should look like "bash Discrete_Simulation.sh" 
This took 5-7 days to run using the research cluster at Montana State University using the settings described but may take longer on a less powerful machine.

## Expected Outputs

When you return to the directory, you will find several new folders. Those folders containing the simulation trees, data, and individual results which will be named based on the matrix with which they were simulated. It will read "__.__", where the first letters represent the model of evolution and the second represents the rates.

- Models
    - "ER"      - This denotes an equal rates model
    - "DR"      - This denotes a direction rates model
    - "DEP"     - This denotes a dependent rates model
    - "Random"  - This denotes a model where data is simulated randomly
- Rates
    - "L"   - This denotes that the base rate is LOW for an equal or dependent model.
    - "M"   - This denotes that the base rate is MEDIUM for an equal or dependent model.
    - "H"   - This denotes that the base rate is HIGH for an equal or dependent model.
    - "LM"  - This denotes that the beta rate is LOW and the alpha rate is MEDIUM for a directional model.
    - "LH"  - This denotes that the beta rate is LOW and the alpha rate is HIGH for a directional model.
    - "MH"  - This denotes that the beta rate is MEDIUM and the alpha rate is HIGH for a directional model.

There will also be a Results folder where you can find the full results table for each test. It will have 33 columns of data for each trial.
They are split into four categories: General information, predictive probabilities for each method, 
    accuracy rates for each method, the Log-Loss scores for each method, and finally the state frequencies for each trial.

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

- "#_00" - The is number of "known" taxa that have a character state of 0,0
- "#_01" - The is number of "known" taxa that have a character state of 0,1
- "#_10" - The is number of "known" taxa that have a character state of 1,0
- "#_11" - The is number of "known" taxa that have a character state of 1,1

## Most common error

If you notice that you are missing the full results, this was most likely due to a lack of variation in the simulated data of a run. BayesTraits needs to see a variation in the trait for which it is predicting (Trait B); otherwise, it will fail to run. Here are steps to fix that issue.

First, find where the simulation stopped. (This is most likely to happen in the trials with low rates)
- Check each of the Results files within each matrix's folder. Do this in the same order as the "stypes" setting in the 'Discrete_Simulation.sh' script.
- If you find only a 'NP.Results' file, at least one BayesTraits run wasn't completed for that matrix's trials.
- Within that folder, locate the error trial. There are several ways to do this, but some tips are that it will lack both of BayesTraits' outputted schedule files and its '.log.txt' files will be around an order of magnitude smaller than the other trials.

Once the problematic file(s) has been found...
- Edit the 'predict_data' and 'edited_data' files to add the necessary variation to allow BayesTraits can run. 
- Be sure you edit the same traits for the same taxa for both files and adjust the non-phylogenetic predictive probabilities as necessary.
- Rerun BayesTraits for that trial, making sure all output files are added to the directory storing all information for that trial (e.g. "ER.L").

From here, you can retry running 'CompileResults' to summarize all results.
