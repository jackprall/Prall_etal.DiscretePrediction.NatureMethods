################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender, and Chris Organ
# October 2024

################################################################################
# Below defines the function called in 'SimInstructions'
# Script uses functions defined in 'DiscreteFunctions'
# Running 'Test_NP_Models' generates tree and trait data files for ONE q matrix across many trees 
#   and runs & writes out non-phylogenetic predictions 
# BayesTraits phylogenetic analyses called in external Bash job script 'Discrete_Simulation'
################################################################################

Test_NP_Models <- function(tree_size, pop_size, tree_scales, trial_amount, qmat, ultrametric = TRUE,
                      symmetric = TRUE, model, sample_type, unknown_n, trial_name = NULL, Random = FALSE) {

  # Get the current working directory
  original <- getwd()
  
  #Start by defining the results table and its column names
  Results <- matrix(nrow = 0, ncol = 9)

  colnames(Results) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                         "Branch_length", "Beta_Bin_Prob", 
                         "Beta_Bin_Accuracy", "Naive_Prob", "Naive_Accuracy")

  # If this run samples more than one tip, we need to save that trait data for later
  if (unknown_n > 1) {
    Trait_B_Data <- matrix(nrow = 0, ncol = as.numeric(unknown_n))
    NB_Predictions <- matrix(nrow = 0, ncol = as.numeric(unknown_n))
  }

  #Create the loop to that will run each trial
  for (i in 1:trial_amount) {

####Generate trees and data
    # Generate a random tree with pop_size tips from the tree function
    full_tree <- generate_tree(tree_size, pop_size, ultrametric, symmetric)
    
    # Next, standardize the tip count and tree length
    full_tree <- standardize_tree(full_tree, i, pop_size, trial_name)
    
    #This generates the data along the tree and stores it in the table "TempMat"
    TempMat <- generate_data(full_tree, tree_scales, pop_size, qmat, Random)
    
    #This function samples and saves a random taxon from the list
    rTaxon <- sample_taxon(full_tree, TempMat, sample_type, unknown_n)
    
    # Save the trait data
    rTaxonValueA <- as.numeric(TempMat[rTaxon, 2])
    rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
    
    # This takes the sampled taxon
    # and returns a tree and data table without that taxon.
    edited_TempMat <- edit_tree(i, rTaxon, full_tree, TempMat, trial_name)
    
    # This function records the number of each taxon for each character state.
    Counts <- countup(i, edited_TempMat, trial_name)
    
    # This function will record the branch length leading to the unknown taxon.
    BL <- unknown_branch_length(rTaxon, full_tree)
    
    #### Non-phylogenetic statistical methods
    # This is where we can use the simulated data
    # and non-phylogenetic methods to predict discrete traits.
    
    # Next, this function will make a Bayesian
    # prediction using a Beta Binomial distribution
    BB.Prob <- Beta_Bin_predict(edited_TempMat, unknown_n)
    
    # Then, this function will use Bayes Theorem to predict the value.
    Naive.Prob <- Naive_Bayes_predict(Counts, rTaxonValueA)

#### These lines will create the files for the phylogenetic runs
    # First determine which model types will be used for BayesTraits
    # Write the MCMC instructions that you need
    if (model %in% c("MCMC", "BOTH")) {
      # Write the settings for an MCMC run of BayesTraits
      #Second, the rate calculation for the MCMC run.
      Multistate_Settings(i, FirstRun = TRUE, RJ = FALSE, trial_name)
      #Last, the prediction settings for the MCMC run.
      Multistate_Settings(i, FirstRun = FALSE, RJ = FALSE, trial_name)
      #First, the rate calculation for the independent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = FALSE, trial_name)
      #Second, the rate calculation for the dependent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = FALSE, trial_name)
      # Next, the prediction settings for the independent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = FALSE, trial_name)
      #Last, the prediction settings for the dependent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = FALSE, trial_name)
    }

    # Write the RJMCMC instructions that you need
    if (model %in% c("RJMCMC", "BOTH")) {
      # Write the settings for an MCMC run of BayesTraits
      #First, the rate calculation for the RJ run.
      Multistate_Settings(i, FirstRun = TRUE, RJ = TRUE, trial_name)
      # Next, the prediction settings for the RJ run.
      Multistate_Settings(i, FirstRun = FALSE, RJ = TRUE, trial_name)      
      #First, the rate calculation for the independent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, trial_name)
      #Second, the rate calculation for the dependent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = TRUE, trial_name)
      # Next, the prediction settings for the independent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = TRUE, trial_name)
      #Last, the prediction settings for the dependent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = TRUE, trial_name)
    }

    #Also, we have to add back in the dataset for predicting
    Replace_taxon_data(i, TempMat, rTaxon, trial_name)

#### Add to results and summary steps

    #This labels each results for each method was correct or not
    BB.acc <- calculate_accuracy(rTaxonValueB, BB.Prob)
    Naive.acc <- calculate_accuracy(rTaxonValueB, Naive.Prob)

    # If this run samples more than one tip, we need to save that trait data for later
    if (unknown_n > 1) {
      # Pad rTaxonValueB with NULL values if it's shorter than unknown_n
      if (length(rTaxonValueB) < unknown_n) {
        rTaxonValueB <- as.numeric(c(rTaxonValueB, rep(NA, unknown_n - length(rTaxonValueB))))
        Naive.Prob <- as.numeric(c(Naive.Prob, rep(NA, unknown_n - length(Naive.Prob))))
      }
  
      # Now bind the padded or correctly sized vector to the dataframe
      Trait_B_Data <- rbind(Trait_B_Data, rTaxonValueB)
      NB_Predictions <- rbind(NB_Predictions, Naive.Prob)
    
      #Store the results in the Results Table
      Trial <- c(i, "Various", "Various", "Various", BL, 
                 BB.Prob, BB.acc, "Various", Naive.acc)
      Results <- rbind(Results, Trial)
    } else {

      #Store the results in the Results Table
      Trial <- c(i, rTaxon, rTaxonValueA, rTaxonValueB, BL, 
                 BB.Prob, BB.acc, Naive.Prob, Naive.acc)
      Results <- rbind(Results, Trial)
    }
  }

  #This will print the full Results Table
  print(Results)

  # Write out results table
  ResultFileName <- paste0("./Results/", trial_name, ".NP.Results.txt")
  write.table(Results, ResultFileName, quote = F, sep = "\t", row.names = T)

  if (unknown_n > 1) {
    Trait_B_dataname <- paste0("./Results/", trial_name, "TraitBInfo.txt")
    write.table(Trait_B_Data, Trait_B_dataname, quote = F, sep = "\t")

    NB_filename <- paste0("./Results/", trial_name, ".NB.Predictions.txt")
    write.table(NB_Predictions, NB_filename, quote = F, sep = "\t")
  }
}