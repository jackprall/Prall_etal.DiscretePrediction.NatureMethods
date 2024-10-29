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

Test_NP_Models <- function(pop_size, trial_amount, qmat, ultrametric = TRUE,
                      symmetric = TRUE, trial_name = NULL, Random = FALSE) {

    # Get the current working directory
    original <- getwd()
  
    #Start by defining the results table and its column names
    Results <- matrix(nrow = 0, ncol = 9)

  colnames(Results) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                         "Branch_length", "Beta_Bin_Prob", 
                         "Beta_Bin_Accuracy", "Naive_Prob", "Naive_Accuracy")

    #Create the loop to that will run each trial
    for (i in 1:trial_amount) {

####Generate trees and data
      # Generate a random tree with pop_size tips from the tree function
      generate_tree(pop_size, ultrametric, symmetric, trial_name)

      # Next, standardize the tip count and tree length
      standardize_tree(full_tree, i, pop_size, trial_name)

      #This generates the data along the tree and stores it in the table "TempMat"
      generate_data(full_tree, pop_size, qmat, Random)

      #This function samples and saves a random taxon from the list
      sample_taxon(TempMat, trial_name)

      # This takes the sampled taxon
      # and returns a tree and data table without that taxon.
      edit_tree(i, rTaxon, full_tree, TempMat, trial_name)

      # This function records the number of each taxon for each character state.
      countup(i, edited_TempMat, trial_name)

      # This function will record the branch length leading to the unknown taxon.
      unknown_branch_length(rTaxon, full_tree)

#### Non-phylogenetic statistical methods
      # This is where we can use the simulated data
      # and non-phylogenetic methods to predict discrete traits.

    # Next, this function will make a Bayesian
    # prediction using a Beta Binomial distribution
    Beta_Bin_predict(edited_TempMat)

    #Finally, this function will use Bayes Theorem to predict the value.
    Naive_Bayes_predict(Counts, rTaxonValueA)

### These lines will create the files for the phylogenetic runs

    #First, the rate calculation for the independent run.
    Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, trial_name)
    #Second, the rate calculation for the dependent run.
    Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, trial_name)
    # Next, the prediction settings for the independent run.
    Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, trial_name)
    #Last, the prediciton settings for the dependent run.
    Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, trial_name)

    #Also, we have to add back in the dataset for predicting
    Replace_taxon_data(i, TempMat, rTaxon, trial_name)

#### Add to results and summaryization steps

    #This labels each results for each method was correct or not
    calculate_NP_accuracy(rTaxonValueB, BB.Prob, Naive.Prob)

    #Store the results in the Results Table
    Trial <- c(i, rTaxon, rTaxonValueA, rTaxonValueB, BL, 
     BB.Prob, BB.acc, Naive.Prob, Naive.acc)
    Results <- rbind(Results, Trial)
  }

  #This will print the full Results Table
  print(Results)

  # Write out results table
  ResultFileName <- paste0("./Results/", trial_name, ".NP.Results.txt")
  write.table(Results, ResultFileName, quote = F, sep = "\t", row.names = T)
}