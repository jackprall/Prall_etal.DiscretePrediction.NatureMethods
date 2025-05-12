################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# October 2024

################################################################################
# Script uses functions defined in 'DiscreteFunctions'
# Calls BayesTraits output files from 'Discrete_Simulation' and 'TestModels'
# Gathers phylogenetic predictions and compiles all results into one of three files, listed below.

# Outputs:
    #'*PH.Results.txt' = table with the results of the phylogenetic prediction methods
    #'*Full.Results.txt = table with the complete results of this simulation for
    #       a given matrix or 'type' 
################################################################################

# This script will be used to compile the full results off all of the multiple tip trials

### Here, we will set up the parameters for the rest of the study
# First, make sure to call the ape library
library(ape)

# Get the current working directory
original_dir <- getwd()

## Get the types and trials arrays
bash_file <- "ExtendedResultsRerun.sh"
types_list <- grep("^types=", readLines(bash_file), value = TRUE)
cleaned_string <- gsub("types=\\(|\\)", "", types_list)
cleaned_string <- gsub("\"", "", cleaned_string)
types <- strsplit(cleaned_string, " ")[[1]]

# Number of iterations when running as batch job/ on HPC 
iter_line <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", iter_line))

# Get the sampling type and number
# Get the sampling type and number
unknown_n <- as.numeric(gsub(".*=(.*)", "\\1", grep("^unknown_n=", readLines(bash_file), value = TRUE)))
sample_type <- gsub(".*=(.*)", "\\1", grep("^sample_type=", readLines(bash_file), value = TRUE))
  
# Get the MCMC type
RJmodel <- grep("^RJmodel=", readLines(bash_file), value = TRUE)
model <- gsub(".*=(.*)", "\\1", RJmodel)

# Pull necessary functions for this script
source(Sys.glob("DiscreteFunctions.V*"))


### This is a series of for-loops that will pull all of the information with the settings provided
# This loop will cycle through each "type", which is determined by the q-matrix used to simulate the data
for (type in types) {

  # Set your directories, and create the data matrix
  original_dir <- getwd()
  setwd(type)

  Results.Full <- matrix(nrow = 0, ncol = 33)
  colnames(Results.Full) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", "Terminal_Branch_Length",
                              "Beta_Bin_Prob", "Naive_Prob", "MS_MCMC_Prob", "Ind_MCMC_Prob", 
                              "Dep_MCMC_Prob", "MS_RJ_Prob", "Ind_RJ_Prob", "Dep_RJ_Prob",
                              "Beta_Bin_Acc", "Naive_Acc", "MS_MCMC_Acc", "Ind_MCMC_Acc",
                              "Dep_MCMC_Acc", "MS_RJ_Acc",  "Ind_RJ_Acc",  "Dep_RJ_Acc",
                              "Beta_Bin_LL", "Naive_LL", "MCMC_MS_LL", "MCMC_Ind_LL", 
                              "MCMC_Dep_LL", "RJ_MS_LL", "RJ_Ind_LL", "RJ_Dep_LL", "#_00", "#_01", 
                              "#_10", "#_11",)

  # Next, we will loop through each trial
  for (i in 1:trial_amount) {

    # We start on this level by pulling the necessary data
    #First call the tree, unknown taxon and trait A
    dataname <- paste0(type, ".", i, ".predict_data.txt")
    data <- read.table(file = dataname)

    tree_name <- paste0(type, ".", i, ".full_tree.tre")
    full_tree <- read.nexus(file = tree_name)
      
    rTaxon <- which(data[, 3] == "?", arr.ind = TRUE)
    rTaxonValueA <- data[rTaxon, 2]
      
    # Then trait B
    if (length(rTaxon) > 1) {
      Bname <- paste0("./Results/", type, "TraitBInfo.txt")
      Bdata <- read.table(file = Bname, row.names = NULL)
      Bdata <- Bdata[,-1]
      rTaxonValueB <- as.vector(unlist(Bdata[i, ]))
    } else {
      Bname <- paste0("./Results/", type, ".NP.Results.txt")
      Bdata <- read.table(file = Bname, skip = 1)
      rTaxonValueB <- Bdata[i, 4]
    }

    # Next, we will pull the Beta Binomial prediciton, which doesn't change for each tip prediction, but instead for each new dataset
    BBname <- paste0("./Results/", type, ".NP.Results.txt")
    BBdata <- read.table(file = BBname, skip = 1)
    BB.Prob <- BBdata[i, 6]
        
    # We also need to get the Naive Bayes predictions
    countname <- paste0(type, ".", i, ".Counts.txt")
    counts <- read.table(file = countname, skip = 1)
    counts <- counts[, 4:5]
        
    # Finally, the Naive predictions
    NB.Predictions <- Naive_Bayes_predict(counts, rTaxonValueA)


    ### This loop will cover any situations where rTaxon is more than one value
    for (j in 1:length(rTaxon)) {
      # First we need to pull the probs
      NB.Prob <- NB.Predictions[j]

      skiplength <- length(rTaxon)
      Calculate_Post_Prob(i, model, skiplength, j, trial_name = type)

      # Second, we will calculate the accuracy of each of these
      BB.acc <- calculate_accuracy(rTaxonValueB[j], BB.Prob)
      NB.acc <- calculate_accuracy(rTaxonValueB[j], NB.Prob)
      MS.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], MS.MCMC.Prob)
      Ind.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Ind.MCMC.Prob)
      Dep.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Dep.MCMC.Prob)
      MS.RJ.acc <- calculate_accuracy(rTaxonValueB[j], MS.RJ.Prob)
      Ind.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Ind.RJ.Prob)
      Dep.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Dep.RJ.Prob)

      # Next, we calculate the LogLoss for each of these values
      BB.LL  <- LogLoss(rTaxonValueB[j], BB.Prob)
      NB.LL  <- LogLoss(rTaxonValueB[j], NB.Prob)
      MS.MCMC.LL  <- LogLoss(rTaxonValueB[j], MS.MCMC.Prob)
      Ind.MCMC.LL  <- LogLoss(rTaxonValueB[j], Ind.MCMC.Prob)
      Dep.MCMC.LL  <- LogLoss(rTaxonValueB[j], Dep.MCMC.Prob)
      MS.RJ.LL  <- LogLoss(rTaxonValueB[j], MS.RJ.Prob)
      Ind.RJ.LL  <- LogLoss(rTaxonValueB[j], Ind.RJ.Prob)
      Dep.RJ.LL  <- LogLoss(rTaxonValueB[j], Dep.RJ.Prob)

      # Add the Terminal Branch Length information
      if (length(rTaxon) > 1) {
        TBL <- unknown_branch_length(rTaxon[j], full_tree)
      } else {
        TBL <- unknown_branch_length(rTaxon, full_tree)
      }

      # Finally, we assign the trial and add it to the matrix
      Trial <- c(i, rTaxon[j], rTaxonValueA[j], rTaxonValueB[j], TBL,
          BB.Prob, NB.Prob, MS.MCMC.Prob, Ind.MCMC.Prob,
          Dep.MCMC.Prob, MS.RJ.Prob, Ind.RJ.Prob, Dep.RJ.Prob,
          BB.acc, NB.acc, MS.MCMC.acc, Ind.MCMC.acc,  
          Dep.MCMC.acc, MS.RJ.acc,  Ind.RJ.acc,  Dep.RJ.acc,
          BB.LL, NB.LL, MS.MCMC.LL, Ind.MCMC.LL, 
          Dep.MCMC.LL, MS.RJ.LL, Ind.RJ.LL, Dep.RJ.LL,
          counts[1,1], counts[1,2], counts[2,1], counts[2,2])

      # Finally, add it to the phylogenetic results table
      Results.Full <- rbind(Results.Full, Trial)

    }

  }

  # Write the table into your results and return to the original directory
  finalname <- paste0("./Results/", type, ".Full.Results.txt")
  write.table(Results.Full, file = finalname, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
  # Clean up some of the additional, superflous scripts
  if (length(list.files(getwd())) > 46000) {
    file.remove(list.files(pattern = "\\.Rates\\."))
    file.remove(list.files(pattern = "\\.In\\.txt"))
    file.remove(list.files(pattern = "\\.ModelFile\\.bin"))
    file.remove(list.files(pattern = "\\.Schedule\\."))
    file.remove(list.files(pattern = "edited"))
  }    

    setwd(original_dir)
}
