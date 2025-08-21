################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender, and Chris Organ
# June 2025

################################################################################
# Below defines the functions called in each other R script in this study
# Script uses functions defined in 'DiscreteFunctions'
# Each function has a unique purpose in the simulation, which is further explained below

################################################################################



# This function is to help the R scripts read T/F in the bash file
read_logical <- function(key, file) {
  line <- grep(paste0("^", key, "="), readLines(file), value = TRUE)
  # Remove any comment and trailing/leading whitespace
  value <- trimws(gsub("#.*", "", sub(".*=", "", line)))
  tolower(value) == "true"
}



# This function will generate the q-matrix used to simulate the data
generate_qmat <- function(gain, loss, scaler, adjuster) {
  qmat <- matrix(data = c(0, gain*adjuster, gain*scaler, 0,
                          loss*adjuster, 0, 0, gain*adjuster,
                          loss*adjuster, 0, 0, gain*adjuster,
                          0, loss*adjuster, loss*scaler, 0), 
                          nrow = 4, ncol = 4, byrow = TRUE)
  return(qmat)
}



# This function simulates trait data and puts it in a temporary matrix.
generate_data <- function(full_tree, tree_scales, pop_size, qmat, Random = FALSE) {
  # First, using the vector in the settings, scale the branch lengths to model transition-rate heterogeneity
  branch_count <- Nedge(full_tree)
  branch_scales <- sample(tree_scales, branch_count, replace = TRUE)
  sim_tree <- full_tree
  sim_tree$edge.length <- full_tree$edge.length * branch_scales
  
  # Run rTraitDisc with the random root state to create a dataset
  TempData <- rTraitDisc(sim_tree, model = qmat, 
                         states = c("00", "01", "10", "11"), 
                         root.value = 1)
  TempMat <- matrix(TempData, nrow = pop_size, ncol = 4)
  colnames(TempMat) <- c("Taxon #", "Trait A", "Trait B", "Class")
  
  # Rework the matrix to make it a little easier to analyze
  for (j in 1:pop_size) {
    if (TempMat[j, 1] == "00") {
      TempMat[j,2] <- 0
      TempMat[j,3] <- 0
    }
    if (TempMat[j, 1] == "01") {
      TempMat[j,2] <- 0
      TempMat[j,3] <- 1
    }
    if (TempMat[j, 1] == "10") {
      TempMat[j,2] <- 1
      TempMat[j,3] <- 0
    }
    if (TempMat[j, 1] == "11") {
      TempMat[j,2] <- 1
      TempMat[j,3] <- 1
    }
    
    # This is where we can ignore the phylogeny, if we should choose
    if (Random == TRUE) {
      TempMat[, 2] <- rbinom(pop_size, size = 1, prob = 0.5)
      TempMat[, 3] <- rbinom(pop_size, size = 1, prob = 0.5)
    }
  }
  
  #Now re-add the taxon numbers to this matrix
  TempMat[, 4] <- TempMat[, 1]
  TempMat[, 1] <- 1:pop_size
  
  #Assign the Temporary Matrix to the global environment
  return(TempMat)
}



# This function is used to randomly sample a clade within the tree
find_random_clade <- function(full_tree, unknown_n) {
  # Get all internal nodes
  internal_nodes <- (length(full_tree$tip.label) + 1):max(full_tree$edge)
  
  # Shuffle nodes to randomize search order
  internal_nodes <- sample(internal_nodes)
  
  # A sub-function to find a clade with a given number of tips
  find_clade_with_n_tips <- function(n) {
    for (node in internal_nodes) {
      clade_tips <- extract.clade(full_tree, node)$tip.label
      if (length(clade_tips) == n) {
        return(clade_tips)
      }
    }
    return(NULL)
  }
  
  # Try finding a clade with exactly unknown_n tips
  result <- find_clade_with_n_tips(unknown_n)
  if (!is.null(result)) return(result)
  
  # Expand search range: try -1, then -2, then -3
  for (offset in 1:50) {
    if (unknown_n - offset > 0) {  # Ensure we don’t look for a negative tip count
      result <- find_clade_with_n_tips(unknown_n - offset)
      if (!is.null(result)) return(result)
    }
  }
  
  # If no clade is found, return a randomly chosen tip
  return(sample(full_tree$tip.label, 1))
}



# This function removes the sampled taxon from the tree and data set.
edit_data<-function(i, rTaxon, TempMat, trial_name = NULL) {   
  #Create a new "edited" matrix, dropping the taxon number row.
  edited_TempMat <- TempMat[, c(2:3)]
  
  # Assign the first row as row names
  rownames(edited_TempMat) <- TempMat[, 1]
  
  # Set the trait value of the sampled taxon to missing for prediction.
  edited_TempMat[rTaxon, 2] <- "-"
  
  # Name and save this new data set
  rates_name <- paste0(trial_name, ".", i, ".rates_data.txt")
  write.table(edited_TempMat, file = rates_name, sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
  
  # Add another matrix that will be saved to predict
  # Save the unknown trait as ? instead of -
  edited_TempMat[rTaxon, 2] <- "?"

  # Name and save the prediction data set
  predict_name <- paste0(trial_name, ".", i, ".predict_data.txt")
  write.table(edited_TempMat, file = predict_name, sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
}



# This function saves a table of the number of taxa for each character state.
countup <- function(i, TempMat, path_name) {
  # Count how many of each class there are
  s00 <- sum(TempMat[, 4] == "00")
  s01 <- sum(TempMat[, 4] == "01")
  s10 <- sum(TempMat[, 4] == "10")
  s11 <- sum(TempMat[, 4] == "11")
  
  # Create a matrix for counts
  counts <- matrix(0, nrow = 2, ncol = 2)
  colnames(counts) <- c("A=0", "A=1")
  rownames(counts) <- c("B=0", "B=1")
  
  # Fill the matrix with counts based on class values
  counts[1, 1] <- s00  # A=0, B=0
  counts[2, 1] <- s01  # A=0, B=1
  counts[1, 2] <- s10  # A=1, B=0
  counts[2, 2] <- s11  # A=1, B=1
  
  # Print to confirm mid-simulation
  print(counts)
  
  # Save the count table
  filename <- paste0(path_name, ".", i, ".Counts.txt")
  write.table(counts, file = filename, sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)

  return(counts)
}



# This function will record the terminal branch length of the unknown taxon.
unknown_branch_length <- function(rTaxon, full_tree) {
  # Identify the correct edges
  leaf_indices <- which(full_tree$tip.label %in% as.character(rTaxon))
  
  # Pull out the branch lengths
  branch_lengths <- full_tree$edge.length[which(full_tree$edge[,2] %in% leaf_indices)]
  
  # If rTaxon is a vector, compute the average branch length
  if (length(branch_lengths) > 1) {
    branch_length <- mean(branch_lengths, na.rm = TRUE)
  } else {
    branch_length <- branch_lengths
  }
  
  # Assign the branch length to be recorded with the rest of the trial
  return(branch_length)
}



# This function will return all descendant tip labels from any node(s)
getAllDescendantTips <- function(tree, nodes) {
  tip_labels <- c()
  
  for (node in nodes) {
    if (node %in% tree$tip.label) {
      # Already a tip
      tip_labels <- c(tip_labels, node)
    } else {
      # Not a tip — get descendants
      desc_nodes <- getDescendants(tree, node)
      
      # Separate tips and internal nodes
      tips <- desc_nodes[desc_nodes <= Ntip(tree)]
      internals <- desc_nodes[desc_nodes > Ntip(tree)]
      
      # Add tip labels
      tip_labels <- c(tip_labels, tree$tip.label[tips])
      
      # Recurse into internal nodes (if any)
      if (length(internals) > 0) {
        tip_labels <- c(tip_labels, getAllDescendantTips(tree, internals))
      }
    }
  }
  return(unique(tip_labels))
}



# This function will record the average value of the unknown taxon's sisters
sisterData <- function(tree, taxon, data, trait = "A") {
  # Get sister nodes of your unknown taxon (use numeric nodes or labels as needed)
  sisters <- getSisters(tree, taxon, mode = "number")  # use mode="number" for numeric nodes
  
  # Get all descendant tip labels of the sister nodes recursively
  descendant_tips_all <- getAllDescendantTips(tree, sisters)
  
  # Extract trait values for those tips from data
  if (!all(descendant_tips_all %in% rownames(data))) {
    stop("Some descendant tips not found in the data.")
  }
  
  # Return the correct trait data
  if (trait == "A") {
    sisters_A <- mean(data[descendant_tips_all, 2], na.rm = TRUE)
    return(sisters_A)
  } else if (trait == "B") {
    sisters_B <- mean(data[descendant_tips_all, 3], na.rm = TRUE)
    return(sisters_B)    
  } else {print("Specified trait not found.")}
}



# This function predicts the unknown trait using the Beta Binomial method.
Beta_Bin_predict <- function(edited_TempMat, unknown_n) {
  
  #Start by calculating y, which is the number of successes (or 1s)
  x <- nrow(edited_TempMat)
  y <- sum(edited_TempMat[, 2] == 1)
  
  #We can now calculate the alpha and beta posterior distribution parameters
  #This is because we assume a prior distribution ~Beta(1,1)
  Shape1 <- y + 1
  Shape2 <- 1 + x - y 
  
  #Then we sample the distribution
  pi.values <- rbeta(1000, Shape1, Shape2)
  
  # Next, we will use rbinom with pi as the probability to create
  # a Bayesian posterior predictive distribution.
  PiDist <- rbinom(1000,1, prob = pi.values)
  
  #Get a frequency of 1s from the rbinom, and that's the posterior probability
  frequency <- (sum(PiDist)/1000)
  
  # Repeat the frequency value to match the length of unknown_n and assign it to the environment
  BB.Prob <- rep(frequency, unknown_n)
  return(frequency)
}



# This function will use Bayes Theorem to predict with the Naive Bayes classifier.
Naive_Bayes_predict <- function(Counts, rTaxonValueA) {
  # This is Bayes Theorum both in general and for our case
  # posterior = prior * likelihood / marginal
  # p(B=1|A) = p(B=1) * p(A|B=1) / P(A)
  # However, because we have the counts for each tree, this can be simplified to P(B=1|A)1
  posterior <- ifelse(
    rTaxonValueA == 0,
    Counts[2,1] / (Counts[1,1] + Counts[2,1]),
    Counts[2,2] / (Counts[1,2] + Counts[2,2]))
  
  #Finally, we assign this to the global environment
  return(posterior)
}



# This function sets the Multistate analysis settings
Multistate_Settings <- function(i, FirstRun = TRUE, RJ = TRUE, output_path = NULL) {
  # First, define model saving/loading behavior
  if (FirstRun == TRUE) {
    a <- "SaveModels "
    c <- ".Rates"
    if (RJ == TRUE) {
      e <- "
rjhp exp 0 5"
    } else {
      e <- "
PriorAll exp 2.5"
    }
  } else {
    a <- "LoadModels "
    c <- ".Predict"
    e <- ""
  }
  
  if (RJ == TRUE) {
    f <- "RJMCMC."
  } else {
    f <- "MCMC."
  }
  
  # Set the name of the run
  runname <- paste0(output_path, ".Multistate.", f, i, c)
  
  # Compile the settings for Multistate analysis
  settings <- paste0("1
2
  
it 1100000
BurnIn 100000
Sample 1000
logFile ", runname, e, "
", a, output_path, ".Multistate.", f, i, ".ModelFile.bin

Run")
  
  # Save the settings file
  writeLines(text = settings, con = paste0(runname, ".In.txt"))
}



# This function defines the settings for BayesTraits runs.
# This is for the Independent and Dependent models.
Bayes_Settings <- function(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, output_path = NULL) {
  #First, the differentiate between the different BayesTraits settings
  if (FirstRun == TRUE) {
    a <- "SaveModels "
    c <- ".Rates"
    if (RJ == TRUE) {
      e <- "
rjhp exp 0 5"
    } else {
      e <- "
PriorAll exp 2.5"
    }
  } else {
    a <- "LoadModels "
    c <- ".Predict"
    e <- ""
  }

  if (RJ == TRUE) {
    f <- "RJMCMC."
  } else {
    f <- "MCMC."
  }

  #Second, set the run to either independent or dependent
  #b is for the BT setting, d is for naming
  if (IndependentCharacters == TRUE) {
    b <- "2" 
    d <- ".Ind."
  } else {
    b <- "3"
    d <- ".Dep."
  }

  #Set the name of the run
  runname <- paste0(output_path, d, f, i, c)

  #Now, we compile the settings
  settings <- paste0(b, "
2

it 1100000
BurnIn 100000
Sample 1000
logFile ", runname, e, "
", a, output_path, d, f, i, ".ModelFile.bin

Run")

  #Finally, we save the file to be used later
  writeLines(text = settings, con = paste0(runname, ".In.txt"))
}



# This function calculates the posterior probability from BayesTraits' output files
Calculate_Post_Prob <- function(i, model, multistate_prediction, unknown_n, j, trial_name = NULL) {
  # A sub-function to calculate posterior probability for a specific column
  compute_posterior <- function(logname, skip_lines, col_offset, j) {
    output <- read.table(logname, header = FALSE, skip = skip_lines, nrows = 1000)
    target_col <- col_offset + j - 1  # Adjust column based on j
    return(mean(output[, target_col]))
  }
  
  MS.skip <- 43 + unknown_n
  Ind.skip <- 47 + unknown_n
  Dep.skip <- 55 + unknown_n
  
  if (model %in% c("RJMCMC", "BOTH")) {
    assign("Ind.RJ.Prob", compute_posterior(paste0("./", trial_name, ".Ind.RJMCMC.", i, ".Predict.Log.txt"), Ind.skip, 9, j), envir = .GlobalEnv)
    assign("Dep.RJ.Prob", compute_posterior(paste0("./", trial_name, ".Dep.RJMCMC.", i, ".Predict.Log.txt"), Dep.skip, 13, j), envir = .GlobalEnv)
    if (isTRUE(multistate_prediction)) {
      assign("MS.RJ.Prob", compute_posterior(paste0("./", trial_name, ".Multistate.RJMCMC.", i, ".Predict.Log.txt"), MS.skip, 7, j), envir = .GlobalEnv)
    }
  }
  
  if (model %in% c("MCMC", "BOTH")) {
    assign("Ind.MCMC.Prob", compute_posterior(paste0("./", trial_name, ".Ind.MCMC.", i, ".Predict.Log.txt"), Ind.skip, 9, j), envir = .GlobalEnv)
    assign("Dep.MCMC.Prob", compute_posterior(paste0("./", trial_name, ".Dep.MCMC.", i, ".Predict.Log.txt"), Dep.skip, 13, j), envir = .GlobalEnv)
    if (isTRUE(multistate_prediction)) {    
      assign("MS.MCMC.Prob", compute_posterior(paste0("./", trial_name, ".Multistate.MCMC.", i, ".Predict.Log.txt"), MS.skip, 7, j), envir = .GlobalEnv)
    }
  }
}



# This function labels if a prediction was accurate or not for that trial.
calculate_accuracy <- function(rTaxonValueB, Predictive_Probability) {
  # Determine correctness based on probability thresholds
  accuracy <- ifelse((rTaxonValueB == 1 & Predictive_Probability > 0.5) | (rTaxonValueB == 0 & Predictive_Probability < 0.5), 100, 0)

  return(accuracy)
}



# This funciton will be used to calculate LogLoss Values 
LogLoss <- function(rTaxonValueB, Predictive_Probability) {
  # Ensure probabilities are within valid range to avoid log(0) errors
  Predictive_Probability <- pmax(pmin(Predictive_Probability, 0.99999), 0.00001)
  
  # Compute log-loss formula
  LL <- -(rTaxonValueB * log(Predictive_Probability) + 
            (1 - rTaxonValueB) * log(1 - Predictive_Probability))
  
  # Return the mean log loss (works for both vectors and single values)
  return(mean(LL, na.rm = TRUE))
}



# This function will run BayesTraits in R
run_bayestraits <- function(bt_exe, tree_file, data_file, instructions) {
  if (.Platform$OS.type == "windows") {
    cmd <- sprintf('%s "%s" "%s" < "%s"', bt_exe, tree_file, data_file, instructions)
    shell(cmd, wait = TRUE)
  } else {
    cmd <- sprintf('"%s" "%s" "%s" < "%s"', bt_exe, tree_file, data_file, instructions)
    system(cmd, wait = TRUE)
  }
}


# This function will run the rates calculations for BayesTraits in R
run_rates <- function(i, trial_name, RJmodel, multistate_prediction, bt_path = "BayesTraitsV5.exe") {
  # Construct input names
  trial_number <- paste0(trial_name, ".", i)
  treename <- paste0("Trees/Full_tree.", i, ".tre")
  dataname <- paste0(trial_number, ".rates_data.txt")

  # Helper function to call BayesTraits with a given instruction file
  run_bt <- function(settings) {
    if (!file.exists(treename)) stop("Tree file not found: ", treename)
    if (!file.exists(dataname)) stop("Data file not found: ", dataname)
    if (!file.exists(settings)) stop("Settings file not found: ", settings)
    run_bayestraits(bt_path, treename, dataname, settings)
  }

  # Run MCMC models
  if (RJmodel %in% c("MCMC", "BOTH")) {
    if (isTRUE(multistate_prediction)) {
      settings <- paste0(trial_name, ".Multistate.MCMC.", i, ".Rates.In.txt")
      run_bt(settings)
    }

    settings <- paste0(trial_name, ".Ind.MCMC.", i, ".Rates.In.txt")
    run_bt(settings)

    settings <- paste0(trial_name, ".Dep.MCMC.", i, ".Rates.In.txt")
    run_bt(settings)
  }

  # Run RJMCMC models
  if (RJmodel %in% c("RJMCMC", "BOTH")) {
    if (isTRUE(multistate_prediction)) {
      settings <- paste0(trial_name, ".Multistate.RJMCMC.", i, ".Rates.In.txt")
      run_bt(settings)
    }

    settings <- paste0(trial_name, ".Ind.RJMCMC.", i, ".Rates.In.txt")
    run_bt(settings)

    settings <- paste0(trial_name, ".Dep.RJMCMC.", i, ".Rates.In.txt")
    run_bt(settings)
  }
}



# This function will run the prediction calculations for BayesTraits
run_prediction <- function(i, trial_name, RJmodel, multistate_prediction, bt_path = "BayesTraitsV5.exe") {
  # Construct input names
  trial_number <- paste0(trial_name, ".", i)
  treename <- paste0("Trees/Full_tree.", i, ".tre")
  dataname <- paste0(trial_number, ".predict_data.txt")

  # Helper function to call BayesTraits with a given instruction file
  run_bt <- function(settings) {
    if (!file.exists(treename)) stop("Tree file not found: ", treename)
    if (!file.exists(dataname)) stop("Data file not found: ", dataname)
    if (!file.exists(settings)) stop("Settings file not found: ", settings)
    run_bayestraits(bt_path, treename, dataname, settings)
  }

  # Run MCMC models
  if (RJmodel %in% c("MCMC", "BOTH")) {
    if (isTRUE(multistate_prediction)) {
      settings <- paste0(trial_name, ".Multistate.MCMC.", i, ".Predict.In.txt")
      run_bt(settings)
    }

    settings <- paste0(trial_name, ".Ind.MCMC.", i, ".Predict.In.txt")
    run_bt(settings)

    settings <- paste0(trial_name, ".Dep.MCMC.", i, ".Predict.In.txt")
    run_bt(settings)
  }

  # Run RJMCMC models
  if (RJmodel %in% c("RJMCMC", "BOTH")) {
    if (isTRUE(multistate_prediction)) {
      settings <- paste0(trial_name, ".Multistate.RJMCMC.", i, ".Predict.In.txt")
      run_bt(settings)
    }

    settings <- paste0(trial_name, ".Ind.RJMCMC.", i, ".Predict.In.txt")
    run_bt(settings)

    settings <- paste0(trial_name, ".Dep.RJMCMC.", i, ".Predict.In.txt")
    run_bt(settings)
  }
}
