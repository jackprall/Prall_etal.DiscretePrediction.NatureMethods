################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender, and Chris Organ
# October 2024

################################################################################
# Below defines the functions called in 'TestModels' and 'CompileResults'
# Script uses functions defined in 'DiscreteFunctions'
# Each function has a unique purpose in the simulation, which is further explained below

################################################################################
# These initial functions simulate data and trees for this experiment.
################################################################################

# This function will generate a random tree with "pop_size" number of tips. 
# There are settings to vary the ultrametricity and symmetry of the tree.
generate_tree <- function(tree_size, pop_size, ultrametric = TRUE, symmetric = TRUE) {
  
  # First, we make sure ultrametric is being named correctly
  if (ultrametric == TRUE) {
    complete <- FALSE
  } else {
    complete <- TRUE
  }
  
  # Make sure that tree_size and pop_size are set correctly
  if (pop_size > tree_size) {
    tree_size <- pop_size
  }
  
  # Then, we simulate the tree
  full_tree <- sim.taxa(1, pop_size, m = tree_size, waitsp = "rexp(1.5)", waitex = "rexp(0.5)", symmetric, complete, gsa = TRUE)
  
  # Selecting single tree from list 
  full_tree <- full_tree[[1]]
  
  #Change the labels to be easier to read if necessary
  full_tree$tip.label <- seq_along(full_tree$tip.label)
  
  # Save the tree to the environment
  return(full_tree)
}

# This function rescales and prunes 'full_tree' to set length and tip number
standardize_tree <- function(full_tree, i, pop_size, trial_name) {
  # First, we cut down the tree to the correct tip number
  # Because sim.taxa stops when ultrametric tip count = pop_size
  
  # To start, get number of tips
  tipcount <- Ntip(full_tree)
  
  # Drop excess tips so Ntip = pop_size
  full_tree <- drop.tip(full_tree, sample.int(tipcount, size = (tipcount - pop_size), replace = FALSE ))
  
  # Change the tip labels to be easier to run the tests
  full_tree$tip.label <- seq_along(full_tree$tip.label)
  
  # Next, we have to standardize the tree length
  # To do this, we first get the tree length
  treelength <- as.numeric(sum(full_tree$edge.length))
  
  # Next, we standardize to a total tree length of 700 units
  std.lengths <- lapply(full_tree$edge.length, function(x) x * (700 / treelength))
  
  #Now convert this list back to a numerical vector and replace the old branch lengths
  std.lengths <- as.numeric(unlist(std.lengths))
  full_tree$edge.length <- std.lengths
  
  # Finally, we need to name, save, and assign this tree to use later
  tree_name <- paste0(trial_name, ".", i, ".full_tree.tre")
  write.nexus(full_tree, file = (tree_name))
  return(full_tree)
}

#One function that creates the data and puts it in a temporary matrix.
generate_data <- function(full_tree, tree_scales, pop_size, qmat, Random = FALSE) {
  # First, using the vector in the settings, scale the branch lengths to model transition-rate heterogeneity
  branch_count <- Nedge(full_tree)
  branch_scales <- sample(tree_scales, branch_count, replace = TRUE)
  sim_tree <- full_tree
  for (i in 1:Nedge(full_tree)) {
    sim_tree$edge.length[i] <- full_tree$edge.length[i] * branch_scales[i]
  }
  
  # Run rTraitDisc with the random root state to create a dataset
  TempData <- rTraitDisc(sim_tree, model = qmat, 
                         states = c("00", "01", "10", "11"), 
                         root.value = 1)
  TempMat <- matrix(TempData, nrow = pop_size, ncol = 3)
  colnames(TempMat) <- c("Taxon #", "Trait A", "Trait B")
  
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
  TempMat[, 1] <- 1:pop_size
  
  #Assign the Temporary Matrix to the global environment
  return(TempMat)
}

# This function is used to randomly sample a clade in the following function
find_random_clade <- function(full_tree, unknown_n) {
  # Get all internal nodes
  internal_nodes <- (length(full_tree$tip.label) + 1):max(full_tree$edge)
  
  # Shuffle nodes to randomize search order
  internal_nodes <- sample(internal_nodes)
  
  # Function to find a clade with a given number of tips
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


# This function randomly samples a taxon from a tree and records its number and character states.
sample_taxon <- function(full_tree, TempMat, sample_type = "Random", unknown_n = 1) {
  # First, if a clade is being sampled, we call another function to find an appropriately
  #   sized clade.
  if (sample_type == "Clade") {
    rTaxon <- find_random_clade(full_tree, unknown_n)
  } else {
    # Randomly select the sampled taxa and record the tip numbers
    rTaxon <- sample(1:nrow(TempMat), unknown_n, replace = F)
  }
  
  # Record its character states for traits A and B
  return(as.numeric(rTaxon))
}


#This function removed the sampled taxon from the tree and data set.
edit_tree<-function(i, rTaxon, full_tree, TempMat, trial_name = NULL) {   
  
  #Snip that branch 
  edited_tree <- drop.tip(full_tree, rTaxon)
  
  #Name this new tree
  tree_name <- paste0(trial_name, ".", i, ".edited_tree.tre")
  
  #Create a new "edited" matrix, dropping the sampled taxon's row.
  edited_TempMat <- TempMat
  
  # Assign the first row as row names
  rownames(edited_TempMat) <- edited_TempMat[, 1]
  edited_TempMat <- edited_TempMat[, -1]
  
  # Remove the sampled taxa from the dataset to create the "edited" data
  edited_TempMat <- edited_TempMat[!(rownames(edited_TempMat) %in% rTaxon), , drop = FALSE]
  
  #Name this new data set
  data_name <- paste0(trial_name, ".", i, ".edited_data.txt")
  
  #Save the "edited" tree and write as a .tre file
  write.nexus(edited_tree, file = (tree_name))
  
  #Save the "edited" data as a .txt file
  write.table(edited_TempMat, file = (data_name), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
  
  #Return the new matrix to the environment for later use
  return(edited_TempMat)
}


# This function saves a table of the number of taxa for each character state.
countup <- function(i, edited_TempMat, trial_name) {
  # First, set each count to 0
  s00 <- 0
  s01 <- 0
  s10 <- 0
  s11 <- 0
  
  # Then, this for-loop goes through the data and records the character state of
  # each taxon in the dataset.
  for (j in 1:nrow(edited_TempMat)) {
    if (edited_TempMat[j, 1] == 0 && edited_TempMat[j, 2] == 0) {
      s00 <- s00 + 1
    } else if (edited_TempMat[j, 1] == 0 && edited_TempMat[j, 2] == 1) {
      s01 <- s01 + 1
    } else if (edited_TempMat[j, 1] == 1 && edited_TempMat[j, 2] == 0) {
      s10 <- s10 + 1
    } else if (edited_TempMat[j, 1] == 1 && edited_TempMat[j, 2] == 1) {
      s11 <- s11 + 1
    } else {
    }
  }
  
  # Create a matrix for counts
  counts <- matrix(0, nrow = 2, ncol = 2)
  colnames(counts) <- c("A = 0", "A = 1")
  rownames(counts) <- c("B = 0", "B = 1")
  
  # Fill the matrix with counts
  counts[1, 1] <- s00
  counts[2, 1] <- s01
  counts[1, 2] <- s10
  counts[2, 2] <- s11
  
  # Print the counts. This will help confirm that there is variation mid-simulation.
  print(counts)
  
  #Name the file
  filename <- paste0(trial_name, ".", i, ".Counts.txt")
  
  # Assign the counts to a global variable
  return(counts)
  write.table(counts, file = filename, row.names = T, col.names = T, quote = F)
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


################################################################################
#These functions test the accuracy of standard statistical prediction methods
################################################################################

#This function calculates a Bayesian, non-phylogenetic predictive probability,
# which will be used to compare to the phylogenetic results.
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

#This function will use Bayes Theorem to predict with a "Naive Bayes" model.
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
Multistate_Settings <- function(i, FirstRun = TRUE, RJ = TRUE, trial_name = NULL) {
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
  runname <- paste0(trial_name, ".Multistate.", f, i, c)
  
  # Compile the settings for Multistate analysis
  settings <- paste0("1
2
  
it 1100000
BurnIn 100000
Sample 1000
logFile ", runname, e, "
", a, trial_name, ".Multistate.", f, i, ".ModelFile.bin

Run")
  
  # Save the settings file
  writeLines(text = settings, con = paste0(runname, ".In.txt"))
}

#This function defines the settings for each BayesTraits run.
Bayes_Settings <- function(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, trial_name = NULL) {
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
  runname <- paste0(trial_name, d, f, i, c)

  #Now, we compile the settings
  settings <- paste0(b, "
2

it 1100000
BurnIn 100000
Sample 1000
logFile ", runname, e, "
", a, trial_name, d, f, i, ".ModelFile.bin

Run")

  #Finally, we save the file to be used later
  writeLines(text = settings, con = paste0(runname, ".In.txt"))
}

################################################################################
#These are the functions will perform phylogenetic analyses to predict the trait
################################################################################

#This function re-adds the taxon to the edited data with a "?" character state.
Replace_taxon_data<-function(i, TempMat, rTaxon, trial_name = NULL) {
  # Ensure rTaxon is numeric
  rTaxon <- as.numeric(rTaxon)

  #Replace the sampled taxon's data with a "?"
  Predict_TempMat <- TempMat
  Predict_TempMat[rTaxon, 3] <- "?"

  # Assign the first row (Taxon #) as row names
  rownames(Predict_TempMat) <- Predict_TempMat[, 1]

  # Remove the first column (Taxon #) as it will become row names
  Predict_TempMat <- Predict_TempMat[, -1]

  #Name and save the table.
  dataname <- paste0(trial_name, ".", i, ".predict_data.txt")
  write.table(Predict_TempMat, file = dataname, sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
}

#This function calculates the posterior probability from the output files
Calculate_Post_Prob <- function(i, model, unknown_n, j, trial_name = NULL) {
  
  # Function to calculate posterior probability for a specific column
  compute_posterior <- function(logname, skip_lines, col_offset, j) {
    output <- read.table(logname, header = FALSE, skip = skip_lines, nrows = 1000)
    target_col <- col_offset + j - 1  # Adjust column based on j
    return(mean(output[, target_col]))
  }
  
  MS.skip <- 40 + unknown_n
  Ind.skip <- 44 + unknown_n
  Dep.skip <- 53 + unknown_n
  
  if (model %in% c("RJMCMC", "BOTH")) {
    assign("Ind.RJ.Prob", compute_posterior(paste0("./", trial_name, ".Ind.RJMCMC.", i, ".Predict.Log.txt"), Ind.skip, 9, j), envir = .GlobalEnv)
    assign("Dep.RJ.Prob", compute_posterior(paste0("./", trial_name, ".Dep.RJMCMC.", i, ".Predict.Log.txt"), Dep.skip, 13, j), envir = .GlobalEnv)
    assign("MS.RJ.Prob", compute_posterior(paste0("./", trial_name, ".Multistate.RJMCMC.", i, ".Predict.Log.txt"), MS.skip, 7, j), envir = .GlobalEnv)
  }
  
  if (model %in% c("MCMC", "BOTH")) {
    assign("Ind.MCMC.Prob", compute_posterior(paste0("./", trial_name, ".Ind.MCMC.", i, ".Predict.Log.txt"), Ind.skip, 9, j), envir = .GlobalEnv)
    assign("Dep.MCMC.Prob", compute_posterior(paste0("./", trial_name, ".Dep.MCMC.", i, ".Predict.Log.txt"), Dep.skip, 13, j), envir = .GlobalEnv)
    assign("MS.MCMC.Prob", compute_posterior(paste0("./", trial_name, ".Multistate.MCMC.", i, ".Predict.Log.txt"), MS.skip, 7, j), envir = .GlobalEnv)
  }
}


################################################################################
#These functions are for cleaning up and presenting the results
################################################################################

#This function labels if a prediction was accurate or not for that trial.
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
