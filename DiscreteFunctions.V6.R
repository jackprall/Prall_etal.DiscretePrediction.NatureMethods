################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender, and Chris Organ
# October 2024

################################################################################
# Below defines the functions called in 'TestModels' and 'CompileResults'
# Script uses functions defined in 'DiscreteFunctions'
# Each function has a unique purpose in the simulation, which is further explained below

################################################################################
# These initial functions simulate data and trees for this experiement.
################################################################################

# This function will generate a random tree with "pop_size" number of tips. 
# There are settings to vary the ultrametricity and symmetry of the tree.
generate_tree <- function(pop_size, ultrametric = FALSE, symmetric = TRUE, trial_name = NULL) {

  # First, we make sure ultrametric is being named correctly
  if (ultrametric == TRUE) {
    complete <- FALSE
  } else {
    complete <- TRUE
  }
  
  # Then, we simulate the tree
  full_tree <- sim.taxa(1, pop_size, m = pop_size, waitsp = "rexp(1.5)", waitex = "rexp(0.5)", symmetric, complete)

  # Selecting single tree from list 
  full_tree <- full_tree[[1]]
  
  #Change the labels to be easier to read if necessary
  full_tree$tip.label <- seq_along(full_tree$tip.label)

  # Save the tree to the environment
  assign("full_tree", full_tree, envir = .GlobalEnv)
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
  assign("full_tree", full_tree, envir = .GlobalEnv)
}

#One function that creates the data and puts it in a temporary matrix.
generate_data <- function(full_tree, pop_size, qmat, Random = FALSE) {
  # Run rTraitDisc with the random root state to create a dataset
  TempData <- rTraitDisc(full_tree, model = qmat, 
                         states = c("00", "01", "10", "11"), 
                         root.value = 1)
  TempMat <- matrix(TempData, nrow = pop_size, ncol = 3)
  colnames(TempMat) <- c("Taxon #", "Trait A", "Trait B")

  #Rework the matrix to make it a little easier to analyze
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
  assign("TempMat", TempMat, envir = .GlobalEnv)
}

#This function randomly samples a taxon from a tree and records its number and character states.
sample_taxon <- function(TempMat, trial_name) {
  #Randomly select one branch and record the tip number
  rTaxon <- sample(1:nrow(TempMat), 1)
  assign("rTaxon", rTaxon, envir = .GlobalEnv)

  #Record its character states for traits A and B
  rTaxonValueA <- TempMat[rTaxon, 2]
  assign("rTaxonValueA", rTaxonValueA, envir = .GlobalEnv)
  rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
  assign("rTaxonValueB", rTaxonValueB, envir = .GlobalEnv)
}

#This function removed the sampled taxon from the tree and data set.
edit_tree<-function(i, rTaxon, full_tree, TempMat, trial_name = NULL) {   

  #Snip that branch 
  edited_tree <- drop.tip(full_tree, rTaxon)

  #Name this new tree
  tree_name <- paste0(trial_name, ".", i, ".edited_tree.tre")

  #Create a new "edited" matrix, dropping the sampled taxon's row.
  edited_TempMat <- TempMat
  edited_TempMat <- edited_TempMat[-rTaxon, , drop = FALSE]

  # Assign the first row as row names
  rownames(edited_TempMat) <- edited_TempMat[, 1]

  # Remove the first column (Taxon #) as it will become row names
  edited_TempMat <- edited_TempMat[, -1]

  #Assign the new matrix to the environment
  assign("edited_TempMat", edited_TempMat, envir = .GlobalEnv)

  #Name this new data set
  data_name <- paste0(trial_name, ".", i, ".edited_data.txt")

  #Save the "edited" tree and write as a .tre file
  write.nexus(edited_tree, file = (tree_name))

  #Save the "edited" data as a .txt file
  write.table(edited_TempMat, file = (data_name), sep = "\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
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
  assign("Counts", counts, envir = .GlobalEnv)
  write.table(counts, file = filename, row.names = T, col.names = T, quote = F)
}

# This function will record the terminal branch length of the unknown taxon.
unknown_branch_length <- function(rTaxon, full_tree) {
  # First, identify the correct edge
  leaf_index <- which(full_tree$tip.label == as.character(rTaxon))

  # Then, pull out the branch length
  branch_length <- full_tree$edge.length[which(full_tree$edge[,2] == leaf_index)]

  # Finally, assign the branch length to be recorded with the rest of the trial
  assign("BL", branch_length, envir = .GlobalEnv)
}


################################################################################
#These functions test the accuracy of standard statistical prediction methods
################################################################################

#This function calculates a bayesian, non-phylogenetic predictive probability,
# which will be used to compare to the phylogenetic results.
Beta_Bin_predict <- function(edited_TempMat) {

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
  # a bayesian posterior predictive distribution.
  PiDist <- rbinom(1000,1, prob = pi.values)

  #Get a frequency of 1s from the rbinom, and that's the posterior probability
  frequency <- (sum(PiDist)/1000)

  assign("BB.Prob", frequency, envir = .GlobalEnv)
}

#This function will use Bayes Theorem to predict with a "Naive Bayes" model.
Naive_Bayes_predict <- function(Counts, rTaxonValueA) {
  # This is Bayes Theorum both in general and for our case
  # posterior = prior * likelihood / marginal
  # p(B=1|A) = p(B=1) * p(A|B=1) / P(A)
  # However, because we have the counts for each tree, this can be simplified to P(B=1|A)

  #Next, the likelihood and marginal, based on rTaxonValueB
  if (rTaxonValueA == 0) {
    posterior <- Counts[2,1] / (Counts[1,1] + Counts[2,1])
  } else {
    posterior <- Counts[2,2] / (Counts[1,2] + Counts[2,2])
  }

  #Finally, we assign this to the global environment
  assign("Naive.Prob", posterior, envir = .GlobalEnv)
}

#This function defines the settings for each BayesTraits run.
Bayes_Settings <- function(i, FirstRun = TRUE, IndependentCharacters = TRUE, trial_name = NULL) {

  #First, the one major change between the first and second run
  #a tells BT to either save or call a model, and c is for naming the files
  if (FirstRun == TRUE) {
    a <- "SaveModels "
    c <- ".Rates"
    e <- "
    rjhp exp 0 5"
  } else {
    a <- "LoadModels "
    c <- ".Predict"
    e <- ""
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
  runname <- paste0(trial_name, d, i, c)

  #Now, we compile the settings
  settings <- paste0(b, "
2

it 1100000
BurnIn 100000
Sample 1000
logFile ", runname, e, "
", a, trial_name, d, i, ".ModelFile.bin

Run")

  #Finally, we save the file to be used later
  writeLines(text = settings, con = paste0(runname, ".In.txt"))
}

################################################################################
#These are the functions will perform phylogenetic analyses to predict the trait
################################################################################

#This funciton re-adds the taxon to the edited data with a "?" character state.
Replace_taxon_data<-function(i, TempMat, rTaxon, trial_name = NULL) {

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
Calculate_Post_Prob <- function(i, rTaxon, trial_name = NULL, IndependentCharacters = TRUE) {

  #These will call and edit the output file to gather the data we need.
  if (IndependentCharacters == TRUE) {
    #This selects the correct file
    logname <- paste0("./", trial_name, ".Ind.", i, ".Predict.Log.txt")
    #This calls the log file into R
    output <- read.table(logname, header = FALSE, skip = 45, nrows = 1000)
    #This takes the average of the predictive distribution
    Post.Prob <- mean(output[, 9])
    #This names the object correctly for other functions.
    name <- "Ind.Post.Prob"
  } else {
    logname <- paste0("./", trial_name, ".Dep.", i, ".Predict.Log.txt")
    output <- read.table(logname, header = FALSE, skip = 53, nrows = 1000)
    Post.Prob <- mean(output[, 13])
    name <- "Dep.Post.Prob"
  }

  #Finally, we assign it to the global environment
  assign(name, Post.Prob, envir = .GlobalEnv)
}

################################################################################
#These functions are for cleaning up and presenting the results
################################################################################

#This function labels if a method was accurate or not for that trial.
calculate_NP_accuracy <- function(rTaxonValueB, BB.Prob, Naive.Prob) {
  if (rTaxonValueB == 1) {
    ifelse(BB.Prob > 0.5, BB.acc <- "Correct", BB.acc <- "Incorrect")
    ifelse(Naive.Prob > 0.5, Naive.acc <- "Correct", Naive.acc <- "Incorrect")
  } else {
    ifelse(BB.Prob < 0.5, BB.acc <- "Correct", BB.acc <- "Incorrect")
    ifelse(Naive.Prob < 0.5, Naive.acc <- "Correct", Naive.acc <- "Incorrect")
  }

  assign("BB.acc", BB.acc, envir = .GlobalEnv)
  assign("Naive.acc", Naive.acc, envir = .GlobalEnv)
}

#This function labels if a method was accurate or not for that trial.
calculate_PH_accuracy <- function(rTaxonValueB, Ind.Post.Prob, Dep.Post.Prob) {
  if (rTaxonValueB == 1) {
    ifelse(Ind.Post.Prob > 0.5, Ind.acc <- "Correct", Ind.acc <- "Incorrect")
    ifelse(Dep.Post.Prob > 0.5, Dep.acc <- "Correct", Dep.acc <- "Incorrect")
  } else {
    ifelse(Ind.Post.Prob < 0.5, Ind.acc <- "Correct", Ind.acc <- "Incorrect")
    ifelse(Dep.Post.Prob < 0.5, Dep.acc <- "Correct", Dep.acc <- "Incorrect")
  }

  assign("Ind.acc", Ind.acc, envir = .GlobalEnv)
  assign("Dep.acc", Dep.acc, envir = .GlobalEnv)
}