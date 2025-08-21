################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# June 2025

################################################################################
# This will produce empty results tables that later scripts will fill in
# One table produced for every combination of data simulation inputs

# The output tables can be found in the 'Results' folder.
# They are named after the data simulation and prediction parameters 
# (e.g. ER.Tip is equal rates and single tip prediction)

################################################################################



# Start by calling the necessary parameters and settings
# Call the bash file and any necessary function
bash_file <- Sys.glob("Discrete_Simulation.V*")
version_str <- sub(".*(V[0-9]+).*", "\\1", bash_file)
source(paste0("Scripts/DiscreteFunctions.", version_str, ".R"))

# Get the types and trials arrays
types <- grep("^types=", readLines(bash_file), value = TRUE)
types <- gsub("types=\\(|\\)", "", types)
types <- gsub("\"", "", types)
types <- strsplit(types, " ")[[1]]

# Get the special cases you may run
variable_rates <- read_logical("variable_rates", bash_file)
multiple_prediction <- read_logical("multiple_prediction", bash_file)
clade_prediction <- read_logical("clade_prediction", bash_file)
multistate_prediction <- read_logical("multistate_prediction", bash_file)

# Get the MCMC type
RJmodel <- grep("^RJmodel=", readLines(bash_file), value = TRUE)
RJmodel <- gsub(".*=(.*)", "\\1", RJmodel)

# Number of data sets/trials for each combination of settings
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))



#####
# Determine the necessary size dataframe, beginning with what is constant
tests <- c("Beta_Binom", "Naive_Bayes")
test_size <- 2

# Expand for the BayesTraits trials
if (RJmodel == "MCMC") {
  BTruns <- c("MCMC_Ind", "MCMC_Dep")
  BTsize <- 2
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("MCMC_Multi", BTruns)
    BTsize <- BTsize + 1
  }
}

if (RJmodel == "RJMCMC") {
  BTruns <- c("RJ_Ind", "RJ_Dep")
  BTsize <- 2
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("RJ_Multi", BTruns)
    BTsize <- BTsize + 1
  }
}

if (RJmodel == "BOTH") {
  BTruns <- c("MCMC_Ind", "MCMC_Dep", "RJ_Ind", "RJ_Dep")
  BTsize <- 4
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("MCMC_Multi", "MCMC_Ind", "MCMC_Dep", "RJ_Multi", "RJ_Ind", "RJ_Dep")
    BTsize <- BTsize + 2
  }
}

# Determine the total number of models being tested and their names
tests <- c(tests, BTruns)
test_size <- test_size + BTsize

# Now, we separate out each into three columns
prob_labels <- paste0(tests, "_Prob")
acc_labels <- paste0(tests, "_acc")
LL_labels <- paste0(tests, "_LL")

# Lastly, we want to include the state frequency counts
state_labels <- c("00", "01", "10", "11", "Avg_Sister_A", "Avg_Sister_B")

# We readjust our new column names and column count to match
new_col_names <- c(prob_labels, acc_labels, LL_labels, state_labels)
new_col_number <- (test_size * 3) + 6



#####
# Now we create the new matrices for the data
# First, for constant rates, single tip prediction
for (type in types) {
  # First, pull the unknown matrix
  unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
  unknown_name <- paste0("ConstantRates/", type, "/Single/", type, ".Unknown_info.txt")
  unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
  row_num <- nrow(unknown_info)
  
  # Next, create the new section of the results matrix
  matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
  colnames(matrix_b) <- new_col_names
  
  # Combine the unknown information with the new space for results
  results <- cbind(unknown_info, matrix_b)
  
  # Save this information for later
  OO_index <- ncol(results) - 5
  Ol_index <- ncol(results) - 4
  lO_index <- ncol(results) - 3
  ll_index <- ncol(results) - 2
  
  # Fill in the counts data that has already by calculated
  for (i in 1:trial_amount) {
    # Find the rows that match the trial number
    trial_i <- which(results[, 1] == i)
    
    # Call the count data
    count_name <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Counts.txt")
    counts <- read.table(count_name, skip = 1)
    
    # Fill in the count data to the main results folder
    results[trial_i, OO_index] <- counts[1, 2]
    results[trial_i, Ol_index] <- counts[2, 2]
    results[trial_i, lO_index] <- counts[1, 3]
    results[trial_i, ll_index] <- counts[2, 3]
  }
  
  # Finally, write the matrix somewhere it can be accessed later
  results_name <- paste0("Results/ConstantRates/Single/", type, ".Single.ResultsFull.txt")
  write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
}



#####
# Now we do this for the Random Data
# First, pull the unknown matrix
unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
unknown_name <- "Random/Single/Random.Unknown_info.txt"
unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
row_num <- nrow(unknown_info)

# Next, create the new section of the results matrix
matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
colnames(matrix_b) <- new_col_names

# Combine the unknown information with the new space for results
results <- cbind(unknown_info, matrix_b)

# Save this information for later
OO_index <- ncol(results) - 5
Ol_index <- ncol(results) - 4
lO_index <- ncol(results) - 3
ll_index <- ncol(results) - 2

# Fill in the counts data that has already by calculated
for (i in 1:trial_amount) {
  # Find the rows that match the trial number
  trial_i <- which(results[, 1] == i)
  
  # Call the count data
  count_name <- paste0("Random/Data/Random.", i, ".Counts.txt")
  counts <- read.table(count_name, skip = 1)
  
  # Fill in the count data to the main results folder
  results[trial_i, OO_index] <- counts[1, 2]
  results[trial_i, Ol_index] <- counts[2, 2]
  results[trial_i, lO_index] <- counts[1, 3]
  results[trial_i, ll_index] <- counts[2, 3]
}

# Finally, write the matrix somewhere it can be accessed later
results_name <- "Results/Random/Random.Single.ResultsFull.txt"
write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)



#####
# Below are the reruns for all of the possible runs
# First of these is single prediction, variable rates
if (variable_rates == TRUE) {
  for (type in types) {
    # First, pull the unknown matrix
    unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
    unknown_name <- paste0("VariableRates/", type, "/Single/", type, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
    row_num <- nrow(unknown_info)
    
    # Next, create the new section of the results matrix
    matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
    colnames(matrix_b) <- new_col_names
    
    # Combine the unknown information with the new space for results
    results <- cbind(unknown_info, matrix_b)
    
    # Save this information for later
    OO_index <- ncol(results) - 5
    Ol_index <- ncol(results) - 4
    lO_index <- ncol(results) - 3
    ll_index <- ncol(results) - 2
    
    # Fill in the counts data that has already by calculated
    for (i in 1:trial_amount) {
      # Find the rows that match the trial number
      trial_i <- which(results[, 1] == i)
      
      # Call the count data
      count_name <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Counts.txt")
      counts <- read.table(count_name, skip = 1)
      
      # Fill in the count data to the main results folder
      results[trial_i, OO_index] <- counts[1, 2]
      results[trial_i, Ol_index] <- counts[2, 2]
      results[trial_i, lO_index] <- counts[1, 3]
      results[trial_i, ll_index] <- counts[2, 3]
    }
    
    # Finally, write the matrix somewhere it can be accessed later
    results_name <- paste0("Results/VariableRates/Single/", type, ".Single.ResultsFull.txt")
    write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}



#####
# Now, check for multiple prediction for constant rates
if (multiple_prediction == TRUE) {
  for (type in types) {
    # First, pull the unknown matrix
    unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
    unknown_name <- paste0("ConstantRates/", type, "/Multiple/", type, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
    row_num <- nrow(unknown_info)
    
    # Next, create the new section of the results matrix
    matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
    colnames(matrix_b) <- new_col_names
    
    # Combine the unknown information with the new space for results
    results <- cbind(unknown_info, matrix_b)
    
    # Save this information for later
    OO_index <- ncol(results) - 5
    Ol_index <- ncol(results) - 4
    lO_index <- ncol(results) - 3
    ll_index <- ncol(results) - 2
    
    # Fill in the counts data that has already by calculated
    for (i in 1:trial_amount) {
      # Find the rows that match the trial number
      trial_i <- which(results[, 1] == i)
      
      # Call the count data
      count_name <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Counts.txt")
      counts <- read.table(count_name, skip = 1)
      
      # Fill in the count data to the main results folder
      results[trial_i, OO_index] <- counts[1, 2]
      results[trial_i, Ol_index] <- counts[2, 2]
      results[trial_i, lO_index] <- counts[1, 3]
      results[trial_i, ll_index] <- counts[2, 3]
    }
    
    # Finally, write the matrix somewhere it can be accessed later
    results_name <- paste0("Results/ConstantRates/Multiple/", type, ".Multiple.ResultsFull.txt")
    write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}



#####
# Check multiple prediction, variable rates
if (variable_rates == TRUE && multiple_prediction == TRUE) {
  for (type in types) {
    # First, pull the unknown matrix
    unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
    unknown_name <- paste0("VariableRates/", type, "/Multiple/", type, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
    row_num <- nrow(unknown_info)
    
    # Next, create the new section of the results matrix
    matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
    colnames(matrix_b) <- new_col_names
    
    # Combine the unknown information with the new space for results
    results <- cbind(unknown_info, matrix_b)
    
    # Save this information for later
    OO_index <- ncol(results) - 5
    Ol_index <- ncol(results) - 4
    lO_index <- ncol(results) - 3
    ll_index <- ncol(results) - 2
    
    # Fill in the counts data that has already by calculated
    for (i in 1:trial_amount) {
      # Find the rows that match the trial number
      trial_i <- which(results[, 1] == i)
      
      # Call the count data
      count_name <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Counts.txt")
      counts <- read.table(count_name, skip = 1)
      
      # Fill in the count data to the main results folder
      results[trial_i, OO_index] <- counts[1, 2]
      results[trial_i, Ol_index] <- counts[2, 2]
      results[trial_i, lO_index] <- counts[1, 3]
      results[trial_i, ll_index] <- counts[2, 3]
    }
    
    # Finally, write the matrix somewhere it can be accessed later
    results_name <- paste0("Results/VariableRates/Multiple/", type, ".Multiple.ResultsFull.txt")
    write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}



#####
# Next, check for clade predictions for constant rates
if (clade_prediction == TRUE) {
  for (type in types) {
    # First, pull the unknown matrix
    unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
    unknown_name <- paste0("ConstantRates/", type, "/Clade/", type, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
    row_num <- nrow(unknown_info)
    
    # Next, create the new section of the results matrix
    matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
    colnames(matrix_b) <- new_col_names
    
    # Combine the unknown information with the new space for results
    results <- cbind(unknown_info, matrix_b)
    
    # Save this information for later
    OO_index <- ncol(results) - 5
    Ol_index <- ncol(results) - 4
    lO_index <- ncol(results) - 3
    ll_index <- ncol(results) - 2
    
    # Fill in the counts data that has already by calculated
    for (i in 1:trial_amount) {
      # Find the rows that match the trial number
      trial_i <- which(results[, 1] == i)
      
      # Call the count data
      count_name <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Counts.txt")
      counts <- read.table(count_name, skip = 1)
      
      # Fill in the count data to the main results folder
      results[trial_i, OO_index] <- counts[1, 2]
      results[trial_i, Ol_index] <- counts[2, 2]
      results[trial_i, lO_index] <- counts[1, 3]
      results[trial_i, ll_index] <- counts[2, 3]
    }
    
    # Finally, write the matrix somewhere it can be accessed later
    results_name <- paste0("Results/ConstantRates/Clade/", type, ".Clade.ResultsFull.txt")
    write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}



#####
# Check clade prediction, variable rates
if (variable_rates == TRUE && clade_prediction == TRUE) {
  for (type in types) {
    # First, pull the unknown matrix
    unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
    unknown_name <- paste0("VariableRates/", type, "/Clade/", type, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
    row_num <- nrow(unknown_info)
    
    # Next, create the new section of the results matrix
    matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
    colnames(matrix_b) <- new_col_names
    
    # Combine the unknown information with the new space for results
    results <- cbind(unknown_info, matrix_b)
    
    # Save this information for later
    OO_index <- ncol(results) - 5
    Ol_index <- ncol(results) - 4
    lO_index <- ncol(results) - 3
    ll_index <- ncol(results) - 2
    
    # Fill in the counts data that has already by calculated
    for (i in 1:trial_amount) {
      # Find the rows that match the trial number
      trial_i <- which(results[, 1] == i)
      
      # Call the count data
      count_name <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Counts.txt")
      counts <- read.table(count_name, skip = 1)
      
      # Fill in the count data to the main results folder
      results[trial_i, OO_index] <- counts[1, 2]
      results[trial_i, Ol_index] <- counts[2, 2]
      results[trial_i, lO_index] <- counts[1, 3]
      results[trial_i, ll_index] <- counts[2, 3]
    }
    
    # Finally, write the matrix somewhere it can be accessed later
    results_name <- paste0("Results/VariableRates/Clade/", type, ".Clade.ResultsFull.txt")
    write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
  }
}



#####
# Next, check if we are running randomly generated multiple prediction
if (multiple_prediction == TRUE) {
  unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
  unknown_name <- "Random/Multiple/Random.Unknown_info.txt"
  unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
  row_num <- nrow(unknown_info)
  
  # Next, create the new section of the results matrix
  matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
  colnames(matrix_b) <- new_col_names
  
  # Combine the unknown information with the new space for results
  results <- cbind(unknown_info, matrix_b)
  
  # Save this information for later
  OO_index <- ncol(results) - 5
  Ol_index <- ncol(results) - 4
  lO_index <- ncol(results) - 3
  ll_index <- ncol(results) - 2
  
  # Fill in the counts data that has already by calculated
  for (i in 1:trial_amount) {
    # Find the rows that match the trial number
    trial_i <- which(results[, 1] == i)
    
    # Call the count data
    count_name <- paste0("Random/Data/Random.", i, ".Counts.txt")
    counts <- read.table(count_name, skip = 1)
    
    # Fill in the count data to the main results folder
    results[trial_i, OO_index] <- counts[1, 2]
    results[trial_i, Ol_index] <- counts[2, 2]
    results[trial_i, lO_index] <- counts[1, 3]
    results[trial_i, ll_index] <- counts[2, 3]
  }
  
  # Finally, write the matrix somewhere it can be accessed later
  results_name <- "Results/Random/Random.Multiple.ResultsFull.txt"
  write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
}



#####
# Finally, check if we are running randomly generated clade prediction
if (clade_prediction == TRUE) {
  unknown_labels <- c("Trial", "rTaxon", "Trait_A", "Trait_B", "Terminal_Branch_Length")
  unknown_name <- "Random/Clade/Random.Unknown_info.txt"
  unknown_info <- read.table(unknown_name, skip = 1, col.names = unknown_labels)
  row_num <- nrow(unknown_info)
  
  # Next, create the new section of the results matrix
  matrix_b <- matrix(data=NA, nrow = row_num, ncol = new_col_number)
  colnames(matrix_b) <- new_col_names
  
  # Combine the unknown information with the new space for results
  results <- cbind(unknown_info, matrix_b)
  
  # Save this information for later
  OO_index <- ncol(results) - 5
  Ol_index <- ncol(results) - 4
  lO_index <- ncol(results) - 3
  ll_index <- ncol(results) - 2
  
  # Fill in the counts data that has already by calculated
  for (i in 1:trial_amount) {
    # Find the rows that match the trial number
    trial_i <- which(results[, 1] == i)
    
    # Call the count data
    count_name <- paste0("Random/Data/Random.", i, ".Counts.txt")
    counts <- read.table(count_name, skip = 1)
    
    # Fill in the count data to the main results folder
    results[trial_i, OO_index] <- counts[1, 2]
    results[trial_i, Ol_index] <- counts[2, 2]
    results[trial_i, lO_index] <- counts[1, 3]
    results[trial_i, ll_index] <- counts[2, 3]
  }
  
  # Finally, write the matrix somewhere it can be accessed later
  results_name <- "Results/Random/Random.Clade.ResultsFull.txt"
  write.table(results, file = results_name, quote = F, sep = "\t", row.names = F, col.names = T)
}

print("Finished generating Results tables. They can be found in the Results folder.")
print("Results tables are currently mostly blank. Following scripts will fill them in.")
