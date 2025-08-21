################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# June 2025

################################################################################
# Running this script will sample n taxa from each dataset to be predicted
# n is the number of sampled taxa, which is specified in the shell script
# It will then store the information in the corresponding folder


# Outputs are stored within the folder that corresponds to the test being run
# (e.g. for single prediction, constant rates, with low, equal q-matrix see /ConstantRates/ER.L/Single)

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

# Get the parameters we need to make the trees
pop_size <- grep("^pop_size=", readLines(bash_file), value = TRUE)
pop_size <- as.numeric(gsub("[^0-9.-]", "", pop_size))

# Get the special cases you may run
variable_rates <- read_logical("variable_rates", bash_file)
multiple_prediction <- read_logical("multiple_prediction", bash_file)
clade_prediction <- read_logical("clade_prediction", bash_file)

# Number of data sets/trials for each combination of settings
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))

# Number of taxa to be sampled
unknown_size <- grep("^unknown_size=", readLines(bash_file), value = TRUE)
unknown_size <- as.numeric(gsub("[^0-9.-]", "", unknown_size))
unknown_size <- pmin((pop_size * .5), pmax(1, unknown_size)) # You should have a value where (1 < n < pop_size/2)



#####
# Install and call the necessary packages
# First, a list of the packages we need
required_packages <- c("ape", "phytools")

# Install any packages that aren't already installed
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed])
}

# Load the libraries
lapply(required_packages, library, character.only = TRUE)



#####
# See if we are running multiple prediction
if (multiple_prediction == TRUE) {
  # Sample taxa for the main trials
  for (type in types) {
    # Create a matrix for this trial's data
    trial_data <- matrix(nrow = 0, ncol = 5)
    colnames(trial_data) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                              "Branch_length")
    
    # Now we need to delve into each type
    for (i in 1:trial_amount) {
      # First, call the full data set
      fullname <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
      TempMat <- read.table(file = fullname, skip = 1)
      
      # Ensure an even sample of A=0 and A=1
      if (i %% 2 == 1) {
        # i is odd: sample where column 2 == 0
        eligible_rows <- which(TempMat[, 2] == 0)
      } else {
        # i is even: sample where column 2 == 1
        eligible_rows <- which(TempMat[, 2] == 1)
      }
      
      # Define ineligible rows (everything else)
      all_rows <- seq_len(nrow(TempMat))
      ineligible_rows <- setdiff(all_rows, eligible_rows)
      
      # Now sample:
      if (length(eligible_rows) >= unknown_size) {
        # Enough eligible rows — sample only from them
        rTaxon <- sample(eligible_rows, unknown_size, replace = FALSE)
      } else {
        # Not enough eligible rows — take all of them
        used_eligible <- eligible_rows
        
        # Check how many more we need
        remainder <- unknown_size - length(used_eligible)
        
        # Sample the rest from ineligible rows
        used_ineligible <- sample(ineligible_rows, remainder, replace = FALSE)
        
        # Combine both parts
        rTaxon <- c(used_eligible, used_ineligible)
      }
      
      # Make sure this is in the correct order then store the character data
      rTaxon <- sort(rTaxon)
      rTaxonValueA <- as.numeric(TempMat[rTaxon, 2])
      rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
      
      # Pull the tree to read the terminal branch length
      treename <- paste0("Trees/Full_tree.", i, ".tre")
      full_tree <- read.nexus(treename)
      
      # Terminal branch length with be calculated as we in put the new data
      # Loop through to store each sampled taxon's data
      for (j in 1:length(rTaxon)) {
        # Store the results in the Trial's Data Table
        trial <- c(i, rTaxon[j], rTaxonValueA[j], rTaxonValueB[j], unknown_branch_length(rTaxon[j], full_tree))
        trial_data <- rbind(trial_data, trial)
      }
    }
    # Name and save the list of unknown taxa
    trial_name <- paste0("ConstantRates/", type, "/Multiple/", type, ".Unknown_info.txt")
    write.table(trial_data, file = trial_name, row.names = F, col.names = T, sep = "\t", quote = F)
  }
  
  
  
  #####
  # Run the random rates next
  # Create a matrix for the random trials' data
  trial_data <- matrix(nrow = 0, ncol = 5)
  colnames(trial_data) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                            "Branch_length")
  
  # Now loop through each trial
  for (i in 1:trial_amount) {
    # First, call the full data set
    fullname <- paste0("Random/Data/Random.", i, ".Full_data.txt")
    TempMat <- read.table(file = fullname, skip = 1)
    
    # Ensure an even sample of A=0 and A=1
    if (i %% 2 == 1) {
      # i is odd: sample where column 2 == 0
      eligible_rows <- which(TempMat[, 2] == 0)
    } else {
      # i is even: sample where column 2 == 1
      eligible_rows <- which(TempMat[, 2] == 1)
    }
    
    # Define ineligible rows (everything else)
    all_rows <- seq_len(nrow(TempMat))
    ineligible_rows <- setdiff(all_rows, eligible_rows)
    
    # Now sample:
    if (length(eligible_rows) >= unknown_size) {
      # Enough eligible rows — sample only from them
      rTaxon <- sample(eligible_rows, unknown_size, replace = FALSE)
    } else {
      # Not enough eligible rows — take all of them
      used_eligible <- eligible_rows
      
      # Check how many more we need
      remainder <- unknown_size - length(used_eligible)
      
      # Sample the rest from ineligible rows
      used_ineligible <- sample(ineligible_rows, remainder, replace = FALSE)
      
      # Combine both parts
      rTaxon <- c(used_eligible, used_ineligible)
    }
    
    # Make sure this is in the correct order then store the character data
    rTaxon <- sort(rTaxon)
    rTaxonValueA <- as.numeric(TempMat[rTaxon, 2])
    rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
    
    # Pull the tree to read the terminal branch length
    treename <- paste0("Trees/Full_tree.", i, ".tre")
    full_tree <- read.nexus(treename)
    
    # Terminal branch length with be calculated as we in put the new data
    # Loop through to store each sampled taxon's data
    for (j in 1:length(rTaxon)) {
      # Store the results in the Trial's Data Table
      trial <- c(i, rTaxon[j], rTaxonValueA[j], rTaxonValueB[j], unknown_branch_length(rTaxon[j], full_tree))
      trial_data <- rbind(trial_data, trial)
    }
  
  # Name and save the list of unknown taxa
  write.table(trial_data, file = "Random/Multiple/Random.Unknown_info.txt", row.names = F, col.names = T, sep = "\t", quote = F)
  
  
  
  #####
  # Sample the variable rates if necessary
  if (variable_rates == TRUE) { 
    for (type in types) {
      # Create a matrix for this trial's data
      trial_data <- matrix(nrow = 0, ncol = 5)
      colnames(trial_data) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                                "Branch_length")
      
      # Now loop through each trial
      for (i in 1:trial_amount) {
        # First, call the full data set
        fullname <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
        TempMat <- read.table(file = fullname, skip = 1)
        
        # Ensure an even sample of A=0 and A=1
        if (i %% 2 == 1) {
          # i is odd: sample where column 2 == 0
          eligible_rows <- which(TempMat[, 2] == 0)
        } else {
          # i is even: sample where column 2 == 1
          eligible_rows <- which(TempMat[, 2] == 1)
        }
        
        # Define ineligible rows (everything else)
        all_rows <- seq_len(nrow(TempMat))
        ineligible_rows <- setdiff(all_rows, eligible_rows)
        
        # Now sample:
        if (length(eligible_rows) >= unknown_size) {
          # Enough eligible rows — sample only from them
          rTaxon <- sample(eligible_rows, unknown_size, replace = FALSE)
        } else {
          # Not enough eligible rows — take all of them
          used_eligible <- eligible_rows
          
          # Check how many more we need
          remainder <- unknown_size - length(used_eligible)
          
          # Sample the rest from ineligible rows
          used_ineligible <- sample(ineligible_rows, remainder, replace = FALSE)
          
          # Combine both parts
          rTaxon <- c(used_eligible, used_ineligible)
        }
        
        # Make sure this is in the correct order then store the character data
        rTaxon <- sort(rTaxon)
        rTaxonValueA <- as.numeric(TempMat[rTaxon, 2])
        rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
        
        # Pull the tree to read the terminal branch length
        treename <- paste0("Trees/Full_tree.", i, ".tre")
        full_tree <- read.nexus(treename)
        
        # Terminal branch length with be calculated as we in put the new data
        # Loop through to store each sampled taxon's data
        for (j in 1:length(rTaxon)) {
          # Store the results in the Trial's Data Table
          trial <- c(i, rTaxon[j], rTaxonValueA[j], rTaxonValueB[j], unknown_branch_length(rTaxon[j], full_tree))
          trial_data <- rbind(trial_data, trial)
        }
      }
      # Name and save the list of unknown taxa
      trial_name <- paste0("VariableRates/", type, "/Multiple/", type, ".Unknown_info.txt")
      write.table(trial_data, file = trial_name, row.names = F, col.names = T, sep = "\t", quote = F)
    }
  }
}



#####
# See if we are running clade prediction
if (clade_prediction == TRUE) {
  # Generate the data for the main trials
  for (type in types) {
    # Create a matrix for this trial's data
    trial_data <- matrix(nrow = 0, ncol = 5)
    colnames(trial_data) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                              "Branch_length")
    
    # Now loop through each trial
    for (i in 1:trial_amount) {
      # Pull the tree to find a random clade
      treename <- paste0("Trees/Full_tree.", i, ".tre")
      full_tree <- read.nexus(treename)
      
      # Use the tree to randomly select a clade of tips to be our rTaxon
      rTaxon <- as.numeric(find_random_clade(full_tree, unknown_size))
      
      # Then, call the full data set
      fullname <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
      TempMat <- read.table(file = fullname, skip = 1)
      
      # Make sure this is in the correct order then store the character data
      rTaxon <- sort(rTaxon)
      rTaxonValueA <- as.numeric(TempMat[rTaxon, 2])
      rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
      
      # Terminal branch length with be calculated as we in put the new data
      # Loop through to store each sampled taxon's data
      for (j in 1:length(rTaxon)) {
        # Store the results in the Trial's Data Table
        trial <- c(i, rTaxon[j], rTaxonValueA[j], rTaxonValueB[j], unknown_branch_length(rTaxon[j], full_tree))
        trial_data <- rbind(trial_data, trial)
      }
    }
    # Name and save the list of unknown taxa
    trial_name <- paste0("ConstantRates/", type, "/Clade/", type, ".Unknown_info.txt")
    write.table(trial_data, file = trial_name, row.names = F, col.names = T, sep = "\t", quote = F)
  }
  
  
  
  #####
  # Run the random rates next
  # Create a matrix for the random trials' data
  trial_data <- matrix(nrow = 0, ncol = 5)
  colnames(trial_data) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                            "Branch_length")
  
  # Now loop through each trial
  for (i in 1:trial_amount) {
    # Pull the tree to find a random clade
    treename <- paste0("Trees/Full_tree.", i, ".tre")
    full_tree <- read.nexus(treename)
    
    # Use the tree to randomly select a clade of tips to be our rTaxon
    rTaxon <- as.numeric(find_random_clade(full_tree, unknown_size))
    
    # Then, call the full data set
    fullname <- paste0("Random/Data/Random.", i, ".Full_data.txt")
    TempMat <- read.table(file = fullname, skip = 1)
    
    # Make sure this is in the correct order then store the character data
    rTaxon <- sort(rTaxon)
    rTaxonValueA <- as.numeric(TempMat[rTaxon, 2])
    rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
    
    # Terminal branch length with be calculated as we in put the new data
    # Loop through to store each sampled taxon's data
    for (j in 1:length(rTaxon)) {
      # Store the results in the Trial's Data Table
      trial <- c(i, rTaxon[j], rTaxonValueA[j], rTaxonValueB[j], unknown_branch_length(rTaxon[j], full_tree))
      trial_data <- rbind(trial_data, trial)
    }
  }
  # Name and save the list of unknown taxa
  write.table(trial_data, file = "Random/Clade/Random.Unknown_info.txt", row.names = F, col.names = T, sep = "\t", quote = F)
  
  
  
  #####
  # Run the variable rates if necessary
  if (variable_rates == TRUE) { 
    for (type in types) {
      # Create a matrix for this trial's data
      trial_data <- matrix(nrow = 0, ncol = 5)
      colnames(trial_data) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                                "Branch_length")
      
      # Now loop through each trial
      for (i in 1:trial_amount) {
        # Pull the tree to find a random clade
        treename <- paste0("Trees/Full_tree.", i, ".tre")
        full_tree <- read.nexus(treename)
        
        # Use the tree to randomly select a clade of tips to be our rTaxon
        rTaxon <- as.numeric(find_random_clade(full_tree, unknown_size))
        
        # Then, call the full data set
        fullname <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
        TempMat <- read.table(file = fullname, skip = 1)
        
        # Make sure this is in the correct order then store the character data
        rTaxon <- sort(rTaxon)
        rTaxonValueA <- as.numeric(TempMat[rTaxon, 2])
        rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
        
        # Terminal branch length with be calculated as we in put the new data
        # Loop through to store each sampled taxon's data
        for (j in 1:length(rTaxon)) {
          # Store the results in the Trial's Data Table
          trial <- c(i, rTaxon[j], rTaxonValueA[j], rTaxonValueB[j], unknown_branch_length(rTaxon[j], full_tree))
          trial_data <- rbind(trial_data, trial)
        }
      }
      # Name and save the list of unknown taxa
      trial_name <- paste0("VariableRates/", type, "/Clade/", type, ".Unknown_info.txt")
      write.table(trial_data, file = trial_name, row.names = F, col.names = T, sep = "\t", quote = F)
    }
  }
}

print("Finished sampling taxa for multi-tip predictions.")
