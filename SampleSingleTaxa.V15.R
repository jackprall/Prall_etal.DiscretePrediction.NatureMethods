################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# June 2025

################################################################################
# Running this script will sample one taxon from each dataset to be predicted
# It will then store the information in the corresponding folder
# This should run in under 2 minutes for 1000 iterations


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

# Number of data sets/trials for each combination of settings
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))

# Keep track of the original directory
original_dir <- getwd()



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
    
    # Now sample only from eligible rows
    if (length(eligible_rows) > 0) {
      rTaxon <- sample(eligible_rows, 1)
    } else {
      rTaxon <- sample(1:nrow(TempMat), 1, replace = F)
    }
    
    # Record the character state(s) for rTaxon
    rTaxonValueA <- as.numeric(TempMat[rTaxon, 2])
    rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
    
    # Pull the tree to read the terminal branch length
    treename <- paste0("Trees/Full_tree.", i, ".tre")
    full_tree <- read.nexus(treename)
    
    # This function will record the branch length leading to the unknown taxon.
    BL <- unknown_branch_length(rTaxon, full_tree)
    
    
    # Store the results in the Trial's Data Table
    trial <- c(i, rTaxon, rTaxonValueA, rTaxonValueB, BL)
    trial_data <- rbind(trial_data, trial)
  }
  # Name and save the list of unknown taxa
  trial_name <- paste0("ConstantRates/", type, "/Single/", type, ".Unknown_info.txt")
  write.table(trial_data, file = trial_name, row.names = F, col.names = T, sep = "\t", quote = F)
}




#####
# Create a matrix for the random trials' unknown information
trial_data <- matrix(nrow = 0, ncol = 5)
colnames(trial_data) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                          "Branch_length")

# Now we need to delve into each type
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
  
  # Now sample only from eligible rows
  if (length(eligible_rows) > 0) {
    rTaxon <- sample(eligible_rows, 1)
  } else {
    rTaxon <- sample(1:nrow(TempMat), 1, replace = F)
  }
  
  # Record the character state(s) for rTaxon
  rTaxonValueA <- as.numeric(TempMat[rTaxon, 2])
  rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
  
  # Pull the tree to read the terminal branch length
  treename <- paste0("Trees/Full_tree.", i, ".tre")
  full_tree <- read.nexus(treename)
  
  # This function will record the branch length leading to the unknown taxon.
  BL <- unknown_branch_length(rTaxon, full_tree)
  
  
  #Store the results in the Trial's Data Table
  trial <- c(i, rTaxon, rTaxonValueA, rTaxonValueB, BL)
  trial_data <- rbind(trial_data, trial)
}
# Name and save the list of unknown taxa
write.table(trial_data, file = "Random/Single/Random.Unknown_info.txt", row.names = F, col.names = T, sep = "\t", quote = F)



##### 
# Run the variable rates if necessary
if (variable_rates == TRUE) { 
  for (type in types) {
    # Create a matrix for this trial's data
    trial_data <- matrix(nrow = 0, ncol = 5)
    colnames(trial_data) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                              "Branch_length")
    
    # Now we need to delve into each type
    for (i in 1:trial_amount) {
      # First, call the full data set
      fullname <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
      TempMat <- read.table(file = fullname, skip = 1)
      
      # Ensure an even sample of A=0 and A=1
      if (i %% 2 == 1) {
        # i is odd: sample where column 3 == 0
        eligible_rows <- which(TempMat[, 2] == 0)
      } else {
        # i is even: sample where column 3 == 1
        eligible_rows <- which(TempMat[, 2] == 1)
      }
      
      # Now sample only from eligible rows
      if (length(eligible_rows) > 0) {
        rTaxon <- sample(eligible_rows, 1)
      } else {
        rTaxon <- sample(1:nrow(TempMat), 1, replace = F)
      }
      
      # Record the character state(s) for rTaxon
      rTaxonValueA <- as.numeric(TempMat[rTaxon, 2])
      rTaxonValueB <- as.numeric(TempMat[rTaxon, 3])
      
      # Pull the tree to read the terminal branch length
      treename <- paste0("Trees/Full_tree.", i, ".tre")
      full_tree <- read.nexus(treename)
      
      # This function will record the branch length leading to the unknown taxon.
      BL <- unknown_branch_length(rTaxon, full_tree)
      
      
      #Store the results in the Trial's Data Table
      trial <- c(i, rTaxon, rTaxonValueA, rTaxonValueB, BL)
      trial_data <- rbind(trial_data, trial)
    }
    # Name and save the list of unknown taxa
    trial_name <- paste0("VariableRates/", type, "/Single/", type, ".Unknown_info.txt")
    write.table(trial_data, file = trial_name, row.names = F, col.names = T, sep = "\t", quote = F)
  }
}

print("Finished sampling taxa for single-tip predictions.")
