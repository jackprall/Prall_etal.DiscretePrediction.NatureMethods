################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# June 2025

################################################################################
# Running the script below will create the data and instruction files for BayesTraits


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

# Get the special cases you may run
variable_rates <- read_logical("variable_rates", bash_file)
multiple_prediction <- read_logical("multiple_prediction", bash_file)
clade_prediction <- read_logical("clade_prediction", bash_file)
multistate_prediction <- read_logical("multistate_prediction", bash_file)

# Number of iterations when running as batch job/ on HPC
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))

# Get the MCMC type
RJmodel <- grep("^RJmodel=", readLines(bash_file), value = TRUE)
RJmodel <- gsub(".*=(.*)", "\\1", RJmodel)



### Start the loop for Single Tip Prediction
for (type in types) {
  # Clarify the paths to the data and outputs
  output_path <- paste0("ConstantRates/", type, "/Single/", type)
  unknown_name <- paste0(output_path, ".Unknown_info.txt")
  unknown_info <- read.table(unknown_name, skip = 1)
  colnames(unknown_info) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", "Branch_length")  
  
  # Loop through the individual trials  
  for (i in 1:trial_amount) {
    # Call the data for this trial
    dataname <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
    full_data <- read.table(dataname, skip = 1)
    colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
    
    # Call the unknown taxon
    rTaxon <- as.numeric(unknown_info[i, 2])
    
    # Create the data tables to use in BayesTraits
    edit_data(i, rTaxon, full_data, output_path)
    
    
    # Now we will write the necessary instructions for BayesTraits
    # These lines write the instructions for non-RJ BayesTraits runs, if necessary.
    if (RJmodel %in% c("MCMC", "BOTH")) {
      # These lines write the multistate prediction instructions, if necessary.
      if (multistate_prediction == TRUE) {
        Multistate_Settings(i, FirstRun = TRUE, RJ = FALSE, output_path)
        Multistate_Settings(i, FirstRun = FALSE, RJ = FALSE, output_path)
      }
      #First, the rate calculation for the independent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
      #Second, the rate calculation for the dependent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
      # Next, the prediction settings for the independent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
      #Last, the prediction settings for the dependent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
    }
    
    # These lines write the RJ BayesTraits runs.
    if (RJmodel %in% c("RJMCMC", "BOTH")) {
      # These lines write the multistate prediction instructions, if necessary.
      if (multistate_prediction == TRUE) {
        Multistate_Settings(i, FirstRun = TRUE, RJ = TRUE, output_path)
        Multistate_Settings(i, FirstRun = FALSE, RJ = TRUE, output_path)
      }     
      #First, the rate calculation for the independent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
      #Second, the rate calculation for the dependent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
      # Next, the prediction settings for the independent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
      #Last, the prediction settings for the dependent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
    }
  }
  


  ##### 
  # Now, we double back if we are testing variable rates
  if (variable_rates == TRUE) {
    # Clarify the paths to the data and outputs
    output_path <- paste0("VariableRates/", type, "/Single/", type)
    unknown_name <- paste0(output_path, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1)
    colnames(unknown_info) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", "Branch_length")  
    
    # Loop through the individual trials  
    for (i in 1:trial_amount) {
      # Call the data for this trial
      dataname <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
      full_data <- read.table(dataname, skip = 1)
      colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
      
      # Call the unknown taxon
      rTaxon <- as.numeric(unknown_info[i, 2])
      
      # Create the data tables to use in BayesTraits
      edit_data(i, rTaxon, full_data, output_path)
      
      
      # Now we will write the necessary instructions for BayesTraits
      # These lines write the instructions for non-RJ BayesTraits runs, if necessary.
      if (RJmodel %in% c("MCMC", "BOTH")) {
        # These lines write the multistate prediction instructions, if necessary.
        if (multistate_prediction == TRUE) {
          Multistate_Settings(i, FirstRun = TRUE, RJ = FALSE, output_path)
          Multistate_Settings(i, FirstRun = FALSE, RJ = FALSE, output_path)
        }
        #First, the rate calculation for the independent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
        #Second, the rate calculation for the dependent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
        # Next, the prediction settings for the independent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
        #Last, the prediction settings for the dependent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
      }
      
      # These lines write the RJ BayesTraits runs.
      if (RJmodel %in% c("RJMCMC", "BOTH")) {
        # These lines write the multistate prediction instructions, if necessary.
        if (multistate_prediction == TRUE) {
          Multistate_Settings(i, FirstRun = TRUE, RJ = TRUE, output_path)
          Multistate_Settings(i, FirstRun = FALSE, RJ = TRUE, output_path)
        }     
        #First, the rate calculation for the independent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
        #Second, the rate calculation for the dependent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
        # Next, the prediction settings for the independent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
        #Last, the prediction settings for the dependent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
      }
    }
  }
}



##### 
# Start the loop for the Random Trials
# Clarify the paths to the data and outputs
output_path <- "Random/Single/Random"
unknown_name <- paste0(output_path, ".Unknown_info.txt")
unknown_info <- read.table(unknown_name, skip = 1)
colnames(unknown_info) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", "Branch_length")  

# Loop through the individual trials  
for (i in 1:trial_amount) {
  # Call the data for this trial
  dataname <- paste0("Random/Data/Random.", i, ".Full_data.txt")
  full_data <- read.table(dataname, skip = 1)
  colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
  
  # Call the unknown taxon
  rTaxon <- as.numeric(unknown_info[i, 2])
  
  # Create the data tables to use in BayesTraits
  edit_data(i, rTaxon, full_data, output_path)
  
  
  # Now we will write the necessary instructions for BayesTraits
  # These lines write the instructions for non-RJ BayesTraits runs, if necessary.
  if (RJmodel %in% c("MCMC", "BOTH")) {
    # These lines write the multistate prediction instructions, if necessary.
    if (multistate_prediction == TRUE) {
      Multistate_Settings(i, FirstRun = TRUE, RJ = FALSE, output_path)
      Multistate_Settings(i, FirstRun = FALSE, RJ = FALSE, output_path)
    }
    #First, the rate calculation for the independent run.
    Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
    #Second, the rate calculation for the dependent run.
    Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
    # Next, the prediction settings for the independent run.
    Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
    #Last, the prediction settings for the dependent run.
    Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
  }
  
  # These lines write the RJ BayesTraits runs.
  if (RJmodel %in% c("RJMCMC", "BOTH")) {
    # These lines write the multistate prediction instructions, if necessary.
    if (multistate_prediction == TRUE) {
      Multistate_Settings(i, FirstRun = TRUE, RJ = TRUE, output_path)
      Multistate_Settings(i, FirstRun = FALSE, RJ = TRUE, output_path)
    }     
    #First, the rate calculation for the independent run.
    Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
    #Second, the rate calculation for the dependent run.
    Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
    # Next, the prediction settings for the independent run.
    Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
    #Last, the prediction settings for the dependent run.
    Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
  }
}



#####
# Start the loop for Multiple Tip Prediction
if (multiple_prediction == TRUE) {
  for (type in types) {
    # Clarify the paths to the data and outputs
    output_path <- paste0("ConstantRates/", type, "/Multiple/", type)
    unknown_name <- paste0(output_path, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1)
    colnames(unknown_info) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", "Branch_length")  
    
    # Loop through the individual trials  
    for (i in 1:trial_amount) {
      # Call the data for this trial
      dataname <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
      full_data <- read.table(dataname, skip = 1)
      colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
      
      # Call the unknown taxon by subsetting the rows where column 1 matches i
      matching_rows <- which(unknown_info[, 1] == i)
      
      # Extract the values from column 2 in those rows
      rTaxon <- unknown_info[matching_rows, 2]
      
      # Create the data tables to use in BayesTraits
      edit_data(i, rTaxon, full_data, output_path)
      
      
      # Now we will write the necessary instructions for BayesTraits
      # These lines write the instructions for non-RJ BayesTraits runs, if necessary.
      if (RJmodel %in% c("MCMC", "BOTH")) {
        # These lines write the multistate prediction instructions, if necessary.
        if (multistate_prediction == TRUE) {
          Multistate_Settings(i, FirstRun = TRUE, RJ = FALSE, output_path)
          Multistate_Settings(i, FirstRun = FALSE, RJ = FALSE, output_path)
        }
        #First, the rate calculation for the independent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
        #Second, the rate calculation for the dependent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
        # Next, the prediction settings for the independent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
        #Last, the prediction settings for the dependent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
      }
      
      # These lines write the RJ BayesTraits runs.
      if (RJmodel %in% c("RJMCMC", "BOTH")) {
        # These lines write the multistate prediction instructions, if necessary.
        if (multistate_prediction == TRUE) {
          Multistate_Settings(i, FirstRun = TRUE, RJ = TRUE, output_path)
          Multistate_Settings(i, FirstRun = FALSE, RJ = TRUE, output_path)
        }     
        #First, the rate calculation for the independent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
        #Second, the rate calculation for the dependent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
        # Next, the prediction settings for the independent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
        #Last, the prediction settings for the dependent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
      }
    }
    
    
    #####
    # Now, we double back if we are testing variable rates
    if (variable_rates == TRUE) {
      # Clarify the paths to the data and outputs
      output_path <- paste0("VariableRates/", type, "/Multiple/", type)
      unknown_name <- paste0(output_path, ".Unknown_info.txt")
      unknown_info <- read.table(unknown_name, skip = 1)
      colnames(unknown_info) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", "Branch_length")  
      
      # Loop through the individual trials  
      for (i in 1:trial_amount) {
        # Call the data for this trial
        dataname <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
        full_data <- read.table(dataname, skip = 1)
        colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
        
        # Call the unknown taxon by subsetting the rows where column 1 matches i
        matching_rows <- which(unknown_info[, 1] == i)
        
        # Extract the values from column 2 in those rows
        rTaxon <- unknown_info[matching_rows, 2]
        
        # Create the data tables to use in BayesTraits
        edit_data(i, rTaxon, full_data, output_path)
        
        
        # Now we will write the necessary instructions for BayesTraits
        # These lines write the instructions for non-RJ BayesTraits runs, if necessary.
        if (RJmodel %in% c("MCMC", "BOTH")) {
          # These lines write the multistate prediction instructions, if necessary.
          if (multistate_prediction == TRUE) {
            Multistate_Settings(i, FirstRun = TRUE, RJ = FALSE, output_path)
            Multistate_Settings(i, FirstRun = FALSE, RJ = FALSE, output_path)
          }
          #First, the rate calculation for the independent run.
          Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
          #Second, the rate calculation for the dependent run.
          Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
          # Next, the prediction settings for the independent run.
          Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
          #Last, the prediction settings for the dependent run.
          Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
        }
        
        # These lines write the RJ BayesTraits runs.
        if (RJmodel %in% c("RJMCMC", "BOTH")) {
          # These lines write the multistate prediction instructions, if necessary.
          if (multistate_prediction == TRUE) {
            Multistate_Settings(i, FirstRun = TRUE, RJ = TRUE, output_path)
            Multistate_Settings(i, FirstRun = FALSE, RJ = TRUE, output_path)
          }     
          #First, the rate calculation for the independent run.
          Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
          #Second, the rate calculation for the dependent run.
          Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
          # Next, the prediction settings for the independent run.
          Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
          #Last, the prediction settings for the dependent run.
          Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
        }
      }
    }
  }
}



#####
# Start the loop for Clade Prediction
if (clade_prediction == TRUE) {
  for (type in types) {
    # Clarify the paths to the data and outputs
    output_path <- paste0("ConstantRates/", type, "/Clade/", type)
    unknown_name <- paste0(output_path, ".Unknown_info.txt")
    unknown_info <- read.table(unknown_name, skip = 1)
    colnames(unknown_info) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", "Branch_length")  
    
    # Loop through the individual trials  
    for (i in 1:trial_amount) {
      # Call the data for this trial
      dataname <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
      full_data <- read.table(dataname, skip = 1)
      colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
      
      # Call the unknown taxon by subsetting the rows where column 1 matches i
      matching_rows <- which(unknown_info[, 1] == i)
      
      # Extract the values from column 2 in those rows
      rTaxon <- unknown_info[matching_rows, 2]
      
      # Create the data tables to use in BayesTraits
      edit_data(i, rTaxon, full_data, output_path)
      
      
      # Now we will write the necessary instructions for BayesTraits
      # These lines write the instructions for non-RJ BayesTraits runs, if necessary.
      if (RJmodel %in% c("MCMC", "BOTH")) {
        # These lines write the multistate prediction instructions, if necessary.
        if (multistate_prediction == TRUE) {
          Multistate_Settings(i, FirstRun = TRUE, RJ = FALSE, output_path)
          Multistate_Settings(i, FirstRun = FALSE, RJ = FALSE, output_path)
        }
        #First, the rate calculation for the independent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
        #Second, the rate calculation for the dependent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
        # Next, the prediction settings for the independent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
        #Last, the prediction settings for the dependent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
      }
      
      # These lines write the RJ BayesTraits runs.
      if (RJmodel %in% c("RJMCMC", "BOTH")) {
        # These lines write the multistate prediction instructions, if necessary.
        if (multistate_prediction == TRUE) {
          Multistate_Settings(i, FirstRun = TRUE, RJ = TRUE, output_path)
          Multistate_Settings(i, FirstRun = FALSE, RJ = TRUE, output_path)
        }     
        #First, the rate calculation for the independent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
        #Second, the rate calculation for the dependent run.
        Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
        # Next, the prediction settings for the independent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
        #Last, the prediction settings for the dependent run.
        Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
      }
    }
    
    

    #####
    # Now, we double back if we are testing variable rates
    if (variable_rates == TRUE) {
      # Clarify the paths to the data and outputs
      output_path <- paste0("VariableRates/", type, "/Clade/", type)
      unknown_name <- paste0(output_path, ".Unknown_info.txt")
      unknown_info <- read.table(unknown_name, skip = 1)
      colnames(unknown_info) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", "Branch_length")  
      
      # Loop through the individual trials  
      for (i in 1:trial_amount) {
        # Call the data for this trial
        dataname <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
        full_data <- read.table(dataname, skip = 1)
        colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
        
        # Call the unknown taxon by subsetting the rows where column 1 matches i
        matching_rows <- which(unknown_info[, 1] == i)
        
        # Extract the values from column 2 in those rows
        rTaxon <- unknown_info[matching_rows, 2]
        
        # Create the data tables to use in BayesTraits
        edit_data(i, rTaxon, full_data, output_path)
        
        
        # Now we will write the necessary instructions for BayesTraits
        # These lines write the instructions for non-RJ BayesTraits runs, if necessary.
        if (RJmodel %in% c("MCMC", "BOTH")) {
          # These lines write the multistate prediction instructions, if necessary.
          if (multistate_prediction == TRUE) {
            Multistate_Settings(i, FirstRun = TRUE, RJ = FALSE, output_path)
            Multistate_Settings(i, FirstRun = FALSE, RJ = FALSE, output_path)
          }
          #First, the rate calculation for the independent run.
          Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
          #Second, the rate calculation for the dependent run.
          Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
          # Next, the prediction settings for the independent run.
          Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
          #Last, the prediction settings for the dependent run.
          Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
        }
        
        # These lines write the RJ BayesTraits runs.
        if (RJmodel %in% c("RJMCMC", "BOTH")) {
          # These lines write the multistate prediction instructions, if necessary.
          if (multistate_prediction == TRUE) {
            Multistate_Settings(i, FirstRun = TRUE, RJ = TRUE, output_path)
            Multistate_Settings(i, FirstRun = FALSE, RJ = TRUE, output_path)
          }     
          #First, the rate calculation for the independent run.
          Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
          #Second, the rate calculation for the dependent run.
          Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
          # Next, the prediction settings for the independent run.
          Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
          #Last, the prediction settings for the dependent run.
          Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
        }
      }
    }
  }
}



##### 
# Reloop if you are doing Multiple Prediction for Random Data
if (multiple_prediction == TRUE) {
  # Clarify the paths to the data and outputs
  output_path <- "Random/Multiple/Random"
  unknown_name <- paste0(output_path, ".Unknown_info.txt")
  unknown_info <- read.table(unknown_name, skip = 1)
  colnames(unknown_info) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", "Branch_length")  
  
  # Loop through the individual trials  
  for (i in 1:trial_amount) {
    # Call the data for this trial
    dataname <- paste0("Random/Data/Random.", i, ".Full_data.txt")
    full_data <- read.table(dataname, skip = 1)
    colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
    
    # Call the unknown taxon by subsetting the rows where column 1 matches i
    matching_rows <- which(unknown_info[, 1] == i)
    
    # Extract the values from column 2 in those rows
    rTaxon <- unknown_info[matching_rows, 2]
    
    # Create the data tables to use in BayesTraits
    edit_data(i, rTaxon, full_data, output_path)
    
    
    # Now we will write the necessary instructions for BayesTraits
    # These lines write the instructions for non-RJ BayesTraits runs, if necessary.
    if (RJmodel %in% c("MCMC", "BOTH")) {
      # These lines write the multistate prediction instructions, if necessary.
      if (multistate_prediction == TRUE) {
        Multistate_Settings(i, FirstRun = TRUE, RJ = FALSE, output_path)
        Multistate_Settings(i, FirstRun = FALSE, RJ = FALSE, output_path)
      }
      #First, the rate calculation for the independent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
      #Second, the rate calculation for the dependent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
      # Next, the prediction settings for the independent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
      #Last, the prediction settings for the dependent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
    }
    
    # These lines write the RJ BayesTraits runs.
    if (RJmodel %in% c("RJMCMC", "BOTH")) {
      # These lines write the multistate prediction instructions, if necessary.
      if (multistate_prediction == TRUE) {
        Multistate_Settings(i, FirstRun = TRUE, RJ = TRUE, output_path)
        Multistate_Settings(i, FirstRun = FALSE, RJ = TRUE, output_path)
      }     
      #First, the rate calculation for the independent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
      #Second, the rate calculation for the dependent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
      # Next, the prediction settings for the independent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
      #Last, the prediction settings for the dependent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
    }
  }
}



##### 
# Reloop if you are doing Clade Prediction for Random Data
if (clade_prediction == TRUE) {
  # Clarify the paths to the data and outputs
  output_path <- "Random/Clade/Random"
  unknown_name <- paste0(output_path, ".Unknown_info.txt")
  unknown_info <- read.table(unknown_name, skip = 1)
  colnames(unknown_info) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", "Branch_length")  
  
  # Loop through the individual trials  
  for (i in 1:trial_amount) {
    # Call the data for this trial
    dataname <- paste0("Random/Data/Random.", i, ".Full_data.txt")
    full_data <- read.table(dataname, skip = 1)
    colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
    
    # Call the unknown taxon by subsetting the rows where column 1 matches i
    matching_rows <- which(unknown_info[, 1] == i)
    
    # Extract the values from column 2 in those rows
    rTaxon <- unknown_info[matching_rows, 2]
    
    # Create the data tables to use in BayesTraits
    edit_data(i, rTaxon, full_data, output_path)
    
    
    # Now we will write the necessary instructions for BayesTraits
    # These lines write the instructions for non-RJ BayesTraits runs, if necessary.
    if (RJmodel %in% c("MCMC", "BOTH")) {
      # These lines write the multistate prediction instructions, if necessary.
      if (multistate_prediction == TRUE) {
        Multistate_Settings(i, FirstRun = TRUE, RJ = FALSE, output_path)
        Multistate_Settings(i, FirstRun = FALSE, RJ = FALSE, output_path)
      }
      #First, the rate calculation for the independent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
      #Second, the rate calculation for the dependent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
      # Next, the prediction settings for the independent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = FALSE, output_path)
      #Last, the prediction settings for the dependent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = FALSE, output_path)
    }
    
    # These lines write the RJ BayesTraits runs.
    if (RJmodel %in% c("RJMCMC", "BOTH")) {
      # These lines write the multistate prediction instructions, if necessary.
      if (multistate_prediction == TRUE) {
        Multistate_Settings(i, FirstRun = TRUE, RJ = TRUE, output_path)
        Multistate_Settings(i, FirstRun = FALSE, RJ = TRUE, output_path)
      }     
      #First, the rate calculation for the independent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
      #Second, the rate calculation for the dependent run.
      Bayes_Settings(i, FirstRun = TRUE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
      # Next, the prediction settings for the independent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = TRUE, RJ = TRUE, output_path)
      #Last, the prediction settings for the dependent run.
      Bayes_Settings(i, FirstRun = FALSE, IndependentCharacters = FALSE, RJ = TRUE, output_path)
    }
  }
}

print("Finished creating all necessary input files for BayesTraits.")
