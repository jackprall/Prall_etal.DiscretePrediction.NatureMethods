################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# June 2025

################################################################################
# Running the script below will run BayesTraits' phylogenetic analyses
# This may take several days, depending on the power of your machine and the number of cores it has

# Outputs are stored within the folder that corresponds to the test being run
# (e.g. for single prediction, constant rates, with low, equal q-matrix see /ConstantRates/ER.L/Single)

################################################################################



#####
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

# Determine the BayesTraits runs you will run
multistate_prediction <- read_logical("multistate_prediction", bash_file)
RJmodel <- grep("^RJmodel=", readLines(bash_file), value = TRUE)
RJmodel <- gsub(".*=(.*)", "\\1", RJmodel)

# Number of iterations when running as batch job/ on HPC
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))



#####
# Install and call the necessary packages
# First, a list of the packages we need
required_packages <- c("ape", "phytools", "parallel")

# Set up a personal library path
user_lib <- file.path(Sys.getenv("HOME"), "Rlibs")
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}
.libPaths(user_lib)  # Prepend to library search path

# Install any packages that aren't already installed
installed <- required_packages %in% rownames(installed.packages())
if (any(!installed)) {
  install.packages(required_packages[!installed], lib = user_lib, repos = "https://cloud.r-project.org")
}

# Load the libraries
lapply(required_packages, library, character.only = TRUE)



#####
# Set the path to BayesTraits and the trial amount
bt_path <- file.path(getwd(), "BayesTraitsV5.exe")

# Define a vector of trial indices for later
trial_indices <- 1:trial_amount

# Create all combinations
combinations <- expand.grid(
  type = types,
  i = trial_indices,
  stringsAsFactors = FALSE
)

# Number of parallel workers (adjust as needed)
n_cores <- detectCores() - 1  # Use all but one core
if (n_cores > trial_amount) {n_cores <- trial_amount}




#####
# Now we do all of the trials, starting with Constant Rates, Single Prediction
print("Starting BayesTraits runs for Constant Rates, Single Prediction.")
# Start cluster
cl <- makeCluster(detectCores() - 1)

# Export needed functions and objects
clusterExport(cl, c("run_rates", "run_prediction", "run_bayestraits"))
clusterExport(cl, c("RJmodel", "multistate_prediction", "combinations", "bt_path"))

# Run the rates for BayesTraits in parallel to maximize speed
results <- parLapply(
  cl,
  seq_len(nrow(combinations)),
  function(idx) {
    row <- combinations[idx, ]
    trial_name <- paste0("ConstantRates/", row$type, "/Single/", row$type)
    run_rates(
      i = row$i,
      trial_name = trial_name,
      RJmodel = RJmodel,
      multistate_prediction = multistate_prediction,
      bt_path = bt_path
    )
  }
)

# Run the prediction for BayesTraits in parallel to maximize speed
results <- parLapply(
  cl,
  seq_len(nrow(combinations)),
  function(idx) {
    row <- combinations[idx, ]
    trial_name <- paste0("ConstantRates/", row$type, "/Single/", row$type)
    run_prediction(
      i = row$i,
      trial_name = trial_name,
      RJmodel = RJmodel,
      multistate_prediction = multistate_prediction,
      bt_path = bt_path
    )
  }
)
# Stop the cluster to let your computer breathe
stopCluster(cl)
print("Finished BayesTraits runs for Constant Rates, Single Prediction.")



#####
# This is for single, random prediction
print("Starting BayesTraits runs for Random Data, Single Prediction.")
# Start the cluster
cl <- makeCluster(n_cores)

# Export needed functions and objects
clusterExport(cl, c("run_rates", "run_prediction", "run_bayestraits"))
clusterExport(cl, c("RJmodel", "multistate_prediction", "bt_path"))

# Run the rates for BayesTraits in parallel to maximize speed
results <- parLapply(
  cl,
  trial_indices,
  function(i) {
    run_rates(
      i = i,
      trial_name = "Random/Single/Random",
      RJmodel = RJmodel,
      multistate_prediction = multistate_prediction,
      bt_path = bt_path
    )
  }
)

# Run the prediction for BayesTraits in parallel to maximize speed
results <- parLapply(
  cl,
  trial_indices,
  function(i) {
    run_prediction(
      i = i,
      trial_name = "Random/Single/Random",
      RJmodel = RJmodel,
      multistate_prediction = multistate_prediction,
      bt_path = bt_path
    )
  }
)
# Stop cluster
stopCluster(cl)
print("Finished BayesTraits runs for Random Data, Single Prediction.")



#####
# Below here is repeated code for all the optional BayesTraits runs
# This checks and runs variable rates, single prediction if necessary
if (isTRUE(variable_rates)) {
  print("Starting BayesTraits runs for Variable Rates, Single Prediction.")
  # Start cluster
  cl <- makeCluster(detectCores() - 1)
  
  # Export needed functions and objects
  clusterExport(cl, c("run_rates", "run_prediction", "run_bayestraits"))
  clusterExport(cl, c("RJmodel", "multistate_prediction", "combinations", "bt_path"))
  
  # Run the rates for BayesTraits in parallel to maximize speed
  results <- parLapply(
    cl,
    seq_len(nrow(combinations)),
    function(idx) {
      row <- combinations[idx, ]
      trial_name <- paste0("VariableRates/", row$type, "/Single/", row$type)
      run_rates(
        i = row$i,
        trial_name = trial_name,
        RJmodel = RJmodel,
        multistate_prediction = multistate_prediction,
        bt_path = bt_path
      )
    }
  )
  
  # Run the prediction for BayesTraits in parallel to maximize speed
  results <- parLapply(
    cl,
    seq_len(nrow(combinations)),
    function(idx) {
      row <- combinations[idx, ]
      trial_name <- paste0("VariableRates/", row$type, "/Single/", row$type)
      run_prediction(
        i = row$i,
        trial_name = trial_name,
        RJmodel = RJmodel,
        multistate_prediction = multistate_prediction,
        bt_path = bt_path
      )
    }
  )
  # Stop the cluster to let your computer breathe
  stopCluster(cl)
  print("Finished BayesTraits runs for Variable Rates, Single Prediction.")
}


  
#####
# This is for Multiple Prediction
if (isTRUE(multiple_prediction)) {
  # Now we do all of the trials, starting with Constant Rates, Single Prediction
  print("Starting BayesTraits runs for Constant Rates, Multiple Prediction.")
  # Start cluster
  cl <- makeCluster(detectCores() - 1)
  
  # Export needed functions and objects
  clusterExport(cl, c("run_rates", "run_prediction", "run_bayestraits"))
  clusterExport(cl, c("RJmodel", "multistate_prediction", "combinations", "bt_path"))
  
  # Run the rates for BayesTraits in parallel to maximize speed
  results <- parLapply(
    cl,
    seq_len(nrow(combinations)),
    function(idx) {
      row <- combinations[idx, ]
      trial_name <- paste0("ConstantRates/", row$type, "/Multiple/", row$type)
      run_rates(
        i = row$i,
        trial_name = trial_name,
        RJmodel = RJmodel,
        multistate_prediction = multistate_prediction,
        bt_path = bt_path
      )
    }
  )
  
  # Run the prediction for BayesTraits in parallel to maximize speed
  results <- parLapply(
    cl,
    seq_len(nrow(combinations)),
    function(idx) {
      row <- combinations[idx, ]
      trial_name <- paste0("ConstantRates/", row$type, "/Multiple/", row$type)
      run_prediction(
        i = row$i,
        trial_name = trial_name,
        RJmodel = RJmodel,
        multistate_prediction = multistate_prediction,
        bt_path = bt_path
      )
    }
  )
  # Stop the cluster to let your computer breathe
  stopCluster(cl)
  print("Finished BayesTraits runs for Constant Rates, Multiple Prediction.")
  
  
  
  #####
  # This is for multiple, random prediction
  print("Starting BayesTraits runs for Random Data, Multiple Prediction.")
  # Start the cluster
  cl <- makeCluster(n_cores)
  
  # Export needed functions and objects
  clusterExport(cl, c("run_rates", "run_prediction", "run_bayestraits"))
  clusterExport(cl, c("RJmodel", "multistate_prediction", "bt_path"))
  
  # Run the rates for BayesTraits in parallel to maximize speed
  results <- parLapply(
    cl,
    trial_indices,
    function(i) {
      run_rates(
        i = i,
        trial_name = "Random/Multiple/Random",
        RJmodel = RJmodel,
        multistate_prediction = multistate_prediction,
        bt_path = bt_path
      )
    }
  )
  
  # Run the prediction for BayesTraits in parallel to maximize speed
  results <- parLapply(
    cl,
    trial_indices,
    function(i) {
      run_prediction(
        i = i,
        trial_name = "Random/Multiple/Random",
        RJmodel = RJmodel,
        multistate_prediction = multistate_prediction,
        bt_path = bt_path
      )
    }
  )
  # Stop cluster
  stopCluster(cl)
  print("Finished BayesTraits runs for Random Data, Multiple Prediction.")
  
  
  
  #####
  # This checks and runs variable rates, multiple prediction if necessary
  if (isTRUE(variable_rates)) {
    print("Starting BayesTraits runs for Variable Rates, Multiple Prediction.")
    # Start cluster
    cl <- makeCluster(detectCores() - 1)
    
    # Export needed functions and objects
    clusterExport(cl, c("run_rates", "run_prediction", "run_bayestraits"))
    clusterExport(cl, c("RJmodel", "multistate_prediction", "combinations", "bt_path"))
    
    # Run the rates for BayesTraits in parallel to maximize speed
    results <- parLapply(
      cl,
      seq_len(nrow(combinations)),
      function(idx) {
        row <- combinations[idx, ]
        trial_name <- paste0("VariableRates/", row$type, "/Multiple/", row$type)
        run_rates(
          i = row$i,
          trial_name = trial_name,
          RJmodel = RJmodel,
          multistate_prediction = multistate_prediction,
          bt_path = bt_path
        )
      }
    )
    
    # Run the prediction for BayesTraits in parallel to maximize speed
    results <- parLapply(
      cl,
      seq_len(nrow(combinations)),
      function(idx) {
        row <- combinations[idx, ]
        trial_name <- paste0("VariableRates/", row$type, "/Multiple/", row$type)
        run_prediction(
          i = row$i,
          trial_name = trial_name,
          RJmodel = RJmodel,
          multistate_prediction = multistate_prediction,
          bt_path = bt_path
        )
      }
    )
    # Stop the cluster to let your computer breathe
    stopCluster(cl)
    print("Finished BayesTraits runs for Variable Rates, Multiple Prediction.")
  }
}



#####
# This is for Clade Prediction
if (isTRUE(clade_prediction)) {
  # Now we do all of the trials, starting with Constant Rates, Clade Prediction
  print("Starting BayesTraits runs for Constant Rates, Clade Prediction.")
  # Start cluster
  cl <- makeCluster(detectCores() - 1)
  
  # Export needed functions and objects
  clusterExport(cl, c("run_rates", "run_prediction", "run_bayestraits"))
  clusterExport(cl, c("RJmodel", "multistate_prediction", "combinations", "bt_path"))
  
  # Run the rates for BayesTraits in parallel to maximize speed
  results <- parLapply(
    cl,
    seq_len(nrow(combinations)),
    function(idx) {
      row <- combinations[idx, ]
      trial_name <- paste0("ConstantRates/", row$type, "/Clade/", row$type)
      run_rates(
        i = row$i,
        trial_name = trial_name,
        RJmodel = RJmodel,
        multistate_prediction = multistate_prediction,
        bt_path = bt_path
      )
    }
  )
  
  # Run the prediction for BayesTraits in parallel to maximize speed
  results <- parLapply(
    cl,
    seq_len(nrow(combinations)),
    function(idx) {
      row <- combinations[idx, ]
      trial_name <- paste0("ConstantRates/", row$type, "/Clade/", row$type)
      run_prediction(
        i = row$i,
        trial_name = trial_name,
        RJmodel = RJmodel,
        multistate_prediction = multistate_prediction,
        bt_path = bt_path
      )
    }
  )
  # Stop the cluster to let your computer breathe
  stopCluster(cl)
  print("Finished BayesTraits runs for Constant Rates, Clade Prediction.")
  
  
  
  #####
  # This is for clade, random prediction
  print("Starting BayesTraits runs for Random Data, Clade Prediction.")
  # Start the cluster
  cl <- makeCluster(n_cores)
  
  # Export needed functions and objects
  clusterExport(cl, c("run_rates", "run_prediction", "run_bayestraits"))
  clusterExport(cl, c("RJmodel", "multistate_prediction", "bt_path"))
  
  # Run the rates for BayesTraits in parallel to maximize speed
  results <- parLapply(
    cl,
    trial_indices,
    function(i) {
      run_rates(
        i = i,
        trial_name = "Random/Clade/Random",
        RJmodel = RJmodel,
        multistate_prediction = multistate_prediction,
        bt_path = bt_path
      )
    }
  )
  
  # Run the prediction for BayesTraits in parallel to maximize speed
  results <- parLapply(
    cl,
    trial_indices,
    function(i) {
      run_prediction(
        i = i,
        trial_name = "Random/Clade/Random",
        RJmodel = RJmodel,
        multistate_prediction = multistate_prediction,
        bt_path = bt_path
      )
    }
  )
  # Stop cluster
  stopCluster(cl)
  print("Finished BayesTraits runs for Random Data, Clade Prediction.")
  
  
  
  #####
  # This checks and runs variable rates, clade prediction if necessary
  if (isTRUE(variable_rates)) {
    print("Starting BayesTraits runs for Variable Rates, Clade Prediction.")
    # Start cluster
    cl <- makeCluster(detectCores() - 1)
    
    # Export needed functions and objects
    clusterExport(cl, c("run_rates", "run_prediction", "run_bayestraits"))
    clusterExport(cl, c("RJmodel", "multistate_prediction", "combinations", "bt_path"))
    
    # Run the rates for BayesTraits in parallel to maximize speed
    results <- parLapply(
      cl,
      seq_len(nrow(combinations)),
      function(idx) {
        row <- combinations[idx, ]
        trial_name <- paste0("VariableRates/", row$type, "/Clade/", row$type)
        run_rates(
          i = row$i,
          trial_name = trial_name,
          RJmodel = RJmodel,
          multistate_prediction = multistate_prediction,
          bt_path = bt_path
        )
      }
    )
    
    # Run the prediction for BayesTraits in parallel to maximize speed
    results <- parLapply(
      cl,
      seq_len(nrow(combinations)),
      function(idx) {
        row <- combinations[idx, ]
        trial_name <- paste0("VariableRates/", row$type, "/Clade/", row$type)
        run_prediction(
          i = row$i,
          trial_name = trial_name,
          RJmodel = RJmodel,
          multistate_prediction = multistate_prediction,
          bt_path = bt_path
        )
      }
    )
    # Stop the cluster to let your computer breathe
    stopCluster(cl)
    print("Finished BayesTraits runs for Variable Rates, Clade Prediction.")
  }
}

print("Finished all BayesTraits runs")
