### This script is for pulling the BayesTraits results into the results file
################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# June 2025

################################################################################
# Running the script below will test multiple phylogenetic and non-phylogenetic prediction models. 
# Calls the files creates by 'TestModels'
# Runs BayesTraits phylogenetic analyses
# Calls and runs 'TestModels' and 'CompileResults' in R
# Written for the Tempest Research Cluster at Montana State University

# Outputs will include any files made by BayesTraits and three Results files
# Results are stored within a 'Results' folder found within the folder for each matrix.

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



### Install and call the necessary packages
# First, a list of the packages we need
required_packages <- c("ape", "phytools")

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
# First, determine the BayesTraits trials which were run by the RJmodel
if (RJmodel == "MCMC") {
  BTruns <- c("MCMC_Ind", "MCMC_Dep")
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("MCMC_Multi", BTruns)
  }
}

if (RJmodel == "RJMCMC") {
  BTruns <- c("RJ_Ind", "RJ_Dep")
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("RJ_Multi", BTruns)
  }
}

if (RJmodel == "BOTH") {
  BTruns <- c("MCMC_Ind", "MCMC_Dep", "RJ_Ind", "RJ_Dep")
  if (isTRUE(multistate_prediction)) {
    BTruns <- c("MCMC_Multi", "MCMC_Ind", "MCMC_Dep", "RJ_Multi", "RJ_Ind", "RJ_Dep")
  }
}

# Now, we separate out each into three columns
prob_labels <- paste0(BTruns, "_Prob")
acc_labels <- paste0(BTruns, "_acc")
LL_labels <- paste0(BTruns, "_LL")
Results_labels <- c(prob_labels, acc_labels, LL_labels)



#####
# Start the loop that will read the trials, beginning with constant rates and single prediction.
for (type in types) {
  print(paste0("Compiling the BayesTraits results from the single-prediction, ", type, " tests with constant rates."))
  # Call the necessary paths
  logfile_path <- paste0("ConstantRates/", type, "/Single/", type)
  results_name <- paste0("Results/ConstantRates/Single/", type, ".Single.ResultsFull.txt")
  
  # Pull the results into R
  colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
  results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
  
  # Start the for-loop that will run each prediction and fill it in
  for (i in 1:trial_amount) {
    # First, call rTaxon by subsetting the rows where column 1 matches i
    matching_rows <- which(results[, 1] == i)
    rTaxon <- results[matching_rows, 2]
    rTaxonValueA <- results[matching_rows, 3]
    rTaxonValueB <- results[matching_rows, 4]

### This loop will cover any situations where rTaxon is more than one value
    for (j in 1:length(rTaxon)) {
      
      # Gather the information from 
      skiplength <- length(rTaxon)
      Calculate_Post_Prob(i, RJmodel, multistate_prediction, skiplength, j, trial_name = logfile_path)
      
      # Use the RJmodel to pull the correct data
      if (RJmodel %in% c("MCMC", "BOTH")) {
        Ind.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Ind.MCMC.Prob)
        Dep.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Dep.MCMC.Prob)
        
        Ind.MCMC.LL  <- LogLoss(rTaxonValueB[j], Ind.MCMC.Prob)
        Dep.MCMC.LL  <- LogLoss(rTaxonValueB[j], Dep.MCMC.Prob)    
        
        if (isTRUE(multistate_prediction)) {
          MS.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], MS.MCMC.Prob)
          MS.MCMC.LL  <- LogLoss(rTaxonValueB[j], MS.MCMC.Prob)
          
        }
      }
      
      if (RJmodel %in% c("RJMCMC", "BOTH")) {
        Ind.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Ind.RJ.Prob)
        Dep.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Dep.RJ.Prob)
        
        Ind.RJ.LL  <- LogLoss(rTaxonValueB[j], Ind.RJ.Prob)
        Dep.RJ.LL  <- LogLoss(rTaxonValueB[j], Dep.RJ.Prob)    
        
        if (isTRUE(multistate_prediction)) {
          MS.RJ.acc <- calculate_accuracy(rTaxonValueB[j], MS.RJ.Prob)
          MS.RJ.LL  <- LogLoss(rTaxonValueB[j], MS.RJ.Prob)
        }
      }
      
      
      # Finally assign the variables to the final table based on the RJ model
      if (RJmodel %in% c("MCMC", "BOTH")) {
        MCMC_prob <- c(Ind.MCMC.Prob, Dep.MCMC.Prob)
        MCMC_acc  <- c(Ind.MCMC.acc, Dep.MCMC.acc)
        MCMC_LL   <- c(Ind.MCMC.LL,  Dep.MCMC.LL)
        
        if (isTRUE(multistate_prediction)) {
          MCMC_prob <- c(MS.MCMC.Prob, MCMC_prob)
          MCMC_acc  <- c(MS.MCMC.acc,  MCMC_acc)
          MCMC_LL   <- c(MS.MCMC.LL,   MCMC_LL)
        }
      }
      
      if (RJmodel %in% c("RJMCMC", "BOTH")) {
        RJ_prob <- c(Ind.RJ.Prob, Dep.RJ.Prob)
        RJ_acc  <- c(Ind.RJ.acc,  Dep.RJ.acc)
        RJ_LL   <- c(Ind.RJ.LL,   Dep.RJ.LL)
        
        if (isTRUE(multistate_prediction)) {
          RJ_prob <- c(MS.RJ.Prob, RJ_prob)
          RJ_acc  <- c(MS.RJ.acc,  RJ_acc)
          RJ_LL   <- c(MS.RJ.LL,   RJ_LL)
        }
      }
      
      if (RJmodel == "MCMC") {
        result_vector <- c(MCMC_prob, MCMC_acc, MCMC_LL)
      } else if (RJmodel == "RJMCMC") {
        result_vector <- c(RJ_prob, RJ_acc, RJ_LL)
      } else if (RJmodel == "BOTH") {
        result_vector <- c(MCMC_prob, RJ_prob, MCMC_acc, RJ_acc, MCMC_LL, RJ_LL)
      }
      
      # Assign results
      results[matching_rows[j], Results_labels] <- result_vector
      
    }
  }
  write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t") 
}



### Now we will do the same for the randomly generated data
print(paste0("Compiling the BayesTraits results from the single-prediction, tests with random data."))
# Call the necessary paths
logfile_path <- "Random/Single/Random"
results_name <- "Results/Random/Random.Single.ResultsFull.txt"

# Pull the results into R
colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)

# Start the for-loop that will run each prediction and fill it in
for (i in 1:trial_amount) {
  # First, call rTaxon by subsetting the rows where column 1 matches i
  matching_rows <- which(results[, 1] == i)
  rTaxon <- results[matching_rows, 2]
  rTaxonValueA <- results[matching_rows, 3]
  rTaxonValueB <- results[matching_rows, 4]
  
  ### This loop will cover any situations where rTaxon is more than one value
  for (j in 1:length(rTaxon)) {
    
    # Gather the information from 
    skiplength <- length(rTaxon)
    Calculate_Post_Prob(i, RJmodel, multistate_prediction, skiplength, j, trial_name = logfile_path)
    
    # Use the RJmodel to pull the correct data
    if (RJmodel %in% c("MCMC", "BOTH")) {
      Ind.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Ind.MCMC.Prob)
      Dep.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Dep.MCMC.Prob)
      
      Ind.MCMC.LL  <- LogLoss(rTaxonValueB[j], Ind.MCMC.Prob)
      Dep.MCMC.LL  <- LogLoss(rTaxonValueB[j], Dep.MCMC.Prob)    
      
      if (isTRUE(multistate_prediction)) {
        MS.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], MS.MCMC.Prob)
        MS.MCMC.LL  <- LogLoss(rTaxonValueB[j], MS.MCMC.Prob)
        
      }
    }
    
    if (RJmodel %in% c("RJMCMC", "BOTH")) {
      Ind.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Ind.RJ.Prob)
      Dep.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Dep.RJ.Prob)
      
      Ind.RJ.LL  <- LogLoss(rTaxonValueB[j], Ind.RJ.Prob)
      Dep.RJ.LL  <- LogLoss(rTaxonValueB[j], Dep.RJ.Prob)    
      
      if (isTRUE(multistate_prediction)) {
        MS.RJ.acc <- calculate_accuracy(rTaxonValueB[j], MS.RJ.Prob)
        MS.RJ.LL  <- LogLoss(rTaxonValueB[j], MS.RJ.Prob)
      }
    }
    
    
    # Finally assign the variables to the final table based on the RJ model
    if (RJmodel %in% c("MCMC", "BOTH")) {
      MCMC_prob <- c(Ind.MCMC.Prob, Dep.MCMC.Prob)
      MCMC_acc  <- c(Ind.MCMC.acc, Dep.MCMC.acc)
      MCMC_LL   <- c(Ind.MCMC.LL,  Dep.MCMC.LL)
      
      if (isTRUE(multistate_prediction)) {
        MCMC_prob <- c(MS.MCMC.Prob, MCMC_prob)
        MCMC_acc  <- c(MS.MCMC.acc,  MCMC_acc)
        MCMC_LL   <- c(MS.MCMC.LL,   MCMC_LL)
      }
    }
    
    if (RJmodel %in% c("RJMCMC", "BOTH")) {
      RJ_prob <- c(Ind.RJ.Prob, Dep.RJ.Prob)
      RJ_acc  <- c(Ind.RJ.acc,  Dep.RJ.acc)
      RJ_LL   <- c(Ind.RJ.LL,   Dep.RJ.LL)
      
      if (isTRUE(multistate_prediction)) {
        RJ_prob <- c(MS.RJ.Prob, RJ_prob)
        RJ_acc  <- c(MS.RJ.acc,  RJ_acc)
        RJ_LL   <- c(MS.RJ.LL,   RJ_LL)
      }
    }
    
    if (RJmodel == "MCMC") {
      result_vector <- c(MCMC_prob, MCMC_acc, MCMC_LL)
    } else if (RJmodel == "RJMCMC") {
      result_vector <- c(RJ_prob, RJ_acc, RJ_LL)
    } else if (RJmodel == "BOTH") {
      result_vector <- c(MCMC_prob, RJ_prob, MCMC_acc, RJ_acc, MCMC_LL, RJ_LL)
    }
    
    # Assign results
    results[matching_rows[j], Results_labels] <- result_vector
    
  }
}
write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t") 



#####
#This is where the repeats for optional runs starts, beginning with variable rates
if (isTRUE(variable_rates)) {
  for (type in types) {
    print(paste0("Compiling the BayesTraits results from the single-prediction, ", type, " tests with variable rates."))
    # Call the necessary paths
    logfile_path <- paste0("VariableRates/", type, "/Single/", type)
    results_name <- paste0("Results/VariableRates/Single/", type, ".Single.ResultsFull.txt")
    
    # Pull the results into R
    colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
    results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
    
    # Start the for-loop that will run each prediction and fill it in
    for (i in 1:trial_amount) {
      # First, call rTaxon by subsetting the rows where column 1 matches i
      matching_rows <- which(results[, 1] == i)
      rTaxon <- results[matching_rows, 2]
      rTaxonValueA <- results[matching_rows, 3]
      rTaxonValueB <- results[matching_rows, 4]
      
      ### This loop will cover any situations where rTaxon is more than one value
      for (j in 1:length(rTaxon)) {
        
        # Gather the information from 
        skiplength <- length(rTaxon)
        Calculate_Post_Prob(i, RJmodel, multistate_prediction, skiplength, j, trial_name = logfile_path)
        
        # Use the RJmodel to pull the correct data
        if (RJmodel %in% c("MCMC", "BOTH")) {
          Ind.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Ind.MCMC.Prob)
          Dep.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Dep.MCMC.Prob)
          
          Ind.MCMC.LL  <- LogLoss(rTaxonValueB[j], Ind.MCMC.Prob)
          Dep.MCMC.LL  <- LogLoss(rTaxonValueB[j], Dep.MCMC.Prob)    
          
          if (isTRUE(multistate_prediction)) {
            MS.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], MS.MCMC.Prob)
            MS.MCMC.LL  <- LogLoss(rTaxonValueB[j], MS.MCMC.Prob)
            
          }
        }
        
        if (RJmodel %in% c("RJMCMC", "BOTH")) {
          Ind.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Ind.RJ.Prob)
          Dep.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Dep.RJ.Prob)
          
          Ind.RJ.LL  <- LogLoss(rTaxonValueB[j], Ind.RJ.Prob)
          Dep.RJ.LL  <- LogLoss(rTaxonValueB[j], Dep.RJ.Prob)    
          
          if (isTRUE(multistate_prediction)) {
            MS.RJ.acc <- calculate_accuracy(rTaxonValueB[j], MS.RJ.Prob)
            MS.RJ.LL  <- LogLoss(rTaxonValueB[j], MS.RJ.Prob)
          }
        }
        
        
        # Finally assign the variables to the final table based on the RJ model
        if (RJmodel %in% c("MCMC", "BOTH")) {
          MCMC_prob <- c(Ind.MCMC.Prob, Dep.MCMC.Prob)
          MCMC_acc  <- c(Ind.MCMC.acc, Dep.MCMC.acc)
          MCMC_LL   <- c(Ind.MCMC.LL,  Dep.MCMC.LL)
          
          if (isTRUE(multistate_prediction)) {
            MCMC_prob <- c(MS.MCMC.Prob, MCMC_prob)
            MCMC_acc  <- c(MS.MCMC.acc,  MCMC_acc)
            MCMC_LL   <- c(MS.MCMC.LL,   MCMC_LL)
          }
        }
        
        if (RJmodel %in% c("RJMCMC", "BOTH")) {
          RJ_prob <- c(Ind.RJ.Prob, Dep.RJ.Prob)
          RJ_acc  <- c(Ind.RJ.acc,  Dep.RJ.acc)
          RJ_LL   <- c(Ind.RJ.LL,   Dep.RJ.LL)
          
          if (isTRUE(multistate_prediction)) {
            RJ_prob <- c(MS.RJ.Prob, RJ_prob)
            RJ_acc  <- c(MS.RJ.acc,  RJ_acc)
            RJ_LL   <- c(MS.RJ.LL,   RJ_LL)
          }
        }
        
        if (RJmodel == "MCMC") {
          result_vector <- c(MCMC_prob, MCMC_acc, MCMC_LL)
        } else if (RJmodel == "RJMCMC") {
          result_vector <- c(RJ_prob, RJ_acc, RJ_LL)
        } else if (RJmodel == "BOTH") {
          result_vector <- c(MCMC_prob, RJ_prob, MCMC_acc, RJ_acc, MCMC_LL, RJ_LL)
        }
        
        # Assign results
        results[matching_rows[j], Results_labels] <- result_vector
        
      }
    }
    write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t") 
  }
}



### Next, a long run for all multiple predictions, beginning with constant rates
if (isTRUE(multiple_prediction)) {
  for (type in types) {
    print(paste0("Compiling the BayesTraits results from the multiple-prediction, ", type, " tests with constant rates."))
    # Call the necessary paths
    logfile_path <- paste0("ConstantRates/", type, "/Multiple/", type)
    results_name <- paste0("Results/ConstantRates/Multiple/", type, ".Multiple.ResultsFull.txt")
    
    # Pull the results into R
    colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
    results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
    
    # Start the for-loop that will run each prediction and fill it in
    for (i in 1:trial_amount) {
      # First, call rTaxon by subsetting the rows where column 1 matches i
      matching_rows <- which(results[, 1] == i)
      rTaxon <- results[matching_rows, 2]
      rTaxonValueA <- results[matching_rows, 3]
      rTaxonValueB <- results[matching_rows, 4]
      
      ### This loop will cover any situations where rTaxon is more than one value
      for (j in 1:length(rTaxon)) {
        
        # Gather the information from 
        skiplength <- length(rTaxon)
        Calculate_Post_Prob(i, RJmodel, multistate_prediction, skiplength, j, trial_name = logfile_path)
        
        # Use the RJmodel to pull the correct data
        if (RJmodel %in% c("MCMC", "BOTH")) {
          Ind.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Ind.MCMC.Prob)
          Dep.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Dep.MCMC.Prob)
          
          Ind.MCMC.LL  <- LogLoss(rTaxonValueB[j], Ind.MCMC.Prob)
          Dep.MCMC.LL  <- LogLoss(rTaxonValueB[j], Dep.MCMC.Prob)    
          
          if (isTRUE(multistate_prediction)) {
            MS.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], MS.MCMC.Prob)
            MS.MCMC.LL  <- LogLoss(rTaxonValueB[j], MS.MCMC.Prob)
            
          }
        }
        
        if (RJmodel %in% c("RJMCMC", "BOTH")) {
          Ind.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Ind.RJ.Prob)
          Dep.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Dep.RJ.Prob)
          
          Ind.RJ.LL  <- LogLoss(rTaxonValueB[j], Ind.RJ.Prob)
          Dep.RJ.LL  <- LogLoss(rTaxonValueB[j], Dep.RJ.Prob)    
          
          if (isTRUE(multistate_prediction)) {
            MS.RJ.acc <- calculate_accuracy(rTaxonValueB[j], MS.RJ.Prob)
            MS.RJ.LL  <- LogLoss(rTaxonValueB[j], MS.RJ.Prob)
          }
        }
        
        
        # Finally assign the variables to the final table based on the RJ model
        if (RJmodel %in% c("MCMC", "BOTH")) {
          MCMC_prob <- c(Ind.MCMC.Prob, Dep.MCMC.Prob)
          MCMC_acc  <- c(Ind.MCMC.acc, Dep.MCMC.acc)
          MCMC_LL   <- c(Ind.MCMC.LL,  Dep.MCMC.LL)
          
          if (isTRUE(multistate_prediction)) {
            MCMC_prob <- c(MS.MCMC.Prob, MCMC_prob)
            MCMC_acc  <- c(MS.MCMC.acc,  MCMC_acc)
            MCMC_LL   <- c(MS.MCMC.LL,   MCMC_LL)
          }
        }
        
        if (RJmodel %in% c("RJMCMC", "BOTH")) {
          RJ_prob <- c(Ind.RJ.Prob, Dep.RJ.Prob)
          RJ_acc  <- c(Ind.RJ.acc,  Dep.RJ.acc)
          RJ_LL   <- c(Ind.RJ.LL,   Dep.RJ.LL)
          
          if (isTRUE(multistate_prediction)) {
            RJ_prob <- c(MS.RJ.Prob, RJ_prob)
            RJ_acc  <- c(MS.RJ.acc,  RJ_acc)
            RJ_LL   <- c(MS.RJ.LL,   RJ_LL)
          }
        }
        
        if (RJmodel == "MCMC") {
          result_vector <- c(MCMC_prob, MCMC_acc, MCMC_LL)
        } else if (RJmodel == "RJMCMC") {
          result_vector <- c(RJ_prob, RJ_acc, RJ_LL)
        } else if (RJmodel == "BOTH") {
          result_vector <- c(MCMC_prob, RJ_prob, MCMC_acc, RJ_acc, MCMC_LL, RJ_LL)
        }
        
        # Assign results
        results[matching_rows[j], Results_labels] <- result_vector
        
      }
    }
    write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t") 
  }
  
  
  
  ### Now we will do the same for the randomly generated data
  print(paste0("Compiling the BayesTraits results from the multiple-prediction, tests with random data."))
  # Call the necessary paths
  logfile_path <- "Random/Multiple/Random"
  results_name <- "Results/Random/Random.Multiple.ResultsFull.txt"
  
  # Pull the results into R
  colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
  results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
  
  # Start the for-loop that will run each prediction and fill it in
  for (i in 1:trial_amount) {
    # First, call rTaxon by subsetting the rows where column 1 matches i
    matching_rows <- which(results[, 1] == i)
    rTaxon <- results[matching_rows, 2]
    rTaxonValueA <- results[matching_rows, 3]
    rTaxonValueB <- results[matching_rows, 4]
    
    ### This loop will cover any situations where rTaxon is more than one value
    for (j in 1:length(rTaxon)) {
      
      # Gather the information from 
      skiplength <- length(rTaxon)
      Calculate_Post_Prob(i, RJmodel, multistate_prediction, skiplength, j, trial_name = logfile_path)
      
      # Use the RJmodel to pull the correct data
      if (RJmodel %in% c("MCMC", "BOTH")) {
        Ind.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Ind.MCMC.Prob)
        Dep.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Dep.MCMC.Prob)
        
        Ind.MCMC.LL  <- LogLoss(rTaxonValueB[j], Ind.MCMC.Prob)
        Dep.MCMC.LL  <- LogLoss(rTaxonValueB[j], Dep.MCMC.Prob)    
        
        if (isTRUE(multistate_prediction)) {
          MS.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], MS.MCMC.Prob)
          MS.MCMC.LL  <- LogLoss(rTaxonValueB[j], MS.MCMC.Prob)
          
        }
      }
      
      if (RJmodel %in% c("RJMCMC", "BOTH")) {
        Ind.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Ind.RJ.Prob)
        Dep.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Dep.RJ.Prob)
        
        Ind.RJ.LL  <- LogLoss(rTaxonValueB[j], Ind.RJ.Prob)
        Dep.RJ.LL  <- LogLoss(rTaxonValueB[j], Dep.RJ.Prob)    
        
        if (isTRUE(multistate_prediction)) {
          MS.RJ.acc <- calculate_accuracy(rTaxonValueB[j], MS.RJ.Prob)
          MS.RJ.LL  <- LogLoss(rTaxonValueB[j], MS.RJ.Prob)
        }
      }
      
      
      # Finally assign the variables to the final table based on the RJ model
      if (RJmodel %in% c("MCMC", "BOTH")) {
        MCMC_prob <- c(Ind.MCMC.Prob, Dep.MCMC.Prob)
        MCMC_acc  <- c(Ind.MCMC.acc, Dep.MCMC.acc)
        MCMC_LL   <- c(Ind.MCMC.LL,  Dep.MCMC.LL)
        
        if (isTRUE(multistate_prediction)) {
          MCMC_prob <- c(MS.MCMC.Prob, MCMC_prob)
          MCMC_acc  <- c(MS.MCMC.acc,  MCMC_acc)
          MCMC_LL   <- c(MS.MCMC.LL,   MCMC_LL)
        }
      }
      
      if (RJmodel %in% c("RJMCMC", "BOTH")) {
        RJ_prob <- c(Ind.RJ.Prob, Dep.RJ.Prob)
        RJ_acc  <- c(Ind.RJ.acc,  Dep.RJ.acc)
        RJ_LL   <- c(Ind.RJ.LL,   Dep.RJ.LL)
        
        if (isTRUE(multistate_prediction)) {
          RJ_prob <- c(MS.RJ.Prob, RJ_prob)
          RJ_acc  <- c(MS.RJ.acc,  RJ_acc)
          RJ_LL   <- c(MS.RJ.LL,   RJ_LL)
        }
      }
      
      if (RJmodel == "MCMC") {
        result_vector <- c(MCMC_prob, MCMC_acc, MCMC_LL)
      } else if (RJmodel == "RJMCMC") {
        result_vector <- c(RJ_prob, RJ_acc, RJ_LL)
      } else if (RJmodel == "BOTH") {
        result_vector <- c(MCMC_prob, RJ_prob, MCMC_acc, RJ_acc, MCMC_LL, RJ_LL)
      }
      
      # Assign results
      results[matching_rows[j], Results_labels] <- result_vector
      
    }
  }
  write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t") 
  
  
  

  ### Now we will repeat the process for the variable rates with multiple prediction
  if (isTRUE(variable_rates)) {
    for (type in types) {
      print(paste0("Compiling the BayesTraits results from the multiple-prediction, ", type, " tests with variable rates."))
      # Call the necessary paths
      logfile_path <- paste0("VariableRates/", type, "/Multiple/", type)
      results_name <- paste0("Results/VariableRates/Multiple/", type, ".Multiple.ResultsFull.txt")
      
      # Pull the results into R
      colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
      results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
      
      # Start the for-loop that will run each prediction and fill it in
      for (i in 1:trial_amount) {
        # First, call rTaxon by subsetting the rows where column 1 matches i
        matching_rows <- which(results[, 1] == i)
        rTaxon <- results[matching_rows, 2]
        rTaxonValueA <- results[matching_rows, 3]
        rTaxonValueB <- results[matching_rows, 4]
        
        ### This loop will cover any situations where rTaxon is more than one value
        for (j in 1:length(rTaxon)) {
          
          # Gather the information from 
          skiplength <- length(rTaxon)
          Calculate_Post_Prob(i, RJmodel, multistate_prediction, skiplength, j, trial_name = logfile_path)
          
          # Use the RJmodel to pull the correct data
          if (RJmodel %in% c("MCMC", "BOTH")) {
            Ind.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Ind.MCMC.Prob)
            Dep.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Dep.MCMC.Prob)
            
            Ind.MCMC.LL  <- LogLoss(rTaxonValueB[j], Ind.MCMC.Prob)
            Dep.MCMC.LL  <- LogLoss(rTaxonValueB[j], Dep.MCMC.Prob)    
            
            if (isTRUE(multistate_prediction)) {
              MS.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], MS.MCMC.Prob)
              MS.MCMC.LL  <- LogLoss(rTaxonValueB[j], MS.MCMC.Prob)
              
            }
          }
          
          if (RJmodel %in% c("RJMCMC", "BOTH")) {
            Ind.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Ind.RJ.Prob)
            Dep.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Dep.RJ.Prob)
            
            Ind.RJ.LL  <- LogLoss(rTaxonValueB[j], Ind.RJ.Prob)
            Dep.RJ.LL  <- LogLoss(rTaxonValueB[j], Dep.RJ.Prob)    
            
            if (isTRUE(multistate_prediction)) {
              MS.RJ.acc <- calculate_accuracy(rTaxonValueB[j], MS.RJ.Prob)
              MS.RJ.LL  <- LogLoss(rTaxonValueB[j], MS.RJ.Prob)
            }
          }
          
          
          # Finally assign the variables to the final table based on the RJ model
          if (RJmodel %in% c("MCMC", "BOTH")) {
            MCMC_prob <- c(Ind.MCMC.Prob, Dep.MCMC.Prob)
            MCMC_acc  <- c(Ind.MCMC.acc, Dep.MCMC.acc)
            MCMC_LL   <- c(Ind.MCMC.LL,  Dep.MCMC.LL)
            
            if (isTRUE(multistate_prediction)) {
              MCMC_prob <- c(MS.MCMC.Prob, MCMC_prob)
              MCMC_acc  <- c(MS.MCMC.acc,  MCMC_acc)
              MCMC_LL   <- c(MS.MCMC.LL,   MCMC_LL)
            }
          }
          
          if (RJmodel %in% c("RJMCMC", "BOTH")) {
            RJ_prob <- c(Ind.RJ.Prob, Dep.RJ.Prob)
            RJ_acc  <- c(Ind.RJ.acc,  Dep.RJ.acc)
            RJ_LL   <- c(Ind.RJ.LL,   Dep.RJ.LL)
            
            if (isTRUE(multistate_prediction)) {
              RJ_prob <- c(MS.RJ.Prob, RJ_prob)
              RJ_acc  <- c(MS.RJ.acc,  RJ_acc)
              RJ_LL   <- c(MS.RJ.LL,   RJ_LL)
            }
          }
          
          if (RJmodel == "MCMC") {
            result_vector <- c(MCMC_prob, MCMC_acc, MCMC_LL)
          } else if (RJmodel == "RJMCMC") {
            result_vector <- c(RJ_prob, RJ_acc, RJ_LL)
          } else if (RJmodel == "BOTH") {
            result_vector <- c(MCMC_prob, RJ_prob, MCMC_acc, RJ_acc, MCMC_LL, RJ_LL)
          }
          
          # Assign results
          results[matching_rows[j], Results_labels] <- result_vector
          
        }
      }
      write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t") 
    }
  }
}



### Finally, a long run for all clade predictions, beginning with constant rates
if (isTRUE(clade_prediction)) {
  for (type in types) {
    print(paste0("Compiling the BayesTraits results from the clade-prediction, ", type, " tests with equal rates."))
    # Call the necessary paths
    logfile_path <- paste0("ConstantRates/", type, "/Clade/", type)
    results_name <- paste0("Results/ConstantRates/Clade/", type, ".Clade.ResultsFull.txt")
    
    # Pull the results into R
    colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
    results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
    
    # Start the for-loop that will run each prediction and fill it in
    for (i in 1:trial_amount) {
      # First, call rTaxon by subsetting the rows where column 1 matches i
      matching_rows <- which(results[, 1] == i)
      rTaxon <- results[matching_rows, 2]
      rTaxonValueA <- results[matching_rows, 3]
      rTaxonValueB <- results[matching_rows, 4]
      
      ### This loop will cover any situations where rTaxon is more than one value
      for (j in 1:length(rTaxon)) {
        
        # Gather the information from 
        skiplength <- length(rTaxon)
        Calculate_Post_Prob(i, RJmodel, multistate_prediction, skiplength, j, trial_name = logfile_path)
        
        # Use the RJmodel to pull the correct data
        if (RJmodel %in% c("MCMC", "BOTH")) {
          Ind.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Ind.MCMC.Prob)
          Dep.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Dep.MCMC.Prob)
          
          Ind.MCMC.LL  <- LogLoss(rTaxonValueB[j], Ind.MCMC.Prob)
          Dep.MCMC.LL  <- LogLoss(rTaxonValueB[j], Dep.MCMC.Prob)    
          
          if (isTRUE(multistate_prediction)) {
            MS.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], MS.MCMC.Prob)
            MS.MCMC.LL  <- LogLoss(rTaxonValueB[j], MS.MCMC.Prob)
            
          }
        }
        
        if (RJmodel %in% c("RJMCMC", "BOTH")) {
          Ind.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Ind.RJ.Prob)
          Dep.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Dep.RJ.Prob)
          
          Ind.RJ.LL  <- LogLoss(rTaxonValueB[j], Ind.RJ.Prob)
          Dep.RJ.LL  <- LogLoss(rTaxonValueB[j], Dep.RJ.Prob)    
          
          if (isTRUE(multistate_prediction)) {
            MS.RJ.acc <- calculate_accuracy(rTaxonValueB[j], MS.RJ.Prob)
            MS.RJ.LL  <- LogLoss(rTaxonValueB[j], MS.RJ.Prob)
          }
        }
        
        
        # Finally assign the variables to the final table based on the RJ model
        if (RJmodel %in% c("MCMC", "BOTH")) {
          MCMC_prob <- c(Ind.MCMC.Prob, Dep.MCMC.Prob)
          MCMC_acc  <- c(Ind.MCMC.acc, Dep.MCMC.acc)
          MCMC_LL   <- c(Ind.MCMC.LL,  Dep.MCMC.LL)
          
          if (isTRUE(multistate_prediction)) {
            MCMC_prob <- c(MS.MCMC.Prob, MCMC_prob)
            MCMC_acc  <- c(MS.MCMC.acc,  MCMC_acc)
            MCMC_LL   <- c(MS.MCMC.LL,   MCMC_LL)
          }
        }
        
        if (RJmodel %in% c("RJMCMC", "BOTH")) {
          RJ_prob <- c(Ind.RJ.Prob, Dep.RJ.Prob)
          RJ_acc  <- c(Ind.RJ.acc,  Dep.RJ.acc)
          RJ_LL   <- c(Ind.RJ.LL,   Dep.RJ.LL)
          
          if (isTRUE(multistate_prediction)) {
            RJ_prob <- c(MS.RJ.Prob, RJ_prob)
            RJ_acc  <- c(MS.RJ.acc,  RJ_acc)
            RJ_LL   <- c(MS.RJ.LL,   RJ_LL)
          }
        }
        
        if (RJmodel == "MCMC") {
          result_vector <- c(MCMC_prob, MCMC_acc, MCMC_LL)
        } else if (RJmodel == "RJMCMC") {
          result_vector <- c(RJ_prob, RJ_acc, RJ_LL)
        } else if (RJmodel == "BOTH") {
          result_vector <- c(MCMC_prob, RJ_prob, MCMC_acc, RJ_acc, MCMC_LL, RJ_LL)
        }
        
        # Assign results
        results[matching_rows[j], Results_labels] <- result_vector
        
      }
    }
    write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t") 
  }
  
  
  
  ### Now we will do the same for the randomly generated data
  print(paste0("Compiling the BayesTraits results from the clade-prediction, tests with random data."))
  # Call the necessary paths
  logfile_path <- "Random/Clade/Random"
  results_name <- "Results/Random/Random.Clade.ResultsFull.txt"
  
  # Pull the results into R
  colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
  results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
  
  # Start the for-loop that will run each prediction and fill it in
  for (i in 1:trial_amount) {
    # First, call rTaxon by subsetting the rows where column 1 matches i
    matching_rows <- which(results[, 1] == i)
    rTaxon <- results[matching_rows, 2]
    rTaxonValueA <- results[matching_rows, 3]
    rTaxonValueB <- results[matching_rows, 4]
    
    ### This loop will cover any situations where rTaxon is more than one value
    for (j in 1:length(rTaxon)) {
      
      # Gather the information from 
      skiplength <- length(rTaxon)
      Calculate_Post_Prob(i, RJmodel, multistate_prediction, skiplength, j, trial_name = logfile_path)
      
      # Use the RJmodel to pull the correct data
      if (RJmodel %in% c("MCMC", "BOTH")) {
        Ind.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Ind.MCMC.Prob)
        Dep.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Dep.MCMC.Prob)
        
        Ind.MCMC.LL  <- LogLoss(rTaxonValueB[j], Ind.MCMC.Prob)
        Dep.MCMC.LL  <- LogLoss(rTaxonValueB[j], Dep.MCMC.Prob)    
        
        if (isTRUE(multistate_prediction)) {
          MS.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], MS.MCMC.Prob)
          MS.MCMC.LL  <- LogLoss(rTaxonValueB[j], MS.MCMC.Prob)
          
        }
      }
      
      if (RJmodel %in% c("RJMCMC", "BOTH")) {
        Ind.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Ind.RJ.Prob)
        Dep.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Dep.RJ.Prob)
        
        Ind.RJ.LL  <- LogLoss(rTaxonValueB[j], Ind.RJ.Prob)
        Dep.RJ.LL  <- LogLoss(rTaxonValueB[j], Dep.RJ.Prob)    
        
        if (isTRUE(multistate_prediction)) {
          MS.RJ.acc <- calculate_accuracy(rTaxonValueB[j], MS.RJ.Prob)
          MS.RJ.LL  <- LogLoss(rTaxonValueB[j], MS.RJ.Prob)
        }
      }
      
      
      # Finally assign the variables to the final table based on the RJ model
      if (RJmodel %in% c("MCMC", "BOTH")) {
        MCMC_prob <- c(Ind.MCMC.Prob, Dep.MCMC.Prob)
        MCMC_acc  <- c(Ind.MCMC.acc, Dep.MCMC.acc)
        MCMC_LL   <- c(Ind.MCMC.LL,  Dep.MCMC.LL)
        
        if (isTRUE(multistate_prediction)) {
          MCMC_prob <- c(MS.MCMC.Prob, MCMC_prob)
          MCMC_acc  <- c(MS.MCMC.acc,  MCMC_acc)
          MCMC_LL   <- c(MS.MCMC.LL,   MCMC_LL)
        }
      }
      
      if (RJmodel %in% c("RJMCMC", "BOTH")) {
        RJ_prob <- c(Ind.RJ.Prob, Dep.RJ.Prob)
        RJ_acc  <- c(Ind.RJ.acc,  Dep.RJ.acc)
        RJ_LL   <- c(Ind.RJ.LL,   Dep.RJ.LL)
        
        if (isTRUE(multistate_prediction)) {
          RJ_prob <- c(MS.RJ.Prob, RJ_prob)
          RJ_acc  <- c(MS.RJ.acc,  RJ_acc)
          RJ_LL   <- c(MS.RJ.LL,   RJ_LL)
        }
      }
      
      if (RJmodel == "MCMC") {
        result_vector <- c(MCMC_prob, MCMC_acc, MCMC_LL)
      } else if (RJmodel == "RJMCMC") {
        result_vector <- c(RJ_prob, RJ_acc, RJ_LL)
      } else if (RJmodel == "BOTH") {
        result_vector <- c(MCMC_prob, RJ_prob, MCMC_acc, RJ_acc, MCMC_LL, RJ_LL)
      }
      
      # Assign results
      results[matching_rows[j], Results_labels] <- result_vector
      
    }
  }
  write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t") 
  
  
  
  ### Finally, we run this for the last time with variable rates and clade-prediction
  if (isTRUE(variable_rates)) {
    for (type in types) {
      print(paste0("Compiling the BayesTraits results from the clade-prediction, ", type, " tests with equal rates."))
      # Call the necessary paths
      logfile_path <- paste0("VariableRates/", type, "/Clade/", type)
      results_name <- paste0("Results/VariableRates/Clade/", type, ".Clade.ResultsFull.txt")
      
      # Pull the results into R
      colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
      results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
      
      # Start the for-loop that will run each prediction and fill it in
      for (i in 1:trial_amount) {
        # First, call rTaxon by subsetting the rows where column 1 matches i
        matching_rows <- which(results[, 1] == i)
        rTaxon <- results[matching_rows, 2]
        rTaxonValueA <- results[matching_rows, 3]
        rTaxonValueB <- results[matching_rows, 4]
        
        ### This loop will cover any situations where rTaxon is more than one value
        for (j in 1:length(rTaxon)) {
          
          # Gather the information from 
          skiplength <- length(rTaxon)
          Calculate_Post_Prob(i, RJmodel, multistate_prediction, skiplength, j, trial_name = logfile_path)
          
          # Use the RJmodel to pull the correct data
          if (RJmodel %in% c("MCMC", "BOTH")) {
            Ind.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Ind.MCMC.Prob)
            Dep.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], Dep.MCMC.Prob)
            
            Ind.MCMC.LL  <- LogLoss(rTaxonValueB[j], Ind.MCMC.Prob)
            Dep.MCMC.LL  <- LogLoss(rTaxonValueB[j], Dep.MCMC.Prob)    
            
            if (isTRUE(multistate_prediction)) {
              MS.MCMC.acc <- calculate_accuracy(rTaxonValueB[j], MS.MCMC.Prob)
              MS.MCMC.LL  <- LogLoss(rTaxonValueB[j], MS.MCMC.Prob)
              
            }
          }
          
          if (RJmodel %in% c("RJMCMC", "BOTH")) {
            Ind.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Ind.RJ.Prob)
            Dep.RJ.acc <- calculate_accuracy(rTaxonValueB[j], Dep.RJ.Prob)
            
            Ind.RJ.LL  <- LogLoss(rTaxonValueB[j], Ind.RJ.Prob)
            Dep.RJ.LL  <- LogLoss(rTaxonValueB[j], Dep.RJ.Prob)    
            
            if (isTRUE(multistate_prediction)) {
              MS.RJ.acc <- calculate_accuracy(rTaxonValueB[j], MS.RJ.Prob)
              MS.RJ.LL  <- LogLoss(rTaxonValueB[j], MS.RJ.Prob)
            }
          }
          
          
          # Finally assign the variables to the final table based on the RJ model
          if (RJmodel %in% c("MCMC", "BOTH")) {
            MCMC_prob <- c(Ind.MCMC.Prob, Dep.MCMC.Prob)
            MCMC_acc  <- c(Ind.MCMC.acc, Dep.MCMC.acc)
            MCMC_LL   <- c(Ind.MCMC.LL,  Dep.MCMC.LL)
            
            if (isTRUE(multistate_prediction)) {
              MCMC_prob <- c(MS.MCMC.Prob, MCMC_prob)
              MCMC_acc  <- c(MS.MCMC.acc,  MCMC_acc)
              MCMC_LL   <- c(MS.MCMC.LL,   MCMC_LL)
            }
          }
          
          if (RJmodel %in% c("RJMCMC", "BOTH")) {
            RJ_prob <- c(Ind.RJ.Prob, Dep.RJ.Prob)
            RJ_acc  <- c(Ind.RJ.acc,  Dep.RJ.acc)
            RJ_LL   <- c(Ind.RJ.LL,   Dep.RJ.LL)
            
            if (isTRUE(multistate_prediction)) {
              RJ_prob <- c(MS.RJ.Prob, RJ_prob)
              RJ_acc  <- c(MS.RJ.acc,  RJ_acc)
              RJ_LL   <- c(MS.RJ.LL,   RJ_LL)
            }
          }
          
          if (RJmodel == "MCMC") {
            result_vector <- c(MCMC_prob, MCMC_acc, MCMC_LL)
          } else if (RJmodel == "RJMCMC") {
            result_vector <- c(RJ_prob, RJ_acc, RJ_LL)
          } else if (RJmodel == "BOTH") {
            result_vector <- c(MCMC_prob, RJ_prob, MCMC_acc, RJ_acc, MCMC_LL, RJ_LL)
          }
          
          # Assign results
          results[matching_rows[j], Results_labels] <- result_vector
          
        }
      }
      write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t") 
    }
  }
}
print("All BayesTraits predictions have been compiled into the Results folder.")