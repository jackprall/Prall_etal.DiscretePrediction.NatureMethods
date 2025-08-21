################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# June 2025

################################################################################
# This will add the average sister taxon's/taxa's trait values to the final results tables.

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

# Number of data sets/trials for each combination of settings
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))



#####
# Install and call the necessary 
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
# We start the loop that will run the trials
for (type in types) {
  # Call the necessary paths
  data_path <- paste0("ConstantRates/", type, "/Data/", type)
  results_name <- paste0("Results/ConstantRates/Single/", type, ".Single.ResultsFull.txt")
  
  # Pull the results into R
  colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
  results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
  
  # Start the for-loop that will run each prediction and fill it in
  for (i in 1:trial_amount) {
    # First, call rTaxon by subsetting the rows where column 1 matches i
    matching_rows <- which(results[, 1] == i)
    rTaxon <- results[matching_rows, 2]
    
    # Call the data file for this iteration
    data_name <- paste0(data_path, ".", i, ".Full_data.txt")
    data <- read.table(data_name, skip = 1, sep = "\t")
    
    # Pull the tree to read the terminal branch length
    treename <- paste0("Trees/Full_tree.", i, ".tre")
    tree_complete <- read.nexus(treename)  
 
    for (j in 1:length(rTaxon)) {
      sisters_A <- sisterData(tree_complete, rTaxon[j], data, "A")
      sisters_B <- sisterData(tree_complete, rTaxon[j], data, "B")
      
      # Accuracies and Log loss scores are calculated as they're added  
      results[matching_rows[j], c("Avg_Sister_A", "Avg_Sister_B")] <-
        c(sisters_A, sisters_B)   
  
    }
  }
  # Save the table with the new information from the sister(s)
  write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
  


  ##### 
  # Check to see if we are doing this for variable rates
  if (variable_rates == TRUE) {
    for (type in types) {
      # Call the necessary paths
      data_path <- paste0("VariableRates/", type, "/Data/", type)
      results_name <- paste0("Results/VariableRates/Single/", type, ".Single.ResultsFull.txt")
      
      # Pull the results into R
      colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
      results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
      
      # Start the for-loop that will run each prediction and fill it in
      for (i in 1:trial_amount) {
        # First, call rTaxon by subsetting the rows where column 1 matches i
        matching_rows <- which(results[, 1] == i)
        rTaxon <- results[matching_rows, 2]
        
        # Call the data file for this iteration
        data_name <- paste0(data_path, ".", i, ".Full_data.txt")
        data <- read.table(data_name, skip = 1, sep = "\t")
        
        # Pull the tree to read the terminal branch length
        treename <- paste0("Trees/Full_tree.", i, ".tre")
        tree_complete <- read.nexus(treename)  
        
        for (j in 1:length(rTaxon)) {
          sisters_A <- sisterData(tree_complete, rTaxon[j], data, "A")
          sisters_B <- sisterData(tree_complete, rTaxon[j], data, "B")
          
          # Accuracies and Log loss scores are calculated as they're added  
          results[matching_rows[j], c("Avg_Sister_A", "Avg_Sister_B")] <-
            c(sisters_A, sisters_B)   
          
        }
      }
      # Save the table with the new information from the sister(s)
      write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
    }
  }
}



#####
# Now we run this for random data
# Call the necessary paths
data_path <- "Random/Data/Random"
results_name <- "Results/Random/Random.Single.ResultsFull.txt"
  
# Pull the results into R
colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)

# Start the for-loop that will run each prediction and fill it in
for (i in 1:trial_amount) {
  # First, call rTaxon by subsetting the rows where column 1 matches i
  matching_rows <- which(results[, 1] == i)
  rTaxon <- results[matching_rows, 2]
  
  # Call the data file for this iteration
  data_name <- paste0(data_path, ".", i, ".Full_data.txt")
  data <- read.table(data_name, skip = 1, sep = "\t")
  
  # Pull the tree to read the terminal branch length
  treename <- paste0("Trees/Full_tree.", i, ".tre")
  tree_complete <- read.nexus(treename)  
  
  for (j in 1:length(rTaxon)) {
    sisters_A <- sisterData(tree_complete, rTaxon[j], data, "A")
    sisters_B <- sisterData(tree_complete, rTaxon[j], data, "B")
    
    # Accuracies and Log loss scores are calculated as they're added  
    results[matching_rows[j], c("Avg_Sister_A", "Avg_Sister_B")] <-
      c(sisters_A, sisters_B)   
    
  }
}
  # Save the table with the new information from the sister(s)
write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
print("Finished gathering single-prediction sister taxa information.")



#####
# Now repeat the proccess for multi-tip prediction if necessary
if (isTRUE(multiple_prediction)) {
  # We start the loop that will run the trials
  for (type in types) {
    # Call the necessary paths
    data_path <- paste0("ConstantRates/", type, "/Data/", type)
    results_name <- paste0("Results/ConstantRates/Multiple/", type, ".Multiple.ResultsFull.txt")
    
    # Pull the results into R
    colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
    results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
    
    # Start the for-loop that will run each prediction and fill it in
    for (i in 1:trial_amount) {
      # First, call rTaxon by subsetting the rows where column 1 matches i
      matching_rows <- which(results[, 1] == i)
      rTaxon <- results[matching_rows, 2]
      
      # Call the data file for this iteration
      data_name <- paste0(data_path, ".", i, ".Full_data.txt")
      data <- read.table(data_name, skip = 1, sep = "\t")
      
      # Pull the tree to read the terminal branch length
      treename <- paste0("Trees/Full_tree.", i, ".tre")
      tree_complete <- read.nexus(treename)  
      
      for (j in 1:length(rTaxon)) {
        # Make sure you aren't using the other unknown taxa to evaluate the sister taxa
        other_unknowns <- setdiff(rTaxon, rTaxon[j])
        tree_edited <- drop.tip(tree_complete, other_unknowns)
        data_edited <- data[-as.numeric(other_unknowns), ]
        
        sisters_A <- sisterData(tree_edited, rTaxon[j], data_edited, "A")
        sisters_B <- sisterData(tree_edited, rTaxon[j], data_edited, "B")
        
        # Accuracies and Log loss scores are calculated as they're added  
        results[matching_rows[j], c("Avg_Sister_A", "Avg_Sister_B")] <-
          c(sisters_A, sisters_B)   
        
      }
    } 
    # Save the table with the new information from the sister(s)
    write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
    
    
    
    ##### 
    # Check to see if we are doing this for variable rates
    if (variable_rates == TRUE) {
      
      # Call the necessary paths
      data_path <- paste0("VariableRates/", type, "/Data/", type)
      results_name <- paste0("Results/VariableRates/Multiple/", type, ".Multiple.ResultsFull.txt")
      
      # Pull the results into R
      colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
      results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
      
      # Start the for-loop that will run each prediction and fill it in
      for (i in 1:trial_amount) {
        # First, call rTaxon by subsetting the rows where column 1 matches i
        matching_rows <- which(results[, 1] == i)
        rTaxon <- results[matching_rows, 2]
        
        # Call the data file for this iteration
        data_name <- paste0(data_path, ".", i, ".Full_data.txt")
        data <- read.table(data_name, skip = 1, sep = "\t")
        
        # Pull the tree to read the terminal branch length
        treename <- paste0("Trees/Full_tree.", i, ".tre")
        tree_complete <- read.nexus(treename)  
        
        for (j in 1:length(rTaxon)) {
          # Make sure you aren't using the other unknown taxa to evaluate the sister taxa
          other_unknowns <- setdiff(rTaxon, rTaxon[j])
          tree_edited <- drop.tip(tree_complete, other_unknowns)
          data_edited <- data[-as.numeric(other_unknowns), ]
          
          sisters_A <- sisterData(tree_edited, rTaxon[j], data_edited, "A")
          sisters_B <- sisterData(tree_edited, rTaxon[j], data_edited, "B")
          
          # Accuracies and Log loss scores are calculated as they're added  
          results[matching_rows[j], c("Avg_Sister_A", "Avg_Sister_B")] <-
            c(sisters_A, sisters_B)   
          
        }
      }
      # Save the table with the new information from the sister(s)
      write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
    }
  }

  
  
  #####
  # Now we run this for random data
  # Call the necessary paths
  data_path <- "Random/Data/Random"
  results_name <- "Results/Random/Random.Multiple.ResultsFull.txt"
  
  # Pull the results into R
  colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
  results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
  
  # Start the for-loop that will run each prediction and fill it in
  for (i in 1:trial_amount) {
    # First, call rTaxon by subsetting the rows where column 1 matches i
    matching_rows <- which(results[, 1] == i)
    rTaxon <- results[matching_rows, 2]
    
    # Call the data file for this iteration
    data_name <- paste0(data_path, ".", i, ".Full_data.txt")
    data <- read.table(data_name, skip = 1, sep = "\t")
    
    # Pull the tree to read the terminal branch length
    treename <- paste0("Trees/Full_tree.", i, ".tre")
    tree_complete <- read.nexus(treename)  
    
    for (j in 1:length(rTaxon)) {
      # Make sure you aren't using the other unknown taxa to evaluate the sister taxa
      other_unknowns <- setdiff(rTaxon, rTaxon[j])
      tree_edited <- drop.tip(tree_complete, other_unknowns)
      data_edited <- data[-as.numeric(other_unknowns), ]
      
      sisters_A <- sisterData(tree_edited, rTaxon[j], data_edited, "A")
      sisters_B <- sisterData(tree_edited, rTaxon[j], data_edited, "B")
      
      # Accuracies and Log loss scores are calculated as they're added  
      results[matching_rows[j], c("Avg_Sister_A", "Avg_Sister_B")] <-
        c(sisters_A, sisters_B)   
      
    }
  }
  # Save the table with the new information from the sister(s)
  write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
  print("Finished gathering multiple-prediction sister taxa information.")
}



#####
# Now repeat this process for clade prediction if necessary
if (isTRUE(clade_prediction)) {
  # We start the loop that will run the trials
  for (type in types) {
    # Call the necessary paths
    data_path <- paste0("ConstantRates/", type, "/Data/", type)
    results_name <- paste0("Results/ConstantRates/Clade/", type, ".Clade.ResultsFull.txt")
    
    # Pull the results into R
    colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
    results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
    
    # Start the for-loop that will run each prediction and fill it in
    for (i in 1:trial_amount) {
      # First, call rTaxon by subsetting the rows where column 1 matches i
      matching_rows <- which(results[, 1] == i)
      rTaxon <- results[matching_rows, 2]
      
      # Call the data file for this iteration
      data_name <- paste0(data_path, ".", i, ".Full_data.txt")
      data <- read.table(data_name, skip = 1, sep = "\t")
      
      # Pull the tree to read the terminal branch length
      treename <- paste0("Trees/Full_tree.", i, ".tre")
      tree_complete <- read.nexus(treename)  
      
      for (j in 1:length(rTaxon)) {
        # Make sure you aren't using the other unknown taxa to evaluate the sister taxa
        other_unknowns <- setdiff(rTaxon, rTaxon[j])
        tree_edited <- drop.tip(tree_complete, other_unknowns)
        data_edited <- data[-as.numeric(other_unknowns), ]
        
        sisters_A <- sisterData(tree_edited, rTaxon[j], data_edited, "A")
        sisters_B <- sisterData(tree_edited, rTaxon[j], data_edited, "B")
        
        # Accuracies and Log loss scores are calculated as they're added  
        results[matching_rows[j], c("Avg_Sister_A", "Avg_Sister_B")] <-
          c(sisters_A, sisters_B)   
        
      }
    }
    # Save the table with the new information from the sister(s)
    write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
    


    #####
    # Check to see if we are doing this for variable rates
    if (variable_rates == TRUE) {
      
      # Call the necessary paths
      data_path <- paste0("VariableRates/", type, "/Data/", type)
      results_name <- paste0("Results/VariableRates/Clade/", type, ".Clade.ResultsFull.txt")
      
      # Pull the results into R
      colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
      results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
      
      # Start the for-loop that will run each prediction and fill it in
      for (i in 1:trial_amount) {
        # First, call rTaxon by subsetting the rows where column 1 matches i
        matching_rows <- which(results[, 1] == i)
        rTaxon <- results[matching_rows, 2]
        
        # Call the data file for this iteration
        data_name <- paste0(data_path, ".", i, ".Full_data.txt")
        data <- read.table(data_name, skip = 1, sep = "\t")
        
        # Pull the tree to read the terminal branch length
        treename <- paste0("Trees/Full_tree.", i, ".tre")
        tree_complete <- read.nexus(treename)  
        
        for (j in 1:length(rTaxon)) {
          # Make sure you aren't using the other unknown taxa to evaluate the sister taxa
          other_unknowns <- setdiff(rTaxon, rTaxon[j])
          tree_edited <- drop.tip(tree_complete, other_unknowns)
          data_edited <- data[-as.numeric(other_unknowns), ]
          
          sisters_A <- sisterData(tree_edited, rTaxon[j], data_edited, "A")
          sisters_B <- sisterData(tree_edited, rTaxon[j], data_edited, "B")
          
          # Accuracies and Log loss scores are calculated as they're added  
          results[matching_rows[j], c("Avg_Sister_A", "Avg_Sister_B")] <-
            c(sisters_A, sisters_B)   
          
        }
      }
      # Save the table with the new information from the sister(s)
      write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
    }
  }
  

  
  #####
  # Now we run this for random data
  # Call the necessary paths
  data_path <- "Random/Data/Random"
  results_name <- "Results/Random/Random.Clade.ResultsFull.txt"
  
  # Pull the results into R
  colnames_results <- strsplit(readLines(results_name, n = 1), "\t")[[1]]
  results <- read.table(results_name, skip = 1, sep = "\t", col.names = colnames_results)
  
  # Start the for-loop that will run each prediction and fill it in
  for (i in 1:trial_amount) {
    # First, call rTaxon by subsetting the rows where column 1 matches i
    matching_rows <- which(results[, 1] == i)
    rTaxon <- results[matching_rows, 2]
    
    # Call the data file for this iteration
    data_name <- paste0(data_path, ".", i, ".Full_data.txt")
    data <- read.table(data_name, skip = 1, sep = "\t")
    
    # Pull the tree to read the terminal branch length
    treename <- paste0("Trees/Full_tree.", i, ".tre")
    tree_complete <- read.nexus(treename)  
    
    for (j in 1:length(rTaxon)) {
      # Make sure you aren't using the other unknown taxa to evaluate the sister taxa
      other_unknowns <- setdiff(rTaxon, rTaxon[j])
      tree_edited <- drop.tip(tree_complete, other_unknowns)
      data_edited <- data[-as.numeric(other_unknowns), ]
      
      sisters_A <- sisterData(tree_edited, rTaxon[j], data_edited, "A")
      sisters_B <- sisterData(tree_edited, rTaxon[j], data_edited, "B")
      
      # Accuracies and Log loss scores are calculated as they're added  
      results[matching_rows[j], c("Avg_Sister_A", "Avg_Sister_B")] <-
        c(sisters_A, sisters_B)   
      
    }
  }
  # Save the table with the new information from the sister(s)
  write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
  print("Finished gathering clade-prediction sister taxa information.")
}
  
print("Stored the average trait values for the sister taxa of the unknown in the results files.")
