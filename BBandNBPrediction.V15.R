# This script is for making prediction with the Beta Binomial and Naive Bayes
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

# Number of iterations when running as batch job/ on HPC
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))



##### 
# Make the constant rates predictions
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
    rTaxonValueA <- results[matching_rows, 3]
    rTaxonValueB <- results[matching_rows, 4]
    
    # Call the data file for this iteration
    data_name <- paste0(data_path, ".", i, ".Full_data.txt")
    data <- read.table(data_name, skip = 1, sep = "\t")
    
    # Now establish some details we need to predict
    prediction_data <- data[-rTaxon,]
    
    OO <- sum(prediction_data[, 4] == 0)  # This saves as 0 instead of 00
    Ol <- sum(prediction_data[, 4] == 1)  # This saves as 1 instead of 01
    lO <- sum(prediction_data[, 4] == 10) # This saves intuitively
    ll <- sum(prediction_data[, 4] == 11) # This saves intuitively
    counts <- matrix(data = c(OO, Ol, lO, ll), nrow = 2, ncol = 2)
    
    prediction_data <- prediction_data[, c(2:3)]
    
    # Now we can make the Beta Binomial and Naive Bayes Predictions
    BB.Prob <- Beta_Bin_predict(prediction_data, length(rTaxon))
    NB.Prob <- Naive_Bayes_predict(counts, rTaxonValueA)
    
    # Assign these values to the results table
    # Accuracies and Log loss scores are calculated as they're added  
    results[i, c("Beta_Binom_Prob", "Beta_Binom_acc", "Beta_Binom_LL",
                 "Naive_Bayes_Prob", "Naive_Bayes_acc", "Naive_Bayes_LL")] <-
      c(BB.Prob, calculate_accuracy(rTaxonValueB, BB.Prob), LogLoss(rTaxonValueB, BB.Prob), NB.Prob,
        calculate_accuracy(rTaxonValueB, NB.Prob), LogLoss(rTaxonValueB, NB.Prob))
  }
  # Rewrite the results file with the new tests
  write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
  
  
  
  ### Check and potentially rerun this for Multiple and Clade prediction
  if (multiple_prediction == TRUE) {
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
        rTaxonValueA <- results[matching_rows, 3]
        rTaxonValueB <- results[matching_rows, 4]
        
        # Call the data file for this iteration
        data_name <- paste0(data_path, ".", i, ".Full_data.txt")
        data <- read.table(data_name, skip = 1, sep = "\t")
        
        # Now establish some details we need to predict
        prediction_data <- data[-rTaxon,]
        
        OO <- sum(prediction_data[, 4] == 0)  # This saves as 0 instead of 00
        Ol <- sum(prediction_data[, 4] == 1)  # This saves as 1 instead of 01
        lO <- sum(prediction_data[, 4] == 10) # This saves intuitively
        ll <- sum(prediction_data[, 4] == 11) # This saves intuitively
        counts <- matrix(data = c(OO, Ol, lO, ll), nrow = 2, ncol = 2)
        
        prediction_data <- prediction_data[, c(2:3)]
        
        # Now we can make the Beta Binomial and Naive Bayes Predictions
        BB.Prob <- Beta_Bin_predict(prediction_data, length(rTaxon))
        NB.Prob <- Naive_Bayes_predict(counts, rTaxonValueA)
        
        # Assign these values to the results table, with a small adjustment for multiple prediction
        # Naive Bayes accuracies and log loss scores are calculated as they're being assigned
        for (j in 1:length(rTaxon)) {
          results[matching_rows[j], c("Beta_Binom_Prob", "Beta_Binom_acc", "Beta_Binom_LL",
                                      "Naive_Bayes_Prob", "Naive_Bayes_acc", "Naive_Bayes_LL")] <-
            c(BB.Prob, calculate_accuracy(rTaxonValueB[j], BB.Prob), LogLoss(rTaxonValueB[j], BB.Prob), 
              NB.Prob[j], calculate_accuracy(rTaxonValueB[j], NB.Prob[j]), LogLoss(rTaxonValueB[j], NB.Prob[j])) 
        }
      }
      # Rewrite the results file with the new tests
      write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
    }
  }
  
  
  # Check and run clade prediction
  if (clade_prediction == TRUE) {
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
        rTaxonValueA <- results[matching_rows, 3]
        rTaxonValueB <- results[matching_rows, 4]
        
        # Call the data file for this iteration
        data_name <- paste0(data_path, ".", i, ".Full_data.txt")
        data <- read.table(data_name, skip = 1, sep = "\t")
        
        # Now establish some details we need to predict
        prediction_data <- data[-rTaxon,]
        
        OO <- sum(prediction_data[, 4] == 0)  # This saves as 0 instead of 00
        Ol <- sum(prediction_data[, 4] == 1)  # This saves as 1 instead of 01
        lO <- sum(prediction_data[, 4] == 10) # This saves intuitively
        ll <- sum(prediction_data[, 4] == 11) # This saves intuitively
        counts <- matrix(data = c(OO, Ol, lO, ll), nrow = 2, ncol = 2)
        
        prediction_data <- prediction_data[, c(2:3)]
        
        # Now we can make the Beta Binomial and Naive Bayes Predictions
        BB.Prob <- Beta_Bin_predict(prediction_data, length(rTaxon))
        NB.Prob <- Naive_Bayes_predict(counts, rTaxonValueA)
        
        # Assign these values to the results table, with a small adjustment for multiple prediction
        # Naive Bayes accuracies and log loss scores are calculated as they're being assigned
        for (j in 1:length(rTaxon)) {
          results[matching_rows[j], c("Beta_Binom_Prob", "Beta_Binom_acc", "Beta_Binom_LL",
                                      "Naive_Bayes_Prob", "Naive_Bayes_acc", "Naive_Bayes_LL")] <-
            c(BB.Prob, calculate_accuracy(rTaxonValueB[j], BB.Prob), LogLoss(rTaxonValueB[j], BB.Prob), 
              NB.Prob[j], calculate_accuracy(rTaxonValueB[j], NB.Prob[j]), LogLoss(rTaxonValueB[j], NB.Prob[j])) 
        }
      }
      # Rewrite the results file with the new tests
      write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
    }
  }
}



##### 
# Run the predictions for random data
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
  rTaxonValueA <- results[matching_rows, 3]
  rTaxonValueB <- results[matching_rows, 4]
  
  # Call the data file for this iteration
  data_name <- paste0(data_path, ".", i, ".Full_data.txt")
  data <- read.table(data_name, skip = 1, sep = "\t")
  
  # Now establish some details we need to predict
  prediction_data <- data[-rTaxon,]
  
  OO <- sum(prediction_data[, 4] == 0)  # This saves as 0 instead of 00
  Ol <- sum(prediction_data[, 4] == 1)  # This saves as 1 instead of 01
  lO <- sum(prediction_data[, 4] == 10) # This saves intuitively
  ll <- sum(prediction_data[, 4] == 11) # This saves intuitively
  counts <- matrix(data = c(OO, Ol, lO, ll), nrow = 2, ncol = 2)
  
  prediction_data <- prediction_data[, c(2:3)]
  
  # Now we can make the Beta Binomial and Naive Bayes Predictions
  BB.Prob <- Beta_Bin_predict(prediction_data, length(rTaxon))
  NB.Prob <- Naive_Bayes_predict(counts, rTaxonValueA)
  
  # Assign these values to the results table
  # Accuracies and Log loss scores are calculated as they're added  
  results[i, c("Beta_Binom_Prob", "Beta_Binom_acc", "Beta_Binom_LL",
               "Naive_Bayes_Prob", "Naive_Bayes_acc", "Naive_Bayes_LL")] <-
    c(BB.Prob, calculate_accuracy(rTaxonValueB, BB.Prob), LogLoss(rTaxonValueB, BB.Prob), NB.Prob,
      calculate_accuracy(rTaxonValueB, NB.Prob), LogLoss(rTaxonValueB, NB.Prob))
}
# Rewrite the results file with the new tests
write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")



### Check and potentially rerun this for Multiple and Clade prediction
if (multiple_prediction == TRUE) {
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
    rTaxonValueA <- results[matching_rows, 3]
    rTaxonValueB <- results[matching_rows, 4]
    
    # Call the data file for this iteration
    data_name <- paste0(data_path, ".", i, ".Full_data.txt")
    data <- read.table(data_name, skip = 1, sep = "\t")
    
    # Now establish some details we need to predict
    prediction_data <- data[-rTaxon,]
    
    OO <- sum(prediction_data[, 4] == 0)  # This saves as 0 instead of 00
    Ol <- sum(prediction_data[, 4] == 1)  # This saves as 1 instead of 01
    lO <- sum(prediction_data[, 4] == 10) # This saves intuitively
    ll <- sum(prediction_data[, 4] == 11) # This saves intuitively
    counts <- matrix(data = c(OO, Ol, lO, ll), nrow = 2, ncol = 2)
    
    prediction_data <- prediction_data[, c(2:3)]
    
    # Now we can make the Beta Binomial and Naive Bayes Predictions
    BB.Prob <- Beta_Bin_predict(prediction_data, length(rTaxon))
    NB.Prob <- Naive_Bayes_predict(counts, rTaxonValueA)
    
    # Now check the accuracy and Log loss scores
    # These will be the same for every prediction, so we only need to do it once
    BB.acc <- 
    
    # Assign these values to the results table, with a small adjustment for multiple prediction
    # Naive Bayes accuracies and log loss scores are calculated as they're being assigned
    for (j in 1:length(rTaxon)) {
      results[matching_rows[j], c("Beta_Binom_Prob", "Beta_Binom_acc", "Beta_Binom_LL",
                                  "Naive_Bayes_Prob", "Naive_Bayes_acc", "Naive_Bayes_LL")] <-
        c(BB.Prob, calculate_accuracy(rTaxonValueB[j], BB.Prob), LogLoss(rTaxonValueB[j], BB.Prob), 
          NB.Prob[j], calculate_accuracy(rTaxonValueB[j], NB.Prob[j]), LogLoss(rTaxonValueB[j], NB.Prob[j])) 
    }
  }
  # Rewrite the results file with the new tests
  write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
}


# Check and run clade prediction
if (clade_prediction == TRUE) {
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
    rTaxonValueA <- results[matching_rows, 3]
    rTaxonValueB <- results[matching_rows, 4]
    
    # Call the data file for this iteration
    data_name <- paste0(data_path, ".", i, ".Full_data.txt")
    data <- read.table(data_name, skip = 1, sep = "\t")
    
    # Now establish some details we need to predict
    prediction_data <- data[-rTaxon,]
    
    OO <- sum(prediction_data[, 4] == 0)  # This saves as 0 instead of 00
    Ol <- sum(prediction_data[, 4] == 1)  # This saves as 1 instead of 01
    lO <- sum(prediction_data[, 4] == 10) # This saves intuitively
    ll <- sum(prediction_data[, 4] == 11) # This saves intuitively
    counts <- matrix(data = c(OO, Ol, lO, ll), nrow = 2, ncol = 2)
    
    prediction_data <- prediction_data[, c(2:3)]
    
    # Now we can make the Beta Binomial and Naive Bayes Predictions
    BB.Prob <- Beta_Bin_predict(prediction_data, length(rTaxon))
    NB.Prob <- Naive_Bayes_predict(counts, rTaxonValueA)
    
    # Now check the accuracy and Log loss scores
    # These will be the same for every prediction, so we only need to do it once
    BB.acc <- calculate_accuracy(rTaxonValueB, BB.Prob)
    BB.LL <- LogLoss(rTaxonValueB, BB.Prob)
    
    # Assign these values to the results table, with a small adjustment for multiple prediction
    # Naive Bayes accuracies and log loss scores are calculated as they're being assigned
    for (j in 1:length(rTaxon)) {
      results[matching_rows[j], c("Beta_Binom_Prob", "Beta_Binom_acc", "Beta_Binom_LL",
                                  "Naive_Bayes_Prob", "Naive_Bayes_acc", "Naive_Bayes_LL")] <-
        c(BB.Prob, calculate_accuracy(rTaxonValueB[j], BB.Prob), LogLoss(rTaxonValueB[j], BB.Prob), 
          NB.Prob[j], calculate_accuracy(rTaxonValueB[j], NB.Prob[j]), LogLoss(rTaxonValueB[j], NB.Prob[j])) 
    }
  }
  # Rewrite the results file with the new tests
  write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
}


##### 
# Check and run the variable rates predictions
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
      rTaxonValueA <- results[matching_rows, 3]
      rTaxonValueB <- results[matching_rows, 4]
      
      # Call the data file for this iteration
      data_name <- paste0(data_path, ".", i, ".Full_data.txt")
      data <- read.table(data_name, skip = 1, sep = "\t")
      
      # Now establish some details we need to predict
      prediction_data <- data[-rTaxon,]
      
      OO <- sum(prediction_data[, 4] == 0)  # This saves as 0 instead of 00
      Ol <- sum(prediction_data[, 4] == 1)  # This saves as 1 instead of 01
      lO <- sum(prediction_data[, 4] == 10) # This saves intuitively
      ll <- sum(prediction_data[, 4] == 11) # This saves intuitively
      counts <- matrix(data = c(OO, Ol, lO, ll), nrow = 2, ncol = 2)
      
      prediction_data <- prediction_data[, c(2:3)]
      
      # Now we can make the Beta Binomial and Naive Bayes Predictions
      BB.Prob <- Beta_Bin_predict(prediction_data, length(rTaxon))
      NB.Prob <- Naive_Bayes_predict(counts, rTaxonValueA)
      
      # Assign these values to the results table
      # Accuracies and Log loss scores are calculated as they're added  
      results[i, c("Beta_Binom_Prob", "Beta_Binom_acc", "Beta_Binom_LL",
                   "Naive_Bayes_Prob", "Naive_Bayes_acc", "Naive_Bayes_LL")] <-
        c(BB.Prob, calculate_accuracy(rTaxonValueB, BB.Prob), LogLoss(rTaxonValueB, BB.Prob), NB.Prob,
          calculate_accuracy(rTaxonValueB, NB.Prob), LogLoss(rTaxonValueB, NB.Prob))
    }
    # Rewrite the results file with the new tests
    write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
    
    
    
    ### Check and potentially rerun this for Multiple and Clade prediction
    if (multiple_prediction == TRUE) {
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
        rTaxonValueA <- results[matching_rows, 3]
        rTaxonValueB <- results[matching_rows, 4]
        
        # Call the data file for this iteration
        data_name <- paste0(data_path, ".", i, ".Full_data.txt")
        data <- read.table(data_name, skip = 1, sep = "\t")
        
        # Now establish some details we need to predict
        prediction_data <- data[-rTaxon,]
        
        OO <- sum(prediction_data[, 4] == 0)  # This saves as 0 instead of 00
        Ol <- sum(prediction_data[, 4] == 1)  # This saves as 1 instead of 01
        lO <- sum(prediction_data[, 4] == 10) # This saves intuitively
        ll <- sum(prediction_data[, 4] == 11) # This saves intuitively
        counts <- matrix(data = c(OO, Ol, lO, ll), nrow = 2, ncol = 2)
        
        prediction_data <- prediction_data[, c(2:3)]
        
        # Now we can make the Beta Binomial and Naive Bayes Predictions
        BB.Prob <- Beta_Bin_predict(prediction_data, length(rTaxon))
        NB.Prob <- Naive_Bayes_predict(counts, rTaxonValueA)
        
        # Now check the accuracy and Log loss scores
        # These will be the same for every prediction, so we only need to do it once
        BB.acc <- calculate_accuracy(rTaxonValueB, BB.Prob)
        BB.LL <- LogLoss(rTaxonValueB, BB.Prob)
        
        # Assign these values to the results table, with a small adjustment for multiple prediction
        # Naive Bayes accuracies and log loss scores are calculated as they're being assigned
        for (j in 1:length(rTaxon)) {
          results[matching_rows[j], c("Beta_Binom_Prob", "Beta_Binom_acc", "Beta_Binom_LL",
                                      "Naive_Bayes_Prob", "Naive_Bayes_acc", "Naive_Bayes_LL")] <-
            c(BB.Prob, calculate_accuracy(rTaxonValueB[j], BB.Prob), LogLoss(rTaxonValueB[j], BB.Prob), 
              NB.Prob[j], calculate_accuracy(rTaxonValueB[j], NB.Prob[j]), LogLoss(rTaxonValueB[j], NB.Prob[j])) 
        }
      }
      # Rewrite the results file with the new tests
      write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
    }
    
    
    # Check and run clade prediction
    if (clade_prediction == TRUE) {
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
        rTaxonValueA <- results[matching_rows, 3]
        rTaxonValueB <- results[matching_rows, 4]
        
        # Call the data file for this iteration
        data_name <- paste0(data_path, ".", i, ".Full_data.txt")
        data <- read.table(data_name, skip = 1, sep = "\t")
        
        # Now establish some details we need to predict
        prediction_data <- data[-rTaxon,]
        
        OO <- sum(prediction_data[, 4] == 0)  # This saves as 0 instead of 00
        Ol <- sum(prediction_data[, 4] == 1)  # This saves as 1 instead of 01
        lO <- sum(prediction_data[, 4] == 10) # This saves intuitively
        ll <- sum(prediction_data[, 4] == 11) # This saves intuitively
        counts <- matrix(data = c(OO, Ol, lO, ll), nrow = 2, ncol = 2)
        
        prediction_data <- prediction_data[, c(2:3)]
        
        # Now we can make the Beta Binomial and Naive Bayes Predictions
        BB.Prob <- Beta_Bin_predict(prediction_data, length(rTaxon))
        NB.Prob <- Naive_Bayes_predict(counts, rTaxonValueA)
        
        # Now check the accuracy and Log loss scores
        # These will be the same for every prediction, so we only need to do it once
        BB.acc <- calculate_accuracy(rTaxonValueB, BB.Prob)
        BB.LL <- LogLoss(rTaxonValueB, BB.Prob)
        
        # Assign these values to the results table, with a small adjustment for multiple prediction
        # Naive Bayes accuracies and log loss scores are calculated as they're being assigned
        for (j in 1:length(rTaxon)) {
          results[matching_rows[j], c("Beta_Binom_Prob", "Beta_Binom_acc", "Beta_Binom_LL",
                                      "Naive_Bayes_Prob", "Naive_Bayes_acc", "Naive_Bayes_LL")] <-
            c(BB.Prob, calculate_accuracy(rTaxonValueB[j], BB.Prob), LogLoss(rTaxonValueB[j], BB.Prob), 
              NB.Prob[j], calculate_accuracy(rTaxonValueB[j], NB.Prob[j]), LogLoss(rTaxonValueB[j], NB.Prob[j]))
        }
      }
      # Rewrite the results file with the new tests
      write.table(results, file = results_name, quote = F, row.names = F, col.names = T, sep = "\t")
    }  
  }
}

print("Finished with all Beta Binomial and Naive Bayes Predictions.")
print("The results can be found in the corresponding files in the Results folder.")