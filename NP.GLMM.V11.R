################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# May 2025

################################################################################
# Script uses functions defined in 'DiscreteFunctions'
# Calls tree and data files from 'TestModels' and 'SimInstructions'
# Gathers these to produce predictions using logistic regression.

# Outputs:
    #'*NP.GLMM.Results.txt' = table with the results of the non-phylogenetic 
    #                           logistic regression prediction method

################################################################################

# This script was written to include non-phylogenetic logistric regression predictions

# First, make sure to call the necessary libraries
library(ape)
library(phytools)
library(MCMCglmm)

# Get the current working directory
original_dir <- getwd()

## Get the types and trials arrays
bash_file <- Sys.glob("Discrete_Simulation.V*")
types_list <- grep("^types=", readLines(bash_file), value = TRUE)
cleaned_string <- gsub("types=\\(|\\)", "", types_list)
cleaned_string <- gsub("\"", "", cleaned_string)
types <- strsplit(cleaned_string, " ")[[1]]

# Number of iterations when running as batch job/ on HPC 
iter_line <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", iter_line))

# Get the sampling type and number
unknown_line <- grep("^unknown_n=", readLines(bash_file), value = TRUE)
unknown_n <- as.numeric(gsub("[^0-9.-]", "", unknown_line))

# Pull necessary functions for this analysis
source(Sys.glob("DiscreteFunctions.V*"))

# Now we can loop through the trials and complete these trials
for (type in types) {

    # First set the directory and create the results matrix
    setwd(type)
    Results.NP.GLMM <- matrix(nrow = 0, ncol = 4)
    colnames(Results.NP.GLMM) <- c("Trial_#", "GLMM_Prob", "GLMM_Acc", "GLMM_LL")

    # Now, loop through the files for each iteration
    for (i in 1:trial_amount) {
        # First, we need to call back the correct data
        # Start by calling the taxon and trait A
        dataname <- paste0(type, ".", i, ".predict_data.txt")
        TempMat <- read.table(file = dataname)
      
        rTaxon <- which(TempMat[, 3] == "?", arr.ind = TRUE)
        rTaxonValueA <- TempMat[rTaxon, 2]

        if (unknown_n > 1) {      
            # Then trait B
            Bname <- paste0("./Results/", type, "TraitBInfo.txt")
            Bdata <- read.table(file = Bname, row.names = NULL)
            Bdata <- Bdata[,-1]
            rTaxonValueB <- as.vector(unlist(Bdata[i, ]))

            # Correct the lengths of rTaxonValueB if it's shorter than unknown_n
            rTaxonValueB <- rTaxonValueB[1:length(rTaxon)]

        } else if (unknown_n == 1) {
            # Now, we need to call the data back
            Bname <- paste0("./Results/", type, ".NP.Results.txt")
            Bdata <- read.table(file = Bname, skip = 1)
            Bdata <- Bdata[, 3:5]
            rTaxonValueB <- Bdata[i, 3]
        }
        # Pull the tree and edit out the unknown(s)
        tree_name <- paste0(type, ".", i, ".full_tree.tre")
        tree_complete <- read.nexus(tree_name)
        tree_train <- drop.tip(tree_complete, rTaxon)        

        # Now, assign the variables to match the language
        edited_TempMat <- TempMat[-rTaxon, ]
        test_species <- rTaxon
        x_test <- rTaxonValueA
        x_train <- edited_TempMat[, 2]
        y_test <- rTaxonValueB
        y_train <- edited_TempMat[, 3]
        trait_test <- data.frame(y = y_test, x = x_test, species = test_species)
        trait_train <- data.frame(
            y = y_train,
            x = x_train,
            species = tree_train$tip.label
        )

        ### "Save" model.
        mcmcglmm_priors <- list(
            B = list(mu = c(0, 0), V = ((1 + (pi^2) / 3) * diag(2))),
            R = list(V = 1, fix = 1)
        )

        mcmcglmm_model <- MCMCglmm(
            fixed = y ~ x,
            family = "categorical",
            data = trait_train,
            prior = mcmcglmm_priors,
            nitt = 8000000,
            thin = 1000,
            burnin = 4000000
            )  # takes about 5 minutes

        mcmcglmm_post <- mcmcglmm_model$Sol

        ### "Load" model.
        e <- rlogis(n = nrow(mcmcglmm_post), location = 0, scale = 1)
        l <- as.vector(c(1, trait_test$x) %*% t(as.matrix(mcmcglmm_post))) + e
        mcmcglmm_ppd <- as.numeric(l > 0)
        mcmcglmm_point_pred_prob <- mean(mcmcglmm_ppd)

        mcmcglmm_log_loss <- LogLoss(
            rTaxonValueB = y_test,
            Predictive_Probability = mcmcglmm_point_pred_prob
        )

        # Calculate the accuracy of this trial
        mcmcglmm_acc <- ifelse((rTaxonValueB == 1 & mcmcglmm_point_pred_prob > 0.5) | (rTaxonValueB == 0 & mcmcglmm_point_pred_prob < 0.5), 1, 0)

        # Save the values as single numbers
        mcmcglmm_point_pred_prob <- mean(mcmcglmm_point_pred_prob)
        mcmcglmm_acc <- mean(mcmcglmm_acc)
        mcmcglmm_log_loss <- median(mcmcglmm_log_loss)

        # Now save this trial
        Trial <- c(i, mcmcglmm_point_pred_prob, mcmcglmm_acc, mcmcglmm_log_loss)
        Results.NP.GLMM <- rbind(Results.NP.GLMM, Trial)
    }

    # Save the matrix as a file for later
    filename <- paste0("./Results/", type, ".NP.GLMM.Results.txt")
    write.table(Results.NP.GLMM, file = filename, quote = F, sep = "\t")

    # Return to the original directory
    setwd(original_dir)
}