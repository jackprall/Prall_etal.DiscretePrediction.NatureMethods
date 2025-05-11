################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# May 2025

################################################################################
# Script uses functions defined in 'DiscreteFunctions'
# Calls tree and data files from 'TestModels' and 'SimInstructions'
# Gathers these to produce predictions using phylogenetic logistic regression.

# Outputs:
    #'*PH.GLMM.Results.txt' = table with the results of the phylogenetic 
    #                           logistic regression prediction method

################################################################################

# This script was written to include phylogenetic logistric regression predictions

# First, make sure to call the ape library
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
    Results.PH.GLMM <- matrix(nrow = 0, ncol = 4)
    colnames(Results.PH.GLMM) <- c("Trial_#", "Phy_GLMM_Prob", "Phy_GLMM_Acc", "Phy_GLMM_LL")

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
        x_test <- as.numeric(rTaxonValueA)
        x_train <- as.numeric(edited_TempMat[, 2])
        y_test <- as.numeric(rTaxonValueB)
        y_train <- as.numeric(edited_TempMat[, 3])
        trait_test <- data.frame(y = y_test, x = x_test, species = test_species)
        trait_train <- data.frame(
            y = y_train,
            x = x_train,
            species = tree_train$tip.label
        )

        ### "Save" model.
        mcmcglmm_phy_priors <- list(
            B = list(mu = c(0, 0), V = ((1 + (pi^2) / 3) * diag(2))),
            R = list(V = 1, fix = 1),
            G = list(G1 = list(V = 1, nu = 1000, alpha.mu = 0, alpha.V = 1))
        )

        mcmcglmm_C_inv <- inverseA(pedigree = tree_train)$Ainv

        mcmcglmm_phy_model <- MCMCglmm(
            fixed = y ~ x,
            random = ~species,
            family = "categorical",
            data = trait_train,
            prior = mcmcglmm_phy_priors,
            ginverse = list(species = mcmcglmm_C_inv),
            nitt = 5000000,
            thin = 2500,
            burnin = 1250000
        )

        mcmcglmm_phy_post <- cbind(mcmcglmm_phy_model$Sol, mcmcglmm_phy_model$VCV[, 1])

    #     ### "Load" model.
        mcmcglmm_phy_post_pred_object <- post_pred_phy_latent(
            tree_full = tree_complete,
            test_species = test_species,
            y_train = trait_train$y,
            x_train = trait_train$x,
            x_test = trait_test$x,
            mcmc_samp = mcmcglmm_phy_post,
            variance = TRUE
        )

        mcmcglmm_phy_ppd <- as.numeric(mcmcglmm_phy_post_pred_object$latent > 0)
        mcmcglmm_phy_point_pred_prob <- mean(mcmcglmm_phy_ppd)

        mcmcglmm_phy_log_loss <- LogLoss(
            rTaxonValueB = y_test,
            Predictive_Probability = mcmcglmm_point_pred_prob
        )

        # Calculate the accuracy of this trial
        mcmcglmm_phy_acc <- ifelse((rTaxonValueB == 1 & mcmcglmm_phy_point_pred_prob > 0.5) | (rTaxonValueB == 0 & mcmcglmm_phy_point_pred_prob < 0.5), 1, 0)

        # Save the values as single numbers
        mcmcglmm_phy_point_pred_prob <- mean(mcmcglmm_phy_point_pred_prob)
        mcmcglmm_phy_acc <- mean(mcmcglmm_phy_acc)
        mcmcglmm_phy_log_loss <- median(mcmcglmm_phy_log_loss)

        # Now save this trial
        Trial <- c(i, mcmcglmm_phy_point_pred_prob, mcmcglmm_phy_acc, mcmcglmm_phy_log_loss)
        Results.PH.GLMM <- rbind(Results.PH.GLMM, Trial)
    }

    # Save the matrix as a file for later
    filename <- paste0("./Results/", type, ".PH.GLMM.Results.txt")
    write.table(Results.PH.GLMM, file = filename, quote = F, sep = "\t")

    # Return to the original directory
    setwd(original_dir)
}
