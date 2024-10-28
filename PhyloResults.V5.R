# This script summarizes the results of all the BayesTraits trials
# writes one table in format .txt for each trial type

# load ape
library(ape)

# Get the current working directory
original_dir <- getwd()

# Set the trial types
Ind.types <- c("ER.L", "ER.M", "ER.H", "DR.LH", "DR.LM", "DR.MH", "Random.Ind")
Dep.types <- c("DEP.M", "DEP.H", "DEP.X", "Random.Dep")

# Number of iterations when running as batch job/ on HPC 
iter_line <- grep("num_iterations=", readLines("Discrete_Simulation.V5.sh.BAK"), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", iter_line))

# Loading in functions: 'Calculate_Post_Prob' and 'calculate_PH_accuracy'
source("DiscreteFunctions.V5.R")


################################################################################
# This loop summarizes the independent trials
for (type in Ind.types) {
    # Change directory to the current model type
    setwd(type)
    
    # Now, we begin compiling the results
    Results.PH <- matrix(nrow = 0, ncol = 5)

    colnames(Results.PH) <- c("Trial_#", "Taxon_#", "Trait_A",  
                              "Phylo_Prob","Phylo_Accuracy")

    # Then go through each iteration of that trial and gather necessary data
    for (i in 1:trial_amount) {
        # Now, we need to call the data back
        dataname <- paste0("./Results/", type, ".NP.Results.txt")
        data <- read.table(file = dataname, skip = 1)

        # Next, we identify the data we need
        rTaxon <- data[i, 3]
        rTaxonValueA <- data[i, 4]

        # Calculate the average prediction from the independent model
        Calculate_Post_Prob(i, rTaxon, type, IndependentCharacters = TRUE)

        # Record the accuracy and precision of this iteration
        calculate_PH_accuracy(rTaxonValueA, IndependentCharacters = TRUE, Ind.Post.Prob)

        # Then we save this into the results file.
        Trial <- c(i, rTaxon, rTaxonValueA, Ind.Post.Prob, Ind.acc)
        Results.PH <- rbind(Results.PH, Trial)
    }

    # Switch to a special results folder
    setwd("./Results")

    # Name and save the file
    savename <- paste0(type, ".PH.Results.txt")
    write.table(Results.PH, file = savename, quote = FALSE, sep = "\t", 
                na = "NA", row.names = FALSE, col.names = TRUE)

    # Now, we combine the Results into one final file
    nname <- paste0(type, ".NP.Results.txt")
    pname <- paste0(type, ".PH.Results.txt")
  
    ndf <- read.table(nname, skip = 1, sep = "\t")
    pdf <- read.table(pname, skip = 1, sep = "\t")

    ndf <- ndf[, 2:7]
    pdf <- pdf[, 4:5]

    Results.Full <- cbind(ndf, pdf)

    colnames(Results.Full) <- c("Trial_#", "Taxon_#", "Trait_A",
                                "Terminal_Branch_Length", "BB_Prob", 
                                "BB_Accuarcy", "Phylo_Prob", "Phylo_Accuracy")
    
  #Don't forget to save the results table.
  #First, name the results file.
  filename <- paste0(type, ".Results.Full.txt")
  
  #Next save the file.
  write.table(Results.Full, file = (filename), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  setwd(original_dir)
}


################################################################################
# Now we tweak the loop to account for the differences of a dependent run
for (type in Dep.types) {
    #Change directory to the type
    setwd(type)
    
    # Now, we begin compiling the results
    Results.PH <- matrix(nrow = 0, ncol = 8)

    colnames(Results.PH) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B", 
                              "Ind_Prob","Ind_Accuracy", "Dep_Prob", 
                              "Dep_Accuracy")

    # Then go through each iteration of that trial and gather necessary data
    for (i in 1:trial_amount) {
        # Now, we need to call the data back
        dataname <- paste0("./Results/", type, ".NP.Results.txt")
        data <- read.table(file = dataname, skip = 1)

        # Next, we identify the data we need
        rTaxon <- data[i, 3]
        rTaxonValueA <- data[i, 4]
        rTaxonValueB <- data[i, 5]

        # Calculate the average prediciton from the phylogenetic methods
        Calculate_Post_Prob(i, rTaxon, type, IndependentCharacters = FALSE)

        # Record the accuracy and precision of this iteration
        calculate_PH_accuracy(rTaxonValueB, IndependentCharacters = FALSE, Ind.Post.Prob, Dep.Post.Prob)

        # Then we save this into the results file.
        Trial <- c(i, rTaxon, rTaxonValueA, rTaxonValueB, Ind.Post.Prob, Ind.acc, Dep.Post.Prob, Dep.acc)
        Results.PH <- rbind(Results.PH, Trial)
    }

    # Switch to a special results folder
    setwd("./Results")

    # Name and save the file
    savename <- paste0(type, ".PH.Results.txt")
    write.table(Results.PH, file = savename, quote = FALSE, sep = "\t", 
                na = "NA", row.names = FALSE, col.names = TRUE)

    # Now, we combine the Results into one final file
    nname <- paste0(type, ".NP.Results.txt")
    pname <- paste0(type, ".PH.Results.txt")
  
    ndf <- read.table(nname, skip = 1, sep = "\t")
    pdf <- read.table(pname, skip = 1, sep = "\t")

    ndf <- ndf[, 2:10]
    pdf <- pdf[, 5:8 ]

    Results.Full <- cbind(ndf, pdf)

    colnames(Results.Full) <- c("Trial_#", "Taxon_#", "Trait_A", "Trait_B",
                                "Terminal_Branch_Length", "Beta_Bin_Prob", 
                                "Beta_Bin_Accuracy", "Naive_Prob", 
                                "Naive_Accuracy", "Ind_Prob", "Ind_Accuracy",
                                "Dep_Prob", "Dep_Accuracy")
    
    #First, name the final results file.
    filename <- paste0(type, ".Results.Full.txt")
    
    #Lastly, save the file and reset the directory.
    write.table(Results.Full, file = (filename), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    setwd(original_dir)
}