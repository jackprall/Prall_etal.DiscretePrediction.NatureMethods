################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# October 2024

################################################################################
# Running below script will create the tree and trait files passed to BayesTraits across models 'types'
# Calls functions defined in 'TestModels' and 'DiscreteFunctions'
# BayesTraits phylogenetic analyses called in external Bash job script 'Discrete_Simulation'

# Outputs:
    #'*edited_tree.tree' = nexus format phylogeny with number of tips = pop_size-1, IE one OTU randomly dropped
    #'*edited_data.txt = table with discrete trait data excluding tip dropped from edited_tree, 
    #         column 1 = taxon names; column 2 = trait A; column 3 = trait B
    #'*full_tree.tree'  = nexus format phylogeny with number of tips = pop_size
    #'*predict_data.txt' = table with discrete trait data; rows match tips of '*full_tree.nex'
    #         rightmost trait column replaced with '?' for taxon dropped from '*edited_tree.nex'
################################################################################

# Get the current working directory
original <- getwd()

# Get simulation bash file name and then parameters: population size, number of trials, rates, and dependency
bash_file <- Sys.glob("Discrete_Simulation.V*")

# Get population size (number of taxa)
pop_line <- grep("^pop_size=", readLines(bash_file), value = TRUE)
pop_size <- as.numeric(gsub("[^0-9.-]", "", pop_line))

## Get the types and trials arrays
types_list <- grep("^types=", readLines(bash_file), value = TRUE)
cleaned_string <- gsub("types=\\(|\\)", "", types_list)
cleaned_string <- gsub("\"", "", cleaned_string)
types <- strsplit(cleaned_string, " ")[[1]]

# Get rates
low <- grep("^low_rate=", readLines(bash_file), value = TRUE)
med <- grep("^medium_rate=", readLines(bash_file), value = TRUE)
hi <- grep("^high_rate=", readLines(bash_file), value = TRUE)
L <- as.numeric(gsub("[^0-9.-]", "", low))
M <- as.numeric(gsub("[^0-9.-]", "", med))
H <- as.numeric(gsub("[^0-9.-]", "", hi))

# Get dependency scale and adjustment so that the sum of the rates remains constant
dep_scale <- grep("^depend_scale=", readLines(bash_file), value = TRUE)
depend_scale <- as.numeric(gsub("[^0-9.-]", "", dep_scale))
dep_adj <- grep("^depend_adj=", readLines(bash_file), value = TRUE)
depend_adj <- as.numeric(gsub("[^0-9.-]", "", dep_adj))

# Number of iterations when running as batch job/ on HPC
iter_line <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", iter_line))

# Set CRAN mirror & install packages
options(repos = list(CRAN = "https://cloud.r-project.org"))
install.packages(c("ape", "ggnewscale", "paletteer", "phytools", "tidytree", "TreeSimGM", "BiocManager"))
BiocManager::install("ggtree")

# Call the other scripts that define the functions
source(Sys.glob("TestModels.V*"))
source(Sys.glob("DiscreteFunctions.V*"))

library(ape)
library(abind)
library(dplyr)
library(ggtree)
library(ggplot2)
library(ggnewscale)
library(paletteer)
library(phytools)
library(tidytree)
library(TreeSimGM)

# Set a for loop for each trial type
for (type in types) {

  # Create the directory for the current type if it doesn't exist
  dir_name <- paste0("./", type)
  if (!dir.exists(dir_name)) {
    dir.create(dir_name, recursive = TRUE)
  }
  
  # Navigate to the current type directory
  setwd(dir_name)
  
  # Create a results folder if it doesn't exist
  if (!dir.exists("Results")) {
    dir.create("Results")
  }

# Check if the directory name contains "ER.L"
if (type == "ER.L") {
    # A check to see where the run is.
    print("Low equal rates model.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, L, L, 0,
                            L, 0, 0, L,
                            L, 0, 0, L,
                            0, L, L, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "ER.L")


# Check if the directory name contains "ER.M"
} else if (type == "ER.M") {
    # A check to see where the run is.
    print("Medium equal rates model.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, M, M, 0,
                            M, 0, 0, M,
                            M, 0, 0, M,
                            0, M, M, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "ER.M")


# Check if the directory name contains "ER.H"
} else if (type == "ER.H") {
    # A check to see where the run is.
    print("High equal rates model.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, H, H, 0,
                            H, 0, 0, H,
                            H, 0, 0, H,
                            0, H, H, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "ER.H")


# Check if the directory name contains "DR.LM"
} else if (type == "DR.LM") {
    # A check to see where the run is.
    print("Directional rates model, with low loss and medium gain.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, M, M, 0,
                            L, 0, 0, M,
                            L, 0, 0, M,
                            0, L, L, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DR.LM")     

                   
# Check if the directory name contains "DR.LH"
} else if (type == "DR.LH") {
    # A check to see where the run is.
    print("Directional rates model, with low loss and high gain.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, H, H, 0,
                            L, 0, 0, H,
                            L, 0, 0, H,
                            0, L, L, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DR.LH")  

                   
# Check if the directory name contains "DR.MH"
} else if (type == "DR.MH") {
    # A check to see where the run is.
    print("Directional rates model, with medium loss and high gain.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, H, H, 0,
                            M, 0, 0, H,
                            M, 0, 0, H,
                            0, M, M, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DR.MH")                 


# Check if the directory name contains "DEP.L"
} else if (type == "DEP.L") {
    # A check to see where the run is.
    print("Low base rate, 1-order magnitude lower into 0/1.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, L*depend_adj, L*depend_scale, 0,
                            L*depend_adj, 0, 0, L*depend_adj,
                            L*depend_adj, 0, 0, L*depend_adj,
                            0, L*depend_adj, L*depend_scale, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DEP.L")


# Check if the directory name contains "DEP.M"
} else if (type == "DEP.M") {
    # A check to see where the run is.
    print("Medium base rate, 1-order magnitude lower into 0/1.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, M*depend_adj, M*depend_scale, 0,
                            M*depend_adj, 0, 0, M*depend_adj,
                            M*depend_adj, 0, 0, M*depend_adj,
                            0, M*depend_adj, M*depend_scale, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DEP.M")
                   
                   
# Check if the directory name contains "DEP.H"
} else if (type == "DEP.H") {
    # A check to see where the run is.
    print("High base rate, 1-order magnitude lower into 0/1.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, H*depend_adj, H*depend_scale, 0,
                            H*depend_adj, 0, 0, H*depend_adj,
                            H*depend_adj, 0, 0, H*depend_adj,
                            0, H*depend_adj, H*depend_scale, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DEP.H")


# Check if the trial type doesn't exist."
} else if (type == "Random") {
    # A check to see where the run is.
    print("Beginning randomized trials.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, 100, 100, 0,
                            100, 0, 0, 100,
                            100, 0, 0, 100,
                            0, 100, 100, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #Then we will run this four times, once for each combination of the following:
    #Ultrametric, non-ultrametric, symmetric, and asymmetric.

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "Random",  Random = TRUE)                    
} else {
    print("This directory does not exist.")
}
    # Then return to the original directory
    setwd(original)
}