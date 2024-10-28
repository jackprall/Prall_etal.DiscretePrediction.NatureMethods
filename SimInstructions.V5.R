################################################################################
# Written by Jack Prall, Liam Feigin, and Chris Organ
# October 2024

################################################################################
# Running below script will create the tree and trait files passed to BayesTraits across models 'types'
# Calls functions defined in 'NonPhyloOneScript' and 'DiscreteFunctions'
# BayesTraits phylogenetic analyses called in external Bash job script 'Discrete_Simulation'

# Outputs:
    #'*edited_tree.tree' = nexus format phylogeny with number of tips = pop_size-1, IE one OTU randomly dropped
    #'*edited_data.txt = table with discrete trait data excluding tip dropped from edited_tree, 
    #         column 1 = taxon names; column 2 = trait A; column 3 = trait B
    #'*full_tree.tree'  = nexus format phylogeny with number of tips = pop_size
    # '*predict_data.txt' = table with discrete trait data; rows match tips of '*full_tree.nex'
    #         rightmost trait column replaced with '?' for taxon dropped from '*edited_tree.nex'
################################################################################


# First, we'll determine the population size and trial amount for the simulation
# Number of tips when running as batch job/ on HPC 
pop_line<- grep("pop_size=", readLines("Discrete_Simulation.V5.Bender.sh.BAK"), value = TRUE)
pop_size <- as.numeric(gsub("[^0-9.-]", "", pop_line))
# Number of iterations when running as batch job/ on HPC 
iter_line <- grep("num_iterations=", readLines("Discrete_Simulation.V5.sh.BAK"), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", iter_line))

# Set CRAN mirror
options(repos = list(CRAN = "https://cloud.r-project.org"))

# Call the other scripts that define the functions
source("NonPhyloOneScript.V5.R")
source("DiscreteFunctions.V5.R")

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

# Get the current working directory
original <- getwd()

## Define the types and trials arrays
types <- c("ER.L", "ER.M", "ER.H", "DR.LH", "DR.LM", "DR.MH", "DEP.M", "DEP.H", "DEP.X", "Random")

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
    print("Beginning low, equal rates trials.")

    #First, we'll examine equal rates models.
    #Beginning with low, equal rates
    qmat <- matrix(data = c(0, 0.02, 
                            0.02, 0), 
                   nrow = 2, ncol = 2, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "ER.L", 
                   IndependentCharacters = TRUE)


# Check if the directory name contains "ER.M"
} else if (type == "ER.M") {
    # A check to see where the run is.
    print("Beginning medium, equal rates trials.")

    #Next, we'll examine medium, equal rates
    qmat <- matrix(data = c(0, 0.05, 0.05, 0), nrow = 2, ncol = 2, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "ER.M", 
                   IndependentCharacters = TRUE)


# Check if the directory name contains "ER.H"
} else if (type == "ER.H") {
    # A check to see where the run is.
    print("Beginning high, Equal rates trials.")

    #Finally, we'll examine high, equal rates
    qmat <- matrix(data = c(0, 0.1, 0.1, 0), nrow = 2, ncol = 2, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "ER.H", 
                   IndependentCharacters = TRUE)


# Check if the directory name contains "DR.LM"
} else if (type == "DR.LM") {
    # A check to see where the run is.
    print("Beginning medium gain, low loss rates trials.")

    #These are the matrices for when gain and loss are inequal
    #First, low loss, medium gain rates.
    qmat <- matrix(data = c(0, 0.05, 0.02, 0), nrow = 2, ncol = 2, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DR.LM", 
                   IndependentCharacters = TRUE)                        


# Check if the directory name contains "DR.LH"
} else if (type == "DR.LH") {
    # A check to see where the run is.
    print("Beginning high gain, low loss rates trials.")

    # Next, low loss, high gain rates.
    qmat <- matrix(data = c(0, 0.1, 0.02, 0), nrow = 2, ncol = 2, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DR.LH", 
                   IndependentCharacters = TRUE)                          


# Check if the directory name contains "DR.MH"
} else if (type == "DR.MH") {
    # A check to see where the run is.
    print("Beginning high gain, medium loss rates trials.")

    # Last, medium loss, high gain rates.
    qmat <- matrix(data = c(0, 0.1, 0.05, 0), nrow = 2, ncol = 2, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DR.MH", 
                   IndependentCharacters = TRUE)  


# Next, we have our dependent matrices
# Check if the directory name contains "DEP.M"
} else if (type == "DEP.M") {
    # A check to see where the run is.
    print("Beginning medium, dependent rate trials.")

    #This is for medium rates, with low dependency
    qmat <- matrix(data = c(0, 0.005, 0.095, 0,
                            0.095, 0, 0, 0.005,
                            0.005, 0, 0, 0.095,
                            0, 0.095, 0.005, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DEP.M", 
                   IndependentCharacters = FALSE)  


# Check if the directory name contains "DEP.H"
} else if (type == "DEP.H") {
    # A check to see where the run is.
    print("Beginning high, dependent rate trials.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, 0.01, 0.19, 0,
                            0.19, 0, 0, 0.01,
                            0.01, 0, 0, 0.19,
                            0, 0.19, 0.01, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DEP.H", 
                   IndependentCharacters = FALSE) 

# Check if the directory name contains "DEP.X"
} else if (type == "DEP.X") {
    # A check to see where the run is.
    print("Beginning extreme, dependent rate trials.")

    #This is for high rates, with low dependency
    qmat <- matrix(data = c(0, 0.04, 0.76, 0,
                            0.76, 0, 0, 0.04,
                            0.04, 0, 0, 0.76,
                            0, 0.76, 0.04, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "DEP.X", 
                   IndependentCharacters = FALSE) 
                   
# Check if the trial type doesn't exist."
} else if (type == "Random") {
    # A check to see where the run is.
    print("Beginning extreme, dependent rate trials.")

    #This is an arbitrary maxtrix, as one is required to run the code.
    #Any data simulated by evolution will be randomly replaced later.
    qmat <- matrix(data = c(0, 100, 100, 0,
                            100, 0, 0, 100,
                            100, 0, 0, 100,
                            0, 100, 100, 0), 
                            nrow = 4, ncol = 4, byrow = TRUE)

    #First, we'll run an ultrametric, symmetric trial
    Test_NP_Models(pop_size, trial_amount, qmat, ultrametric = TRUE, 
                   symmetric = TRUE, trial_name = "Random", 
                   IndependentCharacters = FALSE, Random = TRUE)                    
} else {
    print("This directory does not exist.")
}
    # Then return to the original directory
    setwd(original)
}