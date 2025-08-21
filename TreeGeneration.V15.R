################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# June 2025

################################################################################
# This will generate phylogenies which will be used to simulated data and predict.
# The output tables can be found in the 'Trees' folder.
# It will produce n trees, where n is the trial amount set in the shell script.

################################################################################



# Start by calling the necessary parameters and settings
# Call the bash file and any necessary function
bash_file <- Sys.glob("Discrete_Simulation.V*")
version_str <- sub(".*(V[0-9]+).*", "\\1", bash_file)
source(paste0("Scripts/DiscreteFunctions.", version_str, ".R"))

# Get the parameters we need to make the trees
pop_size <- grep("^pop_size=", readLines(bash_file), value = TRUE)
pop_size <- as.numeric(gsub("[^0-9.-]", "", pop_size))

# Number of data sets/trials for each combination of settings
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))



#####
# Install and call the necessary packages
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
# Simulate your trees
for (i in 1:trial_amount) {
  # First, create a tree
  full_tree <- ape::rphylo(n = pop_size, birth = 1, death = 0)
  
  # Change the tip labels to be easier to run the tests
  full_tree$tip.label <- seq_along(full_tree$tip.label)
  
  # Next, we have to standardize the tree length
  # To do this, we first get the tree length
  treelength <- as.numeric(sum(full_tree$edge.length))
  
  # Next, we standardize to a total tree length of 700 units
  std.lengths <- lapply(full_tree$edge.length, function(x) x * (700 / treelength))
  
  #Now convert this list back to a numerical vector and replace the old branch lengths
  std.lengths <- as.numeric(unlist(std.lengths))
  full_tree$edge.length <- std.lengths
  
  # Finally, we need to name, save, and assign this tree to use later
  tree_name <- paste0("Trees/Full_tree.", i, ".tre")
  write.nexus(full_tree, file = (tree_name))
}

print(paste0("Simulated ", trial_amount, " trees with ", pop_size, " tips. Found in 'Trees' folder."))
