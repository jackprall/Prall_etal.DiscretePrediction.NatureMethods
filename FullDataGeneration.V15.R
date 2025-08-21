################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# June 2025

################################################################################
# Running this script genereates the full data sets to be used in this simulation
# This took ~24 minutes in our simulation for 1000 data sets for all tests


# Outputs are stored within the folder that corresponds to the test being run
# (e.g. for constant rates, with low, equal q-matrix see /ConstantRates/ER.L/Data)

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

# Get the parameters we need to make the trees
pop_size <- grep("^pop_size=", readLines(bash_file), value = TRUE)
pop_size <- as.numeric(gsub("[^0-9.-]", "", pop_size))

# Get the special cases you may run
variable_rates <- read_logical("variable_rates", bash_file)

# Get rates
L <- grep("^low_rate=", readLines(bash_file), value = TRUE)
M <- grep("^medium_rate=", readLines(bash_file), value = TRUE)
H <- grep("^high_rate=", readLines(bash_file), value = TRUE)

L <- as.numeric(gsub("[^0-9.-]", "", L))
M <- as.numeric(gsub("[^0-9.-]", "", M))
H <- as.numeric(gsub("[^0-9.-]", "", H))

# Get dependency scales and adjustments so that the sum of the rates remains constant
depend_scale1 <- as.numeric(gsub(".*=(.*)", "\\1", grep("^depend_scale1=", readLines(bash_file), value = TRUE)))
depend_adj1   <- as.numeric(gsub(".*=(.*)", "\\1", grep("^depend_adj1=", readLines(bash_file), value = TRUE)))

depend_scale2 <- as.numeric(gsub(".*=(.*)", "\\1", grep("^depend_scale2=", readLines(bash_file), value = TRUE)))
depend_adj2   <- as.numeric(gsub(".*=(.*)", "\\1", grep("^depend_adj2=", readLines(bash_file), value = TRUE)))

# Number of data sets/trials for each combination of settings
trial_amount <- grep("^num_iterations=", readLines(bash_file), value = TRUE)
trial_amount <- as.numeric(gsub("[^0-9.-]", "", trial_amount))

# Get the tree scales
CR_tree_scales <- c(1, 1, 1, 1)
VR_tree_scales <- grep("^tree_scales=", readLines(bash_file), value = TRUE)
VR_tree_scales <- gsub(".*=\\((.*)\\)", "\\1", VR_tree_scales)
VR_tree_scales <- as.numeric(unlist(strsplit(VR_tree_scales, "\\s+")))



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
# Set the q-matrices that we will use for these runs
q1 <- generate_qmat(L, L, 1, 1)   # Low, equal rates
q2 <- generate_qmat(M, M, 1, 1)   # Medium, equal rates  
q3 <- generate_qmat(H, H, 1, 1)   # High, equal rates
q4 <- generate_qmat(M, L, 1, 1)   # Low loss, medium gain
q5 <- generate_qmat(H, L, 1, 1)   # Low loss, high gain
q6 <- generate_qmat(H, M, 1, 1)   # Medium loss, high gain
q7 <- generate_qmat(L, L, depend_scale1, depend_adj1)   # Low rates, lower dependency
q8 <- generate_qmat(M, M, depend_scale1, depend_adj1)   # Medium rates, lower dependency
q9 <- generate_qmat(H, H, depend_scale1, depend_adj1)   # High rates, lower dependency
q10 <- generate_qmat(L, L, depend_scale2, depend_adj2)  # Low rates, higher dependency
q11 <- generate_qmat(M, M, depend_scale2, depend_adj2)  # Medium rates, higher dependency
q12 <- generate_qmat(H, H, depend_scale2, depend_adj2)  # High rates, higher dependency
qRandom <- generate_qmat(L, L, 1, 1) #This is just a placeholder
    


#####
# Generate the data for the study
for (i in 1:trial_amount) {
  treename <- paste0("Trees/Full_tree.", i, ".tre")
  full_tree <- read.nexus(treename)
  
  # Now we need to delve into each type
  for (type in types) {
    # First, set the correct Q-matrix
    if(type == "ER.L") {qmat <- q1}
    if(type == "ER.M") {qmat <- q2}
    if(type == "ER.H") {qmat <- q3} 
    if(type == "DR.LM") {qmat <- q4}
    if(type == "DR.LH") {qmat <- q5}
    if(type == "DR.MH") {qmat <- q6} 
    if(type == "DEP1.L") {qmat <- q7}
    if(type == "DEP1.M") {qmat <- q8}
    if(type == "DEP1.H") {qmat <- q9} 
    if(type == "DEP2.L") {qmat <- q10}
    if(type == "DEP2.M") {qmat <- q11}
    if(type == "DEP2.H") {qmat <- q12}
    


    #####
    # Run the variable rates if necessary
    if (variable_rates == TRUE) {
      
      
      # Run the function that will generate data based on the correct settings
      full_data <- generate_data(full_tree, VR_tree_scales, pop_size, qmat, Random = F)
      colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")

      # Save the data to be used later
      full_data_name <- paste0("VariableRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
      write.table(full_data, full_data_name, quote = F, sep = "\t", row.names = F, col.names = T)
          
      # Do a count of each character state across the tree
      countname <- paste0("VariableRates/", type, "/Data/", type)
      countup(i, full_data, countname)
    }
    


    #####
    # Run the function that will generate data based on the correct settings
    full_data <- generate_data(full_tree, CR_tree_scales, pop_size, qmat, Random = F)
    colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
    
    # Save the data to be used later
    full_data_name <- paste0("ConstantRates/", type, "/Data/", type, ".", i, ".Full_data.txt")
    write.table(full_data, full_data_name, quote = F, sep = "\t", row.names = F, col.names = T)
    
    # Do a count of each character state across the tree
    countname <- paste0("ConstantRates/", type, "/Data/", type)
    countup(i, full_data, countname)
  }
  

  
  #####
  # Finally, run the function that will generate data based on the correct settings
  full_data <- generate_data(full_tree, CR_tree_scales, pop_size, qRandom, Random = T)
  colnames(full_data) <- c("Taxon_#", "Trait_A", "Trait_B", "Class")
  
  # Save the data to be used later
  full_data_name <- paste0("Random/Data/Random.", i, ".Full_data.txt")
  write.table(full_data, full_data_name, quote = F, sep = "\t", row.names = F, col.names = T)
  
  # Do a count of each character state across the tree
  countup(i, full_data, path_name = "Random/Data/Random")
}

print("Finished simulating trait data. They can be found in corresponding 'Data' folders.")
print("(e.g. ConstantRates/ER.L/Data  or  Random/Data)")