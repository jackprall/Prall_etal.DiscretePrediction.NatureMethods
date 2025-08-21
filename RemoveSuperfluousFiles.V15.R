################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# June 2025

################################################################################
# This script is to delete any superfluous files at the end of the simulation
# This script is optional, as you may want to review these files

# Deleted files include the following
#   All schedule files from BayesTraits
#   All instruction files for BayesTraits
#   All ModelFile files, which are only readable to a machine (BayesTraits)
#   The data files that were edited for model calculates (Rates)
#   The log files from the first BayesTraits run (Rates)

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

# Set the files you wish to delete by the parts of their names
  # ".Schedule." files are not used in this analysis
  # ".Rates." and ".rates_" files are only used for the initial model building
  # ".ModelFile." is the evolutionary model built, but isn't readable to humans
  # ".In." files are instructions for BayesTraits, you may wish to review these
patterns <- c(".Schedule.", ".Rates.", ".rates_", ".ModelFile.", ".In.")
filesToRemove <- paste(gsub("\\.", "\\\\.", patterns), collapse = "|")

# Finally, set the original directory and file paths
original_dir <- getwd()
constant_path <- file.path(original_dir, "ConstantRates")
random_path   <- file.path(original_dir, "Random")
if (isTRUE(variable_rates)) {variable_path <- file.path(original_dir, "VariableRates")}



#####
# Start with Single prediction, specifically with constant rates
for (type in types) {
  # Get to the constant rates folder
  folder_path <- file.path(constant_path, type, "Single")
  setwd(folder_path)
  
  # Get all files in the working directory
  files <- list.files()
  
  # Filter files that match any of the patterns
  target_files <- files[grepl(filesToRemove, files)]
  
  # Delete the matched files
  file.remove(target_files)
  
  

  #####
  # Now check variable rates
  if (isTRUE(variable_rates)) {
    # Get to the constant rates folder
    folder_path <- file.path(variable_path, type, "Single")
    setwd(folder_path)
    
    # Get all files in the working directory
    files <- list.files()
    
    # Filter files that match any of the patterns
    target_files <- files[grepl(filesToRemove, files)]
    
    # Delete the matched files
    file.remove(target_files)
  }
}

# Get to the random rates folder
folder_path <- file.path(random_path, "Single")
setwd(folder_path)

# Get all files in the working directory
files <- list.files()

# Filter files that match any of the patterns
target_files <- files[grepl(filesToRemove, files)]

# Delete the matched files
file.remove(target_files)



#####
# Now we repeat the same for multiple prediction
if (isTRUE(multiple_prediction)) {
  for (type in types) {
    # Get to the constant rates folder
    folder_path <- file.path(constant_path, type, "Multiple")
    setwd(folder_path)
    
    # Get all files in the working directory
    files <- list.files()
    
    # Filter files that match any of the patterns
    target_files <- files[grepl(filesToRemove, files)]
    
    # Delete the matched files
    file.remove(target_files)
    
    
    #####
    # Now check variable rates
    if (isTRUE(variable_rates)) {
      # Get to the constant rates folder
      folder_path <- file.path(variable_path, type, "Multiple")
      setwd(folder_path)
      
      # Get all files in the working directory
      files <- list.files()
      
      # Filter files that match any of the patterns
      target_files <- files[grepl(filesToRemove, files)]
      
      # Delete the matched files
      file.remove(target_files)
    }
  }
  
  # Get to the random rates folder
  folder_path <- file.path(random_path, "Multiple")
  setwd(folder_path)
  
  # Get all files in the working directory
  files <- list.files()
  
  # Filter files that match any of the patterns
  target_files <- files[grepl(filesToRemove, files)]
  
  # Delete the matched files
  file.remove(target_files)
}



#####
# Now we repeat the same for multiple prediction
if (isTRUE(clade_prediction)) {
  for (type in types) {
    # Get to the constant rates folder
    folder_path <- file.path(constant_path, type, "Clade")
    setwd(folder_path)
    
    # Get all files in the working directory
    files <- list.files()
    
    # Filter files that match any of the patterns
    target_files <- files[grepl(filesToRemove, files)]
    
    # Delete the matched files
    file.remove(target_files)
    
    
    #####
    # Now check variable rates
    if (isTRUE(variable_rates)) {
      # Get to the constant rates folder
      folder_path <- file.path(variable_path, type, "Clade")
      setwd(folder_path)
      
      # Get all files in the working directory
      files <- list.files()
      
      # Filter files that match any of the patterns
      target_files <- files[grepl(filesToRemove, files)]
      
      # Delete the matched files
      file.remove(target_files)
    }
  }
  
  # Get to the random rates folder
  folder_path <- file.path(random_path, "Clade")
  setwd(folder_path)
  
  # Get all files in the working directory
  files <- list.files()
  
  # Filter files that match any of the patterns
  target_files <- files[grepl(filesToRemove, files)]
  
  # Delete the matched files
  file.remove(target_files)
}

# Now tell the user that this process is done
setwd(original_dir)
print("Finished deleting superfluous files.")
