################################################################################
# Written by Jack Prall, Liam Feigin, Emerald Bender and Chris Organ
# October 2024

################################################################################
# Running the script below will test multiple phylogenetic and non-phylogenetic prediction models. 
# Calls the files creates by 'TestModels'
# Runs BayesTraits phylogenetic analyses
# Calls and runs 'TestModels' and 'CompileResults' in R
# Written for the Tempest Research Cluster at Montana State University

# Outputs will include any files made by BayesTraits and three Results files
# Results are stored within a 'Results' folder found within the folder for each matrix.

################################################################################

#!/bin/bash
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.

#SBATCH --account=priority-chrisorgan   # The class account
#SBATCH --partition=priority               # queue partition to run the job in
#SBATCH --ntasks=1                         # number of mpi processes (1 for multi-threaded tasks)
#SBATCH --threads-per-core=2               # Turn off hyperthreading
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G
#SBATCH --time=7-00:00:00                 # maximum job run time in days-hours:minutes:seconds
#SBATCH --job-name=JOB.BayesTraits
#SBATCH --output=BT-%j.out                 # standard output from job
#SBATCH --error=BT-%j.err                  # standard error from job
#SBATCH --mail-user=john.prall@msu.montana.edu
#SBATCH --mail-type=NONE
#SBATCH --no-requeue


## Set the number of iterations
num_iterations=1000
## Set the number of taxa
pop_size=500
## Define the types and trials arrays
types=("ER.L" "ER.M" "ER.H" "DR.LH" "DR.LM", "DR.MH", "D19.M", "D19.H", "D19.V", "D199.M", "D199.H", "D199.V", "Random")

## First, take note of the current directory, so we can return to it later.
original_dir=$(pwd)

## Loop over each type and create their directories
for type in "${types[@]}"; do
  ## Create a directory and name it
  directoryname="${type}"
  mkdir -p "$directoryname"
done

## Load R and BayesTraits into Tempest
module load R
module load BayesTraits/4.1.3-Linux

## Call the two scripts that define the necessary functions
Rscript DiscreteFunctions.V6.R
Rscript TestModels.V6.R

## Call the script that will run the whole simulation
Rscript SimInstructions.V6.R

## Function to run BayesTraits for both independent and dependent tests
run_rates() {
  local i=$1
  local trial_name=$2

  ## Independent test
  trial_number="${trial_name}.${i}"
  treename="${trial_number}.edited_tree.tre"
  dataname="${trial_number}.edited_data.txt"
  settings="${trial_name}.Ind.${i}.Rates.In.txt"
  BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"

  ## Dependent test
  settings="${trial_name}.Dep.${i}.Rates.In.txt"
  BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"
}

## Function to run BayesTraits for both independent and dependent tests
run_prediction() {
  local i=$1
  local trial_name=$2

  ## Independent test
  trial_number="${trial_name}.${i}"
  treename="${trial_number}.full_tree.tre"
  dataname="${trial_number}.predict_data.txt"
  settings="${trial_name}.Ind.$i.Predict.In.txt"
  BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"

  ## Dependent test
  settings="${trial_name}.Dep.${i}.Predict.In.txt"
  BayesTraitsV4 "${treename}" "${dataname}" <"${settings}"
}

export -f run_rates
export -f run_prediction


## Then we do all the phylogenetic stuff
## Loop over each type
for type in "${types[@]}"; do
  directoryname="${type}"
  cd "./$directoryname" || { echo "Failed to change to directory: $directoryname"; exit 1; }

  ## Set the trial name
  trial_name="${type}"

  ## Placeholder for additional operations with $trial_name if needed
  echo "Processing trial: $trial_name"

  ## Run the rate calculation jobs in parallel
  for i in $(seq 1 $num_iterations); do
    run_rates $i $trial_name &
  done

  ## Wait for all background jobs to complete
  wait

  ## Run the prediction jobs in parallel
  for i in $(seq 1 $num_iterations); do
    run_prediction $i $trial_name &
  done

  ## Wait for all background jobs to complete
  wait

  ## Change back to the original directory
  cd "$original_dir"

done

## This script will compile and summarize all our results
Rscript CompileResults.V6.R