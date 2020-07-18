# Covariate Testing PopPK

## Description 
This repository is developed to conduct monte-carlo simulation to compare the emprical power and type I error of Wald, Score and Likelihood Ratio Tests for the covariate testing probelm in population pharmacokinetic (PopPK) model.

- The functions.r contains all the R function that I created to facilitate the simulation process
- The simulation.r contains the simulation scenarios and the main simulation function
- PAT_SIM_DATA.CSV is real subject covariates downloaded from https://wwwn.cdc.gov/nchs/nhanes/
- pilot_sim_results is the pilot simulation results folder. It contains the pilot simulation results when we used sample size 100, intensive sampling design, and age was added on clearance using the power model with the coefficient 0.5. You can easily calculate the emprical power and type I error using the csv files


## Instruction
In order to run the simulation, you need to follow the steps:

1. Download the codes
2. Put nmfe74.bat or nmfe74 (linux) into the folder
3. Go to simulation.r to change the home directory and the simulation scenarios
4. Run simulation.r
 
