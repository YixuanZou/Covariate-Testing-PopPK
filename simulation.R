#-------------------------------------------------------------
# Title: Comparison of Wald, Score and Likelihood Ratio Tests
#        for Covariate Testing in PopPK
# Author: Yixuan Zou
# Create Date: 5/26/2020
# Last Update: 6/13/2020
# CAUTION: You need to put NONMEM batch file into the folder
#-------------------------------------------------------------

#-------------------------------------------------------------
#### r setting
# clear all objects
rm(list=ls(all=TRUE))

# load packages
requiredPackages <- c("ggplot2", "dplyr", "tidyr", "ggthemes")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
#-------------------------------------------------------------

#-------------------------------------------------------------
#### simulation setting (change here)
# number of simulation
num_sim <- 1000

# home directory
home_dir <- "D:/Pharmacy/Projects/Cov_test"
setwd(home_dir)
source("functions.r")

# sample size
sample_size_list <- c(50, 100, 200)

# coefficient of effects of age on clearance
cl_age_list <- c(0.1, 0.25, 0.5)

# coefficient of effects of gender on clearance
cl_sex_list <- c(0.1, 0.25, 0.5)

# sample desgin
sample_design_list <- c('Sparse', 'Intensive')

#-------------------------------------------------------------
# run simulation for each situation
start_time <- Sys.time()
for(sample_size in sample_size_list){
  for(sample_design in sample_design_list){
    for(i in 1:length(cl_age_list)){
      run_simulation(sample_size, sample_design, i,
                     cl_age_list, "AGE", num_sim,
                     home_dir)
      run_simulation(sample_size, sample_design, i,
                     cl_sex_list, "SEX", num_sim,
                     home_dir)
    }
  }
}
run_time <- format(Sys.time() - start_time)
setwd(home_dir)
write.table(run_time, file = "run_time.txt", sep = "", row.names = F, col.names = F)
#-------------------------------------------------------------
# pilot simulation 
# sample_size <- 100
# sample_design <- 'Intensive'
# i <- 2
# num_sim <- 1000
# start_time <- Sys.time()
# run_simulation(sample_size, sample_design, i,
#                cl_age_list, "AGE", num_sim,
#                home_dir)
# run_time <- format(Sys.time() - start_time)
# setwd(home_dir)
# write.table(run_time, file = "run_time.txt", sep = "", row.names = F, col.names = F)
# 
# type_analysis <- c('type_one', 'power')
# for(analysis in type_analysis){
#   run_analysis(sample_size, sample_design, i,
#                cl_age_list, "AGE", num_sim, home_dir,
#                analysis)
# }
