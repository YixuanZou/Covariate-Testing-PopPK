#-------------------------------------------------------------
# Title: Comparison of Wald, Score and Likelihood Ratio Tests
#        for Covariate Testing in PopPK
# Author: Yixuan Zou
# Create Date: 5/26/2020
# Last Update: 6/28/2020
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
num_sim <- 10

# home directory
home_dir <- "D:/Pharmacy/Projects/Cov_test"
setwd(home_dir)
source("functions.r")

# sample size
sample_size_list <- c(50, 100, 200)

# coefficient of effects of continuous variable on clearance
# when coefficient=0 power is calculated
cl_cont_list <- c(0, 0.1, 0.25, 0.5, 0.75, 1)

# coefficient of effects of categorical variable on clearance
cl_cat_list <- c(0, 0.1, 0.25, 0.5, 0.75, 1)

# sample desgin
sample_design_list <- c('Sparse', 'Intensive')

# omega sqaure
omegasq_list <- c(0.03, 0.1)

# sigma sqaure
sigmasq_list <- c(0.03, 0.1)
#-------------------------------------------------------------
# print all the covariate names
demo_df <- read.csv("pat_sim_data.csv")
print(colnames(demo_df))

# run simulation for each situation
start_time <- Sys.time()
for(sample_size in sample_size_list){
  for(sample_design in sample_design_list){
    for(index in 1:length(cl_cat_list)){
      for(omegasq in omegasq_list){
        for(sigmasq in sigmasq_list){
          run_simulation(sample_size, sample_design, index,
                         cl_cont_list, "WT", num_sim,
                         omegasq, sigmasq, home_dir)
          run_simulation(sample_size, sample_design, index,
                         cl_cat_list, "SEX", num_sim,
                         omegasq, sigmasq, home_dir)
        }
      }
    }
  }
}
run_time <- format(Sys.time() - start_time)
setwd(home_dir)
write.table(run_time, file = "run_time.txt", sep = "", 
            row.names = F, col.names = F)
#-------------------------------------------------------------
# analyze the results
for(sample_size in sample_size_list){
  for(sample_design in sample_design_list){
    for(index in 1:length(cl_cat_list)){
      for(omegasq in omegasq_list){
        for(sigmasq in sigmasq_list){
          run_analysis(sample_size, sample_design, index,
                       cl_cont_list, "WT", num_sim,
                       omegasq, sigmasq, home_dir)
          run_analysis(sample_size, sample_design, index,
                       cl_cat_list, "SEX", num_sim,
                       omegasq, sigmasq, home_dir)
        }
      }
    }
  }
}
#-------------------------------------------------------------
# pilot simulation 
# sample_size <- 100
# sample_design <- 'Intensive'
# i <- 0
# num_sim <- 1000
# cov_name <- "WT"
# omegasq <- sigmasq <- 0.1
# start_time <- Sys.time()
# run_simulation(sample_size, sample_design, i,
#                cl_cont_list, "WT", num_sim,
#                omegasq, sigmasq, home_dir)
# run_time <- format(Sys.time() - start_time)
# setwd(home_dir)
# write.table(run_time, file = "run_time.txt", sep = "", row.names = F, col.names = F)
# 
# run_analysis(sample_size, sample_design, i,
#              cl_cont_list, "WT", num_sim,
#              omegasq, sigmasq, home_dir)
