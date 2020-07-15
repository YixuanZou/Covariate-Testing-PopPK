#-------------------------------------------------------------
# Title: Comparison of Wald, Score and Likelihood Ratio Tests
#        for Covariate Testing in PopPK
# Author: Yixuan Zou
# Create Date: 5/26/2020
# Last Update: 7/12/2020
# CAUTION: You need to put NONMEM batch file into the folder
#-------------------------------------------------------------

#-------------------------------------------------------------
#### r setting
# clear all objects
rm(list=ls(all=TRUE))
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
# clearance and volume correlation
# in our dicussion the correlation can be 0, 0.3
cl_v_corr_list <- c(0)
# correlated covariate
cov_corr_list <- list(c("WT","BSA"),
                      c("WT", "HT"),
                      c("WT", "ALP"))
num_sim <- 10
sample_size_list <- 50
cl_cont_list <- 1
cl_cat_list <- 1
sample_design_list <- "Sparse"
omegasq_list <- 0.1
sigmasq_list <- 0.1
cov_corr_list <- list(c("WT","HT"))
#-------------------------------------------------------------

#-------------------------------------------------------------
#### Part 1: single covariate
# run simulation for each scenario
start_time <- Sys.time()
for(sample_size in sample_size_list){
  for(sample_design in sample_design_list){
    for(index in 1:length(cl_cat_list)){
      for(omegasq in omegasq_list){
        for(sigmasq in sigmasq_list){
          for(cl_v_corr in cl_v_corr_list){
            run_simulation(sample_size, sample_design, index,
                           cl_cont_list, "WT", num_sim, cl_v_corr,
                           omegasq, sigmasq, home_dir)
            run_simulation(sample_size, sample_design, index,
                           cl_cat_list, "SEX", num_sim, cl_v_corr,
                           omegasq, sigmasq, home_dir)
          }
        }
      }
    }
  }
}
run_time <- format(Sys.time() - start_time)
setwd(home_dir)
write.table(run_time, file = "run_time.txt", sep = "", 
            row.names = F, col.names = F)

# analyze the results
for(sample_size in sample_size_list){
  for(sample_design in sample_design_list){
    for(index in 1:length(cl_cat_list)){
      for(omegasq in omegasq_list){
        for(sigmasq in sigmasq_list){
          for (cl_v_corr in cl_v_corr_list) {
            run_analysis(sample_size, sample_design, index,
                         cl_cont_list, "WT", num_sim, cl_v_corr,
                         omegasq, sigmasq, home_dir)
            run_analysis(sample_size, sample_design, index,
                         cl_cat_list, "SEX", num_sim, cl_v_corr,
                         omegasq, sigmasq, home_dir) 
          }
        }
      }
    }
  }
}

#-------------------------------------------------------------

#-------------------------------------------------------------
#### Part 2: correlated covariate
start_time <- Sys.time()
for(sample_size in sample_size_list){
  for(sample_design in sample_design_list){
    for(index in 1:length(cl_cat_list)){
      for(omegasq in omegasq_list){
        for(sigmasq in sigmasq_list){
          for(cl_v_corr in cl_v_corr_list){
            for(cov_corr in cov_corr_list){
              run_simulation_corr(sample_size, sample_design, index,
                                  cl_cont_list, cov_corr, num_sim, 
                                  cl_v_corr, omegasq, sigmasq, home_dir)
            }
          }
        }
      }
    }
  }
}
run_time <- format(Sys.time() - start_time)
setwd(home_dir)
write.table(run_time, file = "run_time_corr.txt", sep = "", 
            row.names = F, col.names = F)

# analyze the results
for(sample_size in sample_size_list){
  for(sample_design in sample_design_list){
    for(index in 1:length(cl_cat_list)){
      for(omegasq in omegasq_list){
        for(sigmasq in sigmasq_list){
          for(cl_v_corr in cl_v_corr_list){
            for(cov_corr in cov_corr_list){
              run_analysis_corr(sample_size, sample_design, index,
                                cl_cont_list, cov_corr, num_sim, 
                                cl_v_corr, omegasq, sigmasq, home_dir)
            }
          }
        }
      }
    }
  }
}