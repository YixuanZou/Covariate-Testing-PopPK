#-------------------------------------------------------------
# Title: Analyze the results for Wald, Score and 
#   Likelihood Ratio Tests for Covariate Testing in PopPK
# Author: Yixuan Zou
# Create Date: 6/14/2020
# Last Update: 6/14/2020
#-------------------------------------------------------------

#-------------------------------------------------------------
# clear all objects
# rm(list=ls(all=TRUE))

# load packages 
requiredPackages <- c()
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
#-------------------------------------------------------------
type_analysis <- c('type_one', 'power')
# analyze the results
for(sample_size in sample_size_list){
  for(sample_design in sample_design_list){
    for(i in 1:length(cl_age_list)){
      for(analysis in type_analysis){
        run_analysis(sample_size, sample_design, i,
                     cl_age_list, "AGE", num_sim, home_dir,
                     analysis)
        run_analysis(sample_size, sample_design, i,
                     cl_sex_list, "SEX", num_sim, home_dir,
                     analysis) 
      }
    }
  }
}