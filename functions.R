#-------------------------------------------------------------
# Title: R functions for simulation study
# Author: Yixuan Zou
# Create Date: 5/26/2020
# Last Update: 7/12/2020
#-------------------------------------------------------------

#-------------------------------------------------------------
# load packages
requiredPackages <- c("ggplot2", "dplyr", "tidyr", "ggthemes")
for(p in requiredPackages){
  if(!require(p,character.only = TRUE)) install.packages(p)
  library(p,character.only = TRUE)
}
#-------------------------------------------------------------
# create the template dataset for the simulation
create_data_template <- function(sample_design, sample_size, cov_name){
  ID <- c()
  TIME <- c()
  EVID <- c()
  AMT <- c()
  RATE <- c()
  DV <- c()
  if (startsWith(sample_design, "Intensive")){
    # keep track of time from t=0, 2, 4, 8, 12, 24
    t <- c(0, 2, 4, 8, 12, 24)
    n <- length(t)
    for(i in 1:sample_size){
      ID[seq((i-1)*n+1, i*n)] <- rep(i, n)
      TIME[seq((i-1)*n+1, i*n)] <- c(0, 2, 4, 8, 12, 24)
      EVID[seq((i-1)*n+1, i*n)] <- 0
      EVID[(i-1)*n+1] <- 1
      AMT[seq((i-1)*n+1, i*n)] <- 0
      AMT[(i-1)*n+1] <- 1000
      RATE[seq((i-1)*n+1, i*n)] <- 0
      DV[seq((i-1)*n+1, i*n)] <- 0
    }
  } else{
    # keep track of time from 3 sampling points t
    t <- rep(0, 3)
    n <- length(t)
    for(i in 1:sample_size){
      ID[seq((i-1)*n+1, i*n)] <- rep(i, n)
      TIME[seq((i-1)*n+1, i*n)] <- c(0, sample(1:12, 1), sample(13:24, 1))
      EVID[seq((i-1)*n+1, i*n)] <- 0
      EVID[(i-1)*n+1] <- 1
      AMT[seq((i-1)*n+1, i*n)] <- 0
      AMT[(i-1)*n+1] <- 1000
      RATE[seq((i-1)*n+1, i*n)] <- 0
      DV[seq((i-1)*n+1, i*n)] <- 0
    }
  }
  data <- data.frame(ID, TIME, EVID, AMT, RATE, DV)
  pat_data <- read.csv('PAT_SIM_DATA.CSV')
  pat_sub_data <- sample_n(pat_data, sample_size)
  pat_sub_data$ID <- 1:sample_size
  cov_median <- median(pat_sub_data[, colnames(pat_sub_data) == cov_name])
  final_data <- merge(x = data, y = pat_sub_data, by = "ID", all.x = TRUE)
  write.table(final_data, 'FULL_SIM_DATA.CSV', sep = ',', row.names = FALSE, col.names = FALSE)
  write.table(cov_median, 'COV_MEDIAN.TXT', sep = '', row.names = FALSE, col.names = FALSE)
}

# create the NONMEM script of the base model
create_baseNM <- function(script){
    cov_names <- paste(colnames(read.csv('PAT_SIM_DATA.CSV'))[-1], collapse = " ")
    base_model_ctl <- c("$PROBLEM 1-COMPARTMENT MODEL",
                        gsub("&", cov_names, "$INPUT ID TIME EVID AMT RATE DV &"),
                        "$DATA DATA_SIM.CSV IGNORE=@",
                        "$SUBROUTINES ADVAN1 TRANS2",
                        "$PK",
                        "CLP=THETA(1)",
                        "VP=THETA(2)",
                        "CLT=CLP",
                        "VT=VP",
                        "CL=CLT*EXP(ETA(1))",
                        "V=VT*EXP(ETA(2))",
                        "S1=V",
                        "$ERROR",
                        "IPRED=F",
                        "Y=F*(1+ERR(1))+ERR(2)",
                        "$THETA",
                        "(1E-10, 0.1)",
                        "(1E-10, 1)",
                        "$OMEGA BLOCK(2)",
                        "0.1",
                        "0.025 0.1",
                        "$SIGMA",
                        "0.1",
                        "0 FIXED",
                        "$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5 NOABORT",
                        "$COV UNCONDITIONAL MATRIX=R",
                        gsub("&", cov_names, "$TABLE ID TIME EVID AMT RATE DV & IPRED PRED CWRES"),
                        "FILE=RESULT.CSV NOAPPEND NOPRINT")
  writeLines(paste(base_model_ctl, collapse = "\n"), script, sep="")
}

# create the simulation NONMEM script with the covariate
create_simCovNM <- function(coef, cov_name, cov_median, sim_index, 
                            cl_v_corr, omegasq, sigmasq){
  cov_names <- paste(colnames(read.csv('PAT_SIM_DATA.CSV'))[-1], collapse = " ")
  data_sim_ctl <- c(
    "$PROBLEM SIMULATION DATA",
    gsub("&", cov_names, "$INPUT ID TIME EVID AMT RATE DV &"),
    "$DATA FULL_SIM_DATA.CSV IGNORE=@",
    "$SUBROUTINES  ADVAN1 TRANS2",
    "$PK",
    "CLP=THETA(1)",
    "VP=THETA(2)",
    "CLSEX=EXP(SEX*THETA(3))",
    "CL=CLP*CLSEX*EXP(ETA(1))",
    "V=VP*EXP(ETA(2))",
    "S1=V",
    "$ERROR",
    "IPRED=F",
    "Y=F*(1+ERR(1))+ERR(2)",
    "$THETA",
    "(0.1, FIXED)",
    "(1, FIXED)",
    gsub("&", coef, "(&, FIXED)"),
    "$OMEGA BLOCK(2)",
    gsub("&", omegasq, "&"),
    gsub("&1", cl_v_corr^2*omegasq, gsub("&2", omegasq, "&1 &2 FIXED")),
    "$SIGMA",
    gsub("&", sigmasq, "& FIXED"),
    "0 FIXED",
    gsub("&", sim_index, "$SIMULATION (123456&) ONLYSIMULATION SUBPROBLEM=1"),
    gsub("&", cov_names, "$TABLE ID TIME EVID AMT RATE DV & FILE=DATA_SIM.TAB NOAPPEND NOPRINT ONEHEADER"))
  if(cov_name != "SEX"){
    line1 <- gsub("&2", cov_median, gsub("&1", cov_name, "CL&1=(&1/&2)**THETA(3)"))
    line2 <- gsub("&1", cov_name, "CL=CLP*CL&1*EXP(ETA(1))")
    data_sim_ctl[startsWith(data_sim_ctl, "CLSEX")] = line1
    data_sim_ctl[startsWith(data_sim_ctl, "CL=")] = line2
  }
  writeLines(paste(data_sim_ctl, collapse = "\n"), "data_sim.ctl", sep="")
}

# create the score test NONMEM script
create_scoreNM <- function(cov_name, cov_median, base_result_file, script){
  base_result <- read.table(base_result_file,header = T, sep = "",
                            as.is= T,skip= 1)
  estimates <- base_result[base_result$ITERATION == -1000000000, ][-c(1, ncol(base_result))]
  cov_names <- paste(colnames(read.csv('PAT_SIM_DATA.CSV'))[-1], collapse = " ")
  score_model_ctl <- c(
    "$PROBLEM SIMULATION DATA",
    gsub("&", cov_names, "$INPUT ID TIME EVID AMT RATE DV &"),
    "$DATA DATA_SIM.CSV IGNORE=@",
    "$SUBROUTINES  ADVAN1 TRANS2",
    "$PK",
    "CLP=THETA(1)",
    "VP=THETA(2)",
    "CLSEX=EXP(SEX*(THETA(3)-1.0))",
    "CL=CLP*CLSEX*EXP(ETA(1))",
    "V=VP*EXP(ETA(2))",
    "S1=V",
    "$ERROR",
    "Y=F*(1+ERR(1))+ERR(2)",
    "$THETA",
    "THETA1",
    "THETA2",
    "1.0",
    "$OMEGA BLOCK(2)",
    "OMEGA.1.1.",
    "OMEGA.2.1. OMEGA.2.2.",
    "$SIGMA",
    "SIGMA.1.1.",
    "0 FIXED",
    "$ESTIMATION METHOD=1 INTERACTION MAXEVAL=1 PRINT=5 NOABORT SIGL=10",
    "$COV UNCONDITIONAL MATRIX=R",
    gsub("&", cov_names, "$TABLE ID TIME EVID AMT RATE DV AGE SEX FILE=COV_MODEL.TAB NOAPPEND NOPRINT ONEHEADER"))
  if(cov_name != "SEX"){
    line1 <- gsub("&2", cov_median, gsub("&1", cov_name, "CL&1=(&1/&2)**(THETA(3)-1.0)"))
    line2 <- gsub("&1", cov_name, "CL=CLP*CL&1*EXP(ETA(1))")
    score_model_ctl[startsWith(score_model_ctl, "CLSEX")] = line1
    score_model_ctl[startsWith(score_model_ctl, "CL=")] = line2
  }
  for(para_name in colnames(estimates)){
    if(estimates[para_name] != 0){
      score_model_ctl <- gsub(para_name, estimates[para_name], score_model_ctl)
    }
  }
  writeLines(paste(score_model_ctl, collapse = "\n"), script, sep="")  
}

# create the covariate model NONMEM script for LRT and Wald's test
create_covNM <- function(cov_name, cov_median, script){
  cov_names <- paste(colnames(read.csv('PAT_SIM_DATA.CSV'))[-1], collapse = " ")
  cov_model_ctl <- c(
    "$PROBLEM SIMULATION DATA",
    gsub("&", cov_names, "$INPUT ID TIME EVID AMT RATE DV &"),
    "$DATA DATA_SIM.CSV IGNORE=@",
    "$SUBROUTINES  ADVAN1 TRANS2",
    "$PK",
    "CLP=THETA(1)",
    "VP=THETA(2)",
    "CLSEX=EXP(SEX*THETA(3))",
    "CL=CLP*CLSEX*EXP(ETA(1))",
    "V=VP*EXP(ETA(2))",
    "S1=V",
    "$ERROR",
    "Y=F*(1+ERR(1))+ERR(2)",
    "$THETA",
    "(0, 0.1)",
    "(0, 1)",
    "(0, 1)",
    "$OMEGA BLOCK(2)",
    "0.1",
    "0.025 0.1",
    "$SIGMA",
    "0.1",
    "0 FIXED",
    "$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5 NOABORT",
    "$COV UNCONDITIONAL MATRIX=R",
    gsub("&", cov_names, "$TABLE ID TIME EVID AMT RATE DV AGE SEX FILE=COV_MODEL.TAB NOAPPEND NOPRINT ONEHEADER"))
  if(cov_name != "SEX"){
    line1 <- gsub("&2", cov_median, gsub("&1", cov_name, "CL&1=(&1/&2)**THETA(3)"))
    line2 <- gsub("&1", cov_name, "CL=CLP*CL&1*EXP(ETA(1))")
    cov_model_ctl[startsWith(cov_model_ctl, "CLSEX")] = line1
    cov_model_ctl[startsWith(cov_model_ctl, "CL=")] = line2
  }
  writeLines(paste(cov_model_ctl, collapse = "\n"), script, sep="")
}

# read the object function value(OFV) 
read_obj <- function(result_file){
  result <- read.table(result_file, header = T, sep = "",
                            as.is= T,skip= 1)  
  return(result[result$ITERATION == -1000000000, "OBJ"])
}

# run NONMEM script
run_NM <- function(NM_script){
  NM_result <- gsub("ctl", "out", NM_script)
  cmd_line <- paste("nmfe74.bat", NM_script, NM_result)
  shell(cmd_line)
}

# summarize the output from NONMEM
run_analysis <- function(sample_size, sample_design, coef_index,
                         cl_cov_list, cov_name, num_sim, cl_v_corr,
                         omegasq, sigmasq, home_dir){
  dir_name_0 <-  paste("Sim", sample_size, sample_design,
                       cov_name, "cl_v_corr", cl_v_corr, 
                       "omegasq", omegasq, "sigmasq", 
                       sigmasq, coef_index, sep = "_")
  proj_dir_0 <- paste(home_dir, dir_name_0, sep = "/")
  setwd(proj_dir_0)
  base_result_list <- paste0("base_model_", 1:num_sim, ".ext")
  cov_result_list <- paste0("cov_model_", 1:num_sim, ".ext")
  score_result_list <- paste0("score_model_", 1:num_sim, ".ext")
  score_coi_list <- paste0("score_model_", 1:num_sim, ".coi")
  score_stat_list <- wam_stat_list <- lrt_stat_list <- c()
  for (i in 1:num_sim) {
    tryCatch({
      base_obj <- read_obj(base_result_list[i])
      cov_obj <-  read_obj(cov_result_list[i])
      # likelihood ratio statistic
      lrt_stat_list[i] <- base_obj - cov_obj
      # wald's test
      cov_result <- read.table(cov_result_list[i], header = T, sep = "",
                               as.is= T,skip= 1) 
      estimate <- cov_result[cov_result$ITERATION == -1000000000, "THETA3"]
      sd <- cov_result[cov_result$ITERATION == -1000000001, "THETA3"]
      wam_stat_list[i] <- (estimate/sd)^2
      # score's test 
      score_result <- read.table(score_result_list[i], header = T, sep = "",
                                 as.is= T,skip= 1)
      score_coi <- read.table(score_coi_list[i], header = T, sep = "",
                              as.is= T,skip= 1)
      score <- score_result[score_result$ITERATION == -1000000008, "THETA3"]
      fim  <- score_coi[score_coi$NAME == "THETA3", "THETA3"]
      score_stat_list[i] <- score^2/fim
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
  final_result <- data.frame(sample_size, sample_design, 
                             coef=cl_cov_list[coef_index],
                             cov_name, num_sim, cl_v_corr,
                             omegasq, sigmasq,
                             lrt=lrt_stat_list,
                             wald=wam_stat_list,
                             score=score_stat_list)
  # create a result directory under the home directory
  setwd(home_dir)
  if(!dir.exists("Results")){
    dir.create("Results")
  }
  setwd(paste(home_dir, "Results", sep="/"))
  write.csv(final_result, paste(dir_name_0, ".csv", sep=""), row.names = F)  
}


# summarize the output from NONMEM
run_analysis_corr <- function(sample_size, sample_design, coef_index,
                              cl_cov_list, cov_corr, num_sim, cl_v_corr, 
                              omegasq, sigmasq, home_dir){
  dir_name_0 <-  paste("Sim", sample_size, sample_design,
                       paste(cov_corr, collapse = "_"), 
                       "cl_v_corr", cl_v_corr, 
                       "omegasq", omegasq, 
                       "sigmasq", sigmasq, 
                       coef_index, sep = "_")
  proj_dir_0 <- paste(home_dir, dir_name_0, sep = "/")
  setwd(proj_dir_0)
  base_result_list <- paste0("base_model_", 1:num_sim, ".ext")
  cov_result_list <- paste0("cov_model_", 1:num_sim, ".ext")
  score_result_list <- paste0("score_model_", 1:num_sim, ".ext")
  score_coi_list <- paste0("score_model_", 1:num_sim, ".coi")
  score_stat_list <- wam_stat_list <- lrt_stat_list <- c()
  cov_result_corr_list <- paste0("cov_model_corr", 1:num_sim, ".ext")
  score_result_corr_list <- paste0("score_model_corr", 1:num_sim, ".ext")
  score_coi_corr_list <- paste0("score_model_corr", 1:num_sim, ".coi")
  score_stat_corr_list <- wam_stat_corr_list <- lrt_stat_corr_list <- c()
  for (i in 1:num_sim) {
    tryCatch({
      base_obj <- read_obj(base_result_list[i])
      cov_obj <-  read_obj(cov_result_list[i])
      # likelihood ratio statistic
      lrt_stat_list[i] <- base_obj - cov_obj
      # wald's test
      cov_result <- read.table(cov_result_list[i], header = T, sep = "",
                               as.is= T,skip= 1) 
      estimate <- cov_result[cov_result$ITERATION == -1000000000, "THETA3"]
      sd <- cov_result[cov_result$ITERATION == -1000000001, "THETA3"]
      wam_stat_list[i] <- ifelse(length(sd)>0, (estimate/sd)^2, 0)
      # score's test 
      score_result <- read.table(score_result_list[i], header = T, sep = "",
                                 as.is= T,skip= 1)
      score_coi <- read.table(score_coi_list[i], header = T, sep = "",
                              as.is= T,skip= 1)
      score <- score_result[score_result$ITERATION == -1000000008, "THETA3"]
      fim  <- score_coi[score_coi$NAME == "THETA3", "THETA3"]
      score_stat_list[i] <- ifelse(length(fim)>0, score^2/fim, 0)
      
      cov_corr_obj <-  read_obj(cov_result_corr_list[i])
      # likelihood ratio statistic
      lrt_stat_corr_list[i] <- base_obj - cov_corr_obj
      # wald's test
      cov_corr_result <- read.table(cov_result_corr_list[i], header = T, sep = "",
                                    as.is= T,skip= 1) 
      estimate <- cov_corr_result[cov_corr_result$ITERATION == -1000000000, "THETA3"]
      sd <- cov_corr_result[cov_corr_result$ITERATION == -1000000001, "THETA3"]
      wam_stat_corr_list[i] <- ifelse(length(sd)>0, (estimate/sd)^2, 0)
      # score's test 
      score_corr_result <- read.table(score_result_corr_list[i], header = T, sep = "",
                                      as.is= T,skip= 1)
      score_coi <- read.table(score_coi_corr_list[i], header = T, sep = "",
                              as.is= T,skip= 1)
      score <- score_corr_result[score_corr_result$ITERATION == -1000000008, "THETA3"]
      fim  <- score_coi[score_coi$NAME == "THETA3", "THETA3"]
      score_stat_corr_list[i] <- ifelse(length(fim)>0, score^2/fim, 0)
    }, error=function(e){cat("ERROR : sim #",conditionMessage(e), "\n")})
  }
  final_result_true <- data.frame(sample_size, sample_design, 
                                  coef=cl_cov_list[coef_index],
                                  cov_name=cov_corr[1], num_sim, cl_v_corr,
                                  omegasq, sigmasq,
                                  lrt=lrt_stat_list,
                                  wald=wam_stat_list,
                                  score=score_stat_list)
  final_result_false <- data.frame(sample_size, sample_design, 
                                   coef=cl_cov_list[coef_index],
                                   cov_name=cov_corr[2], num_sim, cl_v_corr,
                                   omegasq, sigmasq,
                                   lrt=lrt_stat_corr_list,
                                   wald=wam_stat_corr_list,
                                   score=score_stat_corr_list)
  # create a result directory under the home directory
  setwd(home_dir)
  if(!dir.exists("Results")){
    dir.create("Results")
  }
  setwd(paste(home_dir, "Results", sep="/"))
  write.csv(rbind(final_result_true, final_result_false), 
            paste(dir_name_0, ".csv", sep=""), row.names = F)  
}

# the main simulation function
run_simulation <- function(sample_size, sample_design, coef_index,
                           cl_cov_list, cov_name, num_sim, cl_v_corr,
                           omegasq, sigmasq, home_dir){
  # set the home directory
  setwd(home_dir)
  coef <- cl_cov_list[coef_index]

  # create a project directory for each situation
  dir_name_0 <-  paste("Sim", sample_size, sample_design,
                       cov_name, "cl_v_corr", cl_v_corr, 
                       "omegasq", omegasq, "sigmasq", 
                       sigmasq, coef_index, sep = "_")
  if(!dir.exists(dir_name_0)){
    dir.create(dir_name_0)
  }
  proj_dir_0 <- paste(home_dir, "/", dir_name_0, sep = "")
  files <- c("nmfe74.bat", "PAT_SIM_DATA.CSV")
  file.copy(from=files, to=proj_dir_0,
            overwrite = TRUE, recursive = FALSE,
            copy.mode = TRUE)
  setwd(proj_dir_0)
  base_model_list <- paste0("base_model_", 1:num_sim, ".ctl")
  base_result_list <- paste0("base_model_", 1:num_sim, ".ext")
  cov_model_list <- paste0("cov_model_", 1:num_sim, ".ctl")
  score_model_list <- paste0("score_model_", 1:num_sim, ".ctl")
  i <- j <- 1
  while(i < (1+num_sim)){
    status <- tryCatch({
      setwd(proj_dir_0)
      # bootstrap the covariate info from NHANES database
      create_data_template(sample_design, sample_size, cov_name)
      # read the median of the covariate
      cov_median <- read.table("COV_MEDIAN.TXT")$V1
      # create NONMEM script for simulation
      create_simCovNM(coef, cov_name, cov_median, i, 
                      cl_v_corr, omegasq, sigmasq)
      # run the sim script
      run_NM("data_sim.ctl")
      # clean the sim output from NONMEM
      sim_df <- read.table(file="DATA_SIM.TAB",
                           header = T, sep = "",
                           as.is= T,skip= 1)
      # delete the negative simulated observation
      sim_df_1 <- sim_df %>% dplyr::filter(DV>=0) 
      write.table(sim_df_1, "DATA_SIM.CSV",  sep=",", row.names = F, col.names = F)
      
      # create the base model and run it
      create_baseNM(base_model_list[i])
      run_NM(base_model_list[i])
      # create score test script and run it
      create_scoreNM(cov_name, cov_median, base_result_list[i], score_model_list[i])
      run_NM(score_model_list[i])
      # create  LRT and Wald's test script and run it
      create_covNM(cov_name, cov_median, cov_model_list[i])
      run_NM(cov_model_list[i])},
      error = function(e) "error")
    if (status == "error"){
      i <- i 
    } else{
      i <- i + 1
    }
    j <- j + 1
  }
}

# the main simulation function for correlated covariate
# cl_cov_list should not contain 0
run_simulation_corr <- function(sample_size, sample_design, coef_index,
                                cl_cov_list, cov_corr, num_sim, cl_v_corr,
                                omegasq, sigmasq, home_dir){
  # set the home directory
  setwd(home_dir)
  coef <- cl_cov_list[coef_index]
  cov_true <- cov_corr[1]
  cov_false <- cov_corr[2]
  
  # create a project directory for each situation
  dir_name_0 <-  paste("Sim", sample_size, sample_design,
                       paste(cov_corr, collapse = "_"), 
                       "cl_v_corr", cl_v_corr, 
                       "omegasq", omegasq, 
                       "sigmasq", sigmasq, 
                       coef_index, sep = "_")
  if(!dir.exists(dir_name_0)){
    dir.create(dir_name_0)
  }
  proj_dir_0 <- paste(home_dir, "/", dir_name_0, sep = "")
  files <- c("nmfe74.bat", "PAT_SIM_DATA.CSV")
  file.copy(from=files, to=proj_dir_0,
            overwrite = TRUE, recursive = FALSE,
            copy.mode = TRUE)
  setwd(proj_dir_0)
  base_model_list <- paste0("base_model_", 1:num_sim, ".ctl")
  base_result_list <- paste0("base_model_", 1:num_sim, ".ext")
  cov_model_list <- paste0("cov_model_", 1:num_sim, ".ctl")
  score_model_list <- paste0("score_model_", 1:num_sim, ".ctl")
  cov_model_corr_list <- paste0("cov_model_corr", 1:num_sim, ".ctl")
  score_model_corr_list <- paste0("score_model_corr", 1:num_sim, ".ctl")
  i <- 1
  j <- 1
  while(i < (1+num_sim)){
    indication <- tryCatch({
      setwd(proj_dir_0)
      # bootstrap the covariate info from NHANES database
      create_data_template(sample_design, sample_size, cov_true)
      # read the median of the covariate
      cov_true_median <- read.table("COV_MEDIAN.TXT")$V1
      # create NONMEM script for simulation
      create_simCovNM(coef, cov_true, cov_true_median, j, 
                      cl_v_corr, omegasq, sigmasq)
      # run the sim script
      run_NM("data_sim.ctl")
      # clean the sim output from NONMEM
      sim_df <- read.table(file="DATA_SIM.TAB",
                           header = T, sep = "",
                           as.is= T,skip= 1)
      # delete the negative simulated observation
      sim_df_1 <- sim_df %>% dplyr::filter(DV>=0) 
      write.table(sim_df_1, "DATA_SIM.CSV",  sep=",", row.names = F, col.names = F)
      # cov_true median
      cov_false_median <- median(sim_df %>% 
                                   group_by(ID) %>% 
                                   filter(row_number() == 1) %>%
                                   ungroup(ID) %>%
                                   select(cov_false) %>%
                                   unlist())
      # create the base model and run it
      create_baseNM(base_model_list[i])
      run_NM(base_model_list[i])
      # create score test script and run it
      create_scoreNM(cov_true, cov_true_median, base_result_list[i], score_model_list[i])
      run_NM(score_model_list[i])
      create_scoreNM(cov_false, cov_false_median, base_result_list[i], score_model_corr_list[i])
      run_NM(score_model_corr_list[i])
      # create  LRT and Wald's test script and run it
      create_covNM(cov_true, cov_true_median, cov_model_list[i])
      run_NM(cov_model_list[i])
      create_covNM(cov_false, cov_false_median, cov_model_corr_list[i])
      run_NM(cov_model_corr_list[i])},
      error = function(e) "error")
    if (indication == "error"){
      i <- i 
    } else{
      i <- i + 1
    }
    j <- j + 1
  }
}
