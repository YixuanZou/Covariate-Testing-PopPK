#-------------------------------------------------------------
# Title: R functions for simulation study
# Author: Yixuan Zou
# Create Date: 5/26/2020
# Last Update: 6/2/2020
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
create_data_template <- function(sample_design, sample_size){
  col_names <- c('CID', 'TIME', 'EVID', 'AMT',
                 'RATE', 'DV', 'AGE', 'SEX')
  ID <- c()
  TIME <- c()
  EVID <- c()
  AMT <- c()
  RATE <- c()
  DV <- c()
  AGE <- c()
  WT <- c()
  SEX <- c()
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
      AGE <- NA
      SEX <- NA
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
      AGE <- NA
      SEX <- NA
    }
  }
  data <- data.frame(ID, TIME, EVID, AMT, RATE, DV, AGE, SEX)
  pat_data <- read.csv('pat_sim_data.csv')
  for (i in 1:sample_size) {
    id_index <- which(data$ID == i)
    for(j in id_index){
      data[j, c('AGE', 'SEX')] <- pat_data[i, c('AGE', 'SEX')]
    }
  }
  colnames(data) <- col_names
  write.table(data, 'FULL_SIM_DATA.CSV', sep = ',', row.names = FALSE, col.names = FALSE)
}

# create the NONMEM script of the base model
create_baseNM <- function(script){
    base_model_ctl <- c("$PROBLEM 1-COMPARTMENT MODEL",
                        "$INPUT ID TIME EVID AMT RATE DV AGE SEX",
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
                        "$OMEGA",
                        "0.1",
                        "0.1",
                        "$SIGMA",
                        "0.1",
                        "0 FIXED",
                        "$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5 NOABORT",
                        "$COV UNCONDITIONAL MATRIX=R",
                        "$TABLE ID TIME EVID AMT RATE DV AGE SEX IPRED PRED CWRES",
                        "FILE=RESULT.CSV NOAPPEND NOPRINT")
  writeLines(paste(base_model_ctl, collapse = "\n"),
             script, sep="")
}

# create the simulation NONMEM script with the covariate
create_simCovNM <- function(coef, cov_name, sim_index){
  data_sim_ctl <- c("$PROBLEM SIMULATION DATA",
                    "$INPUT ID TIME EVID AMT RATE DV AGE SEX",
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
                    "$OMEGA",
                    "0.1 FIXED",
                    "0.1 FIXED",
                    "$SIGMA",
                    "0.1 FIXED",
                    "0 FIXED",
                    gsub("&", sim_index, "$SIMULATION (123456&) ONLYSIMULATION SUBPROBLEM=1"),
                    "$TABLE ID TIME EVID AMT RATE DV AGE SEX FILE=DATA_SIM.TAB NOAPPEND NOPRINT ONEHEADER")
  if(cov_name == "AGE"){
    data_sim_ctl[startsWith(data_sim_ctl, "CLSEX")] = "CLAGE=(AGE/52)**THETA(3)"
    data_sim_ctl[startsWith(data_sim_ctl, "CL=")] = "CL=CLP*CLAGE*EXP(ETA(1))"
  }
  writeLines(paste(data_sim_ctl, collapse = "\n"),
             "data_sim.ctl", sep="")
}

# create the simulation NONMEM script
create_simNM <- function(sim_index){
  data_sim_ctl <- c("$PROBLEM SIMULATION DATA",
                    "$INPUT ID TIME EVID AMT RATE DV AGE SEX",
                    "$DATA FULL_SIM_DATA.CSV IGNORE=@",
                    "$SUBROUTINES  ADVAN1 TRANS2",
                    "$PK",
                    "CLP=THETA(1)",
                    "VP=THETA(2)",
                    "CL=CLP*EXP(ETA(1))",
                    "V=VP*EXP(ETA(2))",
                    "S1=V",
                    "$ERROR",
                    "IPRED=F",
                    "Y=F*(1+ERR(1))+ERR(2)",
                    "$THETA",
                    "(0.1, FIXED)",
                    "(1, FIXED)",
                    "$OMEGA",
                    "0.1 FIXED",
                    "0.1 FIXED",
                    "$SIGMA",
                    "0.1 FIXED",
                    "0 FIXED",
                    gsub("&", sim_index, "$SIMULATION (123456&) ONLYSIMULATION SUBPROBLEM=1"),
                    "$TABLE ID TIME EVID AMT RATE DV AGE SEX FILE=DATA_SIM.TAB NOAPPEND NOPRINT ONEHEADER")
  writeLines(paste(data_sim_ctl, collapse = "\n"),
             "data_sim.ctl", sep="")

}

# create the score test NONMEM script
create_scoreNM <- function(cov_name, base_result_file, script){
  base_result <- read.table(base_result_file,header = T, sep = "",
                            as.is= T,skip= 1)
  estimates <- base_result[base_result$ITERATION == -1000000000, ][-c(1, ncol(base_result))]
  score_model_ctl <- c(
    "$PROBLEM SIMULATION DATA",
    "$INPUT ID TIME EVID AMT RATE DV AGE SEX",
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
    "$OMEGA",
    "OMEGA.1.1.",
    "OMEGA.2.2.",
    "$SIGMA",
    "SIGMA.1.1.",
    "0 FIXED",
    "$ESTIMATION METHOD=1 INTERACTION MAXEVAL=1 PRINT=5 NOABORT SIGL=10",
    "$COV UNCONDITIONAL MATRIX=R",
    "$TABLE ID TIME EVID AMT RATE DV AGE SEX FILE=COV_MODEL.TAB NOAPPEND NOPRINT ONEHEADER")
  if(cov_name == "AGE"){
    score_model_ctl[startsWith(score_model_ctl, "CLSEX")] = "CLAGE=(AGE/52)**(THETA(3)-1.0)"
    score_model_ctl[startsWith(score_model_ctl, "CL=")] = "CL=CLP*CLAGE*EXP(ETA(1))"
  }
  for(para_name in colnames(estimates)){
    if(estimates[para_name] > 0){
      score_model_ctl <- gsub(para_name, estimates[para_name], score_model_ctl)
    }
  }
  writeLines(paste(score_model_ctl, collapse = "\n"),
             script, sep="")  
}

# create the covariate model NONMEM script for LRT and Wald's test
create_covNM <- function(cov_name, script){
  cov_model_ctl <- c(
    "$PROBLEM SIMULATION DATA",
      "$INPUT ID TIME EVID AMT RATE DV AGE SEX",
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
      "$OMEGA",
      "0.1",
      "0.1",
      "$SIGMA",
      "0.1",
      "0 FIXED",
      "$ESTIMATION METHOD=1 INTERACTION MAXEVAL=9999 PRINT=5 NOABORT",
      "$COV UNCONDITIONAL MATRIX=R",
      "$TABLE ID TIME EVID AMT RATE DV AGE SEX FILE=COV_MODEL.TAB NOAPPEND NOPRINT ONEHEADER")
  if(cov_name == "AGE"){
    cov_model_ctl[startsWith(cov_model_ctl, "CLSEX")] = "CLAGE=(AGE/52)**THETA(3)"
    cov_model_ctl[startsWith(cov_model_ctl, "CL=")] = "CL=CLP*CLAGE*EXP(ETA(1))"
  }
  writeLines(paste(cov_model_ctl, collapse = "\n"),
             script, sep="")
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
                         cl_cov_list, cov_name, num_sim, home_dir,
                         analysis){
  dir_name_0 <- paste("Sim", sample_size, sample_design,
                      cov_name, coef_index, sep = "_")
  proj_dir_0 <- paste(home_dir, dir_name_0, analysis, sep = "/")
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
                             cov_name, num_sim, analysis,
                             lrt=lrt_stat_list,
                             wald=wam_stat_list,
                             score=score_stat_list)
  # create a result directory under the home directory
  setwd(home_dir)
  if(!dir.exists("Results")){
    dir.create("Results")
  }
  setwd(paste(home_dir, "Results", sep="/"))
  write.csv(final_result, paste(dir_name_0, "_", analysis, ".csv", sep=""), row.names = F)  
  # create a plot directory under the home directory
  setwd(home_dir)
  if(!dir.exists("Plots")){
    dir.create("Plots")
  }
  setwd(paste(home_dir, "Plots", sep="/"))
  df_long <- final_result %>% select(lrt, wald, score) %>%
    pivot_longer(c(lrt, wald, score), names_to = "Test", values_to = "Statistic")
  p <- ggplot(df_long, aes(Statistic, color = Test)) + geom_density() + 
    xlim(0, 25) + 
    theme_few() + scale_colour_few("Light")
  ggsave(paste(dir_name_0, "_", analysis, ".jpg", sep=""), p, dpi = 300)
}

# the main simulation function
run_simulation <- function(sample_size, sample_design, coef_index,
                           cl_cov_list, cov_name, num_sim, home_dir){
  # set the home directory
  setwd(home_dir)
  coef <- cl_cov_list[coef_index]

  # create a project directory for each situation
  dir_name_0 <- paste("Sim", sample_size, sample_design,
                      cov_name, coef_index, sep = "_")
  if(!dir.exists(dir_name_0)){
    dir.create(dir_name_0)
  }
  proj_dir_0 <- paste(home_dir, "/", dir_name_0, sep = "")
  files <- c("nmfe74.bat", "pat_sim_data.csv")
  file.copy(from=files, to=proj_dir_0,
            overwrite = TRUE, recursive = FALSE,
            copy.mode = TRUE)
  setwd(proj_dir_0)
  type_analysis <- c('type_one', 'power')
  base_model_list <- paste0("base_model_", 1:num_sim, ".ctl")
  base_result_list <- paste0("base_model_", 1:num_sim, ".ext")
  cov_model_list <- paste0("cov_model_", 1:num_sim, ".ctl")
  score_model_list <- paste0("score_model_", 1:num_sim, ".ctl")
  for(i in 1:num_sim){
    setwd(proj_dir_0)
    # bootstrap the covariate info from NHANES database
    create_data_template(sample_design, sample_size)
    # create directory for type 1 and power analysis
    for(analysis in type_analysis){
      setwd(proj_dir_0)
      dir_name_1 <- analysis
      if(!dir.exists(dir_name_1)){
        dir.create(dir_name_1)
      }
      proj_dir_1 <- paste(proj_dir_0, "/", dir_name_1, sep = "")
      files <- c("nmfe74.bat", "FULL_SIM_DATA.CSV")
      file.copy(from=files, to=proj_dir_1,
                overwrite = TRUE, recursive = FALSE,
                copy.mode = TRUE)
      setwd(proj_dir_1)
      # create NONMEM script for simulation
      if(analysis == 'type_one'){
        create_simNM(i)
      } else{
        create_simCovNM(coef, cov_name, i)
      }
      # run the sim script
      run_NM("data_sim.ctl")
      # clean the sim output from NONMEM
      sim_df <- read.table(file="DATA_SIM.TAB",
                           header = T, sep = "",
                           as.is= T,skip= 1)
      write.table(sim_df, "DATA_SIM.CSV",  sep=",", row.names = F, col.names = F)
      
      # create the base model and run it
      create_baseNM(base_model_list[i])
      run_NM(base_model_list[i])
      # create score test script and run it
      create_scoreNM(cov_name, base_result_list[i], score_model_list[i])
      run_NM(score_model_list[i])
      # create  LRT and Wald's test script and run it
      create_covNM(cov_name, cov_model_list[i])
      run_NM(cov_model_list[i])
    }
  }
}
