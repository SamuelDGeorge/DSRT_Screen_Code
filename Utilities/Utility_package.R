if (!require("tidyr")) {
  install.packages("tidyr", dependencies = TRUE)
  library(tidyr)
}

if (!require("stringr")) {
  install.packages("stringr", dependencies = TRUE)
  library(stringr)
}

require(parallel)

## This function is described in "CrisprFunctionScripts.R"
collect_file_vector <- function(FolderPath) {
  files <- list.files(FolderPath) 
  dirs <- list.dirs(FolderPath, full.names = FALSE)
  dirs = dirs[-1]
  list_items = setdiff(files,dirs)
  return(list_items)
}

## A pipeline Function. Takes all files within a folder, performs the function
## and puts output in output_location.
## only supports functions that take filename as input and output example
## function(fileInputName,fileOutputName)
## additional arguments can be added as a vector (***coerced to a string when given as arguments!)

run_parallel <- function(Function_to_perform, input_location, additional_args = c(), cores = 8, environment = lsf.str()) {
  files <- collect_file_vector(input_location)
  inputs <- c()
  outputs <- c()
  for (i in 1:length(files)) {
    InputFile = paste(input_location,"/",sep = "")
    InputFile = paste(InputFile,files[i],sep = "")
    inputs <- c(inputs,InputFile)
    OutputName = str_split(files[i],pattern = "[.]")
    OutputName = OutputName[[1]][1]
    outputs <- c(outputs,OutputName)
  }
  
  core_count = detectCores()
  cores_to_use = 1
  if (core_count < cores){
    cores_to_use = core_count
  } else {
    cores_to_use = cores
  }
  print(cat("Running task on", cores_to_use, "cores\n", sep = " "))
  if(.Platform$OS.type == "windows") {
    c1 = makeCluster(cores_to_use)
    clusterExport(c1,varlist = as.vector(environment))
    result = clusterMap(c1,Function_to_perform, as.list(inputs), as.list(outputs), MoreArgs = as.list(additional_args))
    stopCluster(c1)
  } else {
    result = mcmapply(Function_to_perform,as.list(inputs),as.list(outputs), MoreArgs = as.list(additional_args), mc.cores = cores_to_use)
  }  

  return(result)
}

run_sequential <- function(Function_to_perform, input_location, additional_args = c()) {
  files <- collect_file_vector(input_location)
  result = c()
  for (i in 1:length(files)) {
    InputFile = paste(input_location,"/",sep = "")
    InputFile = paste(InputFile,files[i],sep = "")
    OutputName = str_split(files[i],pattern = "[.]")
    OutputName = OutputName[[1]][1]
    args = c(InputFile,OutputName)
    args = c(args,additional_args)
    result = c(result, do.call(Function_to_perform,as.list(args)))
  }
}

pipeline <- function(Function_to_perform, input_location, output_location, additional_args = c(), parallel_run = FALSE, environment = lsf.str(), cores = 8) {
  wd <- getwd()
  setwd(output_location)
  result = c()
  if (parallel_run) {
    result = run_parallel(Function_to_perform, input_location, additional_args, cores, environment)
  } else {
    result = run_sequential(Function_to_perform, input_location, additional_args)
  }
  setwd(wd) ## set the wd back to the starting wd, had to change to the output_location
  return(result)
}

pipeline_with_finish_function <- function(Function_to_perform, finish_function,input_location, output_location, func_one_additional_args = c(), func_two_args = c(), parallel_run = FALSE, environment = lsf.str()) {
  wd <- getwd()
  setwd(output_location)
  result = c()
  if (parallel_run) {
    result = run_parallel(Function_to_perform, input_location, func_one_additional_args, environment)
  } else {
    result = run_sequential(Function_to_perform, input_location, func_one_additional_args)
  }
  
  result = do.call(finish_function,as.list(func_two_args))
  setwd(wd) ## set the wd back to the starting wd, had to change to the output_location
  return(result)
}



