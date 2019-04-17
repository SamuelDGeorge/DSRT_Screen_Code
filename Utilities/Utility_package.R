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
collect_file_vector <- function(FolderPath, file_type) {
  files <- list_files_with_exts(FolderPath, file_type, full.names = FALSE)
  return(files)
}

## A pipeline Function. Takes all files within a folder, performs the function
## and puts output in output_location.
## only supports functions that take filename as input and output example
## function(fileInputName,fileOutputName)
## additional arguments can be added as a vector (***coerced to a string when given as arguments!)

run_parallel <- function(Function_to_perform, input_location, additional_args = c(), file_types = "xlsx",cores = 8, environment = lsf.str()) {
  files <- collect_file_vector(input_location, file_types)
  if(length(files) == 0){return()}
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

run_sequential <- function(Function_to_perform, input_location, additional_args = c(), file_types = "xlsx") {
  files <- collect_file_vector(input_location, file_types)
  if(length(files) == 0){return()}
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

pipeline <- function(Function_to_perform, input_location, output_location, additional_args = c(), file_types = "xlsx",parallel_run = FALSE, cores = 8, environment = lsf.str()) {
  wd <- getwd()
  setwd(output_location)
  result = c()
  if (parallel_run) {
    result = run_parallel(Function_to_perform, input_location, additional_args, file_types, cores, environment)
  } else {
    result = run_sequential(Function_to_perform, input_location, additional_args, file_types)
  }
  setwd(wd) ## set the wd back to the starting wd, had to change to the output_location
  return(result)
}

pipeline_with_finish_function <- function(Function_to_perform, finish_function,input_location, output_location, func_one_additional_args = c(), func_two_args = c(), file_types = "xlsx", parallel_run = FALSE, cores = 8, environment = lsf.str()) {
  wd <- getwd()
  setwd(output_location)
  result = c()
  if (parallel_run) {
    result = run_parallel(Function_to_perform, input_location, func_one_additional_args, file_types, cores, environment)
  } else {
    result = run_sequential(Function_to_perform, input_location, func_one_additional_args, file_types)
  }
  
  result = do.call(finish_function,as.list(func_two_args))
  setwd(wd) ## set the wd back to the starting wd, had to change to the output_location
  return(result)
}

pipeline_folder_recursive <- function(Function_to_perform, highest_level_input, highest_level_output, additional_args = c(), file_types = ".xlsx", parallel_run = FALSE, cores = 8, environment = lsf.str()){
  
  input_folders <- list.dirs(highest_level_input, full.names = FALSE)
  for (i in input_folders){
    #Have the input directory
    input_directory = paste(highest_level_input,i,sep = "/")
    
    #Make output directory if it doesn't exist
    output_directory = paste(highest_level_output,i,sep="/")
    dir.create(output_directory, showWarnings = FALSE)
    pipeline(Function_to_perform, input_directory, output_directory, additional_args, file_types, parallel_run, cores, environment)
    
  }
  
}

run_parallel_combine <- function(Function_to_perform, input_location, root_directory, additional_args = c(), file_types = "xlsx",cores = 8, environment = lsf.str()) {
  files <- collect_file_vector(input_location, file_types)
  if(length(files) == 0){return()}
  inputs <- c()
  outputs <- c()
  root_append = str_replace(root_directory, "/", "-")
  for (i in 1:length(files)) {
    InputFile = paste(input_location,"/",sep = "")
    InputFile = paste(InputFile,files[i],sep = "")
    inputs <- c(inputs,InputFile)
    OutputName = str_split(files[i],pattern = "[.]")
    OutputName = OutputName[[1]][1]
    if (!identical(root_append,"")){OutputName = paste(root_append,OutputName,sep = "-")}
    OutputName = str_replace(OutputName, "/", "-")
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

run_sequential_combine <- function(Function_to_perform, input_location, root_directory, additional_args = c(), file_types = "xlsx") {
  files <- collect_file_vector(input_location, file_types)
  if(length(files) == 0){return()}
  result = c()
  root_append = str_replace(root_directory, "/", "-")
  for (i in 1:length(files)) {
    InputFile = paste(input_location,"/",sep = "")
    InputFile = paste(InputFile,files[i],sep = "")
    OutputName = str_split(files[i],pattern = "[.]")
    OutputName = OutputName[[1]][1]
    if (!identical(root_append,"")){OutputName = paste(root_append,OutputName,sep = "-")}
    OutputName = str_replace(OutputName, "/", "-")
    print(OutputName)
    args = c(InputFile,OutputName)
    args = c(args,additional_args)
    result = c(result, do.call(Function_to_perform,as.list(args)))
  }
}

pipeline_combine <- function(Function_to_perform, input_location, root_directory,additional_args = c(), file_types = "xlsx",parallel_run = FALSE, cores = 8, environment = lsf.str()) {

  result = c()
  if (parallel_run) {
    result = run_parallel_combine(Function_to_perform, input_location, root_directory, additional_args, file_types, cores, environment)
  } else {
    result = run_sequential_combine(Function_to_perform, input_location, root_directory,additional_args, file_types)
  }
  return(result)
}

pipeline_folder_combine <- function(Function_to_perform, highest_level_input, highest_level_output, additional_args = c(), file_types = ".xlsx", parallel_run = FALSE, cores = 8, environment = lsf.str()){
  wd <- getwd()
  dir.create(highest_level_output, showWarnings = FALSE)
  setwd(highest_level_output)
  
  input_folders <- list.dirs(highest_level_input, full.names = FALSE)
  for (i in input_folders){
    #Have the input directory
    #print(i)
    input_directory = paste(highest_level_input,i,sep = "/")
    pipeline_combine(Function_to_perform, input_directory, i,additional_args, file_types, parallel_run, cores, environment)
    
  }
  setwd(wd) ## set the wd back to the starting wd, had to change to the output_location
}





