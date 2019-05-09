library(dplyr)
library(tidyr)
library(synergyfinder)
library(nplr)
library(RColorBrewer)
library(ggplot2)
library(magrittr)
library(tools)

# 
# if (!require("dplyr")) {
#   install.packages("dplyr", dependencies = TRUE)
#   library(dplyr)
# }
# 
# if (!require("tidyr")) {
#   install.packages("tidyr", dependencies = TRUE)
#   library(tidyr)
# }
# 
# if (!require("synergyfinder")) {
#   if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#   BiocManager::install("synergyfinder", version = "3.8")
#   library(synergyfinder)
# }
# 
# if (!require("nplr")) {
#   install.packages("nplr", dependencies = TRUE)
#   library(nplr)
# }
# 
# if (!require("RColorBrewer")) {
#   install.packages("RColorBrewer", dependencies = TRUE)
#   library(RColorBrewer)
# }
# 
# if (!require("ggplot2")) {
#   install.packages("ggplot2", dependencies = TRUE)
#   library(ggplot2)
# }

#Functions Needed to Run in Parallel
'%>%' <- magrittr::`%>%`


collect_file_vector <- function(FolderPath, file_type) {
  files <- list_files_with_exts(FolderPath, file_type, full.names = FALSE)
  return(files)
}

collect_file_vector_no_ext <- function(FolderPath, file_type){
  files <- list_files_with_exts(FolderPath, file_type, full.names = FALSE)
  return(files)
}

##The following will get the drug concentrations as a vector, which will define the size of the drug matrix
build_data_frame <- function(input_file,range, Viability_Data = TRUE, Normalize = TRUE, Drug_A = "Drug_A", Drug_B = "Drug_B", Cell_Line = "Cell", Conc_A = "uM", Conc_B = "uM", replicate = 1){
  suppressMessages(x <- readxl::read_xlsx(paste(input_file), range = range, col_names = FALSE))
  xl <- as.matrix(x)
  DrugDoseA = na.omit(as.numeric(as.vector(xl[1,2:ncol(xl)])))
  DrugDoseB = na.omit(as.numeric(as.vector(xl[2:nrow(xl),1])))
  matrix_x = length(DrugDoseA)
  matrix_y = length(DrugDoseB)
  
  
  
  metadata = data.frame(c(Drug_A), c(Drug_B), c(Cell_Line), c(Conc_A), c(Conc_B), c(replicate))
  names(metadata) = c("DrugA","DrugB", "Cell_Line","DrugA.conc","DrugB.conc","Replicate")
  y = "inhibition"
  if (Viability_Data){y = "viability"}
  metadata[1,"type"] <- y
  
  ##First isolate the matrix that corresponds to the Response values, and set that to a vector
  
  Response.mat = as.matrix.data.frame(x[2:(matrix_y+1),2:(matrix_x+1)])
  Response = as.numeric(as.vector(Response.mat))
  
  ##Then define the DrugA and DrugB concentration values. Transform those values into a vector of the same length as the Response vector.
  
  DrugDoseArep = as.vector(rep(DrugDoseA, each = matrix_y))
  DrugDoseBrep = as.vector(rep(DrugDoseB, matrix_x))
  DrugConcA = unique(DrugDoseA)
  DrugConcB = unique(DrugDoseB)
  combo_columns = as.numeric(length(DrugConcA))
  combo_rows = as.numeric(length(DrugConcB))
  
  ##Bind the reponse, DrugDoseArep, and DrugDoseBrep vectors. Make sure they are numeric
  
  df = as.data.frame(cbind(Response, DrugDoseArep, DrugDoseBrep))
  df$DrugDoseArep = as.numeric(df$DrugDoseArep)
  df$DrugDoseBrep = as.numeric(df$DrugDoseBrep)
  df$Response = as.numeric(df$Response)
  
  ##Make a new column that will take in to account any repeat treatments (i.e. extra control cells or if single treatments were run in duplicate)
  
  df <- df %>%
    tidyr::unite(AplusB, DrugDoseArep, DrugDoseBrep, remove = T)
  
  ##Then we group the dataframe by repeat treatments and only cary the mean value forward. This will eliminate duplicated treatments
  
  df = df %>% 
    dplyr::group_by(AplusB) %>%
    dplyr::summarise(Response = mean(Response)) %>%
    tidyr::separate(AplusB, c("ConcCol","ConcRow"), sep = "_")
  df$ConcCol = round(as.numeric(df$ConcCol), digits = 2)
  df$ConcRow = round(as.numeric(df$ConcRow), digits = 2)

  ## Find the control well and normalize
  if (Viability_Data) {
    if(Normalize) {
      ctrldf = df %>%
        dplyr::filter(ConcCol==0,ConcRow==0)
      ctrl=ctrldf$Response[1]
      df$Response = df$Response/ctrl*100
    }
  } else {
    if(Normalize) {
      max_response <- max(df$Response)
      min_response <- min(df$Response)
      df$Response = (df$Response-min_response)/(max_response-min_response)*100
    }
  }
  
  ## Now shape the data to the format needed for synergyFinder
  df = df %>%
    dplyr::mutate(BlockID = 1) %>%
    dplyr::mutate(DrugRow = as.character(metadata$DrugB[1])) %>%
    dplyr::mutate(DrugCol = as.character(metadata$DrugA[1])) %>%
    dplyr::mutate(ConcRowUnit = as.character(metadata$DrugB.conc)) %>%
    dplyr::mutate(ConcColUnit = as.character(metadata$DrugA.conc)) %>%
    as.data.frame(df) %>%
    dplyr::arrange(ConcCol) %>%
    dplyr::arrange(ConcRow) %>%
    dplyr::mutate(Col = rep(1:combo_columns, combo_rows)) %>%
    dplyr::mutate(Row = rep(1:combo_rows, each = combo_columns))
  
  final_list <- list(df, metadata)
  
  ## Print data frame and metadata
  return(final_list)
}

get_logistic_model <- function(data_frame, num_params){
  function_ret = nplr::nplr(data_frame$ConcCol, data_frame$Response, silent = TRUE, method = "res", LPweight = 0, npars = num_params) 
}

## inflection points from the n-parameter logistic curve. the n needs to be defined as argument 'y'
print_logistic_function <- function(input_data_frame, Logistic_Output_Folder = "logistic_curves", filename = NULL, num_param_logistic = "all"){
  x = input_data_frame
  fitting_func_df = x[[1]] %>%
    dplyr::select(Response, ConcCol, ConcRow) %>%
    dplyr::filter(ConcCol != 0) %>%
    dplyr::mutate(logConcCol = log10(ConcCol))
  fitting_list = split(fitting_func_df, fitting_func_df$ConcRow)
  
  ##Create overlappying n-parameter logistic curves list and extract inflection points
  
  Models <- lapply(fitting_list, function(tmp){nplr::nplr(tmp$ConcCol, tmp$Response, silent = TRUE, method = "res", LPweight = 0, npars = num_param_logistic)})
  infl_list <- lapply(Models, function(tmp){nplr::getInflexion(tmp)})
  infl_pts <- as.data.frame(unlist(infl_list),use.names= F)
  infl_pts <- infl_pts %>%
    dplyr::mutate(drugBconc = rep(names(infl_list), each = 2)) %>%
    dplyr::mutate(xy = rep(c("x","y"),times = length(infl_list))) %>%
    dplyr::filter(infl_pts, xy == "x")
  infl_pts$logconc<-as.numeric(infl_pts$`unlist(infl_list)`)
  infl_pts$conc <- 10^(infl_pts$logconc)
  
  ## define metadata to us in plots
  cellline <- x[[2]][["Cell_Line"]]
  DrugA <- x[[2]][["DrugA"]]
  DrugB <- x[[2]][["DrugB"]]
  DrugA.conc <- x[[2]][["DrugA.conc"]]
  DrugB.conc <- x[[2]][["DrugB.conc"]]
  Replicate <-x[[2]][["Replicate"]]
  inh_type <- x[[2]][["type"]]
  
  mytitle <- paste(cellline, DrugA, "with", DrugB, "Rep", Replicate, sep = "_")
  xlabel_1 <- paste("[",DrugA,"]","(",DrugA.conc,")")
  xlabel_2 <- paste("[",DrugB,"]","(",DrugB.conc,")")
  ylabel_1 <-paste("Inflection point of",DrugA,"(",DrugA.conc,")")
  ylabel_2 <- if (inh_type == "viability"){"% viability"} else if(inh_type == "inhibition"){"% Inhibition (response - min) / (max - min))"}
  
  ##save plots as tiffs
  
  
  if (is.null(filename)){
    filename = paste(mytitle,"_dose_response", sep = "")
  }
  
  if (is.null(Logistic_Output_Folder)){
    jpeg(file=paste(filename, ".jpeg", sep = ""), quality = 100)
  } else {
    if(!dir.exists(Logistic_Output_Folder)){dir.create(Logistic_Output_Folder); print(paste("Directory", Logistic_Output_Folder, "did not exist and was created"))}
    jpeg(file=paste(Logistic_Output_Folder, "/", filename, ".jpeg", sep = ""), quality = 100)
  }
  device = dev.cur()
  nplr::overlay(Models, xlab = bquote(paste('Log'[10],.(xlabel_1))),
          ylab = ylabel_2, main = mytitle, cex.main=1.5, cex.lab=1.5, font = 2, lwd = 4,
          Cols = RColorBrewer::brewer.pal(n = length(Models),name = "Oranges"), pty = 'm', xlog = T, xaxs = "r")
  dev.off(which = device)
}

## inflection points from the n-parameter logistic curve. the n needs to be defined as argument 'y'
print_logistic_inflection_points <- function(input_data_frame, Inflection_Output_Folder = "inflection_points", filename = NULL, num_param_logistic = "all"){
  x = input_data_frame
  fitting_func_df = x[[1]] %>%
    dplyr::select(Response, ConcCol, ConcRow) %>%
    dplyr::filter(ConcCol != 0) %>%
    dplyr::mutate(logConcCol = log10(ConcCol))
  fitting_list = split(fitting_func_df, fitting_func_df$ConcRow)
  
  ##Create overlappying n-parameter logistic curves list and extract inflection points
  
  Models <- lapply(fitting_list, function(tmp){nplr::nplr(tmp$ConcCol, tmp$Response, silent = TRUE, method = "res", LPweight = 0, npars = num_param_logistic)})
  infl_list <- lapply(Models, function(tmp){nplr::getInflexion(tmp)})
  infl_pts <- as.data.frame(unlist(infl_list),use.names= F)
  infl_pts <- infl_pts %>%
    dplyr::mutate(drugBconc = rep(names(infl_list), each = 2)) %>%
    dplyr::mutate(xy = rep(c("x","y"),times = length(infl_list))) %>%
    dplyr::filter(infl_pts, xy == "x")
  infl_pts$logconc<-as.numeric(infl_pts$`unlist(infl_list)`)
  infl_pts$conc <- 10^(infl_pts$logconc)
  
  ## define metadata to us in plots
  cellline <- x[[2]][["Cell_Line"]]
  DrugA <- x[[2]][["DrugA"]]
  DrugB <- x[[2]][["DrugB"]]
  DrugA.conc <- x[[2]][["DrugA.conc"]]
  DrugB.conc <- x[[2]][["DrugB.conc"]]
  Replicate <-x[[2]][["Replicate"]]
  inh_type <- x[[2]][["type"]]
  
  mytitle <- paste(cellline, DrugA, "with", DrugB, "Rep", Replicate, sep = "-")
  xlabel_1 <- paste("[",DrugA,"]","(",DrugA.conc,")")
  xlabel_2 <- paste("[",DrugB,"]","(",DrugB.conc,")")
  ylabel_1 <-paste("Inflection point of",DrugA,"(",DrugA.conc,")")
  ylabel_2 <- if (inh_type == "viability"){"% viability"} else if(inh_type == "inhibition"){"% Inhibition (response - min) / (max - min))"}
  
  if (is.null(filename)){
    filename = paste(mytitle,"_dose_response", sep = "")
  }
  
  if (is.null(Inflection_Output_Folder)){
    jpeg(file=paste(filename, ".jpeg", sep = ""), quality = 100)
  } else {
    if(!dir.exists(Inflection_Output_Folder)){dir.create(Inflection_Output_Folder); print(paste("Directory", Inflection_Output_Folder, "did not exist and was created"))}
    jpeg(file=paste(Inflection_Output_Folder, "/", filename, ".jpeg", sep = ""), quality = 100)
  }
  device = dev.cur()
  print(ggplot2::ggplot(infl_pts, ggplot2::aes(x=as.numeric(drugBconc),y=conc)) + ggplot2::geom_point(color = "black", size = 4) + ggplot2::geom_line(color = "red", size = 1, linetype = "dashed") + ggplot2::labs(title = paste(num_param_logistic,"parameter logistic curve of", mytitle)) + ggplot2::xlab(paste(xlabel_2)) + ggplot2::ylab(paste(ylabel_1)) + ggplot2::theme_classic())
  dev.off(which = device)
}

synergy_analysis <- function(input_data_frame,synergy_type = "Bliss", Output_Folder = "synergy_folder", filename = NULL){
  x = input_data_frame
  ## create directory for files
  
  
  
  ## define metadata to use in plots
  
  cellline <- x[[2]][["Cell_Line"]]
  DrugA <- x[[2]][["DrugA"]]
  DrugB <- x[[2]][["DrugB"]]
  DrugA.conc <- x[[2]][["DrugA.conc"]]
  DrugB.conc <- x[[2]][["DrugB.conc"]]
  Replicate <-x[[2]][["Replicate"]]
  inh_type <-x[[2]][["type"]]
  
  mytitle <- paste(cellline, DrugA, "with", DrugB, "Rep", Replicate, sep = "_")
  
  ## Calculate and plot synergy
  df <- as.data.frame(x[[1]])
  dose.response.mat <- synergyfinder::ReshapeData(df, data.type = inh_type)
  synergy.score <- synergyfinder::CalculateSynergy(data = dose.response.mat,
                                    method = paste(synergy_type),
                                    correction = F,
                                    Emin = 0, Emax = 100)
  
  if (is.null(filename)){
    filename = paste(mytitle,synergy_type, sep = "_")
  }
  
  if (is.null(Output_Folder)){
    jpeg(file=paste(filename, ".jpeg", sep = ""), quality = 100)
  } else {
    if(!dir.exists(Output_Folder)){dir.create(Output_Folder); print(paste("Directory", Output_Folder, "did not exist and was created"))}
    jpeg(file=paste(Output_Folder, "/", filename, ".jpeg", sep = ""), quality = 100)
  }
  device = dev.cur()
  synergyfinder::PlotSynergy(synergy.score, type = "3D", len = 5, save.file = FALSE)
  dev.off(which = device)
  
  # return the average synergy score
  
  synergy_scores <- synergy.score[["scores"]][[1]]
  synergy_scores <- synergy_scores[2:nrow(synergy_scores),2:ncol(synergy_scores)]
  average_score <- mean(synergy_scores)
  names(average_score) <- paste(mytitle)
  if (is.null(Output_Folder)){
    write.csv(unlist(average_score),file=paste(filename, "_Score" ,".csv", sep = ""))
  } else {
    write.csv(unlist(average_score), file = paste(Output_Folder, "/", filename,"_Score",".csv", sep = ""))
  }
}

Combine_Synergy_Files <- function(Input_Directory){
  files = collect_file_vector(Input_Directory, "csv")
  name_vector = c()
  value_vector = c()
  for (i in files){
    current_file = file_path_sans_ext(i)
    name_vector = c(name_vector, current_file)
    current_data = read.csv(paste(Input_Directory, "/",i, sep = ""),header = FALSE)
    value_vector = c(value_vector, as.character(current_data[2,2]))
  }
  return_frame = data.frame(name_vector, value_vector)
  colnames(return_frame) <- c("Condition", "Synergy_Score")
  return(return_frame)
}

Combine_Synergy_Files_CSV <- function(Input_Directory, Output_Directory, Output_Filename = "Combined_Synergy", allow_overwrite = FALSE ,remove_csv = FALSE){
  files = collect_file_vector(Input_Directory, "csv")
  name_vector = c()
  value_vector = c()
  for (i in files){
    current_file = file_path_sans_ext(i)
    name_vector = c(name_vector, current_file)
    current_data = read.csv(paste(Input_Directory, "/",i, sep = ""),header = FALSE)
    if(remove_csv){file.remove(paste(Input_Directory, "/",i, sep = ""))}
    value_vector = c(value_vector, as.character(current_data[2,2]))
  }
  return_frame = data.frame(name_vector, value_vector)
  if (!dir.exists(Output_Directory)){dir.create(Output_Directory)}
  if (!file.exists(paste(Output_Directory, "/", Output_Filename, ".csv", sep = "")) | allow_overwrite){
    colnames(return_frame) <- c("Condition", "Synergy_Score")
    write.csv(return_frame, file = paste(Output_Directory, "/", Output_Filename, ".csv", sep = ""))
  } else {
    print("File Already Exists. Please Use a different File Name")
  }
  
}

#Functions for the pipeline function to use
pipelined_Logistic_Curves <- function(InputFilePath, OutputFileName, data_location, Normalize = TRUE, Viability_Data = TRUE, num_param = "all") {
  data_frames <- try(build_data_frame(InputFilePath, range=data_location, Viability_Data = Viability_Data, 
                                  Normalize = Normalize, Drug_A = "Drug_A", Drug_B = "Drug_B",
                                  Cell_Line = "Cell_Line", Conc_A = "uM", Conc_B = "uM", replicate = 1))
    
  try(print_logistic_function(data_frames, Logistic_Output_Folder = NULL, num_param_logistic = num_param, filename = paste(OutputFileName, "_Logistic", sep = "")))


}

pipelined_Inflection_Curves <- function(InputFilePath, OutputFileName, data_location, Normalize = TRUE, Viability_Data = TRUE, num_param = "all") {
  data_frames <- try(build_data_frame(InputFilePath, range=data_location, Viability_Data = Viability_Data, 
                                  Normalize = Normalize, Drug_A = "Drug_A", Drug_B = "Drug_B",
                                  Cell_Line = "Cell_Line", Conc_A = "uM", Conc_B = "uM", replicate = 1))
  try(print_logistic_inflection_points(data_frames, Inflection_Output_Folder = NULL, num_param_logistic = num_param, filename = paste(OutputFileName, "_Inflection", sep = "")))
}

pipelined_Synergy_Graphs <- function(InputFilePath, OutputFileName, data_location, synergy_type = "Bliss", Normalize = TRUE, Viability_Data = TRUE) {
  data_frames <- try(build_data_frame(InputFilePath, range=data_location, Viability_Data = Viability_Data, 
                                  Normalize = Normalize, Drug_A = "Drug_A", Drug_B = "Drug_B",
                                  Cell_Line = "Cell_Line", Conc_A = "uM", Conc_B = "uM", replicate = 1))
  try(synergy_analysis(input_data_frame = data_frames, synergy_type = synergy_type, filename = paste(OutputFileName, "_Synergy", sep = ""),Output_Folder = NULL))
  if (file.exists("Rplots.pdf")){file.remove("Rplots.pdf")}
}

