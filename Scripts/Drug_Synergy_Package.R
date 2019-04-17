if (!require("drc")) {
  install.packages("drc", dependencies = TRUE)
  library(drc)
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

if (!require("readxl")) {
  install.packages("tidyverse", dependencies = TRUE)
  library(readxl)
}

library(tools)

read_xlsx <- readxl::read_xlsx
heatmap.2 <- gplots::heatmap.2
file_ext <- tools::file_ext

import_plate <- function(FileName) {
  ext = file_ext(FileName)
  if (ext == "xlsx"){
    df <- as.data.frame(read_xlsx(FileName, col_names = TRUE))
    rows <- as.vector(df[,1])
    df[,1] <- NULL
    rownames(df) <- rows
    return(as.matrix(df))
  } else if (ext == "csv"){
    #range only supported in xlsx files
    df <- read.csv2(FileName, sep = ',', check.names = FALSE)
    rows <- as.vector(df[,1])
    cols <- colnames(df[,2:ncol(df)])
    df[,1] <- NULL
    rownames(df) <- rows
    matrix <- as.matrix(df)
    matrix <- mapply(matrix, FUN=as.numeric)
    matrix <- matrix(matrix, ncol = ncol(df), nrow = nrow(df))
    colnames(matrix) = cols
    rownames(matrix) = rows
    return(matrix)
  } 
  return(NULL)
}

import_plate_range <- function(FileName, ImportRange) {
  ext = file_ext(FileName)
  if (ext == "xlsx"){
    df <- as.data.frame(read_xlsx(FileName, col_names = TRUE, range = ImportRange))
    rows <- as.vector(df[,1])
    df[,1] <- NULL
    rownames(df) <- rows
    return(as.matrix(df))
  } else if (ext == "csv"){
    #range only supported in xlsx files
    df <- read.csv2(FileName, sep = ',', check.names = FALSE)
    rows <- as.vector(df[,1])
    cols <- colnames(df[,2:ncol(df)])
    df[,1] <- NULL
    rownames(df) <- rows
    matrix <- as.matrix(df)
    matrix <- mapply(matrix, FUN=as.numeric)
    matrix <- matrix(matrix, ncol = ncol(df), nrow = nrow(df))
    colnames(matrix) = cols
    rownames(matrix) = rows
    return(matrix)
  } 
  return(NULL)
}

un_craig_plate <- function(Plate) {
  rownames(Plate) = rev(as.vector(rownames(Plate)))
  colnames(Plate) = rev(as.vector(colnames(Plate)))
  
  fixed_plate <- matrix(rep(0, times = ncol(Plate) * nrow(Plate)),ncol = ncol(Plate), nrow = nrow(Plate))
  row_index <- 1
  column_index <- 1
  for ( i in nrow(Plate):1 ) {
    for (j in ncol(Plate):1) {
      fixed_plate[row_index,column_index] <- Plate[i,j]
      column_index = column_index + 1
    }
    column_index = 1
    row_index = row_index + 1
  }
  
  rownames(fixed_plate) <- as.vector(rownames(Plate))
  colnames(fixed_plate) <- as.vector(colnames(Plate))
  return(fixed_plate)
}

growth_effect_calculator <- function(raw_growth, column_to_normalize_by = 1, row_to_normalize_by = 1){
  growth_plate <- raw_growth/raw_growth[row_to_normalize_by,column_to_normalize_by]
  return(growth_plate)
}


kill_effect_matrix <- function(growth_effect_matrix, CraigFactor = FALSE){
  fold_over = growth_effect_matrix
  if (CraigFactor == TRUE) {
    fold_over = un_craig_plate(growth_effect_matrix)
  }
  Kill_effect = 1-fold_over
  return(Kill_effect)
}


bliss_calculator <- function(growth_effect_matrix,CraigFactor = FALSE) {
  bliss <- matrix(ncol = ncol(growth_effect_matrix) - 1,nrow = nrow(growth_effect_matrix) - 1)
  Kill_effect <- kill_effect_matrix(growth_effect_matrix,CraigFactor)
  
  row_index = 1
  column_index = 1
  for(i in 2:nrow(Kill_effect)) {
    for(j in 2:ncol(Kill_effect)) {
      fr <- Kill_effect[i,1]
      fb <- Kill_effect[1,j]
      union = fr * fb
      bliss[row_index,column_index] <- ((fr + fb) - union)/Kill_effect[i,j]
      column_index = column_index + 1
    }
    column_index = 1
    row_index = row_index + 1
    
  }
  
  row_names = as.vector(rownames(Kill_effect))
  col_names = as.vector(colnames(Kill_effect))
  
  row_names <- row_names[-1]
  col_names <- col_names[-1]
  
  rownames(bliss) <- row_names
  colnames(bliss) <- col_names
  
  return(bliss)
}


print_growth_curves <- function(growth_effect_matrix,jpeg_name = "default.jpg", numCurves = 2,CraigFactor = FALSE) {
  df1 <- data.frame(as.numeric(as.vector(colnames(growth_effect_matrix))))
  
  colors <- brewer.pal(nrow(growth_effect_matrix),"Set1")
  
  if(numCurves > nrow(growth_effect_matrix) || numCurves < 0) {
    return()
  }
  
  for(i in 1:numCurves) {
    df1 <- cbind(df1,as.vector(growth_effect_matrix[i,]))
  }
  
  jpeg(filename=jpeg_name,res=600,height = 10,width = 10,units = "in")
  index = 2;
  for(i in 1:numCurves) {
    mL <- drm(df1[,index] ~ df1[,1],data = df1 ,fct = LL.4(fixed = c(NA,NA,NA,NA)), type="continuous")
    if (i == 1) {
      plot(mL,type="all",col = colors[i],pch = i,lwd = 5,xlab = NA,ylab = NA,main = "Combination Points",ylim= c(0,2))
    }else {
      plot(mL,type="all",col = colors[i],add=TRUE,pch = i,lwd = 3,xlab = NA,ylab = NA)
    }
    index = index + 1
  }
  
  legend("topright",lwd = 2,col = colors,pch = c(1:numCurves),legend=as.character(rownames(growth_effect_matrix)[1:numCurves]))
  
  dev.off()
}


print_growth_curves_error <- function(growth_effect_matrix, error_matrix,jpeg_name = "default.jpg", numCurves = 2,CraigFactor = FALSE) {
  df1 <- data.frame(as.numeric(as.vector(colnames(growth_effect_matrix))))
  df2 <- data.frame(as.numeric(as.vector(colnames(error_matrix))))
  df2[1,1] = .001
  
  colors <- brewer.pal(nrow(growth_effect_matrix),"Set1")
  
  if(numCurves > nrow(growth_effect_matrix) || numCurves < 0) {
    return()
  }
  
  for(i in 1:numCurves) {
    df1 <- cbind(df1,as.vector(growth_effect_matrix[i,]))
    df2 <- cbind(df2,as.vector(error_matrix[i,]))
  }
  
  
  jpeg(filename=jpeg_name,res=600,height = 10,width = 10,units = "in")
  index = 2;
  for(i in 1:numCurves) {
    mL <- drm(df1[,index] ~ df1[,1],data = df1 ,fct = LL.4(fixed = c(NA,NA,NA,NA)), type="continuous")
    if (i == 1) {
      plot(mL,type="all",col = colors[i],pch = i,lwd = 5,xlab = NA,ylab = NA,main = "Combination Points",ylim= c(0,2))
      arrows(df2[,1], df1[,index] - df2[,index], df2[,1], df1[,index] + df2[,index], length=0.05, angle=90, code=3, col = colors[i])
    }else {
      plot(mL,type="all",col = colors[i],add=TRUE,pch = i,lwd = 3,xlab = NA,ylab = NA)
      arrows(df2[,1], df1[,index] - df2[,index], df2[,1], df1[,index] + df2[,index], length=0.05, angle=90, code=3, col = colors[i])
    }
    index = index + 1
  }
  dev.off()
}


build_parsed_bliss_map <- function(growth_effect_matrix, max_gi50 = 0.95, min_gi50 = 0.05){  
  data_bliss <- bliss_calculator(growth_effect_matrix)
  cols_to_keep = c()
  row_to_keep = c()
  
  for (item in 1:ncol(growth_effect_matrix)){
    curr_item = growth_effect_matrix[1,item]
    if( curr_item < max_gi50 && curr_item > min_gi50) {
      cols_to_keep = c(cols_to_keep,item)
    }
  }
  
  for (item in 1:nrow(growth_effect_matrix)){
    curr_item = growth_effect_matrix[item, 1]
    if( curr_item < max_gi50 && curr_item > min_gi50) {
      row_to_keep = c(row_to_keep,item)
    }
  }
  row_to_keep = row_to_keep - 1
  cols_to_keep = cols_to_keep - 1
  
  keep <- data_bliss[row_to_keep, cols_to_keep]
  return(keep)
}




print_heatmap_bliss_scaled <- function(growth_effect_matrix,export_name="default_heatmap.jpg",xlabel="x-Axis",ylabel="y-Axis") {
  data_plate = bliss_calculator(growth_effect_matrix)
  
  color_range = 399
  
  jpeg(filename=export_name,res=600,height = 12,width = 12,units = "in")
  
  Synergy <- colorRampPalette(c("red","white","Blue"))(n = color_range)
  colors_pal = sort(c(seq(-1,-0.001, length=200), seq(0,1,length=200)))

  for (i in 1:nrow(data_plate)) {
    for(j in 1:ncol(data_plate)) {
      current_value <- data_plate[i,j]
      if (current_value <= 0) {
        data_plate[i,j] = -1
      } else if (current_value > 2) {
        data_plate[i,j] = 1
      } else {
        data_plate[i,j] = current_value - 1
      } 
    }
  }
  
  map = heatmap.2(data_plate, dendrogram = "none",
            Rowv = NA, 
            Colv = NA, 
            xlab = xlabel, 
            ylab = ylabel,
            col = Synergy,
            trace='none',
            density.info = c("none"),
            keysize = 1,
            symkey = FALSE,
            breaks = colors_pal)

  dev.off()
}

print_heatmap_bliss_outlier_elim <- function(growth_effect_matrix,export_name="default_heatmap.jpg",xlabel="x-Axis",ylabel="y-Axis") {
  data_plate = bliss_calculator(growth_effect_matrix)
  color_range = 399
  jpeg(filename=export_name,res=600,height = 12,width = 12,units = "in")
  
  Synergy <- colorRampPalette(c("red","white","Blue"))(n = color_range)
  colors_pal = sort(c(seq(-1,-0.001, length=200), seq(0,1,length=200)))
  
  for (i in 1:nrow(data_plate)) {
    for(j in 1:ncol(data_plate)) {
      current_value <- data_plate[i,j]
      if (current_value < 0) {
        data_plate[i,j] = 0
      } else if (current_value > 2) {
        data_plate[i,j] = 0
      } else {
        data_plate[i,j] = current_value - 1
      } 
    }
  }
  
  map = heatmap.2(data_plate, dendrogram = "none",
            Rowv = NA, 
            Colv = NA, 
            xlab = xlabel, 
            ylab = ylabel,
            col = Synergy,
            trace='none',
            density.info = c("none"),
            keysize = 1,
            symkey = FALSE,
            breaks = colors_pal)
  
  dev.off()
}

calculate_bliss_sum_scaled <- function(growth_effect_matrix) {
  data_plate = bliss_calculator(growth_effect_matrix)
  for (i in 1:nrow(data_plate)) {
    for(j in 1:ncol(data_plate)) {
      current_value <- data_plate[i,j]
      if (current_value <= 0) {
        data_plate[i,j] = -1
      } else if (current_value > 2) {
        data_plate[i,j] = 1
      } else {
        data_plate[i,j] = current_value - 1
      } 
    }
  }
  
  bliss_sum = 0
  for (i in 1:nrow(data_plate)) {
    for(j in 1:ncol(data_plate)) {
      current_value <- data_plate[i,j]
      bliss_sum = current_value + bliss_sum
    }
  }
  return(bliss_sum)
}

calculate_bliss_sum_outlier_elim <- function(growth_effect_matrix) {
  data_plate = bliss_calculator(growth_effect_matrix)
  for (i in 1:nrow(data_plate)) {
    for(j in 1:ncol(data_plate)) {
      current_value <- data_plate[i,j]
      if (current_value <= 0) {
        data_plate[i,j] = 0
      } else if (current_value > 2) {
        data_plate[i,j] = 0
      } else {
        data_plate[i,j] = current_value - 1
      } 
    }
  }
  
  bliss_sum = 0
  for (i in 1:nrow(data_plate)) {
    for(j in 1:ncol(data_plate)) {
      current_value <- data_plate[i,j]
      bliss_sum = current_value + bliss_sum
    }
  }
  return(bliss_sum)
}


parse_bliss_plate <- function(data_plate) {
  for (i in 1:nrow(data_plate)) {
    for(j in 1:ncol(data_plate)) {
      current_value <- data_plate[i,j]
      if (current_value < 0) {
        data_plate[i,j] = 0
      } else if (current_value > 2) {
        data_plate[i,j] = 1
      } else {
        data_plate[i,j] = current_value - 1
      } 
    }
  }
  return(data_plate)
}



print_heatmap_bliss_scaled_parsed <- function(growth_plate,export_name="default_heatmap.jpg", GI_Max = 0.95, GI_Min = 0.05, xlabel="x-Axis",ylabel="y-Axis") {
  #Get Bliss Scores
  data_plate <- build_parsed_bliss_map(growth_plate, max_gi50 = GI_Max, min_gi50 = GI_Min)
  color_range = 399
  
  jpeg(filename=export_name,res=600,height = 12,width = 12,units = "in")
  
  Synergy <- colorRampPalette(c("red","white","Blue"))(n = color_range)
  colors_pal = sort(c(seq(-1,-0.001, length=200), seq(0,1,length=200)))
  
  for (i in 1:nrow(data_plate)) {
    for(j in 1:ncol(data_plate)) {
      current_value <- data_plate[i,j]
      if (current_value < 0) {
        data_plate[i,j] = 0
      } else if (current_value > 2) {
        data_plate[i,j] = 1
      } else {
        data_plate[i,j] = current_value - 1
      } 
    }
  }
  
  map = heatmap.2(data_plate, dendrogram = "none",
            Rowv = NA, 
            Colv = NA, 
            xlab = xlabel, 
            ylab = ylabel,
            col = Synergy,
            trace='none',
            density.info = c("none"),
            keysize = 1,
            symkey = FALSE,
            breaks = colors_pal)
  
  dev.off()
  return(map)
}

print_heatmap_bliss_outlier_elim_parsed <- function(growth_plate,export_name="default_heatmap.jpg", GI_Max = 0.95, GI_Min = 0.05, xlabel="x-Axis",ylabel="y-Axis") {
  data_plate <- build_parsed_bliss_map(growth_plate, max_gi50 = GI_Max, min_gi50 = GI_Min)
  color_range = 399
  
  jpeg(filename=export_name,res=600,height = 12,width = 12,units = "in")
  
  Synergy <- colorRampPalette(c("red","white","Blue"))(n = color_range)
  colors_pal = sort(c(seq(-1,-0.001, length=200), seq(0,1,length=200)))
  
  for (i in 1:nrow(data_plate)) {
    for(j in 1:ncol(data_plate)) {
      current_value <- data_plate[i,j]
      if (current_value < 0) {
        data_plate[i,j] = 0
      } else if (current_value > 2) {
        data_plate[i,j] = 0
      } else {
        data_plate[i,j] = current_value - 1
      } 
    }
  }
  
  heatmap.2(data_plate, dendrogram = "none",
            Rowv = NA, 
            Colv = NA, 
            xlab = xlabel, 
            ylab = ylabel,
            col = Synergy,
            trace='none',
            density.info = c("none"),
            keysize = 1,
            symkey = FALSE,
            breaks = colors_pal)
  
  dev.off()
  mapped_data = data_plate
  return(mapped_data)
}


pipelined_print_heatmap_growth_scaled <- function(InputFilePath, OutputFileName, data_location) {
  OutputName = paste(OutputFileName, ".jpg",sep = "")
  print_heatmap_bliss_scaled(import_plate_range(InputFilePath, data_location), OutputName)
}

pipelined_print_heatmap_growth_outlier_elim <- function(InputFilePath, OutputFileName, data_location) {
  OutputName = paste(OutputFileName, ".jpg",sep = "")
  print_heatmap_bliss_outlier_elim(import_plate_range(InputFilePath, data_location), OutputName)
}

pipelined_bliss_heatmap_growth_scaled_parsed <- function(InputFilePath, OutputFileName, data_location, GI_Max = 0.95, GI_Min = 0.05) {
  OutputName = paste(OutputFileName, ".jpg",sep = "")
  print_heatmap_bliss_scaled_parsed(import_plate_range(InputFilePath, data_location), OutputName, GI_Max, GI_Min)
}

pipelined_bliss_heatmap_growth_outlier_elim_parsed <- function(InputFilePath, OutputFileName, data_location, GI_Max = 0.95, GI_Min = 0.05) {
  OutputName = paste(OutputFileName, ".jpg",sep = "")
  print_heatmap_bliss_outlier_elim_parsed(import_plate_range(InputFilePath, data_location), OutputName, GI_Max, GI_Min)
}

pipelined_print_bliss_sum_growth_scaled <- function(InputFilePath, OutputFileName, data_location, combined_file="default_file.csv") {
  current_file <- NULL
  
  bliss_sum <- calculate_bliss_sum_scaled(import_plate_range(InputFilePath, data_location))
  if(file.exists(combined_file)){ ## requires that we overwrite our "final_file_name" for each addition
    current_file = read.csv(combined_file, header = TRUE, na.strings = c("",0,"<NA>","NA"), row.names = 1)
  } else { ## this is for processing the first data frame/heatmap column
    current_file = data.frame(row.names = c("Bliss_Score_Sum"))
  }
  current_file$default <- bliss_sum
  names(current_file)[ncol(current_file)]<-OutputFileName
  write.csv(current_file, file = combined_file)
}

pipelined_print_bliss_sum_growth_outlier_elim <- function(InputFilePath, OutputFileName, data_location, combined_file="default_file.csv") {
  current_file <- NULL
  
  bliss_sum <- calculate_bliss_sum_outlier_elim(import_plate_range(InputFilePath, data_location))
  if(file.exists(combined_file)){ ## requires that we overwrite our "final_file_name" for each addition
    current_file = read.csv(combined_file, header = TRUE, na.strings = c("",0,"<NA>","NA"), row.names = 1)
  } else { ## this is for processing the first data frame/heatmap column
    current_file = data.frame(row.names = c("Bliss_Score_Sum"))
  }
  current_file$default <- bliss_sum
  names(current_file)[ncol(current_file)]<-OutputFileName
  write.csv(current_file, file = combined_file)
}

pipelined_uncraig_batch_normalize <- function(InputFilePath, OutputFileName, data_location) {
  OutputName = paste(OutputFileName, ".jpg",sep = "")
  plate <- un_craig_plate(import_plate_range(InputFilePath, data_location))
  plate <- growth_effect_calculator(plate)
  write.csv(plate, file = paste(OutputFileName,".csv", sep = ""))
}

