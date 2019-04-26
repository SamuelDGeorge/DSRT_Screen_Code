##The following will get the drug concentrations as a vector, which will define the size of the drug matrix
final_list_method <- function(x,y,z){
  x <- read.csv(file = paste(x), header = F)
  DrugDoseA = na.omit(as.numeric(as.vector(x[1,2:11])))
  DrugDoseB = na.omit(as.numeric(as.vector(x[2:7,1])))
  matrix_x = length(DrugDoseA)
  matrix_y = length(DrugDoseB)
  
  ##Define metadata such as DrugA (horizontal titration), DrugB (vertical titration), cell line, and concentration units for each DrugA and              DrugB. Also assign the appropriate header names.
  
  metadata = x[2,12:17]
  names(metadata) = c("DrugA","DrugB", "Cell_Line","DrugA.conc","DrugB.conc","Replicate")
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
    unite(AplusB, DrugDoseArep, DrugDoseBrep, remove = T)
  
  ##Then we group the dataframe by repeat treatments and only cary the mean value forward. This will eliminate duplicated treatments
  
  df = df %>% 
    group_by(AplusB) %>%
    summarise(Response = mean(Response)) %>%
    separate(AplusB, c("ConcCol","ConcRow"), sep = "_")
  df$ConcCol = round(as.numeric(df$ConcCol), digits = 2)
  df$ConcRow = round(as.numeric(df$ConcRow), digits = 2)
  
  ## Find the control well and normalize
  if (y == "viability") {
      if(z == "YES") {
    ctrldf = df %>%
      filter(ConcCol==0,ConcRow==0)
    ctrl=ctrldf$Response[1]
    df$Response = df$Response/ctrl*100
      }
  }
  else if (y == "inhibition"){
      if(z == "YES") {
    max_response <- max(df$Response)
    min_response <- min(df$Response)
    df$Response = (df$Response-min_response)/(max_response-min_response)*100
      }
  }
  
  ## Now shape the data to the format needed for synergyFinder
  
  df = df %>%
    mutate(BlockID = 1) %>%
    mutate(DrugRow = as.character(metadata$DrugB[1])) %>%
    mutate(DrugCol = as.character(metadata$DrugA[1])) %>%
    mutate(ConcRowUnit = as.character(metadata$DrugB.conc)) %>%
    mutate(ConcColUnit = as.character(metadata$DrugA.conc)) %>%
    as.data.frame(df) %>%
    arrange(ConcCol) %>%
    arrange(ConcRow) %>%
    mutate(Col = rep(1:combo_columns, combo_rows)) %>%
    mutate(Row = rep(1:combo_rows, each = combo_columns))
  
  final_list <- list(df, metadata)
  
  ## Print data frame and metadata
  return(final_list)
}