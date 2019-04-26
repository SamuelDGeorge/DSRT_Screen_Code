synergy_analysis <- function(x,y){
  
  ## create directory for files
  
  if(!dir.exists("synergy_profiles")){dir.create("synergy_profiles")}
  
  ## define metadata to use in plots
  
  cellline <- x[[2]][["Cell_Line"]]
  DrugA <- x[[2]][["DrugA"]]
  DrugB <- x[[2]][["DrugB"]]
  DrugA.conc <- x[[2]][["DrugA.conc"]]
  DrugB.conc <- x[[2]][["DrugB.conc"]]
  Replicate <-x[[2]][["Replicate"]]
  inh_type <-x[[2]][["type"]]
  
  mytitle <- paste(cellline, DrugA, "combination with", DrugB, "Rep", Replicate)
  
  ## Calculate and plot synergy
df <- as.data.frame(x[[1]])
dose.response.mat <- ReshapeData(df, data.type = inh_type)
synergy.score <- CalculateSynergy(data = dose.response.mat,
                                        method = paste(y),
                                        correction = F,
                                        Emin = 0, Emax = 100)

png(file=paste("synergy_profiles/",mytitle,"Bliss analysis",".png"))
PlotSynergy(synergy.score, type = "3D", len = 5)
dev.off()

# return the average synergy score

synergy_scores <- synergy.score[["scores"]][[1]]
synergy_scores <- synergy_scores[2:nrow(synergy_scores),2:ncol(synergy_scores)]
average_score <- mean(synergy_scores)
names(average_score) <- paste(mytitle)

return(average_score)

}