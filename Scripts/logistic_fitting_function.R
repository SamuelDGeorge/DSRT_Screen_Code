##This function will take a dataframe shaped by final_list_method() to create an overlay of n-parameter logistic curves, and a plot of the

if(!dir.exists("logistic_curves")){dir.create("logistic_curves")}

## inflection points from the n-parameter logistic curve. the n needs to be defined as argument 'y'
logistic_function <- function(x,y,z){
  fitting_func_df = x[[1]] %>%
    select(Response, ConcCol, ConcRow) %>%
    filter(ConcCol != 0) %>%
    mutate(logConcCol = log10(ConcCol))
  fitting_list = split(fitting_func_df, fitting_func_df$ConcRow)

  ##Create overlappying n-parameter logistic curves list and extract inflection points
  
  Models <- lapply(fitting_list, function(tmp){nplr(tmp$ConcCol, tmp$Response, silent = TRUE, method = "res", LPweight = 0, npars = y)})
  infl_list <- lapply(Models, function(tmp){getInflexion(tmp)})
  infl_pts <- as.data.frame(unlist(infl_list),use.names= F)
  infl_pts <- infl_pts %>%
    mutate(drugBconc = rep(names(infl_list), each = 2)) %>%
    mutate(xy = rep(c("x","y"),times = length(infl_list))) %>%
    filter(infl_pts, xy == "x")
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
  
  mytitle <- paste(cellline, DrugA, "combination with", DrugB, "Rep", Replicate)
  xlabel_1 <- paste("[",DrugA,"]","(",DrugA.conc,")")
  xlabel_2 <- paste("[",DrugB,"]","(",DrugB.conc,")")
  ylabel_1 <-paste("Inflection point of",DrugA,"(",DrugA.conc,")")
  ylabel_2 <- if (inh_type == "viability"){"% viability"} else if(inh_type == "inhibition"){"% Inhibition (response - min) / (max - min))"}
 
   ##save plots as tiffs
     tiff(file=paste("logistic_curves/",mytitle,"dose response", ".tif"))
     overlay(Models, xlab = bquote(paste('Log'[10],.(xlabel_1))),
             ylab = ylabel_2, main = mytitle, cex.main=1.5, cex.lab=1.5, font = 2, lwd = 4,
             Cols = brewer.pal(n = length(Models),name = "Oranges"), pty = 'm', xlog = T, xaxs = "r")
     dev.off()
  if(z == "YES"){
    if(!dir.exists("inflection_points")){dir.create("inflection_points")}
    tiff(file=paste("inflection_points/",mytitle,"inflection points",".tif"))
    print(ggplot(infl_pts, aes(x=as.numeric(drugBconc),y=conc)) + geom_point(color = "black", size = 4) + geom_line(color = "red", size = 1, linetype = "dashed") + labs(title = paste(y,"parameter logistic curve of", mytitle)) + xlab(paste(xlabel_2)) + ylab(paste(ylabel_1)) + theme_classic())
    dev.off()}
}