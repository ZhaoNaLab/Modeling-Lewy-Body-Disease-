library("batchtools")

## use the bath tools for RNA contamination detection and doublet removal
setwd("./Step5.FinalClustering")
reg <- makeExperimentRegistry("FindMarker.Registry",
                              packages = c("tidyverse", "Seurat", "data.table"), 
                              seed = 20221005)

reg <- loadRegistry("FindMarker.Registry", writeable = TRUE)

Dir <- "./Step5.FinalClustering/Seurat.SNC.final.rds"

## Define instance
# define instance (data preparation)
subsample <- function(data, job, i) {
  MajorCluster <- c("Ex", "In", "Ast", "Olig", "OPC", "Mic", "Endo", "Peri")
  Cluster <- MajorCluster[i]
  data <- readRDS(data)
  DefaultAssay(data) <- "RNA"
  list(data = data, Cluster = Cluster)  # obtain different parameters that can be parsed to the next step
}

# add the dataset for this problem; this step can either generate a static object or 
addProblem(name = "ClusterMarkers", 
           data = Dir, 
           fun = subsample) 

# Define tasks
FindClusterMarker <- function(data, job, instance) { # instance is inherited from the last step
  data = instance$data
  Cluster = instance$Cluster
  
  ####### FindMarkers #######
  Idents(data) <- data$MajorCluster
  Markers <- FindMarkers(data, slot = 'counts',
                         ident.1 = Cluster,
                         logfc.threshold = 0.25,
                         min.pct = 0.1,
                         latent.vars = c("cngeneson", "decontX_contamination", "percent.mt"),
                         test.use = "MAST",
                         only.pos = FALSE)
  
  write.csv(Markers, paste0("CellMarkers/", Cluster, ".csv"))
}

addAlgorithm(name = "MarkerFinder", fun = FindClusterMarker)

# Creating jobs
# problem design: try two values for the ratio parameter
prob.design <- CJ(i = 1:8)

# final problem design
pdes <- list(ClusterMarkers = prob.design)

# algorithm design: try combinations of kernel and epsilon exhaustively,
# try different number of trees for the forest
ades <- list(MarkerFinder = data.table())

addExperiments(pdes, ades) # pdes set parameter to the data part; ades set parameters to the instance parts

# summarize the experiment
summarizeExperiments()

# Submitting and Collecting Results
options(batchtools.progress = FALSE)
submitJobs() # ids = c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56)
getStatus()
getJobTable()
getErrorMessages()
findExpired() 

