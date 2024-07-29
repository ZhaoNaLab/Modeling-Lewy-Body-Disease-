library("tidyverse")
library("batchtools")
library("Seurat")
library("data.table")
library("future")
library("MAST")
library("SingleCellExperiment")

## use the bath tools for RNA contamination detection and doublet removal
setwd("./Step5.FinalClustering")

# load the seurat object
reg <- makeExperimentRegistry("FindMarker.MAST.Registry",
                              packages = c("tidyverse", "Seurat", "data.table", 
                                           "future", "MAST", "SingleCellExperiment"), 
                              seed = 20221005)

reg <- loadRegistry("FindMarker.MAST.Registry", writeable = TRUE)


Dir <- "./Step5.FinalClustering/scaObjs.rds"

## Define instance
# define instance (data preparation)
subsample <- function(data, job, i) {
  data <- readRDS(data)
  Cluster <- unique(data@colData$MajorCluster)[i]
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
  
  setwd("./Step5.FinalClustering")
 
  ####### FindMarkers #######
  Feature.Cell.Percent <- readRDS("Feature.Cell.Percent.rds") %>% 
    as.data.frame() 
  
  Pos <- grep(Cluster, names(Feature.Cell.Percent))

  ## Select features expressed by at least 25% of the target cells
  FeatureSelected <- rownames(Feature.Cell.Percent)[Feature.Cell.Percent[, Pos] >= 25]
  FeatureSelected <- intersect(FeatureSelected, rownames(data))
  data <- data[FeatureSelected, ]
  
  ## Make a new varible CulsterOfInterest for Marker calculation
  colData(data)$CulsterOfInterest <- "Ref"
  colData(data)$CulsterOfInterest[colData(data)$MajorCluster == Cluster] <- "Target"
  colData(data)$CulsterOfInterest <- factor(colData(data)$CulsterOfInterest,
                                            level = c("Ref", "Target"))
  
  # MAST model for marker identification 
  zlmCond <- zlm(~ CulsterOfInterest + cngeneson + percent.mt + decontX_contamination,
                 parallel = 4,
                 sca = data)
  
  summaryCond <- summary(zlmCond, doLRT= c('CulsterOfInterestTarget'))
  # print out the data.table
  summaryDt <- summaryCond$datatable
  Markers <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                    summaryDt[component=='logFC', .(primerid, coef, contrast)], 
                    by=c('primerid', 'contrast')) #logFC coefficients
  
  Markers$contrast <- Cluster
  
  # add expression level information
  Markers <- Markers %>% 
    left_join(Feature.Cell.Percent %>% rownames_to_column("primerid"), by = "primerid")
  
  # save the result
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
# testJob(id, external = FALSE, reg = getDefaultRegistry())
submitJobs() # ids = c(1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56)
getStatus()
getJobTable()
getErrorMessages()
findExpired() 
# testJob(id = 3)
# killJobs()

