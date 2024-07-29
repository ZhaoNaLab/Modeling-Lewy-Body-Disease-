library("batchtools")

## use the bath tools for RNA contamination detection and doublet removal
setwd("./Step4.IterativeCleaning")
reg <- makeExperimentRegistry("ReclusterMajorCluster.Registry",
                              packages = c("tidyverse", "Seurat", "data.table", "harmony"), 
                              seed = 20221005)

reg <- loadRegistry("ReclusterMajorCluster.Registry", writeable = TRUE)

Dir <- "./Step3.Clustering/Seurat.SNC.Ind.rds"

## Define instance
# define instance (data preparation)
subsample <- function(data, job, i) {
  MajorCluster <- c("Ex", "In", "Ast", "Olig", "OPC", "Mic", "Vascular")
  Cluster <- MajorCluster[i]
  data <- readRDS(data)
  data <- subset(data, MajorCluster == Cluster)
  data <- DietSeurat(data, assays = c("RNA", "decontX"), counts = TRUE)
  DefaultAssay(data) <- "RNA"
  list(data = data, Cluster = Cluster)  # obtain different parameters that can be parsed to the next step
}

# add the dataset for this problem; this step can either generate a static object or 
addProblem(name = "ReclusterMajorCluster", 
           data = Dir, 
           fun = subsample) 

# Define tasks
HarmoneyRecluster <- function(data, job, instance) { # instance is inherited from the last step
  data = instance$data
  Cluster = instance$Cluster
  
  ####### Harmony reclustering #######
  options(repr.plot.height = 2.5, repr.plot.width = 6)
  data <- FindVariableFeatures(data, 
                               selection.method = "vst", 
                               nfeatures = 3000,
                               verbose = F)
  data <- ScaleData(data, 
                    vars.to.regress = c("nFeature_RNA", "percent.mt", "decontX_contamination"), 
                    verbose = FALSE)
  
  data <- RunPCA(data, 
                 npcs = 30, 
                 pc.genes = data@var.genes,
                 verbose = TRUE)
  
  
  ## integration based on individual
  data <- data %>% RunHarmony("orig.ident", plot_convergence = TRUE)
  
  
  ## Identify clusters
  data <- data %>% 
    RunUMAP(reduction = "harmony", dims = 1:30) %>% 
    FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
    FindClusters(resolution = 1.0) %>% 
    identity()
  
  
  # Save the individual Seurat file
  # dir.create(paste0("./Step4.IterativeCleaning/", Cluster),
  #                       showWarnings = FALSE)
  SaveDir <- "./Step4.IterativeCleaning/Recluster.R1"
  
  saveRDS(data, paste0(SaveDir, "/", Cluster, ".rds"))
}

addAlgorithm(name = "HarmoneyRecluster", fun = HarmoneyRecluster)

# Creating jobs
# problem design: try two values for the ratio parameter
prob.design <- CJ(i = 1:7)

# final problem design
pdes <- list(ReclusterMajorCluster = prob.design)

# algorithm design: try combinations of kernel and epsilon exhaustively,
# try different number of trees for the forest
ades <- list(HarmoneyRecluster = data.table())

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

