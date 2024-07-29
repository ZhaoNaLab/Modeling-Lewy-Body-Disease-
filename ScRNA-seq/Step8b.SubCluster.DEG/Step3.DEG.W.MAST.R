# https://mllg.github.io/batchtools/articles/batchtools.html
library("batchtools")
setwd("./Step8b.SubCluster.DEG")
reg <- makeExperimentRegistry("DEGRegistry",
                              packages = c("MAST", "data.table", "Seurat", "lme4"),
                              seed = 22112019)

reg <- loadRegistry("DEGRegistry", writeable = TRUE)


# Step1: function to subset data
# sca file has been filtered and were siplit into chunks
scaObjs <- readRDS('sca.SubClusters.Objs.rds')

# define instance (data preparation)
subsample <- function(data, job, i) {
  Cluster <- names(data)[i]
  data <- data[[i]]
  list(data = data, Cluster = Cluster)  # obtain different parameters that can be parsed to the next step
}

# add the dataset for this problem; this step can either generate a static object or 
addProblem(name = "DEGCal", data = scaObjs, fun = subsample) 

# create an algorithm which applies a support vector machine
MASTDEG <- function(data, job, instance) {           # instance is inherited from the last step
  data = instance$data
  Cluster = instance$Cluster
  data@colData$Dx <- "LBD"
  data@colData$Dx[data@colData$PathDx == "Normal"] <- "Ctrl"
  data@colData$Dx <- relevel(as.factor(data@colData$Dx), ref = "Ctrl")
  zlmCond <- zlm(~ Dx + zAge + cngeneson + percent.mt + decontX_contamination + SeqSat + RIN + (1|orig.ident),
                 method = 'glmer',
                 ebayes = FALSE,
                 fitArgsD = list(nAGQ = 0),
                 sca = data)
  
  summaryCond <- summary(zlmCond, doLRT= c('DxLBD'))
  
  # print out the data.table
  summaryDt <- summaryCond$datatable
  fcHurdle <- merge(summaryDt[component=='H',.(primerid, `Pr(>Chisq)`, contrast)], #hurdle P values
                    summaryDt[component=='logFC', .(primerid, coef, ci.hi, ci.lo, contrast, z)], by=c('primerid', 'contrast')) #logFC coefficients
  
  
  write.csv(fcHurdle, paste0('DEG.random/', Cluster, '.Random.DEG.csv'))
}

addAlgorithm(name = "MASTDEG", fun = MASTDEG)

# Creating jobs
# problem design: try two values for the ratio parameter
prob.design <- CJ(i = 1:50)

# final problem design
pdes <- list(DEGCal = prob.design)

# algorithm design: try combinations of kernel and epsilon exhaustively,
# try different number of trees for the forest
ades <- list(MASTDEG =data.table())

addExperiments(pdes, ades) # pdes set parameter to the data part; ades set parameters to the instance parts

# summarize the experiment
summarizeExperiments()


# Submitting and Collecting Results
options(batchtools.progress = FALSE)
submitJobs()
getStatus()
getErrorMessages()
getJobTable()
# clearRegistry(reg = getDefaultRegistry()
# killJobs(ids = c(1:7))
# removeExperiments(ids = c(9, 12, 19, 20))

