library("tidyverse")
library("data.table")

setwd("./Step6.DEG/DEGs")
FileList <- list.files(pattern = "Random.DEG.csv")
Celltypes <- str_split_fixed(FileList, "\\.", n = 4)[, 1]

DEGFileCombind <- list()

for (i in 1:length(Celltypes)){
  dat <- read.csv(FileList[i], row.names = 1)
  setDT(dat)[, fdr:= p.adjust(`Pr..Chisq.`, 'fdr'), by = "contrast"]
  dat$DEG[dat$fdr < 0.05 & dat$coef >= 0.1] <- "up"
  dat$DEG[dat$fdr < 0.05 & dat$coef <= -0.1] <- "down"
  dat <- arrange(dat, fdr)
  
  ## save each individual file
  write.csv(dat, paste0(Celltypes[i], ".DEG.updated.csv"), row.names = FALSE)
  # save the updated DEG files
  DEGFileCombind[[i]] <- dat
}

names(DEGFileCombind) <- Celltypes
saveRDS(DEGFileCombind, "DEGFileCombind.rds")





