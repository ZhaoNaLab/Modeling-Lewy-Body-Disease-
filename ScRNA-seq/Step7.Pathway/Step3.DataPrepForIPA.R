library("tidyverse")
library("gprofiler2")
library("data.table")
library("BiocParallel")


# load the DEG datasets
setwd("./Projects/SNCSynuclein")

# load the DEG files
DEGFileCombind <- readRDS("Step6.DEG/DEGs/DEGFileCombind.rds")
Cluser <- names(DEGFileCombind)
DEGdf <- list()

for (i in 1:length(DEGFileCombind)){
  dat <- DEGFileCombind[[i]] %>% as.data.frame() %>%
    select(primerid, coef, Pr..Chisq., fdr) 
  names(dat) <- unlist(c("Symbol", paste0(Cluser[i], c(".logFC", ".Pvalue", ".FDR"))))
  DEGdf[[i]] <- dat
}

# Combine all data
DegDfComb <- purrr::reduce(DEGdf, full_join, by = "Symbol")
DegDfComb[is.na(DegDfComb)] <- ""

# Save the data
write.csv(DegDfComb, "Step7.Pathway/DEGDataForIPA.csv", row.names = FALSE)


