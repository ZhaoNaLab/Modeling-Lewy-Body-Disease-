library("tidyverse")
library("data.table")

setwd("./Step6.DEG/DEGs")

# Load the combined DEG files
DEGFileCombind <- readRDS("DEGFileCombind.rds")

DEG.tables <- list()
for (i in 1:length(DEGFileCombind)){
  df <- DEGFileCombind[[i]]
  deg.sum <- setDT(df)[, .(DEG.number = .N), by = "DEG"]
  deg.sum <- drop_na(deg.sum)
  deg.sum$CellType <- names(DEGFileCombind)[i]
  DEG.tables[[names(DEGFileCombind)[i]]] <- deg.sum
}

DEG.table.combined <- do.call(rbind, DEG.tables) %>% dplyr::select(CellType, DEG, DEG.number)
write.csv(DEG.table.combined, "./Step6.DEG/DEG.number.table.csv", row.names = FALSE)



