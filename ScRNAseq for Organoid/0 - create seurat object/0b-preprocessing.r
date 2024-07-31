# ***************
# 0b-preprocessing.r
# Initial preprocessing of the data.
source('Analysis/util.r')
library(Seurat)
library(dplyr)
library(magrittr)

x@meta.data <- x@meta.data
  group_by(orig.ident) %>%
  mutate(keep = ifelse(nCount_RNA < quantile(nCount_RNA, 0.99) & nFeature_RNA < quantile(nFeature_RNA, 0.99), TRUE, FALSE)) %>% 
  data.frame(row.names=colnames(x))
x.filt <- subset(x, cells=colnames(x)[x$keep])
x.filt$keep <- NULL

x.filt@meta.data %>% 
  ggplot(aes(x=nFeature_RNA, y=nCount_RNA, color=percent.mt)) +
  scale_color_viridis_c() + 
  theme_classic() +
  geom_point(alpha=0.3) +
  facet_wrap(~orig.ident)

saveRDS(x.filt, 'Data/Processed/alphasyn-triplication-initial-filtered-seurat.rds')

x.filt <- subset(x.filt, orig.ident != 'BU-SNCA-13')
saveRDS(x.filt, "Data/Processed/alphasyn-triplication-initial-filtered-seurat-wo13.rds")