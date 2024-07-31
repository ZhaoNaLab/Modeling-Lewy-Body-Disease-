# ************* 
# 1b-celltype-annotation.r
# Note: these codes are for exploratory data analysis -- you can ignore this script.

source("Analysis/util.r")
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(data.table)
file.dir = 'Results/Initial/Initial integration/'

# Load result of integrating by subclone.
x <- readRDS("Data/Processed/alphasyn-triplication-initial-integrated.rds")
x$orig.ident <- factor(x$orig.ident, paste0('BU-SNCA-', 1:24))
p <- DimPlot(x, split.by='orig.ident', raster=F, ncol=6) + NoAxes()
cowplot::save_plot(plot = p, filename=file.path(file.dir, 'Na YJ SNCA integration by sub-clone.png', base_asp=1.2*6/4, base_height=10))

p <- DimPlot(x, split.by='orig.ident', raster=F, ncol=6) + NoAxes()
cowplot::save_plot(plot = p, filename=file.path(file.dir, 'Na YJ SNCA integration by sub-clone res 0.1 wo-13.png'), base_asp=1.2*6/4, base_height=10)

p <- DimPlot(x, raster=F, pt.size=0.1, label=T, label.size=4)
cowplot::save_plot(plot=p, filename='Results/Initial/Initial integration/snca-org-trip-umap-wo13-labelled.png', base_asp=1.2, base_height=5)

DefaultAssay(x) <- 'SCT'
p <- FeaturePlot(x, features='SNCA', raster=F, pt.size=0.1, split.by='Genotype', slot = 'data', cols=c('gray95', 'red'))
cowplot::save_plot(plot = p, filename = 'Results/Initial/Initial integration/snca-org-trip-SNCA-feature-umap.png', base_asp=2.4, base_height=5)

p <- DimPlot(x, raster=F, pt.size=0.1, split.by='Genotype')
cowplot::save_plot(plot = p, filename = 'Results/Initial/Initial integration/snca-org-trip-umap-by-genotype-wo13.png', base_asp=2.4, base_height=5)

x$group <- paste0(x$Sex, re$Genotype)
x$id <- paste0(x$line.id, ' ', re$Genotype, '\n', re$orig.ident)
p <- DimPlot(x, split.by='id', ncol=6, raster=F, pt.size = 0.1)
cowplot::save_plot(plot=p, filename='Results/Initial/Initial integration/snca-org-trip-umap-split-wo13.png', base_asp=1.2 * 1.5, base_height=12)

p <- DimPlot(x, split.by='collection.date', ncol=3, raster=F, pt.size=0.1)
cowplot::save_plot(plot=p, filename='Results/Initial/Initial integration/snca-org-trip-umap-split-by-date-wo13.png', base_asp=1.2*1.5, base_height=10)

p <- DimPlot(x, split.by='group', ncol=3, raster=F, pt.size=0.1)
cowplot::save_plot(plot=p, filename='Results/Initial/Initial integration/snca-org-trip-umap-split-by-sex-genotype-wo13.png', base_asp=1.2*3, base_height=5)

x$line.genotype <- paste0(x$line.id, ' ', re$Genotype)
p <- DimPlot(x, split.by='line.genotype', ncol=3, raster=F, pt.size=0.1) + NoAxes()
cowplot::save_plot(plot=p, filename='Results/Initial/Initial integration/snca-org-trip-umap-split-by-line-wo13.png', base_asp=1.2, base_height=10)

p <- DimPlot(x, split.by='Lane', ncol=2, raster=F, pt.size=0.1) + NoAxes()
cowplot::save_plot(plot=p, filename='Results/Initial/Initial integration/snca-org-trip-umap-split-by-Lane-wo13.png', base_asp=1.2, base_height=10)

x@meta.data %>% distinct(orig.ident, .keep_all=T) %>% arrange(orig.ident) %>% tibble::rownames_to_column('ID') %>% select(-ID, -id, -group) %>% write.table(file.path(file.dir, 'snca-org-trip-metrics.txt'), sep='\t', row.names=F,col.names=T,quote=F)

avg <- AverageExpression(x, assays='SCT', slot='data', group.by = c('seurat_clusters', 'line.genotype'))$SCT
write.table(avg, file.path(file.dir, 'snca-org-trip-avg-expression-by-line-genotype.txt'), sep='\t', col.names=T, row.names=T, quote=F)

x@meta.data %>% reshape2::dcast(orig.ident + flow.cell + Lane + collection.date + line.id + Genotype + Sex ~ seurat_clusters) %>% write.table(file.path(file.dir, 'snca-org-trip-cluster-cell-counts.txt'), sep='\t', col.names=T, row.names=F, quote=F)

# Initial annotation of clusters
x@meta.data <- x@meta.data %>% 
  mutate(
    cell_type = case_when(
      seurat_clusters %in% c(9, 10, 29, 1, 0, 15, 3, 4, 14, 2, 8, 11, 16, 7, 21, 26, 28, 13, 20, 25) ~ 'EN',
      seurat_clusters == 12 ~ 'IP',
      seurat_clusters %in% c(5, 18, 6, 27) ~ 'RG',
      seurat_clusters %in% c(22, 23, 17) ~ 'AS',
      seurat_clusters %in% c(19, 24) ~ 'DV'
    )
  ) %>%
  mutate(
    cell_type = factor(cell_type, c('EN', 'RG', 'AS', 'IP', 'DV'))
  )

p <- x@meta.data %>% 
  filter(!is.na(cell_type)) %>% mutate(
    cell_type = case_when(
      cell_type == 'EN' ~ 'EN (72.14%)',
      cell_type == 'RG' ~ 'RG (15.98%)',
      cell_type == 'AS' ~ 'AS (5.49%)',
      cell_type == 'IP' ~ 'IP (3.22%)',
      cell_type == 'DV' ~ 'DV (3.17%)'
    )
  ) %>% mutate(
    cell_type = factor(cell_type, c('EN (72.14%)', 'RG (15.98%)', 'AS (5.49%)', 'IP (3.22%)', 'DV (3.17%)'))
  ) %>%
  rename(`Cell type` = cell_type) %>%
  ggplot(aes(x=UMAP_1, y=UMAP_2, color=`Cell type`)) +
  geom_point(size=0.1) +
  theme_classic() + 
  guides(colour = guide_legend(override.aes = list(size=3))) +
  scale_color_brewer(type='qual', palette='Set1')
cowplot::save_plot(plot = p, filename = 'Results/Initial/Initial integration/snca-org-trip-umap-wo13-celltype.png', base_asp=1.2, base_height=5)

p <- x@meta.data %>% 
  filter(!is.na(cell_type)) %>% 
  group_by(cell_type) %>% 
  mutate(Genotype = ifelse(Genotype == 'Ctrl', 'Ctrl', 'SNCA-T')) %>%
  rename(`Cell type` = cell_type) %>%
  ggplot(aes(x=Genotype, fill = `Cell type`)) + geom_bar(position='fill') + 
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  scale_fill_brewer(type='qual', palette='Set1') +
  ylab('Cell distribution (%)') + 
  geom_segment(aes(x=1,y=1.02,xend=2,yend=1.02)) + 
  geom_text(aes(x=1.5, y=1.05, label='N.S.'))
cowplot::save_plot(plot = p, filename='Results/Initial/Initial integration/snca-org-trip-cell-distribution.png', base_asp=1, base_height=5)
