# Find Markers for round two clusters
source('Analysis/util.r')
library(Seurat)
library(presto)
library(magrittr)
library(dplyr)
library(ggplot2)

# x <- load_integrated()
DefaultAssay(x) <- 'SCT'
Idents(x) <- x$cell.broadtype
x$cngeneson <- scale(x$nFeature_RNA)
x <- PrepSCTFindMarkers(x)
sct.markers <- FindAllMarkers(x, slot='scale.data', assay='SCT', only.pos = TRUE, min.pct=0.1)
sct.markers <- sct.markers %>% mutate(
  cell.broadtype = as.character(cluster),
  cell.broadtype = case_when(
  cell.broadtype == 'CN 3/Photo 6' ~ 'ExN 3/Photo 6',
  cell.broadtype == 'Neuron 5' ~ 'Neuron (InN immature) 5',
  cell.broadtype == 'Neuron 0' ~ 'Neuron (ExN immature) 0',
  cell.broadtype == 'DV 3' ~ 'proRG 3',
  cell.broadtype == 'AS 7' ~ 'Astro 7',
  cell.broadtype == 'CN 1' ~ 'ExN 1',
  cell.broadtype == 'CN 2' ~ 'ExN 2',
  cell.broadtype == 'GPC 4' ~ 'RG 4',
  cell.broadtype == 'PGC 9' ~ 'RG 9',
  cell.broadtype == 'Inhibitory neuron 8' ~ 'InN 8',
  TRUE ~ cell.broadtype
  ),
  cell.broadtype = factor(cell.broadtype, c('Neuron (ExN immature) 0', 'Neuron (InN immature) 5', 'ExN 1', 'ExN 2', 'ExN 3/Photo 6', 'InN 8', 'proRG 3', 'RG 4', 'RG 9', 'Astro 7')))
write.table(sct.markers, file='Results/Markers/YJ wilcox broadtype SCT markers.txt', sep='\t', row.names=T, quote=F, col.names=T)

DefaultAssay(x) <- 'RNA'
rna.markers <- FindAllMarkers(x, slot='data', assay='RNA', only.pos = TRUE, min.pct=0.1)
rna.markers <- rna.markers %>% mutate(
  cell.broadtype = as.character(cluster),
  cell.broadtype = case_when(
  cell.broadtype == 'CN 3/Photo 6' ~ 'ExN 3/Photo 6',
  cell.broadtype == 'Neuron 5' ~ 'Neuron (InN immature) 5',
  cell.broadtype == 'Neuron 0' ~ 'Neuron (ExN immature) 0',
  cell.broadtype == 'DV 3' ~ 'proRG 3',
  cell.broadtype == 'AS 7' ~ 'Astro 7',
  cell.broadtype == 'CN 1' ~ 'ExN 1',
  cell.broadtype == 'CN 2' ~ 'ExN 2',
  cell.broadtype == 'GPC 4' ~ 'RG 4',
  cell.broadtype == 'PGC 9' ~ 'RG 9',
  cell.broadtype == 'Inhibitory neuron 8' ~ 'InN 8',
  TRUE ~ cell.broadtype
  ),
  cell.broadtype = factor(cell.broadtype, c('Neuron (ExN immature) 0', 'Neuron (InN immature) 5', 'ExN 1', 'ExN 2', 'ExN 3/Photo 6', 'InN 8', 'proRG 3', 'RG 4', 'RG 9', 'Astro 7')))
write.table(rna.markers, file='Results/Markers/YJ wilcox broadtype RNA markers.txt', sep='\t', row.names=T, quote=F, col.names=T)

DefaultAssay(x) <- 'SCT'
genes.to.plot = sct.markers %>% 
  filter(gene %in% rownames(x@assays$SCT@scale.data)) %>% 
  group_by(cell.broadtype) %>% 
  slice_head(n=10) %>%
  ungroup %>%
  select(gene) %>%
  unlist

p <- DoHeatmap(x, features=genes.to.plot, group.by='cell.broadtype')
cowplot::save_plot(plot=p, filename='Results/Markers/YJ top10 heatmap.png',
                   base_asp=1.8,
                   base_height=15)

#### Plot cell type marker genes of our choice ####
f <- function(x) {
  DefaultAssay(x) <- 'SCT'
  p.list = FeaturePlot(x, cols=c('gray96', 'red'), features=c('STMN2', 'GAP43', 'SNCA', 'DCX', 'VIM', 'HES1', 'SOX2', 'TBR1', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'EOMES', 'TOP2A', 'MKI67', 'APOE', 'S100B', 'SLC1A3', 'MNS1', 'NPHP1', 'BMP4', 'MSX1', 'MYL1', 'MYH3', 'BGN', 'DCN', 'DDIT3'), combine = FALSE, max.cutoff = 'q95')
  return(p.list)
}
g <- function(x) {
  DefaultAssay(x) <- 'SCT'
  #x@active.ident = x$SCT.seurat_snn_res.0.1
  p = VlnPlot(x, features=c('STMN2', 'GAP43', 'SNCA', 'DCX', 'VIM', 'HES1', 'SOX2', 'TBR1', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'EOMES', 'TOP2A', 'MKI67', 'APOE', 'S100B', 'SLC1A3', 'MNS1', 'NPHP1', 'BMP4', 'MSX1', 'MYL1', 'MYH3', 'BGN', 'DCN', 'DDIT3'), ncol = 5, pt.size=0)
  return(p)
}
p.list = f(x)
cowplot::save_plot(plot = wrap_plots(p.list, ncol=5), filename = 'Results/Markers/YJ cell type marker umap.png', base_height=15, base_asp=1.2)
Idents(x) <- x$celltype
p = g(x)
cowplot::save_plot(plot = p, filename = 'Results/Markers/YJ cell type marker violin-plots.png', base_height=20, base_asp=1.2)


genes = c('STMN2', 'GAP43', 'DCX', 'VIM', 'HES1', 'SOX2', 'TBR1', 'SLC17A7', 'GAD1', 'GAD2', 'SLC32A1', 'EOMES', 'TOP2A', 'MKI67', 'S100B', 'SLC1A3', 'DCN', 'DDIT3')
p = DotPlot(x, assay='SCT', features=genes, group.by='cell.broadtype') + coord_flip()
p = p + theme_classic() + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
cowplot::save_plot(plot=p, filename='Results/Markers/YJ cell type dot plot.png',
                   base_asp=1.2, base_height=8)
