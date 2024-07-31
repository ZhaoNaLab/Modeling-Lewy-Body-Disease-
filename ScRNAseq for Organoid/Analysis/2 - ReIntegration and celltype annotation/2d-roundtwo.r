# ************
# 2d-roundtwo.r
# Another round of integrating data, following removal of contaminating cell clusters
source("Analysis/util.r")
library(Seurat)
library(data.table)
library(magrittr)
library(dplyr)
library(ggplot2)

# -------------- Removal of contaming subclusters -----------------
#re <- readRDS('Results/Integration-RoundOne/roundone-integrated-list.rds')
# Identified the following contaminating subclusters:
contaminating = list(
  g0 = 5, 
  g2 = 9,
  g4 = 4,
  g6 = 5,
  g7 = 4,
  g8 = 4
)
resolutions = list(
  g0 = 0.1,
  g2 = 0.3,
  g4 = 0.1,
  g6 = 0.1,
  g7 = 0.4,
  g8 = 0.1
)
to_remove <- list(
  g0 = 5,
  g2 = 9,
  g4 = 4,
  g6 = 5,
  g7 = 6,
  g8 = 4
)
for (i in c(0,2,4,6,7,8)) {
  re[[paste0("g", i)]] <- RunPCA(re[[paste0("g", i)]]) %>%
    FindNeighbors(dims=1:40) %>%
    FindClusters(resolution=0.1, graph.name = c('SCT_nn', 'SCT_snn')) %>%
    RunUMAP(dims=dims, reduction='pca')
  clusters_to_keep <- unique(re[[paste0("g", i)]]$seurat_clusters)
  clusters_to_keep <- clusters_to_keep[! clusters_to_keep %in% to_remove[[paste0("g", i)]]]
  re <- subset(re, seurat_clusters %in% clusters_to_keep)
}

# Additional splitting of data
re[["g7.0"]] < subset(re[["g7"]], seurat_clusters == 0)
re[["g7.not0"]] < subset(re[["g7"]], seurat_clusters != 0)
re <- re[c(paste0("g", 0:6), "g7.0", "g7.not0", paste0("g", 8:10))]
names(re) <- c(
  'Intermediate cell 0',
  'CN 1',
  'Intermediate cell 2',
  'NEC 3',
  'GPC 4',
  'Neuron (immature) 5',
  'Photoreception 6',
  'AS 0',
  'AS 7',
  'Inhibitory neuron 8',
  'PGC 9' 
)

# Create "celltype" column.
for (n in names(re)) {
  re[[n]]@meta.data$celltype <- n
}
saveRDS(re, 'Data/Processed/YJ SNCA split seurat.rds')

# -------- Re-integration of round two -- final integration ---------
# re <- readRDS("Data/Processed/YJ SNCA split seurat.rds")
x <- purrr::reduce(.x = purrr::map(.x=re, ~{
  DefaultAssay(.x) <- 'RNA'
  .x <- DietSeurat(.x, assays="RNA")
  return(.x)
}), .f=merge)

# apply additional threshold filtering
x <- subset(x, percent.mt < 25 & nFeature_RNA > 250)

x@meta.data[['line.id']] <- factor(x=x@meta.data[['line.id']])
x <- SplitObject(object=x, split.by='line.id')
x <- lapply(x, function(y)  SCTransform(object=y, assay="RNA", vars.to.regress=c('percent.mt', 'percent.rp', 'nfeature_rna'), method='glmGamPoi', vst.flavor='v2'))
for (i in 1:length(x)) {
  levels(x[[i]][['SCT']]) <- names(x)[[i]]
}
# Update -- upon reviewing the code and figures in the manuscript, this last integration was not performed. 
# Although this might lead to some unwanted subclone to subclone heterogeneity,
# this should not significantly affect results, as our cluster labels are quite broad (see below) and the
# subclone to subclone heterogeneity remains even after integration.
#features <- SelectIntegrationFeatures(object.list=x, nfeatures=3000)
#x <- PrepSCTIntegration(object.list=x, anchor.features=features)
#anchors <- FindIntegrationAnchors(object.list=x, anchor.features=features, normalization.method="SCT", reference=1, dims=1:40, scale=FALSE)
#x <- IntegrateData(anchorset=anchors, normalization.method='SCT', dims=1:40)
x <- RunPCA(x) %>%
  FindNeighbors(dims=1:40) %>%
  FindClusters(resolution=0.1, graph.name = c('SCT.seurat_nn', 'SCT.seurat_snn')) %>%
  RunUMAP(dims=1:40, reduction='pca')
saveRDS(x, "Data/Roundtwo-Integration/alphasyn-triplication-final-roundtwo-integrated.rds")

###  ----------- Re-annotate clusters and make plots ------------
# The below steps are also saved inside the load_integrated() method inside Analysis/util.r
#x <- readRDS("Data/Roundtwo-Integration/alphasyn-triplication-final-roundtwo-integrated.rds")
x@meta.data %<>% 
  mutate(cell.broadtype = paste0(celltype, ' ', SCT.seurat_snn_res.0.1),
         cell.broadtype = case_when(
           cell.broadtype == 'CN 6' ~ 'ExN 3/Photo 6',
           cell.subtype == 'AS 0' ~ 'AS 0',
           #stringr::str_detect(cell.broadtype, 'AS') ~ cell.subtype,
           TRUE ~ cell.broadtype),
         cell.broadtype = case_when(
           cell.broadtype == 'Neuron (immature) 5' ~ 'Neuron (InN immature) 5',
           cell.broadtype == 'Intermediate cell 0' ~ 'Neuron (ExN immature) 0',
           cell.broadtype == 'NEC 3' ~ 'proRG 3',
           cell.broadtype == 'AS 7' ~ 'Astro 7',
           cell.broadtype == 'CN 1' ~ 'ExN 1',
           cell.broadtype == 'CN 2' ~ 'ExN 2',
           cell.broadtype == 'GPC 4' ~ 'RG 4',
           cell.broadtype == 'PGC 9' ~ 'RG 9',
           cell.broadtype == 'Inhibitory neuron 8' ~ 'InN 8',
           TRUE ~ cell.broadtype),
         cell.broadtype = factor(cell.broadtype, cell.broadtype.levels) 
         )
x <- subset(x, cell.broadtype != 'AS 0')
x$cell.broadtype <- droplevels(x$cell.broadtype)

x@meta.data %>%
  mutate(orig.ident = factor(orig.ident, paste0('BU-SNCA-', 1:24))) %>%
  reshape2::dcast(orig.ident + flow.cell + Lane + collection.date + line.id + Genotype + Sex ~ cell.broadtype) %>%
  writexl::write_xlsx('YJ celltype cell counts.xlsx')

p <- DimPlot(x, group.by='cell.broadtype', raster=F)
cowplot::save_plot(plot=p, filename='YJ combined umap.png',
                   base_asp=1.8, base_height=5)
p <- DimPlot(x, group.by='cell.broadtype', raster=F, label=TRUE)
cowplot::save_plot(plot=p, filename='YJ combined labelled umap.png',
                   base_asp=1.8, base_height=5)
x$orig.ident = factor(x$orig.ident, paste0('BU-SNCA-', 1:24))
p <- DimPlot(x, group.by='cell.broadtype', split.by='orig.ident', ncol=6, raster=F) + NoAxes()
cowplot::save_plot(plot=p, filename='YJ umap split.png',
                   base_asp=1.3*6/4 + .1, base_height=10)
DefaultAssay(x) <- 'SCT'
p <- VlnPlot(x, 'APOE', group.by='cell.broadtype', split.by='Genotype')
cowplot::save_plot(plot=p, filename='YJ APOE violin.png',
                   base_asp=1.6, base_height=5)
p <- VlnPlot(x, 'SNCA', group.by='cell.broadtype', split.by='Genotype')
cowplot::save_plot(plot=p, filename='YJ SNCA violin.png',
                   base_asp=1.6, base_height=5)
p <- VlnPlot(x, 'MAPT', group.by='cell.broadtype', split.by='Genotype')
cowplot::save_plot(plot=p, filename='YJ MAPT violin.png',
                   base_asp=1.6, base_height=5)
p <- FeaturePlot(x, "APOE", slot='data', cols=c('gray95', 'red'), split.by='Genotype')
cowplot::save_plot(plot=p, filename='YJ APOE umap.png', base_asp=1.2*2, base_height=5)
p <- FeaturePlot(x, "SNCA", slot='data', cols=c('gray95', 'red'), split.by='Genotype')
cowplot::save_plot(plot=p, filename='YJ SNCA umap.png', base_asp=1.2*2, base_height=5)
p <- FeaturePlot(x, "MAPT", slot='data', cols=c('gray95', 'red'), split.by='Genotype')
cowplot::save_plot(plot=p, filename='YJ MAPT umap.png', base_asp=1.2*2, base_height=5)

genes = c('SNCA', 'HERC5', 
          #'PIGY', 
          'HERC3', 'NAP1L5', 
          #'FAM13AOS', 
          'FAM13A', 
          'TIGD2', 'GPRIN3', 'MMRN1') 
          #'LOC644248', 
          #'PYURI')
p = VlnPlot(x, genes, group.by='cell.broadtype', split.by='Genotype', ncol=2)
cowplot::save_plot(plot=p, filename='YJ triplication genes.png',
                   base_asp=1.2, 
                   base_height=20)

p <- DimPlot(x, group.by='orig.ident')
cowplot::save_plot(
  plot = p,
  filename = 'Results/Processed/umap by sample.pdf',
  base_height = 5,
  base_asp = 1.8
)
p <- DimPlot(x, group.by='line.id')
cowplot::save_plot(
  plot = p,
  filename = 'Results/Processed/umap by line id.pdf',
  base_height = 5,
  base_asp = 1.8
)

p <- DimPlot(x, group.by='cell.broadtype', split.by='line.id', ncol=2)
cowplot::save_plot(
  plot = p,
  filename = 'Results/Processed/umap split by line id.pdf',
  base_height = 10,
  base_asp = 1.3
)

source('Analysis/plot/plot.r')
plot_composition_barchart(x, cluster.by='cell.broadtype', group.by='line.id', file.dir = 'Results/Processed')