# ********
# 3-de.r
# Code for conducting cellwise differential expression using the MAST model
source("Analysis/util.r")
source("Analysis/singlecell/differential-genes.R")
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
library(MAST)

# -------- load data  -----------
#x <- load_integrated()
DefaultAssay(x) <- 'RNA'
x <- NormalizeData(x)
x@meta.data %<>% mutate(cell.broadtype = paste0(celltype, ' ', SCT.seurat_snn_res.0.1),
  cell.broadtype = case_when(
    cell.broadtype == 'CN 2' ~ 'Intermediate cell 2',
    cell.broadtype == 'CN 6' ~ 'Photoreception 6',
    cell.subtype == 'AS 0' ~ 'AS 0',
    #stringr::str_detect(cell.broadtype, 'AS') ~ cell.subtype,
    TRUE ~ cell.broadtype
  )
) 
x@meta.data %<>% mutate(
  cell.broadtype = factor(cell.broadtype, 
                          c('CN 1', 
                            'Photoreception 6', 
                            'Inhibitory neuron 8', 
                            'Neuron (immature) 5', 
                            'Intermediate cell 0',
                            'Intermediate cell 2',
                            'NEC 3',                           
                            'AS 7', 
                            'AS 0',                           
                            'GPC 4',
                            'PGC 9'))
)
Idents(x) <- x$cell.broadtype
x$cngeneson <- scale(x$nFeature_RNA)
x$snca <- scale(x@assays$RNA@data['SNCA',])

# -------- Computing co-expression with SNCA gene-expression using MAST model ----------
snca.results <- list()
for (broadtype in unique(x$cell.broadtype)) {
  xs <- subset(x, cell.broadtype == broadtype)
  snca.results[[broadtype]] <- mast_lrt(
    xs,
    latent.vars = c('cngeneson', 'percent.mt', 'line.id', 'snca'),
    cpc = 0.005,
    min.pct = 0.005,
    doLRT = TRUE 
  )
}

gene.list = list()
for (broadtype in c('ExN 1', 'ExN 2')) {
  gene.list[[broadtype]] = snca.results[[broadtype]]$summaryCond$datatable %>% 
    filter(contrast == 'snca' & component == 'D' & z > 0 & primerid %in% rownames(subcluster.list[[broadtype]]@assays$SCT@scale.data)) %>% 
    arrange(`Pr(>Chisq)`, desc(coef)) %>% 
    head(40) %>% 
    select(primerid) %>% 
    unlist
}

out.dir = file.path('Results/Markers/snca.mast/cngeneson+percent.mt+snca+line.id')
dir.create(out.dir, recursive = TRUE)
for (broadtype in c('CN 1', 'Intermediate cell 2')) {
  dat = snca.results[[broadtype]]$summaryCond$datatable %>% filter(contrast == 'snca' & !is.na(z)) %>%
    arrange(`Pr(>Chisq)`,desc(coef)) %>% 
    as.data.table
  for (component in unique(dat$component)) {
    dat[dat$component == component, fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  }
  dat <- dat %>% filter(component %in% c('C', 'D'))
  write.table(dat,
              file = file.path(out.dir, paste('YJ SNCA', broadtype, 'snca coexpression.txt'),
              sep='\t', quote=F, row.names = F, col.names = T))
}

DefaultAssay(x) <- 'SCT'
p.list = list()
subcluster.list <- SplitObject(x, split.by='cell.broadtype')
for (broadtype in c('ExN 1', 'ExN 2')) {
  p.list = list()
  DefaultAssay(subcluster.list[[broadtype]]) <- 'SCT'
  for (g1 in gene.list[[broadtype]][1:10]) {
    for (g2 in gene.list[[broadtype]][1:10]) {
      p.list[[broadtype]][[paste0(g1, '.', g2)]] = FeatureScatter(subcluster.list[[broadtype]], feature1=g1, feature2=g2, group.by='line.id', slot='scale.data') + NoLegend() + NoAxes(keep.text = TRUE)
    }
  }
  cowplot::save_plot(plot=patchwork::wrap_plots(p.list[[broadtype]]),
                     filename = file.path(out.dir, paste('YJ SCNA', broadtype, 'coexpression scatterplots.png')),
                     base_asp=1.2,
                     base_height=30)
}

# ----- Perform cellwise differential expression ------
# Comparison of triplicate vs control groups using the MAST two-part hurdle DE model
de.results <- list()
for (broadtype in unique(x$cell.broadtype)) {
  xs <- subset(x, cell.broadtype == broadtype)
  de.results[[broadtype]] <- mast_lrt(
    xs,
    condition = 'Genotype',
    latent.vars = c('cngeneson', 'percent.mt', 'line.id'),
    cpc = 0.005,
    min.pct = 0.005,
    doLRT = TRUE
  )
}

for (broadtype in names(de.results)) {
  summaryDt = de.results[[broadtype]]$summaryCond$datatable
  cont <- "GenotypeSNCA trip"
  fcHurdle <- merge(
    summaryDt[contrast==cont & component=='H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
    summaryDt[contrast==cont & component=='logFC', .(primerid, coef, ci.hi, ci.lo, z)], by='primerid'
  )
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdle %<>% arrange(`Pr(>Chisq)`, desc(abs(coef)))
  write.table(fcHurdle, 
              file = paste0('Results/Results-for-resubmission/Markers/mast.de/YJ Genotype Trip vs Ctrl ', broadtype, '.txt'),
              sep='\t',
              col.names=T,
              row.names=F,
              quote=F)
}
