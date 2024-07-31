# ********
# 3-pseudode.r
# Code for running pseudobulk DEG analysis.
# In addition, we filter DEG's with the results of 3-de.r, MAST model DEG results.

source("Analysis/util.r")
library(scuttle)
library(edgeR)
library(Seurat)
library(dplyr)
library(magrittr)
library(ggplot2)
#x <- load_integrated()
DefaultAssay(x) <- 'RNA'
x <- NormalizeData(x)

# Create pseudobulk profile
for (broadtype in unique(x$cell.broadtype)) {
  xs <- subset(x, cell.broadtype == broadtype)
  xs = scuttle::summarizeAssayByGroup(x=xs@assays$RNA@counts, 
                                      ids = xs$orig.ident,
                                      statistics='sum')
  meta = xs@meta.data %>% 
    as_tibble %>%
    distinct(orig.ident, .keep_all=TRUE) %>%
    arrange(orig.ident) %>%
    select(orig.ident, Genotype, Line, line.id, flow.cell)
  
  dge = DGEList(counts=assay(xs), samples = meta, group=meta$Genotype)
  keep = filterByExpr(dge, design=model.matrix(~Genotype + flow.cell, data=meta)) # TODO: Fit subclone instead of Genotype
  summary(keep)
  dge <- dge[keep,]
  dge <- calcNormFactors(dge)
  
  design = model.matrix(~Genotype + flow.cell, data=meta)
  dge <- estimateDisp(dge, design=design, robust=TRUE)
  #plotBCV(dge)
  summary(dge$trended.dispersion)
  fit <- glmQLFit(dge, design) 
  #plotQLDisp(fit)
  
  result = glmQLFTest(fit, coef=2)
  write.csv(x=topTags(result, n=3e4)$table %>%
                tibble::rownames_to_column('gene'),
            file=paste0('Results/Results-for-resubmission/PseudobulkResults/', 'YJ ', stringr::str_replace(broadtype, '\\/', '.'), ' pseudobulk DEG.csv'), quote=F)
}

################################################################################
# Gene filtering
p.list = list()
labels = tribble(
  ~pseudo.de, ~cell.de,
  'CN.1', 'CN 1',
  'CN.2', 'Intermediate cell 2',
  'Neuron.0', 'Intermediate cell 0',
  'DV.3', 'NEC 3',
  'GPC.4', 'GPC 4',
  'Neuron.5', 'Neuron (immature) 5',
  'CN.3.Photo 6', 'Photoreception 6',
  'AS.7', 'AS 7',
  'Inhibitory.neuron 8', 'Inhibitory neuron 8',
  'PGC.9', 'PGC 9'
)
for (i in 1:nrow(labels)) {
  celltype = labels[i,'pseudo.de'][[1]]
  pseudo.de = read.csv(paste('Results/PseudobulkResults/YJ', labels[i,'pseudo.de'], 'pseudobulk DEG.csv'))
  cellwise.de =  read.table(paste0('Results/Results-for-resubmission/Markers/mast.de/YJ Genotype Trip vs Ctrl ', labels[i,'cell.de'], '.txt'), sep='\t', header = T)
  colnames(cellwise.de)[1] = 'gene'
  combined <- pseudo.de %>% left_join(cellwise.de, by='gene')
    
  p.list[[celltype]] <- combined %>%
    tidyr::drop_na() %>%
    mutate(
      DEG = ifelse(FDR < 0.05 & fdr < 0.05 & sign(logFC) == sign(coef) & abs(coef) > 0.05 & abs(logFC) > 0.05, 'yes', 'no')) %>%
    ggplot(aes(x=coef, y=logFC, fill=DEG, alpha=ifelse(DEG == 'yes', 1, 0.001))) + 
    geom_point(shape=21, show.legend = FALSE) +
    scale_fill_manual(values=c('black', 'red')) + 
    labs(title=celltype) + 
    xlab('cellwise comparison') +
    ylab('pseudobulk comparison') + 
    facet_wrap(~DEG)
  write.csv(combined %>% mutate(
      DEG = ifelse(FDR < 0.05 & fdr < 0.05 & sign(logFC) == sign(coef) & abs(coef) > 0.05 & abs(logFC) > 0.05, 'yes', 'no')) %>% filter(DEG == 'yes'), file=paste('Results/Results-for-resubmission/Filtered-DEG-Tables/YJ', celltype, 'filtered deg table.csv'), row.names = F, quote = F)
  write.csv(combined %>% mutate(
      DEG = ifelse(FDR < 0.05 & fdr < 0.05 & sign(logFC) == sign(coef) & abs(coef) > 0.05 & abs(logFC) > 0.05, 'yes', 'no')), file=paste('Results/Results-for-resubmission/Filtered-DEG-Tables-including-nonDEGs/YJ', celltype, 'filtered deg table with non-degs.csv'), row.names = F, quote = F)
}

dir.create("Results/Results-for-resubmission/DEGs", recursive = TRUE)
purrr::walk2(.x=p.list, .y=names(p.list),
      .f=~cowplot::save_plot(plot=.x, filename=paste0('Results/Results-for-resubmission/DEGs/DEG filter ', .y, '.png'), base_asp=1.5, base_height=8))
