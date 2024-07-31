source('Analysis/util.r')
x <- load_integrated()
Idents(x) <- x$cell.broadtype
p1 = FeaturePlot(x, 'PRKN', cols=c('gray95', 'red'), min.cutoff = 'q5', max.cutoff='q95', split.by='Genotype', label=TRUE)
p2 = VlnPlot(x, 'PRKN', assay='SCT', group.by='cell.broadtype', split.by='Genotype', pt.size=0)
cowplot::save_plot(plot=p1, filename='Summary/Figures/PRKN umap genotype split.png', base_asp=1.2*2.2, base_height=5)
cowplot::save_plot(plot=p2, filename='Summary/Figures/PRKN violin genotype split.png', base_asp=3, base_height=5)

#load_integrated_metadata() %>% dplyr::count(orig.ident, flow.cell, Lane, collection.date, Genotype, Sex, cell.broadtype) %>% tidyr::pivot_wider(values_from=n, names_from=cell.broadtype, values_fill = 0)
yj_counts = readxl::read_excel('Results/Processed/YJ celltype cell counts.xlsx', range='T29:V39', sheet = 'proportions (%) re-named')
colnames(yj_counts)[1] = 'celltype'
p = yj_counts %>%
  tidyr::pivot_longer(-celltype, names_to='Genotype', values_to='prop') %>% 
  mutate(celltype = factor(celltype, celltype.levels),
         Genotype = factor(Genotype, Genotype.levels)) %>%
  ggplot(aes(y = Genotype, x=prop, fill=fct_rev(celltype))) + 
  geom_col(color='black') + 
  scale_fill_brewer(palette='Set3', direction=-1) + 
  theme_classic() + 
  theme(text=element_text(size=8), legend.position='bottom') + 
  scale_x_continuous(expand = expansion(), breaks = seq(0,100,20)) + 
  guides(fill = guide_legend(reverse=TRUE, title = 'cell type'))
cowplot::save_plot(plot=p, filename='Results/Processed/YJ celltype composition by genotype.png', base_asp=4, base_height=2)
