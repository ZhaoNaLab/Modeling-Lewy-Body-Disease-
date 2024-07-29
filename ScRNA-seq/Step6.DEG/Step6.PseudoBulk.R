library("MAST")
library("edgeR")


setwd("./Step6.DEG")

# load the SCA objects
scaObjs <- readRDS('scaObjs.rds')

for (i in 1:length(scaObjs)) {
  CellType <- names(scaObjs)[i]
  Sca = scaObjs[[i]]
  x = scuttle::summarizeAssayByGroup(x=assay(Sca), 
                                     ids = Sca@colData$orig.ident,
                                     statistics='sum')
  
  meta = Sca@colData %>% 
    as_tibble %>%
    distinct(orig.ident, .keep_all=TRUE) %>%
    arrange(orig.ident) %>%
    select(orig.ident, Sex, Age, APOE, RIN, batch, SeqSat, zAge, Dx, PathDx, TYPE)
  
  dge = DGEList(counts=assay(x), samples = meta, group=meta$Dx)
  plotMDS.DGEList(dge, plot = TRUE)
  
  keep = filterByExpr(dge, design = model.matrix(~ Dx + Age + RIN + SeqSat + as.factor(batch), data=meta))
  summary(keep)
  dge <- dge[keep,]
  dge <- calcNormFactors(dge)
  
  design = model.matrix(~ Dx + Age, data=meta) # Age + RIN + SeqSat + as.factor(batch)
  dge <- estimateDisp(dge, design=design)
  plotBCV(dge)
  summary(dge$trended.dispersion)
  fit <- glmQLFit(dge, design)
  plotQLDisp(fit)
  
  result = glmQLFTest(fit, coef=2)
  write.csv(x=topTags(result, n=Inf)$table %>%
              tibble::rownames_to_column('gene'),
            file=paste0('PseudoBulk.DEG/', CellType, '.pseudobulk DEG.csv'), quote=F)
}
