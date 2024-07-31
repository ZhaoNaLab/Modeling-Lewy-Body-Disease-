#### DoubletFinder ####
run_doubletfinder = function(x, labels=NULL, ncells=ncol(x), data=NA, PCs=1:30, sct=FALSE, GT=FALSE, num.cores=6, pN=0.25, return.metadata=TRUE, ...) {
  # Detect doublets with DoubletFinder
  # Recommended job submission:
  # doubletfinder = CJ(
  #  sct=TRUE,
  #    num.cores=1,
  #      labels='seurat_clusters'
  #  )
  # submitJobs(resources=list(h_vmem='30G', queue='1-hour'))
  stopifnot(!is.null(labels))
  ## pK identification (no ground-truth)
  
  library(DoubletFinder)
  
  sweep.res.list <- paramSweep_v3(x, PCs=PCs, sct=sct, num.cores=num.cores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=GT)
  bcmvn <- find.pK(sweep.stats)
  print(bcmvn)
  pK = bcmvn[which.max(bcmvn$BCmetric),'pK']
  print(pK)
  
  ## Homotypic Doublet Proportion Estimate
  annotations <- x@meta.data[, labels]
  homotypic.prop <- modelHomotypic(annotations)

  # estimated multiplet rate for v3.1 10X chemistry
  # https://kb.10xgenomics.com/hc/en-us/articles/360001378811-What-is-the-maximum-number-of-cells-that-can-be-profiled-
  # https://support.10xgenomics.com/single-cell-gene-expression/library-prep/doc/user-guide-chromium-single-cell-3-reagent-kits-user-guide-v31-chemistry

  dat <- data.frame(
    multiplet_rate=c(0.4, 0.8, 1.6, 2.3, 3.1, 3.9, 4.6, 5.4, 6.1, 6.9, 7.6) * 0.01,
    number_cells_loaded=c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
    number_cells_recovered=c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
  )

  # create a linear model from the above table, assuming the number of cells corresponds to number of cells recovered.
  # should work alright even if we extrapolate
  fit <- lm(multiplet_rate ~ number_cells_recovered, data=dat)
  a = as.numeric(fit$coefficients[1])
  b = as.numeric(fit$coefficients[2])

  print(paste('ncells:', ncells))
  tenX_doublet_prop = a + b*ncells
  print(paste('tenX_doublet_prop:', tenX_doublet_prop))
  nExp_poi <- round(tenX_doublet_prop*nrow(x@meta.data)) ## Assuming the above table of doublet proportions for v3.1
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  ## Run DoubletFinder with varying classification stringencies
  x <- doubletFinder_v3(x, PCs=PCs, pN=pN, pK=pK, nExp=nExp_poi, reuse.pANN=FALSE, sct=sct)
  pANN = colnames(x@meta.data)[stringr::str_detect(colnames(x@meta.data), 'pANN')]
  x <- doubletFinder_v3(x, PCs=PCs, pN=pN, pK=pK, nExp=nExp_poi.adj, reuse.pANN=pANN, sct=sct)

  if (return.metadata) {
    return(x@meta.data)
  } else {
    return(x)
  }
}

# Clean up doublet Finder results.
clean_doubletfinder = function(x) {
  library(stringr)
  pANN = colnames(x@meta.data)[str_detect(colnames(x@meta.data), 'pANN')]
  x$DF.pANN <- x@meta.data[, pANN]
  DF.classifications = colnames(x@meta.data)[str_detect(colnames(x@meta.data), 'DF.classification')] 
  if (length(DF.classifications) == 2) {
    x$DF.classification.upper <- x@meta.data[, DF.classifications[1]]
    x$DF.classification.lower <- x@meta.data[, DF.classifications[2]]
  } else {
    x$DF.classification.upper <- x@meta.data[, DF.classifications]
  }
  x@meta.data <- x@meta.data[, !colnames(x@meta.data) %in% c(pANN, DF.classifications)]
  return(x)
}

# Filter doublet finder results.
filter_doubletfinder = function(re=NULL, meta.data=NULL, group.by = 'decontX_clusters', filter.mt = TRUE, filter.df = TRUE, filter.ct = TRUE, percent.mt.qt = 0.9, percent.mt.th = 50, doublet.prop=0.4, ct.score.qt = 0.9, ct.score.th = 0.01) {
  # Filtering is applied cluster-wise.
  # For filtering mitochondrial expressing cluster, we compute the quantile of percent mitochondria
  # at percent.mt.qt quantile, then filter a cluster if this value is greater than percent.mt.th threshold.
  # For doublet cluster, we remove if greater than doublet.prop proportion of the cluster is doublet.
  # For contaminating cluster, we compute the quantile of the contamination score at ct.qt, then
  # filter a cluster if this value is greater than ct.score.th
 
  if (is.null(re) & is.null(meta.data)) {
    stop('One but not both of re and meta.data must be provided.')
  } 
  stopifnot(!is.null(re) | !is.null(meta.data)) 
  if (!is.null(meta.data)) {
    x = meta.data
  } else {
    x = re@meta.data
  }
  
  barcodes = rownames(x)
  mt.barcodes = ct.barcodes = df.barcodes = c()
  
  # remove low quality clusters
  if (filter.mt) {
    mt.barcodes = x %>% 
      rownames_to_column('barcode') %>%
      group_by(!!sym(group.by)) %>% 
      summarise(mt = quantile(percent.mt, percent.mt.qt), barcode) %>% 
      filter(mt > percent.mt.th) %>% 
      ungroup() %>%
      select(barcode) %>% unlist %>% unname
  }
  
  # remove contaminated clusters
  if (filter.ct) {
    ct.barcodes = x %>% 
      rownames_to_column('barcode') %>%
      group_by(!!sym(group.by)) %>% 
      summarise(ct = quantile(decontX_contamination, ct.score.qt), barcode) %>% 
      filter(ct > ct.score.th) %>% 
      ungroup %>%
      select(barcode) %>% unlist %>% unname
  }
    
  # remove doublet clusters
  if (filter.df) {
    df.barcodes <- x %>% 
      rownames_to_column('barcode') %>%
      group_by(!!sym(group.by), DF.classification.upper, orig.ident) %>% 
      group_by(!!sym(group.by), orig.ident) %>% 
      summarise(
        ncells = n(), 
        pDoublet = sum(DF.classification.upper == 'Doublet') / sum(ncells),
        barcode
      ) %>%
      filter(pDoublet > doublet.prop) %>% 
      ungroup %>%
      select(barcode) %>% unlist %>% unname
  }
  
  barcodes_to_keep = !barcodes %in% c(mt.barcodes, ct.barcodes, df.barcodes)
  summary(barcodes_to_keep)
  if (sum(barcodes_to_keep) < 0.9 * length(barcodes)) {
    warnings('barcodes_to_keep is less than 90% of all barcodes')
  }
    
  if (!is.null(re)) {
    re <- subset(re, cells=barcodes_to_keep)
    return(re)
  } else {
    meta.data = meta.data[barcodes_to_keep,]
    return(meta.data)
  }
}

