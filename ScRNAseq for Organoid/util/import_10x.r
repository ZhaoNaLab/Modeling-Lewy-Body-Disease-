# Given secondary path and meta information (a data frame)
# construct a merged SeuratObject.
create_merged_seurat_object <- function(secondary.path, meta=NULL, raw=FALSE, subsample=FALSE, ...) {
  source('Analysis/util/matrix.r')

  metric.summaries <- gather_metric_summary(secondary.path)
  x <- create_from_10x(secondary.path, raw=raw, subsample=subsample, ...) # subsampling here will make merge_metadata run much faster
  if (!is.null(meta)) {
    meta <- left_join(meta, metric.summaries, by = 'orig.ident')
  } else {
    meta <- metric.summaries
  }
  x <- merge_metadata(x = x, meta = meta)
  x <- subsample(x, subsample=subsample)
  x
}


create_from_secondary <- function(secondary.path, raw=FALSE, subsample=FALSE, ...) {
  source('Analysis/util/matrix.r') # subsample

  secondary.files <- list.files(secondary.path, full.names=FALSE)
  type <- ifelse(raw, 'raw', 'filtered')
  lapply(secondary.files,
         function(n) {
           counts <- Seurat::Read10X(paste0(secondary.path, '/', n, '/outs/', type, '_feature_bc_matrix/'), strip.suffix=TRUE)
           counts <- subsample(counts, subsample=subsample)
           x <- Seurat::CreateSeuratObject(counts = counts, 
                                           project = stringr::str_remove_all(n, 'out_'), 
                                           ...)
           x <- Seurat::RenameCells(x, add.cell.id = levels(x$orig.ident))
           return(x)
         })
}

# Provide a path to secondary analysis folders generated by 10X,
# each of which contains an 'outs' folder.
# I suppose that may change later -- right now we're on
# 10X v3. (I'm hoping this never changes ...)
create_from_10x <- function(secondary.path, raw=FALSE, subsample=FALSE, ...) {
  library(magrittr)
  source("Analysis/util/matrix.r") #subsample

  secondary.files <- list.files(secondary.path, full.names=FALSE)
  type <- ifelse(raw, 'raw', 'filtered')
  lapply(secondary.files,
         function(n) {
           counts <- Seurat::Read10X(paste0(secondary.path, '/', n, '/outs/', type, '_feature_bc_matrix/'), strip.suffix=TRUE)
           counts <- subsample(counts, subsample=subsample)
           x <- Seurat::CreateSeuratObject(counts = counts, 
                                           project = stringr::str_remove_all(n, 'out_'), 
                                           ...)
           x <- Seurat::RenameCells(x, add.cell.id = levels(x$orig.ident))
           return(x)
         })
}

custom_Read10X <- function(data.dir, subsample=FALSE) {
  source("Analysis/util/matrix.r") #subsample

  barcode.loc = list.files(data.dir, '*barcode.*.tsv', full.names=TRUE)
  gene.loc = list.files(data.dir, '*(gene|feature).*.tsv', full.names=TRUE)
  matrix.loc = list.files(data.dir, '*matrix.*.mtx', full.names=TRUE)
  
  stopifnot(length(barcode.loc) == 1)
  stopifnot(length(gene.loc) == 1)
  stopifnot(length(matrix.loc) == 1)

  data = Matrix::readMM(file = matrix.loc)
  cell.names <- readLines(barcode.loc)
  gene.names <- readLines(gene.loc)

  if (all(grepl(pattern = "\\-1$", x = cell.names))) {
    cell.names <- as.vector(x = as.character(x = sapply(X = cell.names, FUN = Seurat:::ExtractField, field = 1, delim = "-")))
  } 
  
  rownames(x = data) <- make.unique(names = as.character(x = sapply(X = gene.names, FUN = Seurat:::ExtractField, field = 2, delim = "\\t")))
  colnames(x = data) <- cell.names
  
  data <- subsample(data, subsample=subsample)
  return(data)
}

# x - a list of matrices imported from count data
# meta - meta information. this can be return value of import_secondary,
#        although you also typically want to include subject-level information.
#        by default, assume that there is an orig.ident column.
merge_metadata <- function(x, meta) {
  library(Seurat)
  if (! 'orig.ident' %in% colnames(meta)) {
    stop('orig.ident not found in columns in meta. Re-organize meta so that orig.ident is a column corresponding to the sample of origin.')
  }
  x.merged <- purrr::reduce(.x=x, .f=merge) #, collapse=TRUE) # 01/22/2024 a layer was being created for each orig.ident -- but collapsing layers not yet supported
  #x.merged <- Reduce(f = merge, x) # if we're running directly on create_from_10x, each cell id should already be unique.
  #x.merged <- RenameCells(x.merged, new.names = paste0(stringr::str_split_fixed(colnames(x.merged), '-', 2)[,1], '-', x.merged$orig.ident))
  if (R.version$minor >= 1) {
    x.merged <- JoinLayers(x.merged) # forms a single counts layer
  }
  x.merged@meta.data <- as.data.frame(dplyr::left_join(x.merged@meta.data, meta, by = 'orig.ident'), row.names=colnames(x.merged))
  x.merged$percent.mt <- PercentageFeatureSet(x.merged, pattern = '^(mt|MT)-')
  x.merged$percent.rp <- PercentageFeatureSet(x.merged, pattern = '^(Rp[s|l]|RP[S|L])')

  assay = DefaultAssay(x.merged)
  x.merged$cngeneson <- scale(log10(x.merged[[paste0('nFeature_', assay)]]))
  return(x.merged)
}

# Given secondary path, create a table of metric summaries for all outs
gather_metric_summary <- function(secondary.path) {
  library(dplyr)
  library(magrittr)
  source('Analysis/util/clean.R')
  # secondary path is a folder containing all your 10X results, which themselves are folders.
  # Example:
  # /path/to/secondary/results/
  #                           /<sample 1>/outs/...
  #                           /<sample 2>/outs/... 
  # etc.
  # Returns a table containing metric summaries for each sample, organized into rows.
  # See the column 'results.folder' to know which sample each row corresponds to.
  secondary.files <- list.files(path = secondary.path, full.names = TRUE)
  metric.summaries <- lapply(secondary.files, function(f) read.csv(paste0(f, '/outs/metrics_summary.csv'), header = T))
  metric.summaries <- Reduce(rbind, metric.summaries)
  secondary.files <- list.files(path=secondary.path)
  metric.summaries <- metric.summaries %>% 
    mutate(across(.cols = Estimated.Number.of.Cells : Number.of.Reads, ~ clean(.x))) %>%
    mutate(across(.cols = Total.Genes.Detected : Median.UMI.Counts.per.Cell, ~ clean(.x))) %>%
    mutate(across(.cols = Valid.Barcodes : Fraction.Reads.in.Cells, ~ stringr::str_remove(.x, '%') %>% as.numeric())) %>%
    mutate(across(.cols = Valid.Barcodes : Fraction.Reads.in.Cells, ~ .x / 100)) %>%
    mutate(orig.ident = stringr::str_remove_all(secondary.files, 'out_'))
  return(metric.summaries)
}