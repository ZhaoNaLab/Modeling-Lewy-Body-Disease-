FindAllMarkers_for_MAST <- function(
  object,
  assay = NULL,
  features = NULL,
  slot = 'data',
  logfc.threshold = 0.25,
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = c('cngeneson'),
  min.cells.feature = 3,
  min.cells.group = 3,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  workers=4,
  return.model=FALSE,
  ...
  )
{
  source('Analysis/util/future.r')

  library(furrr)
  setup_future(object, 'multicore', min(workers, length(levels(object))))
  results <- future_map(levels(Idents(object)), 
    ~FindMarkers_for_MAST(
      object=object,
      ident.1=.x,
      assay=assay,
      latent.vars=latent.vars,
      features=features,
      slot=slot,
      logfc.threshold=logfc.threshold,
      min.pct=min.pct,
      min.diff.pct = min.diff.pct,
      verbose=verbose,
      only.pos=only.pos,
      max.cells.per.ident=max.cells.per.ident,
      random.seed=random.seed,
      min.cells.feature=min.cells.feature,
      min.cells.group=min.cells.group,
      mean.fxn=mean.fxn,
      fc.name=fc.name,
      base=base,
      return.model=return.model),
    ...)
  names(results) <- levels(object)
  return(results)
}

#' FindMarkers, implemented to return MAST summary datatable output
#' alongside the usual FindMarkers output.
#' It only works on one of the slots of an assay, so not
#' on reductions.
#' One additional parameter is return.model (Defult: FALSE)
#' If set to TRUE, then the zlm fit and summary objects will be returned.
#' @export
FindMarkers_for_MAST = function(
  object,
  ident.1 = NULL, 
  ident.2 = NULL,
  return.model=FALSE, # By default, return a table of marker genes
  group.by = NULL,
  subset.ident = NULL,
  assay=NULL,
  slot='data',
  features = NULL,
  logfc.threshold = 0.25,
#test.use = 'MAST',
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  only.pos = FALSE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = 'cngeneson',
  random.effect.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  mean.fxn = NULL,
  fc.name = NULL,
  base = 2,
  ...)
{
  library(MAST)
  
  assay <- assay %||% DefaultAssay(object = object)
  if (!is.null(x = group.by)  ) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }

  if (is.null(ident.2)) {
    ident.2 <- levels(object)[levels(object) != ident.1]
  }

  # select which data to use
  #data.use <- GetAssayData(object = object, slot = slot)
  data.use <- LayerData(object=object, layer=slot)

  #Idents(object) <- object$parental.or.isogenic
  cells <- Seurat:::IdentsToCells(
    object = object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    cellnames.use = colnames(x=data.use)
  )
  
  # fetch latent.vars
  # since we're using MAST, include cngeneson as latent.vars, always
  
  # Honestly not sure if we're to scale nFeature_RNA or nCount_RNA -- 
  # these are more or less the same thing, numerically
  #object$cngeneson <- scale(log(object[[paste0('nFeature_', object@active.assay)]]))  # -- 2024-02-22: again I disagree that we set these kinds of variables within the method. It's opaque. I'd much rather error.
  if (!is.null(x = latent.vars)) {
    ## There isn't any issue with latent.vars having spaces in the names here --
    # but such an issue does arise later.
    # For now, I'm applying a temporary patch where I replace all spaces in latent.vars with underscores, and 
    # then log that here.

    if (!all(latent.vars == stringr::str_replace_all(latent.vars, ' ', '_'))) {
      print('Some latent variables have spaces -- replacing with underscores.')
      latent.vars = stringr::str_replace_all(latent.vars, ' ', '_')
      colnames(object@meta.data) = stringr::str_replace_all(colnames(object@meta.data), ' ', '_')
    }
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(cells$cells.1, cells$cells.2)
    )
  }

  if (!is.null(x = random.effect.vars)) {
    random.effect.vars <- FetchData(
      object = object,
      vars = random.effect.vars,
      cells = c(cells$cells.1, cells$cells.2)
    )
  }

  cells.1 = cells$cells.1
  cells.2 = cells$cells.2
  
  fc.results <- Seurat:::FoldChange.Assay(
    object = object,
    slot = slot,
    cells.1 = cells.1,
    cells.2 = cells.2,
    features = features,
    pseudocount.use = pseudocount.use,
    mean.fxn = mean.fxn,
    fc.name = fc.name,
    base = base
  )
  Seurat:::ValidateCellGroups(
    object = object,
    cells.1 = cells.1,
    cells.2 = cells.2,
    min.cells.group = min.cells.group
  )
  features <- features %||% rownames(x = object)
  # reset parameters so no feature filtering is performed
  #if (test.use %in% DEmethods_noprefilter()) {
  #  features <- rownames(x = object)
  #  min.diff.pct <- -Inf
  #  logfc.threshold <- 0
  #}

  # feature selection (based on percentages)
  alpha.min <- pmax(fc.results$pct.1, fc.results$pct.2)
  names(x = alpha.min) <- rownames(x = fc.results)
  features <- names(x = which(x = alpha.min >= min.pct))
  if (length(x = features) == 0) {
    warning("No features pass min.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
  alpha.diff <- alpha.min - pmin(fc.results$pct.1, fc.results$pct.2)
  features <- names(
    x = which(x = alpha.min >= min.pct & alpha.diff >= min.diff.pct)
  )
  if (length(x = features) == 0) {
    warning("No features pass min.diff.pct threshold; returning empty data.frame")
    return(fc.results[features, ])
  }
  # feature selection (based on logFC)
  if (slot != "scale.data") {
    total.diff <- fc.results[, 1] #first column is logFC
    names(total.diff) <- rownames(fc.results)
    features.diff <- if (only.pos) {
      names(x = which(x = total.diff >= logfc.threshold))
    } else {
      names(x = which(x = abs(x = total.diff) >= logfc.threshold))
    }
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      warning("No features pass logfc.threshold threshold; returning empty data.frame")
      return(fc.results[features, ])
    }
  }
  # subsample cell groups if they are too large
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2), , drop = FALSE]
    }
  }

  #############################################
  # Perform MAST
  
  #de.results <- PerformDE(
  #  object = object,
  #  cells.1 = cells.1,
  #  cells.2 = cells.2,
  #  features = features,
  #  test.use = test.use,
  #  verbose = verbose,
  #  min.cells.feature = min.cells.feature,
  #  latent.vars = latent.vars,
  #  ...
  #)
  #de.results <- cbind(de.results, fc.results[rownames(x = de.results), , drop = FALSE])
 
  data.use = data.use[features, c(cells.1, cells.2), drop = FALSE]
  
  # Check for MAST
  #if (!Seurat:::PackageCheck('MAST', error = FALSE)) {
  #  stop("Please install MAST - learn more at https://github.com/RGLab/MAST")
  #}
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  latent.vars <- latent.vars %||% group.info
  group.info[cells.1, "group"] <- ident.1
  ident.2 = ifelse(length(ident.2) > 1, paste0('not_', ident.1), ident.2)
  group.info[cells.2, "group"] <- ident.2
  group.info[, "group"] <- factor(x = group.info[, "group"])
  latent.vars.names <- c("group", colnames(x = latent.vars))
  if (!is.null(random.effect.vars))  {
    random.effect.vars.names <- colnames(x = random.effect.vars)
  }
  latent.vars <- cbind(latent.vars, group.info)
  latent.vars$wellKey <- rownames(x = latent.vars)
  fdat <- data.frame(rownames(x = data.use))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  if (!is.null(random.effect.vars)) {
    sca <- MAST::FromMatrix(
      exprsArray = as.matrix(x = data.use),
      check_sanity = FALSE,
      cData = cbind(latent.vars, random.effect.vars),
      fData = fdat
    )
  } else {
    sca <- MAST::FromMatrix(
      exprsArray=as.matrix(x = data.use),
      check_sanity=FALSE,
      cData = latent.vars,
      fData=fdat
    )
  }
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = ident.2)
  SummarizedExperiment::colData(sca)$group <- cond
  if (!is.null(random.effect.vars)) {
    fmla <- as.formula( object = paste0(" ~ ", 
                                        paste(c(latent.vars.names, 
                                                paste(paste0('(1|', random.effect.vars.names,')'), collapse='+')),
                                              collapse = "+")
                                        )
    )
  } else {
    fmla <- as.formula( object = paste0(" ~ ", paste(latent.vars.names, collapse='+')))
  }

  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  summaryCond <- MAST::summary(object = zlmCond, doLRT = colnames(zlmCond@coefC)[2])
  summaryDt <- summaryCond$datatable
  cont <- levels(summaryDt$contrast)[1]
  fcHurdle <- merge(
    summaryDt[contrast==cont & component=='H', .(primerid, `Pr(>Chisq)`)], #hurdle P values
    summaryDt[contrast==cont & component=='logFC', .(primerid, coef, ci.hi, ci.lo, z)], by='primerid'
  ) #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  de.results <- cbind(fcHurdle, fc.results[fcHurdle$primerid, , drop = F])
  colnames(de.results)[2] <- 'p_val' 
  
  if (only.pos) {
    de.results <- de.results[de.results$coef > 0, , drop = FALSE]
  }
  #if (test.use %in% DEmethods_nocorrect()) {
  #  de.results <- de.results[order(-de.results$power, -de.results[, 1]), ]
  #} else {
  de.results <- de.results[order(de.results$p_val, decreasing=FALSE), ] 
  de.results$p_val_adj = p.adjust(
    p = de.results$p_val,
    method = "bonferroni",
    n = nrow(x = object)
  )
  #}
  de.results <- tibble::column_to_rownames(.data = as.data.frame(de.results), var = 'primerid')
  if (return.model) {
    return(list(zlmCond=zlmCond, summaryCond=summaryCond, de.results=de.results))
  } else {
    return(de.results)
  }
}

#' mast_lrt
# cpc : counts per cell -- genes with counts per cell smaller than this specified value will not be analyzed. cpc is defined as the ratio between the total counts per gene and number of cells. 
# min.pct : genes that are expressed in less than the minimum percentage in all levels of condition will not be analyzed.
# condition : grouping variable for setting up contrasts
# contrasts : for specifying differential expression - use levels in condition
# latent.vars : covariates to also model, including factor variables
# random.effect.vars : variables for modeling random effect. In general, we model as (1|A) + (1|B) + ... for re variables c('A', 'B', ...)
#' @export
mast_lrt = function(
  object,
  cpc = 0.0,                     
  min.pct=0.1, 
  condition=NULL,
  contrasts=NULL,
  formula=NULL,
  latent.vars=NULL,
  random.effect.vars=NULL,
  return.table.only=FALSE,
  doLRT=FALSE,
  logFC=TRUE,
  return.model=FALSE,
  ...)
{
  library(SingleCellExperiment)
  library(future.apply)
  library(sparseMatrixStats)
  library(scuttle)
  library(Seurat)
  library(MAST)
  plan(sequential)  # set to multicore if multiple contrasts

  if (is.null(condition)) {
    stopifnot("Latent vars must be specified if there is no condition." = !is.null(latent.vars))
    stopifnot("Contrasts specified without condition." = is.null(contrasts))
  }
    
  # Convert to SingleCellExperiment for scuttle later
  if (class(object) == 'Seurat') {
    if (all(dim(GetAssayData(object, 'data')) == c(0,0))) {
        DefaultAssay(object) <- 'RNA'
        object <- NormalizeData(object) # normalize for tpm counts
    }
    sce <- as.SingleCellExperiment(object, assay='RNA')
  } else if (class(object) == 'SingleCellAssay') {
    sce <- SingleCellExperiment(assays(object), colData = colData(object))
  }
  stopifnot(class(sce) == 'SingleCellExperiment')

  if (cpc > 0) {
    sce <- sce[sparseMatrixStats::rowSums2(assay(sce, 'counts')) / ncol(sce) > cpc, ]
  }
  # Create aggregate proportions 
  if (!is.null(condition)) {
    summed <- scuttle::summarizeAssayByGroup(sce, ids = sce[[condition]], statistics = 'prop.detected')
      
    # Keep only the genes that are expressed in min.pct proportion of cells
    # in at least one of the groups
    keep <- rowSums(assay(summed, 'prop.detected') > min.pct) > 0 
    sce <- sce[keep,]
    print(paste('Keeping', nrow(sce), 'genes ...'))
  } else {
    summed = scuttle::summarizeAssayByGroup(sce, ids=rep('cpc', ncol(sce)), statistics = 'prop.detected')
  }
                                        
  sca <- MAST::FromMatrix(exprsArray = list(logcounts=as.matrix(assay(sce, 'logcounts'))),
                          cData = colData(sce))

  if (is.null(formula)) {
    if (!is.null(condition)) {
      if (!is.null(random.effect.vars)) {
          fmla <- formula(paste0("~", condition, '+', paste0(c(latent.vars, paste0(paste0('(1|', random.effect.vars,')'), collapse='+')), collapse = "+")
      ))
      } else if (!is.null(latent.vars)) {
        fmla <- formula(paste0("~", condition, '+', paste0(latent.vars, collapse='+')))
      } else {
        fmla <- formula(paste0("~", condition))
      }
    } else {
      if (!is.null(random.effect.vars)) {
          fmla <- formula(paste0("~", paste0(c(latent.vars, paste0(paste0('(1|', random.effect.vars,')'), collapse='+')), collapse = "+")
      ))
      } else {
        fmla <- formula(paste0("~", paste0(latent.vars, collapse='+')))
      }
    }
  } else {
    fmla <- as.formula(formula)
  }
  # By default, obtain residuals, if random.effect.vars is NULL
  if (is.null(random.effect.vars)) {
    print('Obtaining combined residuals, check zlmCond@hookOut')
    zlmCond <- MAST::zlm(formula = fmla, sca = sca, hook=combined_residuals_hook, ...)
  } else {
    zlmCond <- MAST::zlm(formula = fmla, sca = sca, method='glmer', ebayes=FALSE, strictConvergence=FALSE, ...) 
  }

  if (!is.null(contrasts)) {
    coefnames <- colnames(coef(zlmCond, 'D'))
    #design = model.matrix(fmla, data=colData(sca))
    lfc <- lapply(contrasts, function(cont) {
      contrast1 = limma::makeContrasts(contrasts=cont, levels=coefnames)
      lfc <- logFC(zlmCond, contrast1=contrast1)
    })

    #plan(multisession, workers = min(length(contrasts), 4))
    #options(future.globals.maxSize = object.size(zlmCond) * 1.1)
    lrt.results <- lapply(contrasts, function(cont) { lrTest(zlmCond, hypothesis=limma::makeContrasts(contrasts=cont, levels=coefnames)) }) #Hypothesis(cont, terms=colnames(design)))})
    if (return.table.only) {
      summaryList <- purrr::map2(.x=lrt.results, .y=contrasts, .f= ~{
          lrt <- data.table::as.data.table(.x)
          lrt <- lrt[test.type == 'hurdle' & metric == 'Pr(>Chisq)', .(primerid, p_val = value)]
          lrt[,fdr:=p.adjust(p_val, 'fdr')]
          lrt[, contrast := .y]
          return(lrt[order(p_val)])
        })
      summaryDt <- purrr::map2(lfc, summaryList, ~{cbind(logFC=.x$logFC[.y$primerid,], .y) })
      summaryDt <- dplyr::bind_rows(summaryDt)
      return(summaryDt)
    } else {
      summaryCond <- summary(zlmCond, doLRT=doLRT, logFC=logFC)
      if (return.model) {
        return(list(lrt.results = lrt.results, prop.detected = assay(summed), zlmCond=zlmCond, summaryCond=summaryCond, lfc=lfc))
      } else {
        return(list(lrt.results=lrt.results, prop.detected=assay(summed), summaryCond=summaryCond, lfc=lfc))
      }
    }
  } else {
    summaryCond <- summary(zlmCond, doLRT=doLRT, logFC=logFC)
    if (return.model) {
      return(list(prop.detected = assay(summed), zlmCond=zlmCond, summaryCond=summaryCond))
    } else {
      return(list(summaryCond=summaryCond, prop.detected=assay(summed)))
    }
  }
}