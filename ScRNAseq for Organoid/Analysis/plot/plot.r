# ---- Utility functions for making nice file names ----
prefix_filename <- function(file.name, prefix=NULL) {
  if (!is.null(prefix)) {
    file.name <- paste(prefix, file.name)
  }
  return(file.name)
}

groupby_splitby_filename <- function(file.name, group.by=NULL, split.by=NULL) {
  if (!is.null(group.by)) {
    file.name = paste(file.name, 'grouped by', group.by)
  }
  if (!is.null(split.by)) {
    file.name = paste(file.name, 'split by', split.by)
  }
  return(file.name)
}

# ---- Plotting functions -----
plot_featureplots <- function(
  object,
  features,
  min.cutoff = 'q5',
  max.cutoff = 'q95',
  cols = c('gray90', 'red'),
  order = TRUE,
  split.by = NULL,
  assay = NULL,
  combine = FALSE,
  prefix = NULL,
  file.dir = NULL,
  file.name = 'feature plot',
  base_height = 5,
  base_asp = 1.3,
  raster = TRUE,
  ...
) {
  assay <- assay %||% DefaultAssay(object = object)
  DefaultAssay(object) <- assay
  p <- FeaturePlot(
    object,
    features, 
    cols = cols,
    order = order,
    min.cutoff = min.cutoff,
    max.cutoff = max.cutoff,
    split.by = split.by,
    combine = combine,
    raster = raster,
    ...
  )
  if (is.null(file.dir)) {
    return(p)
  } else {
    file.name = file.name  %>%
      groupby_splitby_filename(split.by=split.by) %>%
      prefix_filename(prefix)
    #pdf(file.path(file.dir, paste0(file.name, '.pdf')))
    #p
    #dev.off()
    if (!combine) {
      N = length(p)
      dims = ggplot2::wrap_dims(N)
      #num_col = ceiling(sqrt(N))
      #num_row = ceiling(N / num_col)
      cowplot::save_plot(
        plot = patchwork::wrap_plots(p, ncol = dims[2]),
        filename = file.path(file.dir, paste0(file.name, '.pdf')),
        base_asp = base_asp * dims[2] / dims[1],
        base_height = base_height * dims[1],
        limitsize = FALSE
      )
    } else {
      cowplot::save_plot(
        plot = p,
        filename = file.path(file.dir, paste0(file.name, '.pdf')),
        base_asp = base_asp,
        base_height = base_height,
        limitsize = FALSE
      )
    }
  }
}

plot_umap <- function(
  object,
  group.by = NULL,
  split.by = NULL,
  prefix = NULL,
  file.dir = NULL,
  file.name = 'umap',
  base_height = 5,
  base_asp = 1.3,
  raster=FALSE,
  ncol=NULL,
  nrow=NULL,
  ...
) {
  if (!is.null(split.by)) {
    N = dplyr::n_distinct(object@meta.data[[split.by]])
    if (is.null(ncol) | is.null(nrow)) {
      dims = ggplot2::wrap_dims(N)
      nrow = dims[1]
      ncol = dims[2]
    }
  } else {
    nrow = ncol = NULL
  }


  p <- DimPlot(object, group.by=group.by, split.by=split.by, ncol = ncol, raster=raster, ...)

  if (is.null(file.dir)) {
    return(p)
  } else {
    file.name <- file.name %>%
      groupby_splitby_filename(group.by, split.by) %>%
      prefix_filename(prefix)
    if (is.null(split.by)) {
      #pdf(file.path(file.dir, paste0(file.name, '.pdf')))
      #p
      #dev.off()
      cowplot::save_plot(
        plot = p,
        filename = file.path(file.dir, paste0(file.name, '.pdf')),
        base_asp = base_asp,
        base_height = base_height
      )
    } else {
      cowplot::save_plot(
        plot = p,
        filename = file.path(file.dir, paste0(file.name, '.pdf')),
        base_asp = base_asp*ncol/nrow,
        base_height = base_height*nrow
      )
    }
  }
}

plot_violins <- function(
  object,
  features,
  group.by = NULL,
  split.by = NULL,
  prefix = NULL,
  file.dir = NULL,
  file.name = 'violins',
  base_height = 3.71,
  base_asp = 1.618, 
  #base_asp = 1.618 + 0.1*dplyr::n_distinct(object@meta.data[[group.by]]),
  stack = FALSE,
  ...
) {
  if (is.null(group.by)) {
    group.by = 'orig.ident'
  }
  num_groups = dplyr::n_distinct(object@meta.data[[group.by]])
  p <- VlnPlot(object, features=features, group.by=group.by, split.by=split.by, stack=stack, ...)
  if (is.null(file.dir)) {
    return(p)
  } else {
    file.name <- file.name %>%
      groupby_splitby_filename(group.by, split.by) %>%
      prefix_filename(prefix)
 
    # NOTE 08/16/2023: adding 0.2 to base_asp to make the qc plots with 5 qc metrics stacked, grouped by cluster, look pretty.
    if (stack) {
      cowplot::save_plot(
        plot = p,
        filename = file.path(file.dir, paste0(file.name, '.pdf')),
        limitsize = FALSE,
        base_height = base_height * length(features),
        base_asp = 0.2 + max(base_asp, num_groups * 0.2) / length(features)
      )
    } else {
      dims = ggplot2::wrap_dims(length(features))
      cowplot::save_plot(
        plot = p,
        filename = file.path(file.dir, paste0(file.name, '.pdf')),
        limitsize = FALSE,
        base_height = base_height * dims[1],
        base_asp = 0.2 + max(base_asp, num_groups*0.2) * dims[2] / dims[1]
      )
    }
  }
}

plot_composition_barchart <- function(
  object,
  cluster.by,  # this is the thing we want to visualize compositions of
  group.by = 'orig.ident',
  split.by = NULL,
  colors = NULL,
  file.dir = NULL,
  base_asp = 1.5,
  base_height = 5,
  prefix = NULL,
  file.name = 'composition barchart',
  seed = 42,
  ...
) {
  set.seed(seed)
  if (is.null(colors)) {
    #colors = sample(grDevices::colorRampPalette(colors = RColorBrewer::brewer.pal(n=12, name='Set3'))(length(unique(object@meta.data[[cluster.by]]))))
    colors = gg_color_hue(length(unique(object@meta.data[[cluster.by]])))
  }
  if (is.null(split.by)) {
    p = object@meta.data %>% 
      #mutate(group.id = paste0(.data[[group.by]], '.', orig.ident)) %>%
      #ggplot(aes(x=group.id, fill=.data[[cluster.by]])) + 
      ggplot(aes(x=.data[[group.by]], fill=.data[[cluster.by]])) + 
      geom_bar(position='fill', color='black') + 
      scale_fill_manual(values = colors) + 
      theme_classic() + 
      theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
      xlab(group.by) + 
      ylab(paste(cluster.by, 'proportion'))

    if (!is.null(file.dir)) {
      file.name = file.name %>%
        groupby_splitby_filename(group.by=group.by, split.by=split.by) %>%
        prefix_filename(prefix = prefix)
      cowplot::save_plot(
        plot=p,
        filename = file.path(file.dir, paste0(file.name, '.pdf')),
        base_height = base_height + 0.1 * dplyr::n_distinct(object@meta.data[[cluster.by]]),
        base_asp = base_asp * dplyr::n_distinct(object@meta.data[[group.by]]) / dplyr::n_distinct(object@meta.data[[cluster.by]]),
        limitsize = FALSE
        #base_asp=1.2*2,
    #    base_height=10
      )
    }
  } else {
    p = object@meta.data %>% 
#      mutate(group.id = paste0(.data[[group.by]], '.', orig.ident)) %>%
#      ggplot(aes(x=orig.ident, fill=.data[[cluster.by]])) + 
      ggplot(aes(x = .data[[group.by]], fill = .data[[cluster.by]])) + 
      geom_bar(position='fill', color='black') + 
      scale_fill_manual(values = colors) + 
      theme_classic() + 
      theme(axis.text.x=element_text(angle=90, hjust=0.5, vjust=0.5)) +
      xlab(group.by) + 
      ylab(paste(cluster.by, 'proportion')) + 
      facet_wrap(split.by, scales='free_x')
    
    if (!is.null(file.dir)) {
      dims = ggplot2::wrap_dims(dplyr::n_distinct(object@meta.data[[split.by]]))
      file.name = file.name %>%
        groupby_splitby_filename(group.by=group.by, split.by=split.by) %>%
        prefix_filename(prefix = prefix)
      cowplot::save_plot(
        plot = p,
        filename = file.path(file.dir, paste0(file.name, '.pdf')),
        base_asp = base_asp * dims[2] / dims[1],
        base_height = (base_height + 0.1*dplyr::n_distinct(object@meta.data[[cluster.by]])) * dims[1],
        limitsize = FALSE
      )
    }
  }
  if (is.null(file.dir)) { return(p) }
}

plot_dotplot <- function(
  object,
  features,
  group.by,
  split.by=NULL,
  cols='RdYlBu',
  cluster.idents=TRUE,
  prefix = NULL,
  file.dir = NULL,
  file.name = 'dotplot',
  base_height = 3.71,
  base_asp = 1.618, 
  ...
) {
  features <- unique(features)
  p <- DotPlot(object, features=features, group.by=group.by, split.by=split.by, cols=cols, cluster.idents=cluster.idents)  +
    coord_flip() + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  if (!is.null(file.dir)) {
    file.name <- prefix_filename(file.name, prefix=prefix) %>%
      groupby_splitby_filename(group.by=group.by, split.by=split.by)
    cowplot::save_plot(
      plot = p,
      filename = file.path(file.dir, paste0(file.name, '.pdf')),
      base_height = base_height + 0.25*length(features),
      base_asp = base_asp
    )
  }
  return(p)
}


plot_metadata_histogram <- function(
  object,
  feature,
  group.by = NULL,
  split.by = NULL,
  scales = 'free_y',
  prefix=NULL,
  file.dir = NULL,
  base_asp = 1.5,
  base_height = 5,
  file.name = 'histogram'
) {
  
  if (! feature %in% colnames(object@meta.data)) {
    rlog::log_error(paste(feature, 'is not in object metadata.'))
    stop()
  }
  p <- ggplot(data = object@meta.data, mapping=aes(x = .data[[feature]]))
  if (is.null(group.by)) {
    p <- p + geom_histogram()
  } else {
    p <- p + geom_histogram(mapping = aes(fill=.data[[group.by]]))
  }
  if (!is.null(split.by)) {
    p <- p + facet_wrap(split.by, scales=scales)
    N = length(unique(object@meta.data[[feature]]))
    dims = ggplot2::wrap_dims(N)
  }

  if (!is.null(file.dir)) {
    file.name <- file.name %>%
      prefix_filename(prefix=prefix) %>%
      groupby_splitby_filename(group.by=group.by, split.by=split.by)

    cowplot::save_plot(
      plot = p,
      filename = file.path(file.dir, paste0(file.name, '.pdf')),
      base_height = base_height * ifelse(is.null(split.by), 1, dims[1]),
      base_asp = base_asp * ifelse(is.null(split.by), 1, dims[2] / dims[1])
    )
  }
  return(p)
}


# ---- Other plotting utilities
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' plot_volcano_template
# this is not an actual function, it's just a template for making nice volcano plots
plot_volcano_template <- function(toptable) {
  EnhancedVolcano::EnhancedVolcano(
    toptable,
    x = 'logFC',
    y = 'PValue', 
    pCutoffCol = 'FDR',
    lab = toptable$gene,
    selectLab = toptable$gene[toptable$FDR < 0.4],
    drawConnectors = TRUE,
    min.segment.length = 0.5,
    widthConnectors=0.1,
    labSize=2,
    title = '',
    subtitle='',
    caption=''
  )
}
