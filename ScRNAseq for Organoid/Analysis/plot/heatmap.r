# These functions contain some default arguments that worked
# in the past for creating heat maps.
# It also helps with some of the preprocessing.

library(ComplexHeatmap)
library(dplyr)
library(tibble)
library(magrittr)
library(RColorBrewer)
library(circlize)


# Input matrix pretty much needs to be exactly the heatmap you want to plot.
# If you want labels, they need to be in colnames or rownames, or added as
# column_split, row_split arguments.
make_heatmap <- function(matrix,
                         name,
                         colors = c('dodgerblue3', '#EEEEEE', 'firebrick2'), # Pretty colors, but maybe there's something prettier...
                         col=colorRamp2(seq(min(matrix), max(matrix), length=3), colors),
                         show_row_dend = F,
                         show_column_dend = F,
                         rect_gp = gpar(col='white', lwd=2),
                         width = ncol(matrix) * unit(5, 'mm'),
                         height = nrow(matrix) * unit(5, 'mm'),
                         heatmap_legend_param = list(
                           title = 'Normalized gene expression',
                           direction = 'vertical',
                           legend_height = unit(3, 'cm')
                         ),
                         ...
                         ) {
  h <- Heatmap(as.matrix(matrix),
               name = name,
               col = col,
               show_row_dend = show_row_dend,
               show_column_dend = show_column_dend,
               rect_gp = rect_gp,
               width = width,
               height = height,
               heatmap_legend_param = heatmap_legend_param,
               ...
               )
  h
}

draw_heatmap <- function(h, heatmap_legend_side = 'right', ...) {
  grid.grabExpr(draw(h, heatmap_legend_side = heatmap_legend_side, ...))
}
