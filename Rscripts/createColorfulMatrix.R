# Create a colourful matrix for figure abstract.

# Load libraries
library(ComplexHeatmap)
library(circlize)

# Create random matrix


create_random_matrix <- function(seed=123){
  set.seed(seed)
  nr1 = 2; nr2 = 3; nr3 = 2; nr = nr1 + nr2 + nr3
  nc1 = 2; nc2 = 3; nc3 = 2; nc = nc1 + nc2 + nc3
  mat = cbind(rbind(matrix(rnorm(nr1*nc1, mean = 1,   sd = 0.5), nr = nr1),
                    matrix(rnorm(nr2*nc1, mean = 0,   sd = 0.5), nr = nr2),
                    matrix(rnorm(nr3*nc1, mean = 0,   sd = 0.5), nr = nr3)),
              rbind(matrix(rnorm(nr1*nc2, mean = 0,   sd = 0.5), nr = nr1),
                    matrix(rnorm(nr2*nc2, mean = 1,   sd = 0.5), nr = nr2),
                    matrix(rnorm(nr3*nc2, mean = 0,   sd = 0.5), nr = nr3)),
              rbind(matrix(rnorm(nr1*nc3, mean = 0.5, sd = 0.5), nr = nr1),
                    matrix(rnorm(nr2*nc3, mean = 0.5, sd = 0.5), nr = nr2),
                    matrix(rnorm(nr3*nc3, mean = 1,   sd = 0.5), nr = nr3))
  )
  mat = mat[sample(nr, nr), sample(nc, nc)] # random shuffle rows and columns
  rownames(mat) = paste0("row", seq_len(nr))
  colnames(mat) = paste0("column", seq_len(nc))
  return(mat)
}

# colors
# Turquois
mat <- create_random_matrix(seed=123)
col_fun = colorRamp2(c(-1, 2), c("white", "#6BBFBD"))

pdf(file="~/src/TissueDeconv/figures/Tcell_turquoise.pdf", height=5, width=5)
Heatmap(mat, 
        col = col_fun, 
        rect_gp = gpar(col = "white", lwd = 2), 
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        width = unit(7, "cm"), 
        height = unit(7, "cm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE)
dev.off()

# Pink
mat <- create_random_matrix(seed=120)
col_fun = colorRamp2(c(-1.5, 1.5), c("white", "#B7477C"))

pdf(file="~/src/TissueDeconv/figures/Bcell_pink.pdf", height=5, width=5)
Heatmap(mat, 
        col = col_fun, 
        rect_gp = gpar(col = "white", lwd = 2), 
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        width = unit(7, "cm"), 
        height = unit(7, "cm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE)
dev.off()

# Bulk - grey
mat <- create_random_matrix(seed=30)
col_fun = colorRamp2(c(-1, 2), c("white", "black"))

pdf(file="~/src/TissueDeconv/figures/Bulk_grey.pdf", height=5, width=5)
Heatmap(mat, 
        col = col_fun, 
        rect_gp = gpar(col = "white", lwd = 2), 
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        width = unit(7, "cm"), 
        height = unit(7, "cm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE)
dev.off()

# Tumor dark violet
mat <- create_random_matrix(seed=11)
col_fun = colorRamp2(c(-1, 1.5), c("white", "#2D306E"))

pdf(file="~/src/TissueDeconv/figures/Tumor_violet.pdf", height=5, width=5)
Heatmap(mat, 
        col = col_fun, 
        rect_gp = gpar(col = "white", lwd = 2), 
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        width = unit(7, "cm"), 
        height = unit(7, "cm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE)
dev.off()

# NK cells green
mat <- create_random_matrix(seed=2)
col_fun = colorRamp2(c(-1.5, 1.5), c("white", "#13673E"))

pdf(file="~/src/TissueDeconv/figures/NK_green.pdf", height=5, width=5)
Heatmap(mat, 
        col = col_fun, 
        rect_gp = gpar(col = "white", lwd = 2), 
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        width = unit(7, "cm"), 
        height = unit(7, "cm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = FALSE)
dev.off()

# Confidence matrix
mat <- create_random_matrix(seed=2)
col_fun = colorRamp2(c(-0.5,2), c("#27AAE1", "#EF4136"))

pdf(file="~/src/TissueDeconv/figures/confidence_mat.pdf", height=5, width=5)
Heatmap(mat, 
        col = col_fun, 
        rect_gp = gpar(col = "white", lwd = 2), 
        show_column_dend = FALSE,
        show_row_dend = FALSE,
        width = unit(7, "cm"), 
        height = unit(7, "cm"),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_heatmap_legend = TRUE)
dev.off()


