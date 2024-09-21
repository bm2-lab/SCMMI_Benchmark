library(ggplot2)

source("plot_scmmib_table.r")
# setting row groups if available
p_tab_rank<-read.table("./test/scmmib_table_demo.txt",sep='\t',header=T)

mosaic_rna_atac_acc_row_info <- data.frame(id = p_tab_rank$method)
mosaic_rna_atac_acc_column_info <- data.frame(id = colnames(p_tab_rank),
                          group = c("Text","Text",
                                    rep("Embedding accuracy", 4),
                                    rep("Cell alignment accuracy",2),
                                    "Embedding accuracy",
                                    "Cell alignment accuracy"), 
                          geom = c(rep("text",2),
                                   rep("circle",6),
                                   rep("bar",2)),
                          width = c(1.5,8.5,rep(2,8)),
                          overlay = F)
mosaic_rna_atac_acc_palettes <- list(   "Text" = "black",
                    "Embedding accuracy" = "Blues",
                   "Cell alignment accuracy" = "Greens"
                   )
p<-plot_scmmib_summary(data = p_tab_rank,
                       row_info = mosaic_rna_atac_acc_row_info,
                       column_info =mosaic_rna_atac_acc_column_info,
                       palettes = mosaic_rna_atac_acc_palettes,
                       rank = TRUE, # T for ranking in rectangles
                       rect_normalize = T, # F if the size of rectangle is not normalized.
                       extend_figure = F # T if the size extend the margin of scmmib summary plots.
                       )

p
