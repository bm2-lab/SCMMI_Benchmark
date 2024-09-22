
## User manual of `plot_scmmib_summary.r`

```R
plot_scmmib_summary(data,row_info,column_info,
palettes,rank=FALSE,rect_normalize=TRUE,
extend_figure=FALSE)
```
- data: data.frame. the method column name must be "method" <br> 
- row_info: data.frame of row information. same `id` col with data  `method` col.
- column_info: data.frame of column information. same colnames with data.
- rank: whether label top 3 method in the rectangle.
- rect_normalize: whetherthe length of rect is normalized.
- extend_figure: whether add extra white space for too wide or too long figure.


```R
# Example

cite_usb_row_info <- data.frame(id = cite_rank2$method)
# if the group is set, user should also set same length of geom, and width of the geom.
cite_usb_column_info <- data.frame(id = colnames(cite_rank2),
                          group = c("Text","Text", "Text", "Text", "Text", 
                                    rep("Code&Paper", 5),
                                    rep("Scalability", 7), 
                                    "Overall score"), 
                          geom = c(rep("text",17),
                                   rep("bar",1)),
                          width = c(1.5,9.5,3.5,3, rep(2,14)),
                          overlay = F)
# For all group in column and row, the color should be listed in palettes list() object
cite_usb_palettes <- list(   "Text" = "black",
                    "Code&Paper" = basecol(3)[1],
                   "Scalability" = basecol(3)[2],
                   "Overall score" = "Oranges"
                   )
p<-plot_scmmib_summary(data = cite_rank2,
                       row_info = cite_usb_row_info,
                       column_info =cite_usb_column_info,
                       palettes = cite_usb_palettes,
                       rank = FALSE
                       )
```