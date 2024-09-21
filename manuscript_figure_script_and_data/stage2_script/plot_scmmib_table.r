library(dplyr)
library(scales)
library(ggimage)
library(cowplot)

add_column_if_missing <- function(df, ...) {
  column_values <- list(...)
  for (column_name in names(column_values)) {
    default_val <- rep(column_values[[column_name]], nrow(df))
    
    if (column_name %in% colnames(df)) {
      df[[column_name]] <- ifelse(is.na(df[[column_name]]), default_val, df[[column_name]])
    } else {
      df[[column_name]] <- default_val
    }
  }
  df
}

label_top_3 <- function(scores, ...) {
  ranks <- rank(scores, ...)
  ifelse(ranks <= 3, as.character(ranks), "")
}

plot_scmmib_summary<-function(data,row_info,column_info,palettes,rank=FALSE,rect_normalize=TRUE,extend_figure=FALSE){
  row_height <- 1.1
  row_space <- .1
  row_bigspace <- .5
  col_width <- 1.1
  col_space <- .2
  col_bigspace <- .5
  segment_data <- NULL
  
  
  if (!"group" %in% colnames(row_info) || all(is.na(row_info$group))) {
    row_info$group <- ""
    row_groups <- tibble(group = "")
    plot_row_annotation <- FALSE
  } else {
    plot_row_annotation <- TRUE
  }
  # 计算行坐标
  row_pos <-
    row_info %>%
    group_by(group) %>%
    dplyr::mutate(group_i = row_number()) %>%
    ungroup() %>%
    dplyr::mutate(
      row_i = row_number(),
      colour_background = group_i %% 2 == 1,
      do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
      ysep = ifelse(do_spacing, row_height + 2 * row_space, row_space),
      y = - (row_i * row_height + cumsum(ysep)),
      ymin = y - row_height / 2,
      ymax = y + row_height / 2
    )
  
  # 计算列坐标
  if (!"group" %in% colnames(column_info) || all(is.na(column_info$group))) {
    column_info$group <- ""
    plot_column_annotation <- FALSE
  } else {
    plot_column_annotation <- TRUE
  }
  
  column_info <-
    column_info %>%
    add_column_if_missing(width = col_width, overlay = FALSE)
  
  
  column_pos <-
    column_info %>%
    dplyr::mutate(
      do_spacing = c(FALSE, diff(as.integer(factor(group))) != 0),
      xsep = case_when(
        overlay ~ c(0, -head(width, -1)),
        do_spacing ~ col_bigspace,
        TRUE ~ col_space
      ),
      xwidth = case_when(
        overlay & width < 0 ~ width - xsep,
        overlay ~ -xsep,
        TRUE ~ width
      ),
      xmax = cumsum(xwidth + xsep),
      xmin = xmax - xwidth,
      x = xmin + xwidth / 2
    )
  
  
  
  ##########################
  #### CREATE GEOM DATA ####
  ##########################
  
  # gather circle data
  ind_circle <- which(column_info$geom == "circle")
  if(length(ind_circle) > 0){
    dat_mat <- as.matrix(data[, ind_circle],drop=F)
    col_palette <- data.frame(metric = colnames(dat_mat), 
                              group = column_info[match(colnames(dat_mat), column_info$id), "group"])
    
    col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])
    
    circle_data <- data.frame(label = unlist(lapply(colnames(dat_mat), 
                                                    function(x) rep(x, nrow(dat_mat)))), 
                              x0 = unlist(lapply(column_pos$x[ind_circle], 
                                                 function(x) rep(x, nrow(dat_mat)))), 
                              y0 = rep(row_pos$y, ncol(dat_mat)),
                              r = row_height/2*as.vector(sqrt(dat_mat))
    )
    for(l in unique(circle_data$label)){
      ind_l <- which(circle_data$label == l)
      circle_data[ind_l, "r"] <- rescale(circle_data[ind_l, "r"], to = c(0.05, 0.55), from = range(circle_data[ind_l, "r"], na.rm = T))
    }
    
    colors <- NULL
    
    
    for(i in 1:ncol(dat_mat)){
      palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
      colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
    }
    
    circle_data$colors <- colors
  }
  
  
  # gather bar data
  ind_bar <- which(column_info$geom == "bar")
  dat_mat <- as.matrix(data[, ind_bar],drop=F)
  
  if (rect_normalize){
    dat_mat <- data.frame(dat_mat) %>%
      dplyr::mutate(across(everything(), ~./max(.,na.rm = T)))
    dat_mat <- as.matrix(dat_mat,drop=F)
  }
  # if (length(ind_bar)==1){
  colnames(dat_mat) = column_info[ind_bar,1] # avoid renames space in column names
  # }
  col_palette <- data.frame(metric = colnames(dat_mat), 
                            group = column_info[match(colnames(dat_mat), column_info$id), "group"])
  
  col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])
  
  
  rect_data <- data.frame(label = unlist(lapply(colnames(dat_mat), 
                                                function(x) rep(x, nrow(dat_mat)))),
                          method = rep(row_info$id, ncol(dat_mat)),
                          value = as.vector(dat_mat),
                          xmin = unlist(lapply(column_pos[ind_bar, "xmin"], 
                                               function(x) rep(x, nrow(dat_mat)))),
                          xmax = unlist(lapply(column_pos[ind_bar, "xmax"], 
                                               function(x) rep(x, nrow(dat_mat)))),
                          ymin = rep(row_pos$ymin, ncol(dat_mat)),
                          ymax = rep(row_pos$ymax, ncol(dat_mat)),
                          xwidth = unlist(lapply(column_pos[ind_bar, "xwidth"], 
                                                 function(x) rep(x, nrow(dat_mat))))
  )
  rect_data <- rect_data %>%
    add_column_if_missing(hjust = 0) %>%
    dplyr::mutate(
      xmin = xmin + (1 - value) * xwidth * hjust,
      xmax = xmax - (1 - value) * xwidth * (1 - hjust)
    )
  
  colors <- NULL
  for(i in 1:ncol(dat_mat)){
    palette <- colorRampPalette(rev(brewer.pal(9, col_palette$name_palette[[i]])))(nrow(data)-sum(is.na(dat_mat[,i])))
    colors <- c(colors, palette[rank(dat_mat[,i], ties.method = "average", na.last = "keep")])
  }
  
  rect_data$colors <- colors
  
  
  
  # gather text data
  ind_text <- which(column_info$geom == "text")
  dat_mat <- as.matrix(data[, ind_text],drop=F)
  
  
  text_data <- data.frame(label_value = as.vector(dat_mat), 
                          group = rep(colnames(dat_mat), each = nrow(dat_mat)),
                          xmin = unlist(lapply(column_pos[ind_text, "xmin"], 
                                               function(x) rep(x, nrow(dat_mat)))),
                          xmax = unlist(lapply(column_pos[ind_text, "xmax"], 
                                               function(x) rep(x, nrow(dat_mat)))),
                          ymin = rep(row_pos$ymin, ncol(dat_mat)),
                          ymax = rep(row_pos$ymax, ncol(dat_mat)),
                          size = 4, fontface = "plain", stringsAsFactors = F)
  col_palette <- data.frame(metric = colnames(dat_mat), 
                            group = column_info[match(colnames(dat_mat), column_info$id), "group"])
  
  col_palette$name_palette <- lapply(col_palette$group, function(x) palettes[[as.character(x)]])
  
  colors <- NULL
  for(i in 1:ncol(dat_mat)){
    colors <- c(colors, rep(col_palette$name_palette[[i]],nrow(data)-sum(is.na(dat_mat[,i]))))
  }
  text_data$colors<-colors
  # text_data$colors <- "black"
  
  text_data[text_data$label_value == "python", "colors"] <- "darkgreen"
  text_data[text_data$label_value == "R", "colors"] <- "darkblue"
  
  # replace yes to √ 
  text_data$label_value <- plyr::mapvalues(text_data$label_value, from = c("yes", "no"), 
                                           to = c("\u221A", ""))
  
  text_data[text_data$label_value == "\u221A" , "size"] <- 4
  text_data[text_data$label_value == "\u221A", "fontface"] <- "bold"
  
  
  # ADD top3 ranking for each bar column
  if (rank){
    cols_bar <- unique(rect_data$label)
    cols_bar <- as.character(cols_bar[!is.na(cols_bar)])
    for(c in cols_bar){
      rect_tmp <- rect_data[rect_data$label == c,]
      rect_tmp <- add_column(rect_tmp, "label_value" = as.character(rank(-rect_tmp$value, ties.method = "min")))
      rect_tmp <- rect_tmp[rect_tmp$label_value %in% c("1", "2", "3"), c("label_value", "xmin", "xmax", "ymin", "ymax")]
      rect_tmp <- add_column(rect_tmp, "size" = 3.5, .after = "ymax")
      rect_tmp <- add_column(rect_tmp, "colors" = "black", .after = "size")
      rect_tmp <- add_column(rect_tmp, "fontface" = "plain", .after = "colors")
      rect_tmp <- add_column(rect_tmp, "group" = "top3", .after = "fontface")
      text_data <- bind_rows(text_data, rect_tmp) 
    }
  }
  
  # ADD COLUMN NAMES, 不加的单独列出来。
  df <- column_pos %>% filter(id != "method") # %>% filter(id != "Rank")
  
  if (nrow(df) > 0) {
    segment_data <- segment_data %>% bind_rows(
      df %>% transmute(x = x, xend = x, y = -.3, yend = -.1, size = .5)
    )
    text_data <-
      bind_rows(
        text_data,
        df %>% transmute(
          xmin = x, xmax = x, ymin = 0, ymax = -0.5,
          angle = 30, vjust = 0, hjust = 0,
          label_value = id, 
          size = 3
        )
      )
  }
  
  
  # GENERATE ROW ANNOTATION，应该没有
  if (plot_row_annotation) {
    row_annotation <-
      row_pos %>% 
      select(group, ymin, ymax) %>%
      group_by(group) %>%
      summarise(
        ymin = min(ymin),
        ymax = max(ymax),
        y = (ymin + ymax) / 2
      ) %>%
      ungroup() %>%
      dplyr::mutate(xmin = -.5, xmax = 5) %>%
      filter(!is.na(group), group != "")
    
    text_data <- text_data %>% bind_rows(
      row_annotation %>%
        transmute(xmin, xmax, ymin = ymax + row_space, label_value = group %>% gsub("\n", " ", .), 
                  hjust = 0, vjust = .5, fontface = "bold", size = 4) %>%
        dplyr::mutate(ymax = ymin + row_height)
    )
  }
  # segment_data
  suppressWarnings({
    minimum_x <- min(column_pos$xmin, segment_data$x, segment_data$xend, 
                     text_data$xmin, na.rm = TRUE)
    maximum_x <- max(column_pos$xmax, segment_data$x, segment_data$xend, 
                     text_data$xmax, na.rm = TRUE)
    minimum_y <- min(row_pos$ymin, segment_data$y, segment_data$yend,  
                     text_data$ymin, na.rm = TRUE)
    maximum_y <- max(row_pos$ymax, segment_data$y, segment_data$yend, 
                     text_data$ymax, na.rm = TRUE)
  })
  
  ####################################
  ###   CREATE HARDCODED LEGENDS   ### 
  ####################################
  atac=FALSE
  
  x_min_output <- minimum_x+0.5
  x_min_scaling <- minimum_x + 5.5
  x_min_ranking <- ifelse(atac, minimum_x + 5.5, minimum_x + 10.5)
  x_min_score <-  ifelse(atac, minimum_x + 11, minimum_x + 17)
  
  leg_max_y <- minimum_y - .5
  
  # Create legend for Output(不需要了)
  leg_min_x <- x_min_output
  
  # CREATE LEGEND for ranking colors
  leg_min_x <- x_min_ranking
  rank_groups <- as.character(column_info[column_info$geom == "bar", "group"])
  
  j = 0
  tmp2 = palettes
  item_names = names(tmp2)
  rank_minimum_x = list()
  for (i in seq_along(tmp2)){
    if (i >1){
      rank_minimum_x[[item_names[i]]]<-leg_min_x+j
      j = j+1
    }
  }
  
  leg_max_x <- leg_min_x+j
  
  
  rank_title_data <- data.frame(xmin = leg_min_x, 
                                xmax = leg_min_x+ 2, 
                                ymin = leg_max_y - 1, 
                                ymax = leg_max_y, 
                                label_value = "Ranking", 
                                hjust = 0, vjust = 0, 
                                fontface = "bold")
  
  for(rg in rank_groups){
    rank_palette <- colorRampPalette(rev(brewer.pal(9, palettes[[rg]])))(5)
    
    
    rank_data <- tibble(xmin = rank_minimum_x[[rg]],
                        xmax = rank_minimum_x[[rg]] + .8,
                        ymin = seq(leg_max_y-4, leg_max_y - 2, by = .5),
                        ymax = seq(leg_max_y-3.5, leg_max_y -1.5, by = .5),
                        border = TRUE,
                        colors = rank_palette
    )
    rect_data <- bind_rows(rect_data, rank_data)
    
  }
  
  # add text next to the arrow
  arrow_data <- data.frame(x = leg_max_x + 1.5, 
                           xend = leg_max_x +1.5, 
                           y = leg_max_y-4, 
                           yend = leg_max_y -1.5)
  
  
  arrow_text <- data.frame(xmin = leg_max_x +2, 
                           xmax = leg_max_x +2.5, 
                           ymin = c(leg_max_y-2, leg_max_y-4), 
                           ymax = c(leg_max_y-1.5, leg_max_y-3.5 ), 
                           label_value = c("1", as.character(nrow(data))), 
                           hjust = 0, vjust = 0, size = 2.5)
  
  
  text_data <- bind_rows(text_data, rank_title_data, arrow_text)
  
  
  # circle legend
  if(length(ind_circle) > 0){
    cir_minimum_x <- x_min_score
    
    cir_legend_size <- 1
    cir_legend_space <- .1
    
    cir_legend_dat <-
      data.frame(
        value = seq(0, 1, by = .2),
        r = row_height/2*seq(0, 1, by = .2)
      )
    cir_legend_dat$r <- rescale(cir_legend_dat$r, to = c(0.05, 0.55), from = range(cir_legend_dat$r, na.rm = T))
    
    x0 <- vector("integer", nrow(cir_legend_dat))
    for(i in 1:length(x0)){
      if(i == 1){
        x0[i] <- cir_minimum_x + cir_legend_space + cir_legend_dat$r[i]
      }
      else {
        x0[i] <- x0[i-1] + cir_legend_dat$r[i-1] + cir_legend_space + cir_legend_dat$r[i]
      }
    }
    
    cir_legend_dat$x0 <- x0
    cir_legend_min_y <- leg_max_y-4
    cir_legend_dat$y0 <- cir_legend_min_y + 1 + cir_legend_dat$r
    
    cir_legend_dat$colors <- NULL
    cir_maximum_x <- max(cir_legend_dat$x0)
    
    cir_title_data <- data.frame(xmin = cir_minimum_x, 
                                 xmax = cir_maximum_x, 
                                 ymin = leg_max_y -1, 
                                 ymax = leg_max_y,
                                 label_value = "Score", 
                                 hjust = 0, vjust = 0, fontface = "bold")
    
    cir_value_data <- data.frame(xmin = cir_legend_dat$x0 - cir_legend_dat$r,
                                 xmax = cir_legend_dat$x0 + cir_legend_dat$r,
                                 ymin = cir_legend_min_y,
                                 ymax = cir_legend_min_y +3,
                                 hjust = .5, vjust = 0, size = 2.5,
                                 label_value = ifelse(cir_legend_dat$value %in% c(0, 1), 
                                                      paste0(cir_legend_dat$value*100, "%"), ""))
    
    circle_data <- bind_rows(circle_data, cir_legend_dat)
    text_data <- bind_rows(text_data, cir_title_data, cir_value_data)
  }
  
  minimum_y <- min(minimum_y, min(text_data$ymin, na.rm = TRUE))
  
  
  ########################################################
  ######  plot according to row_pos,col_pos,rect_data
  ########################################################
  g <-
    ggplot() +
    coord_equal(expand = FALSE) +
    scale_alpha_identity() +
    scale_colour_identity() +
    scale_fill_identity() +
    scale_size_identity() +
    scale_linetype_identity() +
    cowplot::theme_nothing() +
    theme(text=element_text(family = "serif"))
  
  # PLOT ROW BACKGROUNDS
  df <- row_pos %>% filter(colour_background)
  if (nrow(df) > 0) {
    g <- g + geom_rect(aes(xmin = min(column_pos$xmin)-.25, xmax = max(column_pos$xmax)+.25, ymin = ymin - (row_space / 2), ymax = ymax + (row_space / 2)), df, fill = "#DDDDDD")
  } 
  
  # PLOT CIRCLES
  if (length(ind_circle) > 0) {
    g <- g + ggforce::geom_circle(aes(x0 = x0, y0 = y0, fill= colors, r = r), circle_data, size=.25)
  }
  
  
  # PLOT RECTANGLES
  if (nrow(rect_data) > 0) {
    # add defaults for optional values
    rect_data <- rect_data %>%
      add_column_if_missing(alpha = 1, border = TRUE, border_colour = "black") %>%
      dplyr::mutate(border_colour = ifelse(border, border_colour, NA))
    
    g <- g + geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = colors, colour = border_colour, alpha = alpha), rect_data, linewidth = .25)
  }
  # add defaults for optional values
  text_data <- text_data %>%
    add_column_if_missing(
      hjust = .5,
      vjust = .5,
      size = 3,
      fontface = "plain",
      colors = "black",
      lineheight = 1,
      angle = 0
    ) %>%
    dplyr::mutate(
      angle2 = angle / 360 * 2 * pi,
      cosa = cos(angle2) %>% round(2),
      sina = sin(angle2) %>% round(2),
      alphax = ifelse(cosa < 0, 1 - hjust, hjust) * abs(cosa) + ifelse(sina > 0, 1 - vjust, vjust) * abs(sina),
      alphay = ifelse(sina < 0, 1 - hjust, hjust) * abs(sina) + ifelse(cosa < 0, 1 - vjust, vjust) * abs(cosa),
      x = (1 - alphax) * xmin + alphax * xmax,
      y = (1 - alphay) * ymin + alphay * ymax
    ) %>%
    filter(label_value != "")
  
  # text_data[rownames(subset(text_data,text_data$group=="method")),"hjust"]=0.3 # too long algorithm names
  # Set fontface for legend bold
  # text_data[text_data$label_value == "Ranking", "fontface"] <- "bold"
  text_data[text_data$label_value == "Rank", "fontface"] <- "bold" 
  # text_data[text_data$label_value == "method", "fontface"] <- "bold" # 
  # Set fontface for ranking numbers and software as bold，
  # text_data[1:nrow(data), "fontface"] <- "bold"
  
  # subset text_data to left-aligned rows, 把算法名字和排名分出来，单独画
  text_data_left <- text_data[which(text_data$group == "method" | text_data$group == "top3"), ]
  text_data <- text_data[-which(text_data$group == "method" | text_data$group == "top3"), ]
  
  g <- g + geom_text(aes(x = x, y = y, label = label_value, colour = colors, hjust = hjust, vjust = vjust, size = size, fontface = fontface, angle = angle), data = text_data)
  
  text_data_left[text_data_left$group == "method", "x"] <- text_data_left[text_data_left$group == "method", "x"] - 3
  
  text_data_left[text_data_left$group == "top3", "x"] <- text_data_left[text_data_left$group == "top3", "xmin"] + .3
  # text_data_left[text_data_left$group == "method", "x"] <- text_data_left[text_data_left$group == "method", "x"] + .5
  
  g <- g + geom_text(aes(x = x, y = y, label = label_value, colour = colors, hjust = 'left', vjust = vjust, size = size, fontface = fontface, angle = angle), data = text_data_left)
  
  # PLOT SEGMENTS
  if (nrow(segment_data) > 0) {
    # add defaults for optional values
    segment_data <- segment_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")
    
    g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), segment_data)
  }
  
  # PLOT ARROW RANKING
  if (nrow(arrow_data) > 0) {
    # add defaults for optional values
    arrow_data <- arrow_data %>% add_column_if_missing(size = .5, colour = "black", linetype = "solid")
    
    g <- g + geom_segment(aes(x = x, xend = xend, y = y, yend = yend, size = size, colour = colour, linetype = linetype), arrow_data, arrow = arrow(length = unit(0.1, "cm")), lineend = "round", linejoin = "bevel")
  }
  
  # ADD SIZE
  # reserve a bit more room for text that wants to go outside the frame
  minimum_x <- minimum_x - 3
  maximum_x <- maximum_x + 6
  
  minimum_y <- minimum_y - 3
  maximum_y <- maximum_y + 5
  if (extend_figure){
    maximum_x <- maximum_x + 5
    maximum_y <- maximum_y + 10
  }
  g$width <- maximum_x - minimum_x
  g$height <- maximum_y - minimum_y
  
  g <- g + expand_limits(x = c(minimum_x, maximum_x), y = c(minimum_y, maximum_y))
  return(g)
}