##### Plotting Function for Epistasis output ####
#                                              ##
#      Author: Matthias Arnold, PhD            ##
#                                              ##
#               Nov 29, 2024                   ##
#                                              ##
#################################################

require(pheatmap)
require(ggpubr)
require(ggplotify)

starManhattan <-
  function(data,
           x = "xcord",
           trait,
           color = "ac",
           traitmax,
           xticks = "xticks") {
    x <- ggscatter(
      data,
      x = x,
      y = trait,
      color = color,
      palette = rep(RColorBrewer::brewer.pal(9, "Blues")[c(7, 9)], 7),
      parse = F,
      label = "label",
      label.select = list(criteria = eval(
        paste0("`", trait, "` > -log10(0.05/11) & `", traitmax, "`")
      )),
      repel = T,
      label.rectangle = T,
      star.plot = T,
      star.plot.lty = 3,
      lineheight = 0.2
    ) +
      scale_x_continuous(breaks = unique(data[, xticks]),
                         labels = unique(as.character(data[, color]))) +
      xlab("Genetic model") + ylab(expression(paste('-', log[10](P)))) + ylim(c(0, 5.9)) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(
          angle = 45,
          vjust = 1,
          hjust = 1
        ),
        plot.margin = unit(c(0.02, 0, 0, .01), units = "npc")
      )
    x$layers[[3]]$aes_params$fill <- rgb(1, 1, 1, .75)
    x$layers[[3]]$aes_params$alpha <- NULL
    x$layers[[3]]$aes_params$size <- 3
    x$layers[[3]]$geom_params$min.segment.length <- 0
    x$layers[[3]]$geom_params$force <- 6
    x$layers[[3]]$geom_params$force_pull <- 0
    x
  }

starManhattanGeno <-
  function(data,
           x = "xcord",
           trait,
           color = "ac",
           traitmax,
           xticks = "xticks") {
    x <- ggscatter(
      data,
      x = x,
      y = trait,
      color = color,
      palette = rep(RColorBrewer::brewer.pal(9, "Blues")[c(7, 9)], 7),
      parse = F,
      label = "label",
      label.select = list(criteria = eval(
        paste0("`", trait, "` > -log10(0.01) & `", traitmax, "`")
      )),
      repel = T,
      label.rectangle = T,
      star.plot = T,
      star.plot.lty = 3,
      lineheight = 0.2
    ) +
      scale_x_continuous(breaks = unique(data[, xticks]),
                         labels = unique(as.character(data[, color]))) +
      xlab("Genetic model") + ylab(expression(paste('-', log[10](P)))) + ylim(c(0, 3.7)) +
      theme(legend.position = "none",
            axis.text.x = element_text(
              angle = 45,
              vjust = 1,
              hjust = 1
            ))
    x$layers[[3]]$aes_params$fill <- rgb(1, 1, 1, .75)
    x$layers[[3]]$aes_params$alpha <- NULL
    x$layers[[3]]$aes_params$size <- 3
    x$layers[[3]]$geom_params$min.segment.length <- 0
    x$layers[[3]]$geom_params$force <- 6
    x$layers[[3]]$geom_params$force_pull <- 0
    x
  }

epiHeatmap <- function(plot.data) {
  lbls <-
    c(0, seq(floor(min(plot.data[plot.data > 1.5])), round(max(plot.data), digits = 1), by = .5))
  plot.data[plot.data > 0] <- plot.data[plot.data > 0] - 1.5
  brks <- c(seq(0, round(max(plot.data), digits = 1), by = .5))
  brks[-c(1, 2, 4)] <- brks[-c(1, 2, 4)] - 0.1
  brks[2] <- brks[2] - 0.05
  
  x <-
    pheatmap(
      t(plot.data[, colSums(plot.data) > 0]),
      cluster_cols = F,
      gaps_col = c(5, 10),
      color = c(rep("white", 30), colorRampPalette(
        c("white", RColorBrewer::brewer.pal(9, "Blues")[3:9])
      )(150)),
      #legend_breaks = c(0, seq(floor(min(plot.data[plot.data>0])), round(max(plot.data), digits = 1), by = .5)),
      legend_breaks = brks,
      legend_labels = lbls,
      angle_col = 45,
      silent = T,
      treeheight_row = 20,
      fontsize_row = 12,
      fontsize_col = 12,
      labels_col = c(
        expression(paste("CSF A", beta[1 - 42])),
        "CSF p-tau",
        "FDG-PET",
        "ADAS-Cog. 13",
        "Diagnosis",
        "Global cognition",
        "Amyloid load",
        "PHF tangle load",
        "Global pathology",
        "Diagnosis",
        "Diagnosis [staged]",
        "Diagnosis [binary]"
      )
    )
  
  xcg <- x$gtable$grobs[[5]]$children[[1]]$x
  ycg <- x$gtable$grobs[[5]]$children[[1]]$y
  wcg <- x$gtable$grobs[[5]]$children[[1]]$width
  hcg <- x$gtable$grobs[[5]]$children[[1]]$height
  
  x$gtable$layout[5, c(1:4)] <- c(5, 4, 5, 4)
  x$gtable$grobs[[5]]$children[[1]]$x <- ycg
  x$gtable$grobs[[5]]$children[[1]]$y <-
    xcg + unit(0.3, units = "npc")
  x$gtable$grobs[[5]]$children[[1]]$width <- hcg
  x$gtable$grobs[[5]]$children[[1]]$height <- wcg * 2
  
  xcg <- x$gtable$grobs[[5]]$children[[2]]$x
  ycg <- x$gtable$grobs[[5]]$children[[2]]$y
  wcg <- x$gtable$grobs[[5]]$children[[2]]$width
  hcg <- x$gtable$grobs[[5]]$children[[2]]$height
  
  x$gtable$grobs[[5]]$children[[2]]$x <- ycg
  x$gtable$grobs[[5]]$children[[2]]$y <-
    xcg + unit(0.05, units = "npc")
  x$gtable$grobs[[5]]$children[[2]]$width <- hcg
  x$gtable$grobs[[5]]$children[[2]]$height <- wcg
  x$gtable$grobs[[5]]$children[[2]]$gp$fontsize <- 12
  
  x <- x %>% as.ggplot +
    ylab("Genetic model") +
    annotate(
      geom = "text",
      x = .77,
      y = 0.15,
      vjust = 0,
      label = "`-` * log[10](P)",
      parse = T,
      size = 4.33333
    ) +
    annotate(
      geom = "text",
      x = .18,
      y = 1,
      vjust = 0,
      label = "ADNI",
      parse = F,
      size = 4.33333,
      fontface = "bold"
    ) +
    annotate(
      geom = "text",
      x = .415,
      y = 1,
      vjust = 0,
      label = "ROS/MAP",
      parse = F,
      size = 4.33333,
      fontface = "bold"
    ) +
    annotate(
      geom = "text",
      x = .575,
      y = 1,
      vjust = 0,
      label = "Mayo",
      parse = F,
      size = 4.33333,
      fontface = "bold"
    ) +
    labs(tag = "a") + coord_cartesian(clip = "off") +
    theme(
      axis.title.y = element_text(
        face = "plain",
        angle = 270,
        vjust = 135,
        hjust = 0.32,
        size = 12
      ),
      plot.tag = element_text(size = rel(1.235), face = "bold"),
      plot.margin = unit(c(.018, 0, 0, .016), units = "npc")
    )
  x
}