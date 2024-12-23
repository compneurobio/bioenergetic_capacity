##### Plotting Functions for SGI output ####
#                                         ##
#      Author: Matthias Arnold, PhD       ##
#                                         ##
#               Nov 28, 2024              ##
#                                         ##
############################################

require(ggpubr)
require(ggstatsplot)
require(dplyr)

## Wrapper function for ggboxplot with jitter for SGI output
SgiBoxScatter <-
  function(data,
           xvar,
           yvar,
           xlbl = "Cluster",
           ylbl = NULL,
           tag = " ") {
    if (is.null(ylbl)) {
      ylbl <- yvar
    }
    p <- ggboxplot(
      data,
      x = xvar,
      y = yvar,
      add = "jitter",
      notch = T,
      fill = xvar,
      add.params = list(
        alpha = .5,
        shape = 21,
        color = "black",
        fill = xvar
      )
    ) +
      stat_compare_means(
        label = "p.format",
        comparisons = list((
          data[, xvar] %>% na.omit %>% unique %>% as.character
        )),
        method = "t.test",
        method.args = list(var.equal = T),
        label.x = 1.35,
        size = rel(4.09)
      ) + xlab(xlbl) +
      ylab(ylbl) +
      labs(tag = tag) +
      coord_cartesian(clip = "off") +
      theme(
        legend.position = "none",
        plot.tag = element_text(face = "bold"),
        plot.margin = unit(c(.05, .03, .03, .03), units = "npc"),
        plot.tag.position = c(.01, 1)
      )
    
    if (any(is.na(data[, xvar]))) {
      p <-
        p + scale_x_discrete(labels = c(unique(data[, xvar]) %>% 
                                          sort(na.last = NA), "other"))
    }
    return(p)
  }

## Wrapper function for ggbarstats for SGI output (phenotype associations)
SgiStackedBar <-
  function(data,
           xvar,
           yvar,
           xlbl = "Cluster",
           ylbl = NULL,
           level,
           lbls = NULL,
           as) {
    if (is.null(ylbl)) {
      ylbl = yvar
    }
    if (is.null(lbls)) {
      ifelse(is.factor(data[, xvar]), lbls <-
               rev(levels(data[, xvar])), lbls <-
               rev(levels(factor(data[, xvar]))))
    }
    x <- dplyr::enquo(xvar)
    y <- dplyr::enquo(yvar)
    p <- data %>% ggbarstats(
      data = .,
      x = !!x,
      y = !!y,
      sampling.plan = "hypergeom",
      xlab = xlbl,
      legend.title = xvar,
      palette = "Set2",
      results.subtitle = F,
      subtitle = format(as$results[[xvar]][level, 4], digits = 2),
      proportion.test = F,
      sample.size.label = "",
      label.args = list(
        alpha = 1,
        fill = "white",
        size = rel(4.09)
      ),
      messages = FALSE
    ) + theme(
      axis.text = element_text(colour = "black", size = rel(1.19)),
      axis.title.x = element_text(face = "plain", size = rel(1.19)),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      plot.margin = unit(c(0, .03, .03, .03), units = "npc"),
      plot.subtitle = element_text(
        size = rel(1.19),
        hjust = .5,
        vjust = -4.5
      ),
      legend.box.margin = margin(l = -.1, unit = "npc")
    ) +
      scale_fill_brewer(name = xvar,
                        palette = "Set2",
                        labels = lbls) +
      geom_segment(aes(
        x = 1,
        y = 1.04,
        xend = 2,
        yend = 1.04
      ), size = rel(.1)) +
      geom_segment(aes(
        x = 1,
        y = 1.04,
        xend = 1,
        yend = 1.02
      ), size = rel(.1)) +
      geom_segment(aes(
        x = 2,
        y = 1.04,
        xend = 2,
        yend = 1.02
      ), size = rel(.1))
    p$layers[[3]] <- NULL
    return(p)
  }