custom_dotplot <- function(object, x = "Cluster", colorBy = "p.adjust",
                                         showCategory = 5, by = "geneRatio", size = "geneRatio",
                                         split = NULL, includeAll = TRUE,
                                         font.size = 12, title = "", label_format = 30,
                                         group = FALSE, shape = FALSE, facet = NULL, strip_width = 15) {
  color <- NULL
  if (is.null(size)) size <- by ## by may be deprecated in future release
  
  if (!is.null(facet) && facet == "intersect") {
    object <- append_intersect(object)
  }
  
  df <- fortify(object, showCategory = showCategory, by = size,
                includeAll = includeAll, split = split)
  
  # Ensure the 'Description' factor levels match the order in 'showCategory'
  if (is.character(showCategory) || is.list(showCategory)) {
    df$Description <- factor(df$Description, levels = showCategory)
  } else {
    # Sort based on the original logic if 'showCategory' is a number
    idx <- order(df[[by]], decreasing = TRUE)
    df$Description <- factor(df$Description, levels = rev(unique(df$Description[idx])))
  }
  
  label_func <- default_labeller(label_format)
  if (is.function(label_format)) {
    label_func <- label_format
  }
  
  if (size %in% c("rowPercentage", "count", "geneRatio")) {
    by2 <- switch(size, rowPercentage = "Percentage",
                  count = "Count",
                  geneRatio = "GeneRatio")
  } else {
    by2 <- size
  }
  
  p <- ggplot(df, aes_string(x = x, y = "Description", size = by2)) +
    scale_y_discrete(labels = label_func)
  
  if (group) {
    p <- p + geom_line(aes_string(color = "Cluster", group = "Cluster"), size = .3) +
      ggnewscale::new_scale_colour()
  }
  
  if (shape) {
    check_installed('ggstar', 'for `dotplot()` with `shape = TRUE`.')
    ggstar <- "ggstar"
    require(ggstar, character.only = TRUE)
    p <- p + ggstar::geom_star(aes_string(starshape = "Cluster", fill = colorBy)) +
      set_enrichplot_color(type = "fill")
  } else {
    p <- p + geom_point(aes_string(fill = colorBy)) +
      aes(shape = I(enrichplot_point_shape))
  }
  
  p <- p + set_enrichplot_color(type = "fill") +
    ylab(NULL) + ggtitle(title) + DOSE::theme_dose(font.size) +
    scale_size_continuous(range = c(3, 8))
  
  if (!is.null(facet)) {
    p <- p + facet_grid(.data[[facet]] ~ .,
                        scales = "free", space = 'free',
                        switch = 'y',
                        labeller = ggplot2::label_wrap_gen(strip_width)) +
      theme(strip.text = element_text(size = 14))
  }
  
  class(p) <- c("enrichplotDot", class(p))
  
  return(p)
}
default_labeller <- function(n) {
  fun <- function(str){
    str <- gsub("_", " ", str)
    yulab.utils::str_wrap(str, n)
  }
  
  structure(fun, class = "labeller")
}
enrichplot_point_shape <- ggfun:::enrichplot_point_shape

set_enrichplot_color <- function(colors = get_enrichplot_color(2), 
                                 type = "color", name = NULL, .fun = NULL, reverse=TRUE, ...) {
  
  type <- match.arg(type, c("color", "colour", "fill"))
  if (!reverse) colors = rev(colors)
  n <- length(colors)
  if (n < 2) {
    stop("'colors' should be of length >= 2")
  } else if (n == 2) {
    params <- list(low = colors[1], high = colors[2])
    fn_suffix <- "continuous"
  } else if (n == 3) {
    params <- list(low = colors[1], mid = colors[2], high = colors[3])
    fn_suffix <- "gradient2"
  } else {
    params <- list(colors = colors) 
    fn_suffix <- "gradientn"   
  }
  
  if (!is.null(.fun)) {
    if (n == 3) { 
      # should determine parameter for user selected functions: 'gradient2' or 'gradientn'
      fn_type <- which_scale_fun(.fun)
      if (fn_type == "gradientn") {
        params <- list(colors = colors) 
      } else {
        params <- list(low = colors[1], mid = colors[2], high = colors[3])
      }
    }
  } else {
    fn <- sprintf("scale_%s_%s", type, fn_suffix)
    .fun <- getFromNamespace(fn, "ggplot2")
  }
  
  params$guide <- guide_colorbar(reverse=reverse, order=1)
  params$name <- name # no legend name setting by default as 'name = NULL'
  
  params <- modifyList(params, list(...))
  
  do.call(.fun, params)
}

get_enrichplot_color <- function(n = 2) {
  colors <- getOption("enrichplot.colours")
  if (!is.null(colors)) return(colors)
  
  if (n != 2 && n != 3) stop("'n' should be 2 or 3")
  
  colors = c("#e06663", "#327eba")
  if (n == 2) return(colors)
  
  if (n == 3) return(c(colors[1], "white", colors[2]))
}
