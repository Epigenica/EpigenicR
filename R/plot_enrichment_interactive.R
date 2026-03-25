#' Interactive enrichment profile plot
#'
#' Generates an interactive multi-sample enrichment profile plot from binned
#' signal values. Lines are colored by sample class (control vs non-control),
#' and hovered samples are highlighted while non-hovered samples are dimmed.
#'
#' @param df A data frame containing at least the columns
#'   \code{index}, \code{mean}, and \code{sample}. Required unless
#'   \code{epk}, \code{marker}, and \code{loci} are all provided.
#' @param epk Optional EPK object. When supplied together with \code{marker}
#'   and \code{loci}, the profile data frame is extracted from
#'   \code{epk$enrichment_results$enrichment_profile[[loci]][[marker]]}.
#' @param marker Character; marker name used to look up the profile inside
#'   \code{epk$enrichment_results$enrichment_profile[[loci]]}.
#' @param loci Character; loci key (e.g. \code{"protein_coding"},
#'   \code{"CpG_islands"}) used to look up the profile list inside
#'   \code{epk$enrichment_results$enrichment_profile}.
#' @param bin_size Numeric; bin size in base pairs. Default: \code{100}.
#' @param window_bp Numeric; half window size in base pairs. Default:
#'   \code{2500} (i.e. +/- 2.5 kb).
#' @param highlight_opacity Numeric in [0,1]; opacity applied to non-hovered
#'   traces. Default: \code{0.5}.
#' @param mid_coord Character; reference point used for center label and x-axis
#'   title. One of \code{"center"} or \code{"start"}. Default:
#'   \code{"center"}.
#' @param plot_title Character; title for the plot. Default:
#'   \code{"Enrichment Profile"}.
#' @param control_samples Optional character vector of sample names to treat as
#'   controls. If provided, takes precedence over \code{control_pattern}.
#' @param control_pattern Optional regex pattern used to classify controls when
#'   \code{control_samples} is \code{NULL}.
#' @param control_color Character color for control traces. Default:
#'   \code{"rgba(0,0,0,0.90)"}.
#' @param other_color Character color for non-control traces. Default:
#'   \code{"rgba(120,180,255,0.55)"}.
#' @param hover_color Character color used for highlighted trace on hover.
#'   Default: \code{"rgba(0,120,255,1)"}.
#' @param group_by Optional character scalar naming a column used for color
#'   grouping (e.g. \code{"condition"}, \code{"batch"}). If provided,
#'   samples are colored by this grouping instead of control vs non-control.
#' @param metadata Optional data frame containing sample annotations to join
#'   into \code{df} before plotting (for example, condition assignments).
#' @param metadata_sample_col Character scalar naming the sample-id column in
#'   \code{metadata}. Default: \code{"sample"}.
#' @param group_palette Optional named character vector mapping group names to
#'   colors. When named, names should match the observed \code{group_by}
#'   levels. Missing levels are auto-filled (with warning); extra names are
#'   ignored (with warning).
#' @param control_overrides_group Logical; when \code{TRUE} and controls are
#'   defined, control samples use \code{control_color} even when
#'   \code{group_by} is provided. Default: \code{FALSE}.
#'
#' @return A \code{plotly} htmlwidget.
#'
#' @details
#' The x-axis is computed as genomic distance in bp from the mid-window and
#' labeled as \code{-window}, \code{mid_coord}, and \code{+window}. A dashed
#' vertical line marks \code{0 bp}. Hover highlighting remains sample-specific,
#' even when base colors are assigned by group.
#'
#' @examples
#' \dontrun{
#' # From a data frame directly
#' set.seed(1)
#' ex <- expand.grid(index = 1:51, sample = c("INPUT_A1", "H3K4me3_A1"))
#' ex$mean <- c(dnorm(seq(-2.5, 2.5, length.out = 51), sd = 0.8),
#'              dnorm(seq(-2.5, 2.5, length.out = 51), sd = 0.5))
#' p <- plot_enrichment_interactive(df = ex, bin_size = 100, window_bp = 2500)
#' p
#'
#' # From an EPK object
#' p <- plot_enrichment_interactive(
#'   epk    = epk,
#'   marker = "H3K4me3",
#'   loci   = "protein_coding"
#' )
#' p
#' }
#'
#' @importFrom dplyr mutate distinct left_join
#' @importFrom plotly plot_ly layout highlight highlight_key
#' @importFrom grDevices hcl.colors rgb
#' @importFrom stats setNames
#' @export
plot_enrichment_interactive <- function(
    df = NULL,
    epk = NULL,
    marker = NULL,
    loci = NULL,
    bin_size = 100,
    window_bp = 2500,
    highlight_opacity = 0.5,
    mid_coord = c("center", "start"),
    plot_title = "Enrichment Profile",
    control_samples = NULL,
    control_pattern = NULL,
    control_color = "rgba(0,0,0,0.90)",
    other_color = "rgba(120,180,255,0.55)",
    hover_color = "rgba(0,120,255,1)",
    group_by = NULL,
    metadata = NULL,
    metadata_sample_col = "sample",
    group_palette = NULL,
    control_overrides_group = FALSE
) {
  # ── Resolve df from epk + marker + loci ────────────────────────────────────
  if (is.null(df)) {
    if (is.null(epk) || is.null(marker) || is.null(loci)) {
      stop(
        "Provide either 'df' directly, or 'epk', 'marker', and 'loci' together."
      )
    }
    profile_list <- epk$enrichment_results$enrichment_profile[[loci]]
    if (is.null(profile_list)) {
      stop(sprintf(
        "No enrichment_profile found for loci '%s' in epk. ",
        loci
      ))
    }
    df <- profile_list[[marker]]
    if (is.null(df)) {
      stop(sprintf(
        "No profile data found for marker '%s' under loci '%s'.", marker, loci
      ))
    }
    if (identical(plot_title, "Enrichment Profile")) {
      plot_title <- paste0(marker, " — ", loci)
    }
  }

  # ── Validation ──────────────────────────────────────────────────────────────
  if (!is.data.frame(df)) {
    stop("'df' must be a data.frame. Got: ", paste(class(df), collapse = "/"))
  }

  required_cols <- c("index", "mean", "sample")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("'df' is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.numeric(bin_size) || length(bin_size) != 1 || is.na(bin_size) || bin_size <= 0) {
    stop("'bin_size' must be a single positive numeric value.")
  }
  if (!is.numeric(window_bp) || length(window_bp) != 1 || is.na(window_bp) || window_bp <= 0) {
    stop("'window_bp' must be a single positive numeric value.")
  }
  if (!is.numeric(highlight_opacity) || length(highlight_opacity) != 1 ||
      is.na(highlight_opacity) || highlight_opacity < 0 || highlight_opacity > 1) {
    stop("'highlight_opacity' must be a single numeric value in [0, 1].")
  }
  if (!is.null(group_by) && (!is.character(group_by) || length(group_by) != 1)) {
    stop("'group_by' must be NULL or a single column name.")
  }
  if (!is.character(metadata_sample_col) || length(metadata_sample_col) != 1) {
    stop("'metadata_sample_col' must be a single column name.")
  }
  if (!is.logical(control_overrides_group) || length(control_overrides_group) != 1) {
    stop("'control_overrides_group' must be TRUE or FALSE.")
  }

  mid_coord <- match.arg(mid_coord)

  .normalize_plotly_color <- function(x) {
    x <- as.character(x)
    out <- x
    is_rgba <- grepl("^rgba\\(", x)

    if (any(is_rgba)) {
      parsed <- regexec(
        "^rgba\\((\\d+),(\\d+),(\\d+),([0-9]*\\.?[0-9]+)\\)$",
        gsub("\\s+", "", x[is_rgba])
      )
      m <- regmatches(gsub("\\s+", "", x[is_rgba]), parsed)

      out_rgba <- vapply(m, function(v) {
        if (length(v) != 5) return(NA_character_)
        r <- as.numeric(v[2])
        g <- as.numeric(v[3])
        b <- as.numeric(v[4])
        a <- as.numeric(v[5])
        if (any(is.na(c(r, g, b, a)))) return(NA_character_)
        grDevices::rgb(r, g, b, alpha = a * 255, maxColorValue = 255)
      }, character(1))

      out[is_rgba] <- out_rgba
    }

    out
  }

  control_color <- .normalize_plotly_color(control_color)
  other_color   <- .normalize_plotly_color(other_color)
  hover_color   <- .normalize_plotly_color(hover_color)

  df2 <- dplyr::mutate(
    df,
    index  = suppressWarnings(as.numeric(index)),
    mean   = suppressWarnings(as.numeric(mean)),
    sample = trimws(as.character(sample))
  )

  if (anyNA(df2$index)) stop("'index' contains non-numeric values (or NA after coercion).")
  if (anyNA(df2$mean))  stop("'mean' contains non-numeric values (or NA after coercion).")
  if (any(df2$sample == "")) stop("'sample' contains empty values.")

  if (!is.null(metadata)) {
    if (!is.data.frame(metadata)) stop("'metadata' must be a data.frame when provided.")
    if (!metadata_sample_col %in% names(metadata)) {
      stop("'metadata' is missing join column: ", metadata_sample_col)
    }

    md <- as.data.frame(metadata, stringsAsFactors = FALSE)
    md[[metadata_sample_col]] <- trimws(as.character(md[[metadata_sample_col]]))

    if (any(md[[metadata_sample_col]] == "")) {
      stop("'metadata' contains empty sample identifiers in column '", metadata_sample_col, "'.")
    }
    if (anyDuplicated(md[[metadata_sample_col]]) > 0) {
      stop("'metadata' contains duplicated sample identifiers in column '", metadata_sample_col, "'.")
    }

    idx <- match(df2$sample, md[[metadata_sample_col]])
    md_cols <- setdiff(names(md), metadata_sample_col)
    for (cn in md_cols) df2[[cn]] <- md[[cn]][idx]
  }

  is_control_vec <- rep(FALSE, nrow(df2))
  if (!is.null(control_samples)) {
    is_control_vec <- df2$sample %in% as.character(control_samples)
  } else if (!is.null(control_pattern)) {
    is_control_vec <- grepl(control_pattern, df2$sample)
  }

  df2 <- dplyr::mutate(
    df2,
    is_control = is_control_vec,
    pos_bp     = -window_bp + (index - 1) * bin_size,
    sample     = factor(sample)
  )

  sample_info <- dplyr::distinct(df2, sample, is_control)

  if (!is.null(group_by)) {
    if (!group_by %in% names(df2)) stop("'group_by' column not found in data: ", group_by)

    sample_group <- dplyr::distinct(df2, sample, group = .data[[group_by]])
    sample_group$group <- as.character(sample_group$group)
    sample_group$group[is.na(sample_group$group) | trimws(sample_group$group) == ""] <- "Unknown"

    if (any(duplicated(sample_group$sample))) {
      stop("Each sample must map to exactly one group in column '", group_by, "'.")
    }

    grp_levels <- sort(unique(sample_group$group))

    if (is.null(group_palette)) {
      grp_colors <- grDevices::hcl.colors(length(grp_levels), palette = "Set 2")
      group_color_map <- stats::setNames(grp_colors, grp_levels)
    } else {
      if (!is.character(group_palette)) stop("'group_palette' must be a character vector of colors.")
      group_palette <- .normalize_plotly_color(group_palette)
      if (is.null(names(group_palette))) {
        if (length(group_palette) < length(grp_levels)) {
          stop("Unnamed 'group_palette' must have at least as many colors as unique groups.")
        }
        group_color_map <- stats::setNames(group_palette[seq_along(grp_levels)], grp_levels)
      } else {
        group_color_map <- group_palette
        missing_groups <- setdiff(grp_levels, names(group_color_map))
        extra_groups   <- setdiff(names(group_color_map), grp_levels)

        if (length(extra_groups) > 0) {
          warning(
            "'group_palette' contains names not present in observed groups: ",
            paste(extra_groups, collapse = ", "), ". They will be ignored."
          )
          group_color_map <- group_color_map[setdiff(names(group_color_map), extra_groups)]
        }
        if (length(missing_groups) > 0) {
          warning(
            "'group_palette' is missing colors for group level(s): ",
            paste(missing_groups, collapse = ", "), ". Auto-generating colors."
          )
          auto_cols       <- grDevices::hcl.colors(length(missing_groups), palette = "Set 2")
          group_color_map <- c(group_color_map, stats::setNames(auto_cols, missing_groups))
        }
      }
    }

    sample_color_map <- stats::setNames(
      group_color_map[sample_group$group],
      as.character(sample_group$sample)
    )
    sample_info <- dplyr::left_join(sample_info, sample_group, by = "sample")

    if (isTRUE(control_overrides_group)) {
      sample_color_map[as.character(sample_info$sample[sample_info$is_control])] <- control_color
    }

    color_map <- sample_color_map
  } else {
    color_map <- stats::setNames(
      ifelse(sample_info$is_control, control_color, other_color),
      as.character(sample_info$sample)
    )
  }

  hover_template <- "<b>%{fullData.name}</b><br>mean: %{y:.3f}<extra></extra>"
  if (!is.null(group_by)) {
    hover_template <- paste0(
      "<b>%{fullData.name}</b><br>", group_by,
      ": %{customdata}<br>mean: %{y:.3f}<extra></extra>"
    )
    df2$..group_plotly <- as.character(df2[[group_by]])
    df2$..group_plotly[is.na(df2$..group_plotly) | trimws(df2$..group_plotly) == ""] <- "Unknown"
  }

  dfk <- plotly::highlight_key(df2, ~sample)

  window_kb   <- formatC(window_bp / 1000, format = "fg", digits = 3)
  x_tick_vals <- c(-window_bp, 0, window_bp)
  x_tick_text <- c(paste0("-", window_kb, " kb"), mid_coord, paste0("+", window_kb, " kb"))

  if (!is.null(group_by)) {
    p <- plotly::plot_ly(
      data = dfk, x = ~pos_bp, y = ~mean,
      type = "scatter", mode = "lines",
      color = ~sample, colors = color_map,
      line = list(width = 1),
      customdata = ~..group_plotly,
      hovertemplate = hover_template
    )
  } else {
    p <- plotly::plot_ly(
      data = dfk, x = ~pos_bp, y = ~mean,
      type = "scatter", mode = "lines",
      color = ~sample, colors = color_map,
      line = list(width = 1),
      hovertemplate = hover_template
    )
  }

  p <- p |>
    plotly::layout(
      showlegend = FALSE,
      title = plot_title,
      xaxis = list(
        title    = paste0("Distance from ", mid_coord, " (bp)"),
        tickmode = "array",
        tickvals = x_tick_vals,
        ticktext = x_tick_text
      ),
      yaxis  = list(title = "Mean signal"),
      shapes = list(list(
        type = "line", x0 = 0, x1 = 0, y0 = 0, y1 = 1,
        yref = "paper",
        line = list(color = "rgba(100,100,100,0.8)", dash = "dash", width = 1)
      ))
    ) |>
    plotly::highlight(
      on         = "plotly_hover",
      off        = "plotly_doubleclick",
      persistent = FALSE,
      opacityDim = highlight_opacity,
      selected   = list(line = list(color = hover_color, width = 3))
    )

  p
}
