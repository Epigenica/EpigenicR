#' Interactive ChromHMM chromatin state heatmap
#'
#' Generates an interactive heatmap of mean RPGC per chromatin state across
#' samples for a given marker and loci set.
#'
#' @param epk An EPK object containing
#'   \code{epk$enrichment_results$chromatin_states[[loci]][[marker]]}.
#' @param marker Character; marker name (e.g. \code{"5mC"}, \code{"H3K4me3"}).
#' @param loci Character; loci key matching a name in
#'   \code{epk$enrichment_results$chromatin_states}
#'   (e.g. \code{"protein_coding"}, \code{"cpg_islands"}).
#'
#' @return A \code{plotly} htmlwidget heatmap.
#'
#' @examples
#' \dontrun{
#' p <- interactive_heatmap_chromhmm(epk, marker = "5mC", loci = "cpg_islands")
#' p
#' }
#'
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom plotly plot_ly layout
#' @export
interactive_heatmap_chromhmm <- function(epk, marker, loci) {
  if (is.null(epk$enrichment_results$chromatin_states[[loci]])) {
    stop(sprintf("No chromatin_states data found for loci '%s' in epk.", loci))
  }

  df_cs <- epk$enrichment_results$chromatin_states[[loci]][[marker]]

  if (is.null(df_cs)) {
    stop(sprintf(
      "No chromatin_states data found for marker '%s' under loci '%s'.", marker, loci
    ))
  }

  heat_df <- df_cs |>
    dplyr::select(sample_id_rep, Chromatin_State, mean_rpgc_val) |>
    tidyr::pivot_wider(names_from = Chromatin_State, values_from = mean_rpgc_val)

  heat_mat <- as.matrix(
    data.frame(heat_df, row.names = heat_df$sample_id_rep)[, -1, drop = FALSE]
  )

  plotly::plot_ly(
    x = colnames(heat_mat),
    y = rownames(heat_mat),
    z = heat_mat,
    type = "heatmap"
  ) |>
    plotly::layout(
      title = paste0("Chromatin State Distribution - ", marker),
      yaxis = list(autorange = "reversed")
    )
}
