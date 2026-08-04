#' Interactive ChromHMM chromatin state heatmap
#'
#' Generates an interactive heatmap of mean RPGC per chromatin state across
#' samples for a given marker and ChromHMM reference annotation.
#'
#' @param epk An EPK object containing
#'   \code{epk$enrichment_results$chromatin_states[[chromhmm_ref]][[marker]]}.
#' @param marker Character; marker name (e.g. \code{"5mC"}, \code{"H3K4me3"}).
#' @param chromhmm_ref Character; ChromHMM annotation key matching a name in
#'   \code{epk$enrichment_results$chromatin_states}. Corresponds to the
#'   ChromHMM reference BED filename (without \code{.bed}), e.g.
#'   \code{"E107_15_coreMarks_hg38lift_mnemonics"}.
#'
#' @return A \code{plotly} htmlwidget heatmap.
#'
#' @examples
#' \dontrun{
#' p <- interactive_heatmap_chromhmm(
#'   epk,
#'   marker      = "H3K4me3",
#'   chromhmm_ref = "E107_15_coreMarks_hg38lift_mnemonics"
#' )
#' p
#' }
#'
#' @importFrom dplyr select
#' @importFrom tidyr pivot_wider
#' @importFrom plotly plot_ly layout
#' @export
interactive_heatmap_chromhmm <- function(epk, marker, chromhmm_ref, show_pooled = FALSE) {
  if (is.null(epk$enrichment_results$chromatin_states[[chromhmm_ref]])) {
    stop(sprintf(
      "No chromatin_states data found for chromhmm_ref '%s' in epk.", chromhmm_ref
    ))
  }

  df_cs <- epk$enrichment_results$chromatin_states[[chromhmm_ref]][[marker]]
  df_cs <- df_cs |> dplyr::filter(show_pooled | !grepl("pooled", sample_id_rep))

  if (is.null(df_cs)) {
    stop(sprintf(
      "No chromatin_states data found for marker '%s' under chromhmm_ref '%s'.",
      marker, chromhmm_ref
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
    type = "heatmap",
    hovertemplate = "Chromatin State: %{x}<br>Sample: %{y}<br>Mean RPGC: %{z}<extra></extra>"
  ) |>
    plotly::layout(
      title = paste0("Chromatin State Distribution - ", marker),
      yaxis = list(autorange = "reversed")
    )
}
