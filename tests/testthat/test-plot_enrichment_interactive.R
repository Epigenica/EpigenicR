test_that("plot_enrichment_interactive sets x-axis label and ticks for center", {
  df <- expand.grid(index = 1:5, sample = c("S1", "S2"), stringsAsFactors = FALSE)
  df$mean <- seq_len(nrow(df)) / 10

  p <- plot_enrichment_interactive(
    df = df,
    bin_size = 100,
    window_bp = 200,
    mid_coord = "center"
  )
  pb <- plotly::plotly_build(p)

  expect_s3_class(p, "plotly")
  expect_identical(pb$x$layout$xaxis$title, "Distance from center (bp)")
  expect_equal(pb$x$layout$xaxis$tickvals, c(-200, 0, 200))
  expect_match(pb$x$layout$xaxis$ticktext[[1]], "0\\.2 kb")
  expect_identical(pb$x$layout$xaxis$ticktext[[2]], "center")
  expect_match(pb$x$layout$xaxis$ticktext[[3]], "0\\.2 kb")
})

test_that("plot_enrichment_interactive sets x-axis label for start", {
  df <- expand.grid(index = 1:5, sample = c("S1", "S2"), stringsAsFactors = FALSE)
  df$mean <- seq_len(nrow(df)) / 10

  p <- plot_enrichment_interactive(
    df = df,
    bin_size = 100,
    window_bp = 200,
    mid_coord = "start"
  )
  pb <- plotly::plotly_build(p)

  expect_identical(pb$x$layout$xaxis$title, "Distance from start (bp)")
  expect_match(pb$x$layout$xaxis$ticktext[[1]], "0\\.2 kb")
  expect_identical(pb$x$layout$xaxis$ticktext[[2]], "start")
  expect_match(pb$x$layout$xaxis$ticktext[[3]], "0\\.2 kb")
})

test_that("plot_enrichment_interactive validates required columns", {
  bad_df <- data.frame(index = 1:3, mean = c(1, 2, 3))
  expect_error(
    plot_enrichment_interactive(bad_df),
    "missing required columns",
    fixed = FALSE
  )
})

test_that("plot_enrichment_interactive supports group_by coloring", {
  df <- expand.grid(index = 1:5, sample = c("S1", "S2"), stringsAsFactors = FALSE)
  df$mean <- seq_len(nrow(df)) / 10
  df$condition <- ifelse(df$sample == "S1", "Ctrl", "Treat")

  p <- plot_enrichment_interactive(
    df = df,
    group_by = "condition",
    group_palette = c(Ctrl = "#000000", Treat = "#ff0000")
  )
  pb <- plotly::plotly_build(p)

  trace_colors <- vapply(
    pb$x$data,
    function(tr) {
      if (!is.null(tr$line$color)) tr$line$color else tr$marker$color
    },
    character(1)
  )
  expect_true(any(grepl("rgba\\(0,0,0", trace_colors)))
  expect_true(any(grepl("rgba\\(255,0,0", trace_colors)))
})

test_that("plot_enrichment_interactive can resolve group from metadata", {
  df <- expand.grid(index = 1:5, sample = c("S1", "S2"), stringsAsFactors = FALSE)
  df$mean <- seq_len(nrow(df)) / 10

  metadata <- data.frame(
    sample_id = c("S1", "S2"),
    batch = c("A", "B"),
    stringsAsFactors = FALSE
  )

  p <- plot_enrichment_interactive(
    df = df,
    metadata = metadata,
    metadata_sample_col = "sample_id",
    group_by = "batch"
  )

  expect_s3_class(p, "plotly")
})
