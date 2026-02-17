#!/usr/bin/env Rscript

library(ggplot2)
library(optparse)
library(glue)
library(jsonlite)
library(fs)
library(stringr)

option_list <- list(
  make_option(
    c("-f", "--file"),
    type = "character",
    default = "../ENCSR230ORT_prostate_ReMap/CTCF_ENCSR230ORT_prostate/CTCF_ENCSR230ORT_prostate_matrix.gz",
    help = "Path to matrix gz file"
  ),
  make_option(
    c("-s", "--span"),
    type = "double",
    default = 0.1,
    help = "Base span value used to derive three spans"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

file <- opts$file
span_val <- opts$span

spans_in_plot = c(span_val, span_val * 2, span_val * 3)

if(length(spans_in_plot) == 3){
  span_labels = c(paste0("span = ", spans_in_plot[1]),
                  paste0("span = ", spans_in_plot[2]),
                  paste0("span = ", spans_in_plot[3])) |> 
    paste0(collapse = ", ")
} else {
  stop("Please provide exactly three span values for the plot.")
}

read_matrix_gz <- function(path) {

  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)

  first_line <- readLines(con, n = 1)
  if (length(first_line) == 0) {
    stop("File is empty: ", path)
  }

  header_json <- sub("^@", "", first_line)
  metadata <- jsonlite::fromJSON(header_json, simplifyVector = TRUE)

  raw <- read.table(
    con,
    header = FALSE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  segment_info <- raw[, 1:6, drop = FALSE]
  segment_info[[1]] <- as.character(segment_info[[1]])
  segment_info[[2]] <- as.numeric(segment_info[[2]])
  segment_info[[3]] <- as.numeric(segment_info[[3]])
  segment_info[[4]] <- as.character(segment_info[[4]])
  segment_info[[5]] <- as.character(segment_info[[5]])
  segment_info[[6]] <- as.character(segment_info[[6]])

  matrix_data <- raw[, 7:ncol(raw), drop = FALSE]
  matrix_data <- apply(matrix_data, 2, as.numeric)
  matrix_data <- as.matrix(matrix_data)

  list(
    metadata = metadata,
    segment_info = segment_info,
    data = matrix_data
  )
}

bin_positions <- function(metadata) {
  upstream <- metadata$upstream
  downstream <- metadata$downstream
  body <- metadata$body
  bin_size <- metadata$`bin size`

  total_span <- upstream + body + downstream
  if (total_span %% bin_size != 0) {
    stop("Total span must be divisible by bin size.")
  }

  start <- -upstream + bin_size / 2
  end <- downstream + body - bin_size / 2

  seq.int(start, end, by = bin_size)
}

matrix_data <- read_matrix_gz(file)

target_name <- file |> fs::path_file() |> stringr::str_extract("[:graph:]+(?=_matrix.gz)")
sample_name <- matrix_data$metadata$sample_labels
title_text <- paste0(target_name, " - ", sample_name)

matrix_data$data |> apply(2, sum, na.rm = TRUE) -> ctotals
matrix_data$data |> apply(2, mean, na.rm = TRUE) -> cmeans
data.frame(
  bp = bin_positions(matrix_data$metadata),
  total = ctotals,
  mean = cmeans
) -> vert_agg_df

loess(mean ~ bp, data = vert_agg_df, span = span_val) -> loess_fit
vert_agg_df$loess_fitted_mean <- loess_fit$fitted

havg = vert_agg_df$mean |> mean(na.rm = TRUE)
zero_predicted = predict(loess_fit, c(0))
score = zero_predicted / havg

mean_caption = paste0("Overall Mean: ", round(havg, 2))

plot_spans <- vert_agg_df |> 
  ggplot(aes(x = bp, y = mean)) +
  geom_line() +
  theme_classic() +
  theme(text = element_text(size = 8)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = title_text,
       x = "Distance to center (bp)",
       caption = glue::glue("{mean_caption}\n{span_labels}"),
       y = "Average signal") + 
  geom_smooth(method = "loess", span = spans_in_plot[1], se = FALSE, color = "blue", linetype = "dashed") +
  geom_smooth(method = "loess", span = spans_in_plot[2], se = FALSE, color = "red", linetype = "dashed") +
  geom_smooth(method = "loess", span = spans_in_plot[3], se = FALSE, color = "purple", linetype = "dashed")

plot_loess <- vert_agg_df |> 
  ggplot(aes(x = bp, y = mean)) +
  geom_line() +
  theme_classic() +
  theme(text = element_text(size = 8), plot.title = element_text(size = 9)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = glue::glue("{title_text}\nscore: {round(score, 2)}"),
       x = "Distance to center (bp)",
       caption = glue::glue("{mean_caption}"),
       y = "Average signal") + 
  scale_x_continuous(breaks = seq(-2000, 2000, by = 1000),
        limits = c(-2000, 2000)) +
  geom_hline(yintercept = havg, linetype = "dotted", color = "gray") +
  geom_line(aes(y = loess_fitted_mean), color = "blue", size = 1)

output_stub <- glue::glue("{target_name}_{sample_name}")

plot_spans_pdf <- paste0(output_stub, "_spans.pdf")
plot_spans_png <- paste0(output_stub, "_spans.png")
plot_loess_pdf <- paste0(output_stub, "_loess.pdf")
plot_loess_png <- paste0(output_stub, "_loess.png")

ggsave(plot_spans_pdf, plot = plot_spans, width = 2, height = 3, units = "in")
ggsave(plot_spans_png, plot = plot_spans, width = 2, height = 3, units = "in", dpi = 300)
ggsave(plot_loess_pdf, plot = plot_loess, width = 2, height = 3, units = "in")
ggsave(plot_loess_png, plot = plot_loess, width = 2, height = 3, units = "in", dpi = 300)

output_json <- paste0(output_stub, "_loess_score.json")
output_payload <- list(
  `loess score` = as.numeric(score),
  score = as.numeric(score),
  file = fs::path_abs(file),
  target = target_name,
  sample = sample_name,
  spans = as.numeric(spans_in_plot),
  span_base = as.numeric(span_val),
  score_at_bp = 0,
  overall_mean = as.numeric(havg),
  predicted_at_center = as.numeric(zero_predicted),
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  metadata = matrix_data$metadata
)

jsonlite::write_json(output_payload, output_json, auto_unbox = TRUE, pretty = TRUE)
