#!/usr/bin/env Rscript

# Reproducible bulk RNA-seq analysis for CRISPR-targeted macrophages
#
# Run from the repository root:
#   Rscript analysis.R
#
# Expected input:
#   data/gene_count.xls
#
# Despite its extension, gene_count.xls is expected to be a tab-delimited text
# file with gene_id, gene_name, and the 28 sample columns listed below.

options(stringsAsFactors = FALSE)
set.seed(123)

# -----------------------------------------------------------------------------
# 1. Configuration and package checks
# -----------------------------------------------------------------------------

input_file <- file.path("data", "gene_count.xls")
results_dir <- "results"
figures_dir <- "figures"

de_dir <- file.path(results_dir, "deseq2")
pca_table_dir <- file.path(results_dir, "pca")
gsea_dir <- file.path(results_dir, "hallmark_gsea")
progeny_dir <- file.path(results_dir, "progeny")
pca_figure_dir <- file.path(figures_dir, "pca")
gsea_figure_dir <- file.path(figures_dir, "hallmark_gsea")
progeny_figure_dir <- file.path(figures_dir, "progeny")

output_dirs <- c(
  de_dir, pca_table_dir, gsea_dir, progeny_dir,
  pca_figure_dir, gsea_figure_dir, progeny_figure_dir
)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

required_packages <- c(
  "DESeq2", "limma", "edgeR", "decoupleR", "clusterProfiler",
  "msigdbr", "dplyr", "tidyr", "tibble", "ggplot2", "ggrepel",
  "stringr", "readr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required packages: ", paste(missing_packages, collapse = ", "),
    ". See README.md for installation instructions.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

sample_cols <- c(
  "HD164_NOEP", "HD164_NTC", "HD164_CLEC4E", "HD164_CLEC6A",
  "HD164_CD209", "HD164_SELL", "HD164_SIG9",
  "HD180_NOEP", "HD180_NTC", "HD180_CLEC4E", "HD180_CLEC6A",
  "HD180_CD209", "HD180_SELL", "HD180_SIG9",
  "HD183_NOEP", "HD183_NTC", "HD183_CLEC4E", "HD183_CLEC6A",
  "HD183_CD209", "HD183_SELL", "HD183_SIG9",
  "HD181_NOEP", "HD181_NTC", "HD181_CLEC4E", "HD181_CLEC6A",
  "HD181_CD209", "HD181_SELL", "HD181_SIG9"
)

condition_levels <- c(
  "NOEP", "NTC", "CLEC4E", "CLEC6A", "CD209", "SELL", "SIGLEC9"
)
ko_conditions <- c("CLEC4E", "CLEC6A", "CD209", "SELL", "SIGLEC9")

# Raw count columns use SIG9; biological labels and all reported results use
# SIGLEC9. This map is the only place where that naming difference is resolved.
sample_metadata <- tibble(
  sample = sample_cols,
  donor = sub("_.*$", "", sample_cols),
  raw_condition = sub("^[^_]+_", "", sample_cols),
  condition = if_else(raw_condition == "SIG9", "SIGLEC9", raw_condition)
) %>%
  mutate(
    donor = factor(donor, levels = c("HD164", "HD180", "HD183", "HD181")),
    condition = factor(condition, levels = condition_levels)
  )

coldata <- as.data.frame(sample_metadata)
rownames(coldata) <- coldata$sample

comparison_table <- tribble(
  ~comparison,       ~numerator, ~denominator,
  "NTC_vs_NOEP",     "NTC",      "NOEP",
  "CLEC4E_vs_NTC",   "CLEC4E",   "NTC",
  "CLEC6A_vs_NTC",   "CLEC6A",   "NTC",
  "CD209_vs_NTC",    "CD209",    "NTC",
  "SELL_vs_NTC",     "SELL",     "NTC",
  "SIGLEC9_vs_NTC",  "SIGLEC9",  "NTC"
)

# -----------------------------------------------------------------------------
# 2. Count import and validation
# -----------------------------------------------------------------------------

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file, call. = FALSE)
}

counts_raw <- read.delim(
  input_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

required_columns <- c("gene_id", "gene_name", sample_cols)
missing_columns <- setdiff(required_columns, colnames(counts_raw))
if (length(missing_columns) > 0L) {
  stop(
    "The count file is missing required columns: ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (anyNA(counts_raw[, sample_cols, drop = FALSE])) {
  stop("The count matrix contains missing values.", call. = FALSE)
}
if (anyDuplicated(counts_raw$gene_id)) {
  stop("gene_id values must be unique.", call. = FALSE)
}

count_mat <- as.matrix(counts_raw[, sample_cols, drop = FALSE])
storage.mode(count_mat) <- "numeric"
rownames(count_mat) <- counts_raw$gene_id

if (any(!is.finite(count_mat)) || any(count_mat < 0)) {
  stop("Counts must be finite and non-negative.", call. = FALSE)
}
if (any(abs(count_mat - round(count_mat)) > sqrt(.Machine$double.eps))) {
  stop("DESeq2 requires integer-valued raw counts.", call. = FALSE)
}
storage.mode(count_mat) <- "integer"

stopifnot(
  identical(colnames(count_mat), rownames(coldata)),
  ncol(count_mat) == 28L,
  all(table(coldata$donor) == 7L),
  all(table(coldata$condition) == 4L)
)

gene_annotation <- counts_raw %>%
  select(gene_id, gene_name) %>%
  distinct(gene_id, .keep_all = TRUE)

readr::write_csv(sample_metadata, file.path(results_dir, "sample_metadata.csv"))

# -----------------------------------------------------------------------------
# 3. Paired differential expression with DESeq2
# -----------------------------------------------------------------------------

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = coldata,
  design = ~ donor + condition
)

# Retain genes with at least 10 counts in at least four samples (one donor set).
keep_deseq <- rowSums(DESeq2::counts(dds) >= 10L) >= 4L
dds <- dds[keep_deseq, ]
if (nrow(dds) == 0L) stop("No genes passed the count filter.", call. = FALSE)

dds <- DESeq2::DESeq(dds, quiet = TRUE)

extract_deseq_results <- function(numerator, denominator, comparison) {
  result <- DESeq2::results(
    dds,
    contrast = c("condition", numerator, denominator),
    alpha = 0.10
  )

  table <- as.data.frame(result) %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_annotation, by = "gene_id") %>%
    relocate(gene_name, gene_id) %>%
    mutate(
      comparison = comparison,
      numerator = numerator,
      denominator = denominator,
      neglog10_pvalue = -log10(pmax(pvalue, .Machine$double.xmin)),
      neglog10_padj = -log10(pmax(padj, .Machine$double.xmin)),
      .before = 1
    ) %>%
    arrange(padj, pvalue)

  readr::write_csv(table, file.path(de_dir, paste0(comparison, ".csv")))
  table
}

deseq_tables <- lapply(seq_len(nrow(comparison_table)), function(i) {
  extract_deseq_results(
    comparison_table$numerator[i],
    comparison_table$denominator[i],
    comparison_table$comparison[i]
  )
})
names(deseq_tables) <- comparison_table$comparison

deseq_all <- bind_rows(deseq_tables)
readr::write_csv(deseq_all, file.path(de_dir, "all_comparisons.csv"))

# -----------------------------------------------------------------------------
# 4. Global PCA and donor-adjusted pairwise PCA
# -----------------------------------------------------------------------------

vsd <- DESeq2::vst(dds, blind = FALSE)
vst_matrix <- SummarizedExperiment::assay(vsd)

# plotPCA uses the 500 most variable genes by default; reproduce that selection
# explicitly so the coordinates and variance explained can also be exported.
top_variable_rows <- function(matrix, n = 500L) {
  gene_variance <- apply(matrix, 1L, stats::var)
  order(gene_variance, decreasing = TRUE)[seq_len(min(n, length(gene_variance)))]
}

global_selected <- top_variable_rows(vst_matrix, 500L)
global_pca <- prcomp(t(vst_matrix[global_selected, , drop = FALSE]))
global_variance <- 100 * global_pca$sdev^2 / sum(global_pca$sdev^2)

global_pca_data <- as.data.frame(global_pca$x[, 1:2, drop = FALSE]) %>%
  rownames_to_column("sample") %>%
  left_join(sample_metadata %>% mutate(across(where(is.factor), as.character)), by = "sample")

readr::write_csv(global_pca_data, file.path(pca_table_dir, "global_pca_coordinates.csv"))

global_pca_plot <- ggplot(
  global_pca_data,
  aes(PC1, PC2, color = donor, shape = condition, label = sample)
) +
  geom_point(size = 3.5) +
  ggrepel::geom_text_repel(size = 2.8, show.legend = FALSE, max.overlaps = Inf) +
  scale_shape_manual(values = setNames(0:6, condition_levels)) +
  labs(
    x = sprintf("PC1: %.1f%% variance", global_variance[1]),
    y = sprintf("PC2: %.1f%% variance", global_variance[2]),
    color = "Donor", shape = "Condition"
  ) +
  theme_classic()

ggsave(
  file.path(pca_figure_dir, "global_PCA.pdf"),
  global_pca_plot, width = 7.2, height = 5.8, units = "in"
)
ggsave(
  file.path(pca_figure_dir, "global_PCA.png"),
  global_pca_plot, width = 7.2, height = 5.8, units = "in", dpi = 600
)

run_pairwise_pca <- function(target_condition, ntop = 500L) {
  keep_samples <- coldata$condition %in% c("NTC", target_condition)
  pair_matrix <- vst_matrix[, keep_samples, drop = FALSE]
  pair_meta <- droplevels(coldata[keep_samples, , drop = FALSE])
  selected <- top_variable_rows(pair_matrix, ntop)

  # Remove donor-associated variation for visualization while preserving the
  # condition term. Differential expression itself remains donor-adjusted in
  # the DESeq2 model above.
  condition_design <- model.matrix(~ condition, data = pair_meta)
  adjusted_matrix <- limma::removeBatchEffect(
    pair_matrix[selected, , drop = FALSE],
    batch = pair_meta$donor,
    design = condition_design
  )

  pca <- prcomp(t(adjusted_matrix))
  variance <- 100 * pca$sdev^2 / sum(pca$sdev^2)
  pca_data <- as.data.frame(pca$x[, 1:2, drop = FALSE]) %>%
    rownames_to_column("sample") %>%
    left_join(
      sample_metadata %>% mutate(across(where(is.factor), as.character)),
      by = "sample"
    )

  readr::write_csv(
    pca_data,
    file.path(pca_table_dir, paste0(target_condition, "_vs_NTC_coordinates.csv"))
  )

  plot <- ggplot(
    pca_data,
    aes(PC1, PC2, color = donor, shape = condition)
  ) +
    stat_ellipse(
      aes(group = condition, fill = condition),
      geom = "polygon", type = "norm", level = 0.95,
      alpha = 0.12, colour = NA, show.legend = FALSE
    ) +
    geom_point(size = 4) +
    ggrepel::geom_text_repel(
      aes(label = sample), size = 3.2, show.legend = FALSE, max.overlaps = Inf
    ) +
    labs(
      title = target_condition,
      x = sprintf("PC1: %.1f%% variance", variance[1]),
      y = sprintf("PC2: %.1f%% variance", variance[2]),
      color = "Donor", shape = "Condition"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(face = "bold.italic", hjust = 0.5),
      legend.position = "none"
    )

  ggsave(
    file.path(pca_figure_dir, paste0(target_condition, "_vs_NTC_donor_adjusted_PCA.pdf")),
    plot, width = 6, height = 6.2, units = "in"
  )
  ggsave(
    file.path(pca_figure_dir, paste0(target_condition, "_vs_NTC_donor_adjusted_PCA.png")),
    plot, width = 6, height = 6.2, units = "in", dpi = 600
  )

  list(data = pca_data, plot = plot)
}

pairwise_pca <- setNames(lapply(ko_conditions, run_pairwise_pca), ko_conditions)

# -----------------------------------------------------------------------------
# 5. Hallmark gene-set enrichment analysis
# -----------------------------------------------------------------------------

msig_hallmark <- tryCatch(
  msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
  error = function(e) msigdbr::msigdbr(species = "Homo sapiens", category = "H")
)

symbol_column <- if ("gene_symbol" %in% colnames(msig_hallmark)) {
  "gene_symbol"
} else if ("human_gene_symbol" %in% colnames(msig_hallmark)) {
  "human_gene_symbol"
} else {
  stop("No gene-symbol column was found in the msigdbr Hallmark data.", call. = FALSE)
}

hallmark_term_to_gene <- msig_hallmark %>%
  select(gs_name, gene_symbol = all_of(symbol_column)) %>%
  distinct()

condition_label <- function(x) {
  ifelse(x %in% ko_conditions, paste0(x, " KO"), x)
}

run_hallmark_gsea <- function(numerator, denominator, comparison) {
  result <- DESeq2::results(
    dds,
    contrast = c("condition", numerator, denominator),
    alpha = 0.10
  )

  rank_table <- as.data.frame(result) %>%
    rownames_to_column("gene_id") %>%
    left_join(gene_annotation, by = "gene_id") %>%
    filter(!is.na(gene_name), gene_name != "", gene_name != "-", !is.na(stat)) %>%
    arrange(desc(abs(stat))) %>%
    distinct(gene_name, .keep_all = TRUE)

  gene_rank <- rank_table$stat
  names(gene_rank) <- rank_table$gene_name
  gene_rank <- sort(gene_rank, decreasing = TRUE)

  set.seed(123)
  enrichment <- clusterProfiler::GSEA(
    geneList = gene_rank,
    TERM2GENE = hallmark_term_to_gene,
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    eps = 0,
    verbose = FALSE,
    seed = TRUE
  )

  enrichment_table <- as.data.frame(enrichment) %>%
    mutate(
      comparison = comparison,
      numerator = numerator,
      denominator = denominator,
      direction = case_when(
        NES > 0 ~ paste0("Enriched in ", condition_label(numerator)),
        NES < 0 ~ paste0("Enriched in ", condition_label(denominator)),
        TRUE ~ "No direction"
      ),
      description_clean = Description %>%
        stringr::str_replace("^HALLMARK_", "") %>%
        stringr::str_replace_all("_", " "),
      .before = 1
    ) %>%
    arrange(p.adjust)

  readr::write_csv(
    rank_table,
    file.path(gsea_dir, paste0(comparison, "_ranked_genes.csv"))
  )
  readr::write_csv(
    enrichment_table,
    file.path(gsea_dir, paste0(comparison, "_hallmark_GSEA.csv"))
  )

  enrichment_table
}

hallmark_tables <- lapply(seq_len(nrow(comparison_table)), function(i) {
  run_hallmark_gsea(
    comparison_table$numerator[i],
    comparison_table$denominator[i],
    comparison_table$comparison[i]
  )
})
names(hallmark_tables) <- comparison_table$comparison

hallmark_all <- bind_rows(hallmark_tables)
hallmark_significant <- hallmark_all %>%
  filter(p.adjust < 0.05) %>%
  arrange(comparison, p.adjust)

readr::write_csv(hallmark_all, file.path(gsea_dir, "all_hallmark_results.csv"))
readr::write_csv(
  hallmark_significant,
  file.path(gsea_dir, "FDR_0.05_hallmark_results.csv")
)

# Manuscript dot plot: retain pathways significant in at least two KO
# comparisons or FDR-significant with an absolute NES in the top decile.
ko_comparison_names <- paste0(ko_conditions, "_vs_NTC")
ko_hallmark <- hallmark_all %>% filter(comparison %in% ko_comparison_names)

nes_threshold <- quantile(abs(ko_hallmark$NES), probs = 0.90, na.rm = TRUE)
recurrent_pathways <- ko_hallmark %>%
  filter(p.adjust < 0.05) %>%
  group_by(ID) %>%
  summarise(n_significant_conditions = n_distinct(comparison), .groups = "drop") %>%
  filter(n_significant_conditions >= 2L) %>%
  pull(ID)

strong_pathways <- ko_hallmark %>%
  filter(p.adjust < 0.05, abs(NES) >= nes_threshold) %>%
  distinct(ID) %>%
  pull(ID)

selected_pathways <- union(recurrent_pathways, strong_pathways)
if (length(selected_pathways) == 0L) {
  warning("No Hallmark pathways met the manuscript dot-plot selection criteria.")
} else {
  hallmark_plot_data <- ko_hallmark %>%
    filter(ID %in% selected_pathways) %>%
    mutate(
      condition = factor(sub("_vs_NTC$", "", comparison), levels = ko_conditions),
      pathway = stringr::str_to_title(description_clean),
      significant = p.adjust < 0.05,
      neglog10_padj = pmin(-log10(pmax(p.adjust, .Machine$double.xmin)), 10)
    )

  pathway_order <- hallmark_plot_data %>%
    group_by(pathway) %>%
    summarise(
      n_significant_conditions = sum(significant, na.rm = TRUE),
      max_abs_NES = max(abs(NES), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(n_significant_conditions, max_abs_NES) %>%
    pull(pathway)

  hallmark_plot_data <- hallmark_plot_data %>%
    mutate(pathway = factor(pathway, levels = pathway_order))

  readr::write_csv(
    hallmark_plot_data,
    file.path(gsea_dir, "manuscript_dotplot_data.csv")
  )

  nes_limit <- max(abs(hallmark_plot_data$NES), na.rm = TRUE)
  hallmark_dotplot <- ggplot(hallmark_plot_data, aes(condition, pathway)) +
    geom_point(aes(colour = NES, size = neglog10_padj, alpha = significant)) +
    scale_colour_gradient2(
      low = "#5B6BD5", mid = "white", high = "#872B2A",
      midpoint = 0, limits = c(-nes_limit, nes_limit), name = "NES"
    ) +
    scale_size_continuous(
      name = expression(-log[10]("adjusted p-value")), range = c(1.5, 6)
    ) +
    scale_alpha_manual(
      values = c(`TRUE` = 1, `FALSE` = 0.20), name = "FDR < 0.05"
    ) +
    labs(x = NULL, y = NULL) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(
    file.path(gsea_figure_dir, "Hallmark_GSEA_combined_dotplot.pdf"),
    hallmark_dotplot, width = 10, height = 9, units = "in"
  )
  ggsave(
    file.path(gsea_figure_dir, "Hallmark_GSEA_combined_dotplot.png"),
    hallmark_dotplot, width = 10, height = 9, units = "in", dpi = 600
  )
}

# -----------------------------------------------------------------------------
# 6. PROGENy pathway activity and paired limma testing
# -----------------------------------------------------------------------------

# Reproduce the voom-scale expression matrix used as PROGENy input. edgeR is
# used here only for filtering/TMM normalization before voom, not as a second
# differential-expression method.
progeny_design <- model.matrix(~ donor + condition, data = coldata)
colnames(progeny_design)[1] <- "Intercept"

dge <- edgeR::DGEList(counts = count_mat)
keep_voom <- edgeR::filterByExpr(dge, design = progeny_design)
dge <- dge[keep_voom, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")
voom_data <- limma::voom(dge, design = progeny_design, plot = FALSE)

mean_voom_expression <- rowMeans(voom_data$E)
progeny_gene_map <- gene_annotation %>%
  filter(
    gene_id %in% rownames(voom_data$E),
    !is.na(gene_name), gene_name != "", gene_name != "-"
  ) %>%
  mutate(mean_voom_expression = mean_voom_expression[gene_id]) %>%
  arrange(desc(mean_voom_expression)) %>%
  distinct(gene_name, .keep_all = TRUE)

progeny_expression <- voom_data$E[
  match(progeny_gene_map$gene_id, rownames(voom_data$E)),
  , drop = FALSE
]
rownames(progeny_expression) <- progeny_gene_map$gene_name

stopifnot(
  !anyDuplicated(rownames(progeny_expression)),
  identical(colnames(progeny_expression), rownames(coldata))
)

progeny_network <- decoupleR::get_progeny(organism = "human", top = 500)
readr::write_csv(progeny_network, file.path(progeny_dir, "progeny_network.csv"))

progeny_coverage <- progeny_network %>%
  group_by(source) %>%
  summarise(
    model_genes = n_distinct(target),
    detected_genes = n_distinct(target[target %in% rownames(progeny_expression)]),
    coverage_percent = 100 * detected_genes / model_genes,
    .groups = "drop"
  ) %>%
  arrange(coverage_percent)
readr::write_csv(progeny_coverage, file.path(progeny_dir, "model_coverage.csv"))

progeny_activity_raw <- decoupleR::run_mlm(
  mat = progeny_expression,
  network = progeny_network,
  .source = "source",
  .target = "target",
  .mor = "weight",
  minsize = 5
)

progeny_activity_long <- progeny_activity_raw %>%
  transmute(
    pathway = source,
    sample = condition,
    activity = score,
    footprint_p_value = p_value
  ) %>%
  left_join(
    sample_metadata %>%
      transmute(
        sample,
        donor = as.character(donor),
        condition = as.character(condition)
      ),
    by = "sample"
  )

if (anyNA(progeny_activity_long$donor)) {
  stop("PROGENy returned sample names that do not match the metadata.", call. = FALSE)
}
readr::write_csv(
  progeny_activity_long,
  file.path(progeny_dir, "sample_pathway_activities.csv")
)

progeny_activity_matrix <- progeny_activity_long %>%
  select(pathway, sample, activity) %>%
  pivot_wider(names_from = sample, values_from = activity) %>%
  column_to_rownames("pathway") %>%
  as.matrix()
progeny_activity_matrix <- progeny_activity_matrix[, rownames(coldata), drop = FALSE]

stopifnot(identical(colnames(progeny_activity_matrix), rownames(coldata)))

progeny_contrasts <- limma::makeContrasts(
  NTC_vs_NOEP = conditionNTC,
  CLEC4E_vs_NTC = conditionCLEC4E - conditionNTC,
  CLEC6A_vs_NTC = conditionCLEC6A - conditionNTC,
  CD209_vs_NTC = conditionCD209 - conditionNTC,
  SELL_vs_NTC = conditionSELL - conditionNTC,
  SIGLEC9_vs_NTC = conditionSIGLEC9 - conditionNTC,
  levels = progeny_design
)

progeny_fit <- limma::lmFit(progeny_activity_matrix, progeny_design) %>%
  limma::contrasts.fit(progeny_contrasts) %>%
  limma::eBayes()

progeny_tables <- lapply(seq_len(ncol(progeny_contrasts)), function(i) {
  limma::topTable(
    progeny_fit,
    coef = i,
    number = Inf,
    sort.by = "P",
    adjust.method = "BH"
  ) %>%
    rownames_to_column("pathway") %>%
    mutate(comparison = colnames(progeny_contrasts)[i], .before = 1)
})
names(progeny_tables) <- colnames(progeny_contrasts)

progeny_results <- bind_rows(progeny_tables) %>%
  left_join(comparison_table, by = "comparison") %>%
  mutate(
    direction = case_when(
      logFC > 0 ~ paste0("Higher inferred activity in ", condition_label(numerator)),
      logFC < 0 ~ paste0("Higher inferred activity in ", condition_label(denominator)),
      TRUE ~ "No direction"
    ),
    global_FDR = p.adjust(P.Value, method = "BH")
  ) %>%
  relocate(comparison, numerator, denominator, direction, pathway) %>%
  arrange(comparison, adj.P.Val)

readr::write_csv(progeny_results, file.path(progeny_dir, "limma_results.csv"))
readr::write_csv(
  progeny_results %>% filter(P.Value < 0.05),
  file.path(progeny_dir, "nominal_p_0.05_results.csv")
)
readr::write_csv(
  progeny_results %>% filter(adj.P.Val < 0.05),
  file.path(progeny_dir, "FDR_0.05_results.csv")
)

pathway_order <- progeny_results %>%
  group_by(pathway) %>%
  summarise(max_abs_change = max(abs(logFC), na.rm = TRUE), .groups = "drop") %>%
  arrange(max_abs_change) %>%
  pull(pathway)

progeny_plot_data <- progeny_results %>%
  mutate(
    comparison = factor(comparison, levels = comparison_table$comparison),
    pathway = factor(pathway, levels = pathway_order),
    significance_label = if_else(P.Value < 0.05, "*", "")
  )

activity_limit <- max(abs(progeny_plot_data$logFC), na.rm = TRUE)
progeny_heatmap <- ggplot(
  progeny_plot_data,
  aes(comparison, pathway, fill = logFC)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = significance_label), size = 5) +
  scale_fill_gradient2(
    low = "#5B6BD5", mid = "white", high = "#872B2A",
    midpoint = 0, limits = c(-activity_limit, activity_limit)
  ) +
  labs(
    x = NULL, y = NULL, fill = "Activity\nchange",
    caption = "* nominal p < 0.05 (paired limma model)"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(progeny_figure_dir, "PROGENy_pathway_activity.pdf"),
  progeny_heatmap, width = 9, height = 6, units = "in"
)
ggsave(
  file.path(progeny_figure_dir, "PROGENy_pathway_activity.png"),
  progeny_heatmap, width = 9, height = 6, units = "in", dpi = 600
)

# -----------------------------------------------------------------------------
# 7. Reproducibility record
# -----------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()), file.path(results_dir, "sessionInfo.txt"))
saveRDS(
  list(
    dds = dds,
    vst = vsd,
    comparison_table = comparison_table,
    hallmark_results = hallmark_all,
    progeny_results = progeny_results
  ),
  file.path(results_dir, "analysis_objects.rds"),
  compress = "xz"
)

message("Analysis complete. Tables are in 'results/' and figures are in 'figures/'.")
