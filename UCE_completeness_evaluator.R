# ***************************************************************
# UCE COMPLETENESS THRESHOLD EVALUATOR ----
# ***************************************************************
#
#   Compare rooted phylogenetic trees inferred from datasets generated
#   under different UCE completeness thresholds.
# 
# Citation:
# Tyts, V., Schwartz, M.D., Chérot, F., Hoey-Chamberlain, R.,
# Rosas Sánchez, C., Stukel, M., Schroeder, S., Konstantinov, F.,
# Namyatova, A., Forthman, M. and Weirauch, C. (in preparation).
# The first UCE-based phylogenies and divergence-time estimates for the
# superdiverse Mirinae (Heteroptera: Miridae): implications for data
# completeness and dating strategies.
#
# Input:
#   Already-rooted trees placed in the folders below and named as:
#     iqtree/iqtree_50p.tre
#     iqtree/iqtree_55p.tre
#     iqtree/iqtree_60p.tre
#     ...
#
#   Optional ASTRAL trees:
#     astral/astral_50p.tre
#     astral/astral_55p.tre
#     astral/astral_60p.tre
#     ...
#
# Notes:
#   - Trees MUST be rooted before running this script. I avoided coding this
#   as some of the ape versions lose tip names while rooting.
#   - Any number of completeness thresholds can be used.
#   - Threshold values are extracted automatically from filenames, make sure to
#   name them using the pattern above.
#   - All results are written to the results/ folder.
#   - This script has titles for navigation and collapsing finished steps.
#   - The script can be run comfortably in the terminal as well as RStudio.


# ***************************************************************
# 1. USER SETTINGS ----
# ***************************************************************
# Those are general settings for the analyses, change them accordingly.

# Set to TRUE if ASTRAL trees should also be analyzed. If FALSE, the script
# will skip ASTRAL-related lines automatically.
include_astral <- TRUE

# Node-support threshold you consider strongly supported. Is used for summary
# graphs and stats similar to the paper.
strong_support <- 95

# Optional: add one selected IQ-TREE topology (completeness threshold for it) 
# to the ASTRAL RF comparison. This reproduces the comparison of ASTRAL trees 
# against a preferred IQ-TREE topology. Set to NA to compare only ASTRAL trees 
# with one another.
iqtree_reference_threshold <- 50

# Ask the user to confirm the automatically detected input files before trees
# are read. This prompt is used only in an interactive R session.
confirm_input_files <- TRUE

# Input/output folders. Paths are relative to the working directory.
iqtree_dir <- "iqtree"
astral_dir <- "astral"
results_dir <- "results"


# ***************************************************************
# 2. PACKAGES ----
# ***************************************************************
# Load all the packages and check if any need installation

required_packages <- c(
  "ggplot2",
  "treeio",
  "dplyr",
  "tidyr",
  "FSA",
  "ape",
  "phangorn"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "The following required R packages are not installed: ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running the script."
  )
}

library(ggplot2)
library(treeio)
library(dplyr)
library(tidyr)
library(FSA)
library(ape)
library(phangorn)


# ***************************************************************
# 3. HELPER FUNCTIONS ----
# ***************************************************************
# Load all the functions needed for the comparison before proceeding.

# Find tree files and extract completeness thresholds from filenames.
detect_tree_files <- function(folder, prefix) {

  if (!dir.exists(folder)) {
    stop("Input folder not found: ", folder)
  }

  pattern <- paste0("^", prefix, "_[0-9]+p\\.(tre|tree|treefile)$")

  files <- list.files(
    folder,
    pattern = pattern,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(files) < 2) {
    stop(
      "Fewer than two valid tree files were detected in '", folder, "'.\n",
      "Expected filenames such as ", prefix, "_50p.tre, ",
      prefix, "_75p.tre, etc."
    )
  }

  thresholds <- as.numeric(
    sub(
      paste0(".*", prefix, "_([0-9]+)p\\.(tre|tree|treefile)$"),
      "\\1",
      basename(files),
      ignore.case = TRUE
    )
  )

  if (any(is.na(thresholds))) {
    stop("At least one completeness threshold could not be extracted from filenames.")
  }

  if (anyDuplicated(thresholds)) {
    stop(
      "Duplicate completeness thresholds were detected in '", folder, "'.\n",
      "Each threshold should have exactly one tree."
    )
  }

  data.frame(
    threshold = thresholds,
    label = paste0(thresholds, "p"),
    file = files,
    filename = basename(files),
    stringsAsFactors = FALSE
  ) |>
    arrange(threshold)
}


# Print detected files in an easy-to-check format.
print_input_table <- function(x, title) {
  cat("\n", title, "\n", sep = "")
  cat(strrep("-", nchar(title)), "\n", sep = "")
  print(x[, c("threshold", "filename")], row.names = FALSE)
}


# Read each tree once with treeio::read.iqtree().
read_tree_set <- function(input_table) {
  trees <- lapply(input_table$file, treeio::read.iqtree)
  names(trees) <- as.character(input_table$threshold)
  trees
}


# Convert treeio objects to phylo objects without rereading files.
to_phylo_list <- function(tree_list) {
  lapply(tree_list, ape::as.phylo)
}


# Confirm that every tree is rooted.
check_rooting <- function(tree_list, tree_set_name) {
  rooted <- vapply(tree_list, ape::is.rooted, logical(1))

  if (!all(rooted)) {
    stop(
      tree_set_name,
      " contains one or more unrooted trees: ",
      paste(names(rooted)[!rooted], collapse = ", "),
      ".\nRoot all trees before running this workflow."
    )
  }
}


# Confirm that all trees contain the same taxa before RF comparisons.
check_tip_sets <- function(tree_list, tree_set_name) {
  reference <- sort(tree_list[[1]]$tip.label)

  same_tips <- vapply(
    tree_list,
    function(x) identical(sort(x$tip.label), reference),
    logical(1)
  )

  if (!all(same_tips)) {
    stop(
      tree_set_name,
      " trees do not all contain the same set of tip labels.\n",
      "RF comparisons in this workflow require identical taxon sampling."
    )
  }
}


# Extract a named support column from all treeio objects.
extract_support <- function(tree_list, column, support_name, multiplier = 1) {

  out <- lapply(names(tree_list), function(thr) {
    dat <- tree_list[[thr]]@data

    if (!column %in% names(dat)) {
      stop(
        "Support column '", column, "' was not found for threshold ", thr, "p.\n",
        "Available columns are: ", paste(names(dat), collapse = ", ")
      )
    }

    data.frame(
      threshold = as.numeric(thr),
      topology = thr,
      support = support_name,
      value = as.numeric(dat[[column]]) * multiplier
    )
  })

  bind_rows(out) |>
    drop_na(value)
}


# Find an ASTRAL support column automatically.
# read.iqtree() often stores ASTRAL local posterior probabilities in SH_aLRT.
extract_astral_support <- function(tree_list) {

  candidate_columns <- c(
    "LPP", "lpp", "localPP", "localpp",
    "posterior", "SH_aLRT"
  )

  available <- names(tree_list[[1]]@data)
  selected <- candidate_columns[candidate_columns %in% available][1]

  if (is.na(selected)) {
    stop(
      "ASTRAL support values could not be identified automatically.\n",
      "Available columns are: ", paste(available, collapse = ", ")
    )
  }

  raw_values <- unlist(lapply(tree_list, function(x) x@data[[selected]]))

  # ASTRAL LPP is commonly stored from 0 to 1. Convert to percentages if needed.
  multiplier <- if (max(raw_values, na.rm = TRUE) <= 1) 100 else 1

  cat(
    "\nASTRAL support column detected: ", selected,
    if (multiplier == 100) " (converted from 0-1 to 0-100)" else "",
    "\n",
    sep = ""
  )

  extract_support(
    tree_list = tree_list,
    column = selected,
    support_name = "LPP",
    multiplier = multiplier
  )
}


# Create an R-style summary table equivalent to summary() for every threshold
# and support metric.
make_support_summary <- function(df) {
  df |>
    group_by(support, threshold) |>
    summarise(
      n = n(),
      minimum = min(value),
      q1 = unname(quantile(value, 0.25)),
      median = median(value),
      mean = mean(value),
      q3 = unname(quantile(value, 0.75)),
      maximum = max(value),
      .groups = "drop"
    ) |>
    arrange(support, threshold)
}


# Kruskal-Wallis, Dunn post-hoc, and chi-squared tests for one support metric.
run_support_statistics <- function(df, metric, output_prefix) {

  sub <- df |>
    filter(support == metric) |>
    mutate(topology = factor(topology, levels = as.character(sort(unique(threshold)))))

  kw <- kruskal.test(value ~ topology, data = sub)

  dunn <- FSA::dunnTest(
    value ~ topology,
    data = sub,
    method = "bh"
  )$res |>
    tidyr::separate(
      Comparison,
      into = c("tree1", "tree2"),
      sep = " - ",
      remove = FALSE
    )

  chi <- chisq.test(table(sub$topology, sub$value >= strong_support))

  write.csv(
    dunn,
    file.path(results_dir, paste0(output_prefix, "_Dunn_BH.csv")),
    row.names = FALSE
  )

  capture.output(
    kw,
    file = file.path(results_dir, paste0(output_prefix, "_Kruskal_Wallis.txt"))
  )

  capture.output(
    chi,
    file = file.path(results_dir, paste0(output_prefix, "_Chi_squared_strong_support.txt"))
  )

  list(kruskal = kw, dunn = dunn, chisq = chi)
}


# Dunn-test heatmap using the same viridis option as the original workflow.
plot_dunn <- function(dunn_df, filename) {

  tree_levels <- sort(unique(c(
    as.numeric(dunn_df$tree1),
    as.numeric(dunn_df$tree2)
  )))

  p <- dunn_df |>
    mutate(
      tree1 = factor(tree1, levels = as.character(tree_levels)),
      tree2 = factor(tree2, levels = as.character(tree_levels))
    ) |>
    ggplot(aes(tree1, tree2, fill = P.adj)) +
    geom_tile() +
    geom_text(
      aes(label = round(signif(P.adj), digits = 2)),
      color = "white",
      size = 5,
      fontface = 1
    ) +
    scale_fill_viridis_c(option = "G") +
    theme_classic() +
    coord_equal() +
    labs(x = NULL, y = NULL, fill = "Adjusted p")

  n_thr <- length(tree_levels)
  size_cm <- max(12, 2.5 * n_thr)

  ggsave(
    file.path(results_dir, filename),
    p,
    device = "pdf",
    width = size_cm,
    height = size_cm,
    units = "cm"
  )

  p
}


# Support histogram, matching the original log1p count transformation.
plot_support_histogram <- function(df, filename, use_fill = TRUE) {

  df <- df |>
    mutate(
      topology = factor(
        topology,
        levels = as.character(sort(unique(threshold)))
      )
    )

  if (use_fill) {
    p <- ggplot(df, aes(value, fill = support)) +
      geom_histogram(binwidth = 2)
  } else {
    p <- ggplot(df, aes(value)) +
      geom_histogram(binwidth = 2)
  }

  p <- p +
    scale_y_continuous(trans = "log1p") +
    facet_grid(support ~ topology) +
    theme_bw() +
    labs(
      x = "Node support",
      y = "Count (log1p scale)"
    )

  n_thr <- length(unique(df$threshold))

  ggsave(
    file.path(results_dir, filename),
    p,
    device = "pdf",
    width = max(18, 6 * n_thr),
    height = if (length(unique(df$support)) > 1) 12 else 8,
    units = "cm"
  )

  p
}


# Proportion of nodes above/below the strong-support threshold.
plot_strong_support <- function(df, filename) {

  plot_df <- df |>
    mutate(
      topology = factor(
        topology,
        levels = as.character(sort(unique(threshold)))
      ),
      high = value >= strong_support
    )

  p <- ggplot(plot_df, aes(topology, fill = high)) +
    geom_bar(position = "fill") +
    facet_wrap(~support) +
    theme_bw() +
    ylab(paste0("Proportion of >=", strong_support, " nodes")) +
    xlab("Completeness threshold") +
    labs(fill = paste0("Support >=", strong_support))

  n_thr <- length(unique(plot_df$threshold))

  ggsave(
    file.path(results_dir, filename),
    p,
    device = "pdf",
    width = max(13, 3.5 * n_thr),
    height = 9,
    units = "cm"
  )

  p
}


# Calculate tree heights (maximum root-to-tip distance).
calculate_tree_heights <- function(tree_list) {
  heights <- vapply(
    tree_list,
    function(tr) max(ape::node.depth.edgelength(tr)),
    numeric(1)
  )

  data.frame(
    threshold = as.numeric(names(heights)),
    height = as.numeric(heights)
  ) |>
    arrange(threshold)
}


# Plot tree heights using the original steelblue color scheme.
plot_tree_heights <- function(df_height, filename) {

  p <- ggplot(df_height, aes(threshold, height)) +
    geom_line(linewidth = 1, color = "steelblue") +
    geom_point(size = 4, color = "steelblue") +
    theme_bw() +
    ylab("Tree height (root-to-tip distance)") +
    xlab("Completeness threshold")

  ggsave(
    file.path(results_dir, filename),
    p,
    device = "pdf",
    width = 15,
    height = 10,
    units = "cm"
  )

  p
}


# Calculate normalized RF distance and corresponding topology similarity.
calculate_rf <- function(tree_list) {
  
  # Convert the ordinary list of phylo objects to a multiPhylo object
  trees_multi <- tree_list
  class(trees_multi) <- "multiPhylo"
  
  rf_matrix <- phangorn::RF.dist(
    trees_multi,
    normalize = TRUE
  )
  
  rf_df <- as.data.frame(as.table(as.matrix(rf_matrix)))
  names(rf_df) <- c("tree1", "tree2", "RF")
  
  rf_df |>
    dplyr::mutate(
      tree1 = as.character(tree1),
      tree2 = as.character(tree2),
      RF = as.numeric(RF),
      RF_percent = RF * 100,
      similarity_percent = (1 - RF) * 100
    )
}


# Plot topology similarity heatmap.
# This uses 100 for identical trees and 0 for maximally different trees.
plot_rf_similarity <- function(rf_df, tree_order, filename) {

  plot_df <- rf_df |>
    mutate(
      tree1 = factor(tree1, levels = tree_order),
      tree2 = factor(tree2, levels = tree_order),
      shared = round(similarity_percent)
    ) |>
    filter(as.numeric(tree1) <= as.numeric(tree2))

  p <- ggplot(plot_df, aes(tree1, tree2, fill = shared)) +
    geom_tile() +
    geom_text(aes(label = shared), color = "white", size = 5) +
    scale_fill_gradient(
      low = "white",
      high = "black",
      limits = c(0, 100)
    ) +
    coord_equal() +
    labs(x = NULL, y = NULL, fill = "% similarity") +
    theme_classic() +
    theme(axis.ticks = element_blank())

  n_trees <- length(tree_order)
  size_cm <- max(12, 2.5 * n_trees)

  ggsave(
    file.path(results_dir, filename),
    p,
    device = "pdf",
    width = size_cm,
    height = size_cm,
    units = "cm"
  )

  p
}


# ***************************************************************
# 4. CREATE RESULTS FOLDER ----
# ***************************************************************

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}


# ***************************************************************
# 5. DETECT AND CHECK INPUT FILES ----
# ***************************************************************

iq_input <- detect_tree_files(iqtree_dir, "iqtree")
print_input_table(iq_input, "IQ-TREE files detected")

if (include_astral) {
  astral_input <- detect_tree_files(astral_dir, "astral")
  print_input_table(astral_input, "ASTRAL files detected")
}

# Save the detected file list before reading trees.
write.csv(
  iq_input[, c("threshold", "filename")],
  file.path(results_dir, "detected_IQ_TREE_files.csv"),
  row.names = FALSE
)

if (include_astral) {
  write.csv(
    astral_input[, c("threshold", "filename")],
    file.path(results_dir, "detected_ASTRAL_files.csv"),
    row.names = FALSE
  )
}

# Confirmation pause for the user to check that the trees were loaded correctly.
# Introduced for running in terminal or generally the whole script at once.
if (confirm_input_files && interactive()) {
  answer <- readline(
    "\nCheck the files listed above. Press ENTER to continue or type STOP to quit: "
  )

  if (toupper(trimws(answer)) == "STOP") {
    stop("Analysis stopped before reading trees.")
  }
}


# ***************************************************************
# 6. READ TREES ----
# ***************************************************************

iq_trees <- read_tree_set(iq_input)
iq_phylo <- to_phylo_list(iq_trees)

check_rooting(iq_phylo, "IQ-TREE")
check_tip_sets(iq_phylo, "IQ-TREE")

if (include_astral) {
  astral_trees <- read_tree_set(astral_input)
  astral_phylo <- to_phylo_list(astral_trees)

  check_rooting(astral_phylo, "ASTRAL")
  check_tip_sets(astral_phylo, "ASTRAL")
}


# ***************************************************************
# 7. IQ-TREE SUPPORT ANALYSES ----
# ***************************************************************

SH_df <- extract_support(
  iq_trees,
  column = "SH_aLRT",
  support_name = "SH-aLRT"
)

UF_df <- extract_support(
  iq_trees,
  column = "UFboot",
  support_name = "UFboot"
)

iq_support <- bind_rows(SH_df, UF_df) |>
  mutate(topology = as.character(threshold))

write.csv(
  iq_support,
  file.path(results_dir, "IQ_TREE_support_values.csv"),
  row.names = FALSE
)

write.csv(
  make_support_summary(iq_support),
  file.path(results_dir, "IQ_TREE_support_summary.csv"),
  row.names = FALSE
)

# 7a. Kruskal-Wallis, Dunn, and chi-squared tests
SH_stats <- run_support_statistics(
  iq_support,
  metric = "SH-aLRT",
  output_prefix = "IQ_TREE_SH_aLRT"
)

UF_stats <- run_support_statistics(
  iq_support,
  metric = "UFboot",
  output_prefix = "IQ_TREE_UFboot"
)

# 7b. Dunn-test heatmaps
plot_dunn(
  SH_stats$dunn,
  "IQ_TREE_Dunn_SH_aLRT.pdf"
)

plot_dunn(
  UF_stats$dunn,
  "IQ_TREE_Dunn_UFboot.pdf"
)

# 7c. Faceted support histograms
plot_support_histogram(
  iq_support,
  "IQ_TREE_support_histograms_log1p.pdf",
  use_fill = TRUE
)

# 7d. Proportion of strongly supported nodes
plot_strong_support(
  iq_support,
  "IQ_TREE_strong_support_proportions.pdf"
)


# ***************************************************************
# 8. ASTRAL SUPPORT ANALYSES ----
# ***************************************************************

if (include_astral) {

  astral_support <- extract_astral_support(astral_trees) |>
    mutate(topology = as.character(threshold))

  write.csv(
    astral_support,
    file.path(results_dir, "ASTRAL_support_values.csv"),
    row.names = FALSE
  )

  write.csv(
    make_support_summary(astral_support),
    file.path(results_dir, "ASTRAL_support_summary.csv"),
    row.names = FALSE
  )

  # 8a. Kruskal-Wallis, Dunn, and chi-squared tests
  LPP_stats <- run_support_statistics(
    astral_support,
    metric = "LPP",
    output_prefix = "ASTRAL_LPP"
  )

  # 8b. Dunn-test heatmap
  plot_dunn(
    LPP_stats$dunn,
    "ASTRAL_Dunn_LPP.pdf"
  )

  # 8c. Faceted support histograms
  plot_support_histogram(
    astral_support,
    "ASTRAL_support_histograms_log1p.pdf",
    use_fill = FALSE
  )

  # 8d. Proportion of strongly supported nodes
  plot_strong_support(
    astral_support,
    "ASTRAL_strong_support_proportions.pdf"
  )
}


# ***************************************************************
# 9. IQ-TREE TREE-HEIGHT COMPARISON ----
# ***************************************************************

iq_height <- calculate_tree_heights(iq_phylo)

write.csv(
  iq_height,
  file.path(results_dir, "IQ_TREE_tree_heights.csv"),
  row.names = FALSE
)

plot_tree_heights(
  iq_height,
  "IQ_TREE_tree_height.pdf"
)


# ***************************************************************
# 10. IQ-TREE ROBINSON-FOULDS COMPARISON ----
# ***************************************************************

iq_rf <- calculate_rf(iq_phylo)

write.csv(
  iq_rf,
  file.path(results_dir, "IQ_TREE_RF_distances.csv"),
  row.names = FALSE
)

# Keep threshold order numerical rather than alphabetical.
iq_order <- as.character(sort(iq_input$threshold))

plot_rf_similarity(
  iq_rf,
  tree_order = iq_order,
  filename = "IQ_TREE_RF_similarity.pdf"
)


# ***************************************************************
# 11. ASTRAL ROBINSON-FOULDS COMPARISON ----
# ***************************************************************

if (include_astral) {

  astral_rf_trees <- astral_phylo
  names(astral_rf_trees) <- paste0(names(astral_rf_trees), "a")

  # Optionally add the preferred IQ-TREE topology to the ASTRAL comparison,
  # reproducing the comparison used in the original workflow.
  if (!is.na(iqtree_reference_threshold)) {

    reference_name <- as.character(iqtree_reference_threshold)

    if (!reference_name %in% names(iq_phylo)) {
      stop(
        "iqtree_reference_threshold = ", iqtree_reference_threshold,
        " was requested, but no matching IQ-TREE file was detected."
      )
    }

    # Check that ASTRAL and reference IQ-TREE have the same taxa.
    if (!identical(
      sort(astral_rf_trees[[1]]$tip.label),
      sort(iq_phylo[[reference_name]]$tip.label)
    )) {
      stop(
        "The selected IQ-TREE reference and ASTRAL trees do not contain ",
        "the same set of taxa."
      )
    }

    astral_rf_trees[[reference_name]] <- iq_phylo[[reference_name]]
  }

  astral_rf <- calculate_rf(astral_rf_trees)

  write.csv(
    astral_rf,
    file.path(results_dir, "ASTRAL_RF_distances.csv"),
    row.names = FALSE
  )

  astral_order <- paste0(sort(astral_input$threshold), "a")

  if (!is.na(iqtree_reference_threshold)) {
    astral_order <- c(astral_order, as.character(iqtree_reference_threshold))
  }

  plot_rf_similarity(
    astral_rf,
    tree_order = astral_order,
    filename = "ASTRAL_RF_similarity.pdf"
  )
}


# ***************************************************************
# 12. SAVE A COMBINED TEXT SUMMARY OF STATISTICAL TESTS ----
# ***************************************************************

summary_file <- file.path(results_dir, "statistical_tests_summary.txt")

capture.output({

  cat("UCE COMPLETENESS THRESHOLD EVALUATOR\n")
  cat("====================================\n\n")

  cat("Strong-support cutoff:", strong_support, "\n\n")

  cat("IQ-TREE: SH-aLRT\n")
  cat("------------------\n")
  print(SH_stats$kruskal)
  cat("\n")
  print(SH_stats$chisq)
  cat("\n\n")

  cat("IQ-TREE: UFBoot\n")
  cat("----------------\n")
  print(UF_stats$kruskal)
  cat("\n")
  print(UF_stats$chisq)
  cat("\n\n")

  if (include_astral) {
    cat("ASTRAL: LPP\n")
    cat("-----------\n")
    print(LPP_stats$kruskal)
    cat("\n")
    print(LPP_stats$chisq)
    cat("\n")
  }

}, file = summary_file)


# ***************************************************************
# 13. OPTIONAL: VISUAL COMPARISON OF TWO SPECIFIC TOPOLOGIES ----
# ***************************************************************
#
# This section is intentionally manual. It is useful when the RF results show
# that two thresholds differ and the user wants to identify the exact
# topological changes.
#
# Replace the threshold values below with any two thresholds present in the
# corresponding tree list, then uncomment the lines (Ctrl+Shift+C).
#
# IQ-TREE example:
# comparePhylo(
#   iq_phylo[["75"]],
#   iq_phylo[["80"]],
#   plot = TRUE,
#   rooted = TRUE,
#   cex = 0.5
# )
#
# IQ-TREE vs ASTRAL example:
# comparePhylo(
#   iq_phylo[["80"]],
#   astral_phylo[["85"]],
#   plot = TRUE,
#   rooted = TRUE,
#   cex = 0.5
# )


# ***************************************************************
# 14. FINISHED ----
# ***************************************************************

cat(
  "\nAnalysis complete.\n",
  "All tables, statistical-test outputs, and plots have been saved to:\n",
  normalizePath(results_dir),
  "\n",
  sep = ""
)
