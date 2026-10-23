# UCE Completeness Threshold Evaluator

An R workflow for comparing phylogenetic trees inferred from datasets generated under different ultraconserved element (UCE) completeness thresholds.

The workflow was developed to provide a simple, evidence-based way to evaluate completeness-threshold selection by examining how **tree topology, branch lengths, and node support change across thresholds**. Rather than selecting a completeness threshold a priori, users can compare a series of thresholds and identify the region in which these properties become relatively stable.

The script was developed for the analyses presented in Tyts et al. (in preparation).

## What the script does

For a set of rooted IQ-TREE topologies generated under different UCE completeness thresholds, the script automatically:

- detects the completeness thresholds from tree filenames;
- reads and checks the input trees;
- summarizes SH-aLRT and UFBoot support;
- compares support distributions among thresholds using Kruskal–Wallis tests and Dunn's post-hoc tests with Benjamini–Hochberg correction;
- compares the proportion of strongly supported nodes using Pearson's chi-squared tests;
- produces support-distribution histograms;
- plots the proportion of strongly supported nodes;
- calculates and plots maximum root-to-tip tree height;
- calculates normalized Robinson–Foulds (RF) distances among topologies;
- produces a topology-similarity heatmap; and
- exports all plots, summary tables, and statistical results to a `results/` folder.

ASTRAL trees can optionally be analyzed using the same workflow. Local posterior probability (LPP) support is summarized and compared among thresholds, and RF distances among ASTRAL topologies can also be calculated.

The script additionally contains an optional `comparePhylo()` section for examining the exact topological differences between selected trees.

## Input requirements

**All trees must be rooted before running the script.**

Trees compared using Robinson–Foulds distances must also contain the same set of taxa.

### IQ-TREE files

Place rooted IQ-TREE files in a folder called:

```text
iqtree/
```

Files must be named using the following pattern:

```text
iqtree_50p.tre
iqtree_55p.tre
iqtree_60p.tre
iqtree_65p.tre
...
```

where the number indicates the completeness threshold.

The extensions `.tre`, `.tree`, and `.treefile` are accepted.

Any number or combination of completeness thresholds can be analyzed. Thresholds do not need to occur at regular intervals.

### Optional ASTRAL files

To include ASTRAL analyses, place rooted ASTRAL trees in:

```text
astral/
```

using the same naming convention:

```text
astral_50p.tre
astral_55p.tre
astral_60p.tre
...
```

Set:

```r
include_astral <- TRUE
```

near the beginning of the script.

If ASTRAL analyses are not required, use:

```r
include_astral <- FALSE
```

## Folder structure

A minimal analysis directory should look like this:

```text
UCE-completeness-evaluator/
│
├── UCE_completeness_evaluator.R
│
├── iqtree/
│   ├── iqtree_50p.tre
│   ├── iqtree_55p.tre
│   ├── iqtree_60p.tre
│   └── ...
│
├── astral/                    # optional
│   ├── astral_50p.tre
│   ├── astral_55p.tre
│   ├── astral_60p.tre
│   └── ...
│
└── results/                   # created automatically if absent
```

## R packages

The workflow requires the following R packages:

```text
ggplot2
treeio
dplyr
tidyr
FSA
ape
phangorn
```

Most can be installed from CRAN:

```r
install.packages(c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "FSA",
  "ape",
  "phangorn"
))
```

`treeio` is available through Bioconductor:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("treeio")
```

## Running the analysis

1. Download or clone this repository.
2. Place the rooted trees into the appropriate `iqtree/` and, if applicable, `astral/` folders.
3. Make sure the filenames follow the patterns described above.
4. Open `UCE_completeness_evaluator.R` in RStudio or run it using R.
5. Check the list of detected tree files printed near the beginning of the analysis.
6. Run the complete script.

In an interactive R session, the script displays the detected thresholds and filenames before reading the trees and allows the analysis to be stopped if the files were detected incorrectly.

Threshold values are obtained automatically from filenames, so **the script does not need to be edited when a different number of completeness thresholds is analyzed**.

Overall, the script is thoroughly annotated, allowing you to jump straight to it.

## User settings

A small number of settings can be changed near the beginning of the script.

### Include ASTRAL analyses

```r
include_astral <- TRUE
```

Change to `FALSE` if only IQ-TREE topologies are being compared.

### Strong-support cutoff

```r
strong_support <- 95
```

The default workflow classifies nodes with support ≥95 as strongly supported. This cutoff can be changed if desired.

### IQ-TREE reference for ASTRAL comparison

```r
iqtree_reference_threshold <- 80
```

If ASTRAL analyses are included, one IQ-TREE topology can optionally be added to the ASTRAL RF-distance comparison. For example, setting the value to `80` adds the IQ-TREE topology generated from the 80% completeness dataset.

To compare only ASTRAL topologies, use:

```r
iqtree_reference_threshold <- NA
```

## Output

All output is written automatically to:

```text
results/
```

### Support analyses

For IQ-TREE, outputs include:

- raw SH-aLRT and UFBoot support values;
- support summary statistics;
- Kruskal–Wallis test results;
- Dunn pairwise comparisons with Benjamini–Hochberg correction;
- Dunn-test heatmaps;
- log1p-scaled support histograms;
- Pearson's chi-squared tests comparing strongly and weakly supported nodes; and
- plots showing the proportion of strongly supported nodes.

Equivalent outputs are produced for ASTRAL LPP when ASTRAL analysis is enabled.

### Branch-length comparison

Maximum root-to-tip distance is calculated for each IQ-TREE topology and exported as both:

```text
IQ_TREE_tree_heights.csv
IQ_TREE_tree_height.pdf
```

The plot can be used to evaluate whether changes in overall tree height stabilize as completeness thresholds are lowered.

### Topological comparison

Normalized Robinson–Foulds distances are calculated among trees and exported in:

```text
IQ_TREE_RF_distances.csv
```

and, when applicable:

```text
ASTRAL_RF_distances.csv
```

The output includes both normalized RF distance and the corresponding topology-similarity percentage:

```text
similarity = (1 - normalized RF distance) × 100
```

Thus, **100% indicates identical topologies**, whereas lower values indicate increasing topological differences.

Topology-similarity heatmaps are saved as PDF files.

## Interpreting completeness-threshold comparisons

The workflow is intended to assist with **empirical threshold selection rather than prescribe a universal completeness value**.

We recommend examining the metrics together, particularly:

1. changes in topology across completeness thresholds;
2. changes in maximum root-to-tip tree height; and
3. changes in node-support distributions.

A useful threshold may occur within the region where topology and tree height become relatively stable with the inclusion of additional loci while node support remains sufficiently high. The appropriate threshold is expected to be dataset-specific.

## Optional examination of individual topology changes

Normalized RF distances summarize the magnitude of differences between trees but do not identify which relationships changed.

The final section of the script therefore provides an optional `comparePhylo()` command. Users can manually specify two thresholds of interest, for example:

```r
comparePhylo(
  iq_phylo[["75"]],
  iq_phylo[["80"]],
  plot = TRUE,
  rooted = TRUE,
  cex = 0.5
)
```

This section is optional and is not part of the automated workflow.

IQ-TREE and ASTRAL topologies from the same threshold can also be compared:

```r
comparePhylo(
  iq_phylo[["80"]],
  astral_phylo[["85"]],
  plot = TRUE,
  rooted = TRUE,
  cex = 0.5
)
```

## Citation

If you use this workflow, please cite:

Tyts, V., Schwartz, M.D., Chérot, F., Hoey-Chamberlain, R., Rosas Sánchez, C., Stukel, M., Schroeder, S., Konstantinov, F., Namyatova, A., Forthman, M. and Weirauch, C. (in preparation). *The first UCE-based phylogenies and divergence-time estimates for the superdiverse Mirinae (Heteroptera: Miridae): implications for data completeness and dating strategies*.

Citation information is also provided in `CITATION.cff`.
