# Case study: per-gene EC expression, all cells, two conditions

A worked example of the `run_expression_workflow.R` pattern.

## Question

In endothelial cells (EC), compare **prox vs dist** expression of a small panel:

`PFKFB3`, `SLC2A4`, `SLC2A3`, `HK2`, `PKM`, `LDHA`, `VEGFA`

Requirements:

- use **all EC** (do not split by EC subtype);
- one figure per gene, not a faceted multi-gene panel;
- patient-paired inference, because conditions come from the same donors.

## SCOP functions

```r
FeatureStatPlot(srt, stat.by = gene, group.by = "condition",
                plot_type = "violin",
                comparisons = list(c("prox", "dist")),
                add_stat = "none", sig_label = "p.format")

FeatureDimPlot(srt, features = gene, reduction = "umap")

GroupHeatmap(srt, features = genes, group.by = "condition",
             exp_method = "zscore", limits = c(-2, 2))
```

## Statistical design

- 3 patients × 2 conditions = 6 pseudobulk samples;
- `design = ~ patient + condition`;
- test `dist` vs `prox`;
- primary inference at the **sample** level, not the cell level.

With `pseudobulk_deseq2()`:

```r
stats <- pseudobulk_deseq2(
  srt, genes = genes,
  group_col = "condition", group1 = "dist", group2 = "prox",
  sample_col = "sample_id", covariates = "patient"
)
```

## Typical findings (illustrative)

With only three donors, most genes do not reach significance. The pattern that
recurred in one project:

| gene | direction | interpretation |
|---|---|---|
| PFKFB3 | no difference | bulk EC activity is stable |
| SLC2A3 | lower in dist (trend) | glucose uptake may be reduced |
| VEGFA | lower in dist (trend) | angiogenic output may be reduced |
| LDHA | slightly lower in dist | mild glycolytic shift |
| PKM | no difference | |
| SLC2A4 | almost undetectable (3/983 cells) | GLUT4 is not the EC transporter |
| HK2 | 10/983 cells | too sparse to interpret |

Single-cell Wilcoxon p-values can be very small (e.g. SLC2A3 p = 0.0018,
VEGFA p = 0.0021) but are **descriptive**; pseudobulk p-values were ~0.05–0.2.
Report both, and do not call the cell-level test significant inference.

## Practical notes

1. Filter non-EC contaminants before the analysis.
2. Keep `SLC2A4` in the table but annotate low detection; do not drop silently.
3. Use one SCOP figure per gene (`FeatureStatPlot`, `FeatureDimPlot`) and one
   `GroupHeatmap` for the panel.
4. Save results as PDF + PNG for each figure; keep a `gene_summary.csv` and
   `pseudobulk_DESeq2.csv`.
5. If the user asks for subtype-level results, run a separate workflow with
   `group.by = "ec_subtype"`; do not mix it into the all-EC figures.
