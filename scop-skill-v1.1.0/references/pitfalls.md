# SCOP pitfalls and fixes

These are real failures encountered while using SCOP 0.9.2 on macOS (R 4.6),
with fixes that are known to work.

## 1. Network plots look blank

**Symptom:** `CCCNetworkPlot(plot_type = "circle")` saved via `pdf(); print(p);
dev.off()` gives a blank page or a file that only contains text.

**Cause:** SCOP's circlize/base network plots return a **`recordedplot`**, not a
ggplot. `print()` on a recordedplot outside its original device does nothing.

**Fix:** evaluate the plotting call **inside** an open device, once per device:

```r
save_scop_baseplot(
  CCCNetworkPlot(srt, plot_type = "circle"),
  "results/ccc/01_circle", 7, 7
)
```

`save_scop_baseplot()` uses `substitute()` + `eval()` so the expression is
re-evaluated for the PDF and the PNG. Do not use `ragg::agg_png` for base
graphics — it can produce a blank/missing PNG; use `grDevices::png()`.

## 2. A helper evaluates a plot only once (R promise caching)

**Symptom:** PDF has content, PNG is blank; or the second device has no plot.

**Cause:** passing a bare expression to a helper evaluates it once in the caller;
the value (a recordedplot) is reused for the second device.

**Fix:** capture with `substitute(expr)` and `eval(expr, envir = parent.frame())`
inside each device. See `save_scop_baseplot()`.

## 3. SCENIC: `cannot unserialize ALTVEC object ... arrow`

**Symptom:** `RunSCENIC(backend = "cpp")` reports hundreds of failed motif
results and `No regulons remain after filtering`.

**Cause:** `cistarget2()` parallelises with psock workers; arrow feather data
contains ALTREP vectors that cannot be serialized to workers.

**Fix:** use fork on macOS/Linux:

```r
RunSCENIC(..., backend = "cpp", cores = 4, parallel_backend = "fork")
```

On Windows, fork is unavailable; use `cores = 1` (slower but works).

## 4. "SCENIC needs Python" misconception

`RunSCENIC(backend = "cpp")` is a **pure R/C++** implementation:
GRNBoost2 → C++ cisTarget → C++ AUCell. It does **not** call `PrepareEnv()` or
`check_python()`. Only `backend = "python"` requires `pyscenic==0.12.1`,
`arboreto==0.1.6`, `ctxcore==0.2.0`, etc.

If you do install pySCENIC, a workspace venv works; `check_python()` may still be
called by the Python backend, so install into the env SCOP resolves or point
`RETICULATE_PYTHON` at the venv.

## 5. RSS is always 1

**Symptom:** every cell type's top regulon has RSS = 1, and the ranking is
identical across cell types.

**Cause:** hand-written Jensen–Shannon code comparing a single regulon's scalar
mean to the global scalar mean: both normalise to 1, JS = 0, RSS = 1.

**Fix:** use SCOP's internal RSS:

```r
r <- SCENICPlot(srt, group.by = "celltype", plot_type = "rss_heatmap",
                top_n = 30, return_data = TRUE)
r$rss_matrix   # regulon x group
r$rank_table   # group, regulon, specificity_score, rank
```

## 6. `GroupHeatmap` clips bottom rows / top labels

**Symptom:** the last row (e.g. `TNFA_SIGNALING_VIA_NFKB`) or the rotated column
names are cut off.

**Cause:** SCOP heatmaps are fixed-size. `width`/`height` describe the **heatmap
body**, not the whole plot; annotations, column names, and legends need extra
device space.

**Fix:** make the device larger than the body, e.g. body `6 x 7.5 in` → device
`15 x 12 in`. If you pass `width`/`height`, set them to fit inside the device.

## 7. `FeatureStatPlot` argument errors

- `add_stat` only accepts `"none"`, `"mean"`, `"median"`; `"wilcox"` errors with
  `'arg' should be one of “none”, “mean”, “median”`.
- Use `comparisons = list(c("A", "B"))`, `pairwise_method = "wilcox.test"`,
  `sig_label = "p.format"` or `"p.signif"` for tests.
- `CellStatPlot(stat.by = gene)` errors: `stat.by` must be a **metadata column**.
  Create `srt$GENE_pos <- as.integer(expr > 0)` first; note `stat_level = "1"`
  may still show a 0/1 stack — use `FeatureStatPlot` for percent-positive unless
  you need `CellStatPlot`'s composition mode.

## 8. `CellStatPlot` percent stack shows 0/1

For percent-positive of a gene, `CellStatPlot(stat.by = "GENE_pos",
stat_type = "percent", position = "stack")` stacks the 0 and 1 levels. If you
need a single positive-percent bar, compute the summary manually or use
`FeatureStatPlot(plot_type = "bar")` with a metadata indicator.

## 9. Sample-level `RunDEtest` failures

- `Design and contrast are supported only for sample-level edgeR testing.` —
  SCOP only accepts explicit `design`/`contrast` for `test.use = "edgeR"`.
- `Explicit edgeR models need at least 2 biological samples per condition` —
  check `sample_col`; it must identify **independent samples** (donor × condition),
  not the condition column. Set `group.by = NULL` for a global comparison.
- `Design metadata 'position' must be constant within each aggregate sample` —
  use a sample ID that is unique per donor × condition, and include donor in
  `design`.
- Contrast must be a **named numeric vector matching every design column**:

```r
con <- c("(Intercept)" = 0, patientPatient2 = 0,
         patientPatient3 = 0, positiondist = 1)
RunDEtest(..., test.use = "edgeR", sample_col = "sample_id",
          condition_col = "position", design = ~ patient + position,
          contrast = con)
```

Low-count genes are filtered by edgeR (`SLC2A4`, `HK2` may disappear). Report
that explicitly rather than treating absent genes as non-significant.

## 10. `RunSCENIC` regulon naming mismatch

SCOP names regulons `TF(+)`; some pipelines (e.g. SCENIC R) use plain `TF`.
Before comparing two runs, strip the suffix:

```r
strip_reg <- function(x) sub("\\(\\+\\)$", "", x)
```

## 11. CellChat installation

If `pak::pkg_install("jinworks/CellChat")` fails with GitHub auth/HTTP2 errors,
install from a tarball:

```bash
curl -L --http1.1 https://github.com/jinworks/CellChat/archive/refs/heads/master.tar.gz -o CellChat.tar.gz
R CMD INSTALL --library=Rlibs CellChat.tar.gz
```

Dependencies may include `ggalluvial`, `sna`, `network`, `ggnetwork`, `collapse`.

## 12. Cell-level p-values are not inference

`FeatureStatPlot(comparisons=...)` and `RunDEtest(test.use="wilcox")` test cells
as if independent. With thousands of cells, p-values are tiny and misleading.
For formal inference use patient/sample-level pseudobulk (3+ donors), or state
clearly that the test is descriptive.

## 13. `presto` missing

`RunCellChat(do.fast = TRUE)` and Seurat's Wilcoxon fast path need `presto`.
Use `do.fast = FALSE` or install `presto` (`devtools::install_github(
"immunogenomics/presto")`).

## 14. `sctransform` / layer names

SCOP's `layer` argument expects Seurat v5 layers (`counts`, `data`) or v4 slots.
Use `scop_get_layer()` if your Seurat version differs.
