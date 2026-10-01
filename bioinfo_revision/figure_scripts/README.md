# `figure_scripts/` — two main-text figures for the revision

Standalone figure scripts written for the v9 revision. Each runs on its own from the
repository root and writes to [`../reports/figures/`](../reports/figures/) in both PDF (for
the manuscript) and PNG at 300 dpi (for quick viewing and drafts).

```
Rscript bioinfo_revision/figure_scripts/make_gtex_noperm_figure.R
Rscript bioinfo_revision/figure_scripts/make_table2_extended_figure.R
```

| script | figure | what it plots | source data |
| --- | --- | --- | --- |
| `make_gtex_noperm_figure.R` | `fig_gtex_noperm.{pdf,png}` | GTEx Whole Blood, MRGN **without** the permutation test — selected PCs, model calls, T1–T2 edge, mediation type | `GTEx/data/WholeBloodmaster_table.RData`, `Permutation_test_analysis/results/no.perm.all.trios.WB.RData`, `GTEx/results/GMAC/gmac.results.tables.combined.RData` |
| `make_table2_extended_figure.R` | `fig_table2_extended_csq.{pdf,png}` | MRGN with CS-q across the four covariate structures — recall and precision, five models plus the T1–T2 edge | `Final_tables/Final_rec_prec_tables.xlsx`, sheet `Table2 Extended` |

Both use the legacy manuscript figure style — `ggthemes::theme_hc()`, `ggpubr::ggarrange`,
bold size-16 panel letters, a shared legend harvested as a grob — so they sit beside
`Manuscript/figures/` output rather than beside the `theme_bw` figures in
`../simulation_results/results_scripts/`.

---

## `make_gtex_noperm_figure.R`

Rebuilds the manuscript's **Figure 3**
(`Manuscript/figures/MF3_GTEx.model.and.t1t2.edge.bargraphs.pdf`, from
`Manuscript/scripts/create_GTEx_figs.R`) panel for panel on the master table's
`MRGN.Inferred.Model.no.perm` column. Figure 3 plots
`MRGN.Inferred.Model.libconf.alpha05`; this column had never been plotted anywhere.

Same palette, same 12 × 10 in page. MRPC-ADDIS and GMAC have no perm / no-perm variants, so
their bars are identical to Figure 3's; they are kept so the two figures can be read side
by side.

**The panel grid is not Figure 3's.** It reads:

| | |
| --- | --- |
| **A** number of selected PCs | **B** the six model classes |
| **C** T1–T2 edge | **D** type of mediation |

The plot objects therefore go into `ggarrange` as `D, A, C, B` while `labels` stays
`c("A","B","C","D")`, so the displayed letters run in reading order. This follows
`Manuscript/scripts/create_main_figs.R:162`, which likewise passes its panels out of order
and keeps the letters sequential. If the letters should instead stay glued to Figure 3's
panel names, change `labels` — not the argument order.

Three things worth knowing before editing it:

- **Panel A uses `MRGN.number.of.PCs`, not `MRGN.libconf.alpha05.number.of.PCs`.** The
  no-perm call was made on the standard confounder set (median 30 PCs, range 6–51), not the
  liberal α = 0.05 set Figure 3's PC histogram describes (median 51). It has to describe
  the set the other three panels were actually run on.
- **Panel C pairs the model label with indicators from the same run.** Figure 3 pairs its
  `libconf.alpha05` labels with the indicator matrix in
  `reg.res.wo.pseudo.list.all.tissues.RData`, which comes from the *permutation* run — the
  label and the six indicators `get.adj.from.class()` refines the adjacency with come from
  two different runs. Here both come from `no.perm.all.trios.WB.RData`: row 14 is the model
  call, rows 1–6 are its indicators. The script asserts row 14 equals the master table
  column, so the two files cannot drift apart unnoticed.
- **It does not rewrite `ST_GTEx_all_trios_master.csv`.** `create_GTEx_figs.R:26` writes
  that supplementary table as a side effect; a figure script should not.

What it should print:

```
  A | selected PCs       : median=30 range=6-51
  B | MRGN no.perm models: M0=712 M1=76 M2=11 M3=2278 M4=152 Other=19
  C | MRGN T1-T2 edge    : absent=3007 present=241
  D | MRGN M1 split      : V->Cis->Trans=43 V->Trans->Cis=33
```

Different numbers mean an input was regenerated, not that the script drifted.

---

## `make_table2_extended_figure.R`

Turns the `Table2 Extended` sheet into two line panels: **A** recall, **B** precision, the
four covariate structures on x, six lines (M0–M4 and the T1–T2 edge). CS-q arm only.

**It reads the workbook rather than re-deriving from `confusion_structures.R`.** The
workbook is what ships with the manuscript, so sourcing the figure from it means figure and
table cannot disagree. The cost is parsing `"0.683 (0.060)"` strings; the assertions cover
it — an unparseable cell is an error, not an `NA`, because a silent `NA` here would produce
a plausible-looking wrong figure.

- **The x order is not the sheet order.** The sheet runs U, U+Z, U+W, U+W+Z; the figure runs
  **U → U+W → U+Z → U+W+Z**, so the two hazards are introduced one at a time and the drop at
  U+Z reads left to right. Blocks are matched **by title substring**, never by row position,
  so reordering the sheet cannot silently mis-label the figure.
- **The T1–T2 edge is drawn black, dashed and thicker.** It scores the edge-present class
  rather than a model label, so it must not read as a sixth model.
- **Both panels are fixed to y ∈ [0, 1]** with identical breaks, so recall and precision are
  directly comparable across the pair. This leaves the lower half of each panel empty; that
  is deliberate — truncating a proportion axis would overstate how much the structure moves
  performance, which is the one thing this figure is for.
- Error bars are ±1 Wald SE, clamped to [0, 1], on a `position_dodge(0.4)` so six series at
  four x positions stay legible.

The result the figure exists to show, from its own run log:

```
  CS-q recall, Null (M0)         : 0.683 -> 0.783 -> 0.367 -> 0.383
  CS-q recall, Cond. indep. (M3) : 0.650 -> 0.650 -> 0.300 -> 0.300
  CS-q precision, T1 - T2 edge   : 0.879 -> 0.893 -> 0.680 -> 0.710
```

Adding the intermediate (`U + W`) costs nothing; adding the common child (`U + Z`) halves
M0 and M3 recall and takes 20 points off edge precision. The common child is a collider, and
`get.conf.trios()` still selects the trio's own `Z` in 65.7% of trios — see
[`../reports/CONFOUNDER_STRUCTURE.md`](../reports/CONFOUNDER_STRUCTURE.md).

### Cross-check

The T1–T2 edge recall row can be checked against a file built by a completely separate path,
`../simulation_results/tables/structure_comparison.csv` (`edge_present_recall`, `arm == "CSq"`):

| structure | CSV | sheet |
| --- | --- | --- |
| `u_only` | 0.8500 | 0.850 |
| `u_w` | 0.8833 | 0.883 |
| `u_z` | 0.7778 | 0.778 |
| `u_w_z` | 0.7889 | 0.789 |

This also confirms the block-to-scenario mapping, since `u_w` and `u_z` are the two that
would be swapped by a mis-parse.
