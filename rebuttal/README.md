# Re-analysis of Figs 1–2 with pair-respecting permutations

This directory holds a Snakemake workflow that re-runs the simulations behind
Figs 1 and 2 (and S1–S3) of Dos Santos, Murariu *et al.* (2026, PLoS Comput
Biol) under different permutation schemes. The rest of the repository is left
untouched; the workflow reads `data/` and reuses `code/get_confusion.R`.

## Why

TCGA (tumour/normal from the same patient), the PD1 immunotherapy dataset
(pre/on-treatment from the same patient) and the single-cell pseudobulk
(CTRL/STIM from the same donor × cell type) are paired. The paper's null data
come from `seqgendiff::thin_2group()`, which assigns a fresh random 50/50 split
ignoring pairs, and all tools are fitted with unpaired models.

| dataset | samples | complete pairs | unpaired |
|---|---|---|---|
| TCGA BRCA / KIRC / LIHC / LUAD / PRAD / THCA | 224 / 144 / 100 / 114 / 104 / 116 | all | 0 |
| PD1 immunotherapy (corrected) | 108 | 43 | 22 |
| single-cell pseudobulk (8 donors × 13 cell types) | 206 | 102 | 2 |
| yeast | 86 | 0 | 86 |

## Schemes (`config/config.yaml`)

| scheme | data | permutation | model |
|---|---|---|---|
| `published` | as in the paper's scripts | `thin_2group` (paper) | `~ group` |
| `pairsplit` | corrected | each pair split across groups (coin flip per pair) | `~ group` |
| `pairsplit_balanced` | corrected | pairs split *and* exactly half flipped, so each new group holds n/2 of each original condition | `~ group` |
| `paired_model` | corrected | pairsplit | `~ pair + group`, complete pairs only |

`run_schemes` selects which run (default: `published`, `pairsplit`).

For a given iteration the same genes get the same simulated effects in every
scheme (the signal is drawn in the same RNG order as `thin_2group`); only the
grouping differs, so schemes can be compared iteration by iteration.

**Data variants.** `published` reproduces the paper's loading exactly,
including two problems with the PD1 data: `read.table(skip = 35)` on a file
whose header is on line 1 (drops 35 genes and uses gene 36 as the header —
the paper's confusion matrices, 20,478 genes, match this), and Pt109's two
"On" libraries. `corrected` reads the file properly and drops the second
Pt109 library.

**Tools.** `deseq`, `edger` (paper: `glmQLFit` on the raw matrix — this is
already edgeR v4 QL because `legacy = FALSE` is the default in edgeR ≥ 4, but
without TMM), `edger_ql` (new: `DGEList` → `normLibSizes` → `glmQLFit(legacy =
FALSE)` → `glmQLFTest`), `limma` (voom), `aldex3` (γ = 0.001, 0.1 … 0.5). ALDEx2
is off by default (`enable_aldex2: true` to add it; ~23 min and ~45 GB per
BRCA call).

As in the paper, every iteration *i* runs `set.seed(20 + i)` → thinning → tool
on the original data with permuted labels (`.r`) → tool on the thinned data
(`.t`), and one run on the original labels after `set.seed(2025)` (`.u`). All
three are kept. Figs 1/2 use only `.t`.

## Running

Build the container (once):

```bash
apptainer build rebuttal/container/usri-r.sif rebuttal/container/usri-r.def
```

Then, **from the repository root**:

```bash
uv run --project rebuttal snakemake -s rebuttal/Snakefile -n
```

```bash
uv run --project rebuttal snakemake -s rebuttal/Snakefile --profile rebuttal/profiles/slurm
```

Edit `profiles/slurm/config.yaml` (partition/account) first. Running from the
repository root matters: Apptainer binds the current directory, and all paths
are relative to it.

Before the full run, check the ALDEx3 cost on your cluster with
`n_iter: 1` (e.g. `--config n_iter=1 chunk_size=1`); the
`results/benchmarks/` files record wall time and max RSS per job, and
`resources:` in the config should be adjusted to match.

## Outputs (`results/`)

- `data/{published,corrected}/<dataset>.rds` — filtered counts, labels, pairs
- `<scheme>/runs/<dataset>/<tool>/chunk_*.rds` — per-iteration `.r`/`.t` tool
  tables plus `coefmat`/`designmat` (the thinned counts are not stored);
  `unpermuted.rds` — the `.u` run. `load_run()` in `R/common.R` reassembles
  these into the paper's object format.
- `<scheme>/confusion/<dataset>.<tool>.rds` — `get_confusion()` output
- `<scheme>/tpr_fdr.tsv`, `<scheme>/per_iteration.tsv` — summaries (`cFDR0` /
  `cTPR0` are what Figs 1/2 plot)
- `published/check_vs_paper.tsv` — iteration-by-iteration comparison with
  `analysis/confusionMats/`. DESeq2/edgeR/limma should match exactly; ALDEx3
  only if the same package commit is used.
- `figures/` — `fdr_<scheme>.png`, `tpr_<scheme>.png`, `compare_fdr.png`,
  `compare_tpr.png` (schemes overlaid, including the paper's own numbers),
  `fdr_per_iteration.png`

## Compute (one core per job, measured)

Per call, `.t` or `.r`: limma/edgeR 2–5 s, DESeq2 30–50 s (≤ 2 GB) on the
20k-gene datasets. With the defaults (2 schemes, 9 datasets, 100 iterations,
`.r` + `.t`): DESeq2 ≈ 30 core-hours, edgeR/limma a few core-hours. ALDEx3 is
not yet measured (`resources.aldex3` is a guess). ALDEx2 would add ≈ 2,000
core-hours per scheme.

## Caveats

- Package versions can only be approximated (the paper used Bioconductor-devel
  DESeq2 1.49.1 / edgeR 4.7.2 and an unrecorded ALDEx3 commit).
- `get_confusion()` counts a thinned gene whose |simulated LFC| is below the
  threshold as a false positive when detected, so FDR rises with the
  threshold by construction.
- For the `paired_model` scheme with ALDEx2, `paired.test = TRUE` is used;
  ALDEx3 uses fixed pair effects (`~ pair + group`, `lm`), not `lme4`.
