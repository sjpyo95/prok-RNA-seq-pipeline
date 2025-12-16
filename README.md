# prok-rnaseq-pipeline

Snakemake workflow for archaeal/prokaryotic RNA-seq: raw QC → trimming → alignment → counting (CDS/tRNA/rRNA) → DE → downstream plots. Designed to be config-driven and reproducible (single mamba/conda env).

## Requirements
- mamba/conda
- Snakemake (used via the env below)
- Reference files: genome FASTA, annotation GTF (CDS), optional tRNA/rRNA GTFs, hisat2 index (or enable index build)

## Install env
```bash
mamba env create -f envs/workflow.yaml        # first time
# or update
mamba env update -n prok-rnaseq-workflow -f envs/workflow.yaml
conda activate prok-rnaseq-workflow
```

## Configure
Edit a config YAML (e.g., `config/fm_config.yaml`):
- `paths`: samples table, results/logs/scratch roots
- `references`: fasta, gtf, tRNA/rRNA gtf, hisat2 index (set `build_*` booleans if you want the workflow to build them)
- `qc.extra_params`: extra fastp flags (passed through)
- Threads for qc/hisat2/featureCounts/samtools

Samples table (`sample_id	condition	fastq1	fastq2	library	strandedness	batch`), tab-delimited.

## Run
Dry-run first:
```bash
snakemake -n -p --cores 8 --configfile config/fm_config.yaml
```

Execute (adjust cores):
```bash
snakemake --cores 32 --configfile config/fm_config.yaml
```

Main rules/outputs (under `results/<project>/`):
- `fastp/` trimmed FASTQs + `qc/fastp/*.html/.json`
- `qc/raw` and `qc/fastp` FastQC, aggregated `qc/multiqc/multiqc_report.html` (includes hisat2 logs)
- `align/*.bam` (hisat2 aligned, sorted)
- `counts/raw/*.featurecounts.txt` (+ summary), combined `counts/featurecounts_matrix.tsv` and `counts/featurecounts_summary.tsv`

Optional strandedness inference (if enabled in config): subsample → RSeQC, summary in `qc/strandedness/`.

## Downstream DE and plots (R scripts)
Run these inside an R-capable env (same env works; DE packages are listed in `envs/workflow.yaml`):

- `workflow/scripts/deseq2_contrast.R`  
  Inputs: counts TSV, sample metadata TSV. Args: `-n` numerator condition, `-d` denominator, `-o` outdir. Outputs: DE results, significant table, normalized counts, MA plot; with ≥3 conditions also PCA and raw-vs-norm boxplot; heatmaps for all genes and DEGs.

- `workflow/scripts/deseq2_threeway.R`  
  Convenience wrapper for three contrasts (MMF_vs_MMM, MMFM_vs_MMM, MMFM_vs_MMF). Produces per-contrast DE results.

- `workflow/scripts/integrate_deseq2_results.R`  
  Merge `*_vs_*` contrast folders’ `deseq2_results.tsv` into `deseq2_integrated.tsv` (log2FC/padj/deg flags per contrast).

- `workflow/scripts/cluster_multi_heatmap_bar.R`  
  Combined heatmap + barplots across contrasts. Input: integrated DE table or single annotated DESeq2 file; optional normalized counts + metadata. Supports per-sample heatmap or `--condition_means` (average by condition). Cluster labels on the right of the barplots.

- `workflow/scripts/multi_anno_plot.R`  
  arCOG/Cluster bubble charts sized by %DEGs (|log2FC|>=1 & padj<alpha) per contrast.

- `workflow/scripts/counts_biotype_stats.py`  
  Summarize assigned reads to CDS/tRNA/rRNA via featureCounts; outputs `count_stat.tsv` and stacked % plot.

- `workflow/scripts/gene_group_heatmap.R` / `gene_heatmap_all.R`  
  Heatmaps for selected genes or all genes/DEGs, using normalized counts + metadata.

Example DE + integrate:
```bash
Rscript workflow/scripts/deseq2_contrast.R \
  -c results/fm/de/rawCounts_filt.tsv \
  -s results/fm/de/sample_metadata.tsv \
  -n MMF -d MMM -o results/fm/de/MMF_vs_MMM --alpha 0.05

Rscript workflow/scripts/integrate_deseq2_results.R \
  -i results/fm/de \
  -o results/fm/de/integrate
```

Example multi-contrast heatmap+bars:
```bash
Rscript workflow/scripts/cluster_multi_heatmap_bar.R \
  -i results/fm/de/integrate/deseq2_integrated.tsv \
  -m results/fm/de/normalized_counts.tsv \
  --metadata results/fm/de/sample_metadata.tsv \
  --condition_means \
  -o results/fm/de/plots
```

## Tips
- Keep paths in config absolute or relative to repo root.
- Update `qc.extra_params` to tune fastp (e.g., poly-G trimming).
- If you want per-condition heatmaps, provide metadata (`sample_id`, `condition`) and add `--condition_means`.
- MultiQC scans both results and logs; hisat2 alignment stats appear there.***
