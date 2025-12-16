#!/usr/bin/python3

configfile: "config/fm_config.yaml"

from pathlib import Path

from workflow.lib.sample_utils import read_samples

# Read samples
samples = read_samples(config["paths"]["samples_table"])

# Paths and sample list
SAMPLES = list(samples.keys())
RESULTS = config["paths"]["results_dir"].rstrip("/")
LOGS = config["paths"]["logs_dir"].rstrip("/")
SCRATCH = config["paths"]["scratch_dir"].rstrip("/")
if config["infer_strandedness"]["enabled"]:
    SAMPLES_TO_INFER = [
        s
        for s in SAMPLES
        if samples[s]["strandedness"] in {"", "unknown", "unstranded", "na", "none"}
    ]
else:
    SAMPLES_TO_INFER = []
STRAND_REPORTS = (
    expand(f"{RESULTS}/qc/strandedness/{{sample}}.infer_experiment.txt", sample=SAMPLES_TO_INFER)
    if SAMPLES_TO_INFER
    else []
)
STRAND_SUMMARY = (
    [f"{RESULTS}/qc/strandedness/strandedness_summary.tsv"] if SAMPLES_TO_INFER else []
)
INFER_SUMMARY_PATH = STRAND_SUMMARY[0] if STRAND_SUMMARY else None
_INFER_CACHE = None


def resolved_strand(sample: str) -> str:
    """
    Return the strandedness call for a sample, preferring the sample sheet value.
    For unknown entries, fall back to the inferred summary if available; otherwise unstranded.
    """
    val = samples[sample]["strandedness"].strip().lower()
    if val not in {"", "unknown", "unstranded", "na", "none"}:
        return val
    if not INFER_SUMMARY_PATH:
        return "unstranded"
    global _INFER_CACHE
    if _INFER_CACHE is None:
        _INFER_CACHE = {}
        summary_path = Path(INFER_SUMMARY_PATH)
        if summary_path.exists():
            with summary_path.open() as handle:
                next(handle, None)  # skip header
                for line in handle:
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        _INFER_CACHE[parts[0]] = parts[1].lower()
    return _INFER_CACHE.get(sample, "unstranded")


def hisat2_strand_flag(sample: str) -> str:
    strand = resolved_strand(sample).upper()
    if strand in {"FR", "RF"}:
        return f"--rna-strandness {strand}"
    return ""


rule all:
    input:
        # Linked raw FASTQs (symlinks named by sample)
        expand(f"{RESULTS}/raw/{{sample}}_R1.fastq.gz", sample=SAMPLES),
        expand(f"{RESULTS}/raw/{{sample}}_R2.fastq.gz", sample=SAMPLES),
        # FastQC raw (renamed to sample IDs)
        expand(f"{RESULTS}/qc/raw/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        expand(f"{RESULTS}/qc/raw/{{sample}}_R2_fastqc.zip", sample=SAMPLES),
        # fastp outputs + reports
        expand(f"{RESULTS}/fastp/{{sample}}_R1.fastq.gz", sample=SAMPLES),
        expand(f"{RESULTS}/fastp/{{sample}}_R2.fastq.gz", sample=SAMPLES),
        expand(f"{RESULTS}/qc/fastp/{{sample}}.html", sample=SAMPLES),
        expand(f"{RESULTS}/qc/fastp/{{sample}}.json", sample=SAMPLES),
        # FastQC trimmed
        expand(f"{RESULTS}/qc/fastp/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        expand(f"{RESULTS}/qc/fastp/{{sample}}_R2_fastqc.zip", sample=SAMPLES),
        # Alignment
        expand(f"{RESULTS}/align/{{sample}}.bam", sample=SAMPLES),
        # Counts and summaries
        expand(f"{RESULTS}/counts/raw/{{sample}}.featurecounts.txt", sample=SAMPLES),
        expand(f"{RESULTS}/counts/raw/{{sample}}.featurecounts.txt.summary", sample=SAMPLES),
        f"{RESULTS}/counts/featurecounts_matrix.tsv",
        f"{RESULTS}/counts/featurecounts_summary.tsv",
        # MultiQC summary
        f"{RESULTS}/qc/multiqc/multiqc_report.html",
        # Strandedness inference (only if needed)
        STRAND_REPORTS,
        STRAND_SUMMARY,


################# QC #################
rule link_raw_fastq:
    input:
        fq1=lambda wc: samples[wc.sample]["fastq1"],
        fq2=lambda wc: samples[wc.sample]["fastq2"],
    output:
        r1=f"{RESULTS}/raw/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/raw/{{sample}}_R2.fastq.gz",
    threads: 1
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})"
        ln -sf {input.fq1} {output.r1}
        ln -sf {input.fq2} {output.r2}
        """


rule fastqc_raw:
    input:
        fq1=f"{RESULTS}/raw/{{sample}}_R1.fastq.gz",
        fq2=f"{RESULTS}/raw/{{sample}}_R2.fastq.gz",
    output:
        r1=f"{RESULTS}/qc/raw/{{sample}}_R1_fastqc.zip",
        r2=f"{RESULTS}/qc/raw/{{sample}}_R2_fastqc.zip",
    threads: config["qc"]["threads"]
    log: f"{LOGS}/fastqc/raw/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})" "$(dirname {output.r2})"
        fastqc -t {threads} {input.fq1} {input.fq2} -o "$(dirname {output.r1})" > {log} 2>&1
        """


rule fastp:
    input:
        fq1=f"{RESULTS}/raw/{{sample}}_R1.fastq.gz",
        fq2=f"{RESULTS}/raw/{{sample}}_R2.fastq.gz",
    output:
        r1=f"{RESULTS}/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/fastp/{{sample}}_R2.fastq.gz",
        html=f"{RESULTS}/qc/fastp/{{sample}}.html",
        json=f"{RESULTS}/qc/fastp/{{sample}}.json",
    threads: config["qc"]["threads"]
    params:
        extra=config["qc"].get("extra_params", "")
    log: f"{LOGS}/fastp/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})" "$(dirname {output.html})"
        fastp -i {input.fq1} -I {input.fq2} -o {output.r1} -O {output.r2} \
              -h {output.html} -j {output.json} {params.extra} > {log} 2>&1
        """


rule fastqc_trimmed:
    input:
        r1=f"{RESULTS}/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/fastp/{{sample}}_R2.fastq.gz",
    output:
        r1=f"{RESULTS}/qc/fastp/{{sample}}_R1_fastqc.zip",
        r2=f"{RESULTS}/qc/fastp/{{sample}}_R2_fastqc.zip",
    threads: config["qc"]["threads"]
    log: f"{LOGS}/fastqc/fastp/{{sample}}.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.r1})" "$(dirname {output.r2})"
        fastqc -t {threads} {input.r1} {input.r2} -o "$(dirname {output.r1})" > {log} 2>&1
        """


rule multiqc:
    input:
        raw_r1=expand(f"{RESULTS}/qc/raw/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        raw_r2=expand(f"{RESULTS}/qc/raw/{{sample}}_R2_fastqc.zip", sample=SAMPLES),
        trimmed_r1=expand(f"{RESULTS}/qc/fastp/{{sample}}_R1_fastqc.zip", sample=SAMPLES),
        trimmed_r2=expand(f"{RESULTS}/qc/fastp/{{sample}}_R2_fastqc.zip", sample=SAMPLES),
        fastp_html=expand(f"{RESULTS}/qc/fastp/{{sample}}.html", sample=SAMPLES),
        fastp_json=expand(f"{RESULTS}/qc/fastp/{{sample}}.json", sample=SAMPLES),
    output:
        html=f"{RESULTS}/qc/multiqc/multiqc_report.html",
    threads: 1
    log: f"{LOGS}/multiqc/multiqc.log"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.html})"
        multiqc {RESULTS} {LOGS} -o {RESULTS}/qc/multiqc -n multiqc_report.html > {log} 2>&1
        """


################# Strandedness inference #################
rule subsample_for_strand:
    input:
        fq1=lambda wc: samples[wc.sample]["fastq1"],
        fq2=lambda wc: samples[wc.sample]["fastq2"],
    output:
        r1=f"{SCRATCH}/strandedness/{{sample}}_R1.subsample.fastq.gz",
        r2=f"{SCRATCH}/strandedness/{{sample}}_R2.subsample.fastq.gz",
    threads: 1
    log: f"{LOGS}/strandedness/subsample/{{sample}}.log"
    params:
        subsample_target=config["infer_strandedness"]["subsample_reads"],
    shell:
        r"""
        python - <<'PY' > {log} 2>&1
from workflow.lib.sample_utils import determine_subsample_count, subsample_reads
target = "{params.subsample_target}"
count = determine_subsample_count("{input.fq1}", target)
subsample_reads("{input.fq1}", "{output.r1}", count)
subsample_reads("{input.fq2}", "{output.r2}", count)
print(f"Subsampled {count} reads to {output.r1} and {output.r2}")
PY
        """


rule infer_strand_align:
    input:
        r1=f"{SCRATCH}/strandedness/{{sample}}_R1.subsample.fastq.gz",
        r2=f"{SCRATCH}/strandedness/{{sample}}_R2.subsample.fastq.gz",
    output:
        bam=f"{SCRATCH}/strandedness/{{sample}}.subsample.bam",
    threads: min(4, int(config["hisat2"]["align_threads"]))
    log: f"{LOGS}/strandedness/hisat2/{{sample}}.log"
    params:
        index=config["references"]["hisat2_index_prefix"],
        extra=config["hisat2"]["extra_params"],
    shell:
        "hisat2 -x {params.index} -1 {input.r1} -2 {input.r2} -p {threads} {params.extra} "
        "| samtools sort -@ {threads} -o {output.bam} - 2> {log}"


rule infer_strand_report:
    input:
        bam=f"{SCRATCH}/strandedness/{{sample}}.subsample.bam",
    output:
        report=f"{RESULTS}/qc/strandedness/{{sample}}.infer_experiment.txt",
    threads: 1
    log: f"{LOGS}/strandedness/infer_experiment/{{sample}}.log"
    params:
        bed=config["references"]["annotation_bed"],
    shell:
        "infer_experiment.py -r {params.bed} -i {input.bam} > {output.report} 2> {log}"


rule infer_strand_summary:
    input:
        reports=STRAND_REPORTS,
    output:
        STRAND_SUMMARY
    threads: 1
    log: f"{LOGS}/strandedness/summary.log"
    shell:
        r"""
        python - <<'PY'
from pathlib import Path
from workflow.lib.sample_utils import infer_strand_from_rseqc

reports = {input.reports!r}
out_path = Path("{output[0]}")
out_path.parent.mkdir(parents=True, exist_ok=True)
with out_path.open("w") as handle:
    handle.write("sample\tcall\n")
    for report in reports:
        sample = Path(report).stem.split(".")[0]
        call = infer_strand_from_rseqc(report, min_fraction={config["infer_strandedness"]["min_fraction"]})
        handle.write(f"{sample}\t{call}\n")
PY
        """


################# Alignment #################
rule hisat2_align:
    input:
        r1=f"{RESULTS}/fastp/{{sample}}_R1.fastq.gz",
        r2=f"{RESULTS}/fastp/{{sample}}_R2.fastq.gz",
        strand_summary=INFER_SUMMARY_PATH if SAMPLES_TO_INFER else []
    output:
        bam=f"{RESULTS}/align/{{sample}}.bam",
    threads: config["hisat2"]["align_threads"]
    log: f"{LOGS}/hisat2/{{sample}}.log"
    params:
        index=config["references"]["hisat2_index_prefix"],
        extra=config["hisat2"]["extra_params"],
        strand=lambda wc: hisat2_strand_flag(wc.sample),
    shell:
        r"""
        set -euo pipefail
        mkdir -p "$(dirname {output.bam})" "$(dirname {log})"
        hisat2 -x {params.index} -1 {input.r1} -2 {input.r2} -p {threads} {params.strand} {params.extra} \
            2> {log} \
        | samtools sort -@ {threads} -o {output.bam} - 2>> {log}
        """


################# Counting #################
def strand_to_fc(sample: str) -> str:
    return {
        "fr": "1",
        "rf": "2",
        "unstranded": "0",
    }.get(resolved_strand(sample).lower(), "0")


def fc_paired_flag(sample: str) -> str:
    """
    Return '-p' for paired-end libraries; otherwise ''.
    """
    return "-p" if samples[sample]["library"].upper() == "PE" else ""


rule featurecounts:
    input:
        bam=f"{RESULTS}/align/{{sample}}.bam",
    output:
        counts=f"{RESULTS}/counts/raw/{{sample}}.featurecounts.txt",
        summary=f"{RESULTS}/counts/raw/{{sample}}.featurecounts.txt.summary",
    threads: config["featurecounts"]["threads"]
    log: f"{LOGS}/featurecounts/{{sample}}.log"
    params:
        annotation=config["references"]["annotation_gtf"],
        gtf_feature=config["featurecounts"]["annotation_type"],
        gtf_attr=config["featurecounts"]["gene_id_attribute"],
        extra=config["featurecounts"]["extra_params"],
        strand=lambda wc: strand_to_fc(wc.sample),
        paired=lambda wc: fc_paired_flag(wc.sample),
    shell:
        "featureCounts -T {threads} -a {params.annotation} -t {params.gtf_feature} -g {params.gtf_attr} "
        "-s {params.strand} {params.paired} {params.extra} -o {output.counts} {input.bam} > {log} 2>&1"


rule combine_featurecounts:
    input:
        counts=expand(f"{RESULTS}/counts/raw/{{sample}}.featurecounts.txt", sample=SAMPLES),
    output:
        matrix=f"{RESULTS}/counts/featurecounts_matrix.tsv",
    threads: 1
    params:
        samples=SAMPLES,
    script:
        "workflow/scripts/combine_featurecounts.py"


rule summarize_featurecounts:
    input:
        summaries=expand(f"{RESULTS}/counts/raw/{{sample}}.featurecounts.txt.summary", sample=SAMPLES),
    output:
        table=f"{RESULTS}/counts/featurecounts_summary.tsv",
    threads: 1
    params:
        samples=SAMPLES,
        biotype="all",
    script:
        "workflow/scripts/summarize_featurecounts.py"

rule featurecounts_trna:
    input:
        bam=f"{RESULTS}/align/{{sample}}.bam",
        trna_gtf=lambda wc: config["references"]["trna_gtf"],
    output:
        counts=f"{RESULTS}/counts/trna/{{sample}}.featurecounts.txt",
        summary=f"{RESULTS}/counts/trna/{{sample}}.featurecounts.txt.summary",
    threads: config["featurecounts"]["threads"]
    log: f"{LOGS}/featurecounts/trna/{{sample}}.log"
    params:
        gtf_feature="tRNA",
        gtf_attr=config["featurecounts"]["gene_id_attribute"],
        extra=config["featurecounts"]["extra_params"],
        strand=lambda wc: strand_to_fc(wc.sample),
        paired=lambda wc: fc_paired_flag(wc.sample),
    shell:
        "featureCounts -T {threads} -a {input.trna_gtf} -t {params.gtf_feature} -g {params.gtf_attr} "
        "-s {params.strand} {params.paired} {params.extra} -o {output.counts} {input.bam} > {log} 2>&1"


rule featurecounts_rrna:
    input:
        bam=f"{RESULTS}/align/{{sample}}.bam",
        rrna_gtf=lambda wc: config["references"]["rrna_gtf"],
    output:
        counts=f"{RESULTS}/counts/rrna/{{sample}}.featurecounts.txt",
        summary=f"{RESULTS}/counts/rrna/{{sample}}.featurecounts.txt.summary",
    threads: config["featurecounts"]["threads"]
    log: f"{LOGS}/featurecounts/rrna/{{sample}}.log"
    params:
        gtf_feature="rRNA",
        gtf_attr=config["featurecounts"]["gene_id_attribute"],
        extra=config["featurecounts"]["extra_params"],
        strand=lambda wc: strand_to_fc(wc.sample),
        paired=lambda wc: fc_paired_flag(wc.sample),
    shell:
        "featureCounts -T {threads} -a {input.rrna_gtf} -t {params.gtf_feature} -g {params.gtf_attr} "
        "-s {params.strand} {params.paired} {params.extra} -o {output.counts} {input.bam} > {log} 2>&1"


rule summarize_trna_rrna:
    input:
        cds=f"{RESULTS}/counts/featurecounts_summary.tsv",
        trna=expand(f"{RESULTS}/counts/trna/{{sample}}.featurecounts.txt.summary", sample=SAMPLES),
        rrna=expand(f"{RESULTS}/counts/rrna/{{sample}}.featurecounts.txt.summary", sample=SAMPLES),
    output:
        f"{RESULTS}/counts/trna_rrna_summary.tsv"
    threads: 1
    params:
        samples=SAMPLES,
    script:
        "workflow/scripts/summarize_trna_rrna.py"
