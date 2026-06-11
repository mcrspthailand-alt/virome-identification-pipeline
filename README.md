# Virome Identification Pipeline

Reusable, configurable pipeline template for identifying putative viral contigs from metagenomic assemblies and paired-end reads.

This repository intentionally contains the reusable template only:

```text
virome_pipeline_template/
```

Project-specific scripts, sample-specific configs, sequencing reads, assemblies, databases, and pipeline outputs should not be committed to this repository.

## Quick Start

Copy the template folder into a new project directory:

```bash
cp -r virome_pipeline_template my_project_virome_pipeline
cd my_project_virome_pipeline
```

Edit:

```text
config.env
samples.tsv
```

Run:

```bash
bash run_virome_pipeline.sh
```

or:

```bash
bash run_virome_pipeline.sh /path/to/config.env
```

## Required Inputs

`samples.tsv` must contain one row per sample:

```text
sample	assembly	r1	r2
sample01	/path/to/sample01.contigs.fa.gz	/path/to/sample01_R1.fastq.gz	/path/to/sample01_R2.fastq.gz
```

Required columns:

- `sample`: sample ID
- `assembly`: assembled contig FASTA
- `r1`: paired-end read 1 FASTQ
- `r2`: paired-end read 2 FASTQ

Assemblies can be uncompressed or gzipped FASTA:

```text
.fa
.fasta
.fa.gz
.fasta.gz
```

## Configuration

Set paths and parameters in `config.env`:

```bash
OUT_BASE="/path/to/project/viral_id_pipeline"
SAMPLES_TSV="./samples.tsv"

THREADS=32
MIN_CONTIG_LEN=500
DEPTH_CUTOFF=2
VS2_MIN_LEN=500

DIAMOND_DB="/path/to/databases/diamond/virus_nr.dmnd"
GENOMAD_DB="/path/to/databases/genomad/genomad_db"
VS2_DB_DIR="/path/to/databases/virsorter2_db"
CHECKV_DB="/path/to/databases/checkv_db/checkv-db-v1.5"
```

DIAMOND thresholds:

```bash
DIAMOND_EVALUE="1e-10"
DIAMOND_BITSCORE=100
DIAMOND_ALNLEN=60
DIAMOND_MAX_TARGET_SEQS=5
```

## Workflow

The pipeline runs:

1. Contig filtering and length binning
2. Read mapping to contigs
3. Read count and depth calculation
4. VirSorter2
5. geNomad
6. CheckV
7. Prodigal
8. DIAMOND `blastx`
9. DIAMOND `blastp`
10. Evidence-table merge
11. NCBI taxonomy assignment
12. Final viral-contig tables
13. Output README generation

## Why Both DIAMOND blastx And blastp?

The template uses both homology approaches:

- `blastx`: nucleotide contigs are translated in six frames and searched against a viral protein database. This is useful for fragmented or short environmental contigs.
- `blastp`: Prodigal-predicted proteins are searched against the same viral protein database. This provides ORF-level supporting evidence.

They are reported separately:

```text
blastx_hit
blastp_hit
```

## Viral Calling Rule

Putative viral contigs are called by evidence:

```text
putative_viral = vs2_hit OR genomad_hit OR blastx_hit OR blastp_hit
```

Confidence:

```text
High    = viral_evidence_count >= 2
Medium  = viral_evidence_count == 1
None    = viral_evidence_count == 0
```

Evidence sources counted:

```text
VirSorter2
geNomad
DIAMOND blastx
DIAMOND blastp
```

CheckV is included as annotation only and is not used to discard contigs.

## Output Structure

Outputs are written under `${OUT_BASE}`:

```text
viral_id_pipeline/
  logs/
  samples/
  merged/
  taxonomy/
  README.md
```

Per-sample outputs:

```text
samples/<sample>/
  00_input/
  01_contigs/
  02_mapping/
  03_depth_filtered_contigs/
  04_virsorter2/
  05_genomad/
  06_checkv/
  07_diamond_ge1000/
  07_diamond_len500_999/
  08_coverage/
  run_summary.tsv
```

Important merged outputs:

```text
merged/all_contig_read_counts.tsv
merged/diamond_blastx_best_hits.tsv
merged/diamond_blastp_best_hits.tsv
merged/all_viral_evidence.tsv
merged/all_viral_evidence.parsed_organisms.tsv
merged/virus_taxonomy_map.tsv
merged/virus_reads_matrix.tsv
merged/virus_reads_matrix_with_taxonomy.tsv
merged/virus_presence_absence_matrix.tsv
merged/virus_reads_by_sample_long.tsv
merged/all_contigs_with_annotations.tsv
merged/final_virus_contigs.tsv
merged/final_virus_contigs_bacteriophage.tsv
merged/final_virus_contigs_eukaryotic.tsv
merged/final_virus_contigs_archaeal.tsv
merged/final_virus_contigs_unidentified.tsv
merged/final_virus_summary.txt
```

The pipeline writes `${OUT_BASE}/README.md`, which explains output files and columns based on the actual generated TSV headers.

## Required Tools

Command-line tools expected in `PATH`:

```text
bash
python3
awk
gzip
seqkit
bbduk.sh
minimap2
samtools
mosdepth
virsorter
genomad
checkv
prodigal
diamond
```

Python package:

```text
pandas
```
