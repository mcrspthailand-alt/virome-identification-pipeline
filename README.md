# Virome Identification Pipeline Scripts

## Quick Context For New Chat

This project contains scripts for identifying putative viral contigs from metagenomic assemblies and paired-end reads. The work started from sediment/air-specific scripts, but the current reusable version is the generic template in:

```text
virome_pipeline_template/
```

For new projects, use `virome_pipeline_template/` and edit only:

```text
virome_pipeline_template/config.env
virome_pipeline_template/samples.tsv
```

The older air/sediment scripts remain in the repository as project-specific history/reference. They should not be the first choice for new analyses.

Current local git state at the time this README was updated:

```text
local commit exists: d7388cb Add configurable virome identification pipeline
GitHub remote: not configured yet
```

To push later, add a GitHub remote and push:

```bash
git remote add origin https://github.com/<user>/<repo>.git
git branch -M main
git push -u origin main
```

## What The Template Does

The template runs a configurable virome-identification workflow:

1. Reads sample metadata from `samples.tsv`
2. Processes each sample assembly and paired FASTQ files
3. Maps reads back to contigs
4. Calculates read counts and depth
5. Runs viral prediction tools
6. Runs homology searches
7. Merges all evidence
8. Assigns taxonomy from NCBI taxdump
9. Writes final viral-contig tables
10. Generates an output README explaining result files and columns

## Main Command

```bash
cd virome_pipeline_template
bash run_virome_pipeline.sh
```

Or:

```bash
bash run_virome_pipeline.sh /path/to/config.env
```

## Required Input Files

`samples.tsv` must contain:

```text
sample	assembly	r1	r2
sample01	/path/to/sample01.contigs.fa.gz	/path/to/sample01_R1.fastq.gz	/path/to/sample01_R2.fastq.gz
```

Required columns:

- `sample`: sample ID
- `assembly`: contig FASTA from assembler such as MEGAHIT
- `r1`: paired-end read 1 FASTQ
- `r2`: paired-end read 2 FASTQ

Assemblies can be:

```text
.fa
.fasta
.fa.gz
.fasta.gz
```

## Config File

Edit `config.env`:

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
CHECKV_DB="/path/to/databases/checkv/checkv-db-v1.5"
```

Homology thresholds:

```bash
DIAMOND_EVALUE="1e-10"
DIAMOND_BITSCORE=100
DIAMOND_ALNLEN=60
DIAMOND_MAX_TARGET_SEQS=5
```

## Pipeline Stages

### 1. `scripts/01_process_samples.sh`

Per-sample processing:

- decompress/copy assembly FASTA
- filter low-complexity contigs with `bbduk.sh`
- split contigs into length bins:
  - `len500_999`
  - `ge1000`
- map reads to contigs with `minimap2`
- sort/index BAM with `samtools`
- calculate depth with `mosdepth`
- run `VirSorter2`
- run `geNomad`
- run `CheckV`
- run `Prodigal`
- run DIAMOND `blastx`
- run DIAMOND `blastp`

### 2. `scripts/02_build_evidence_table.py`

Builds the main evidence table by merging:

- read counts
- DIAMOND `blastx`
- DIAMOND `blastp`
- VirSorter2
- geNomad
- CheckV
- Prodigal

It also creates:

```text
putative_viral
viral_evidence
viral_evidence_count
viral_confidence
```

### 3. `scripts/03_assign_taxonomy.py`

Taxonomy stage:

- parses virus/organism names from DIAMOND hit titles
- maps names to NCBI taxid using `nodes.dmp` and `names.dmp`
- downloads NCBI taxdump if needed
- creates taxonomy map and read matrices

### 4. `scripts/04_finalize_results.py`

Final result stage:

- merges taxonomy into the evidence table
- classifies broad virus groups
- writes final viral-contig files

### 5. `scripts/05_write_output_readme.py`

Writes:

```text
${OUT_BASE}/README.md
```

This generated README describes every output file and lists columns found in the real output TSV headers.

## Tools Used

Required command-line tools:

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

## Why Both DIAMOND blastx And blastp?

The template intentionally runs both:

```text
DIAMOND blastx
DIAMOND blastp
```

`blastx` uses nucleotide contigs as query and translates them in six frames. It is more sensitive for fragmented or short environmental contigs.

`blastp` uses Prodigal-predicted proteins as query. It is useful as ORF-level supporting evidence and is often cleaner/easier to interpret.

They are treated as separate evidence sources:

```text
blastx_hit
blastp_hit
```

## Viral Calling Criteria

Main evidence rule:

```text
putative_viral = vs2_hit OR genomad_hit OR blastx_hit OR blastp_hit
```

Confidence rule:

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

CheckV is included as annotation only. It is not used to discard contigs.

## Output Structure

Template outputs:

```text
${OUT_BASE}/
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

Merged cohort-level outputs:

```text
merged/
  all_contig_read_counts.tsv
  diamond_blastx_best_hits.tsv
  diamond_blastp_best_hits.tsv
  all_viral_evidence.tsv
  all_viral_evidence.parsed_organisms.tsv
  virus_names_with_taxid.tsv
  taxid_ranks.tsv
  virus_taxonomy_map.tsv
  virus_reads_matrix.tsv
  virus_reads_matrix_with_taxonomy.tsv
  virus_presence_absence_matrix.tsv
  virus_reads_by_sample_long.tsv
  all_contigs_with_annotations.tsv
  final_virus_contigs.tsv
  final_virus_contigs_bacteriophage.tsv
  final_virus_contigs_eukaryotic.tsv
  final_virus_contigs_archaeal.tsv
  final_virus_contigs_unidentified.tsv
  final_virus_summary.txt
```

## Most Important Output Files

Use this for checking all evidence before taxonomy filtering/finalization:

```text
merged/all_viral_evidence.tsv
```

Use this as the main final viral-contig table:

```text
merged/final_virus_contigs.tsv
```

Use these for broad group subsets:

```text
merged/final_virus_contigs_bacteriophage.tsv
merged/final_virus_contigs_eukaryotic.tsv
merged/final_virus_contigs_archaeal.tsv
merged/final_virus_contigs_unidentified.tsv
```

Use these for sample-by-virus abundance/presence summaries:

```text
merged/virus_reads_matrix.tsv
merged/virus_reads_matrix_with_taxonomy.tsv
merged/virus_presence_absence_matrix.tsv
merged/virus_reads_by_sample_long.tsv
```

## Notes For Future ChatGPT/Codex Sessions

If starting a new chat, the fastest context is:

```text
We have a configurable virome identification pipeline template in virome_pipeline_template/.
It uses samples.tsv + config.env, runs VirSorter2/geNomad/CheckV/Prodigal/DIAMOND blastx/DIAMOND blastp, then builds evidence and taxonomy tables.
putative_viral = vs2_hit OR genomad_hit OR blastx_hit OR blastp_hit.
viral_confidence is High for >=2 evidence sources, Medium for 1.
CheckV is annotation only.
Main output is merged/final_virus_contigs.tsv.
README.md and virome_pipeline_template/README.md describe the project.
```

## Project-Specific Scripts

Root-level air/sediment scripts are project-specific historical scripts. They are kept for reference, but for new projects use:

```text
virome_pipeline_template/
```
