# Virome Identification Pipeline

Reusable metagenomic virome identification pipeline for environmental, animal, clinical, or other metagenomic projects.

This repository is designed to store the reusable pipeline template only. It should not contain sample-specific scripts, sequencing reads, assembly files, database files, or analysis outputs.

```text
virome_pipeline_template/
```

For each new project, copy the template folder, edit the configuration and sample sheet, then run the pipeline in that project directory.

## What This Pipeline Does

The pipeline identifies putative viral contigs from metagenomic assemblies by combining multiple evidence types:

- tool-based viral prediction from VirSorter2
- tool-based viral prediction from geNomad
- nucleotide-to-protein homology from DIAMOND blastx
- predicted-protein homology from DIAMOND blastp
- genome quality annotation from CheckV
- read support from mapping host-removed reads back to contigs
- taxonomy assignment from viral hit names and NCBI taxonomy files

The final viral call is evidence-based:

```text
putative_viral = vs2_hit OR genomad_hit OR blastx_hit OR blastp_hit
```

CheckV results are reported in the final table, but CheckV is not used to discard contigs. This is intentional because environmental viral contigs can be fragmented, novel, or low-completeness while still having useful viral evidence from VirSorter2, geNomad, or protein homology.

## Repository Layout

```text
virome-identification-pipeline/
  README.md
  .gitignore
  virome_pipeline_template/
    README.md
    config.env
    samples.tsv
    run_virome_pipeline.sh
    scripts/
      01_process_samples.sh
      02_build_evidence_table.py
      03_assign_taxonomy.py
      04_finalize_results.py
      05_write_output_readme.py
```

Only the template should be pushed to GitHub. Project-specific folders such as air, sediment, bat, or other sample-specific runs should remain local or be stored in a separate analysis repository if needed.

## Quick Start

Copy the template into a new project folder:

```bash
cp -r virome_pipeline_template my_project_virome_pipeline
cd my_project_virome_pipeline
```

Edit these two files:

```text
config.env
samples.tsv
```

Run the full pipeline:

```bash
bash run_virome_pipeline.sh
```

Or run with an explicit config file:

```bash
bash run_virome_pipeline.sh /path/to/config.env
```

## Required Input Files

The pipeline starts from two input types for each sample:

1. assembled contigs in FASTA format
2. paired-end host-removed reads in FASTQ format

Assemblies may be compressed or uncompressed:

```text
.fa
.fasta
.fa.gz
.fasta.gz
```

Reads should be paired FASTQ files:

```text
.fastq.gz
.fq.gz
```

The pipeline does not perform host removal or assembly. It assumes those steps were already completed upstream.

## Sample Sheet

Samples are defined in `samples.tsv`.

Required columns:

```text
sample	assembly	r1	r2
```

Example:

```text
sample	assembly	r1	r2
sample01	/path/to/MEGAHIT-sample01.contigs.fa.gz	/path/to/sample01.unmapped_1.fastq.gz	/path/to/sample01.unmapped_2.fastq.gz
sample02	/path/to/MEGAHIT-sample02.contigs.fa.gz	/path/to/sample02.unmapped_1.fastq.gz	/path/to/sample02.unmapped_2.fastq.gz
```

Column meaning:

| Column | Meaning |
|---|---|
| `sample` | Sample ID used in output folder names and merged tables |
| `assembly` | Path to assembled contig FASTA |
| `r1` | Path to paired-end read 1 FASTQ |
| `r2` | Path to paired-end read 2 FASTQ |

Use simple sample names without spaces. Hyphens and underscores are fine.

## Configuration

Main settings are stored in `config.env`.

Example:

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

Important parameters:

| Parameter | Meaning |
|---|---|
| `OUT_BASE` | Main output directory for the run |
| `SAMPLES_TSV` | Sample sheet path |
| `THREADS` | Number of CPU threads used by supported tools |
| `MIN_CONTIG_LEN` | Minimum contig length retained for analysis |
| `DEPTH_CUTOFF` | Depth threshold used for depth-filtered contig outputs |
| `VS2_MIN_LEN` | Minimum contig length passed to VirSorter2 |
| `DIAMOND_DB` | DIAMOND-formatted viral protein database |
| `GENOMAD_DB` | geNomad database directory |
| `VS2_DB_DIR` | VirSorter2 database directory |
| `CHECKV_DB` | CheckV database directory |

DIAMOND thresholds:

```bash
DIAMOND_EVALUE="1e-10"
DIAMOND_BITSCORE=100
DIAMOND_ALNLEN=60
DIAMOND_MAX_TARGET_SEQS=5
```

These thresholds are applied when selecting useful DIAMOND hits for downstream evidence tables.

## Workflow Overview

The pipeline runs from contig/read input to final merged viral tables:

```text
samples.tsv + config.env
        |
        v
per-sample contig preparation
        |
        v
read mapping, read counts, depth calculation
        |
        v
VirSorter2 + geNomad + CheckV
        |
        v
Prodigal gene prediction
        |
        v
DIAMOND blastx + DIAMOND blastp
        |
        v
merged viral evidence table
        |
        v
taxonomy assignment
        |
        v
final viral contig tables and abundance matrices
```

## Step-By-Step Explanation

### 1. Read Config And Sample Sheet

`run_virome_pipeline.sh` loads `config.env`, checks required variables, then reads `samples.tsv`.

Each row in `samples.tsv` is treated as one sample. The sample ID becomes the output folder name under:

```text
${OUT_BASE}/samples/<sample>/
```

### 2. Prepare Contigs

For each sample, the pipeline prepares assembly FASTA files for downstream tools.

This step:

- accepts gzipped or uncompressed FASTA
- standardizes input locations inside the per-sample output folder
- filters contigs using `MIN_CONTIG_LEN`
- separates contigs into length groups where needed

The template supports both short environmental contigs and longer assembled contigs. Contigs from 500 bp upward are kept by default.

### 3. Map Reads Back To Contigs

Host-removed paired-end reads are mapped back to assembled contigs.

This step uses:

- `minimap2` for read alignment
- `samtools` for BAM processing
- `mosdepth` for depth calculation

The mapping outputs provide read support for each contig. This is useful for:

- estimating which viral contigs are supported by reads
- building read-count matrices across samples
- comparing viral signal among samples
- separating detected viral evidence from unsupported contigs

### 4. Run VirSorter2

VirSorter2 predicts viral sequences from assembled contigs using viral hallmark and sequence composition signals.

The pipeline records whether each contig has VirSorter2 evidence:

```text
vs2_hit
```

VirSorter2 is one of the tool-based viral prediction sources used in the final evidence rule.

### 5. Run geNomad

geNomad predicts viruses and plasmids from metagenomic contigs.

The pipeline records whether each contig has geNomad viral evidence:

```text
genomad_hit
```

geNomad is another tool-based viral prediction source used in the final evidence rule.

### 6. Run CheckV

CheckV estimates viral genome quality and completeness where possible.

Important: CheckV is annotation only in this pipeline.

The pipeline does not remove low-quality CheckV contigs automatically. Instead, CheckV fields are retained in the final output table so users can filter later depending on the research question.

This is useful because:

- environmental viral contigs may be partial
- novel viruses may not receive high completeness estimates
- discarding low-completeness contigs too early can remove biologically interesting viral candidates

### 7. Predict Proteins With Prodigal

Prodigal predicts protein-coding genes from contigs.

The predicted proteins are used for DIAMOND blastp.

This gives ORF-level homology evidence and complements blastx.

### 8. Run DIAMOND blastx

DIAMOND blastx translates nucleotide contigs in all reading frames and searches them against a viral protein database.

This is useful when:

- contigs are short or fragmented
- gene prediction is uncertain
- a contig contains partial viral protein similarity
- environmental viruses are too divergent for nucleotide-level searches

The pipeline reports best blastx hits in:

```text
merged/diamond_blastx_best_hits.tsv
```

The final evidence table includes:

```text
blastx_hit
```

### 9. Run DIAMOND blastp

DIAMOND blastp searches Prodigal-predicted proteins against the same viral protein database.

This is useful when:

- predicted ORFs are reliable
- protein-level hits need to be interpreted per gene
- multiple genes on the same contig support viral origin

The pipeline reports best blastp hits in:

```text
merged/diamond_blastp_best_hits.tsv
```

The final evidence table includes:

```text
blastp_hit
```

### 10. Merge Viral Evidence

The pipeline merges evidence from:

- VirSorter2
- geNomad
- DIAMOND blastx
- DIAMOND blastp
- CheckV
- read mapping
- contig metadata

The main merged evidence table is:

```text
merged/all_viral_evidence.tsv
```

The central viral calling rule is:

```text
putative_viral = vs2_hit OR genomad_hit OR blastx_hit OR blastp_hit
```

The evidence count is:

```text
viral_evidence_count = number of positive evidence sources among:
  VirSorter2
  geNomad
  blastx
  blastp
```

Confidence labels:

```text
High    = viral_evidence_count >= 2
Medium  = viral_evidence_count == 1
None    = viral_evidence_count == 0
```

This means a contig can be called putatively viral if it is predicted by at least one viral prediction tool or supported by protein homology.

### 11. Assign Taxonomy

The taxonomy scripts parse viral organism names from DIAMOND hits and map them to NCBI taxonomic IDs and ranks.

Taxonomy outputs include:

```text
merged/all_viral_evidence.parsed_organisms.tsv
merged/virus_names_with_taxid.tsv
merged/taxid_ranks.tsv
merged/virus_taxonomy_map.tsv
```

The taxonomy assignment is best interpreted as hit-derived taxonomy, not definitive species identification. Environmental viral contigs often match conserved proteins shared across related viruses.

### 12. Build Final Tables

The final result table is:

```text
merged/final_virus_contigs.tsv
```

This table contains viral candidates with combined evidence, read support, homology information, taxonomy fields, and CheckV annotation.

The pipeline also separates final viral candidates into broad groups:

```text
merged/final_virus_contigs_bacteriophage.tsv
merged/final_virus_contigs_eukaryotic.tsv
merged/final_virus_contigs_archaeal.tsv
merged/final_virus_contigs_unidentified.tsv
```

These grouped files are intended for easier downstream review.

### 13. Write Output README

At the end of the run, the pipeline writes:

```text
${OUT_BASE}/README.md
```

That file describes the actual output files generated by the run and summarizes the columns found in the final TSV files.

## Output Structure

Main output directory:

```text
${OUT_BASE}/
  logs/
  samples/
  merged/
  taxonomy/
  README.md
```

### Per-Sample Output

Each sample gets its own folder:

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

Folder meaning:

| Folder | Meaning |
|---|---|
| `00_input/` | Standardized links or copies of input files |
| `01_contigs/` | Filtered and prepared contig FASTA files |
| `02_mapping/` | Read mapping BAM files and mapping summaries |
| `03_depth_filtered_contigs/` | Contigs passing depth-related filters |
| `04_virsorter2/` | VirSorter2 output |
| `05_genomad/` | geNomad output |
| `06_checkv/` | CheckV output |
| `07_diamond_ge1000/` | DIAMOND outputs for longer contigs |
| `07_diamond_len500_999/` | DIAMOND outputs for 500-999 bp contigs |
| `08_coverage/` | Coverage and read-support summaries |
| `run_summary.tsv` | Per-sample run summary |

### Merged Output

Important merged files:

| File | Meaning |
|---|---|
| `all_contig_read_counts.tsv` | Read counts and support information per contig |
| `diamond_blastx_best_hits.tsv` | Best DIAMOND blastx hits per contig |
| `diamond_blastp_best_hits.tsv` | Best DIAMOND blastp hits per contig or predicted protein |
| `all_viral_evidence.tsv` | Main evidence table before final taxonomy/final filtering |
| `all_viral_evidence.parsed_organisms.tsv` | Evidence table with parsed organism names from hits |
| `virus_names_with_taxid.tsv` | Parsed virus names mapped to NCBI taxids |
| `taxid_ranks.tsv` | NCBI rank information for mapped taxids |
| `virus_taxonomy_map.tsv` | Consolidated taxonomy map used for final tables |
| `virus_reads_matrix.tsv` | Virus read-count matrix across samples |
| `virus_reads_matrix_with_taxonomy.tsv` | Read-count matrix with taxonomy annotation |
| `virus_presence_absence_matrix.tsv` | Presence/absence matrix across samples |
| `virus_reads_by_sample_long.tsv` | Long-format read table for plotting or statistics |
| `all_contigs_with_annotations.tsv` | All contigs with merged annotations |
| `final_virus_contigs.tsv` | Final putative viral contig table |
| `final_virus_contigs_bacteriophage.tsv` | Final candidates classified as bacteriophage-like |
| `final_virus_contigs_eukaryotic.tsv` | Final candidates classified as eukaryotic-virus-like |
| `final_virus_contigs_archaeal.tsv` | Final candidates classified as archaeal-virus-like |
| `final_virus_contigs_unidentified.tsv` | Final candidates without clear broad host group |
| `final_virus_summary.txt` | Human-readable summary of final counts |

## How To Interpret The Final Table

The most important file for manual review is:

```text
merged/final_virus_contigs.tsv
```

Key column groups usually include:

| Column group | Meaning |
|---|---|
| sample and contig columns | Sample ID, contig ID, contig length, and contig metadata |
| read support columns | Mapped read counts, depth, and coverage-related values |
| VirSorter2 columns | Whether VirSorter2 predicted viral origin and related output fields |
| geNomad columns | Whether geNomad predicted viral origin and related output fields |
| CheckV columns | Completeness, quality, contamination, and related CheckV annotations |
| DIAMOND blastx columns | Best translated nucleotide-to-protein hit information |
| DIAMOND blastp columns | Best predicted-protein hit information |
| taxonomy columns | Parsed organism, taxid, lineage, and broad virus group |
| evidence columns | `vs2_hit`, `genomad_hit`, `blastx_hit`, `blastp_hit`, `viral_evidence_count`, `putative_viral`, and confidence |

Recommended interpretation:

- use `putative_viral` to identify candidates called by the pipeline
- use `viral_evidence_count` and confidence to prioritize candidates
- use CheckV quality and completeness as supporting quality information
- use read counts and depth to assess sample support
- manually inspect important contigs, especially when only one evidence source is present

## Why Both blastx And blastp Are Included

Using both DIAMOND blastx and blastp is not strange; they answer related but different questions.

`blastx` asks:

```text
Does this nucleotide contig encode any protein-like region similar to known viral proteins?
```

`blastp` asks:

```text
Do the proteins predicted by Prodigal match known viral proteins?
```

In environmental metagenomics, this is useful because:

- blastx can detect partial coding regions even when gene prediction is imperfect
- blastp gives cleaner ORF-level evidence when Prodigal predicts good proteins
- agreement between blastx and blastp increases confidence
- disagreement can still be informative for fragmented or novel viral contigs

The two results are kept as separate evidence columns so users can judge them independently.

## Required Tools

The following command-line tools should be available in `PATH`:

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

Databases required:

| Database | Used by |
|---|---|
| VirSorter2 database | VirSorter2 |
| geNomad database | geNomad |
| CheckV database | CheckV |
| DIAMOND viral protein database | blastx and blastp |
| NCBI taxonomy dump | taxonomy assignment |

## Common Problems

### Missing CheckV Database

Example error:

```text
[ERROR] missing CHECKV_DB: /path/to/checkv-db-v1.5
```

Fix `CHECKV_DB` in `config.env`:

```bash
CHECKV_DB="/correct/path/to/checkv_db/checkv-db-v1.5"
```

### Command Not Found

Example:

```text
diamond: command not found
```

Load the correct conda environment, module, or software path before running the pipeline.

Example:

```bash
conda activate virome
bash run_virome_pipeline.sh
```

### Wrong Sample Paths

If the pipeline reports missing assembly or FASTQ files, check `samples.tsv`.

Each path should exist on the machine where the pipeline is run:

```bash
ls -lh /path/to/file
```

### Accidentally Adding Large Files To Git

Do not commit:

- FASTQ files
- FASTA assemblies
- BAM files
- database folders
- pipeline output folders
- project-specific copied pipeline folders

The GitHub repository should contain the reusable template, not raw data or finished analysis outputs.

## Updating This GitHub Repository From PowerShell

Run these commands from Windows PowerShell.

Go to the local repository:

```powershell
cd "C:\Users\HP\Documents\KKU\งาน sequencing\eDNA\air_script"
```

Check current remote:

```powershell
git remote -v
```

Expected remote:

```text
https://github.com/mcrspthailand-alt/virome-identification-pipeline.git
```

Check branch and changed files:

```powershell
git status --short --branch
```

Update your local copy from GitHub before editing further:

```powershell
git pull
```

Add only reusable pipeline files:

```powershell
git add README.md .gitignore virome_pipeline_template
```

Check what will be committed:

```powershell
git status --short
```

Commit:

```powershell
git commit -m "Update pipeline documentation"
```

Push to GitHub:

```powershell
git push
```

After pushing, check status again:

```powershell
git status --short --branch
```

If everything is clean, it should show that `main` is up to date with `origin/main`.

## Recommended Git Rule For This Repository

Use this repository for:

- reusable pipeline template
- documentation
- example config structure
- example sample sheet structure
- helper scripts needed by the reusable pipeline

Do not use this repository for:

- one-off air, sediment, bat, or other sample-specific scripts
- raw sequencing files
- assemblies
- output folders
- database folders
- temporary HPC logs

This keeps the GitHub repository clean and makes it useful for other people who want to adapt the pipeline to their own projects.
