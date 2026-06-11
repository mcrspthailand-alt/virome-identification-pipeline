# Virome Identification Scripts

This repository contains scripts for viral-contig identification from metagenomic assemblies and paired-end reads. It includes a reusable pipeline template plus project-specific scripts that were developed for air/environmental samples.

## Recommended Template

Use the generic template in:

```text
virome_pipeline_template/
```

This template is intended for reuse across projects. It is configured through:

```text
virome_pipeline_template/config.env
virome_pipeline_template/samples.tsv
```

Run the full workflow with:

```bash
cd virome_pipeline_template
bash run_virome_pipeline.sh
```

or with an explicit config file:

```bash
bash run_virome_pipeline.sh /path/to/config.env
```

## Template Inputs

`samples.tsv` must contain one row per sample:

```text
sample	assembly	r1	r2
sample01	/path/to/sample01.contigs.fa.gz	/path/to/sample01_R1.fastq.gz	/path/to/sample01_R2.fastq.gz
```

Assemblies may be `.fa`, `.fasta`, `.fa.gz`, or `.fasta.gz`.

Edit `config.env` to set:

```bash
OUT_BASE="/path/to/project/viral_id_pipeline"
DIAMOND_DB="/path/to/databases/diamond/virus_nr.dmnd"
GENOMAD_DB="/path/to/databases/genomad/genomad_db"
VS2_DB_DIR="/path/to/databases/virsorter2_db"
CHECKV_DB="/path/to/databases/checkv/checkv-db-v1.5"
```

## Workflow

The template runs these stages:

1. `01_process_samples.sh`
   - filter contigs
   - map reads to contigs with `minimap2`
   - calculate read counts and depth with `samtools` and `mosdepth`
   - run `VirSorter2`
   - run `geNomad`
   - run `CheckV`
   - run `Prodigal`
   - run DIAMOND `blastx`
   - run DIAMOND `blastp`

2. `02_build_evidence_table.py`
   - merge read counts, BLASTx, BLASTp, VirSorter2, geNomad, CheckV, and Prodigal outputs
   - calculate viral evidence columns

3. `03_assign_taxonomy.py`
   - parse virus names from DIAMOND hit titles
   - map names to NCBI taxonomy
   - generate taxonomy and read-count matrices

4. `04_finalize_results.py`
   - merge taxonomy into the evidence table
   - classify broad virus groups
   - write final viral-contig tables

5. `05_write_output_readme.py`
   - generate an output README describing result files and columns

## Viral Evidence Rule

The template uses an evidence-based call:

```text
putative_viral = vs2_hit OR genomad_hit OR blastx_hit OR blastp_hit
```

Confidence is assigned as:

```text
High    = viral_evidence_count >= 2
Medium  = viral_evidence_count == 1
None    = viral_evidence_count == 0
```

CheckV results are reported as annotation columns but are not used to discard contigs.

## Main Outputs

Project-level results are written to:

```text
${OUT_BASE}/merged/
```

Important files:

```text
all_contig_read_counts.tsv
diamond_blastx_best_hits.tsv
diamond_blastp_best_hits.tsv
all_viral_evidence.tsv
all_viral_evidence.parsed_organisms.tsv
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

The pipeline also writes:

```text
${OUT_BASE}/README.md
```

This generated README explains each output file and lists the columns found in the actual output tables.

## Tools Required

The pipeline expects these tools in `PATH`:

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

Python requires:

```text
pandas
```

## Project-Specific Air Scripts

The root folder also contains project-specific scripts for the air/non-hybrid analysis. For new projects, prefer using `virome_pipeline_template/` instead of editing those air-specific scripts.

