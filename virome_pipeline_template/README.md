# Configurable Virome Identification Pipeline Template

This folder is a reusable template for viral-contig identification from metagenomic assemblies and paired reads.

The intended workflow is:

1. Edit `config.env`.
2. Fill `samples.tsv`.
3. Run `bash run_virome_pipeline.sh`.

The pipeline uses:

- VirSorter2
- geNomad
- CheckV
- Prodigal
- DIAMOND `blastx`
- DIAMOND `blastp`
- minimap2, samtools, mosdepth, seqkit, bbduk.sh

Viral calling is evidence based:

```text
putative_viral = vs2_hit OR genomad_hit OR blastx_hit OR blastp_hit
viral_confidence =
  High    if evidence_count >= 2
  Medium  if evidence_count == 1
  None    if evidence_count == 0
```

CheckV is reported in the final tables but is not used to remove contigs.

## Required Input

Edit `samples.tsv` with one row per sample:

```text
sample	assembly	r1	r2
sample01	/path/to/MEGAHIT-sample01.contigs.fa.gz	/path/to/sample01_R1.fastq.gz	/path/to/sample01_R2.fastq.gz
```

Assemblies may be `.fa`, `.fasta`, `.fa.gz`, or `.fasta.gz`.

## Main Command

```bash
bash run_virome_pipeline.sh
```

## Main Outputs

All project-level outputs are written to:

```text
${OUT_BASE}/merged/
```

Important files:

- `all_contig_read_counts.tsv`
- `diamond_blastx_best_hits.tsv`
- `diamond_blastp_best_hits.tsv`
- `all_viral_evidence.tsv`
- `virus_taxonomy_map.tsv`
- `virus_reads_matrix.tsv`
- `virus_reads_matrix_with_taxonomy.tsv`
- `virus_presence_absence_matrix.tsv`
- `virus_reads_by_sample_long.tsv`
- `final_virus_contigs.tsv`
- `final_virus_contigs_bacteriophage.tsv`
- `final_virus_contigs_eukaryotic.tsv`
- `final_virus_contigs_archaeal.tsv`
- `final_virus_contigs_unidentified.tsv`
- `final_virus_summary.txt`

The pipeline writes an output README to:

```text
${OUT_BASE}/README.md
```

That generated README describes each output file and lists columns found in the actual TSV files.
