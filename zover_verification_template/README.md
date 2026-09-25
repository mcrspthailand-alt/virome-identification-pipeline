# ZOVER verification template

This stage follows the repository's virome-identification workflow. It checks
candidate contigs against a host-associated ZOVER nucleotide reference with
BLASTN and a protein reference with BLASTX. It retains the upstream call
columns and reports every candidate, including no-hits.

**New to this stage? Follow the [step-by-step run guide](RUN_GUIDE.md).** It
covers setup, input checks, the complete verification run, output review, and
troubleshooting. Steps 1-7 are sufficient for virus verification. NCBI
annotation, amino-acid tree inputs, and reference-DB inventory are optional
follow-up tasks; they do not run as part of the core verification command.

The scripts are reusable. Do not commit reference FASTAs, BLAST databases,
sample FASTAs, NCBI metadata, or result directories to this repository.

## Requirements

- Bash, Python 3.9+, NCBI BLAST+ (`makeblastdb`, `blastn`, `blastx`), gzip,
  `sort`, and `awk`.
- A tab-delimited candidate table with `sample` and `contig` columns (column
  names can be overridden). A final table from an upstream screen is suitable.
- A directory of per-sample FASTAs, one subdirectory per sample. The default
  search order is the PANDA/PANDASIA `01_contigs/` layout; override it for other
  layouts with `FASTA_CANDIDATES`.
- Host-associated ZOVER nucleotide and protein FASTAs (`.gz` accepted). An
  NCBI-title-enriched protein FASTA is useful for `stitle` and gene/ORF labels,
  but is not required for BLASTX.

## Run across samples

Keep input and output paths outside this repository. For example:

```bash
export MERGED_TSV=/data/project/final_virus_contigs.tsv
export SAMPLES_ROOT=/data/project/results/samples
export ZOVER_NUCL=/data/databases/Rodent-associated_viruses_nucleotide.fasta.gz
export ZOVER_PROT=/data/databases/Rodent-associated_viruses_protein.fasta.gz
export OUTDIR=/data/project/zover_verified
export OUTPUT_PREFIX=ALL_rat_pools
export THREADS=16
bash zover_verification_template/run_zover_verification.sh
```

`SAMPLE_COLUMN` and `CONTIG_COLUMN` default to `sample` and `contig`.
`FASTA_CANDIDATES` is a comma-separated list of paths relative to each sample
directory; its default is
`01_contigs/contigs.filtered.ge1000.fa,01_contigs/contigs.filtered.len500_999.fa,01_contigs/contigs.filtered.ge500.fa,01_contigs/contigs.raw.fa`.
Set `VERIFY_MAX_EVALUE`, `VERIFY_MIN_BLASTX_PIDENT`, and
`VERIFY_MIN_BLASTX_QCOVS` to adjust the BLASTX-based screening thresholds.
The defaults are `1e-5`, `25`, and `30`. Set `CONTINUE_ON_ERROR=1` only if a
partial multi-sample run is acceptable; otherwise the driver stops at an error.

For one query FASTA, run `run_zover_local_blast.sh QUERY NUCLEOTIDE PROTEIN
OUTDIR [THREADS]` directly. BLASTN supports close matches; BLASTX is the
primary divergent-sequence screen. Neither alone establishes a final viral
identification.

## Main results

`OUTDIR/query_extraction_report.tsv` checks requested, found, and missing
contigs per sample. `OUTDIR/run_status.tsv` records sample-level success or
failure. Detailed BLAST and summary outputs are under `OUTDIR/pools/`.

The `OUTDIR/combined/` directory contains:

- `<prefix>_summary_all_contigs.with_original.tsv`: one row per candidate,
  BLASTN/BLASTX best hits, `stitle`, gene/ORF, e-value, bitscore, and upstream
  fields prefixed `original_`.
- `<prefix>_verified_blastx_contigs.with_original.tsv`: BLASTX-screened
  contigs with upstream fields.
- `<prefix>_verified_virome_profile.tsv`: profile summary (not an independent
  manually confirmed species list).
- `<prefix>_tree_candidates.with_original.tsv`: candidate contigs for review.
- `<prefix>_summary_stats.tsv`: per-sample verification statistics.

Files with `.zover_only.tsv` omit the upstream columns. These results are
evidence for verification, not a substitute for checking false positives,
nonviral homologs, read support, genome context, and phylogeny.

## Optional: NCBI gene/CDS annotation

NCBI E-utilities access is required for the metadata fetch. Set `NCBI_EMAIL`
and optionally `NCBI_API_KEY` in the environment. This step may take time for
large reference sets and should be run once per reference release. If batches
fail, inspect `zover_metadata_fetch_errors.tsv` and rerun before interpreting
the inventory.

```bash
python3 zover_verification_template/fetch_ncbi_metadata.py \
  --fasta "$ZOVER_NUCL" --outdir "$OUTDIR/ncbi_metadata"
python3 zover_verification_template/build_verified_tree_groups.py \
  --candidates "$OUTDIR/combined/${OUTPUT_PREFIX}_tree_candidates.with_original.tsv" \
  --queries "$OUTDIR/queries" --protein-fasta "$ZOVER_PROT" \
  --metadata "$OUTDIR/ncbi_metadata/zover_nuccore_metadata.tsv" \
  --features "$OUTDIR/ncbi_metadata/zover_nuccore_features.tsv" \
  --outdir "$OUTDIR/postprocess"
```

The postprocessor writes
`postprocess/verified_contigs_gene_CDS_annotation.tsv`, a manual curation
template, and `postprocess/tree_groups/<family__gene>/` with candidate
nucleotide contigs and reference proteins. Review annotations before final
virus calls, especially when no `protein_id` match is available.

## Optional: amino-acid tree inputs

The grouped nucleotide contigs above are an extraction/audit step, not the
sequences for a protein tree. This command reruns BLASTX only for selected
candidates to recover their translated query HSPs and outputs matching query
peptides plus reference proteins grouped by family and gene/ORF:

```bash
python3 zover_verification_template/prepare_protein_tree_fastas.py \
  --root "$OUTDIR" \
  --annotation "$OUTDIR/postprocess/verified_contigs_gene_CDS_annotation.tsv" \
  --protein-reference "$ZOVER_PROT" \
  --output /data/project/zover_protein_tree_inputs --threads 16
```

By default this includes all automated candidates. For a final tree set, use
`--curation decisions.tsv` with **every** `sample,contig` row and a `decision`
of `include` or `exclude`. Each group contains `tree_input_unaligned.faa`,
`query_blastx_hsp_proteins.faa`, and `reference_proteins.faa`. Check homologous
regions and reference sampling, then align, trim, and infer/support trees with
the phylogenetic tools chosen for the virus family. The script does not infer
trees or assign species names.

## Optional: reference-DB inventory

Skip this step if you only need virus verification. If the NCBI metadata TSVs
were fetched as above, the command below summarizes families and annotated
gene/CDS marker categories in the downloaded ZOVER reference set:

```bash
python3 zover_verification_template/summarize_zover_db.py \
  --metadata "$OUTDIR/ncbi_metadata/zover_nuccore_metadata.tsv" \
  --features "$OUTDIR/ncbi_metadata/zover_nuccore_features.tsv" \
  --outdir "$OUTDIR/zover_db_summary"
```

This writes family and marker/gene tables plus SVG figures. Marker classes
are keyword summaries, **not** a curated hallmark-gene database.

## Scope and provenance

The verification logic was adapted from the PANDA rodent ZOVER analysis,
including shared reference databases across samples and the `-query_gencode`
BLASTX option for translated HSP extraction. Paths and output names here are
configurable; no host-specific data is bundled.

