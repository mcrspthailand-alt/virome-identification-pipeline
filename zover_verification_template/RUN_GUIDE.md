# Run ZOVER virus verification: step by step

This guide starts **after** an upstream virus-identification screen. You need a
merged candidate TSV and the assembled contig FASTAs for those candidates.
The essential workflow ends at Step 7. NCBI annotation, protein-tree inputs,
and reference-DB inventory are separate optional tasks.

Commands below are for a Linux/HPC Bash shell. Run them from the root of this
repository unless a step says otherwise. Paths in the example are placeholders;
replace them with paths on your system.

## What you need

| Input | Purpose |
| --- | --- |
| Candidate TSV | One row or more per candidate, with `sample` and `contig` columns. Other columns are retained in the combined result with an `original_` prefix. |
| Sample FASTAs | One directory per sample under a common root. Each directory must contain a FASTA with the contig IDs named in the TSV. |
| ZOVER nucleotide FASTA | Local BLASTN reference, plain FASTA or `.gz`. |
| ZOVER protein FASTA | Local BLASTX reference, plain FASTA or `.gz`. |
| Output directory | A new directory outside the Git repository. Have enough space for decompressed references, BLAST databases, and results. |

BLASTN searches nucleotide against nucleotide. BLASTX translates the candidate
contig and searches the protein reference. **The default screening decision is
BLASTX-primary.** This is a homology screen, not automatic confirmation of a
virus species or proof of infection.

## Step 1: get the branch on the HPC

For a new checkout:

```bash
git clone --branch pandasia/zover-virus-verification \
  https://github.com/mcrspthailand-alt/virome-identification-pipeline.git
cd virome-identification-pipeline
```

For an existing checkout of this same repository:

```bash
git fetch origin
git switch --track origin/pandasia/zover-virus-verification
```

If you already created a local branch with that name, switch to it and pull
the latest version instead of using `--track` again. Verify the script exists:

```bash
ls -l zover_verification_template/run_zover_verification.sh
```

The ZOVER step is separate from `run_virome_pipeline.sh`. If your merged
candidate TSV already exists, you do **not** need to rerun identification.

## Step 2: activate your environment and check tools

Use the environment where NCBI BLAST+ and Python are installed. For example:

```bash
conda activate virus_metagenomic
for tool in python3 makeblastdb blastn blastx gzip sort awk; do
  command -v "$tool" || exit 1
done
blastx -version
```

If `conda` is unavailable, use your site's module system or another
environment containing the same commands. The scripts need Python 3.9+.

## Step 3: check the TSV and FASTA layout

The TSV must be tab-delimited. The default headers are `sample` and `contig`:

```text
sample	contig	other_upstream_columns...
Rat-pool-1	k141_12345	...
Rat-pool-2	k141_67890	...
```

The first line above shows **tab-separated fields**; `\t` is shown only to
make the separators visible. Contig IDs must match the first token after `>`
in the FASTA header, including case. The default FASTA layout is:

```text
SAMPLES_ROOT/
  Rat-pool-1/01_contigs/contigs.filtered.ge1000.fa
  Rat-pool-1/01_contigs/contigs.filtered.ge500.fa
  Rat-pool-2/01_contigs/contigs.filtered.ge1000.fa
  ...
```

The driver searches these files, in order, within each sample directory:

```text
01_contigs/contigs.filtered.ge1000.fa
01_contigs/contigs.filtered.len500_999.fa
01_contigs/contigs.filtered.ge500.fa
01_contigs/contigs.raw.fa
```

You may use another layout by setting `FASTA_CANDIDATES` in Step 5. The
current extraction code reads **uncompressed sample FASTAs**. The ZOVER
reference FASTAs, unlike sample FASTAs, may be `.gz`.

## Step 4: choose input and output paths

Set every required path explicitly. Quote paths so spaces are handled.
`OUTDIR` should be new, outside the checkout, to avoid mixing runs or
committing results by accident.

```bash
export MERGED_TSV="/data/my_project/final_virus_contigs.tsv"
export SAMPLES_ROOT="/data/my_project/results/samples"
export ZOVER_NUCL="/data/references/Host-associated_viruses_nucleotide.fasta.gz"
export ZOVER_PROT="/data/references/Host-associated_viruses_protein.fasta.gz"
export OUTDIR="/data/my_project/zover_verified_run1"
export OUTPUT_PREFIX="ALL_samples"
export THREADS=16
```

For the PANDA rodent layout discussed in this project, the corresponding
values are:

```bash
export BASE="/home/panda/workspace/projects/nextflow_metagenomic/PANDA/rodent"
export MERGED_TSV="$BASE/virus_verification/final_virus_contigs.tsv"
export SAMPLES_ROOT="$BASE/05_virome_identification/results/samples"
export ZOVER_NUCL="$BASE/virus_verification/Rodent-associated_viruses_nucleotide.fasta.gz"
export ZOVER_PROT="$BASE/virus_verification/Rodent-associated_viruses_protein.fasta.gz"
export OUTDIR="$BASE/virus_verification/rodent_zover_verified_new_run"
export OUTPUT_PREFIX="ALL_rat_pools"
export THREADS=16
```

Only use the rodent block if those files actually exist at those paths. For
bat or another host, replace the table, sample directory, and both references.
**Do not use rodent references for a bat analysis by accident.**

## Step 5: preflight checks

Run these before the long BLAST step:

```bash
ls -lh "$MERGED_TSV" "$ZOVER_NUCL" "$ZOVER_PROT"
ls -ld "$SAMPLES_ROOT"
head -n 3 "$MERGED_TSV"
find "$SAMPLES_ROOT" -maxdepth 1 -mindepth 1 -type d | head
```

For gzipped references, also check integrity:

```bash
gzip -t "$ZOVER_NUCL" "$ZOVER_PROT"
```

Skip the `gzip -t` command for any reference that is not `.gz`.

If your TSV uses different headers, set them before running, for example:

```bash
export SAMPLE_COLUMN="pool"
export CONTIG_COLUMN="contig_id"
```

If your FASTA is in another location within each sample directory, set a
comma-separated list of **relative paths**, for example:

```bash
export FASTA_CANDIDATES="contigs/viral_candidates.fa,assembly/contigs.fa"
```

Do not set either override unless your files require it. `OUTPUT_PREFIX` is
used only in combined filenames and may contain letters, digits, `_`, or `-`.

## Step 6: run virus verification

From the repository root:

```bash
bash zover_verification_template/run_zover_verification.sh \
  > "$OUTDIR.run.log" 2>&1
```

Follow the log in another terminal:

```bash
tail -f "$OUTDIR.run.log"
```

The driver extracts only TSV-listed contigs, builds shared BLAST databases
once, runs BLASTN and BLASTX for each sample, and combines per-sample outputs.
It stops on a sample error by default. For a deliberately partial run, set
`CONTINUE_ON_ERROR=1`; inspect `run_status.tsv` before using combined results.
When rerunning, prefer a **new** `OUTDIR` so previous results are preserved.

For a **single already-extracted FASTA**, use the lower-level command instead:

```bash
bash zover_verification_template/run_zover_local_blast.sh \
  /data/query_contigs.fa "$ZOVER_NUCL" "$ZOVER_PROT" \
  /data/one_sample_zover 16
```

## Step 7: check the completed run

First check sample status and sequence extraction:

```bash
head -n 15 "$OUTDIR/run_status.tsv"
head -n 15 "$OUTDIR/query_extraction_report.tsv"
```

Every intended sample should show `OK`; `ERROR` means its log needs review.
For complete verification, `missing` should be `0` for each sample. A missing
contig means its ID was absent from the FASTAs searched for that sample; it
has **not** been verified. Check the report's `missing_contigs_file` path.

Open the main combined outputs:

```bash
ls -lh "$OUTDIR/combined"
head -n 3 "$OUTDIR/combined/${OUTPUT_PREFIX}_summary_all_contigs.with_original.tsv"
head -n 3 "$OUTDIR/combined/${OUTPUT_PREFIX}_verified_blastx_contigs.with_original.tsv"
head -n 3 "$OUTDIR/combined/${OUTPUT_PREFIX}_summary_stats.tsv"
```

| Output | Use |
| --- | --- |
| `combined/<prefix>_summary_all_contigs.with_original.tsv` | Main contig-level audit table, including no-hits, BLASTN/BLASTX best hit fields and original upstream fields. |
| `combined/<prefix>_verified_blastx_contigs.with_original.tsv` | Contigs passing the BLASTX screen. Review before finalizing. |
| `combined/<prefix>_verified_virome_profile.tsv` | Sample/profile summary of screened contigs. Do not treat as a final species list. |
| `combined/<prefix>_tree_candidates.with_original.tsv` | Candidate list for manual review and possible phylogeny. |
| `combined/<prefix>_summary_stats.tsv` | Per-sample counts and screening statistics. |
| `pools/<sample>_zover_verified/blastn_results.tsv` and `blastx_results.tsv` | Raw per-sample hits when a best-hit summary is not enough. |

In the per-contig table, inspect `blastx_hit`, `blastx_title`,
`blastx_gene_orf`, `blastx_pident`, `blastx_qcovs`, `blastx_evalue`, and
`blastx_bitscore`. The script labels a BLASTX hit as `VERIFIED` when it passes
its configured thresholds; `HIGH` adds BLASTN support and `MEDIUM` lacks it.
These are **pipeline evidence labels**, not biological confirmation. Look for
nonviral/host homologs, taxonomic ambiguity, weak or short matches, read
support, and genome context before reporting a virus. BLASTN-only hits are
still visible in the all-contigs table but do not pass the BLASTX-primary call.

The default thresholds are BLASTX e-value `<= 1e-5`, identity `>= 25%`, and
query coverage (`qcovs`) `>= 30%`. If changing thresholds, export
`VERIFY_MAX_EVALUE`, `VERIFY_MIN_BLASTX_PIDENT`, and/or
`VERIFY_MIN_BLASTX_QCOVS` before Step 6, and record the values in your
methods. A hit title may only contain a protein product; enriched NCBI titles
or the optional annotation step below can help resolve gene/CDS identity.

**Stop here if you only need virus verification.** No NCBI metadata fetch or
reference-DB summary is needed for Steps 1-7.

## Optional A: NCBI gene/CDS annotation and protein-tree FASTAs

Use this only when you need clearer gene/CDS labels or phylogenetic inputs.
It requires internet access to NCBI and may be slow for a large reference
set. The metadata fetch reads nucleotide accessions from `ZOVER_NUCL` and
writes reusable TSVs under `OUTDIR/ncbi_metadata`.

```bash
export NCBI_EMAIL="your.name@institute.org"
python3 zover_verification_template/fetch_ncbi_metadata.py \
  --fasta "$ZOVER_NUCL" --outdir "$OUTDIR/ncbi_metadata"
```

Inspect `zover_metadata_fetch_errors.tsv` and
`zover_accessions_not_returned.tsv` in that directory. An incomplete NCBI
fetch can leave gene/CDS assignments unresolved. Then annotate the candidate
contigs and create family/gene review groups:

```bash
python3 zover_verification_template/build_verified_tree_groups.py \
  --candidates "$OUTDIR/combined/${OUTPUT_PREFIX}_tree_candidates.with_original.tsv" \
  --queries "$OUTDIR/queries" --protein-fasta "$ZOVER_PROT" \
  --metadata "$OUTDIR/ncbi_metadata/zover_nuccore_metadata.tsv" \
  --features "$OUTDIR/ncbi_metadata/zover_nuccore_features.tsv" \
  --outdir "$OUTDIR/postprocess"
```

Review `OUTDIR/postprocess/verified_contigs_gene_CDS_annotation.tsv`,
`postprocess_query_sequence_qc.tsv`, and
`tree_groups/manual_curation_template.tsv`. The grouped `query_contigs.fna`
files are **nucleotide** audit/extraction files. For amino-acid tree inputs,
run:

```bash
python3 zover_verification_template/prepare_protein_tree_fastas.py \
  --root "$OUTDIR" \
  --annotation "$OUTDIR/postprocess/verified_contigs_gene_CDS_annotation.tsv" \
  --protein-reference "$ZOVER_PROT" \
  --output "${OUTDIR}_protein_tree_inputs" --threads "$THREADS"
```

This reruns BLASTX for selected contigs to recover translated query HSPs. The
new `groups/<family__gene>/tree_input_unaligned.faa` files contain query
peptides and reference proteins. They are **unaligned inputs**, not inferred
trees. By default all automated candidates are included. After reviewing
false positives, pass `--curation /path/to/decisions.tsv`; it must contain
`sample`, `contig`, and `decision` (`include` or `exclude`) for **every** row.
The tool writes a `curation_template.tsv` in its output for a subsequent run;
use a different output directory when rerunning after curation. Align and
inspect homologous regions, add suitable references/outgroups, and infer a
tree using an appropriate phylogenetic method.

## Optional B: summarize the ZOVER reference DB

This is for methods, tables, or figures about the **reference set itself**.
It is **not part of virus verification** and can be skipped even when doing
Optional A. It requires the NCBI metadata TSVs from the first command in
Optional A, but not the candidate-annotation or tree commands.

```bash
python3 zover_verification_template/summarize_zover_db.py \
  --metadata "$OUTDIR/ncbi_metadata/zover_nuccore_metadata.tsv" \
  --features "$OUTDIR/ncbi_metadata/zover_nuccore_features.tsv" \
  --outdir "$OUTDIR/zover_db_summary"
```

This writes family counts, gene/CDS and marker-category tables, plus SVG
figures. These are annotations of the downloaded references, **not** a
manually curated hallmark-gene database. Interpret family counts in the
context of the exact reference release used for the run.

## Common problems

| Symptom | Check |
| --- | --- |
| `Set MERGED_TSV` or another required-variable error | Export all four required paths in Step 4 in the same shell where you run Step 6. |
| `no contigs were extracted` or nonzero `missing` | Compare TSV `sample` with subdirectory names and TSV `contig` with FASTA IDs; set `FASTA_CANDIDATES` if the layout differs. |
| `makeblastdb` failure | Inspect `OUTDIR/logs/makeblastdb_*.log`; confirm reference FASTA integrity and write permission. |
| A sample shows `ERROR` | Read `OUTDIR/logs/<sample>_zover.log`; combined tables may be incomplete. |
| `blastx_title` lacks organism or gene | The original protein FASTA may have only a product label. Use enriched titles or Optional A; do not guess the species from the product. |
| Many `NO_VERIFIED_HIT` or `CHECK` rows | Examine raw BLASTX hits, thresholds, reference coverage, and sequence quality; absence of a ZOVER match is not proof of absence of a virus. |
| `--output already exists` for protein FASTAs | Pick a fresh output directory to preserve the previous run. |

For reproducibility, keep the input TSV, reference FASTA release/checksums,
BLAST+ version, thresholds, `run_status.tsv`, `query_extraction_report.tsv`,
and manual curation decisions together with the analysis results.

