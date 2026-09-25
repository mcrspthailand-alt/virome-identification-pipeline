#!/bin/bash
# Run ZOVER verification across samples in a merged candidate TSV.
#
# What it does:
#   1) Reads sample + contig from the candidate TSV.
#   2) Extracts candidate contigs from per-sample FASTA files.
#   3) Runs run_zover_local_blast.sh for each pool.
#   4) Combines all pool outputs.
#   5) Adds original virus identification columns back to contig-level outputs
#      with prefix original_ so original calls can be compared to ZOVER calls.
#
# Required: MERGED_TSV, SAMPLES_ROOT, ZOVER_NUCL, ZOVER_PROT.
# See README.md for configuration and outputs.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKDIR="${WORKDIR:-$PWD}"
MERGED_TSV="${MERGED_TSV:?Set MERGED_TSV to the merged candidate TSV}"
SAMPLES_ROOT="${SAMPLES_ROOT:?Set SAMPLES_ROOT to the per-sample results directory}"
ZOVER_SCRIPT="${ZOVER_SCRIPT:-$SCRIPT_DIR/run_zover_local_blast.sh}"
ZOVER_NUCL="${ZOVER_NUCL:?Set ZOVER_NUCL to nucleotide FASTA or FASTA.gz}"
ZOVER_PROT="${ZOVER_PROT:?Set ZOVER_PROT to protein FASTA or FASTA.gz}"
OUTDIR="${OUTDIR:-$WORKDIR/zover_verified}"
THREADS="${THREADS:-16}"
CONTINUE_ON_ERROR="${CONTINUE_ON_ERROR:-0}"
FASTA_CANDIDATES="${FASTA_CANDIDATES:-01_contigs/contigs.filtered.ge1000.fa,01_contigs/contigs.filtered.len500_999.fa,01_contigs/contigs.filtered.ge500.fa,01_contigs/contigs.raw.fa}"
SAMPLE_COLUMN="${SAMPLE_COLUMN:-sample}"
CONTIG_COLUMN="${CONTIG_COLUMN:-contig}"
OUTPUT_PREFIX="${OUTPUT_PREFIX:-ALL_samples}"
if [[ ! "$OUTPUT_PREFIX" =~ ^[A-Za-z0-9_-]+$ ]]; then
    echo "ERROR: OUTPUT_PREFIX must contain only letters, digits, underscores or hyphens" >&2
    exit 1
fi

mkdir -p "$OUTDIR/queries" "$OUTDIR/logs" "$OUTDIR/pools" "$OUTDIR/combined" "$OUTDIR/reference_db"

MANIFEST_TSV="$OUTDIR/query_manifest.tsv"
EXTRACT_REPORT_TSV="$OUTDIR/query_extraction_report.tsv"
RUN_STATUS_TSV="$OUTDIR/run_status.tsv"

check_file() {
    local path="$1"
    local label="$2"
    if [[ ! -s "$path" ]]; then
        echo "ERROR: $label not found or empty: $path" >&2
        exit 1
    fi
}

check_file "$MERGED_TSV" "final_virus_contigs.tsv"
check_file "$ZOVER_SCRIPT" "run_zover_local_blast.sh"
check_file "$ZOVER_NUCL" "ZOVER nucleotide FASTA"
check_file "$ZOVER_PROT" "ZOVER protein FASTA"
if [[ ! -d "$SAMPLES_ROOT" ]]; then
    echo "ERROR: sample FASTA directory not found: $SAMPLES_ROOT" >&2
    exit 1
fi

echo "Working directory: $WORKDIR"
echo "Candidate TSV:     $MERGED_TSV"
echo "Samples root:      $SAMPLES_ROOT"
echo "ZOVER nucleotide:  $ZOVER_NUCL"
echo "ZOVER protein:     $ZOVER_PROT"
echo "Output directory:  $OUTDIR"
echo "Threads:           $THREADS"
echo ""

REFERENCE_DB_DIR="$OUTDIR/reference_db"
ZOVER_NUCL_FASTA="$REFERENCE_DB_DIR/zover_nucl_ref.fasta"
ZOVER_PROT_FASTA="$REFERENCE_DB_DIR/zover_prot_ref.fasta"
ZOVER_NUCL_DB="$REFERENCE_DB_DIR/zover_nucl_db"
ZOVER_PROT_DB="$REFERENCE_DB_DIR/zover_prot_db"

prepare_reference() {
    local source="$1"
    local target="$2"
    if [[ ! -s "$target" || "$source" -nt "$target" ]]; then
        if [[ "$source" == *.gz ]]; then
            gzip -cd "$source" > "$target"
        else
            cp "$source" "$target"
        fi
    fi
}

prepare_reference "$ZOVER_NUCL" "$ZOVER_NUCL_FASTA"
prepare_reference "$ZOVER_PROT" "$ZOVER_PROT_FASTA"

if [[ ! -s "${ZOVER_NUCL_DB}.nin" || ! -s "${ZOVER_NUCL_DB}.nsq" || "$ZOVER_NUCL_FASTA" -nt "${ZOVER_NUCL_DB}.nsq" ]]; then
    echo "Building shared nucleotide BLAST database ..."
    if ! makeblastdb -in "$ZOVER_NUCL_FASTA" -dbtype nucl -out "$ZOVER_NUCL_DB" > "$OUTDIR/logs/makeblastdb_nucl.log" 2>&1; then
        tail -n 40 "$OUTDIR/logs/makeblastdb_nucl.log" >&2
        exit 1
    fi
fi

if [[ ! -s "${ZOVER_PROT_DB}.pin" || ! -s "${ZOVER_PROT_DB}.psq" || "$ZOVER_PROT_FASTA" -nt "${ZOVER_PROT_DB}.psq" ]]; then
    echo "Building shared protein BLAST database ..."
    if ! makeblastdb -in "$ZOVER_PROT_FASTA" -dbtype prot -out "$ZOVER_PROT_DB" > "$OUTDIR/logs/makeblastdb_prot.log" 2>&1; then
        tail -n 40 "$OUTDIR/logs/makeblastdb_prot.log" >&2
        exit 1
    fi
fi

echo "[1/4] Extracting per-pool query FASTA from final candidate TSV + sample FASTA files ..."
python3 - \
    "$MERGED_TSV" \
    "$SAMPLES_ROOT" \
    "$OUTDIR/queries" \
    "$MANIFEST_TSV" \
    "$EXTRACT_REPORT_TSV" \
    "$FASTA_CANDIDATES" \
    "$SAMPLE_COLUMN" \
    "$CONTIG_COLUMN" << 'PYEOF'
import csv
import os
import re
import sys
from collections import defaultdict

merged_tsv, samples_root, query_dir, manifest_tsv, report_tsv, fasta_candidates_s, sample_column, contig_column = sys.argv[1:9]
fasta_candidates = [x.strip() for x in fasta_candidates_s.split(",") if x.strip()]

wanted = defaultdict(list)
wanted_seen = defaultdict(set)
with open(merged_tsv, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter="\t")
    required = {contig_column, sample_column}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"ERROR: merged TSV lacks required columns: {', '.join(sorted(missing))}")
    for row in reader:
        contig = (row.get(contig_column) or "").strip()
        sample = (row.get(sample_column) or "").strip()
        if sample and not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", sample):
            raise SystemExit(f"ERROR: unsafe sample name in merged TSV: {sample!r}")
        if contig and sample and contig not in wanted_seen[sample]:
            wanted[sample].append(contig)
            wanted_seen[sample].add(contig)

def sample_key(sample):
    m = re.match(r"^(.*?)(\d+)$", sample)
    return (m.group(1), int(m.group(2))) if m else (sample, 0)

def open_fasta(path):
    name = None
    seq = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if name is not None:
            yield name, "".join(seq)

def first_token(header):
    return header.split()[0]

def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

def find_fastas(sample):
    sample_dir = os.path.join(samples_root, sample)
    found = []
    for rel in fasta_candidates:
        path = os.path.join(sample_dir, rel)
        if os.path.exists(path) and os.path.getsize(path) > 0:
            found.append(path)
    return found

manifest_rows = []
report_rows = []
os.makedirs(query_dir, exist_ok=True)

for sample in sorted(wanted, key=sample_key):
    ids = set(wanted[sample])
    fasta_paths = find_fastas(sample)
    out_fasta = os.path.join(query_dir, f"{sample}_merged_contigs.fa")
    missing_path = os.path.join(query_dir, f"{sample}_missing_contigs.txt")
    found = {}

    for fasta_path in fasta_paths:
        for header, seq in open_fasta(fasta_path):
            contig_id = first_token(header)
            if contig_id in ids and contig_id not in found:
                found[contig_id] = (header, seq, fasta_path)
        if len(found) == len(ids):
            break

    with open(out_fasta, "w") as out:
        for contig_id in wanted[sample]:
            record = found.get(contig_id)
            if record:
                header, seq, _source = record
                out.write(f">{header}\n{wrap(seq)}\n")

    missing_ids = sorted(ids - set(found))
    with open(missing_path, "w") as out:
        for contig in missing_ids:
            out.write(contig + "\n")

    used_fastas = sorted({record[2] for record in found.values()})
    source_label = ";".join(used_fastas) if used_fastas else ("FASTA_NOT_FOUND" if not fasta_paths else "NO_MATCHING_CONTIGS")
    manifest_rows.append([sample, out_fasta, source_label, len(ids), len(found), len(missing_ids), missing_path])
    report_rows.append([sample, len(ids), len(found), len(missing_ids), source_label, out_fasta, missing_path])

with open(manifest_tsv, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample", "query_fasta", "source_fasta", "n_requested_contigs", "n_found_contigs", "n_missing_contigs", "missing_contigs_file"])
    writer.writerows(manifest_rows)

with open(report_tsv, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["sample", "requested", "found", "missing", "source_fasta", "query_fasta", "missing_contigs_file"])
    writer.writerows(report_rows)

if sum(row[4] for row in manifest_rows) == 0:
    raise SystemExit("ERROR: no contigs were extracted; check FASTA paths and contig IDs")

print(f"Samples in merged TSV: {len(wanted)}")
print(f"Manifest: {manifest_tsv}")
print(f"Extraction report: {report_tsv}")
PYEOF

echo "[2/4] Running ZOVER verification per pool ..."
echo -e "sample\tstatus\tquery_fasta\toutdir\tlog" > "$RUN_STATUS_TSV"

tail -n +2 "$MANIFEST_TSV" | while IFS=$'\t' read -r sample query_fasta source_fasta n_requested n_found n_missing missing_file; do
    pool_out="$OUTDIR/pools/${sample}_zover_verified"
    log="$OUTDIR/logs/${sample}_zover.log"

    if [[ "${n_found:-0}" -eq 0 ]]; then
        echo "  [$sample] SKIP: no contigs found in source FASTA"
        echo -e "$sample\tSKIP_NO_CONTIGS\t$query_fasta\t$pool_out\t$log" >> "$RUN_STATUS_TSV"
        continue
    fi

    echo "  [$sample] Running ZOVER verification on $n_found contigs ..."
    if ZOVER_NUCL_DB="$ZOVER_NUCL_DB" ZOVER_PROT_DB="$ZOVER_PROT_DB" \
        bash "$ZOVER_SCRIPT" "$query_fasta" "$ZOVER_NUCL_FASTA" "$ZOVER_PROT_FASTA" "$pool_out" "$THREADS" > "$log" 2>&1; then
        echo -e "$sample\tOK\t$query_fasta\t$pool_out\t$log" >> "$RUN_STATUS_TSV"
    else
        echo "  [$sample] ERROR. See log: $log" >&2
        echo -e "$sample\tERROR\t$query_fasta\t$pool_out\t$log" >> "$RUN_STATUS_TSV"
        if [[ "$CONTINUE_ON_ERROR" != "1" ]]; then
            exit 1
        fi
    fi
done

echo "[3/4] Combining per-pool ZOVER outputs ..."
python3 - "$RUN_STATUS_TSV" "$OUTDIR/combined" "$OUTPUT_PREFIX" << 'PYEOF'
import csv
import os
import sys

run_status_tsv, combined_dir, prefix = sys.argv[1:4]
os.makedirs(combined_dir, exist_ok=True)

targets = [
    ("verified_blastx_contigs.tsv", f"{prefix}_verified_blastx_contigs.zover_only.tsv"),
    ("verified_virome_profile.tsv", f"{prefix}_verified_virome_profile.tsv"),
    ("summary_stats.tsv", f"{prefix}_summary_stats.tsv"),
    ("tree_candidates.tsv", f"{prefix}_tree_candidates.zover_only.tsv"),
    ("summary_all_contigs.tsv", f"{prefix}_summary_all_contigs.zover_only.tsv"),
]

with open(run_status_tsv, newline="") as f:
    statuses = list(csv.DictReader(f, delimiter="\t"))

for per_pool_name, combined_name in targets:
    combined_path = os.path.join(combined_dir, combined_name)
    wrote_header = False
    with open(combined_path, "w", newline="") as out:
        writer = None
        for row in statuses:
            if row.get("status") != "OK":
                continue
            sample = row["sample"]
            pool_out = row["outdir"]
            path = os.path.join(pool_out, per_pool_name)
            if not os.path.exists(path) or os.path.getsize(path) == 0:
                continue
            with open(path, newline="") as f:
                reader = csv.DictReader(f, delimiter="\t")
                fieldnames = ["sample"] + (reader.fieldnames or [])
                if not wrote_header:
                    writer = csv.DictWriter(out, delimiter="\t", lineterminator="\n", fieldnames=fieldnames)
                    writer.writeheader()
                    wrote_header = True
                for rec in reader:
                    writer.writerow({"sample": sample, **rec})

if not any(row.get("status") == "OK" for row in statuses):
    raise SystemExit("ERROR: no sample completed successfully; inspect run_status.tsv and logs")

print(f"Combined ZOVER-only outputs written to: {combined_dir}")
PYEOF

echo "[4/4] Adding original virus identification columns to contig-level outputs ..."
python3 - \
    "$MERGED_TSV" \
    "$OUTDIR/combined/${OUTPUT_PREFIX}_verified_blastx_contigs.zover_only.tsv" \
    "$OUTDIR/combined/${OUTPUT_PREFIX}_tree_candidates.zover_only.tsv" \
    "$OUTDIR/combined/${OUTPUT_PREFIX}_summary_all_contigs.zover_only.tsv" \
    "$OUTDIR/combined/${OUTPUT_PREFIX}_verified_blastx_contigs.with_original.tsv" \
    "$OUTDIR/combined/${OUTPUT_PREFIX}_tree_candidates.with_original.tsv" \
    "$OUTDIR/combined/${OUTPUT_PREFIX}_summary_all_contigs.with_original.tsv" \
    "$SAMPLE_COLUMN" "$CONTIG_COLUMN" << 'PYEOF'
import csv
import os
import sys

(
    merged_tsv,
    verified_zover_only,
    tree_zover_only,
    summary_zover_only,
    verified_with_original,
    tree_with_original,
    summary_with_original,
    sample_column,
    contig_column,
) = sys.argv[1:10]

original_by_key = {}
original_fields = []
with open(merged_tsv, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f, delimiter="\t")
    original_fields = reader.fieldnames or []
    for row in reader:
        sample = (row.get(sample_column) or "").strip()
        contig = (row.get(contig_column) or "").strip()
        if sample and contig and (sample, contig) not in original_by_key:
            original_by_key[(sample, contig)] = row

# Keep original columns in the same order, but avoid duplicating join fields.
original_payload_fields = [x for x in original_fields if x not in {sample_column, contig_column}]
original_prefixed_fields = ["original_" + x for x in original_payload_fields]

def merge_file(in_path, out_path):
    if not os.path.exists(in_path) or os.path.getsize(in_path) == 0:
        with open(out_path, "w") as out:
            out.write("")
        return

    with open(in_path, newline="") as f, open(out_path, "w", newline="") as out:
        reader = csv.DictReader(f, delimiter="\t")
        in_fields = reader.fieldnames or []
        fieldnames = in_fields + original_prefixed_fields
        writer = csv.DictWriter(out, delimiter="\t", lineterminator="\n", fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            sample = (row.get("sample") or "").strip()
            contig = (row.get("contig") or row.get("contig_id") or "").strip()
            original = original_by_key.get((sample, contig), {})
            merged = dict(row)
            for field in original_payload_fields:
                merged["original_" + field] = original.get(field, "")
            writer.writerow(merged)

merge_file(verified_zover_only, verified_with_original)
merge_file(tree_zover_only, tree_with_original)
merge_file(summary_zover_only, summary_with_original)

print(f"Original records loaded: {len(original_by_key)}")
print(f"Verified contigs with original columns: {verified_with_original}")
print(f"Tree candidates with original columns: {tree_with_original}")
print(f"All contigs with original columns: {summary_with_original}")
PYEOF

echo ""
echo "Done. Outputs are in: $OUTDIR"
echo ""
echo "Run QC files:"
echo "  $EXTRACT_REPORT_TSV"
echo "  $RUN_STATUS_TSV"
echo ""
echo "Main combined files with original virus identification columns included:"
echo "  $OUTDIR/combined/${OUTPUT_PREFIX}_verified_blastx_contigs.with_original.tsv"
echo "  $OUTDIR/combined/${OUTPUT_PREFIX}_tree_candidates.with_original.tsv"
echo "  $OUTDIR/combined/${OUTPUT_PREFIX}_summary_all_contigs.with_original.tsv"
echo ""
echo "Virome profile summary:"
echo "  $OUTDIR/combined/${OUTPUT_PREFIX}_verified_virome_profile.tsv"

