#!/bin/bash
# Script: run_zover_local_blast.sh
# Reusable local BLASTN/BLASTX verification against ZOVER references.
# Purpose: BLAST extracted viral contigs against ZOVER references, summarize all
#          hits, make a blastx-primary verified virome profile, and prepare
#          candidate files for downstream phylogenetic tree confirmation.
#
# Usage:
#   bash run_zover_local_blast.sh <query_fasta> <zover_nucl_fasta[.gz]> <zover_prot_fasta_or_enriched_fasta[.gz]> <outdir> [threads]
#
# An NCBI-enriched protein FASTA adds organism/title text to BLASTX stitle.
#
# Verification thresholds can be changed without editing the script:
#   VERIFY_MAX_EVALUE=1e-10 VERIFY_MIN_BLASTX_PIDENT=30 VERIFY_MIN_BLASTX_QCOVS=50 bash run_zover_local_blast.sh ...
# Existing shared BLAST databases can be supplied to avoid rebuilding them per sample:
#   ZOVER_NUCL_DB=/path/to/zover_nucl_db ZOVER_PROT_DB=/path/to/zover_prot_db bash ...
set -euo pipefail

if [[ $# -lt 4 ]]; then
    echo "Usage: bash $0 <query_fasta> <zover_nucl_fasta[.gz]> <zover_prot_fasta_or_enriched_fasta[.gz]> <outdir> [threads]" >&2
    exit 1
fi

QUERY_FASTA="$1"
ZOVER_NUCL_IN="$2"
ZOVER_PROT_IN="$3"
OUTDIR="$4"
THREADS="${5:-8}"
VERIFY_MAX_EVALUE="${VERIFY_MAX_EVALUE:-1e-5}"
VERIFY_MIN_BLASTX_PIDENT="${VERIFY_MIN_BLASTX_PIDENT:-25}"
VERIFY_MIN_BLASTX_QCOVS="${VERIFY_MIN_BLASTX_QCOVS:-30}"

for input in "$QUERY_FASTA" "$ZOVER_NUCL_IN" "$ZOVER_PROT_IN"; do
    if [[ ! -s "$input" ]]; then
        echo "ERROR: FASTA is missing or empty: $input" >&2
        exit 1
    fi
done
for command in python3 blastn blastx makeblastdb; do
    if ! command -v "$command" >/dev/null 2>&1; then
        echo "ERROR: required command not found: $command" >&2
        exit 1
    fi
done

mkdir -p "$OUTDIR"

prep_ref() {
    local infile="$1"
    local label="$2"
    if [[ "$infile" == *.gz ]]; then
        local out="$OUTDIR/${label}.fasta"
        echo "  Decompressing $infile -> $out" >&2
        gunzip -c "$infile" > "$out"
        echo "$out"
    else
        echo "$infile"
    fi
}

run_makeblastdb() {
    local fasta="$1"
    local dbtype="$2"
    local dbout="$3"
    local log="$4"

    if ! makeblastdb -in "$fasta" -dbtype "$dbtype" -out "$dbout" > "$log" 2>&1; then
        echo "ERROR: makeblastdb failed for $dbtype reference. Last log lines:" >&2
        tail -n 40 "$log" >&2
        exit 1
    fi
}

echo "[1/7] Preparing reference FASTA files ..."
ZOVER_NUCL=$(prep_ref "$ZOVER_NUCL_IN" "zover_nucl_ref")
ZOVER_PROT=$(prep_ref "$ZOVER_PROT_IN" "zover_prot_ref")

if grep -m 1 '^>' "$ZOVER_PROT" | grep -q 'ncbi_title='; then
    echo "  Protein reference appears enriched with NCBI titles: $ZOVER_PROT"
else
    echo "  WARNING: Protein reference does not appear to contain ncbi_title= in the FASTA header." >&2
    echo "  blastx stitle will only show the original ZOVER protein header unless you use the enriched FASTA." >&2
fi

echo "  Verification thresholds: blastx evalue <= $VERIFY_MAX_EVALUE, pident >= $VERIFY_MIN_BLASTX_PIDENT, qcovs >= $VERIFY_MIN_BLASTX_QCOVS"

NUCL_DB="${ZOVER_NUCL_DB:-$OUTDIR/zover_nucl_db}"
PROT_DB="${ZOVER_PROT_DB:-$OUTDIR/zover_prot_db}"
FMT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle"

if [[ -n "${ZOVER_NUCL_DB:-}" ]]; then
    if [[ ! -s "${NUCL_DB}.nin" || ! -s "${NUCL_DB}.nsq" ]]; then
        echo "ERROR: shared nucleotide BLAST database is incomplete: $NUCL_DB" >&2
        exit 1
    fi
    echo "[2/7] Using shared nucleotide BLAST db: $NUCL_DB"
else
    echo "[2/7] Building nucleotide BLAST db from $ZOVER_NUCL ..."
    run_makeblastdb "$ZOVER_NUCL" nucl "$NUCL_DB" "$OUTDIR/makeblastdb_nucl.log"
fi

if [[ -n "${ZOVER_PROT_DB:-}" ]]; then
    if [[ ! -s "${PROT_DB}.pin" || ! -s "${PROT_DB}.psq" ]]; then
        echo "ERROR: shared protein BLAST database is incomplete: $PROT_DB" >&2
        exit 1
    fi
    echo "[3/7] Using shared protein BLAST db: $PROT_DB"
else
    echo "[3/7] Building protein BLAST db from $ZOVER_PROT ..."
    run_makeblastdb "$ZOVER_PROT" prot "$PROT_DB" "$OUTDIR/makeblastdb_prot.log"
fi

echo "[4/7] Running blastn (close-match search, nucl vs nucl) ..."
blastn -query "$QUERY_FASTA" -db "$NUCL_DB" -evalue 1e-5 -max_target_seqs 5 \
    -num_threads "$THREADS" -outfmt "$FMT" -out "$OUTDIR/blastn_results.tsv"

echo "[5/7] Running blastx (divergent-match search, translated nucl vs protein) ..."
blastx -query "$QUERY_FASTA" -db "$PROT_DB" -evalue 1e-5 -max_target_seqs 5 \
    -num_threads "$THREADS" -outfmt "$FMT" -out "$OUTDIR/blastx_results.tsv"

echo "[6/7] Selecting best hit per contig per method, merging, flagging no-hit contigs ..."
# outfmt 6 field order: 1 qseqid 2 sseqid 3 pident 4 length 5 mismatch
# 6 gapopen 7 qstart 8 qend 9 sstart 10 send 11 evalue 12 bitscore
# 13 qcovs 14 stitle
for method in blastn blastx; do
    sort -t$'\t' -k1,1 -k11,11g -k12,12gr "$OUTDIR/${method}_results.tsv" | \
        awk -F'\t' -v OFS='\t' -v m="$method" '!seen[$1]++ {print $1, m, $2, $3, $13, $11, $12, $14}' \
        > "$OUTDIR/${method}_best_hit.tsv"
done
# columns now: contig  method  sseqid  pident  qcovs  evalue  bitscore  stitle

{
    echo -e "contig\tblastn_hit\tblastn_accessions\tblastn_pident\tblastn_qcovs\tblastn_evalue\tblastn_bitscore\tblastn_title\tblastn_gene_orf\tblastn_ref_title_by_accession\tblastx_hit\tblastx_accessions\tblastx_pident\tblastx_qcovs\tblastx_evalue\tblastx_bitscore\tblastx_title\tblastx_gene_orf\tblastx_ref_title_by_accession\tany_hit"
    python3 - "$QUERY_FASTA" "$OUTDIR/blastn_best_hit.tsv" "$OUTDIR/blastx_best_hit.tsv" "$ZOVER_NUCL" "$ZOVER_PROT" << 'PYEOF'
import re
import sys

query_fasta, blastn_tsv, blastx_tsv, zover_nucl, zover_prot = sys.argv[1:6]

accession_re = re.compile(r"\b[A-Z]{1,4}_?\d{5,}(?:\.\d+)?\b")

def accession_candidates(text):
    found = []
    for m in re.finditer(r"\|([^|\s()]+)\|", text):
        found.extend(accession_re.findall(m.group(1)))
    for m in re.finditer(r"\(([^)]*)\)", text):
        first = m.group(1).split(",", 1)[0].strip()
        found.extend(accession_re.findall(first))
    found.extend(accession_re.findall(text))

    unique = []
    seen = set()
    for acc in found:
        if acc not in seen:
            unique.append(acc)
            seen.add(acc)
    return unique

def gene_orf_from_title(title):
    text = title or ""
    # Enriched headers append metadata after pipes; product/gene is in the original ZOVER part.
    base = re.split(r"\s+\|\s+(?:ncbi_nuccore|organism|ncbi_title)=", text, maxsplit=1)[0].strip()
    if ")" in base:
        product = base.split(")", 1)[1].strip()
    elif "|" in base:
        product = base.rsplit("|", 1)[-1].strip()
    else:
        product = base
    product = re.sub(r"\s+", " ", product).strip(" ;|")
    return product

def load_ref_titles(paths):
    mapping = {}
    for source, path in paths:
        with open(path) as f:
            for line in f:
                if not line.startswith(">"):
                    continue
                title = line[1:].strip()
                for acc in accession_candidates(title):
                    old = mapping.get(acc)
                    if old is None or (source == "nucl" and old[0] != "nucl"):
                        mapping[acc] = (source, title)
    return {acc: title for acc, (_source, title) in mapping.items()}

ref_titles = load_ref_titles([("prot", zover_prot), ("nucl", zover_nucl)])

all_contigs = []
with open(query_fasta) as f:
    for line in f:
        if line.startswith(">"):
            all_contigs.append(line[1:].strip().split()[0])

def load_best_hits(path):
    d = {}
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            c, method, sseqid, pident, qcovs, evalue, bitscore = parts[:7]
            title = "\t".join(parts[7:])
            accs = accession_candidates(" ".join([sseqid, title]))
            ref_matches = []
            for acc in accs:
                title_match = ref_titles.get(acc)
                if title_match and title_match not in ref_matches:
                    ref_matches.append(title_match)
            d[c] = {
                "hit": sseqid,
                "accs": ";".join(accs),
                "pident": pident,
                "qcovs": qcovs,
                "evalue": evalue,
                "bitscore": bitscore,
                "title": title,
                "gene_orf": gene_orf_from_title(title),
                "ref_title": " | ".join(ref_matches),
            }
    return d

bn = load_best_hits(blastn_tsv)
bx = load_best_hits(blastx_tsv)

def fields(hit):
    if not hit:
        return ["", "", "", "", "", "", "", "", ""]
    return [hit["hit"], hit["accs"], hit["pident"], hit["qcovs"], hit["evalue"], hit["bitscore"], hit["title"], hit["gene_orf"], hit["ref_title"]]

for c in all_contigs:
    b = bn.get(c)
    x = bx.get(c)
    any_hit = "YES" if (b or x) else "NO_HIT"
    print("\t".join([c] + fields(b) + fields(x) + [any_hit]))
PYEOF
} > "$OUTDIR/summary_all_contigs.tsv"

STATS_TSV="$OUTDIR/summary_stats.tsv"
METHOD_STATS_TSV="$OUTDIR/method_hit_stats.tsv"
TOP_TITLES_TSV="$OUTDIR/top_hit_titles.tsv"
TOP_ORGANISMS_TSV="$OUTDIR/top_hit_organisms.tsv"
VERIFIED_CONTIGS_TSV="$OUTDIR/verified_blastx_contigs.tsv"
VIROME_PROFILE_TSV="$OUTDIR/verified_virome_profile.tsv"
TREE_CANDIDATES_TSV="$OUTDIR/tree_candidates.tsv"
TREE_CONTIGS_FASTA="$OUTDIR/tree_candidate_contigs.fna"
TREE_REFS_FASTA="$OUTDIR/tree_candidate_reference_proteins.faa"
FULL_REPORT_TSV="$OUTDIR/zover_blast_full_report.tsv"

echo "[7/7] Writing summary statistics, verified virome profile, and tree candidate files ..."
python3 - \
    "$OUTDIR/summary_all_contigs.tsv" \
    "$STATS_TSV" \
    "$METHOD_STATS_TSV" \
    "$TOP_TITLES_TSV" \
    "$TOP_ORGANISMS_TSV" \
    "$FULL_REPORT_TSV" \
    "$VERIFIED_CONTIGS_TSV" \
    "$VIROME_PROFILE_TSV" \
    "$TREE_CANDIDATES_TSV" \
    "$TREE_CONTIGS_FASTA" \
    "$TREE_REFS_FASTA" \
    "$QUERY_FASTA" \
    "$ZOVER_PROT" \
    "$VERIFY_MAX_EVALUE" \
    "$VERIFY_MIN_BLASTX_PIDENT" \
    "$VERIFY_MIN_BLASTX_QCOVS" << 'PYEOF'
import csv
import re
import statistics
import sys
from collections import Counter, defaultdict

(
    summary_tsv,
    stats_tsv,
    method_stats_tsv,
    top_titles_tsv,
    top_organisms_tsv,
    full_report_tsv,
    verified_contigs_tsv,
    virome_profile_tsv,
    tree_candidates_tsv,
    tree_contigs_fasta,
    tree_refs_fasta,
    query_fasta,
    zover_prot,
    max_evalue_s,
    min_blastx_pident_s,
    min_blastx_qcovs_s,
) = sys.argv[1:17]

max_evalue = float(max_evalue_s)
min_blastx_pident = float(min_blastx_pident_s)
min_blastx_qcovs = float(min_blastx_qcovs_s)

with open(summary_tsv, newline="") as f:
    rows = list(csv.DictReader(f, delimiter="\t"))

def has(row, key):
    return bool(row.get(key, "").strip())

def num(row, key):
    value = row.get(key, "").strip()
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None

def avg(values):
    values = [v for v in values if v is not None]
    return "" if not values else f"{statistics.mean(values):.3f}"

def median(values):
    values = [v for v in values if v is not None]
    return "" if not values else f"{statistics.median(values):.3f}"

def pct(n, d):
    return "0.00" if d == 0 else f"{(100.0 * n / d):.2f}"

def clean_text(text):
    return " ".join((text or "").replace("\t", " ").split())

def organism_from_text(text):
    m = re.search(r"(?:^|\|)\s*organism=([^|]+)", text or "")
    if m:
        return clean_text(m.group(1))
    return ""

def ncbi_title_from_text(text):
    m = re.search(r"(?:^|\|)\s*ncbi_title=([^|]+)", text or "")
    if m:
        return clean_text(m.group(1))
    return ""

def profile_name(row):
    text = " | ".join([row.get("blastx_title", ""), row.get("blastx_ref_title_by_accession", "")])
    return organism_from_text(text) or ncbi_title_from_text(text) or clean_text(row.get("blastx_ref_title_by_accession")) or clean_text(row.get("blastx_title")) or row.get("blastx_hit", "")

def verification(row):
    bx_evalue = num(row, "blastx_evalue")
    bx_pident = num(row, "blastx_pident")
    bx_qcovs = num(row, "blastx_qcovs")
    if not has(row, "blastx_hit"):
        return "NO_VERIFIED_HIT", "NO_BLASTX", "no blastx hit"
    failed = []
    if bx_evalue is None or bx_evalue > max_evalue:
        failed.append("evalue")
    if bx_pident is None or bx_pident < min_blastx_pident:
        failed.append("pident")
    if bx_qcovs is None or bx_qcovs < min_blastx_qcovs:
        failed.append("qcovs")
    if failed:
        return "CHECK", "LOW", "blastx present but below threshold: " + ",".join(failed)
    if has(row, "blastn_hit"):
        return "VERIFIED", "HIGH", "blastx passes thresholds; blastn support present"
    return "VERIFIED", "MEDIUM", "blastx passes thresholds; no blastn support"

def read_fasta(path):
    records = {}
    order = []
    name = None
    chunks = []
    def flush():
        if name is None:
            return
        seq = "".join(chunks)
        records[name] = seq
        order.append(name)
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush()
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        flush()
    return records, order

def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))

def safe_id(text):
    return re.sub(r"[^A-Za-z0-9_.|:-]+", "_", text).strip("_") or "seq"

def first_token(header):
    return header.split()[0]

for row in rows:
    status, confidence, reason = verification(row)
    row["verification_status"] = status
    row["confidence"] = confidence
    row["verification_reason"] = reason
    row["primary_method"] = "blastx" if has(row, "blastx_hit") else ""
    row["profile_name"] = profile_name(row)

total = len(rows)
blastn_hits = sum(1 for r in rows if has(r, "blastn_hit"))
blastx_hits = sum(1 for r in rows if has(r, "blastx_hit"))
both_hits = sum(1 for r in rows if has(r, "blastn_hit") and has(r, "blastx_hit"))
blastn_only = sum(1 for r in rows if has(r, "blastn_hit") and not has(r, "blastx_hit"))
blastx_only = sum(1 for r in rows if has(r, "blastx_hit") and not has(r, "blastn_hit"))
any_hits = sum(1 for r in rows if r.get("any_hit") == "YES")
no_hits = sum(1 for r in rows if r.get("any_hit") == "NO_HIT")
verified = [r for r in rows if r["verification_status"] == "VERIFIED"]
check = [r for r in rows if r["verification_status"] == "CHECK"]
high = [r for r in verified if r["confidence"] == "HIGH"]
medium = [r for r in verified if r["confidence"] == "MEDIUM"]

stats = [
    ("total_contigs", total),
    ("verification_primary_method", "blastx"),
    ("verify_max_evalue", max_evalue_s),
    ("verify_min_blastx_pident", min_blastx_pident_s),
    ("verify_min_blastx_qcovs", min_blastx_qcovs_s),
    ("verified_contigs", len(verified)),
    ("verified_contigs_percent", pct(len(verified), total)),
    ("verified_high_confidence_contigs", len(high)),
    ("verified_medium_confidence_contigs", len(medium)),
    ("check_contigs_blastx_below_threshold", len(check)),
    ("contigs_with_any_zover_hit", any_hits),
    ("contigs_with_any_zover_hit_percent", pct(any_hits, total)),
    ("contigs_with_no_zover_hit", no_hits),
    ("contigs_with_no_zover_hit_percent", pct(no_hits, total)),
    ("contigs_with_blastn_hit", blastn_hits),
    ("contigs_with_blastn_hit_percent", pct(blastn_hits, total)),
    ("contigs_with_blastx_hit", blastx_hits),
    ("contigs_with_blastx_hit_percent", pct(blastx_hits, total)),
    ("contigs_with_both_blastn_and_blastx_hit", both_hits),
    ("contigs_with_blastn_only_hit", blastn_only),
    ("contigs_with_blastx_only_hit", blastx_only),
    ("blastn_mean_pident", avg(num(r, "blastn_pident") for r in rows)),
    ("blastn_median_pident", median(num(r, "blastn_pident") for r in rows)),
    ("blastn_mean_qcovs", avg(num(r, "blastn_qcovs") for r in rows)),
    ("blastn_median_qcovs", median(num(r, "blastn_qcovs") for r in rows)),
    ("blastn_mean_bitscore", avg(num(r, "blastn_bitscore") for r in rows)),
    ("blastn_median_bitscore", median(num(r, "blastn_bitscore") for r in rows)),
    ("blastx_mean_pident", avg(num(r, "blastx_pident") for r in rows)),
    ("blastx_median_pident", median(num(r, "blastx_pident") for r in rows)),
    ("blastx_mean_qcovs", avg(num(r, "blastx_qcovs") for r in rows)),
    ("blastx_median_qcovs", median(num(r, "blastx_qcovs") for r in rows)),
    ("blastx_mean_bitscore", avg(num(r, "blastx_bitscore") for r in rows)),
    ("blastx_median_bitscore", median(num(r, "blastx_bitscore") for r in rows)),
]

with open(stats_tsv, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["stat", "value"])
    writer.writerows(stats)

method_rows = [
    ["blastn", blastn_hits, pct(blastn_hits, total), blastn_only, both_hits, avg(num(r, "blastn_pident") for r in rows), avg(num(r, "blastn_qcovs") for r in rows), avg(num(r, "blastn_bitscore") for r in rows)],
    ["blastx", blastx_hits, pct(blastx_hits, total), blastx_only, both_hits, avg(num(r, "blastx_pident") for r in rows), avg(num(r, "blastx_qcovs") for r in rows), avg(num(r, "blastx_bitscore") for r in rows)],
]
with open(method_stats_tsv, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["method", "hit_contigs", "hit_percent", "method_only_contigs", "shared_with_other_method", "mean_pident", "mean_qcovs", "mean_bitscore"])
    writer.writerows(method_rows)

title_counts = []
for method in ("blastn", "blastx"):
    counter = Counter()
    for r in rows:
        title = clean_text(r.get(f"{method}_ref_title_by_accession")) or clean_text(r.get(f"{method}_title"))
        if title:
            counter[title] += 1
    title_counts.extend([method, count, title] for title, count in counter.most_common(30))
with open(top_titles_tsv, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["method", "count", "title"])
    writer.writerows(title_counts)

organism_counts = []
for method in ("blastn", "blastx"):
    counter = Counter()
    for r in rows:
        text = " | ".join([r.get(f"{method}_title", ""), r.get(f"{method}_ref_title_by_accession", "")])
        organism = organism_from_text(text)
        if organism:
            counter[organism] += 1
    organism_counts.extend([method, count, organism] for organism, count in counter.most_common(30))
with open(top_organisms_tsv, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(["method", "count", "organism"])
    writer.writerows(organism_counts)

summary_fieldnames = list(rows[0].keys()) if rows else []
with open(summary_tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\t", lineterminator="\n", fieldnames=summary_fieldnames)
    writer.writeheader()
    writer.writerows(rows)

with open(verified_contigs_tsv, "w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\t", lineterminator="\n", fieldnames=summary_fieldnames)
    writer.writeheader()
    writer.writerows(verified)

profile_groups = defaultdict(list)
for r in verified:
    profile_groups[r["profile_name"]].append(r)

profile_rows = []
for name, group in sorted(profile_groups.items(), key=lambda kv: (-len(kv[1]), kv[0])):
    confidence_counter = Counter(r["confidence"] for r in group)
    representative = min(group, key=lambda r: num(r, "blastx_evalue") if num(r, "blastx_evalue") is not None else float("inf"))
    profile_rows.append([
        name,
        len(group),
        confidence_counter.get("HIGH", 0),
        confidence_counter.get("MEDIUM", 0),
        avg(num(r, "blastx_pident") for r in group),
        avg(num(r, "blastx_qcovs") for r in group),
        representative.get("blastx_evalue", ""),
        representative.get("blastx_bitscore", ""),
        representative.get("blastx_hit", ""),
        representative.get("blastx_accessions", ""),
        clean_text(representative.get("blastx_gene_orf", "")),
        clean_text(representative.get("blastx_title", "")),
        "tree recommended before final species-level claim",
    ])

profile_header = [
    "profile_name",
    "verified_contig_count",
    "high_confidence_contigs",
    "medium_confidence_contigs",
    "mean_blastx_pident",
    "mean_blastx_qcovs",
    "best_blastx_evalue",
    "representative_blastx_bitscore",
    "representative_blastx_hit",
    "representative_accessions",
    "representative_gene_orf",
    "representative_title",
    "finalization_note",
]
with open(virome_profile_tsv, "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(profile_header)
    writer.writerows(profile_rows)

query_records, query_order = read_fasta(query_fasta)
query_by_id = {first_token(h): (h, seq) for h, seq in query_records.items()}
prot_records, prot_order = read_fasta(zover_prot)
prot_by_id = {first_token(h): (h, seq) for h, seq in prot_records.items()}

with open(tree_candidates_tsv, "w", newline="") as out:
    fieldnames = [
        "contig", "profile_name", "confidence", "blastx_hit", "blastx_accessions",
        "blastx_pident", "blastx_qcovs", "blastx_evalue", "blastx_bitscore", "blastx_gene_orf", "blastx_title", "tree_note"
    ]
    writer = csv.DictWriter(out, delimiter="\t", lineterminator="\n", fieldnames=fieldnames)
    writer.writeheader()
    for r in verified:
        writer.writerow({
            "contig": r.get("contig", ""),
            "profile_name": r.get("profile_name", ""),
            "confidence": r.get("confidence", ""),
            "blastx_hit": r.get("blastx_hit", ""),
            "blastx_accessions": r.get("blastx_accessions", ""),
            "blastx_pident": r.get("blastx_pident", ""),
            "blastx_qcovs": r.get("blastx_qcovs", ""),
            "blastx_evalue": r.get("blastx_evalue", ""),
            "blastx_bitscore": r.get("blastx_bitscore", ""),
            "blastx_gene_orf": clean_text(r.get("blastx_gene_orf", "")),
            "blastx_title": clean_text(r.get("blastx_title", "")),
            "tree_note": "build tree per gene/virus group; do not mix unrelated proteins",
        })

with open(tree_contigs_fasta, "w") as out:
    for r in verified:
        contig = r.get("contig", "")
        record = query_by_id.get(contig)
        if not record:
            continue
        header, seq = record
        out.write(f">{safe_id(contig)} profile={safe_id(r.get('profile_name', ''))} confidence={r.get('confidence', '')}\n{wrap(seq)}\n")

seen_ref_ids = set()
with open(tree_refs_fasta, "w") as out:
    for r in verified:
        hit = r.get("blastx_hit", "")
        record = prot_by_id.get(hit)
        if not record or hit in seen_ref_ids:
            continue
        header, seq = record
        seen_ref_ids.add(hit)
        out.write(f">{safe_id(hit)} {clean_text(header)}\n{wrap(seq)}\n")

def write_section(out, section_name, header, data_rows):
    out.write(f"## {section_name}\n")
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(header)
    writer.writerows(data_rows)
    out.write("\n")

with open(full_report_tsv, "w", newline="") as out:
    write_section(out, "SUMMARY_STATS", ["stat", "value"], stats)
    write_section(out, "METHOD_HIT_STATS", ["method", "hit_contigs", "hit_percent", "method_only_contigs", "shared_with_other_method", "mean_pident", "mean_qcovs", "mean_bitscore"], method_rows)
    write_section(out, "VERIFIED_VIROME_PROFILE", profile_header, profile_rows)
    write_section(out, "TOP_HIT_TITLES", ["method", "count", "title"], title_counts)
    write_section(out, "TOP_HIT_ORGANISMS", ["method", "count", "organism"], organism_counts)
    out.write("## VERIFIED_BLASTX_CONTIGS\n")
    writer = csv.DictWriter(out, delimiter="\t", lineterminator="\n", fieldnames=summary_fieldnames)
    writer.writeheader()
    writer.writerows(verified)
    out.write("\n## ALL_CONTIGS\n")
    writer = csv.DictWriter(out, delimiter="\t", lineterminator="\n", fieldnames=summary_fieldnames)
    writer.writeheader()
    writer.writerows(rows)

print(f"Total contigs tested: {total}")
print(f"Verified blastx-primary contigs: {len(verified)} ({pct(len(verified), total)}%)")
print(f"High confidence verified contigs: {len(high)}; medium confidence verified contigs: {len(medium)}")
print(f"Check manually: {len(check)} blastx hits below threshold")
print(f"Contigs with NO hit to ZOVER at all: {no_hits} ({pct(no_hits, total)}%)")
print(f"blastn hits: {blastn_hits} ({pct(blastn_hits, total)}%); blastx hits: {blastx_hits} ({pct(blastx_hits, total)}%)")
print(f"Virome profile groups: {len(profile_rows)}")
PYEOF

echo ""
echo "Done. Main outputs:"
echo "  Full report:              $FULL_REPORT_TSV"
echo "  Summary table:            $OUTDIR/summary_all_contigs.tsv"
echo "  Verified blastx contigs:  $VERIFIED_CONTIGS_TSV"
echo "  Verified virome profile:  $VIROME_PROFILE_TSV"
echo "  Tree candidates table:    $TREE_CANDIDATES_TSV"
echo "  Tree candidate contigs:   $TREE_CONTIGS_FASTA"
echo "  Tree reference proteins:  $TREE_REFS_FASTA"
echo "  Summary stats:            $STATS_TSV"
echo "  Method stats:             $METHOD_STATS_TSV"
echo "  Top hit titles:           $TOP_TITLES_TSV"
echo "  Top hit organisms:        $TOP_ORGANISMS_TSV"

