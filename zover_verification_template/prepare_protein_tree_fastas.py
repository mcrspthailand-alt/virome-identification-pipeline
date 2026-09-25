#!/usr/bin/env python3
"""Build amino-acid tree inputs from a completed ZOVER run.

BLASTX is repeated only for selected candidate contigs so its translated query
HSP (`qseq`) and reading frame are available. The original verification calls
and confidence labels are read from the supplied annotation table, not re-made.
"""

import argparse
import csv
import gzip
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


ACCESSION_RE = re.compile(r"(?<![A-Za-z0-9_])([A-Za-z]{1,5}_?\d{5,}(?:\.\d+)?)(?![A-Za-z0-9_.])")
SAFE_RE = re.compile(r"[^A-Za-z0-9_.-]+")
BLAST_FIELDS = [
    "qseqid", "sseqid", "pident", "length", "qstart", "qend", "sstart",
    "send", "evalue", "bitscore", "qframe", "qseq",
]
MANIFEST_FIELDS = [
    "sample", "contig", "query_id", "tree_group_id", "virus_family",
    "gene_orf", "ncbi_gene", "ncbi_cds_product", "virus_organism",
    "reference_protein_accession", "annotation_source", "confidence",
    "curation_decision", "status", "original_blastx_pident",
    "original_blastx_qcovs", "original_blastx_evalue", "original_blastx_bitscore",
    "hsp_pident", "hsp_aa_length", "hsp_qstart", "hsp_qend", "hsp_sstart",
    "hsp_send", "hsp_evalue", "hsp_bitscore", "hsp_qframe",
    "hsp_stop_count", "hsp_ambiguous_count", "hsp_note",
]


def read_tsv(path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Empty TSV: {path}")
        return reader.fieldnames, list(reader)


def write_tsv(path, fieldnames, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fasta_records(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8") as handle:
        header, pieces = None, []
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(pieces).upper()
                header, pieces = line[1:], []
            elif line:
                if header is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                pieces.append(line)
        if header is not None:
            yield header, "".join(pieces).upper()


def write_fasta(path, records):
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for pos in range(0, len(seq), 80):
                handle.write(seq[pos:pos + 80] + "\n")


def slug(text):
    return SAFE_RE.sub("_", text).strip("_") or "unknown"


def sample_sort_key(row):
    sample = row["sample"]
    match = re.search(r"(\d+)$", sample)
    return (sample[:match.start()] if match else sample,
            int(match.group(1)) if match else 0, row["contig"])


def load_candidates(path):
    fields, rows = read_tsv(path)
    required = {"sample", "contig", "virus_family", "gene_orf",
                "reference_protein_accession", "tree_group_id"}
    missing = required - set(fields)
    if missing:
        raise ValueError(f"Annotation table lacks: {', '.join(sorted(missing))}")
    seen = set()
    for row in rows:
        key = (row["sample"].strip(), row["contig"].strip())
        if not all(key) or key in seen:
            raise ValueError(f"Missing or duplicate sample/contig: {key}")
        seen.add(key)
        for field in ("virus_family", "gene_orf", "reference_protein_accession", "tree_group_id"):
            if not row[field].strip():
                raise ValueError(f"{field} is blank for {key}")
    if not rows:
        raise ValueError("Annotation table has no candidates")
    return sorted(rows, key=sample_sort_key)


def load_curation(path, candidates):
    if path is None:
        return {(row["sample"], row["contig"]): "include" for row in candidates}
    fields, rows = read_tsv(path)
    if not {"sample", "contig", "decision"}.issubset(fields):
        raise ValueError("Curation TSV needs sample, contig and decision columns")
    decisions = {}
    for row in rows:
        key = (row["sample"].strip(), row["contig"].strip())
        decision = row["decision"].strip().lower()
        if key in decisions or decision not in {"include", "exclude"}:
            raise ValueError(f"Duplicate or invalid curation decision for {key}: {decision}")
        decisions[key] = decision
    expected = {(row["sample"], row["contig"]) for row in candidates}
    if decisions.keys() != expected:
        missing = sorted(expected - decisions.keys())
        extra = sorted(decisions.keys() - expected)
        raise ValueError(f"Curation keys differ from annotation; missing={missing[:5]}, extra={extra[:5]}")
    return decisions


def load_queries(query_root, selected):
    by_sample = defaultdict(set)
    for row in selected:
        by_sample[row["sample"]].add(row["contig"])
    found = {}
    for sample, wanted in by_sample.items():
        path = query_root / "queries" / f"{sample}_merged_contigs.fa"
        if not path.is_file():
            raise FileNotFoundError(f"Missing query FASTA: {path}")
        for header, seq in fasta_records(path):
            contig = header.split()[0]
            if contig in wanted:
                key = (sample, contig)
                if key in found:
                    raise ValueError(f"Duplicate query FASTA ID: {key}")
                if not seq or re.search(r"[^ACGTRYSWKMBDHVN]", seq):
                    raise ValueError(f"Invalid nucleotide sequence: {key}")
                found[key] = seq
    missing = sorted((row["sample"], row["contig"]) for row in selected
                     if (row["sample"], row["contig"]) not in found)
    if missing:
        raise ValueError(f"Selected contigs missing from query FASTA: {missing[:10]}")
    return found


def load_references(path, selected):
    wanted = {row["reference_protein_accession"] for row in selected}
    found = {}
    for header, seq in fasta_records(path):
        token = header.split()[0]
        accessions = set(ACCESSION_RE.findall(token)) & wanted
        for accession in accessions:
            if accession in found and found[accession][1] != seq:
                raise ValueError(f"Reference accession has conflicting sequences: {accession}")
            found[accession] = (header, seq)
    missing = sorted(wanted - found.keys())
    if missing:
        raise ValueError(f"Reference proteins missing from {path}: {', '.join(missing[:20])}")
    return found


def run_logged(command, log_path):
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
    if process.returncode:
        tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-20:])
        raise RuntimeError(f"Command failed ({process.returncode}): {' '.join(command)}\n{tail}")


def load_matching_hsps(path, candidate_by_id, reference_by_id):
    best = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for line_no, line in enumerate(handle, 1):
            values = line.rstrip("\r\n").split("\t")
            if len(values) != len(BLAST_FIELDS):
                raise ValueError(f"BLASTX line {line_no} has {len(values)} fields; expected {len(BLAST_FIELDS)}")
            hit = dict(zip(BLAST_FIELDS, values))
            query_id = hit["qseqid"]
            if query_id.startswith("lcl|"):
                query_id = query_id[4:]
            subject_id = hit["sseqid"]
            if subject_id.startswith("lcl|"):
                subject_id = subject_id[4:]
            if query_id not in candidate_by_id:
                continue
            expected_subject = reference_by_id[candidate_by_id[query_id]["reference_protein_accession"]]
            if subject_id != expected_subject:
                continue
            score = (float(hit["bitscore"]), -float(hit["evalue"]))
            if query_id not in best or score > best[query_id][0]:
                best[query_id] = (score, hit)
    return {query_id: item[1] for query_id, item in best.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True, help="Existing ZOVER verification result directory")
    parser.add_argument("--annotation", type=Path, required=True, help="verified_contigs_gene_CDS_annotation.tsv")
    parser.add_argument("--protein-reference", type=Path, required=True, help="ZOVER-associated protein FASTA or FASTA.gz")
    parser.add_argument("--output", type=Path, required=True, help="New directory for protein tree inputs")
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--curation", type=Path, help="Optional complete include/exclude decision table")
    parser.add_argument("--blastx-tsv", type=Path, help="Existing 12-field BLASTX result for testing or resuming")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    for path in (args.root, args.annotation, args.protein_reference):
        if not path.exists():
            parser.error(f"Input does not exist: {path}")
    if args.output.exists():
        parser.error(f"Output already exists; choose a new directory: {args.output}")

    candidates = load_candidates(args.annotation)
    decisions = load_curation(args.curation, candidates)
    selected = [row for row in candidates if decisions[(row["sample"], row["contig"])] == "include"]
    if not selected:
        parser.error("No candidates selected for protein tree FASTA")
    queries = load_queries(args.root, selected)
    references = load_references(args.protein_reference, selected)
    print(f"Loaded {len(selected)} selected contigs and {len(references)} matched protein references", flush=True)

    output = args.output
    work = output / "audit"
    groups_dir = output / "groups"
    work.mkdir(parents=True)
    groups_dir.mkdir()

    id_by_key = {}
    used_ids = set()
    for row in candidates:
        key = (row["sample"], row["contig"])
        query_id = f"Q__{slug(row['sample'])}__{slug(row['contig'])}"
        if query_id in used_ids:
            raise ValueError(f"Query ID collision after sanitizing names: {key}")
        id_by_key[key] = query_id
        used_ids.add(query_id)
    reference_by_accession = {accession: f"REF__{slug(accession)}" for accession in references}
    candidate_by_id = {id_by_key[(row["sample"], row["contig"])]: row for row in selected}

    write_fasta(work / "candidate_nucleotide_queries.fna", [
        (id_by_key[(row["sample"], row["contig"])], queries[(row["sample"], row["contig"])])
        for row in selected
    ])
    write_fasta(work / "reference_protein_subset.faa", [
        (f"lcl|{reference_by_accession[accession]}", references[accession][1])
        for accession in sorted(references)
    ])
    write_tsv(output / "curation_template.tsv", ["sample", "contig", "virus_family", "gene_orf", "decision", "reason"], [
        {"sample": row["sample"], "contig": row["contig"], "virus_family": row["virus_family"],
         "gene_orf": row["gene_orf"], "decision": "", "reason": ""}
        for row in candidates
    ])

    blastx_path = work / "blastx_translated_hsps.tsv"
    if args.blastx_tsv:
        if not args.blastx_tsv.is_file():
            parser.error(f"--blastx-tsv does not exist: {args.blastx_tsv}")
        blastx_path.write_bytes(args.blastx_tsv.read_bytes())
    else:
        db_prefix = work / "candidate_reference_db"
        print("Building a small BLAST database of the matched protein references ...", flush=True)
        run_logged(["makeblastdb", "-in", str(work / "reference_protein_subset.faa"),
                    "-dbtype", "prot", "-parse_seqids", "-blastdb_version", "4",
                    "-out", str(db_prefix)], work / "makeblastdb.log")
        print(f"Running BLASTX for {len(selected)} selected contigs with {args.threads} threads ...", flush=True)
        run_logged(["blastx", "-query", str(work / "candidate_nucleotide_queries.fna"),
                    "-db", str(db_prefix), "-query_gencode", "1", "-evalue", "100",
                    "-max_target_seqs", str(len(references)), "-max_hsps", "1",
                    "-num_threads", str(args.threads), "-outfmt", "6 " + " ".join(BLAST_FIELDS),
                    "-out", str(blastx_path)], work / "blastx.log")

    hsps = load_matching_hsps(blastx_path, candidate_by_id, reference_by_accession)
    missing_hsps = sorted(set(candidate_by_id) - set(hsps))
    if missing_hsps:
        write_tsv(output / "missing_matching_hsps.tsv", ["query_id", "reference_accession"], [
            {"query_id": query_id,
             "reference_accession": candidate_by_id[query_id]["reference_protein_accession"]}
            for query_id in missing_hsps
        ])
        raise ValueError(f"No BLASTX HSP against the recorded reference for {len(missing_hsps)} candidates; see missing_matching_hsps.tsv")

    manifest = []
    group_queries = defaultdict(list)
    group_references = defaultdict(set)
    group_labels = {}
    for row in candidates:
        key = (row["sample"], row["contig"])
        query_id = id_by_key[key]
        group_id = row["tree_group_id"]
        label = (row["virus_family"], row["gene_orf"])
        if group_id in group_labels and group_labels[group_id] != label:
            raise ValueError(f"Tree group ID has inconsistent family/gene: {group_id}")
        group_labels[group_id] = label
        decision = decisions[key]
        result = {
            "sample": row["sample"], "contig": row["contig"], "query_id": query_id,
            "tree_group_id": group_id, "virus_family": row["virus_family"],
            "gene_orf": row["gene_orf"], "ncbi_gene": row.get("ncbi_gene", ""),
            "ncbi_cds_product": row.get("ncbi_cds_product", ""),
            "virus_organism": row.get("virus_organism", ""),
            "reference_protein_accession": row["reference_protein_accession"],
            "annotation_source": row.get("annotation_source", ""),
            "confidence": row.get("confidence", ""), "curation_decision": decision,
            "status": "excluded" if decision == "exclude" else "created",
            "original_blastx_pident": row.get("blastx_pident", ""),
            "original_blastx_qcovs": row.get("blastx_qcovs", ""),
            "original_blastx_evalue": row.get("blastx_evalue", ""),
            "original_blastx_bitscore": row.get("blastx_bitscore", ""),
        }
        if decision == "include":
            hit = hsps[query_id]
            aa_original = hit["qseq"].replace("-", "").upper()
            if not aa_original or re.search(r"[^A-Z*]", aa_original):
                raise ValueError(f"Invalid translated BLASTX qseq for {query_id}: {aa_original[:40]}")
            stop_count = aa_original.count("*")
            aa = aa_original.replace("*", "X")
            ambiguous_count = sum(letter not in "ACDEFGHIKLMNPQRSTVWY" for letter in aa)
            notes = []
            if stop_count:
                notes.append("stop_replaced_with_X")
            if len(aa) < 50:
                notes.append("short_hsp_under_50_aa")
            if ambiguous_count:
                notes.append("ambiguous_aa")
            original_score = row.get("blastx_bitscore", "").strip()
            if original_score:
                try:
                    if abs(float(hit["bitscore"]) - float(original_score)) > 0.1:
                        notes.append("hsp_bitscore_differs_from_original")
                except ValueError:
                    notes.append("original_bitscore_unparseable")
            result.update({
                "hsp_pident": hit["pident"], "hsp_aa_length": len(aa),
                "hsp_qstart": hit["qstart"], "hsp_qend": hit["qend"],
                "hsp_sstart": hit["sstart"], "hsp_send": hit["send"],
                "hsp_evalue": hit["evalue"], "hsp_bitscore": hit["bitscore"],
                "hsp_qframe": hit["qframe"], "hsp_stop_count": stop_count,
                "hsp_ambiguous_count": ambiguous_count, "hsp_note": ";".join(notes),
            })
            group_queries[group_id].append((query_id, aa))
            group_references[group_id].add(row["reference_protein_accession"])
        manifest.append(result)

    write_tsv(output / "protein_tree_manifest.tsv", MANIFEST_FIELDS, manifest)
    group_summaries = []
    used_folders = set()
    for group_id in sorted(group_queries):
        folder_name = slug(group_id)
        if folder_name in used_folders:
            raise ValueError(f"Tree group folder collision: {group_id}")
        used_folders.add(folder_name)
        group_folder = groups_dir / folder_name
        group_folder.mkdir()
        query_records = sorted(group_queries[group_id])
        ref_records = [
            (reference_by_accession[accession], references[accession][1])
            for accession in sorted(group_references[group_id])
        ]
        write_fasta(group_folder / "query_blastx_hsp_proteins.faa", query_records)
        write_fasta(group_folder / "reference_proteins.faa", ref_records)
        write_fasta(group_folder / "tree_input_unaligned.faa", query_records + ref_records)
        selected_rows = [row for row in manifest if row["tree_group_id"] == group_id and row["status"] == "created"]
        lengths = [int(row["hsp_aa_length"]) for row in selected_rows]
        group_summaries.append({
            "tree_group_id": group_id, "virus_family": group_labels[group_id][0],
            "gene_orf": group_labels[group_id][1], "query_contigs": len(query_records),
            "reference_proteins": len(ref_records),
            "pools_n": len({row["sample"] for row in selected_rows}),
            "pools": ";".join(sorted({row["sample"] for row in selected_rows})),
            "min_query_hsp_aa": min(lengths), "max_query_hsp_aa": max(lengths),
            "query_with_stop": sum(int(row["hsp_stop_count"]) > 0 for row in selected_rows),
            "query_under_50_aa": sum(length < 50 for length in lengths),
            "tree_input_fasta": str(Path("groups") / folder_name / "tree_input_unaligned.faa"),
        })
    group_fields = ["tree_group_id", "virus_family", "gene_orf", "query_contigs",
                    "reference_proteins", "pools_n", "pools", "min_query_hsp_aa",
                    "max_query_hsp_aa", "query_with_stop", "query_under_50_aa", "tree_input_fasta"]
    write_tsv(output / "protein_tree_group_summary.tsv", group_fields, group_summaries)
    write_tsv(output / "protein_tree_run_qc.tsv", ["metric", "value"], [
        {"metric": "automated_candidates", "value": len(candidates)},
        {"metric": "curation_supplied", "value": "YES" if args.curation else "NO"},
        {"metric": "selected_for_tree", "value": len(selected)},
        {"metric": "excluded_by_curation", "value": len(candidates) - len(selected)},
        {"metric": "translated_query_hsps_written", "value": sum(len(rows) for rows in group_queries.values())},
        {"metric": "family_gene_groups", "value": len(group_summaries)},
        {"metric": "unique_reference_proteins", "value": len(references)},
        {"metric": "query_with_stop", "value": sum(int(row.get("hsp_stop_count", 0)) > 0 for row in manifest)},
        {"metric": "query_under_50_aa", "value": sum(0 < int(row.get("hsp_aa_length", 0)) < 50 for row in manifest)},
        {"metric": "hsp_bitscore_differs_from_original", "value": sum(
            "hsp_bitscore_differs_from_original" in row.get("hsp_note", "") for row in manifest)},
    ])
    print(f"Protein tree FASTAs: {len(selected)} queries in {len(group_summaries)} groups -> {output}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, RuntimeError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)

