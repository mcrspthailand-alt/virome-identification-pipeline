#!/usr/bin/env python3
"""Fetch nuccore taxonomy and gene/CDS feature annotations for a FASTA reference."""

import argparse
import csv
import gzip
import hashlib
import os
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path


ACCESSION_RE = re.compile(r"^(?:[A-Z]{1,4}_?\d{5,}|[A-Z]{2}\d{6,})(?:\.\d+)?$", re.I)
DB_TAG_RE = re.compile(r"(?:gb|ref|emb|dbj|tpg|tpe|tpd)\|([^|\s]+)\|", re.I)
EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

METADATA_FIELDS = [
    "input_accession", "accession_version", "organism", "definition",
    "taxonomy", "family", "record_status",
]
FEATURE_FIELDS = [
    "input_accession", "nuccore_accession", "organism", "family",
    "feature_key", "location", "gene", "product", "protein_id",
    "locus_tag", "note", "db_xref",
]


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def accessions_from_header(header):
    found = []
    for match in DB_TAG_RE.finditer(header):
        value = match.group(1).strip().rstrip(".,;")
        if ACCESSION_RE.fullmatch(value):
            found.append(value)
    if not found:
        token = header.split()[0] if header.split() else ""
        if ACCESSION_RE.fullmatch(token):
            found.append(token)
    return found


def fasta_accessions(path):
    seen = set()
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                for accession in accessions_from_header(line[1:].strip()):
                    if accession not in seen:
                        seen.add(accession)
                        yield accession


def norm_acc(value):
    return (value or "").strip().upper()


def base_acc(value):
    return re.sub(r"\.\d+$", "", norm_acc(value))


def child_text(parent, tag):
    node = parent.find(tag)
    return "" if node is None or node.text is None else " ".join(node.text.split())


def family_from_taxonomy(taxonomy):
    for item in (taxonomy or "").split(";"):
        item = item.strip()
        if item.lower().endswith("viridae"):
            return item
    return "Unclassified"


def parse_record(record):
    accession = child_text(record, "GBSeq_accession-version")
    primary = child_text(record, "GBSeq_primary-accession")
    taxonomy = child_text(record, "GBSeq_taxonomy")
    metadata = {
        "accession_version": accession or primary,
        "organism": child_text(record, "GBSeq_organism"),
        "definition": child_text(record, "GBSeq_definition"),
        "taxonomy": taxonomy,
        "family": family_from_taxonomy(taxonomy),
    }
    aliases = {norm_acc(accession), norm_acc(primary)} - {""}
    other_ids = record.find("GBSeq_other-seqids")
    if other_ids is not None:
        for item in other_ids.findall("GBSeqid"):
            value = (item.text or "").strip()
            for match in DB_TAG_RE.finditer(value + "|"):
                aliases.add(norm_acc(match.group(1)))
            if ACCESSION_RE.fullmatch(value):
                aliases.add(norm_acc(value))

    features = []
    table = record.find("GBSeq_feature-table")
    if table is not None:
        for feature in table.findall("GBFeature"):
            key = child_text(feature, "GBFeature_key")
            if key not in {"CDS", "gene"}:
                continue
            qualifiers = {}
            qroot = feature.find("GBFeature_quals")
            if qroot is not None:
                for qualifier in qroot.findall("GBQualifier"):
                    qname = child_text(qualifier, "GBQualifier_name")
                    qvalue = child_text(qualifier, "GBQualifier_value")
                    if qname and qvalue:
                        qualifiers.setdefault(qname, []).append(qvalue)

            def q(name):
                return "; ".join(qualifiers.get(name, []))

            features.append({
                "nuccore_accession": accession or primary,
                "organism": metadata["organism"],
                "family": metadata["family"],
                "feature_key": key,
                "location": child_text(feature, "GBFeature_location"),
                "gene": q("gene"),
                "product": q("product"),
                "protein_id": q("protein_id"),
                "locus_tag": q("locus_tag"),
                "note": q("note"),
                "db_xref": q("db_xref"),
            })
    return metadata, aliases, features


def fetch_batch(accessions, email, api_key, tool, timeout, retries, delay):
    params = {
        "db": "nuccore",
        "id": ",".join(accessions),
        "rettype": "gb",
        "retmode": "xml",
        "tool": tool,
    }
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    data = urllib.parse.urlencode(params).encode("ascii")
    request = urllib.request.Request(
        EFETCH_URL,
        data=data,
        headers={"User-Agent": f"{tool}/1.0 ({email or 'research contact not supplied'})"},
        method="POST",
    )
    last_error = None
    for attempt in range(retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                payload = response.read()
            time.sleep(delay)
            return payload
        except Exception as exc:
            last_error = exc
            if attempt < retries:
                time.sleep(min(60, 2 ** attempt))
    raise RuntimeError(str(last_error))


def write_tsv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, help="ZOVER-associated nucleotide FASTA, plain or .gz")
    parser.add_argument("--outdir", required=True, help="Directory for metadata, features, and request cache")
    parser.add_argument("--email", default=os.environ.get("NCBI_EMAIL", ""))
    parser.add_argument("--api-key", default=os.environ.get("NCBI_API_KEY", ""))
    parser.add_argument("--tool", default="ZOVER_Verification_Metadata")
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--delay", type=float, default=None)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--retries", type=int, default=4)
    args = parser.parse_args()

    if args.batch_size < 1 or args.batch_size > 200:
        parser.error("--batch-size must be from 1 to 200")
    delay = args.delay if args.delay is not None else (0.12 if args.api_key else 0.36)
    accessions = list(fasta_accessions(args.fasta))
    if not accessions:
        raise SystemExit("ERROR: no GenBank/RefSeq accessions were recognized in the FASTA headers")

    outdir = Path(args.outdir)
    cache_dir = outdir / "ncbi_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    records_by_request = {}
    features_by_request = {}
    failed_batches = []
    failed_accessions = set()
    batches = [accessions[i:i + args.batch_size] for i in range(0, len(accessions), args.batch_size)]

    for batch_no, batch in enumerate(batches, 1):
        digest = hashlib.sha256(",".join(batch).encode("utf-8")).hexdigest()[:20]
        cache_path = cache_dir / f"batch_{digest}.xml"
        try:
            if cache_path.exists() and cache_path.stat().st_size:
                payload = cache_path.read_bytes()
            else:
                payload = fetch_batch(batch, args.email, args.api_key, args.tool,
                                      args.timeout, args.retries, delay)
                cache_path.write_bytes(payload)
            root = ET.fromstring(payload)
            parsed = [parse_record(node) for node in root.findall(".//GBSeq")]
            by_alias = {}
            by_base = {}
            for metadata, aliases, feature_rows in parsed:
                for alias in aliases:
                    by_alias[alias] = (metadata, feature_rows)
                    by_base.setdefault(base_acc(alias), (metadata, feature_rows))
            for accession in batch:
                match = by_alias.get(norm_acc(accession)) or by_base.get(base_acc(accession))
                if match:
                    metadata, feature_rows = match
                    records_by_request[accession] = metadata
                    features_by_request[accession] = feature_rows
                else:
                    records_by_request[accession] = None
                    features_by_request[accession] = []
        except Exception as exc:
            failed_batches.append((batch_no, ",".join(batch), str(exc)))
            failed_accessions.update(batch)
        if batch_no % 10 == 0 or batch_no == len(batches):
            print(f"Fetched metadata batches: {batch_no}/{len(batches)}", file=sys.stderr)

    metadata_rows = []
    feature_rows = []
    missing_rows = []
    for accession in accessions:
        record = records_by_request.get(accession)
        if record is None:
            status = "FETCH_ERROR" if accession in failed_accessions else "NOT_FOUND"
            metadata_rows.append({
                "input_accession": accession, "accession_version": "", "organism": "",
                "definition": "", "taxonomy": "", "family": "Unclassified", "record_status": status,
            })
            missing_rows.append({"input_accession": accession, "record_status": status})
            continue
        metadata_rows.append({"input_accession": accession, **record, "record_status": "OK"})
        for feature in features_by_request.get(accession, []):
            feature_rows.append({"input_accession": accession, **feature})

    write_tsv(outdir / "zover_nuccore_metadata.tsv", METADATA_FIELDS, metadata_rows)
    write_tsv(outdir / "zover_nuccore_features.tsv", FEATURE_FIELDS, feature_rows)
    write_tsv(outdir / "zover_accessions_not_returned.tsv", ["input_accession", "record_status"], missing_rows)
    write_tsv(outdir / "zover_metadata_fetch_errors.tsv", ["batch", "accessions", "error"],
              [{"batch": n, "accessions": accs, "error": error} for n, accs, error in failed_batches])

    ok_count = sum(1 for row in metadata_rows if row["record_status"] == "OK")
    print(f"Input accessions: {len(accessions)}")
    print(f"NCBI records matched: {ok_count}")
    print(f"NCBI gene/CDS features: {len(feature_rows)}")
    print(f"Not returned: {len(missing_rows)}")
    print(f"Metadata output: {outdir / 'zover_nuccore_metadata.tsv'}")
    print(f"Feature output: {outdir / 'zover_nuccore_features.tsv'}")
    if failed_batches:
        print(f"WARNING: {len(failed_batches)} batches failed; see zover_metadata_fetch_errors.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()

