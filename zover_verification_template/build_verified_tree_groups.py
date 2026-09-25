#!/usr/bin/env python3
"""Annotate completed ZOVER hits and prepare family/gene tree groups.

This script does not run BLAST. It consumes the completed combined tree-candidate
table, extracted query FASTAs, the ZOVER protein FASTA and NCBI annotation tables.
"""

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path


ACC_RE = re.compile(r"\b[A-Z]{1,5}_?\d{4,}(?:\.\d+)?\b", re.I)
DB_ACC_RE = re.compile(r"(?:gb|ref|emb|dbj|tpg|tpe|tpd)\|([^|\s]+)\|", re.I)
ORF_RE = re.compile(r"\b(orf\s*[0-9]+[a-z]?)\b", re.I)
NSP_RE = re.compile(r"\b(nsp\s*[0-9]+)\b", re.I)
VP_RE = re.compile(r"\b(vp\s*[0-9]+)\b", re.I)


def read_tsv(path):
    with Path(path).open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: no TSV header in {path}")
        return reader.fieldnames, list(reader)


def write_tsv(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "r", encoding="utf-8", errors="replace")


def fasta_records(path):
    name, chunks = None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            elif name is not None:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)


def wrap(seq, width=80):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def clean(value):
    return " ".join((value or "").replace("\t", " ").split())


def norm_acc(value):
    return clean(value).upper()


def base_acc(value):
    return re.sub(r"\.\d+$", "", norm_acc(value))


def accession_candidates(text):
    text = text or ""
    found = []
    for match in DB_ACC_RE.finditer(text):
        found.append(match.group(1).strip().rstrip(".,;"))
    found.extend(ACC_RE.findall(text))
    result, seen = [], set()
    for value in found:
        if not ACC_RE.fullmatch(value):
            continue
        key = norm_acc(value)
        if key not in seen:
            result.append(value)
            seen.add(key)
    return result


def id_token(header):
    return header.split()[0] if header.split() else ""


def slug(value):
    value = clean(value)
    value = re.sub(r"[^A-Za-z0-9._+-]+", "_", value).strip("._-")
    return value[:120] or "unknown"


def first_nonempty(row, names):
    for name in names:
        value = clean(row.get(name, ""))
        if value:
            return value
    return ""


def family_from_meta(meta):
    family = clean(meta.get("family", ""))
    if family and family.lower() != "unclassified":
        return family
    for item in (meta.get("taxonomy", "") or "").split(";"):
        if item.strip().lower().endswith("viridae"):
            return item.strip()
    return ""


def canonical_gene(gene, product, blast_gene, title):
    text = " ".join(x for x in [gene, product, blast_gene, title] if clean(x))
    low = text.lower()
    match = ORF_RE.search(text)
    if match:
        return re.sub(r"\s+", "", match.group(1)).upper()
    match = NSP_RE.search(text)
    if match:
        return re.sub(r"\s+", "", match.group(1)).upper()
    match = VP_RE.search(text)
    if match:
        return re.sub(r"\s+", "", match.group(1)).upper()
    exact_gene = clean(gene).lower()
    if exact_gene in {"pol", "polymerase", "rna-dependent rna polymerase", "rdrp"}:
        return "RdRp/Pol"
    if re.search(r"rna[- ]dependent rna polymerase|rna[- ]directed rna polymerase|\brdrp\b", low):
        return "RdRp"
    if re.search(r"dna polymerase|\bdpol\b", low):
        return "DNA_polymerase"
    if re.search(r"reverse transcriptase|integrase|retroviral protease|retroviral pol|\bpol\b", low):
        return "Pol"
    if re.search(r"\bgag\b", low):
        return "Gag"
    if re.search(r"\benv\b|envelope glycoprotein|env[- ]protein", low):
        return "Env"
    if re.search(r"spike|\bs protein\b", low):
        return "Spike"
    if re.search(r"nucleocapsid|nucleoprotein|\bn protein\b", low):
        return "Nucleocapsid"
    if re.search(r"capsid|coat protein|hexon", low):
        return "Capsid"
    if re.search(r"terminase|packaging protein|portal protein", low):
        return "Terminase"
    if re.search(r"helicase", low):
        return "Helicase"
    if re.search(r"polyprotein", low):
        return "Polyprotein"
    chosen = clean(gene) or clean(product) or clean(blast_gene) or clean(title)
    chosen = re.sub(r"\s+", " ", chosen).strip(" ;|.-")
    return chosen[:100] or "Unannotated"


def parse_protein_header(header):
    accs = accession_candidates(header)
    protein_acc = ""
    match = DB_ACC_RE.search(header)
    if match:
        protein_acc = match.group(1).strip()
    if not protein_acc and accs:
        protein_acc = accs[0]
    nc_candidates = []
    for paren in re.findall(r"\(([^()]*)\)", header):
        candidate = paren.split(",", 1)[0].strip()
        if ACC_RE.fullmatch(candidate):
            nc_candidates.append(candidate)
    for m in re.finditer(r"(?:ncbi_nuccore|nuccore_accession)=([^|\s]+)", header, re.I):
        if ACC_RE.fullmatch(m.group(1)):
            nc_candidates.append(m.group(1))
    return protein_acc, accs, nc_candidates


def load_annotation(metadata_path, features_path):
    _mfields, metadata = read_tsv(metadata_path)
    _ffields, features = read_tsv(features_path)
    meta_by_acc, feature_by_protein, features_by_nuccore = {}, defaultdict(list), defaultdict(list)
    for row in metadata:
        for value in (row.get("input_accession", ""), row.get("accession_version", "")):
            if value:
                meta_by_acc[base_acc(value)] = row
    for row in features:
        pid = row.get("protein_id", "")
        if pid:
            feature_by_protein[base_acc(pid)].append(row)
        for value in (row.get("input_accession", ""), row.get("nuccore_accession", "")):
            if value:
                features_by_nuccore[base_acc(value)].append(row)
    return metadata, features, meta_by_acc, feature_by_protein, features_by_nuccore


def select_feature(accessions, feature_by_protein, features_by_nuccore, gene_hint=""):
    for accession in accessions:
        rows = feature_by_protein.get(base_acc(accession), [])
        cds = [r for r in rows if r.get("feature_key") == "CDS"]
        if cds:
            return cds[0], "NCBI CDS matched by protein_id"
    possible = []
    for accession in accessions:
        possible.extend(r for r in features_by_nuccore.get(base_acc(accession), []) if r.get("feature_key") == "CDS")
    if not possible:
        return {}, ""
    hint_words = set(re.findall(r"[a-z0-9]+", (gene_hint or "").lower()))
    def score(row):
        annotation = " ".join([row.get("gene", ""), row.get("product", ""), row.get("note", "")]).lower()
        words = set(re.findall(r"[a-z0-9]+", annotation))
        return len(hint_words & words)
    possible.sort(key=score, reverse=True)
    return possible[0], "NCBI CDS matched by nuccore accession; closest title terms"


def original_family(row):
    for key, value in row.items():
        if key.startswith("original_") and any(term in key.lower() for term in ("family", "taxon")) and clean(value):
            return clean(value)
    return ""


def original_organism(row):
    for key, value in row.items():
        if key.startswith("original_") and any(term in key.lower() for term in ("virus", "organism", "species")) and clean(value):
            return clean(value)
    return ""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", required=True, help="Combined tree_candidates.with_original.tsv")
    parser.add_argument("--queries", required=True, help="Directory containing <sample>_merged_contigs.fa")
    parser.add_argument("--protein-fasta", required=True, help="ZOVER-associated protein FASTA, plain or .gz")
    parser.add_argument("--metadata", required=True, help="zover_nuccore_metadata.tsv")
    parser.add_argument("--features", required=True, help="zover_nuccore_features.tsv")
    parser.add_argument("--outdir", required=True, help="Postprocess output directory")
    args = parser.parse_args()

    candidate_fields, candidates = read_tsv(args.candidates)
    _metadata, _features, meta_by_acc, feature_by_protein, features_by_nuccore = load_annotation(args.metadata, args.features)
    query_dir, outdir = Path(args.queries), Path(args.outdir)
    groups_dir = outdir / "tree_groups"
    groups_dir.mkdir(parents=True, exist_ok=True)

    query_sequences = {}
    samples = sorted({clean(r.get("sample", "")) for r in candidates if clean(r.get("sample", ""))})
    for sample in samples:
        fasta = query_dir / f"{sample}_merged_contigs.fa"
        if not fasta.is_file():
            continue
        query_sequences[sample] = {id_token(h): (h, seq) for h, seq in fasta_records(fasta)}

    annotated = []
    for row in candidates:
        bx_hit = row.get("blastx_hit", "")
        bx_accessions = row.get("blastx_accessions", "")
        bx_title = row.get("blastx_title", "")
        ref_title = row.get("blastx_ref_title_by_accession", "")
        accessions = accession_candidates(" ".join([bx_hit, bx_accessions, bx_title, ref_title]))
        feature, feature_source = select_feature(accessions, feature_by_protein, features_by_nuccore,
                                                 " ".join([row.get("blastx_gene_orf", ""), bx_title, ref_title]))
        nuc_acc = feature.get("nuccore_accession", "") or next((a for a in accessions if base_acc(a) in meta_by_acc), "")
        meta = meta_by_acc.get(base_acc(nuc_acc), {})
        family = clean(feature.get("family", ""))
        if family.lower() == "unclassified":
            family = ""
        family = family or family_from_meta(meta)
        fam_source = "NCBI GenBank taxonomy" if family else ""
        organism = clean(feature.get("organism", "")) or clean(meta.get("organism", ""))
        if not family:
            family = original_family(row)
            if family:
                fam_source = "original candidate annotation"
        if not organism:
            organism = original_organism(row)
        if not family:
            # Keep taxa without family-level NCBI annotation separate rather than pooling unrelated viruses.
            fallback = organism or (nuc_acc or (accessions[0] if accessions else bx_hit) or "unknown reference")
            family = "Unclassified_" + fallback
            fam_source = "fallback taxon/reference label; verify manually"

        gene = clean(feature.get("gene", ""))
        product = clean(feature.get("product", ""))
        blast_gene = clean(row.get("blastx_gene_orf", ""))
        group_gene = canonical_gene(gene, product, blast_gene, " ".join([bx_title, ref_title]))
        group_id = slug(family) + "__" + slug(group_gene)
        annotation_source = feature_source or ("NCBI record metadata only" if meta else "BLAST title fallback")

        item = dict(row)
        item.update({
            "virus_family": family,
            "family_annotation_source": fam_source,
            "virus_organism": organism,
            "reference_protein_accession": next((a for a in accessions if base_acc(a) in feature_by_protein), accessions[0] if accessions else ""),
            "reference_nuccore_accession": nuc_acc,
            "ncbi_gene": gene,
            "ncbi_cds_product": product,
            "ncbi_cds_note": clean(feature.get("note", "")),
            "ncbi_cds_location": clean(feature.get("location", "")),
            "gene_orf": group_gene,
            "annotation_source": annotation_source,
            "tree_group_id": group_id,
            "query_sequence_found": "YES" if clean(row.get("contig", "")) in query_sequences.get(clean(row.get("sample", "")), {}) else "NO",
            "tree_note": "Candidate only; review taxon, homologous region, alignment quality and biological plausibility before tree inclusion",
        })
        annotated.append(item)

    # Stream the protein reference and retain only records matching a candidate group.
    wanted_groups = {(r["virus_family"], r["gene_orf"]) for r in annotated}
    refs_by_group = defaultdict(dict)
    hit_refs_by_group = defaultdict(dict)
    hit_to_group = {}
    for row in annotated:
        group = (row["virus_family"], row["gene_orf"])
        hit_accessions = accession_candidates(row.get("blastx_hit", ""))
        if hit_accessions:
            hit_to_group.setdefault(base_acc(hit_accessions[0]), set()).add(group)

    for header, seq in fasta_records(args.protein_fasta):
        protein_acc, all_accs, nuccore_accs = parse_protein_header(header)
        feature = None
        if protein_acc:
            rows = feature_by_protein.get(base_acc(protein_acc), [])
            feature = next((r for r in rows if r.get("feature_key") == "CDS"), None)
        if feature is None:
            for acc in all_accs:
                rows = feature_by_protein.get(base_acc(acc), [])
                feature = next((r for r in rows if r.get("feature_key") == "CDS"), None)
                if feature:
                    break
        if feature:
            family = clean(feature.get("family", ""))
            if family.lower() == "unclassified":
                family = ""
            family = family or family_from_meta(meta_by_acc.get(base_acc(feature.get("nuccore_accession", "")), {}))
            gene_group = canonical_gene(feature.get("gene", ""), feature.get("product", ""), "", header)
            group = (family, gene_group)
            if group in wanted_groups:
                key = protein_acc or id_token(header)
                refs_by_group[group].setdefault(base_acc(key), (key, header, seq, feature))
        subject_accession = protein_acc or (all_accs[0] if all_accs else "")
        if subject_accession:
            for group in hit_to_group.get(base_acc(subject_accession), set()):
                hit_refs_by_group[group].setdefault(
                    base_acc(subject_accession), (subject_accession, header, seq, feature or {})
                )

    all_annotated_fields = list(dict.fromkeys(candidate_fields + [
        "virus_family", "family_annotation_source", "virus_organism", "reference_protein_accession",
        "reference_nuccore_accession", "ncbi_gene", "ncbi_cds_product", "ncbi_cds_note",
        "ncbi_cds_location", "gene_orf", "annotation_source", "tree_group_id",
        "query_sequence_found", "tree_note",
    ]))
    annotated_path = outdir / "verified_contigs_gene_CDS_annotation.tsv"
    write_tsv(annotated_path, all_annotated_fields, annotated)

    profile_groups = defaultdict(list)
    for row in annotated:
        key = (row.get("sample", ""), row.get("profile_name", ""), row["virus_family"], row["gene_orf"])
        profile_groups[key].append(row)
    profile_rows = []
    for (sample, profile_name, family, gene_group), members in sorted(profile_groups.items()):
        profile_rows.append({
            "sample": sample, "profile_name": profile_name, "virus_family": family, "gene_orf": gene_group,
            "verified_contig_count": len({r.get("contig", "") for r in members}),
            "high_confidence_contigs": sum(r.get("confidence") == "HIGH" for r in members),
            "medium_confidence_contigs": sum(r.get("confidence") == "MEDIUM" for r in members),
            "contigs": ";".join(sorted({r.get("contig", "") for r in members if r.get("contig", "")})),
            "ncbi_gene": ";".join(sorted({r.get("ncbi_gene", "") for r in members if r.get("ncbi_gene", "")})),
            "ncbi_cds_product": ";".join(sorted({r.get("ncbi_cds_product", "") for r in members if r.get("ncbi_cds_product", "")})),
            "reference_accessions": ";".join(sorted({r.get("blastx_accessions", "") for r in members if r.get("blastx_accessions", "")})),
            "tree_group_id": members[0]["tree_group_id"],
        })
    profile_fields = ["sample", "profile_name", "virus_family", "gene_orf", "verified_contig_count",
                      "high_confidence_contigs", "medium_confidence_contigs", "contigs", "ncbi_gene",
                      "ncbi_cds_product", "reference_accessions", "tree_group_id"]
    write_tsv(outdir / "verified_profile_gene_CDS_annotation.tsv", profile_fields, profile_rows)

    group_members = defaultdict(list)
    for row in annotated:
        group_members[row["tree_group_id"]].append(row)
    group_index, curation_rows = [], []
    for group_id, members in sorted(group_members.items()):
        family, gene_group = members[0]["virus_family"], members[0]["gene_orf"]
        folder = groups_dir / group_id
        folder.mkdir(parents=True, exist_ok=True)
        metadata_fields = list(dict.fromkeys(all_annotated_fields + ["include_in_tree", "curation_note"]))
        for row in members:
            row["include_in_tree"] = ""
            row["curation_note"] = ""
            curation_rows.append({
                "tree_group_id": group_id, "sample": row.get("sample", ""), "contig": row.get("contig", ""),
                "virus_family": family, "gene_orf": gene_group, "virus_organism": row.get("virus_organism", ""),
                "ncbi_gene": row.get("ncbi_gene", ""), "ncbi_cds_product": row.get("ncbi_cds_product", ""),
                "confidence": row.get("confidence", ""), "blastx_pident": row.get("blastx_pident", ""),
                "blastx_qcovs": row.get("blastx_qcovs", ""), "blastx_evalue": row.get("blastx_evalue", ""),
                "blastx_bitscore": row.get("blastx_bitscore", ""), "include_in_tree": "", "curation_note": "",
            })
        write_tsv(folder / "metadata.tsv", metadata_fields, members)

        seen_contigs = set()
        with (folder / "query_contigs.fna").open("w", encoding="utf-8") as out:
            for row in members:
                sample, contig = clean(row.get("sample", "")), clean(row.get("contig", ""))
                rec = query_sequences.get(sample, {}).get(contig)
                key = (sample, contig)
                if rec and key not in seen_contigs:
                    header, seq = rec
                    out.write(f">{slug(sample)}|{slug(contig)}\n{wrap(seq)}\n")
                    seen_contigs.add(key)

        refs = dict(refs_by_group.get((family, gene_group), {}))
        refs.update({k: v for k, v in hit_refs_by_group.get((family, gene_group), {}).items() if k not in refs})
        with (folder / "reference_proteins.faa").open("w", encoding="utf-8") as out:
            for key, header, seq, _feature in sorted(refs.values(), key=lambda x: x[0]):
                out.write(f">{slug(key)} {clean(header)}\n{wrap(seq)}\n")
        ref_annotation_rows = []
        for key, header, _seq, feature in sorted(refs.values(), key=lambda x: x[0]):
            ref_annotation_rows.append({
                "protein_accession": key, "reference_header": clean(header),
                "nuccore_accession": feature.get("nuccore_accession", ""),
                "family": feature.get("family", ""), "organism": feature.get("organism", ""),
                "gene": feature.get("gene", ""), "product": feature.get("product", ""),
                "protein_id": feature.get("protein_id", ""), "note": feature.get("note", ""),
            })
        write_tsv(folder / "reference_annotations.tsv",
                  ["protein_accession", "reference_header", "nuccore_accession", "family", "organism", "gene", "product", "protein_id", "note"],
                  ref_annotation_rows)
        (folder / "README.md").write_text(
            "# Candidate tree group\n\n"
            f"Family: `{family}`\n\nGene/product group: `{gene_group}`\n\n"
            f"Candidate contig rows: {len(members)}\n\n"
            f"Reference protein sequences included: {len(refs)}\n\n"
            "These are candidate inputs, not a finalized tree. Curate `manual_curation_template.tsv` in the parent `tree_groups` directory. "
            "Do not mix different virus families or non-homologous genes/ORFs. The query FASTA contains nucleotide contigs; the tree runner translates the BLASTX-aligned query region.\n",
            encoding="utf-8")
        group_index.append({
            "tree_group_id": group_id, "virus_family": family, "gene_orf": gene_group,
            "candidate_contig_rows": len(members), "unique_contigs": len({r.get("contig", "") for r in members}),
            "reference_protein_sequences": len(refs), "query_sequences_written": len(seen_contigs),
            "directory": str(folder), "family_annotation_sources": ";".join(sorted({r["family_annotation_source"] for r in members})),
            "review_before_tree": "YES",
        })

    curation_fields = ["tree_group_id", "sample", "contig", "virus_family", "gene_orf", "virus_organism",
                       "ncbi_gene", "ncbi_cds_product", "confidence", "blastx_pident", "blastx_qcovs",
                       "blastx_evalue", "blastx_bitscore", "include_in_tree", "curation_note"]
    write_tsv(groups_dir / "manual_curation_template.tsv", curation_fields, curation_rows)
    write_tsv(groups_dir / "group_index.tsv",
              ["tree_group_id", "virus_family", "gene_orf", "candidate_contig_rows", "unique_contigs",
               "reference_protein_sequences", "query_sequences_written", "directory", "family_annotation_sources", "review_before_tree"],
              group_index)

    qc_rows = []
    for sample in samples:
        sample_rows = [r for r in annotated if r.get("sample") == sample]
        missing = [r.get("contig", "") for r in sample_rows if r.get("query_sequence_found") != "YES"]
        qc_rows.append({"sample": sample, "verified_candidate_rows": len(sample_rows),
                        "candidate_rows_missing_query_sequence": len(missing),
                        "missing_contigs": ";".join(missing)})
    write_tsv(outdir / "postprocess_query_sequence_qc.tsv",
              ["sample", "verified_candidate_rows", "candidate_rows_missing_query_sequence", "missing_contigs"], qc_rows)
    print(f"Annotated verified rows: {len(annotated)}")
    print(f"Profile annotation rows: {len(profile_rows)}")
    print(f"Family+gene candidate groups: {len(group_index)}")
    print(f"Annotated contigs: {annotated_path}")
    print(f"Tree groups and manual curation template: {groups_dir}")


if __name__ == "__main__":
    main()

