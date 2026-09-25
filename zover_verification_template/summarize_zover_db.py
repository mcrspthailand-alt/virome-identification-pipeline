#!/usr/bin/env python3
"""Summarize families and annotated genes/CDS in a ZOVER reference set.

Counts describe records and annotations present in the downloaded reference DB;
they are not a curated or experimentally validated hallmark-gene catalogue.
"""

import argparse
import csv
import html
import re
from collections import Counter, defaultdict
from pathlib import Path


MARKERS = [
    ("RdRp / RNA polymerase", r"rdrp|rna[- ]dependent rna polymerase|rna[- ]directed rna polymerase|rna polymerase"),
    ("ORF1ab / replicase", r"orf1ab|replicase|nonstructural polyprotein|nsp\s*[0-9]+"),
    ("DNA polymerase", r"dna polymerase|\bdpol\b"),
    ("Retroviral Pol / RT / integrase", r"\bpol\b|reverse transcriptase|integrase|retroviral protease"),
    ("Gag", r"\bgag\b"),
    ("Env / glycoprotein / spike", r"envelope|glycoprotein|spike|\bs protein\b"),
    ("Capsid / coat / hexon", r"capsid|coat protein|hexon"),
    ("Nucleocapsid / N protein", r"nucleocapsid|nucleoprotein|\bn protein\b"),
    ("Terminase / packaging", r"terminase|packaging protein|portal protein"),
    ("Helicase", r"helicase"),
    ("VP proteins", r"\bvp\s*[0-9]+\b|viral protein [0-9]+"),
    ("NSP / nonstructural", r"\bnsp\s*[0-9]+\b|nonstructural|non-structural"),
    ("Polyprotein", r"polyprotein"),
    ("ORF protein", r"\borf[0-9a-z-]*\b|open reading frame"),
    ("Matrix / M protein", r"\bmatrix\b|\bm protein\b|membrane protein"),
    ("Phosphoprotein / P protein", r"phosphoprotein|\bp protein\b"),
    ("Fusion / F protein", r"\bfusion\b|\bf protein\b"),
    ("Large / L protein", r"\bl protein\b|large protein"),
    ("Replication-associated / Rep", r"replication-associated protein|\brep\b"),
]


def read_tsv(path):
    with Path(path).open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_delimited(path, fields, rows, delimiter):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter,
                                extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def unique_join(values, limit=30):
    vals = sorted({" ".join((v or "").split()) for v in values if (v or "").strip()})
    if len(vals) > limit:
        return "; ".join(vals[:limit]) + f"; ... (+{len(vals) - limit} more)"
    return "; ".join(vals)


def marker_class(text):
    text = (text or "").lower()
    for label, pattern in MARKERS:
        if re.search(pattern, text, flags=re.I):
            return label
    return "Other / unparsed"


def family_of(row):
    family = (row.get("family") or "").strip()
    if family:
        return family
    for item in (row.get("taxonomy") or "").split(";"):
        if item.strip().lower().endswith("viridae"):
            return item.strip()
    return "Unclassified"


def svg_text(x, y, text, size=12, fill="#20262e", anchor="start", weight="normal", rotate=None):
    transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
    return (f'<text x="{x}" y="{y}" font-family="Arial, sans-serif" font-size="{size}px" '
            f'font-weight="{weight}" fill="{fill}" text-anchor="{anchor}"{transform}>'
            f'{html.escape(str(text))}</text>')


def write_family_bar(path, family_rows):
    data = sorted(family_rows, key=lambda r: (-int(r["nuccore_records"]), r["family"]))[:50]
    data.reverse()
    width, row_h, left, right, top, bottom = 1300, 25, 300, 90, 90, 60
    height = top + bottom + row_h * max(1, len(data))
    plot_w = width - left - right
    max_v = max([int(r["nuccore_records"]) for r in data] or [1])
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             svg_text(36, 38, "ZOVER reference records by virus family", 22, weight="bold"),
             svg_text(36, 63, "Top 50 families by number of nucleotide reference records", 12, fill="#59636e")]
    ticks = min(max_v, 5)
    for k in range(ticks + 1):
        value = round(max_v * k / max(1, ticks))
        x = left + plot_w * value / max_v
        parts.append(f'<line x1="{x:.1f}" y1="{top-4}" x2="{x:.1f}" y2="{height-bottom}" stroke="#e5e7eb"/>')
        parts.append(svg_text(x, height - 24, value, 10, fill="#59636e", anchor="middle"))
    for i, row in enumerate(data):
        y = top + i * row_h
        value = int(row["nuccore_records"])
        bar_w = plot_w * value / max_v
        parts.append(svg_text(left - 12, y + 16, row["family"], 11, anchor="end"))
        parts.append(f'<rect x="{left}" y="{y+3}" width="{bar_w:.1f}" height="17" fill="#3676a8"/>')
        parts.append(svg_text(left + bar_w + 6, y + 16, value, 10, fill="#35404a"))
    parts.append(svg_text(left + plot_w / 2, height - 5, "Nucleotide records", 11, fill="#59636e", anchor="middle"))
    parts.append("</svg>")
    Path(path).write_text("\n".join(parts), encoding="utf-8")


def write_marker_bar(path, marker_rows):
    data = sorted(marker_rows, key=lambda r: (-int(r["annotated_features"]), r["marker_class"]))
    data.reverse()
    width, row_h, left, right, top, bottom = 1300, 30, 360, 90, 90, 60
    height = top + bottom + row_h * max(1, len(data))
    plot_w = width - left - right
    max_v = max([int(r["annotated_features"]) for r in data] or [1])
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             svg_text(36, 38, "Annotated gene/CDS marker classes in ZOVER", 22, weight="bold"),
             svg_text(36, 63, "Classification is keyword-derived from NCBI gene, product, note and record definition", 12, fill="#59636e")]
    ticks = min(max_v, 5)
    for k in range(ticks + 1):
        value = round(max_v * k / max(1, ticks))
        x = left + plot_w * value / max_v
        parts.append(f'<line x1="{x:.1f}" y1="{top-4}" x2="{x:.1f}" y2="{height-bottom}" stroke="#e5e7eb"/>')
        parts.append(svg_text(x, height - 24, value, 10, fill="#59636e", anchor="middle"))
    for i, row in enumerate(data):
        y = top + i * row_h
        value = int(row["annotated_features"])
        bar_w = plot_w * value / max_v
        parts.append(svg_text(left - 12, y + 18, row["marker_class"], 11, anchor="end"))
        parts.append(f'<rect x="{left}" y="{y+4}" width="{bar_w:.1f}" height="19" fill="#218777"/>')
        parts.append(svg_text(left + bar_w + 6, y + 18, value, 10))
    parts.append(svg_text(left + plot_w / 2, height - 5, "Annotated CDS/gene features", 11, fill="#59636e", anchor="middle"))
    parts.append("</svg>")
    Path(path).write_text("\n".join(parts), encoding="utf-8")


def write_heatmap(path, family_marker_rows):
    totals = Counter()
    cells = Counter()
    marker_totals = Counter()
    for row in family_marker_rows:
        family, marker = row["family"], row["marker_class"]
        value = int(row["annotated_features"])
        totals[family] += value
        cells[(family, marker)] += value
        marker_totals[marker] += value
    families = [f for f, _ in totals.most_common(40)]
    markers = [m for m, _ in marker_totals.most_common()]
    if not markers:
        markers = ["No annotated features"]
    width = 220 + 150 * len(markers)
    row_h, top, left, bottom = 27, 125, 310, 75
    height = top + bottom + row_h * max(1, len(families))
    max_v = max(cells.values() or [1])
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             svg_text(28, 36, "Virus family by annotated marker class", 22, weight="bold"),
             svg_text(28, 61, "Cell values are annotated CDS/gene features; showing up to 40 families", 12, fill="#59636e")]
    cell_w = 150
    for j, marker in enumerate(markers):
        x = left + j * cell_w + cell_w / 2
        parts.append(svg_text(x, top - 12, marker, 10, anchor="end", rotate=-38))
    for i, family in enumerate(families):
        y = top + i * row_h
        parts.append(svg_text(left - 10, y + 18, family, 11, anchor="end"))
        for j, marker in enumerate(markers):
            value = cells[(family, marker)]
            t = (value / max_v) ** 0.5 if value else 0
            low = (242, 248, 247)
            high = (27, 112, 100)
            rgb = tuple(round(low[k] * (1 - t) + high[k] * t) for k in range(3))
            color = "#%02x%02x%02x" % rgb
            x = left + j * cell_w
            parts.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{row_h}" fill="{color}" stroke="#ffffff"/>')
            if value:
                fg = "white" if t > 0.57 else "#20262e"
                parts.append(svg_text(x + cell_w / 2, y + 18, value, 10, fill=fg, anchor="middle", weight="bold"))
    parts.append("</svg>")
    Path(path).write_text("\n".join(parts), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, help="zover_nuccore_metadata.tsv from the NCBI fetch script")
    parser.add_argument("--features", required=True, help="zover_nuccore_features.tsv from the NCBI fetch script")
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    metadata_rows = read_tsv(args.metadata)
    feature_rows = read_tsv(args.features)
    metadata_by_input = {r.get("input_accession", ""): r for r in metadata_rows}
    metadata_by_acc = {}
    for row in metadata_rows:
        for key in (row.get("input_accession"), row.get("accession_version")):
            if key:
                metadata_by_acc[key.upper().split(".")[0]] = row

    enriched_features = []
    for feature in feature_rows:
        meta = metadata_by_input.get(feature.get("input_accession", ""), {})
        family = feature.get("family") or family_of(meta)
        definition = meta.get("definition", "")
        text = " ".join([feature.get("gene", ""), feature.get("product", ""),
                         feature.get("note", ""), definition])
        marker = marker_class(text)
        specific = feature.get("gene") or feature.get("product") or definition or "Unannotated feature"
        enriched_features.append({
            "family": family,
            "nuccore_accession": feature.get("nuccore_accession", "") or meta.get("accession_version", ""),
            "input_accession": feature.get("input_accession", ""),
            "organism": feature.get("organism", "") or meta.get("organism", ""),
            "feature_key": feature.get("feature_key", ""),
            "gene": feature.get("gene", ""),
            "product": feature.get("product", ""),
            "protein_id": feature.get("protein_id", ""),
            "locus_tag": feature.get("locus_tag", ""),
            "location": feature.get("location", ""),
            "note": feature.get("note", ""),
            "marker_class": marker,
            "specific_gene_product": " ".join(specific.split()),
            "annotation_source": "NCBI GenBank gene/CDS feature",
        })

    # If a nucleotide record has no gene/CDS rows, preserve its definition as an explicitly weaker fallback.
    feature_accs = {r["input_accession"].upper().split(".")[0] for r in enriched_features if r["input_accession"]}
    for meta in metadata_rows:
        input_acc = meta.get("input_accession", "")
        if not input_acc or input_acc.upper().split(".")[0] in feature_accs:
            continue
        definition = meta.get("definition", "")
        enriched_features.append({
            "family": family_of(meta), "nuccore_accession": meta.get("accession_version", ""),
            "input_accession": input_acc, "organism": meta.get("organism", ""),
            "feature_key": "record_definition", "gene": "", "product": definition,
            "protein_id": "", "locus_tag": "", "location": "", "note": "",
            "marker_class": marker_class(definition), "specific_gene_product": definition or "Unannotated",
            "annotation_source": "NCBI record definition only; no gene/CDS feature returned",
        })

    family_records = defaultdict(list)
    for row in metadata_rows:
        family_records[family_of(row)].append(row)
    feature_groups = defaultdict(list)
    for row in enriched_features:
        feature_groups[(row["family"], row["marker_class"])].append(row)

    family_summary = []
    for family, records in sorted(family_records.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        feats = [r for r in enriched_features if r["family"] == family]
        family_summary.append({
            "family": family,
            "nuccore_records": len(records),
            "unique_organisms": len({r.get("organism", "") for r in records if r.get("organism", "")}),
            "annotated_features": sum(1 for r in feats if r["feature_key"] != "record_definition"),
            "marker_classes": unique_join([r["marker_class"] for r in feats]),
            "gene_product_labels": unique_join([r["specific_gene_product"] for r in feats], 60),
        })

    family_marker = []
    for (family, marker), rows in sorted(feature_groups.items()):
        family_marker.append({
            "family": family, "marker_class": marker,
            "annotated_features": len(rows),
            "unique_nuccore_records": len({r["input_accession"] for r in rows if r["input_accession"]}),
            "unique_organisms": len({r["organism"] for r in rows if r["organism"]}),
            "genes_products": unique_join([r["specific_gene_product"] for r in rows], 30),
            "example_accessions": unique_join([r["nuccore_accession"] for r in rows], 20),
        })

    marker_totals = defaultdict(list)
    for row in enriched_features:
        marker_totals[row["marker_class"]].append(row)
    marker_summary = []
    for marker, rows in sorted(marker_totals.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        marker_summary.append({
            "marker_class": marker,
            "annotated_features": len(rows),
            "unique_nuccore_records": len({r["input_accession"] for r in rows if r["input_accession"]}),
            "unique_families": len({r["family"] for r in rows if r["family"] != "Unclassified"}),
            "genes_products": unique_join([r["specific_gene_product"] for r in rows], 50),
        })

    outdir = Path(args.outdir)
    table_dir = outdir / "tables"
    figure_dir = outdir / "figures"
    table_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)
    overview = [
        {"metric": "Input nucleotide reference records", "value": len(metadata_rows)},
        {"metric": "NCBI records successfully fetched", "value": sum(r.get("record_status") == "OK" for r in metadata_rows)},
        {"metric": "Virus families represented (excluding Unclassified)", "value": len({family_of(r) for r in metadata_rows if family_of(r) != "Unclassified"})},
        {"metric": "NCBI gene/CDS features", "value": sum(r["feature_key"] != "record_definition" for r in enriched_features)},
        {"metric": "Reference records using definition-only fallback", "value": sum(r["feature_key"] == "record_definition" for r in enriched_features)},
        {"metric": "Marker classes (including Other/unparsed)", "value": len(marker_summary)},
        {"metric": "Interpretation", "value": "Database annotation inventory; not a curated hallmark-gene catalogue"},
    ]

    tables = [
        ("zover_db_overview", ["metric", "value"], overview),
        ("zover_db_family_summary", list(family_summary[0]) if family_summary else ["family", "nuccore_records", "unique_organisms", "annotated_features", "marker_classes", "gene_product_labels"], family_summary),
        ("zover_db_family_marker_long", list(family_marker[0]) if family_marker else ["family", "marker_class", "annotated_features", "unique_nuccore_records", "unique_organisms", "genes_products", "example_accessions"], family_marker),
        ("zover_db_marker_summary", list(marker_summary[0]) if marker_summary else ["marker_class", "annotated_features", "unique_nuccore_records", "unique_families", "genes_products"], marker_summary),
        ("zover_db_gene_cds_annotations", list(enriched_features[0]) if enriched_features else ["family", "nuccore_accession", "input_accession", "organism", "feature_key", "gene", "product", "protein_id", "locus_tag", "location", "note", "marker_class", "specific_gene_product", "annotation_source"], enriched_features),
    ]
    for name, fields, rows in tables:
        write_delimited(table_dir / f"{name}.tsv", fields, rows, "\t")
        write_delimited(table_dir / f"{name}.csv", fields, rows, ",")

    write_family_bar(figure_dir / "figure_01_zover_family_records.svg", family_summary)
    write_marker_bar(figure_dir / "figure_02_zover_marker_classes.svg", marker_summary)
    write_heatmap(figure_dir / "figure_03_zover_family_marker_heatmap.svg", family_marker)
    print(f"Virus families represented: {overview[2]['value']}")
    print(f"NCBI gene/CDS features: {overview[3]['value']}")
    print(f"Tables: {table_dir}")
    print(f"Figures (SVG): {figure_dir}")


if __name__ == "__main__":
    main()

