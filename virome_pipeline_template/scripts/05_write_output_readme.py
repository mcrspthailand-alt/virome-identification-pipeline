#!/usr/bin/env python3
import argparse
import shlex
import subprocess
from pathlib import Path


FILE_DESCRIPTIONS = {
    "all_contig_read_counts.tsv": "Read mapping counts for every contig in every sample and length bin.",
    "diamond_blastx_best_hits.tsv": "Best DIAMOND blastx hit per contig. Query is nucleotide contig FASTA translated in six frames.",
    "diamond_blastp_best_hits.tsv": "Best DIAMOND blastp hit per contig. Query is Prodigal-predicted proteins.",
    "all_viral_evidence.tsv": "Main pre-taxonomy evidence table combining read counts, blastx, blastp, VirSorter2, geNomad, CheckV, and Prodigal.",
    "all_viral_evidence.parsed_organisms.tsv": "Evidence table with parsed virus_name and virus_name_norm columns from DIAMOND hit titles.",
    "virus_names_with_taxid.tsv": "Parsed virus names mapped to NCBI taxid candidates.",
    "taxid_ranks.tsv": "Taxonomic rank table for mapped taxids.",
    "virus_taxonomy_map.tsv": "Lookup table from virus_name to taxid and lineage.",
    "virus_reads_matrix.tsv": "Virus/organism-by-sample mapped-read matrix.",
    "virus_reads_matrix_with_taxonomy.tsv": "Read matrix with taxonomy metadata columns.",
    "virus_presence_absence_matrix.tsv": "Presence/absence matrix derived from mapped reads.",
    "virus_reads_by_sample_long.tsv": "Long-format table of virus_name, sample, and reads.",
    "taxonomy_assignment_report.txt": "Summary of taxonomy assignment rate.",
    "all_contigs_with_annotations.tsv": "All contigs with evidence and taxonomy annotations, including non-viral contigs.",
    "final_virus_contigs.tsv": "Final putative viral contigs only: putative_viral is TRUE.",
    "final_virus_contigs_bacteriophage.tsv": "Subset of final_virus_contigs.tsv classified as bacteriophage.",
    "final_virus_contigs_eukaryotic.tsv": "Subset of final_virus_contigs.tsv classified as eukaryotic viruses.",
    "final_virus_contigs_archaeal.tsv": "Subset of final_virus_contigs.tsv classified as archaeal viruses.",
    "final_virus_contigs_unidentified.tsv": "Subset of final_virus_contigs.tsv without broad group assignment.",
    "final_virus_summary.txt": "Final row counts and group/confidence summaries.",
}

COLUMN_NOTES = {
    "sample": "Sample ID.",
    "contig": "Contig identifier.",
    "len_bin": "Length bin: len500_999 or ge1000.",
    "contig_length": "Contig length in bp from samtools idxstats.",
    "mapped_reads": "Reads mapped back to this contig.",
    "unmapped_reads": "Unmapped count reported by samtools idxstats for this reference entry.",
    "bam": "Path to the sorted BAM used for mapping/counting.",
    "blastx_hit": "TRUE if contig has a passing DIAMOND blastx hit.",
    "blastp_hit": "TRUE if contig has a passing Prodigal ORF DIAMOND blastp hit.",
    "blastx_sseqid": "Subject sequence ID for best blastx hit.",
    "blastp_sseqid": "Subject sequence ID for best blastp hit.",
    "blastx_evalue": "E-value for best blastx hit.",
    "blastp_evalue": "E-value for best blastp hit.",
    "blastx_bitscore": "Bit score for best blastx hit.",
    "blastp_bitscore": "Bit score for best blastp hit.",
    "blastx_stitle": "Subject title for best blastx hit.",
    "blastp_stitle": "Subject title for best blastp hit.",
    "blastp_n_orfs_strict": "Number of Prodigal ORFs on the contig passing strict blastp criteria.",
    "vs2_hit": "TRUE if VirSorter2 reported this contig.",
    "genomad_hit": "TRUE if geNomad reported this contig.",
    "putative_viral": "TRUE if any of vs2_hit, genomad_hit, blastx_hit, or blastp_hit is TRUE.",
    "viral_evidence_count": "Number of independent evidence sources supporting the contig.",
    "viral_confidence": "High for evidence_count >=2, Medium for evidence_count ==1, None for 0.",
    "viral_evidence": "Semicolon-separated evidence sources.",
    "virus_name": "Parsed virus/organism name from DIAMOND titles.",
    "taxid": "NCBI taxonomy ID assigned from virus_name.",
    "family_or_best": "Family rank if available, otherwise best higher rank.",
    "best_rank": "Rank represented by family_or_best.",
    "lineage_ranks": "Semicolon-separated lineage.",
    "virus_group": "Broad group: Bacteriophage, Eukaryotic viruses, Archaeal viruses, or Unidentified.",
    "prodigal_total_orfs": "Total ORFs predicted by Prodigal on the contig.",
}

PREFIX_NOTES = {
    "vs2_": "VirSorter2 annotation column copied from selected VirSorter2 TSV output.",
    "genomad_": "geNomad annotation column copied from selected geNomad TSV output.",
    "checkv_": "CheckV quality/completeness annotation. CheckV is reported, not used as a filter.",
    "blastx_": "DIAMOND blastx column using nucleotide contigs as query.",
    "blastp_": "DIAMOND blastp column using Prodigal proteins as query.",
}


def load_config(path: str):
    cmd = f"set -a; source {shlex.quote(path)}; env"
    proc = subprocess.run(["bash", "-lc", cmd], check=True, text=True, capture_output=True)
    cfg = {}
    for line in proc.stdout.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            cfg[k] = v
    return cfg


def columns(path: Path):
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        header = handle.readline().rstrip("\n")
    return header.split("\t") if header else []


def describe_col(col):
    if col in COLUMN_NOTES:
        return COLUMN_NOTES[col]
    for prefix, note in PREFIX_NOTES.items():
        if col.startswith(prefix):
            return note
    return "Pipeline output column; see the upstream tool output for exact interpretation."


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    cfg = load_config(args.config)

    out_base = Path(cfg["OUT_BASE"])
    merged = out_base / "merged"
    readme = out_base / "README.md"

    lines = [
        "# Virome Pipeline Outputs",
        "",
        "This README was generated automatically after the pipeline finished.",
        "",
        "## How To Run",
        "",
        "```bash",
        "bash run_virome_pipeline.sh /path/to/config.env",
        "```",
        "",
        "## Evidence Rule",
        "",
        "```text",
        "putative_viral = vs2_hit OR genomad_hit OR blastx_hit OR blastp_hit",
        "viral_confidence = High if evidence_count >= 2, Medium if evidence_count == 1, None if evidence_count == 0",
        "```",
        "",
        "CheckV is included as annotation only and is not used to discard contigs.",
        "",
        "## Output Structure",
        "",
        "```text",
        "viral_id_pipeline/",
        "  logs/",
        "  samples/",
        "  merged/",
        "  taxonomy/",
        "```",
        "",
        "Per-sample outputs are under `samples/<sample>/`. Project-level outputs are under `merged/`.",
        "",
        "## Merged Files",
        "",
    ]

    for filename, desc in FILE_DESCRIPTIONS.items():
        path = merged / filename
        lines.extend([f"### `{filename}`", "", desc, ""])
        cols = columns(path)
        if cols:
            lines.extend(["Columns:", ""])
            for col in cols:
                lines.append(f"- `{col}`: {describe_col(col)}")
            lines.append("")
        else:
            lines.extend(["Columns: file was not present when this README was generated.", ""])

    readme.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(f"[OK] wrote {readme}")


if __name__ == "__main__":
    main()
