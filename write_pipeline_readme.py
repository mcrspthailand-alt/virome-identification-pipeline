#!/usr/bin/env python3
import os
from pathlib import Path


PROJECT_ROOT = Path("/home/panda/workspace/projects/nextflow_metagenomic/eDNA/air")
OUT_BASE = PROJECT_ROOT / "non_hybrid" / "viral_id_pipeline"
MERGED_DIR = OUT_BASE / "merged"
README = OUT_BASE / "README.md"

SAMPLES = ["air13", "air14", "air15", "air17", "air19", "air20", "air21", "air23", "air24"]

FILE_DESCRIPTIONS = {
    "all_contig_read_counts.tsv": "Read mapping counts for every contig in every sample and length bin.",
    "diamond_best_hits_initial.tsv": "Initial per-contig DIAMOND best-hit table generated before the full alignment-field rerun.",
    "diamond_best_hits_with_read_counts.tsv": "Initial DIAMOND best-hit table after adding mapped-read counts.",
    "virus_reads_matrix_initial.tsv": "Initial organism-by-sample read matrix before taxonomy cleaning.",
    "diamond_best_hits_full.tsv": "DIAMOND best-hit table with full alignment fields such as qstart, qend, sstart, send, qcovhsp.",
    "all_viral_evidence.tsv": "Main pre-taxonomy evidence table combining read counts, DIAMOND, VirSorter2, geNomad, CheckV, and Prodigal.",
    "all_viral_evidence.parsed_organisms.tsv": "Evidence table with parsed and normalized organism names from DIAMOND stitle.",
    "virus_taxonomy_map.tsv": "Lookup table mapping parsed virus/organism names to taxid and lineage summary.",
    "taxid.ranks.tsv": "Taxonomy rank table for mapped taxids.",
    "virus_reads_matrix.tsv": "Cleaned virus/organism-by-sample mapped-read matrix.",
    "virus_reads_matrix_with_taxonomy.tsv": "Read matrix with taxid and lineage metadata columns.",
    "virus_presence_absence_matrix.tsv": "Presence/absence matrix derived from the read matrix.",
    "virus_reads_by_sample_long.tsv": "Long-format table of virus/organism, sample, and reads for plotting or downstream summaries.",
    "taxonomy_assignment_report.txt": "Summary of taxonomy assignment counts and mapping rate.",
    "final_virus_contigs.tsv": "Main final table with evidence, taxonomy, lineage, virus group, read counts, CheckV, and Prodigal annotations.",
    "final_virus_contigs_bacteriophage.tsv": "Subset of final_virus_contigs.tsv classified as bacteriophage.",
    "final_virus_contigs_eukaryotic.tsv": "Subset of final_virus_contigs.tsv classified as eukaryotic viruses.",
    "final_virus_contigs_archaeal.tsv": "Subset of final_virus_contigs.tsv classified as archaeal viruses.",
    "final_virus_contigs_unidentified.tsv": "Subset of final_virus_contigs.tsv without a confident group classification.",
    "final_virus_summary.txt": "Final summary report with row counts and virus-group counts.",
}

COLUMN_NOTES = {
    "sample": "Sample ID.",
    "contig": "Contig identifier.",
    "len_bin": "Length bin: len500_999 or ge1000.",
    "contig_length": "Contig length in bp from mapping/index statistics.",
    "length": "Best available contig length column used in final tables.",
    "mapped_reads": "Number of reads mapped back to the contig.",
    "unmapped_reads": "Unmapped-read count reported by samtools idxstats for this reference entry.",
    "bam": "Path to the sorted BAM used for read-count and coverage calculation.",
    "blastx_hit": "TRUE if the contig has a passing DIAMOND/BLASTx hit.",
    "sseqid": "Subject/reference sequence ID from DIAMOND.",
    "pident": "Percent identity of the best DIAMOND alignment.",
    "aln_len": "Alignment length of the best DIAMOND hit.",
    "evalue": "E-value of the best DIAMOND hit.",
    "bitscore": "Bit score of the best DIAMOND hit.",
    "qcovhsp": "Query coverage of the high-scoring pair from DIAMOND.",
    "stitle": "Subject title/description from the DIAMOND database.",
    "n_orfs_strict": "Number of Prodigal ORFs on the contig passing strict DIAMOND criteria.",
    "blastx_organism": "Parsed organism name from DIAMOND stitle.",
    "vs2_hit": "TRUE if VirSorter2 reported this contig in the selected score/summary table.",
    "vs2_source": "VirSorter2 output file used for the evidence row.",
    "genomad_hit": "TRUE if geNomad reported this contig in the selected virus summary table.",
    "genomad_source": "geNomad output file used for the evidence row.",
    "putative_viral": "TRUE if vs2_hit OR genomad_hit OR blastx_hit is TRUE.",
    "viral_evidence": "Semicolon-separated list of evidence tools: VirSorter2, geNomad, BLASTx.",
    "prodigal_total_orfs": "Total number of ORFs predicted by Prodigal on the contig.",
    "prodigal_proteins_faa": "Path to Prodigal protein FASTA for the sample/bin.",
    "prodigal_gff": "Path to Prodigal GFF for the sample/bin.",
    "taxid": "NCBI taxonomy ID assigned from parsed organism name.",
    "family_or_best": "Family-level taxonomy when available, otherwise the best higher rank.",
    "best_rank": "Rank used in family_or_best.",
    "lineage_ranks": "Semicolon-separated lineage from kingdom to species where available.",
    "kingdom": "Taxonomic kingdom/superkingdom field.",
    "phylum": "Taxonomic phylum.",
    "class": "Taxonomic class.",
    "order": "Taxonomic order.",
    "family": "Taxonomic family.",
    "genus": "Taxonomic genus.",
    "species": "Taxonomic species.",
    "lineage_full": "Full lineage string used in the final table.",
    "virus_group": "Broad group classification: Bacteriophage, Eukaryotic viruses, Archaeal viruses, or Unidentified.",
}

PREFIX_NOTES = {
    "vs2_": "VirSorter2 annotation column copied from the selected VirSorter2 TSV output.",
    "genomad_": "geNomad annotation column copied from the selected geNomad TSV output.",
    "checkv_": "CheckV quality/completeness annotation column copied from quality_summary.tsv. CheckV is reported, not used as a filter.",
}


def read_columns(path: Path):
    if not path.exists() or not path.is_file():
        return []
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n")
    except OSError:
        return []
    if not header:
        return []
    return header.split("\t")


def describe_column(col: str) -> str:
    if col in COLUMN_NOTES:
        return COLUMN_NOTES[col]
    for prefix, note in PREFIX_NOTES.items():
        if col.startswith(prefix):
            return note
    return "Pipeline output column; see the source tool documentation or upstream table for exact interpretation."


def file_section(filename: str) -> str:
    path = MERGED_DIR / filename
    lines = [f"### `{filename}`", "", FILE_DESCRIPTIONS.get(filename, "Pipeline output file."), ""]
    cols = read_columns(path)
    if cols:
        lines.append("Columns:")
        lines.append("")
        for col in cols:
            lines.append(f"- `{col}`: {describe_column(col)}")
        lines.append("")
    else:
        lines.append("Columns: not listed because the file was not present when this README was generated.")
        lines.append("")
    return "\n".join(lines)


def main():
    OUT_BASE.mkdir(parents=True, exist_ok=True)
    MERGED_DIR.mkdir(parents=True, exist_ok=True)

    lines = [
        "# Air Non-Hybrid Virome Pipeline Outputs",
        "",
        "This README is generated automatically by `write_pipeline_readme.py` after the pipeline finishes.",
        "",
        "## How To Run",
        "",
        "Run the complete workflow with:",
        "",
        "```bash",
        "bash run_air_virome_pipeline.sh",
        "```",
        "",
        "Sankey outputs are intentionally not generated in this workflow.",
        "",
        "## Samples",
        "",
        ", ".join(SAMPLES),
        "",
        "## Output Structure",
        "",
        "```text",
        "viral_id_pipeline/",
        "  logs/      pipeline and per-sample log files",
        "  samples/   per-sample outputs and intermediate files",
        "  merged/    cohort-level merged tables and final tables",
        "  taxonomy/  NCBI taxdump and taxonomy intermediate files",
        "```",
        "",
        "Per-sample directories under `samples/<sample>/` contain:",
        "",
        "- `00_input/`: symlinks to input reads and assembly.",
        "- `01_filter_lenbin/`: filtered contigs and length-bin FASTA files.",
        "- `02_mapping/`: read-to-contig BAM files.",
        "- `03_depth_tiers/`: depth-filtered contig FASTA files.",
        "- `04_virsorter2/`: VirSorter2 outputs.",
        "- `05_genomad/`: geNomad outputs.",
        "- `06_checkv/`: CheckV outputs.",
        "- `07_diamond_ge1000/`: Prodigal and DIAMOND outputs for contigs >=1000 bp.",
        "- `07_diamond_len500_999/`: Prodigal and DIAMOND outputs for 500-999 bp contigs.",
        "- `08_coverage/`: coverage and read-count tables.",
        "",
        "## Viral Evidence Rule",
        "",
        "`putative_viral` is TRUE when at least one of these is TRUE:",
        "",
        "- `vs2_hit`",
        "- `genomad_hit`",
        "- `blastx_hit`",
        "",
        "CheckV results are included as `checkv_*` columns for interpretation, but CheckV quality/completeness is not used to remove contigs.",
        "",
        "## Merged Files",
        "",
    ]

    for filename in FILE_DESCRIPTIONS:
        lines.append(file_section(filename))

    README.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")
    print(f"[OK] wrote {README}")


if __name__ == "__main__":
    main()
