# Cynopterus Bat Virome Pipeline

This folder is a ready-to-run virome identification pipeline for the Cynopterus bat pooled samples.

## Samples

```text
Bat-pool-1
Bat-pool-2
Bat-pool-5
Bat-pool-6
```

## Input Paths

Host-removed reads:

```text
/home/panda/workspace/projects/nextflow_metagenomic/PANDA/bat/results_Cynopterus/QC_shortreads/remove_host
```

MEGAHIT assemblies:

```text
/home/panda/workspace/projects/nextflow_metagenomic/PANDA/bat/results_Cynopterus/Assembly/MEGAHIT
```

The exact input files are listed in:

```text
samples.tsv
```

## Output Path

Results will be written to:

```text
/home/panda/workspace/projects/nextflow_metagenomic/PANDA/bat/results_Cynopterus/viral_id_pipeline
```

## Run

On the HPC:

```bash
cd /path/to/bat_cynopterus_virome_pipeline
bash run_virome_pipeline.sh
```

or explicitly:

```bash
bash run_virome_pipeline.sh config.env
```

## Workflow

The pipeline runs:

1. Read/assembly preparation and contig filtering
2. Read mapping to contigs with `minimap2`
3. Read counts and depth with `samtools` and `mosdepth`
4. VirSorter2
5. geNomad
6. CheckV
7. Prodigal
8. DIAMOND `blastx`
9. DIAMOND `blastp`
10. Evidence-table merge
11. NCBI taxonomy assignment
12. Final viral-contig tables
13. Output README generation

## Viral Calling Rule

```text
putative_viral = vs2_hit OR genomad_hit OR blastx_hit OR blastp_hit
```

Confidence:

```text
High    = viral_evidence_count >= 2
Medium  = viral_evidence_count == 1
None    = viral_evidence_count == 0
```

CheckV is included as annotation only and is not used to discard contigs.

## Important Outputs

After the run, check:

```text
/home/panda/workspace/projects/nextflow_metagenomic/PANDA/bat/results_Cynopterus/viral_id_pipeline/merged/all_viral_evidence.tsv
/home/panda/workspace/projects/nextflow_metagenomic/PANDA/bat/results_Cynopterus/viral_id_pipeline/merged/final_virus_contigs.tsv
/home/panda/workspace/projects/nextflow_metagenomic/PANDA/bat/results_Cynopterus/viral_id_pipeline/merged/final_virus_summary.txt
```

The generated output documentation will be:

```text
/home/panda/workspace/projects/nextflow_metagenomic/PANDA/bat/results_Cynopterus/viral_id_pipeline/README.md
```

