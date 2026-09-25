import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


MODULE = Path(__file__).resolve().parents[1]


def tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class ZoverTemplateTest(unittest.TestCase):
    def test_translated_hsp_group_and_curation(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "verified"
            queries = root / "queries"
            queries.mkdir(parents=True)
            (queries / "sample-1_merged_contigs.fa").write_text(">contig_A\nATGAAAACC\n", encoding="utf-8")
            annotation = Path(tmp) / "annotation.tsv"
            fields = ["sample", "contig", "virus_family", "gene_orf", "reference_protein_accession", "tree_group_id", "blastx_bitscore"]
            tsv(annotation, fields, [{"sample": "sample-1", "contig": "contig_A", "virus_family": "Testviridae", "gene_orf": "RdRp", "reference_protein_accession": "ABC12345.1", "tree_group_id": "Testviridae__RdRp", "blastx_bitscore": "80"}])
            protein = Path(tmp) / "reference.faa"
            protein.write_text(">gb|ABC12345.1| polymerase\nMKTQQ\n", encoding="utf-8")
            blastx = Path(tmp) / "blastx.tsv"
            blastx.write_text("Q__sample-1__contig_A\tREF__ABC12345.1\t100\t3\t1\t9\t1\t3\t1e-10\t80\t1\tMKT\n", encoding="utf-8")
            output = Path(tmp) / "tree"
            result = subprocess.run([
                sys.executable, str(MODULE / "prepare_protein_tree_fastas.py"),
                "--root", str(root), "--annotation", str(annotation),
                "--protein-reference", str(protein), "--output", str(output),
                "--blastx-tsv", str(blastx),
            ], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            group = output / "groups" / "Testviridae__RdRp"
            self.assertIn("MKT", (group / "query_blastx_hsp_proteins.faa").read_text())
            self.assertIn("MKTQQ", (group / "reference_proteins.faa").read_text())
            with (output / "protein_tree_manifest.tsv").open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(row["hsp_qframe"], "1")
            self.assertEqual(row["reference_protein_accession"], "ABC12345.1")


if __name__ == "__main__":
    unittest.main()

