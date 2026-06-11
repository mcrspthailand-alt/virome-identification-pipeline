#!/usr/bin/env python3
import argparse
import os
import re
import pandas as pd
import numpy as np


def parse_organism_series(stitle: pd.Series) -> pd.Series:
    s = stitle.fillna("").astype(str).str.strip()
    os_name = s.str.extract(r'OS=([^=]+?)(?:\s+OX=|\s+GN=|\s+PE=|\s+SV=|$)', expand=False)
    bracket = s.str.extract(r'\[([^\[\]]+)\]\s*$', expand=False)
    out = os_name.fillna("")
    out = out.mask(out.eq(""), bracket.fillna(""))
    out = out.mask(out.eq(""), s)
    out = out.str.replace(r'\s+', ' ', regex=True).str.strip()
    return out


def split_lineage(df: pd.DataFrame) -> pd.DataFrame:
    if 'lineage_ranks' not in df.columns:
        df['lineage_ranks'] = ''
    parts = df['lineage_ranks'].fillna('').astype(str).str.split(';', expand=True)
    cols = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for i, c in enumerate(cols):
        if i in parts.columns:
            df[c] = parts[i].fillna('').astype(str).str.strip()
        else:
            df[c] = ''
    df['lineage_full'] = df[cols].replace('', np.nan).apply(
        lambda r: ';'.join([x for x in r.tolist() if isinstance(x, str) and x.strip()]), axis=1
    )
    return df


def classify(df: pd.DataFrame) -> pd.DataFrame:
    lineage = df.get('lineage_full', pd.Series('', index=df.index)).fillna('').astype(str)
    family_best = df.get('family_or_best', pd.Series('', index=df.index)).fillna('').astype(str)
    organism = df.get('blastx_organism', pd.Series('', index=df.index)).fillna('').astype(str)
    stitle = df.get('stitle', pd.Series('', index=df.index)).fillna('').astype(str)
    text = (lineage + ' || ' + family_best + ' || ' + organism + ' || ' + stitle).str.lower()

    phage_pat = r'caudovir|ackermannviridae|autographiviridae|chaseviridae|chimalliviridae|demerecviridae|drexlerviridae|fiersviridae|herelleviridae|inoviridae|leviviridae|microviridae|myoviridae|podoviridae|siphoviridae|schitoviridae|tectiviridae|phage|bacteriophage|cyanophage'
    archaeal_pat = r'archaeal|archaea|rudiviridae|fuselloviridae|bicaudaviridae|ampullaviridae|guttaviridae|lipothrixviridae|clavaviridae|tristromaviridae|hafunaviridae'
    euk_pat = r'alloherpesviridae|amalgaviridae|ascoviridae|baculoviridae|botourmiaviridae|chuviridae|circoviridae|endornaviridae|eupolintoviridae|flaviviridae|hepeviridae|hypoviridae|iridoviridae|marseilleviridae|mimiviridae|mitoviridae|nairoviridae|nimaviridae|orthoherpesviridae|orthomyxoviridae|orthototiviridae|partitiviridae|parvoviridae|phenuiviridae|phycodnaviridae|pithoviridae|pneumoviridae|polydnaviriformidae|polymycoviridae|poxviridae|prasinovirus|pseudoviridae|reovir|retroviridae|sedoreoviridae|solemoviridae|togaviridae|tombusviridae|pandoravirus|herpesvirus|poxvirus|retrovirus|partitivirus|totivirus|mitovirus|polymycovirus|adintovirus|circovirus|chlorella virus|ostreococcus|emiliania huxleyi virus|large algae virus|ranavirus|bracovirus|nuclear polyhedrosis virus|papillom|adenovirus|coronavirus|orthomyxo'

    phage = text.str.contains(phage_pat, regex=True, na=False)
    archaeal = (~phage) & text.str.contains(archaeal_pat, regex=True, na=False)
    euk = (~phage) & (~archaeal) & text.str.contains(euk_pat, regex=True, na=False)

    df['virus_group'] = np.select(
        [phage, archaeal, euk],
        ['Bacteriophage', 'Archaeal viruses', 'Eukaryotic viruses'],
        default='Unidentified'
    )
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--merged-dir', required=True)
    args = ap.parse_args()

    merged_dir = args.merged_dir
    master_path = os.path.join(merged_dir, 'ALL_samples.master_option1.combined.fixed.with_mapped_reads.tsv')
    tax_path = os.path.join(merged_dir, 'virus_taxonomy_map.tsv')

    if not os.path.exists(master_path):
        raise SystemExit(f'[ERROR] missing {master_path}')
    if not os.path.exists(tax_path):
        raise SystemExit(f'[ERROR] missing {tax_path}')

    master = pd.read_csv(master_path, sep='\t', dtype=str, low_memory=False)
    tax = pd.read_csv(tax_path, sep='\t', dtype=str, low_memory=False)

    for c in master.columns:
        master[c] = master[c].fillna('')

    if 'length' not in master.columns:
        master['length'] = ''
    if 'contig_length' in master.columns:
        master['length'] = np.where(
            master['length'].fillna('').astype(str).str.strip().eq(''),
            master['contig_length'].fillna('').astype(str).str.strip(),
            master['length'].fillna('').astype(str).str.strip()
        )

    if 'stitle' not in master.columns:
        master['stitle'] = ''
    if 'blastx_organism' not in master.columns:
        master['blastx_organism'] = ''
    master['blastx_organism_raw'] = master['blastx_organism'].fillna('').astype(str)
    master['blastx_organism'] = parse_organism_series(master['stitle'])

    for c in tax.columns:
        tax[c] = tax[c].fillna('').astype(str)
    if 'blastx_organism' not in tax.columns:
        raise SystemExit('[ERROR] tax map has no blastx_organism column')
    tax = tax.drop_duplicates(subset=['blastx_organism'], keep='first')

    out = master.merge(tax, on='blastx_organism', how='left', suffixes=('', '_tax'))
    for c in ['taxid', 'family_or_best', 'best_rank', 'lineage_ranks']:
        if c not in out.columns:
            out[c] = ''
        out[c] = out[c].fillna('')

    out = split_lineage(out)
    out = classify(out)

    front = [c for c in [
        'sample', 'contig', 'len_bin', 'length', 'contig_length', 'mapped_reads', 'unmapped_reads',
        'pident', 'aln_len', 'evalue', 'bitscore', 'n_orfs_strict',
        'blastx_organism', 'blastx_organism_raw', 'stitle', 'sseqid',
        'taxid', 'family_or_best', 'best_rank',
        'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'lineage_ranks', 'lineage_full',
        'virus_group'
    ] if c in out.columns]
    rest = [c for c in out.columns if c not in front]
    out = out[front + rest]

    out_path = os.path.join(merged_dir, 'final_annotated_with_taxonomy.tsv')
    phage_path = os.path.join(merged_dir, 'final_phage_only.tsv')
    euk_path = os.path.join(merged_dir, 'final_eukaryotic_only.tsv')
    arch_path = os.path.join(merged_dir, 'final_archaeal_only.tsv')
    unk_path = os.path.join(merged_dir, 'final_unidentified_only.tsv')
    report_path = os.path.join(merged_dir, 'final_group_report.txt')

    out.to_csv(out_path, sep='\t', index=False)
    out[out['virus_group'] == 'Bacteriophage'].to_csv(phage_path, sep='\t', index=False)
    out[out['virus_group'] == 'Eukaryotic viruses'].to_csv(euk_path, sep='\t', index=False)
    out[out['virus_group'] == 'Archaeal viruses'].to_csv(arch_path, sep='\t', index=False)
    out[out['virus_group'] == 'Unidentified'].to_csv(unk_path, sep='\t', index=False)

    with open(report_path, 'w') as fh:
        fh.write(f'rows\t{len(out)}\n')
        fh.write('group_counts\n')
        fh.write(out['virus_group'].value_counts(dropna=False).to_string() + '\n')
        fh.write(f'nonempty_length\t{int(out["length"].fillna("").astype(str).str.strip().ne("").sum())}\n')
        fh.write(f'nonempty_lineage\t{int(out["lineage_full"].fillna("").astype(str).str.strip().ne("").sum())}\n')

    print('[OK] wrote', out_path)
    print('[OK] wrote', phage_path)
    print('[OK] wrote', euk_path)
    print('[OK] wrote', arch_path)
    print('[OK] wrote', unk_path)
    print('[OK] wrote', report_path)


if __name__ == '__main__':
    main()
