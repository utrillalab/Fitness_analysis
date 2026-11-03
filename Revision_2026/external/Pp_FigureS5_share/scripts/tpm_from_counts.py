"""Compute TPM for all GSE85475 samples from raw counts.

Input:    GEO processed file GSE85475_KBojanovic_Processed_data_Stress.xlsx.
Samples:  every 'Raw Counts <sample>' column (21 samples).
Gene set: PP_ locus tags, excluding tRNA and rRNA (from the NCBI feature table).
Length:   CDS length from the Excel, |Translation Stop - Translation Start| + 1
          (transcript coordinates are missing for ~57% of genes).
Output:   TSV with locus_tag, length_bp and one TPM_<sample> column per sample.

Usage (from the project root):
    python scripts/tpm_from_counts.py
"""
import argparse

import numpy as np
import pandas as pd

COUNT_PREFIX = 'Raw Counts '


def load_counts(xlsx, sheet):
    """Excel -> (DataFrame with locus_tag, CDS coordinates and raw counts, sample names)."""
    df = pd.read_excel(xlsx, sheet_name=sheet)
    samples = [c[len(COUNT_PREFIX):] for c in df.columns if c.startswith(COUNT_PREFIX)]
    cols = ['Synonym', 'Translation Start', 'Translation Stop'] + \
           [COUNT_PREFIX + s for s in samples]
    df = df[cols].rename(columns={'Synonym': 'locus_tag'})
    return df, samples


def load_rna_tags(feature_table):
    """NCBI feature table -> set of tRNA/rRNA locus tags."""
    ft = pd.read_table(feature_table, usecols=['# feature', 'locus_tag'])
    return set(ft.loc[ft['# feature'].isin(['tRNA', 'rRNA']), 'locus_tag'])


def filter_genes(df, rna_tags):
    """Keep PP_ locus tags that are not tRNA/rRNA."""
    is_pp = df['locus_tag'].astype(str).str.match(r'^PP_')
    return df[is_pp & ~df['locus_tag'].isin(rna_tags)].copy()


def cds_length(df):
    """Add length_bp; drop rows with missing coordinates or length <= 0."""
    length = (df['Translation Stop'] - df['Translation Start']).abs() + 1
    df = df.assign(length_bp=length)
    bad = df['length_bp'].isna() | (df['length_bp'] <= 0)
    if bad.any():
        print(f'Dropped {bad.sum()} genes without valid coordinates: '
              f'{list(df.loc[bad, "locus_tag"])}')
    return df[~bad]


def compute_tpm(counts, length_bp):
    """Raw counts (genes x samples) + length -> TPM; each column sums to 1e6."""
    rate = counts.div(length_bp / 1e3, axis=0)
    tpm = rate.div(rate.sum(axis=0), axis=1) * 1e6
    assert np.allclose(tpm.sum(axis=0), 1e6), 'TPM columns do not sum to 1e6'
    return tpm


def build_table(df, tpm, samples):
    """Gene table + TPM -> output table (locus_tag, length_bp, TPM_<sample>)."""
    out = df[['locus_tag', 'length_bp']].copy()
    out[[f'TPM_{s}' for s in samples]] = tpm.values
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--xlsx', default='data/raw/GSE85475_KBojanovic_Processed_data_Stress.xlsx')
    p.add_argument('--sheet', default='KBojanovic_Processed_data_Stres')
    p.add_argument('--feature-table', default='data/raw/GCA_000007565_feature_table.txt')
    p.add_argument('--out', default='data/processed/GSE85475_TPM.tsv')
    a = p.parse_args()

    df, samples = load_counts(a.xlsx, a.sheet)
    df = filter_genes(df, load_rna_tags(a.feature_table))
    df = cds_length(df)
    tpm = compute_tpm(df[[COUNT_PREFIX + s for s in samples]], df['length_bp'])
    out = build_table(df, tpm, samples)
    out.to_csv(a.out, sep='\t', index=False)

    print(f'Samples: {len(samples)}  Genes: {len(out)}  ->  {a.out}')


if __name__ == '__main__':
    main()
