#!/usr/bin/env python3
"""
Bulk DNA/RNA DARLIN Array Processing Script
"""

import os
import sys
import gzip
import subprocess
import hashlib
import argparse
import logging
import shutil
from collections import Counter
from typing import Optional, Union, Sequence
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from fuzzysearch import find_near_matches
from tqdm import tqdm
from umi_tools import UMIClusterer
from darlinpy import analyze_sequences
from darlinpy.config.amplicon_configs import load_carlin_config_by_locus


# ============================================================================
# Argument Parsing and Logging Setup
# ============================================================================

def setup_logging(log_file, log_level=logging.INFO):
    """Setup logging configuration."""
    logger = logging.getLogger(__name__)
    logger.setLevel(log_level)
    
    # Remove existing handlers to avoid duplicates
    logger.handlers = []
    
    # Create formatter - simplified format with time to seconds only
    formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s', 
                                  datefmt='%Y-%m-%d %H:%M:%S')
    
    # File handler
    file_handler = logging.FileHandler(log_file, mode='w')
    file_handler.setLevel(log_level)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    
    # Console handler
    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(log_level)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    
    return logger


def parse_arguments():
    """Parse command line arguments."""
    # Custom formatter that preserves formatting and shows defaults
    class CustomFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass
    
    example_text = '''
Examples:
  Basic usage:
    python run_DARLIN_bulk.py \\
        --sample-id LL583_bulk_DNA_CA \\
        --fq1 ./fastq/GSM6924444_CC-DNA_LL583/SRR23027633_1.fastq.gz \\
        --fq2 ./fastq/GSM6924444_CC-DNA_LL583/SRR23027633_2.fastq.gz \\
        --umi-len 12 \\
        --locus Col1a1 \\
        --output-dir ./output \\
        --reads-cutoff 1
        

  Skip PEAR assembly (use existing assembled file):
    python run_DARLIN_bulk.py \\
        --sample-id LL583_bulk_DNA_CA \\
        --fq1 ./fastq/GSM6924444_CC-DNA_LL583/SRR23027633_1.fastq.gz \\
        --fq2 ./fastq/GSM6924444_CC-DNA_LL583/SRR23027633_2.fastq.gz \\
        --locus Col1a1 \\
        --skip-pear

  Keep PEAR output directory:
    python run_DARLIN_bulk.py \\
        --sample-id LL583_bulk_DNA_CA \\
        --fq1 ./fastq/GSM6924444_CC-DNA_LL583/SRR23027633_1.fastq.gz \\
        --fq2 ./fastq/GSM6924444_CC-DNA_LL583/SRR23027633_2.fastq.gz \\
        --locus Col1a1 \\
        --keep-pear
'''
    
    parser = argparse.ArgumentParser(
        description='Bulk DNA/RNA DARLIN Array Processing Script',
        formatter_class=CustomFormatter,
        epilog=example_text
    )
    
    # Required arguments
    parser.add_argument('--sample-id', type=str, required=True,
                        help='Sample ID for output directory naming')
    parser.add_argument('--fq1', type=str, required=True,
                        help='Path to forward reads FASTQ file')
    parser.add_argument('--fq2', type=str, required=True,
                        help='Path to reverse reads FASTQ file')
    # Optional arguments
    parser.add_argument('--locus', type=str, default='Col1a1',
                        help='Locus name, options: Col1a1, Rosa, Tigre')
    parser.add_argument('--umi-len', type=int, default=12,
                        help='UMI length')
    parser.add_argument('--output-dir', type=str, default='./output',
                        help='Base output directory')
    parser.add_argument('--pear-path', type=str, default='pear',
                        help='Path to PEAR executable')
    parser.add_argument('--threads', type=int, default=8,
                        help='Number of threads for PEAR')
    parser.add_argument('--min-bc-len', type=int, default=20,
                        help='Minimum barcode length')
    parser.add_argument('--reads-cutoff', type=int, default=1,
                        help='Reads cutoff threshold')
    parser.add_argument('--denoise-iter', type=int, default=1,
                        help='Number of denoising iterations')
    parser.add_argument('--umi-ld', type=int, nargs='+', default=[1],
                        help='List of maximum edit/Hamming distances for UMI clustering (e.g. 1 2 3)')
    parser.add_argument('--lb-hd-relative', type=float, nargs='+', default=[0.01],
                        help='List of relative Hamming-distance thresholds for lineage barcode clustering (e.g. 0.01 0.02)')
    parser.add_argument('--skip-pear', action='store_true',
                        help='Skip PEAR assembly (use existing assembled file)')
    parser.add_argument('--keep-pear', action='store_true',
                        help='Keep PEAR output directory after completion')
    parser.add_argument('--log-level', type=str, default='INFO',
                        choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'],
                        help='Logging level')
    
    return parser.parse_args()


# ============================================================================
# Utility Functions
# ============================================================================

def open_fastq_file(file_path):
    """Open text or gz FASTQ in text mode."""
    if file_path.endswith(".gz"):
        return gzip.open(file_path, "rt")
    return open(file_path, "r")


def iter_fastq_raw(handle):
    """
    Minimal FASTQ iterator: read 4 lines per record.
    Returns (read_id, seq, qual) as strings. No conversion on quality.
    """
    while True:
        id_line = handle.readline()
        if not id_line:
            break
        seq_line = handle.readline()
        plus_line = handle.readline()
        qual_line = handle.readline()

        if not (seq_line and plus_line and qual_line):
            raise ValueError("Incomplete FASTQ record encountered.")

        if not id_line.startswith("@") or not plus_line.startswith("+"):
            raise ValueError("Invalid FASTQ structure (missing @ or + line).")

        read_id = id_line[1:].strip()
        seq = seq_line.strip()
        qual = qual_line.strip()
        if len(seq) != len(qual):
            raise ValueError(f"Length mismatch (seq {len(seq)} vs qual {len(qual)}) at read {read_id}")
        yield read_id, seq, qual


def get_mm_dist(seq, rate=0.05, n=2):
    # allow 5% mismatch, at least n
    return max(int(len(seq) * rate), n)


def find_exact_matches(seq_str, patterns):
    """
    Fast exact matching using Python's built-in string.find()
    Returns: dict with pattern_name -> [positions] mapping
    """
    matches = {}
    for pattern_name, pattern_seq in patterns.items():
        positions, start = [], 0
        while True:
            pos = seq_str.find(pattern_seq, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1  # allow overlapping
        matches[pattern_name] = positions
    return matches


def find_fuzzy_matches(seq_str, pattern, max_errors):
    """Fallback fuzzy matching (Levenshtein) when exact matching fails."""
    return find_near_matches(pattern, seq_str, max_l_dist=max_errors)


class MatchResult:
    """Container for match results"""
    def __init__(self):
        self.p3_matches = []
        self.p5_matches = []
        self.match_method = ""  # "exact", "fuzzy", "mixed", or "failed"


def find_all_matches(seq_str, p3_seq, p5_seq, p3_mm, p5_mm):
    """
    Try exact matches first; fall back to fuzzy per element if needed.
    Require exactly one hit for each element.
    """
    result = MatchResult()
    patterns = {'p3': p3_seq, 'p5': p5_seq}
    exact = find_exact_matches(seq_str, patterns)

    exact_success = (len(exact['p3']) == 1 and
                     len(exact['p5']) == 1)
    if exact_success:
        result.p3_matches = [type('Match', (), {'start': exact['p3'][0], 'end': exact['p3'][0] + len(p3_seq)})()]
        result.p5_matches = [type('Match', (), {'start': exact['p5'][0], 'end': exact['p5'][0] + len(p5_seq)})()]
        result.match_method = "exact"
        return result

    mixed_success = True
    methods_used = []

    # p3
    if len(exact['p3']) == 1:
        result.p3_matches = [type('Match', (), {'start': exact['p3'][0], 'end': exact['p3'][0] + len(p3_seq)})()]
        methods_used.append('exact')
    else:
        result.p3_matches = find_fuzzy_matches(seq_str, p3_seq, p3_mm)
        methods_used.append('fuzzy')
        if len(result.p3_matches) != 1:
            mixed_success = False

    # p5
    if mixed_success:
        if len(exact['p5']) == 1:
            result.p5_matches = [type('Match', (), {'start': exact['p5'][0], 'end': exact['p5'][0] + len(p5_seq)})()]
            methods_used.append('exact')
        else:
            result.p5_matches = find_fuzzy_matches(seq_str, p5_seq, p5_mm)
            methods_used.append('fuzzy')
            if len(result.p5_matches) != 1:
                mixed_success = False

    if mixed_success:
        result.match_method = "mixed" if 'fuzzy' in methods_used else "exact"
    else:
        result.match_method = "failed"
        result.p3_matches = []
        result.p5_matches = []

    return result


def hamming1(a, b):
    if len(a) != len(b):
        return False
    diff = 0
    for ca, cb in zip(a, b):
        if ca != cb:
            diff += 1
            if diff > 1:
                return False
    return True


def hamming_dist(a: str, b: str) -> int:
    """Return Hamming distance (requires equal length)."""
    if len(a) != len(b):
        return max(len(a), len(b))  # treat unequal length as very large
    return sum(c1 != c2 for c1, c2 in zip(a, b))


def neighbors_hd1(s):
    bases = ('A', 'C', 'G', 'T')
    out = []
    for i, ch in enumerate(s):
        for b in bases:
            if b != ch:
                out.append(s[:i] + b + s[i+1:])
    return out


def collapse_directional(items):
    """
    items: iterable of (seq, count)
    Rule: parent absorbs neighbor if count_parent >= 2*count_child - 1
    """
    counts = dict(items)
    ordered = sorted(counts, key=lambda s: counts[s], reverse=True)
    parent = {s: s for s in ordered}
    index = set(ordered)

    for s in ordered:
        if parent[s] != s:  # already absorbed
            continue
        c_hi = counts[s]
        for nb in neighbors_hd1(s):
            if nb in index and parent[nb] == nb:
                c_lo = counts[nb]
                if c_hi >= 2*c_lo - 1:
                    parent[nb] = s
    return parent


def collapse_within_hd(items, max_hd: int):
    """
    Like collapse_directional but allows Hamming distance <= max_hd (>=2 typical).
    Naive O(n^2) within a group, which is OK because we operate per-UMI buckets.
    """
    counts = dict(items)
    seqs = sorted(counts, key=lambda s: counts[s], reverse=True)
    parent = {s: s for s in seqs}

    for i, s in enumerate(seqs):
        if parent[s] != s:
            continue
        c_hi = counts[s]
        # only consider candidates not yet absorbed
        for t in seqs[i+1:]:
            if parent[t] != t:
                continue
            if len(t) != len(s):
                continue
            if hamming_dist(s, t) <= max_hd:
                c_lo = counts[t]
                if c_hi >= 2 * c_lo - 1:
                    parent[t] = s
    return parent


def _to_df(
    data: Union[pd.DataFrame, Sequence, "np.ndarray"],
    bc_col: str,
    umi_col: str,
    count_col: Optional[str]
) -> pd.DataFrame:
    """Coerce various inputs to a clean DataFrame with required columns."""
    if isinstance(data, pd.DataFrame):
        df = data.copy()
    else:
        # list[dict], list[tuple], np.ndarray → DataFrame
        df = pd.DataFrame(data)

    # try to infer columns if unnamed tuples/arrays were provided
    if bc_col not in df.columns or umi_col not in df.columns:
        # common fallback for 2-3 column structures
        if df.shape[1] >= 2 and bc_col not in df.columns and umi_col not in df.columns:
            # name first three columns defensively
            cols = list(df.columns)
            rename_map = {}
            if len(cols) >= 1:
                rename_map[cols[0]] = bc_col
            if len(cols) >= 2:
                rename_map[cols[1]] = umi_col
            if len(cols) >= 3 and count_col and count_col not in df.columns:
                rename_map[cols[2]] = count_col
            df = df.rename(columns=rename_map)

    # minimal checks
    missing = [c for c in [bc_col, umi_col] if c not in df.columns]
    if missing:
        raise ValueError(f"Input is missing required column(s): {missing}. "
                         f"Available: {list(df.columns)}")

    # clean up
    for c in [bc_col, umi_col]:
        df[c] = df[c].astype(str).str.upper()
    df = df.dropna(subset=[bc_col, umi_col])

    if count_col and count_col in df.columns:
        df[count_col] = pd.to_numeric(df[count_col], errors="coerce").fillna(1).astype(int)
    else:
        df["__count__"] = 1
        count_col = "__count__"

    return df, count_col


def correct_lineage_and_umi(
    data: Union[pd.DataFrame, Sequence],
    umi_col: str = "UMI",
    bc_col: str = "lineage_bc",
    count_col: Optional[str] = None,
    n_iter: int = 2,
    umi_ld: int = 2,
    lb_hd_relative: float = 0.01,
    logger=None,
):
    """
    Accepts DataFrame, list[dict], list[tuple], or ndarray.
    Returns:
      agg: deduplicated (lineage_bc_corr, umi_corr, n_reads)
      mapping: original→corrected pairs
      stats: dict with merge counts
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    df, count_col = _to_df(data, bc_col=bc_col, umi_col=umi_col, count_col=count_col)
    out = df[[bc_col, umi_col, count_col]].copy()
    out["umi_corr"] = out[umi_col]
    out["lineage_bc_corr"] = out[bc_col]

    total_umi_merges = 0
    total_bc_merges = 0

    # umi_tools clusterer (directional strategy, distance threshold controlled by umi_ld)
    clusterer = UMIClusterer(cluster_method="directional")

    for i in range(n_iter):
        logger.info(f"Iteration {i+1}/{n_iter}")
        # 1) collapse lineage_bc globally by length, using length-aware Hamming threshold
        out["__bc_len__"] = out["lineage_bc_corr"].str.len()
        bc_parent_total = {}
        for blen, sub in tqdm(out.groupby(["__bc_len__"]), desc="Collapsing barcodes (length-aware HD global)", leave=True):
            cnt = Counter(dict(sub.groupby("lineage_bc_corr")[count_col].sum()))
            if blen == 0:
                continue
            # threshold scales with length; ensure at least distance 1
            hd_thresh_len = max(int(round(lb_hd_relative * blen)), 1)
            parent = collapse_within_hd(cnt.items(), max_hd=hd_thresh_len)
            bc_parent_total.update(parent)
        before = out["lineage_bc_corr"].ne(out["lineage_bc_corr"].map(lambda b: bc_parent_total.get(b, b))).sum()
        out["lineage_bc_corr"] = out["lineage_bc_corr"].map(lambda b: bc_parent_total.get(b, b))
        total_bc_merges += int(before)

        # 2) collapse UMIs within each lineage_bc using umi_tools (directional clustering)
        parent_all = {}
        for bc_val, sub in tqdm(out.groupby("lineage_bc_corr"), desc="Collapsing UMIs with umi_tools", leave=True):
            cnt_series = sub.groupby("umi_corr")[count_col].sum()
            if cnt_series.empty:
                continue
            # umi_tools expects dict[bytes, int]
            umi_counts_bytes = {umi.encode(): int(c) for umi, c in cnt_series.items()}
            umi_groups = clusterer(umi_counts_bytes, threshold=umi_ld)

            for group in umi_groups:
                # representative: UMI with max count within the group
                representative_bytes = max(group, key=lambda u: umi_counts_bytes.get(u, 0))
                representative_str = representative_bytes.decode()
                for umi_bytes in group:
                    umi_str = umi_bytes.decode()
                    parent_all[umi_str] = representative_str

        before = out["umi_corr"].ne(out["umi_corr"].map(lambda u: parent_all.get(u, u))).sum()
        out["umi_corr"] = out["umi_corr"].map(lambda u: parent_all.get(u, u))
        total_umi_merges += int(before)

    agg = (out.groupby(["lineage_bc_corr", "umi_corr"], as_index=False)[count_col]
           .sum()
           .rename(columns={count_col: "n_reads"}))

    mapping = out[[bc_col, umi_col, "lineage_bc_corr", "umi_corr"]].copy().drop_duplicates()

    stats = {
        "n_input_rows": int(len(df)),
        "n_unique_pairs_before": int(df.groupby([bc_col, umi_col]).size().shape[0]),
        "n_unique_pairs_after": int(agg.shape[0]),
        "umi_merges": total_umi_merges,
        "barcode_merges": total_bc_merges,
    }
    return agg, mapping, stats


# ============================================================================
# Main Processing Functions
# ============================================================================

def assemble_pe_reads(in_fq1, in_fq2, out_fq, log, pear_path="~/software/pear", threads=8, logger=None):
    """
    Assemble paired-end reads using PEAR.
    
    Args:
        in_fq1: Path to forward reads
        in_fq2: Path to reverse reads
        out_fq: Output prefix for assembled reads
        log: Path to log file
        pear_path: Path to PEAR executable
        threads: Number of threads to use
        logger: Logger instance
    """
    if logger is None:
        logger = logging.getLogger(__name__)
    
    pear_path = os.path.expanduser(pear_path)
    cmd = [
        pear_path,
        "-f", in_fq1,
        "-r", in_fq2,
        "-o", out_fq,
        "-j", str(threads)
    ]
    
    logger.info(f"Running PEAR assembly: {' '.join(cmd)}")
    with open(log, 'w') as log_file:
        subprocess.run(cmd, stdout=log_file, stderr=subprocess.STDOUT, check=True)
    logger.info(f"PEAR assembly completed. Log saved to: {log}")


def extract_lineage_barcode_and_umi(fq_file, umi_len, p3_seq, p5_seq):
    """
    Extract lineage barcode and UMI from assembled FASTQ file.
    
    Returns:
        List of tuples (lineage_bc, umi)
    """
    results = []
    p3_mm = get_mm_dist(p3_seq)
    p5_mm = get_mm_dist(p5_seq)
    
    with open_fastq_file(fq_file) as fq_handle:
        for (read_id, seq, qual) in tqdm(iter_fastq_raw(fq_handle), desc="Processing reads", unit_scale=True, unit=' reads'):
            umi = seq[:umi_len]
            if 'N' in umi:
                continue
            match_result = find_all_matches(seq, p3_seq, p5_seq, p3_mm, p5_mm)
            if (len(match_result.p3_matches) == 1 and len(match_result.p5_matches) == 1):
                s = match_result.p3_matches[0].end
                e = match_result.p5_matches[0].start
                lineage_bc = seq[s:e]
                results.append((lineage_bc, umi))
    
    return results


def filter_reads(results_df, min_bc_len=20, reads_cutoff=1):
    """
    Filter reads by minimum barcode length and reads cutoff.
    
    Returns:
        Filtered DataFrame
    """
    results_df = results_df.copy()
    results_df['bc_len'] = results_df['lineage_bc'].str.len()
    results_df = results_df[results_df['bc_len'] >= min_bc_len]
    
    # Count reads per unique RNA molecule (lineage barcode + UMI)
    results_df = results_df.groupby(['lineage_bc', 'UMI']).size().reset_index(name="reads")
    results_df.sort_values(by='reads', ascending=False, inplace=True)
    
    # Apply reads cutoff
    results_clean = results_df[results_df['reads'] >= reads_cutoff].copy()
    results_clean['bc_len'] = results_clean['lineage_bc'].str.len()
    
    return results_clean


def concat_and_md5(row):
    """Calculate MD5 hash of concatenated aligned_query and aligned_ref."""
    concat_str = str(row['aligned_query']) + str(row['aligned_ref'])
    md5_hash = hashlib.md5(concat_str.encode('utf-8')).hexdigest()
    return md5_hash


# ============================================================================
# Diagnostic Plotting Functions
# ============================================================================

def plot_barcode_length_distribution(results_df, output_dir, title_suffix="", ylabel="Number of reads", unedited_bc_len=None):
    """Plot distribution of lineage barcode lengths by reads."""
    plt.figure(figsize=(4, 2))
    plt.hist(results_df['bc_len'], bins=range(1, 300, 1), edgecolor='black')
    if unedited_bc_len is not None:
        plt.axvline(unedited_bc_len, color='red', linestyle='--', linewidth=.6)
    plt.xlabel('Sequence Length')
    plt.ylabel(ylabel)
    plt.title(f'Distribution of DARLIN Array Sequence\nLengths By Reads{title_suffix}')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'barcode_length_distribution_by_reads.png'), dpi=150)
    plt.close()


def plot_reads_cutoff_distribution(results_df, reads_cutoff, output_dir):
    """Plot distribution of read counts with cutoff line."""
    plt.figure(figsize=(4, 2))
    plt.hist(results_df['reads'], bins=20, edgecolor=None, color='skyblue')
    plt.xlabel('Number of Reads')
    plt.ylabel('Frequency')
    plt.axvline(reads_cutoff, color='red', linestyle='--', linewidth=.6)
    plt.title('Distribution of Read Counts')
    plt.grid(True, alpha=0.3)
    plt.yscale('log')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'reads_counts_distribution.png'), dpi=150)
    plt.close()


def plot_cutoff_vs_num_umis(results_df, reads_cutoff, output_dir):
    """Plot reads cutoff vs number of UMIs."""
    cutoff_values = np.arange(1, results_df['reads'].max()+1, 1)
    num_umis = [(results_df[results_df['reads'] >= c]['UMI'].nunique()) for c in cutoff_values]
    
    plt.figure(figsize=(4, 2))
    plt.plot(cutoff_values, num_umis, marker='o', markersize=2, linewidth=1)
    plt.axvline(reads_cutoff, color='red', linestyle='--', linewidth=.6)
    plt.xlabel('Reads Cutoff')
    plt.ylabel('Number of UMIs')
    plt.yscale('log')
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'UMIs_by_reads_cutoff.png'), dpi=150)
    plt.close()


def plot_cutoff_vs_fraction_reads_retained(results_df, reads_cutoff, output_dir):
    """Plot reads cutoff vs fraction of reads retained."""
    cutoff_values = np.arange(1, results_df['reads'].max()+1, 1)
    fraction_reads_retained = [
        results_df[results_df['reads'] >= c]['reads'].sum() / results_df['reads'].sum()
        for c in cutoff_values
    ]
    
    plt.figure(figsize=(4, 2))
    plt.plot(cutoff_values, fraction_reads_retained, marker='o', markersize=2, linewidth=1)
    plt.axvline(reads_cutoff, color='red', linestyle='--', linewidth=.6)
    plt.xlabel('Reads Cutoff')
    plt.ylabel('Frac. of Reads\nRetained')
    plt.ylim(0, 1.05)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'fraction_retained_reads_by_reads_cutoff.png'), dpi=150)
    plt.close()


def plot_barcode_length_by_umi(agg, output_dir, unedited_bc_len=276):
    """Plot distribution of DARLIN array sequence lengths by UMI after denoising."""
    plt.figure(figsize=(4, 2))
    plt.hist(agg['bc_len'], bins=range(1, 300, 1), edgecolor='black')
    plt.axvline(unedited_bc_len, color='red', linestyle='--', linewidth=.6)
    plt.xlabel('Sequence Length')
    plt.ylabel('Number of UMIs')
    plt.title('Distribution of DARLIN Array Sequence\nLengths By UMI (After denoising)')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'barcode_length_distribution_by_UMI_after_denoising.png'), dpi=150)
    plt.close()


def plot_clone_size_distribution(final, output_dir):
    """Plot clone size distribution."""
    plt.figure(figsize=(5, 2))
    plt.hist(final[final['UMIs'] < 500]['UMIs'], bins=range(0, 501, 1), color='steelblue', edgecolor=None)
    plt.yscale('log')
    plt.xlabel('UMIs (clone size)')
    plt.ylabel('Alleles (Clones)')
    plt.title('Clone Size Distribution')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'clone_size_distribution.png'), dpi=150)
    plt.close()


def plot_barcode_length_by_editing_events(final, output_dir, unedited_bc_len=276):
    """Plot distribution of DARLIN array sequence lengths by editing events."""
    plt.figure(figsize=(4, 2))
    plt.hist(final['bc_len'], bins=range(1, 300, 1), edgecolor='black')
    plt.axvline(unedited_bc_len, color='red', linestyle='--', linewidth=.6)
    plt.xlabel('Sequence Length')
    plt.ylabel('Number of Alleles')
    plt.title('Distribution of DARLIN Array Sequence\nLengths By Editing Events')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'barcode_length_distribution_by_editing_events.png'), dpi=150)
    plt.close()


# ============================================================================
# Main Execution
# ============================================================================

def main():
    """Main execution function."""
    
    # Parse arguments
    args = parse_arguments()
    
    # Setup output directories
    sample_output_dir = os.path.join(args.output_dir, args.sample_id)
    pear_output_dir = os.path.join(sample_output_dir, 'pear')
    log_file = os.path.join(sample_output_dir, f'{args.sample_id}.log')
    
    # Create output directories
    os.makedirs(sample_output_dir, exist_ok=True)
    os.makedirs(pear_output_dir, exist_ok=True)
    
    # Setup logging
    log_level = getattr(logging, args.log_level.upper())
    logger = setup_logging(log_file, log_level)
    logger.info("--------------------------------")
    logger.info(f"Starting bulk DNA/RNA DARLIN Array processing for sample: {args.sample_id}")
    logger.info(f"Output directory: {sample_output_dir}")
    logger.info(f"Arguments: {vars(args)}")
    logger.info("--------------------------------")

    # Library structure of bulk DNA/RNA (CA/RA/TA)
    # R1-adapter--UMI--3primer-rc--lineageBC--5primer-rc--R2-adapter
    # So, in assembled FQ file, we got:
    # UMI--3primer-rc--lineageBC--5primer-rc
    CONFIG_AMPLICON = load_carlin_config_by_locus(locus=args.locus)
    UNEDITED_BC_LEN = len(CONFIG_AMPLICON.carlin_sequence)
    P5_SEQ = CONFIG_AMPLICON.sequence.primer5
    P3_SEQ = CONFIG_AMPLICON.sequence.secondary_sequence + CONFIG_AMPLICON.sequence.primer3
    P3_SEQ_RC = str(Seq(P3_SEQ).reverse_complement())
    P5_SEQ_RC = str(Seq(P5_SEQ).reverse_complement())
    # P3_SEQ_RC = P3_SEQ_RC[-10:] if len(P3_SEQ_RC) > 10 else P3_SEQ_RC
    # P5_SEQ_RC = P5_SEQ_RC[:10] if len(P5_SEQ_RC) > 10 else P5_SEQ_RC
    
    try:
        # Use primer sequences (already reverse complemented at module level)
        current_umi_len = args.umi_len
        current_p3_seq = P3_SEQ_RC
        current_p5_seq = P5_SEQ_RC
        
        #########################################################
        #### Step 1: Assemble PE reads
        #########################################################
        pear_out_prefix = os.path.join(pear_output_dir, 'pear')
        pear_log = os.path.join(pear_output_dir, 'pear.log')
        assembled_fq_file = os.path.join(pear_output_dir, 'pear.assembled.fastq')
        
        if not args.skip_pear:
            logger.info("#### Step 1: Assembling paired-end reads with PEAR...")
            assemble_pe_reads(
                args.fq1, args.fq2, pear_out_prefix, pear_log,
                pear_path=args.pear_path, threads=args.threads, logger=logger
            )
        else:
            logger.info("#### Step 1: Skipping PEAR assembly (using existing assembled file)")
            if not os.path.exists(assembled_fq_file):
                logger.error(f"Assembled FASTQ file not found: {assembled_fq_file}")
                sys.exit(1)
        
        #########################################################
        #### Step 2: Extract lineage barcode and UMI
        #########################################################
        logger.info("#### Step 2: Extracting lineage barcode and UMI...")
        results = extract_lineage_barcode_and_umi(assembled_fq_file, current_umi_len, current_p3_seq, current_p5_seq)
        logger.info(f"Valid reads: {len(results)}")
        
        # Build base dataframe and pre-filter; this part is shared across all parameter combinations
        results_df = pd.DataFrame(results, columns=['lineage_bc', 'UMI'])
        results_df['bc_len'] = results_df['lineage_bc'].str.len()
        
        # Filter by minimum barcode length
        results_df = results_df[results_df['bc_len'] >= args.min_bc_len]
        
        # Count reads per unique RNA molecule (lineage barcode + UMI)
        results_df = results_df.groupby(['lineage_bc', 'UMI']).size().reset_index(name="reads")
        results_df.sort_values(by='reads', ascending=False, inplace=True)
        
        # Apply reads cutoff
        results_clean = results_df[results_df['reads'] >= args.reads_cutoff].copy()
        results_clean['bc_len'] = results_clean['lineage_bc'].str.len()
        
        logger.info(f"Reads after filtering: {len(results_clean)}")
        logger.info(f"Proportion of valid reads: {results_clean['reads'].sum() / results_df['reads'].sum():.4f}")

        # Basic sanity check: lists must be non-empty
        if not args.umi_ld or not args.lb_hd_relative:
            logger.error("umi_ld and lb_hd_relative must both be non-empty lists.")
            sys.exit(1)

        # Loop over full Cartesian product of (umi_ld, lb_hd_relative)
        for umi_ld in args.umi_ld:
            for lb_rel in args.lb_hd_relative:
                combo_tag = f"u_{umi_ld}_l_{lb_rel}"
            combo_output_dir = os.path.join(sample_output_dir, combo_tag)
            os.makedirs(combo_output_dir, exist_ok=True)
            logger.info(
                f"Running denoising/annotation with umi_ld={umi_ld}, "
                f"lb_hd_relative={lb_rel} -> output dir: {combo_output_dir}"
            )

            #########################################################
            #### Diagnostic plots (pre-denoising, per-parameter set)
            #########################################################
            logger.info("Generating diagnostic plots (pre-denoising)...")
            plot_barcode_length_distribution(
                results_df, combo_output_dir, title_suffix="", unedited_bc_len=UNEDITED_BC_LEN
            )
            plot_reads_cutoff_distribution(results_df, args.reads_cutoff, combo_output_dir)
            plot_cutoff_vs_num_umis(results_df, args.reads_cutoff, combo_output_dir)
            plot_cutoff_vs_fraction_reads_retained(results_df, args.reads_cutoff, combo_output_dir)

            #########################################################
            #### Step 3: Denoise lineage barcode and UMI
            #########################################################
            logger.info("#### Step 3: Denoising lineage barcode and UMI...")
            agg, mapping, stats = correct_lineage_and_umi(
                results_clean,
                umi_col="UMI",
                bc_col="lineage_bc",
                n_iter=args.denoise_iter,
                umi_ld=umi_ld,
                lb_hd_relative=lb_rel,
                logger=logger,
            )
            logger.info(f"Denoising stats (umi_ld={umi_ld}, lb_hd_relative={lb_rel}): {stats}")
            
            agg['bc_len'] = agg['lineage_bc_corr'].str.len()
            agg.sort_values(by='n_reads', ascending=False, inplace=True)
            
            # Diagnostic plot 5: Distribution of lineage barcode lengths (after denoising)
            plot_barcode_length_by_umi(agg, combo_output_dir, unedited_bc_len=UNEDITED_BC_LEN)
            
            agg2 = agg.groupby('lineage_bc_corr').size().reset_index(name="UMIs")
            agg2.sort_values(by='UMIs', ascending=False, inplace=True)
            agg2['bc_len'] = agg2['lineage_bc_corr'].str.len()
            logger.info(f"Unique lineage barcodes (umi_ld={umi_ld}, lb_hd_relative={lb_rel}): {len(agg2)}")
            
            #########################################################
            #### Step 4: Annotate alleles
            #########################################################
            logger.info("#### Step 4: Annotating alleles...")
            sequences = agg2['lineage_bc_corr'].tolist()
            # For bulk DNA/RNA, we should perform reverse_complement on the sequences
            sequences = [str(Seq(s).reverse_complement()) for s in sequences]
            agg2['query'] = sequences
            
            results_allele = analyze_sequences(
                sequences, config=args.locus,
                min_sequence_length=args.min_bc_len, verbose=False
            )
            results_allele = results_allele.to_df()

            #########################################################
            #### Step 5: Merge and process final results
            #########################################################
            logger.info("#### Step 5: Processing final results...")
            final = agg2.merge(results_allele, on='query', how='left')
            final = final[['query', 'bc_len', 'UMIs', 'mutations', 'confidence', 'aligned_query', 'aligned_ref']]
            
            # Calculate MD5
            final['md5'] = final.apply(concat_and_md5, axis=1)
            
            # Group by MD5
            final2 = final.drop('query', axis=1).groupby('md5', as_index=False).agg({
                'bc_len': 'first',
                'UMIs': 'sum',
                'mutations': 'first',
                'confidence': 'first',
                'aligned_query': 'first',
                'aligned_ref': 'first'
            }).sort_values(by='UMIs', ascending=False).reset_index(drop=True)
            
            # Diagnostic plot 7: Clone size distribution
            plot_clone_size_distribution(final, combo_output_dir)
            
            # Diagnostic plot 8: Distribution of lineage barcode lengths by editing events
            plot_barcode_length_by_editing_events(final, combo_output_dir, unedited_bc_len=UNEDITED_BC_LEN)
            
            # Save output for this parameter combination
            output_file = os.path.join(combo_output_dir, f'{args.sample_id}_alleles.csv')
            final2.to_csv(output_file, index=False)
            logger.info(f"[{combo_tag}] Results saved to: {output_file}")
            logger.info(f"[{combo_tag}] Total unique alleles: {len(final2)}")
            logger.info(f"[{combo_tag}] Diagnostic plots saved to: {combo_output_dir}")
        
        # Clean up PEAR output directory if requested
        if not args.keep_pear:
            logger.info("Cleaning up PEAR output directory...")
            if os.path.exists(pear_output_dir):
                shutil.rmtree(pear_output_dir)
                logger.info(f"PEAR output directory removed: {pear_output_dir}")
        
        logger.info("Processing completed successfully!\n\n")
        
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()

