#!/usr/bin/env python3
"""
Pipeline: extract one representative protein per protein_coding gene from consensus GFF3
using gffread (assembly + GFF3), tblastn vs genome, then report unknown_likely_coding
genes that get hits not overlapping protein_coding but overlapping themselves, with best
matching protein_coding gene.

Usage:
  python scripts/unknown_coding_tblastn_pipeline.py \
    --consensus-dir panprimate_output/consensus_gene_set \
    --genome-dir panprimate_output/genome_files \
    --work-dir panprimate_output/unknown_coding_tblastn \
    [--genome <PPG...>]  # optional: run for one genome only
"""

import argparse
import os
import re
import subprocess
import sys
from collections import defaultdict

def parse_attrs(attr_str):
    d = {}
    for part in attr_str.split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            d[k.strip()] = v.strip()
    return d

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--consensus-dir", default="panprimate_output/consensus_gene_set", help="Directory with *_consensus.gff3 and *_consensus.gp_info")
    ap.add_argument("--genome-dir", default="panprimate_output/genome_files", help="Directory with PPG*~*.pri~verkko.fa genomes")
    ap.add_argument("--work-dir", default="panprimate_output/unknown_coding_tblastn", help="Working directory for filtered GFF3, proteins, BLAST DB, results")
    ap.add_argument("--genome", default=None, help="Optional: process only this genome base (e.g. PPG00525~Lagothrix_lagotricha.pri~verkko)")
    ap.add_argument("--gffread", default="gffread", help="Path to gffread")
    ap.add_argument("--makeblastdb", default="makeblastdb", help="Path to makeblastdb")
    ap.add_argument("--tblastn", default="tblastn", help="Path to tblastn")
    ap.add_argument("--evalue", default=1e-5, type=float, help="tblastn E-value threshold")
    ap.add_argument("--threads", default=4, type=int, help="BLAST threads")
    ap.add_argument("--stop-after", type=int, default=None, choices=[1, 2, 3, 4, 5], help="Stop after step N (1=filter GFF3, 2=gffread, 3=makeblastdb, 4=tblastn, 5=report). Default: run all.")
    ap.add_argument("--skip-blast", action="store_true", help="Skip makeblastdb and tblastn; use existing tblastn.out (run from step 6).")
    args = ap.parse_args()

    consensus_dir = os.path.abspath(args.consensus_dir)
    genome_dir = os.path.abspath(args.genome_dir)
    work_dir = os.path.abspath(args.work_dir)
    os.makedirs(work_dir, exist_ok=True)

    # Discover GFF3 files
    if args.genome:
        base = args.genome.rstrip("/")
        if not base.endswith("_consensus"):
            base_consensus = base + "_consensus"
        else:
            base_consensus = base
        gff3_path = os.path.join(consensus_dir, base_consensus + ".gff3")
        if not os.path.isfile(gff3_path):
            # try with _consensus already in name
            gff3_path = os.path.join(consensus_dir, base + ".gff3")
        if not os.path.isfile(gff3_path):
            print(f"Error: GFF3 not found: {gff3_path}", file=sys.stderr)
            sys.exit(1)
        gff3_files = [(base.replace("_consensus", "").rstrip("."), gff3_path)]
    else:
        gff3_files = []
        for f in os.listdir(consensus_dir):
            if f.endswith("_consensus.gff3"):
                base = f.replace("_consensus.gff3", "")
                gff3_files.append((base, os.path.join(consensus_dir, f)))

    if not gff3_files:
        print("No consensus GFF3 files found.", file=sys.stderr)
        sys.exit(1)

    for genome_base, gff3_path in sorted(gff3_files):
        gp_info_path = gff3_path.replace("_consensus.gff3", "_consensus.gp_info")
        gp_path = gff3_path.replace("_consensus.gff3", "_consensus.gp")
        genome_fa = os.path.join(genome_dir, genome_base + ".fa")
        if not os.path.isfile(genome_fa):
            print(f"Skipping {genome_base}: genome FASTA not found: {genome_fa}", file=sys.stderr)
            continue
        if not os.path.isfile(gp_info_path):
            print(f"Skipping {genome_base}: gp_info not found: {gp_info_path}", file=sys.stderr)
            continue
        if not os.path.isfile(gp_path):
            print(f"Skipping {genome_base}: gp not found: {gp_path}", file=sys.stderr)
            continue

        run_one(
            genome_base=genome_base,
            gff3_path=gff3_path,
            gp_info_path=gp_info_path,
            gp_path=gp_path,
            genome_fa=genome_fa,
            work_dir=work_dir,
            gffread=args.gffread,
            makeblastdb=args.makeblastdb,
            tblastn=args.tblastn,
            evalue=args.evalue,
            threads=args.threads,
            stop_after=args.stop_after,
            skip_blast=args.skip_blast,
        )

def run_one(genome_base, gff3_path, gp_info_path, gp_path, genome_fa, work_dir, gffread, makeblastdb, tblastn, evalue, threads, stop_after=None, skip_blast=False):
    subdir = os.path.join(work_dir, genome_base)
    os.makedirs(subdir, exist_ok=True)
    filtered_gff3 = os.path.join(subdir, "rep_protein_coding.gff3")
    proteins_fa = os.path.join(subdir, "rep_proteins.fa")
    blast_db_prefix = os.path.join(subdir, "genome_db")
    blast_out = os.path.join(subdir, "tblastn.out")
    report_tsv = os.path.join(subdir, "unknown_likely_coding_best_pc_match.tsv")

    # 1) From gp_info: one representative transcript per protein_coding gene (first with source_transcript != 'N/A', i.e. has a real reference transcript ID)
    rep_tx = {}   # gene_id -> transcript_id
    gene_biotype = {}  # gene_id -> gene_biotype
    with open(gp_info_path) as f:
        header = next(f).strip().split("\t")
        try:
            gid_ix = header.index("gene_id")
            tid_ix = header.index("transcript_id")
            st_ix = header.index("source_transcript")
            gb_ix = header.index("gene_biotype")
        except ValueError as e:
            print(f"Missing column in gp_info: {e}", file=sys.stderr)
            return
        for line in f:
            row = line.strip().split("\t")
            if len(row) <= max(gid_ix, tid_ix, st_ix, gb_ix):
                continue
            gid = row[gid_ix]
            tid = row[tid_ix]
            st = row[st_ix] if st_ix < len(row) else "N/A"
            gb = row[gb_ix] if gb_ix < len(row) else ""
            gene_biotype[gid] = gb
            if gb != "protein_coding":
                continue
            if gid not in rep_tx and st and st != "N/A":
                rep_tx[gid] = tid

    rep_tx_ids = set(rep_tx.values())
    rep_gene_ids = set(rep_tx.keys())
    tx_to_pc_gene = {tid: gid for gid, tid in rep_tx.items()}
    print(f"[{genome_base}] Representative transcripts for {len(rep_tx)} protein_coding genes.")

    # 2) Filter GFF3: keep gene, transcript, exon, CDS for rep transcripts only
    def keep_line(line):
        if line.startswith("#"):
            return True
        parts = line.strip().split("\t")
        if len(parts) < 9:
            return False
        typ = parts[2]
        attrs = parse_attrs(parts[8])
        oid = attrs.get("ID", "")
        parent = attrs.get("Parent", "")
        if typ == "gene":
            return oid in rep_gene_ids
        if typ in ("transcript", "mRNA"):
            return oid in rep_tx_ids
        if typ in ("exon", "CDS"):
            return parent in rep_tx_ids
        return False

    with open(filtered_gff3, "w") as out:
        with open(gff3_path) as inp:
            for line in inp:
                if keep_line(line):
                    out.write(line)
    print(f"[{genome_base}] Wrote filtered GFF3: {filtered_gff3}")
    if stop_after == 1:
        return

    # 3) gffread -y from assembly
    cmd = [gffread, "-g", genome_fa, "-y", proteins_fa, filtered_gff3]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"gffread failed: {r.stderr}", file=sys.stderr)
        return
    if not os.path.isfile(proteins_fa) or os.path.getsize(proteins_fa) == 0:
        print(f"No proteins extracted. Check CDS in filtered GFF3.", file=sys.stderr)
        return
    print(f"[{genome_base}] Extracted proteins: {proteins_fa}")
    if stop_after == 2:
        return

    # 4) makeblastdb on genome
    db_path = blast_db_prefix + ".nhr"
    if not os.path.isfile(db_path):
        cmd = [makeblastdb, "-in", genome_fa, "-dbtype", "nucl", "-out", blast_db_prefix]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"makeblastdb failed: {r.stderr}", file=sys.stderr)
            return
    print(f"[{genome_base}] BLAST DB ready.")
    if stop_after == 3:
        return

    # 5) tblastn
    if not skip_blast:
        cmd = [tblastn, "-query", proteins_fa, "-db", blast_db_prefix, "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore", "-evalue", str(evalue), "-num_threads", str(threads), "-out", blast_out]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            print(f"tblastn failed: {r.stderr}", file=sys.stderr)
            return
        print(f"[{genome_base}] tblastn done: {blast_out}")
    elif not os.path.isfile(blast_out):
        print(f"Error: --skip-blast but {blast_out} not found", file=sys.stderr)
        return
    if stop_after == 4:
        return

    # 6) Gene coordinates from GP: transcript_id -> (chrom, start, end), then gene_id -> (chrom, start, end)
    tx_coords = {}
    with open(gp_path) as f:
        for line in f:
            row = line.strip().split("\t")
            if len(row) < 12:
                continue
            tid, chrom, _, tx_start, tx_end = row[0], row[1], row[2], int(row[3]), int(row[4])
            tx_coords[tid] = (chrom, tx_start, tx_end)
    gene_coords = defaultdict(lambda: {"chrom": None, "start": None, "end": None})
    with open(gp_info_path) as f:
        next(f)  # header
        for line in f:
            row = line.strip().split("\t")
            if len(row) <= gid_ix:
                continue
            gid = row[gid_ix]
            tid = row[tid_ix]
            if tid not in tx_coords:
                continue
            chrom, s, e = tx_coords[tid]
            g = gene_coords[gid]
            if g["chrom"] is None:
                g["chrom"], g["start"], g["end"] = chrom, s, e
            else:
                if g["chrom"] != chrom:
                    continue
                g["start"] = min(g["start"], s)
                g["end"] = max(g["end"], e)

    def overlap(a_start, a_end, b_start, b_end):
        """All coordinates 0-based inclusive."""
        return not (a_end < b_start or b_end < a_start)

    # Gene coords from GP are 0-based start, 1-based end; convert to 0-based inclusive for overlap
    pc_genes = [gid for gid, gb in gene_biotype.items() if gb == "protein_coding"]
    unk_genes = [gid for gid, gb in gene_biotype.items() if gb == "unknown_likely_coding"]
    pc_intervals = {}
    for gid in pc_genes:
        g = gene_coords[gid]
        if g["chrom"] is None:
            continue
        pc_intervals[gid] = (g["chrom"], g["start"], g["end"] - 1 if g["end"] is not None else g["start"])
    unk_intervals = {}
    for gid in unk_genes:
        g = gene_coords[gid]
        if g["chrom"] is None:
            continue
        unk_intervals[gid] = (g["chrom"], g["start"], g["end"] - 1 if g["end"] is not None else g["start"])

    # 7) For each tblastn hit: if overlaps unknown_likely_coding and not protein_coding, record (unk_gene, pc_gene, location, bitscore).
    # Keep all matches so we can flag unknown genes that match multiple genes or multiple locations.
    unk_to_matches = defaultdict(list)  # unk_gid -> list of {pc_gene, bitscore, evalue, hit_chrom, hit_start, hit_end}
    with open(blast_out) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, ev, bitscore = parts[:12]
            pc_gene = tx_to_pc_gene.get(qseqid)
            if not pc_gene:
                continue
            try:
                s1, s2 = int(sstart), int(send)
                hit_start = min(s1, s2) - 1
                hit_end = max(s1, s2) - 1
                hit_chrom = sseqid
            except ValueError:
                continue
            for unk_gid, (uchrom, ustart, uend) in unk_intervals.items():
                if uchrom != hit_chrom or not overlap(hit_start, hit_end, ustart, uend):
                    continue
                overlaps_pc = any(
                    pchrom == hit_chrom and overlap(hit_start, hit_end, pstart, pend)
                    for pgid, (pchrom, pstart, pend) in pc_intervals.items()
                )
                if overlaps_pc:
                    continue
                try:
                    bs = float(bitscore)
                    ev = float(ev)
                except ValueError:
                    continue
                unk_to_matches[unk_gid].append({
                    "pc_gene": pc_gene,
                    "bitscore": bs,
                    "evalue": ev,
                    "hit_chrom": hit_chrom,
                    "hit_start": hit_start,
                    "hit_end": hit_end,
                })

    # 8) Write report: one row per (unknown_gene, pc_gene, location); add multi_locus when unknown matches >1 gene or >1 location
    with open(report_tsv, "w") as out:
        out.write("unknown_likely_coding_gene\tmatching_protein_coding_gene\thit_chrom\thit_start\thit_end\tbitscore\tevalue\tmulti_locus\n")
        for unk_gid in sorted(unk_to_matches.keys()):
            matches = unk_to_matches[unk_gid]
            distinct_pc = len(set(m["pc_gene"] for m in matches))
            distinct_loci = len(set((m["hit_chrom"], m["hit_start"], m["hit_end"]) for m in matches))
            multi_locus = distinct_pc > 1 or distinct_loci > 1
            for m in sorted(matches, key=lambda x: (-x["bitscore"], x["evalue"])):
                loc_str = "yes" if multi_locus else "no"
                out.write(f"{unk_gid}\t{m['pc_gene']}\t{m['hit_chrom']}\t{m['hit_start']}\t{m['hit_end']}\t{m['bitscore']}\t{m['evalue']}\t{loc_str}\n")
    n_unk_with_match = sum(1 for u, ms in unk_to_matches.items() if ms)
    n_multi = sum(1 for u, ms in unk_to_matches.items() if ms and (len(set(m["pc_gene"] for m in ms)) > 1 or len(set((m["hit_chrom"], m["hit_start"], m["hit_end"]) for m in ms)) > 1))
    print(f"[{genome_base}] Report: {report_tsv} ({n_unk_with_match} unknown_likely_coding genes with ≥1 match; {n_multi} match multiple genes/locations).")

if __name__ == "__main__":
    main()
