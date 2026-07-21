#!/usr/bin/env python3
"""
Pipeline: For each genome, identify unknown_likely_coding genes that match protein_coding genes.

Steps per genome:
  1) Pick longest protein per protein_coding gene from consensus protein FASTA
  2) Build GFF3 of unknown_likely_coding genes (excluding those overlapping protein_coding)
  3) Extract unknown_likely_coding cDNA with gffread
  4) makeblastdb on rep proteins, blastx unknown cDNA vs rep proteins
  5) Generate report: all matches, top match per gene, merge nearby clusters, rank by confidence

Usage:
  python scripts/unknown_coding_blastx_pipeline.py \
    --consensus-dir panprimate_output/consensus_gene_set \
    --genome-dir panprimate_output/genome_files \
    --work-dir panprimate_output/unknown_coding_tblastn \
    [--genome PPG00239~Pithecia_pithecia.pri~verkko] \
    [--threads 64] [--evalue 1e-5]
"""

import argparse
import csv
import os
import subprocess
import sys
from collections import defaultdict


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--consensus-dir", default="panprimate_output/consensus_gene_set")
    ap.add_argument("--genome-dir", default="panprimate_output/genome_files")
    ap.add_argument("--work-dir", default="panprimate_output/unknown_coding_tblastn")
    ap.add_argument("--genome", default=None, help="Process only this genome (e.g. PPG00239~Pithecia_pithecia.pri~verkko)")
    ap.add_argument("--threads", default=4, type=int)
    ap.add_argument("--evalue", default=1e-5, type=float)
    ap.add_argument("--merge-distance", default=0, type=int, help="Max gap (bp) to merge nearby unknown genes into clusters (0 = overlap only)")
    return ap.parse_args()


def discover_genomes(consensus_dir, genome_dir, single_genome=None):
    """Find all genomes with required files."""
    genomes = []
    if single_genome:
        bases = [single_genome]
    else:
        bases = []
        for f in os.listdir(consensus_dir):
            if f.endswith("_consensus.gff3"):
                bases.append(f.replace("_consensus.gff3", ""))

    for base in sorted(bases):
        paths = {
            "gff3": os.path.join(consensus_dir, f"{base}_consensus.gff3"),
            "gp_info": os.path.join(consensus_dir, f"{base}_consensus.gp_info"),
            "gp": os.path.join(consensus_dir, f"{base}_consensus.gp"),
            "protein_fa": os.path.join(consensus_dir, f"{base}_consensus_protein.fasta"),
            "genome_fa": os.path.join(genome_dir, f"{base}.fa"),
        }
        missing = [k for k, v in paths.items() if not os.path.isfile(v)]
        if missing:
            print(f"[{base}] Skipping: missing {', '.join(missing)}", file=sys.stderr)
            continue
        genomes.append((base, paths))
    return genomes


def run_cmd(cmd, desc=""):
    """Run a shell command, exit on failure."""
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  FAILED: {desc}\n  {r.stderr}", file=sys.stderr)
        return False
    return True


def step1_rep_proteins(paths, subdir):
    """Pick longest protein per protein_coding gene from consensus protein FASTA."""
    out_fa = os.path.join(subdir, "rep_proteins.fa")
    gp_info = paths["gp_info"]
    protein_fa = paths["protein_fa"]

    # Read gp_info
    tx_to_gene = {}
    gene_biotype = {}
    with open(gp_info) as f:
        header = next(f).strip().split("\t")
        gid_ix = header.index("gene_id")
        tid_ix = header.index("transcript_id")
        gb_ix = header.index("gene_biotype")
        for line in f:
            row = line.strip().split("\t")
            tx_to_gene[row[tid_ix]] = row[gid_ix]
            gene_biotype[row[gid_ix]] = row[gb_ix]

    pc_genes = {g for g, b in gene_biotype.items() if b == "protein_coding"}

    # Read all proteins
    proteins = {}
    cur_id, cur_seq = None, []
    with open(protein_fa) as f:
        for line in f:
            if line.startswith(">"):
                if cur_id:
                    proteins[cur_id] = "".join(cur_seq)
                cur_id = line[1:].strip().split()[0]
                cur_seq = []
            else:
                cur_seq.append(line.strip())
        if cur_id:
            proteins[cur_id] = "".join(cur_seq)

    # Longest protein per pc gene
    gene_best = {}
    for tid, seq in proteins.items():
        gid = tx_to_gene.get(tid)
        if not gid or gid not in pc_genes:
            continue
        seqlen = len(seq.rstrip("*"))
        if gid not in gene_best or seqlen > len(gene_best[gid][1].rstrip("*")):
            gene_best[gid] = (tid, seq)

    with open(out_fa, "w") as out:
        for gid in sorted(gene_best):
            tid, seq = gene_best[gid]
            out.write(f">{tid}\n{seq}\n")

    print(f"  Step 1: {len(pc_genes)} protein_coding genes, {len(gene_best)} have proteins")
    return out_fa, tx_to_gene, gene_biotype, pc_genes


def step2_unknown_gff3(paths, subdir, tx_to_gene, gene_biotype):
    """Filter GFF3 to unknown_likely_coding genes not overlapping protein_coding."""
    gp_info = paths["gp_info"]
    gp_path = paths["gp"]
    gff3 = paths["gff3"]
    out_gff3 = os.path.join(subdir, "unknown_likely_coding.gff3")

    # Read gp_info header indices
    with open(gp_info) as f:
        header = next(f).strip().split("\t")
        gid_ix = header.index("gene_id")
        tid_ix = header.index("transcript_id")

    # Gene coords from GP
    tx_coords = {}
    with open(gp_path) as f:
        for line in f:
            row = line.strip().split("\t")
            if len(row) < 12:
                continue
            tx_coords[row[0]] = (row[1], int(row[3]), int(row[4]))

    gene_coords = defaultdict(lambda: {"chrom": None, "start": None, "end": None})
    for tid, (chrom, s, e) in tx_coords.items():
        gid = tx_to_gene.get(tid)
        if not gid:
            continue
        g = gene_coords[gid]
        if g["chrom"] is None:
            g["chrom"], g["start"], g["end"] = chrom, s, e
        elif g["chrom"] == chrom:
            g["start"] = min(g["start"], s)
            g["end"] = max(g["end"], e)

    # PC intervals
    pc_intervals = []
    for gid, gb in gene_biotype.items():
        if gb != "protein_coding":
            continue
        g = gene_coords[gid]
        if g["chrom"] is not None:
            pc_intervals.append((g["chrom"], g["start"], g["end"]))

    def overlaps_any_pc(chrom, start, end):
        for pc_chrom, pc_start, pc_end in pc_intervals:
            if pc_chrom == chrom and not (end < pc_start or pc_end < start):
                return True
        return False

    # Find non-overlapping unknown genes
    unk_keep = set()
    unk_skip = 0
    for gid, gb in gene_biotype.items():
        if gb != "unknown_likely_coding":
            continue
        g = gene_coords[gid]
        if g["chrom"] is None or overlaps_any_pc(g["chrom"], g["start"], g["end"]):
            unk_skip += 1
        else:
            unk_keep.add(gid)

    # Get transcript IDs for kept genes
    tx_keep = set()
    with open(gp_info) as f:
        next(f)
        for line in f:
            row = line.strip().split("\t")
            if row[gid_ix] in unk_keep:
                tx_keep.add(row[tid_ix])

    # Filter GFF3
    written = 0
    with open(out_gff3, "w") as out:
        out.write("##gff-version 3\n")
        with open(gff3) as inp:
            for line in inp:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                oid, parent = "", ""
                for a in parts[8].split(";"):
                    if a.startswith("ID="):
                        oid = a[3:].strip()
                    elif a.startswith("Parent="):
                        parent = a[7:].strip()
                typ = parts[2]
                if (typ == "gene" and oid in unk_keep) or \
                   (typ in ("transcript", "mRNA") and oid in tx_keep) or \
                   (typ in ("exon", "CDS", "intron") and parent in tx_keep):
                    out.write(line)
                    written += 1

    print(f"  Step 2: {len(unk_keep)} unknown_likely_coding genes kept, {unk_skip} excluded (overlap pc), {written} GFF3 lines")
    return out_gff3, gene_coords, unk_keep


def step3_blastx(paths, subdir, rep_proteins_fa, unk_gff3, threads, evalue):
    """gffread cDNA extraction, makeblastdb, blastx."""
    genome_fa = paths["genome_fa"]
    cdna_fa = os.path.join(subdir, "unknown_likely_coding_cdna.fa")
    db_prefix = os.path.join(subdir, "rep_proteins_db")
    blastx_out = os.path.join(subdir, "blastx.out")

    # gffread -w
    if not run_cmd(["gffread", "-g", genome_fa, "-w", cdna_fa, unk_gff3], "gffread"):
        return None
    print(f"  Step 3a: Extracted cDNA: {cdna_fa}")

    # makeblastdb
    if not os.path.isfile(db_prefix + ".phr"):
        if not run_cmd(["makeblastdb", "-in", rep_proteins_fa, "-dbtype", "prot", "-out", db_prefix], "makeblastdb"):
            return None
    print(f"  Step 3b: Protein DB ready")

    # blastx
    if not run_cmd([
        "blastx", "-query", cdna_fa, "-db", db_prefix,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-evalue", str(evalue), "-num_threads", str(threads), "-out", blastx_out
    ], "blastx"):
        return None
    n_hits = sum(1 for _ in open(blastx_out))
    print(f"  Step 3c: blastx done: {n_hits} hits")
    return blastx_out


def step4_report(paths, subdir, blastx_out, tx_to_gene, gene_biotype, unk_keep, gene_coords):
    """Generate all-matches, top-match-per-gene, and no-match reports."""
    gp_info = paths["gp_info"]
    all_match_tsv = os.path.join(subdir, "unknown_likely_coding_best_pc_match.tsv")
    top_match_tsv = os.path.join(subdir, "unknown_likely_coding_top_match.tsv")
    no_match_tsv = os.path.join(subdir, "unknown_likely_coding_no_match.tsv")

    # Gene common names
    gene_common_name = {}
    with open(gp_info) as f:
        header = next(f).strip().split("\t")
        gid_ix = header.index("gene_id")
        sgcn_ix = header.index("source_gene_common_name")
        for line in f:
            row = line.strip().split("\t")
            gid = row[gid_ix]
            name = row[sgcn_ix] if sgcn_ix < len(row) else "N/A"
            if gid not in gene_common_name or gene_common_name[gid] in ("N/A", "None", ""):
                gene_common_name[gid] = name

    # Parse blastx
    unk_gene_hits = defaultdict(list)
    with open(blastx_out) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            q_tx, s_tx, pident = parts[0], parts[1], float(parts[2])
            ev, bs = float(parts[10]), float(parts[11])
            unk_gene = tx_to_gene.get(q_tx)
            pc_gene = tx_to_gene.get(s_tx)
            if not unk_gene or not pc_gene:
                continue
            if gene_biotype.get(unk_gene) != "unknown_likely_coding":
                continue
            if gene_biotype.get(pc_gene) != "protein_coding":
                continue
            unk_gene_hits[unk_gene].append({
                "pc_gene": pc_gene,
                "pc_name": gene_common_name.get(pc_gene, "N/A"),
                "bitscore": bs, "evalue": ev, "pident": pident,
            })

    # All matches (best per pc gene per unknown gene)
    with open(all_match_tsv, "w") as out:
        out.write("unknown_likely_coding_gene\tmatching_protein_coding_gene\tpc_gene_name\tpident\tbitscore\tevalue\tnum_distinct_pc_genes\n")
        for unk_gid in sorted(unk_gene_hits):
            hits = unk_gene_hits[unk_gid]
            best_per_pc = {}
            for h in hits:
                pg = h["pc_gene"]
                if pg not in best_per_pc or h["bitscore"] > best_per_pc[pg]["bitscore"]:
                    best_per_pc[pg] = h
            n_distinct = len(best_per_pc)
            for h in sorted(best_per_pc.values(), key=lambda x: (-x["bitscore"], x["evalue"])):
                out.write(f"{unk_gid}\t{h['pc_gene']}\t{h['pc_name']}\t{h['pident']:.1f}\t{h['bitscore']:.1f}\t{h['evalue']:.2e}\t{n_distinct}\n")

    # Top match per unknown gene
    top_rows = []
    with open(all_match_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        seen = set()
        for row in reader:
            gid = row["unknown_likely_coding_gene"]
            if gid not in seen:
                seen.add(gid)
                top_rows.append(row)

    with open(top_match_tsv, "w") as out:
        out.write("unknown_likely_coding_gene\tbest_matching_pc_gene\tpc_gene_name\tpident\tbitscore\tevalue\tnum_distinct_pc_genes\n")
        for row in top_rows:
            out.write(f"{row['unknown_likely_coding_gene']}\t{row['matching_protein_coding_gene']}\t{row['pc_gene_name']}\t{row['pident']}\t{row['bitscore']}\t{row['evalue']}\t{row['num_distinct_pc_genes']}\n")

    # No-match genes: in unk_keep but not in unk_gene_hits
    no_match_genes = sorted(unk_keep - set(unk_gene_hits.keys()))
    with open(no_match_tsv, "w") as out:
        out.write("unknown_likely_coding_gene\tchrom\tstart\tend\n")
        for gid in no_match_genes:
            g = gene_coords.get(gid)
            if g and g["chrom"]:
                out.write(f"{gid}\t{g['chrom']}\t{g['start']}\t{g['end']}\n")
            else:
                out.write(f"{gid}\tNA\t0\t0\n")

    print(f"  Step 4: {len(unk_gene_hits)} unknown genes with hits, {len(no_match_genes)} with no match, {len(top_rows)} top matches")
    return top_match_tsv


def step5_merge_and_rank(paths, subdir, top_match_tsv, tx_to_gene, gene_coords, merge_distance, unk_keep):
    """Merge overlapping unknown genes into clusters (both with and without hits) and rank."""
    merged_tsv = os.path.join(subdir, "unknown_likely_coding_top_match_merged.tsv")
    ranked_tsv = os.path.join(subdir, "unknown_likely_coding_top_match_merged_ranked.tsv")

    # Read top matches (genes WITH hits)
    matched_genes = {}  # gid -> row dict
    with open(top_match_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            matched_genes[row["unknown_likely_coding_gene"]] = row

    # Build unified list: all unknown genes (with and without hits)
    all_genes = []
    for gid in sorted(unk_keep):
        g = gene_coords.get(gid)
        chrom = g["chrom"] if g and g["chrom"] else "NA"
        start = g["start"] if g and g["start"] is not None else 0
        end = g["end"] if g and g["end"] is not None else 0
        if gid in matched_genes:
            row = matched_genes[gid]
            all_genes.append({
                "gene_id": gid, "chrom": chrom, "start": start, "end": end,
                "has_match": True,
                "best_matching_pc_gene": row["best_matching_pc_gene"],
                "pc_gene_name": row["pc_gene_name"],
                "pident": row["pident"], "bitscore": row["bitscore"],
                "evalue": row["evalue"], "num_distinct_pc_genes": row["num_distinct_pc_genes"],
            })
        else:
            all_genes.append({
                "gene_id": gid, "chrom": chrom, "start": start, "end": end,
                "has_match": False,
                "best_matching_pc_gene": "N/A", "pc_gene_name": "N/A",
                "pident": "0", "bitscore": "0", "evalue": "N/A",
                "num_distinct_pc_genes": "0",
            })

    all_genes.sort(key=lambda x: (x["chrom"], x["start"]))

    # Cluster by overlap
    clusters = []
    for g in all_genes:
        if not clusters:
            clusters.append([g])
            continue
        last_cluster = clusters[-1]
        cluster_end = max(d["end"] for d in last_cluster)
        cluster_chrom = last_cluster[0]["chrom"]
        if cluster_chrom == g["chrom"] and g["start"] <= cluster_end + merge_distance:
            last_cluster.append(g)
        else:
            clusters.append([g])

    # Write merged
    fieldnames = ["cluster_id", "unknown_likely_coding_genes", "chrom", "cluster_start", "cluster_end",
                  "num_genes_in_cluster", "best_matching_pc_gene", "pc_gene_name",
                  "pident", "bitscore", "evalue", "num_distinct_pc_genes"]
    with open(merged_tsv, "w") as out:
        out.write("\t".join(fieldnames) + "\n")
        for i, cluster in enumerate(clusters):
            gene_ids = ",".join(d["gene_id"] for d in cluster)
            chrom = cluster[0]["chrom"]
            cstart = min(d["start"] for d in cluster)
            cend = max(d["end"] for d in cluster)
            matched_in_cluster = [d for d in cluster if d["has_match"]]
            if matched_in_cluster:
                best = max(matched_in_cluster, key=lambda d: float(d["bitscore"]))
            else:
                best = cluster[0]  # no match; use defaults (N/A)
            out.write(f"cluster_{i+1}\t{gene_ids}\t{chrom}\t{cstart}\t{cend}\t{len(cluster)}\t{best['best_matching_pc_gene']}\t{best['pc_gene_name']}\t{best['pident']}\t{best['bitscore']}\t{best['evalue']}\t{best['num_distinct_pc_genes']}\n")

    # Rank by confidence
    rows = []
    with open(merged_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            try:
                pident = float(row["pident"])
                bitscore = float(row["bitscore"])
            except ValueError:
                pident, bitscore = 0, 0
            if row["best_matching_pc_gene"] == "N/A":
                tier, tier_rank = "no_match", 3
            elif pident >= 70 and bitscore >= 200:
                tier, tier_rank = "high", 0
            elif pident >= 50 and bitscore >= 100:
                tier, tier_rank = "medium", 1
            else:
                tier, tier_rank = "low", 2
            row["confidence"] = tier
            row["_rank"] = (tier_rank, -bitscore)
            rows.append(row)

    rows.sort(key=lambda r: r["_rank"])

    with open(ranked_tsv, "w") as out:
        out.write("\t".join(fieldnames + ["confidence"]) + "\n")
        for i, row in enumerate(rows):
            row["cluster_id"] = f"cluster_{i+1}"
            out.write("\t".join(row[f] for f in fieldnames + ["confidence"]) + "\n")

    n_merged = sum(1 for c in clusters if len(c) > 1)
    named = sum(1 for r in rows if r["confidence"] != "no_match")
    unnamed = sum(1 for r in rows if r["confidence"] == "no_match")
    high = sum(1 for r in rows if r["confidence"] == "high")
    med = sum(1 for r in rows if r["confidence"] == "medium")
    low = sum(1 for r in rows if r["confidence"] == "low")
    print(f"  Step 5: {len(clusters)} total clusters ({n_merged} with overlapping genes merged)")
    print(f"    Named (has pc match):   {named}  (high={high}, medium={med}, low={low})")
    print(f"    Unnamed (no pc match):  {unnamed}")
    return ranked_tsv


def process_genome(base, paths, work_dir, threads, evalue, merge_distance):
    """Run full pipeline for one genome."""
    print(f"\n{'='*60}")
    print(f"Processing: {base}")
    print(f"{'='*60}")

    subdir = os.path.join(work_dir, base)
    os.makedirs(subdir, exist_ok=True)

    rep_fa, tx_to_gene, gene_biotype, pc_genes = step1_rep_proteins(paths, subdir)
    unk_gff3, gene_coords, unk_keep = step2_unknown_gff3(paths, subdir, tx_to_gene, gene_biotype)
    blastx_out = step3_blastx(paths, subdir, rep_fa, unk_gff3, threads, evalue)
    if not blastx_out:
        print(f"  FAILED at blastx step", file=sys.stderr)
        return
    top_match = step4_report(paths, subdir, blastx_out, tx_to_gene, gene_biotype, unk_keep, gene_coords)
    ranked = step5_merge_and_rank(paths, subdir, top_match, tx_to_gene, gene_coords, merge_distance, unk_keep)
    print(f"  DONE: {ranked}")


def main():
    args = parse_args()
    consensus_dir = os.path.abspath(args.consensus_dir)
    genome_dir = os.path.abspath(args.genome_dir)
    work_dir = os.path.abspath(args.work_dir)
    os.makedirs(work_dir, exist_ok=True)

    genomes = discover_genomes(consensus_dir, genome_dir, args.genome)
    if not genomes:
        print("No genomes found.", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(genomes)} genome(s) to process")
    for base, paths in genomes:
        process_genome(base, paths, work_dir, args.threads, args.evalue, args.merge_distance)

    print(f"\n{'='*60}")
    print(f"ALL DONE")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
