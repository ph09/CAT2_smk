#!/usr/bin/env python3
"""
Gene family / ortholog / paralog analysis for a CAT2 run.

CAT2 annotates every target genome by *projecting* the reference annotation, so
orthology is already encoded in the output: every consensus transcript carries the
reference gene it came from (``source_gene`` in ``*_consensus_novel_annotated.gp_info``). This lets
us reconstruct orthologs and paralogs directly from the final gene sets without any
extra sequence clustering:

  * ORTHOLOGS  : consensus loci in different genomes that share the same ``source_gene``.
  * PARALOGS   : two or more distinct consensus loci (``gene_id``) in the *same* genome
                 that share the same ``source_gene`` (i.e. the reference gene was found
                 in >1 place -> a duplication / gene-copy-number gain).
  * COPY NUMBER: number of distinct consensus loci per ``source_gene`` per genome.

From copy number we call, per reference gene and per (heuristic) gene family:

  * lost        : 0 copies in a genome (contraction to zero).
  * single-copy : exactly 1 copy (conserved).
  * expanded    : more copies than the reference has (duplication).
  * contracted  : (family level) fewer total copies than the reference family size.

The analysis is restricted to protein-coding reference genes by default
(``source_gene_biotype == 'protein_coding'``), which is what you almost always want
for expansion/contraction work.

Two grouping levels are produced:

  1. ORTHOLOG GROUPS  -- one per reference gene (``source_gene``). Rigorous: this is
     exactly CAT2's projection-based orthology.
  2. GENE FAMILIES    -- ortholog groups bucketed by a heuristic symbol root (e.g.
     ``KLHL4`` -> ``KLHL``, ``OR4A1`` -> ``OR``). This is a convenience aggregation to
     surface tandem-array / multigene families; genes with no informative symbol
     (``LOC*`` / unnamed) stay as singleton families. Tune with --family-regex.

Inputs are read straight from a finished (or partially finished) work_dir:
  {work_dir}/databases/{ref_genome}.db          reference annotation table
  {work_dir}/consensus_gene_set/{genome}_consensus_novel_annotated.gp_info
  {work_dir}/consensus_gene_set/{genome}_consensus_novel_annotated.gp

Outputs (into --out-dir, default {work_dir}/gene_family_analysis/):
  ortholog_copy_matrix.tsv        genes x genomes, values = # loci (copy number)
  ortholog_transcript_matrix.tsv  genes x genomes, values = # transcripts
  paralogs.tsv                    one row per locus for multi-copy genes (the paralogs)
  gene_family_matrix.tsv          families x genomes, total copies + call
  expansion_contraction.tsv       long form per (gene|family, genome) with call
  novel_genes.tsv                 genome-specific / de-novo PC loci with no ortholog
  summary.md                      human-readable report (phylogenetically ordered)
  copy_number_heatmap.png         heatmap of the most copy-number-variable families
                                  (skipped if matplotlib is unavailable)
"""
import argparse
import csv
import os
import re
import sqlite3
import sys
from collections import defaultdict

# gp_info is small enough that stdlib csv is plenty; pandas only used for the
# matrix pivot / heatmap where it is genuinely convenient.
try:
    import pandas as pd
except Exception:  # pragma: no cover
    pd = None


# ─── config / genome discovery ────────────────────────────────────────────────
def load_yaml_config(path):
    try:
        import yaml
    except Exception:
        sys.exit("PyYAML is required to read a --config file; install pyyaml or pass "
                 "--ref-genome/--genomes explicitly.")
    with open(path) as fh:
        return yaml.safe_load(fh)


def discover_genomes(consensus_dir):
    """All genomes that have a *_consensus_novel_annotated.gp_info in the consensus dir."""
    suffix = "_consensus_novel_annotated.gp_info"
    out = []
    for fn in sorted(os.listdir(consensus_dir)):
        if fn.endswith(suffix):
            out.append(fn[: -len(suffix)])
    return out


def order_genomes(genomes, hal, ref_genome):
    """Order genomes by phylogenetic distance from ref if halStats works, else input order."""
    if not hal or not os.path.exists(hal):
        return list(genomes)
    try:
        # local import so the script still runs without ete3/halStats installed
        sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from tools.hal import build_genome_order
        # build_genome_order is @memoize'd, so every arg must be hashable -> frozenset.
        ordered = build_genome_order(hal, ref_genome, genome_subset=frozenset(genomes),
                                     include_ancestors=True)
        # build_genome_order drops any genome not in the HAL / not a descendant; keep
        # the rest appended in input order so nothing silently disappears.
        ordered = [g for g in ordered if g in genomes]
        ordered += [g for g in genomes if g not in ordered]
        return ordered
    except Exception as e:
        sys.stderr.write(f"[warn] could not order genomes via HAL ({e}); using input order\n")
        return list(genomes)


# ─── reference annotation ──────────────────────────────────────────────────────
def load_reference_pc(ref_db, biotype="protein_coding"):
    """
    Return dict keyed by reference gene id ->
        {name, biotype, n_ref_tx}
    restricted to the requested gene biotype.
    """
    if not os.path.exists(ref_db):
        sys.exit(f"reference db not found: {ref_db}")
    con = sqlite3.connect(ref_db)
    cur = con.execute(
        "SELECT GeneId, GeneName, GeneBiotype, TranscriptId "
        "FROM annotation WHERE GeneBiotype = ?", (biotype,))
    genes = {}
    for gene_id, gene_name, gbiotype, tx_id in cur:
        gid = str(gene_id)
        rec = genes.setdefault(gid, {"name": gene_name or gid,
                                     "biotype": gbiotype, "n_ref_tx": 0})
        rec["n_ref_tx"] += 1
    con.close()
    return genes


# ─── consensus parsing ─────────────────────────────────────────────────────────
def load_gp_coords(gp_path):
    """transcript_id -> (chrom, strand, txStart, txEnd) from a consensus genePred."""
    coords = {}
    if not os.path.exists(gp_path):
        return coords
    with open(gp_path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            # genePred: name, chrom, strand, txStart, txEnd, cdsStart, cdsEnd, ...
            coords[f[0]] = (f[1], f[2], int(f[3]), int(f[4]))
    return coords


def parse_consensus(gp_info_path, gp_path, select_biotype="protein_coding",
                    biotype_field="source_gene_biotype"):
    """
    Parse one genome's consensus gp_info.

    Returns:
      loci        : {source_gene: {gene_id: {"tx": set(tx_ids), "modes": set(),
                                             "classes": set()}}}
      novel       : list of dicts for PC loci with no reference source_gene
      gene_meta   : {source_gene: {"common_name": str}}
    coords are attached per-locus (min start / max end over its transcripts).
    """
    coords = load_gp_coords(gp_path)
    loci = defaultdict(lambda: defaultdict(lambda: {"tx": set(), "modes": set(),
                                                    "classes": set()}))
    gene_meta = {}
    novel = []
    with open(gp_info_path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            src_biotype = (row.get(biotype_field) or "").strip()
            gene_biotype = (row.get("gene_biotype") or "").strip()
            source_gene = (row.get("source_gene") or "").strip()
            gene_id = (row.get("gene_id") or "").strip()
            tx_id = (row.get("transcript_id") or "").strip()
            tclass = (row.get("transcript_class") or "").strip()
            mode = (row.get("alignment_mode") or "").strip()
            common = (row.get("source_gene_common_name") or "").strip()

            # Novel / de novo PC loci: no reference source_gene, but predicted coding.
            # This covers de-novo (augPB/strg) novel genes AND protein-only novel
            # genes (augMP, source_gene "MP-NOVEL-*", transcript_class putative_novel)
            # promoted by consensus so lineage-specific genes are reported as novel.
            is_novel = (not source_gene or src_biotype in ("", "N/A")
                        or tclass == "putative_novel"
                        or str(source_gene).startswith("MP-NOVEL-"))
            if is_novel:
                if gene_biotype == select_biotype:
                    c = coords.get(tx_id, ("NA", ".", 0, 0))
                    novel.append({"gene_id": gene_id, "transcript_id": tx_id,
                                  "chrom": c[0], "start": c[2], "end": c[3],
                                  "transcript_class": tclass, "alignment_mode": mode})
                continue

            # Track fate of reference protein-coding genes (by reference biotype).
            if src_biotype != select_biotype:
                continue

            rec = loci[source_gene][gene_id]
            rec["tx"].add(tx_id)
            rec["modes"].add(mode)
            rec["classes"].add(tclass)
            if source_gene not in gene_meta:
                gene_meta[source_gene] = {"common_name": common or source_gene}

    # attach coordinates per locus
    for source_gene, gmap in loci.items():
        for gene_id, rec in gmap.items():
            chroms, starts, ends, strands = [], [], [], []
            for tx_id in rec["tx"]:
                c = coords.get(tx_id)
                if c:
                    chroms.append(c[0]); strands.append(c[1])
                    starts.append(c[2]); ends.append(c[3])
            rec["chrom"] = chroms[0] if chroms else "NA"
            rec["strand"] = strands[0] if strands else "."
            rec["start"] = min(starts) if starts else 0
            rec["end"] = max(ends) if ends else 0
    return loci, novel, gene_meta


# ─── gene family heuristic ─────────────────────────────────────────────────────
DEFAULT_FAMILY_REGEX = r"^([A-Za-z][A-Za-z]+?)\d"
LOC_RE = re.compile(r"^(LOC)?\d+$", re.IGNORECASE)


def family_key(common_name, gene_id, family_re):
    """
    Heuristic gene-family bucket from a gene symbol.
      KLHL4 -> KLHL ; OR4A1 -> OR ; DOCK11 -> DOCK
    Unnamed / LOC##### / purely numeric ids get their own singleton family so we
    never merge unrelated unnamed genes.
    """
    name = (common_name or "").strip()
    if not name or LOC_RE.match(name):
        return f"__singleton__:{gene_id}"
    m = family_re.match(name)
    if m and len(m.group(1)) >= 2:
        return m.group(1).upper()
    # no trailing digits (e.g. "TP53BP1"): the whole symbol is its own family
    return name.upper()


# ─── classification ────────────────────────────────────────────────────────────
def call_gene(copies):
    if copies == 0:
        return "lost"
    if copies == 1:
        return "single_copy"
    return "expanded"


def call_family(obs, ref_size):
    if obs == 0:
        return "lost"
    if obs > ref_size:
        return "expanded"
    if obs < ref_size:
        return "contracted"
    return "conserved"


def is_significant(parent, child, fold, min_delta):
    """
    Simple significance flag for a branch copy-number change: the change must be at
    least ``min_delta`` copies AND at least a ``fold``-fold change vs the ancestral
    (parent) size. A birth from 0 (parent==0) is significant once it clears min_delta.
    """
    delta = abs(child - parent)
    if delta < min_delta:
        return False
    if parent == 0:
        return child >= min_delta
    if child == 0:
        return parent >= min_delta          # full loss of a previously multi-copy family
    ratio = max(child / parent, parent / child)
    return ratio >= fold


# ─── phylogeny / ancestral reconstruction ──────────────────────────────────────
def load_tree(tree_arg, hal, genomes_with_data):
    """
    Build an ete3 tree for the ancestral reconstruction, restricted to genomes we
    have data for. Source precedence: explicit --tree (newick file or literal) then
    the HAL (via halStats --tree). Returns (tree, None) or (None, reason).
    """
    try:
        import ete3
    except Exception:
        return None, "ete3 not installed"

    newick = None
    if tree_arg:
        if os.path.exists(tree_arg):
            with open(tree_arg) as fh:
                newick = fh.read().strip()
        else:
            newick = tree_arg.strip()
        tree = None
        for fmt in (1, 0, 3, 5):
            try:
                tree = ete3.Tree(newick, format=fmt)
                break
            except Exception:
                continue
        if tree is None:
            return None, "could not parse --tree newick"
    elif hal and os.path.exists(hal):
        try:
            sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            from tools.hal import get_tree
            tree = get_tree(hal)
        except Exception as e:
            return None, f"halStats/HAL tree unavailable ({e})"
    else:
        return None, "no --tree and no readable HAL"

    # Prune to the leaves we actually annotated (keep internal structure).
    keep = [l for l in tree.get_leaf_names() if l in genomes_with_data]
    if len(keep) < 2:
        return None, "fewer than 2 annotated genomes overlap the tree leaves"
    try:
        tree.prune(keep, preserve_branch_length=True)
    except Exception as e:
        return None, f"tree prune failed ({e})"

    # Ensure every node has a unique, stable name (HAL gives ancestors AncN; unnamed
    # internal nodes get nodeN in preorder so branch rows are addressable).
    counter = [0]
    for node in tree.traverse("preorder"):
        if not node.name:
            counter[0] += 1
            node.name = f"node{counter[0]}"
    return tree, None


def sankoff_reconstruct(tree, fixed_values, max_state):
    """
    Wagner/linear-cost (Sankoff) parsimony reconstruction of integer copy number.

    fixed_values : {node_name: int} observed counts. Always includes leaves; also
                   includes internal nodes for ancestral genomes that were themselves
                   annotated (their observed count is used instead of being inferred).
    Returns {node_name: int} for all nodes (observed nodes keep their value).
    """
    INF = float("inf")
    states = range(max_state + 1)
    # postorder: cost vectors
    for node in tree.traverse("postorder"):
        if node.name in fixed_values:
            v = fixed_values[node.name]
            node.add_feature("_cost", [0.0 if s == v else INF for s in states])
        elif node.is_leaf():
            # leaf without data (shouldn't happen after pruning): free
            node.add_feature("_cost", [0.0 for _ in states])
        else:
            cost = []
            for s in states:
                total = 0.0
                for ch in node.children:
                    total += min(ch._cost[t] + abs(s - t) for t in states)
                cost.append(total)
            node.add_feature("_cost", cost)
    # root assignment
    root = tree
    best_root = min(states, key=lambda s: root._cost[s])
    assign = {}

    def descend(node, parent_state):
        if node.name in fixed_values:
            state = fixed_values[node.name]
        elif node.is_root():
            state = best_root
        else:
            state = min(states, key=lambda t: node._cost[t] + abs(parent_state - t))
        assign[node.name] = state
        for ch in node.children:
            descend(ch, state)

    descend(root, best_root)
    return assign


def reconstruct_branch_changes(tree, family_leaf_counts, observed_internal, family_name):
    """
    Reconstruct ancestral family sizes and emit per-branch changes.
    Returns list of dicts (one per non-root branch).
    """
    fixed = dict(family_leaf_counts)
    fixed.update({k: v for k, v in observed_internal.items()})
    max_state = max(list(fixed.values()) + [0])
    if max_state == 0:
        return [], {}
    assign = sankoff_reconstruct(tree, fixed, max_state)
    rows = []
    for node in tree.traverse("preorder"):
        if node.is_root():
            continue
        parent = node.up
        pc = assign[parent.name]
        cc = assign[node.name]
        if pc == cc:
            direction = "stable"
        elif cc > pc:
            direction = "expansion"
        else:
            direction = "contraction"
        leaves = node.get_leaf_names()
        rows.append({
            "family": family_name,
            "branch": f"{parent.name}->{node.name}",
            "child_node": node.name,
            "child_is_leaf": node.is_leaf(),
            "parent_copies": pc,
            "child_copies": cc,
            "delta": cc - pc,
            "direction": direction,
            "observed": node.name in observed_internal or node.is_leaf(),
            "descendant_leaves": ",".join(sorted(leaves)[:8]),
        })
    return rows, assign


# ─── main ──────────────────────────────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--work-dir", required=True, help="CAT2 work_dir (has consensus_gene_set/ and databases/)")
    ap.add_argument("--config", help="CAT2 config yaml (to read ref_genome/genomes/hal). Optional.")
    ap.add_argument("--ref-genome", help="Reference genome name (overrides config).")
    ap.add_argument("--genomes", nargs="+", help="Target genomes (overrides config/auto-detect).")
    ap.add_argument("--hal", help="HAL file for phylogenetic ordering (overrides config).")
    ap.add_argument("--ref-db", help="Explicit reference sqlite db path (overrides work-dir/databases/{ref}.db).")
    ap.add_argument("--out-dir", help="Output dir (default {work_dir}/gene_family_analysis).")
    ap.add_argument("--biotype", default="protein_coding", help="Reference gene biotype to analyze (default protein_coding).")
    ap.add_argument("--family-regex", default=DEFAULT_FAMILY_REGEX,
                    help="Regex whose group(1) is the family root from a gene symbol.")
    ap.add_argument("--heatmap-top", type=int, default=40, help="How many most-variable families in the heatmap (0=off).")
    ap.add_argument("--tree", help="Newick file or literal string for ancestral reconstruction. "
                    "If omitted, the HAL tree (halStats --tree) is used.")
    ap.add_argument("--no-tree", action="store_true", help="Skip the phylogenetic ancestral reconstruction entirely.")
    ap.add_argument("--sig-fold", type=float, default=2.0,
                    help="Fold-change vs the reconstructed ancestral (parent) copy number for a "
                         "branch change to be flagged 'significant' (default 2.0 = doubling/halving).")
    ap.add_argument("--sig-min-delta", type=int, default=2,
                    help="Minimum absolute copy-number change for the significance flag (default 2). "
                         "Suppresses trivial 1->2 style calls.")
    args = ap.parse_args()

    cfg = load_yaml_config(args.config) if args.config else {}
    ref_genome = args.ref_genome or cfg.get("ref_genome")
    if not ref_genome:
        sys.exit("need --ref-genome (or --config with ref_genome).")
    hal = args.hal or cfg.get("hal")
    if hal and not os.path.isabs(hal):
        # config hal paths are relative to the repo/work_dir parent; try a couple of bases
        for base in (os.getcwd(), os.path.dirname(os.path.abspath(args.work_dir))):
            cand = os.path.join(base, hal)
            if os.path.exists(cand):
                hal = cand
                break

    consensus_dir = os.path.join(args.work_dir, "consensus_gene_set")
    if not os.path.isdir(consensus_dir):
        sys.exit(f"no consensus_gene_set/ under {args.work_dir}")

    genomes = args.genomes or cfg.get("genomes") or discover_genomes(consensus_dir)
    # never analyze the reference against itself; keep only genomes that have output
    have = set(discover_genomes(consensus_dir))
    genomes = [g for g in genomes if g != ref_genome and g in have]
    if not genomes:
        sys.exit("no target genomes with consensus output found.")
    genomes = order_genomes(genomes, hal, ref_genome)

    ref_db = args.ref_db or os.path.join(args.work_dir, "databases", f"{ref_genome}.db")
    ref_genes = load_reference_pc(ref_db, args.biotype)
    family_re = re.compile(args.family_regex)

    out_dir = args.out_dir or os.path.join(args.work_dir, "gene_family_analysis")
    os.makedirs(out_dir, exist_ok=True)

    # copy_number[source_gene][genome] = n_loci ; tx_count[source_gene][genome] = n_tx
    copy_number = defaultdict(lambda: defaultdict(int))
    tx_count = defaultdict(lambda: defaultdict(int))
    common_names = {}
    paralog_rows = []      # multi-copy loci
    novel_rows = []
    per_genome_loci = {}

    for genome in genomes:
        gp_info = os.path.join(consensus_dir, f"{genome}_consensus_novel_annotated.gp_info")
        gp = os.path.join(consensus_dir, f"{genome}_consensus_novel_annotated.gp")
        loci, novel, gene_meta = parse_consensus(gp_info, gp, args.biotype)
        per_genome_loci[genome] = loci
        for source_gene, gmap in loci.items():
            copy_number[source_gene][genome] = len(gmap)
            tx_count[source_gene][genome] = sum(len(r["tx"]) for r in gmap.values())
            if source_gene not in common_names:
                cn = gene_meta.get(source_gene, {}).get("common_name")
                common_names[source_gene] = cn or ref_genes.get(source_gene, {}).get("name", source_gene)
            if len(gmap) > 1:  # paralogs in this genome
                for gene_id, rec in sorted(gmap.items(),
                                           key=lambda kv: (kv[1]["chrom"], kv[1]["start"])):
                    paralog_rows.append({
                        "genome": genome, "source_gene": source_gene,
                        "common_name": common_names[source_gene],
                        "locus_gene_id": gene_id, "n_copies_in_genome": len(gmap),
                        "chrom": rec["chrom"], "start": rec["start"], "end": rec["end"],
                        "strand": rec["strand"], "n_transcripts": len(rec["tx"]),
                        "alignment_modes": ",".join(sorted(rec["modes"])),
                        "transcript_classes": ",".join(sorted(rec["classes"])),
                    })
        for n in novel:
            n2 = dict(n); n2["genome"] = genome
            novel_rows.append(n2)

    # universe of reference genes = ref PC genes (so losses are visible) unioned with
    # any source_gene actually seen (guards against ref-db/consensus mismatches).
    all_genes = set(ref_genes) | set(copy_number)
    for gid in all_genes:
        common_names.setdefault(gid, ref_genes.get(gid, {}).get("name", gid))

    # ─── ortholog copy-number matrix ───────────────────────────────────────────
    def write_matrix(path, valuemap):
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t")
            w.writerow(["source_gene", "common_name", "ref_biotype"] + list(genomes)
                       + ["n_genomes_present", "total_copies"])
            for gid in sorted(all_genes, key=lambda g: common_names[g].upper()):
                vals = [valuemap[gid].get(g, 0) for g in genomes]
                present = sum(1 for v in vals if v > 0)
                w.writerow([gid, common_names[gid],
                            ref_genes.get(gid, {}).get("biotype", "NA")]
                           + vals + [present, sum(vals)])

    write_matrix(os.path.join(out_dir, "ortholog_copy_matrix.tsv"), copy_number)
    write_matrix(os.path.join(out_dir, "ortholog_transcript_matrix.tsv"), tx_count)

    # ─── paralogs / novel ──────────────────────────────────────────────────────
    def dump(path, rows, fields):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow(r)

    dump(os.path.join(out_dir, "paralogs.tsv"),
         sorted(paralog_rows, key=lambda r: (-r["n_copies_in_genome"], r["common_name"], r["genome"])),
         ["genome", "source_gene", "common_name", "n_copies_in_genome", "locus_gene_id",
          "chrom", "start", "end", "strand", "n_transcripts", "alignment_modes",
          "transcript_classes"])
    dump(os.path.join(out_dir, "novel_genes.tsv"),
         sorted(novel_rows, key=lambda r: (r["genome"], r["chrom"], r["start"])),
         ["genome", "gene_id", "transcript_id", "chrom", "start", "end",
          "transcript_class", "alignment_mode"])

    # ─── gene families ─────────────────────────────────────────────────────────
    # family -> reference gene ids
    fam_ref = defaultdict(set)
    for gid in ref_genes:
        fam_ref[family_key(common_names[gid], gid, family_re)].add(gid)
    # include families that only appear via consensus (rare) so nothing is dropped
    for gid in copy_number:
        fam_ref[family_key(common_names[gid], gid, family_re)].add(gid)

    fam_rows = []
    ec_rows = []  # long-form expansion/contraction
    for fam, gids in fam_ref.items():
        ref_size = sum(1 for g in gids if g in ref_genes)
        genome_copies = {g: sum(copy_number[gid].get(g, 0) for gid in gids) for g in genomes}
        calls = {g: call_family(genome_copies[g], max(ref_size, 1)) for g in genomes}
        fam_rows.append({
            "family": fam if not fam.startswith("__singleton__") else common_names.get(next(iter(gids)), fam),
            "family_key": fam, "n_ref_genes": ref_size,
            "member_genes": ",".join(sorted(common_names[g] for g in gids)[:25]),
            **{g: genome_copies[g] for g in genomes},
            **{f"{g}__call": calls[g] for g in genomes},
        })
        for g in genomes:
            ec_rows.append({"level": "family", "id": fam,
                            "name": fam if not fam.startswith("__singleton__") else common_names.get(next(iter(gids)), fam),
                            "genome": g, "ref_copies": ref_size,
                            "obs_copies": genome_copies[g],
                            "delta": genome_copies[g] - ref_size, "call": calls[g]})

    # per-gene expansion/contraction long form (ref copy number is 1 per ref gene)
    for gid in sorted(all_genes, key=lambda g: common_names[g].upper()):
        ref_copies = 1 if gid in ref_genes else 0
        for g in genomes:
            obs = copy_number[gid].get(g, 0)
            ec_rows.append({"level": "gene", "id": gid, "name": common_names[gid],
                            "genome": g, "ref_copies": ref_copies, "obs_copies": obs,
                            "delta": obs - ref_copies, "call": call_gene(obs)})

    fam_cols = (["family", "family_key", "n_ref_genes", "member_genes"]
                + list(genomes) + [f"{g}__call" for g in genomes])
    with open(os.path.join(out_dir, "gene_family_matrix.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fam_cols, delimiter="\t")
        w.writeheader()
        for r in sorted(fam_rows, key=lambda r: (-r["n_ref_genes"], str(r["family"]))):
            w.writerow(r)
    dump(os.path.join(out_dir, "expansion_contraction.tsv"), ec_rows,
         ["level", "id", "name", "genome", "ref_copies", "obs_copies", "delta", "call"])

    # ─── presence / absence + lost / core genes ────────────────────────────────
    pa_rows = []
    lost_rows = []
    core_present = []      # gene present (>=1) in every genome
    core_single = []       # present as exactly 1 copy in every genome
    lost_in_all = []       # ref gene recovered in no genome
    ref_gene_ids = sorted(ref_genes, key=lambda g: common_names[g].upper())
    for gid in ref_gene_ids:
        pres = {g: (1 if copy_number[gid].get(g, 0) > 0 else 0) for g in genomes}
        n_present = sum(pres.values())
        pa_rows.append({"source_gene": gid, "common_name": common_names[gid],
                        **{g: pres[g] for g in genomes},
                        "n_present": n_present, "n_absent": len(genomes) - n_present})
        if n_present == 0:
            lost_in_all.append(gid)
        if n_present == len(genomes):
            core_present.append(gid)
            if all(copy_number[gid].get(g, 0) == 1 for g in genomes):
                core_single.append(gid)
        # lineage-specific / partial losses: absent in >=1 but present in >=1
        if 0 < n_present < len(genomes):
            absent = [g for g in genomes if pres[g] == 0]
            lost_rows.append({"source_gene": gid, "common_name": common_names[gid],
                              "n_absent": len(absent),
                              "absent_in": ",".join(absent),
                              "present_in": ",".join(g for g in genomes if pres[g] == 1)})

    with open(os.path.join(out_dir, "presence_absence_matrix.tsv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["source_gene", "common_name"] + list(genomes)
                           + ["n_present", "n_absent"], delimiter="\t")
        w.writeheader()
        for r in pa_rows:
            w.writerow(r)
    dump(os.path.join(out_dir, "lost_genes.tsv"),
         sorted(lost_rows, key=lambda r: -r["n_absent"]),
         ["source_gene", "common_name", "n_absent", "absent_in", "present_in"])
    with open(os.path.join(out_dir, "core_genes.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["source_gene", "common_name", "core_class"])
        for gid in core_present:
            w.writerow([gid, common_names[gid],
                        "single_copy_core" if gid in set(core_single) else "core"])
    pa_stats = {"core_present": len(core_present), "core_single": len(core_single),
                "lost_in_all": len(lost_in_all), "n_ref": len(ref_genes),
                "lineage_losses": len(lost_rows)}

    # ─── phylogenetic ancestral reconstruction (CAFE-style) ────────────────────
    tree_info = {"used": False, "reason": "disabled"}
    if not args.no_tree:
        # observed internal nodes = annotated ancestor genomes present in consensus dir
        all_annotated = set(discover_genomes(consensus_dir))
        tree, reason = load_tree(args.tree, hal, set(genomes) | all_annotated)
        if tree is None:
            tree_info = {"used": False, "reason": reason}
            sys.stderr.write(f"[warn] ancestral reconstruction skipped: {reason}\n")
        else:
            # For each family: leaf counts + observed internal (annotated ancestors).
            branch_rows = []
            branch_summary = defaultdict(lambda: {"expansions": 0, "contractions": 0,
                                                  "sig_expansions": 0, "sig_contractions": 0,
                                                  "genes_gained": 0, "genes_lost": 0})
            tree_leaves = set(tree.get_leaf_names())
            internal_names = {n.name for n in tree.traverse() if not n.is_leaf()}
            observed_internal_genomes = internal_names & all_annotated
            # parse observed internal ancestor copy numbers (family sums) lazily
            anc_copy = defaultdict(lambda: defaultdict(int))  # anc_genome -> source_gene -> loci
            for anc in sorted(observed_internal_genomes):
                gp_info = os.path.join(consensus_dir, f"{anc}_consensus_novel_annotated.gp_info")
                gp = os.path.join(consensus_dir, f"{anc}_consensus_novel_annotated.gp")
                loci, _, _ = parse_consensus(gp_info, gp, args.biotype)
                for sg, gmap in loci.items():
                    anc_copy[anc][sg] = len(gmap)

            real_families = [f for f in fam_rows
                             if not f["family_key"].startswith("__singleton__")]
            for f in real_families:
                gids = [g for g in ref_genes if family_key(common_names[g], g, family_re) == f["family_key"]]
                gids += [g for g in copy_number if family_key(common_names[g], g, family_re) == f["family_key"]]
                gids = set(gids)
                leaf_counts = {g: sum(copy_number[gid].get(g, 0) for gid in gids)
                               for g in tree_leaves}
                observed_internal = {anc: sum(anc_copy[anc].get(gid, 0) for gid in gids)
                                     for anc in observed_internal_genomes}
                rows, _ = reconstruct_branch_changes(tree, leaf_counts, observed_internal,
                                                     f["family"])
                for r in rows:
                    r["significant"] = is_significant(r["parent_copies"], r["child_copies"],
                                                      args.sig_fold, args.sig_min_delta)
                    if r["direction"] == "expansion":
                        branch_summary[r["branch"]]["expansions"] += 1
                        branch_summary[r["branch"]]["genes_gained"] += r["delta"]
                        if r["significant"]:
                            branch_summary[r["branch"]]["sig_expansions"] += 1
                    elif r["direction"] == "contraction":
                        branch_summary[r["branch"]]["contractions"] += 1
                        branch_summary[r["branch"]]["genes_lost"] += -r["delta"]
                        if r["significant"]:
                            branch_summary[r["branch"]]["sig_contractions"] += 1
                    if r["direction"] != "stable":
                        branch_rows.append(r)

            dump(os.path.join(out_dir, "branch_changes.tsv"),
                 sorted(branch_rows, key=lambda r: (-int(r["significant"]), -abs(r["delta"]))),
                 ["family", "branch", "child_node", "child_is_leaf", "parent_copies",
                  "child_copies", "delta", "direction", "significant", "observed",
                  "descendant_leaves"])
            dump(os.path.join(out_dir, "significant_changes.tsv"),
                 sorted((r for r in branch_rows if r["significant"]),
                        key=lambda r: -abs(r["delta"])),
                 ["family", "branch", "child_node", "child_is_leaf", "parent_copies",
                  "child_copies", "delta", "direction", "significant", "observed",
                  "descendant_leaves"])
            with open(os.path.join(out_dir, "branch_summary.tsv"), "w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["branch", "expansions", "contractions", "sig_expansions",
                            "sig_contractions", "families_gained_copies", "families_lost_copies"])
                for br, s in sorted(branch_summary.items(),
                                    key=lambda kv: -(kv[1]["sig_expansions"] + kv[1]["sig_contractions"])):
                    w.writerow([br, s["expansions"], s["contractions"], s["sig_expansions"],
                                s["sig_contractions"], s["genes_gained"], s["genes_lost"]])
            try:
                newick_out = tree.write(format=1)
                with open(os.path.join(out_dir, "tree.nwk"), "w") as fh:
                    fh.write(newick_out + "\n")
            except Exception:
                pass
            n_sig = sum(1 for r in branch_rows if r["significant"])
            tree_info = {"used": True, "n_leaves": len(tree_leaves),
                         "n_observed_ancestors": len(observed_internal_genomes),
                         "n_significant": n_sig, "sig_fold": args.sig_fold,
                         "sig_min_delta": args.sig_min_delta,
                         "sig_changes": sorted((r for r in branch_rows if r["significant"]),
                                               key=lambda r: -abs(r["delta"]))[:15],
                         "branch_summary": dict(branch_summary)}

    # ─── summary.md ────────────────────────────────────────────────────────────
    write_summary(out_dir, genomes, ref_genome, ref_genes, copy_number, common_names,
                  paralog_rows, novel_rows, fam_rows, hal, pa_stats, tree_info)

    # ─── heatmap ───────────────────────────────────────────────────────────────
    if args.heatmap_top and pd is not None:
        try:
            make_heatmap(out_dir, genomes, fam_rows, args.heatmap_top)
        except Exception as e:
            sys.stderr.write(f"[warn] heatmap skipped: {e}\n")

    print(f"Wrote gene-family analysis to {out_dir}")
    print(f"  genomes analyzed  : {len(genomes)}")
    print(f"  reference PC genes: {len(ref_genes)}")
    print(f"  gene families     : {sum(1 for f in fam_rows if not f['family_key'].startswith('__singleton__'))} "
          f"(+{sum(1 for f in fam_rows if f['family_key'].startswith('__singleton__'))} singletons)")
    print(f"  core (all genomes): {pa_stats['core_present']}  |  lineage losses: {pa_stats['lineage_losses']}"
          f"  |  lost in all: {pa_stats['lost_in_all']}")
    print(f"  ancestral recon   : {'on' if tree_info.get('used') else 'skipped (' + str(tree_info.get('reason')) + ')'}"
          + (f"  |  significant branch changes: {tree_info.get('n_significant')}" if tree_info.get('used') else ""))


def write_summary(out_dir, genomes, ref_genome, ref_genes, copy_number, common_names,
                  paralog_rows, novel_rows, fam_rows, hal, pa_stats=None, tree_info=None):
    n_ref = len(ref_genes)
    lines = []
    lines.append(f"# Protein-coding gene family / ortholog analysis")
    lines.append("")
    lines.append(f"- Reference: **{ref_genome}** ({n_ref} protein-coding genes)")
    lines.append(f"- Target genomes: {len(genomes)}"
                 + (" (ordered by phylogenetic distance from reference)" if hal else ""))
    lines.append(f"- Orthology: reference `source_gene` shared across genomes. "
                 f"Copy number = distinct consensus loci per gene per genome.")
    lines.append("")

    # per-genome table
    lines.append("## Per-genome summary")
    lines.append("")
    lines.append("| genome | genes found | single-copy | expanded (>1) | lost | extra copies | novel PC loci |")
    lines.append("|---|---|---|---|---|---|---|")
    novel_by_genome = defaultdict(int)
    for r in novel_rows:
        novel_by_genome[r["genome"]] += 1
    novel_loci = defaultdict(set)
    for r in novel_rows:
        novel_loci[r["genome"]].add(r["gene_id"])
    for g in genomes:
        found = single = expanded = extra = 0
        for gid in ref_genes:
            c = copy_number[gid].get(g, 0)
            if c > 0:
                found += 1
            if c == 1:
                single += 1
            elif c > 1:
                expanded += 1
                extra += c - 1
        lost = n_ref - found
        lines.append(f"| {g} | {found} ({100*found/max(n_ref,1):.1f}%) | {single} | "
                     f"{expanded} | {lost} | {extra} | {len(novel_loci[g])} |")
    lines.append("")

    # biggest expansions (per gene, max copies across genomes)
    lines.append("## Largest single-gene expansions")
    lines.append("")
    lines.append("| gene | genome | copies | (chrom loci) |")
    lines.append("|---|---|---|---|")
    top_par = sorted(paralog_rows, key=lambda r: -r["n_copies_in_genome"])
    seen = set()
    shown = 0
    # collapse to one row per (gene,genome)
    perkey = {}
    for r in paralog_rows:
        k = (r["common_name"], r["genome"])
        perkey.setdefault(k, {"copies": r["n_copies_in_genome"], "chroms": []})
        perkey[k]["chroms"].append(r["chrom"])
    for (name, genome), v in sorted(perkey.items(), key=lambda kv: -kv[1]["copies"]):
        if shown >= 25:
            break
        chroms = ",".join(sorted(set(v["chroms"]))[:5])
        lines.append(f"| {name} | {genome} | {v['copies']} | {chroms} |")
        shown += 1
    lines.append("")

    # family-level expansion / contraction highlights (exclude singletons)
    real_fams = [f for f in fam_rows if not f["family_key"].startswith("__singleton__")
                 and f["n_ref_genes"] >= 1]
    lines.append("## Gene families with the most variable size")
    lines.append("")
    lines.append("Copy number per genome vs the reference family size "
                 "(number of protein-coding genes in that family in the reference).")
    lines.append("")
    def variability(f):
        vals = [f[g] for g in genomes]
        return max(vals) - min(vals) if vals else 0
    top_fams = sorted(real_fams, key=lambda f: (-(variability(f)), -f["n_ref_genes"]))[:20]
    header = "| family | ref genes | " + " | ".join(genomes) + " |"
    lines.append(header)
    lines.append("|---|---|" + "|".join(["---"] * len(genomes)) + "|")
    for f in top_fams:
        row = f"| {f['family']} | {f['n_ref_genes']} | " + " | ".join(str(f[g]) for g in genomes) + " |"
        lines.append(row)
    lines.append("")

    # presence / absence + core / lost
    if pa_stats:
        lines.append("## Presence / absence (gene conservation)")
        lines.append("")
        lines.append(f"- **Core** (present in all {len(genomes)} genomes): "
                     f"{pa_stats['core_present']} / {pa_stats['n_ref']} "
                     f"({100*pa_stats['core_present']/max(pa_stats['n_ref'],1):.1f}%)")
        lines.append(f"- **Single-copy core** (exactly 1 copy in every genome): "
                     f"{pa_stats['core_single']}")
        lines.append(f"- **Lineage / partial losses** (absent in ≥1 but not all): "
                     f"{pa_stats['lineage_losses']}")
        lines.append(f"- **Lost in all** (recovered in no genome): {pa_stats['lost_in_all']}")
        lines.append("")

    # phylogenetic branch analysis
    if tree_info and tree_info.get("used"):
        lines.append("## Phylogenetic expansion / contraction (ancestral reconstruction)")
        lines.append("")
        lines.append(f"Wagner (linear-cost) parsimony over {tree_info['n_leaves']} genomes"
                     + (f" plus {tree_info['n_observed_ancestors']} annotated ancestral nodes "
                        "(their observed copy numbers are used directly instead of inferred)"
                        if tree_info.get("n_observed_ancestors") else "")
                     + ". Per-branch counts of families gaining/losing copies:")
        lines.append("")
        lines.append("| branch (parent→child) | families expanded | families contracted | significant | net copies gained | net copies lost |")
        lines.append("|---|---|---|---|---|---|")
        bs = tree_info["branch_summary"]
        for br, s in sorted(bs.items(),
                            key=lambda kv: -(kv[1]["expansions"] + kv[1]["contractions"]))[:15]:
            sig = s.get("sig_expansions", 0) + s.get("sig_contractions", 0)
            lines.append(f"| {br} | {s['expansions']} | {s['contractions']} | {sig} | "
                         f"{s['genes_gained']} | {s['genes_lost']} |")
        lines.append("")
        lines.append(f"### Significant changes "
                     f"(≥{tree_info['sig_min_delta']} copies and ≥{tree_info['sig_fold']}×"
                     f" vs ancestor): {tree_info['n_significant']}")
        lines.append("")
        if tree_info.get("sig_changes"):
            lines.append("| family | branch | direction | parent→child | Δ | observed |")
            lines.append("|---|---|---|---|---|---|")
            for r in tree_info["sig_changes"]:
                lines.append(f"| {r['family']} | {r['branch']} | {r['direction']} | "
                             f"{r['parent_copies']}→{r['child_copies']} | {r['delta']:+d} | "
                             f"{'yes' if r['observed'] else 'inferred'} |")
            lines.append("")
    elif tree_info:
        lines.append("## Phylogenetic expansion / contraction")
        lines.append("")
        lines.append(f"_Skipped: {tree_info.get('reason')}._ Provide `--tree <newick>` or run "
                     "where `halStats` is available to enable ancestral reconstruction.")
        lines.append("")

    lines.append("## Files")
    lines.append("")
    for fn, desc in [
        ("ortholog_copy_matrix.tsv", "gene x genome copy-number (loci) matrix"),
        ("ortholog_transcript_matrix.tsv", "gene x genome transcript-count matrix"),
        ("presence_absence_matrix.tsv", "gene x genome 1/0 presence matrix"),
        ("lost_genes.tsv", "genes absent in some (but not all) genomes -> lineage-specific losses"),
        ("core_genes.tsv", "genes present in every genome (core / single-copy core)"),
        ("paralogs.tsv", "one row per locus for every multi-copy gene (the paralogs)"),
        ("gene_family_matrix.tsv", "family x genome total copies + per-genome call"),
        ("expansion_contraction.tsv", "long-form calls (lost/single/expanded/contracted) per gene and family"),
        ("branch_changes.tsv", "per-family, per-branch ancestral copy-number changes (CAFE-style), with significance flag"),
        ("significant_changes.tsv", "only the branch changes flagged significant (fold + min-delta thresholds)"),
        ("branch_summary.tsv", "per-branch totals of family expansions / contractions (incl. significant)"),
        ("tree.nwk", "pruned tree used for reconstruction"),
        ("novel_genes.tsv", "genome-specific / de-novo protein-coding loci with no reference ortholog"),
        ("copy_number_heatmap.png", "heatmap of the most copy-number-variable families"),
    ]:
        lines.append(f"- `{fn}` — {desc}")
    lines.append("")
    with open(os.path.join(out_dir, "summary.md"), "w") as fh:
        fh.write("\n".join(lines))


def make_heatmap(out_dir, genomes, fam_rows, top):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    real = [f for f in fam_rows if not f["family_key"].startswith("__singleton__")]
    def variability(f):
        vals = [f[g] for g in genomes]
        return (max(vals) - min(vals)) if vals else 0
    real = sorted(real, key=lambda f: -variability(f))[:top]
    if not real:
        return
    mat = [[f[g] for g in genomes] for f in real]
    labels = [f["family"] for f in real]
    df = pd.DataFrame(mat, index=labels, columns=genomes)
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * len(genomes) + 3), max(4, 0.3 * len(real))))
    im = ax.imshow(df.values, aspect="auto", cmap="viridis")
    ax.set_xticks(range(len(genomes)))
    ax.set_xticklabels(genomes, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    for i in range(len(labels)):
        for j in range(len(genomes)):
            ax.text(j, i, int(df.values[i, j]), ha="center", va="center",
                    color="white", fontsize=6)
    fig.colorbar(im, ax=ax, label="copy number (loci)")
    ax.set_title("Most copy-number-variable protein-coding gene families")
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "copy_number_heatmap.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
