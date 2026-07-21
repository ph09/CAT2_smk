import argparse
import collections
import itertools
import json
import logging
import os
import sys
import warnings
import tempfile
import shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from collections import OrderedDict, defaultdict, Counter
from matplotlib.backends.backend_pdf import PdfPages
import tools.nameConversions
import tools.sqlInterface

# set matplotlib backend and style
matplotlib.rcParams['pdf.fonttype'] = 42

# set backend to avoid display errors
sns.set_style('ticks')
warnings.filterwarnings('ignore')
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO)
BAR_WIDTH = 0.45
BOXPLOT_SATURATION = 0.7

class atomic_file:
    def __init__(self, path):
        self.path = path
        self.temp_path = None
        self.file = None
        
        dirname = os.path.dirname(self.path)
        if dirname:
            os.makedirs(dirname, exist_ok=True)
        else:
            dirname = '.'
        
        fd, self.temp_path = tempfile.mkstemp(dir=dirname)
        os.close(fd)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            shutil.move(self.temp_path, self.path)
        else:
            if self.temp_path and os.path.exists(self.temp_path):
                os.unlink(self.temp_path)

    def move_to_final_destination(self):
        """Move temp file to final path (for use outside context manager)."""
        if self.temp_path and os.path.exists(self.temp_path):
            shutil.move(self.temp_path, self.path)

def main():
    """
    Main entry point for the plotting script. Parses arguments and generates all plots.
    """
    parser = argparse.ArgumentParser(description="Generates a series of PDF plots visualizing metrics from the CAT pipeline.")
    # Input Data
    parser.add_argument("--tm-jsons", nargs=2, action='append', required=True, metavar=("GENOME", "PATH"),
                        help="transMap metrics JSON file. Specify once for each genome.")
    parser.add_argument("--metrics-jsons", nargs=2, action='append', required=True, metavar=("GENOME", "PATH"),
                        help="Consensus metrics JSON file. Specify once for each genome.")
    parser.add_argument("--dbs", nargs=2, action='append', required=True, metavar=("GENOME", "PATH"),
                        help="Genome database file. Specify once for each genome.")
    parser.add_argument("--annotation-db", required=True, help="Path to the reference annotation database.")
    # Plotting Parameters
    parser.add_argument("--ordered-genomes", nargs='+', required=True, help="Ordered list of genomes for plots.")
    parser.add_argument("--pb-genomes", nargs='*', default=[], help="List of genomes with PacBio IsoSeq data.")
    # Output Directory
    parser.add_argument("--out-dir", required=True, help="Directory to write all plot PDFs to.")

    args = parser.parse_args()
    
    os.makedirs(args.out_dir, exist_ok=True)
    generate_plots(args)

def generate_plots(args):
    """
    Generates all plots based on the provided arguments.
    
    :param args: An argparse.Namespace object containing all paths and parameters.
    """
    # Load and structure the input data from the provided paths
    tm_data = collections.OrderedDict(args.tm_jsons)
    consensus_data = collections.OrderedDict(args.metrics_jsons)
    dbs = collections.OrderedDict(args.dbs)
    
    # Load transMap data, but skip if file doesn't exist (e.g., for reference genome)
    for genome, path in list(tm_data.items()):
        if not os.path.exists(path):
            logger.warning(f"TransMap metrics JSON file not found: {path}. Skipping genome {genome} for transMap metrics.")
            # Remove this genome from tm_data since we can't load its metrics
            del tm_data[genome]
            continue
        with open(path) as f:
            tm_data[genome] = json.load(f)
    
    for genome, path in list(consensus_data.items()):
        if not os.path.exists(path):
            logger.warning(f"Consensus metrics JSON file not found: {path}. Skipping genome {genome}.")
            del consensus_data[genome]
            continue
        with open(path) as f:
            consensus_data[genome] = json.load(f)

    # Load data from databases
    tm_metrics = load_tm_metrics(dbs)
    transcript_biotype_map = tools.sqlInterface.get_transcript_biotype_map(args.annotation_db)
    gene_biotype_map = tools.sqlInterface.get_gene_biotype_map(args.annotation_db)
    biotypes = sorted(tools.sqlInterface.get_transcript_biotypes(args.annotation_db))

    if 'protein_coding' in biotypes:
        biotypes.insert(0, biotypes.pop(biotypes.index('protein_coding')))

    # Define output paths based on the output directory
    plot_paths = {
        name: os.path.join(args.out_dir, f"{name}.pdf")
        for name in [
            "transcript_modes", "transmap_coverage", "transmap_identity", "paralogy", "unfiltered_paralogy",
            "gene_family_collapse", "coverage", "identity", "missing_genes_transcripts",
            "consensus_annotation_support", "consensus_extrinsic_support", "completeness", "coding_indels",
            "split_genes", "denovo", "IsoSeq_isoform_validation", "augustus_improvement"
        ]
    }

    def _has_key(data, key):
        """Check if any genome's data has a non-empty value for the given key."""
        for d in data.values():
            if not isinstance(d, dict):
                continue
            val = d.get(key)
            if val is None:
                continue
            if isinstance(val, dict) and not val:
                continue
            if isinstance(val, (list, str)) and len(val) == 0:
                continue
            return True
        return False

    # Generate each plot (all conditional on data availability)
    tx_modes_plot(consensus_data, args.ordered_genomes, plot_paths["transcript_modes"])

    if tm_metrics.get('transMap Coverage') and tm_metrics.get('transMap Identity'):
        tm_metrics_plot(tm_metrics, args.ordered_genomes, biotypes, transcript_biotype_map,
                        plot_paths["transmap_coverage"], plot_paths["transmap_identity"])

    if _has_key(tm_data, 'Paralogy'):
        tm_para_plot(tm_data, args.ordered_genomes, biotypes, plot_paths["paralogy"], plot_paths["unfiltered_paralogy"])

    if _has_key(tm_data, 'Gene Family Collapse'):
        tm_gene_family_plot(tm_data, args.ordered_genomes, biotypes, plot_paths["gene_family_collapse"])

    if _has_key(consensus_data, 'Coverage') and _has_key(consensus_data, 'Identity'):
        consensus_metrics_plot(consensus_data, args.ordered_genomes, biotypes, plot_paths["coverage"], plot_paths["identity"])

    if _has_key(consensus_data, 'Gene Missing') and _has_key(consensus_data, 'Transcript Missing'):
        missing_rate_plot(consensus_data, args.ordered_genomes, biotypes, plot_paths["missing_genes_transcripts"])

    consensus_support_plot(consensus_data, args.ordered_genomes, biotypes,
                           modes=['Splice Annotation Support', 'Exon Annotation Support', 'Original Introns'],
                           title='Reference annotation support', tgt=plot_paths["consensus_annotation_support"])
    consensus_support_plot(consensus_data, args.ordered_genomes, biotypes,
                           modes=['Splice Support', 'Exon Support'],
                           title='Extrinsic support', tgt=plot_paths["consensus_extrinsic_support"])

    if _has_key(consensus_data, 'Completeness'):
        completeness_plot(consensus_data, args.ordered_genomes, biotypes, plot_paths["completeness"],
                          gene_biotype_map, transcript_biotype_map)

    if _has_key(consensus_data, 'transMap Indels') and _has_key(consensus_data, 'Consensus Indels'):
        indel_plot(consensus_data, args.ordered_genomes, plot_paths["coding_indels"])

    if _has_key(tm_data, 'Split Genes'):
        split_genes_plot(tm_data, args.ordered_genomes, plot_paths["split_genes"])
    if _has_key(consensus_data, 'denovo'):
        denovo_plot(consensus_data, args.ordered_genomes, plot_paths["denovo"])
    if _has_key(consensus_data, 'IsoSeq Transcript Validation') and args.pb_genomes:
        pb_support_plot(consensus_data, args.ordered_genomes, args.pb_genomes, plot_paths["IsoSeq_isoform_validation"])
    if _has_key(consensus_data, 'Evaluation Improvement'):
        improvement_plot(consensus_data, args.ordered_genomes, plot_paths["augustus_improvement"])

def load_tm_metrics(dbs):
    """Loads transMap data from PSLs"""
    tm_metrics = {'transMap Coverage': OrderedDict(), 'transMap Identity': OrderedDict()}
    tm_name_map = {'TransMapCoverage': 'transMap Coverage', 'TransMapIdentity': 'transMap Identity'}
    for genome, db_path in dbs.items():
        session = tools.sqlInterface.start_session(db_path)
        table = tools.sqlInterface.TmEval
        for classifier in ['TransMapCoverage', 'TransMapIdentity']:
            query = session.query(table.AlignmentId, table.value).filter(table.classifier == classifier)
            tm_metrics[tm_name_map[classifier]][genome] = dict(query.all())
    return tm_metrics


###
# Plots
###


def tm_metrics_plot(tm_metrics, ordered_genomes, biotypes, transcript_biotype_map, tm_coverage_tgt, tm_identity_tgt):
    """plots for transMap coverage, identity"""
    tm_iter = list(zip(*[['transMap Coverage', 'transMap Identity'],
                    [tm_coverage_tgt, tm_identity_tgt]]))
    for mode, tgt in tm_iter:
        df = dict_to_df_with_biotype(tm_metrics[mode], transcript_biotype_map)
        df = pd.melt(df, id_vars='biotype', value_vars=ordered_genomes).dropna()
        df.columns = ['biotype', 'genome', mode]
        cov_ident_plot(biotypes, ordered_genomes, mode, tgt, df, x=mode, y='genome')


def consensus_metrics_plot(consensus_data, ordered_genomes, biotypes, coverage_tgt, identity_tgt):
    """plots for consensus coverage, identity, score"""
    cons_iter = list(zip(*[['Coverage', 'Identity'],
                      [coverage_tgt, identity_tgt]]))
    for mode, tgt in cons_iter:
        df = json_to_df_with_biotype(consensus_data, mode)
        cov_ident_plot(biotypes, ordered_genomes, mode, tgt, df, x=mode, y='genome')


def consensus_support_plot(consensus_data, ordered_genomes, biotypes, modes, title, tgt):
    """grouped violin plots of original intron / intron annotation / exon annotation support"""
    available_modes = [m for m in modes if any(m in d for d in consensus_data.values())]
    if not available_modes:
        return

    # Skip if all data values are zero (would produce empty-looking violins)
    has_nonzero = False
    for mode in available_modes:
        for d in consensus_data.values():
            if mode in d and isinstance(d[mode], dict):
                for vals in d[mode].values():
                    if isinstance(vals, list) and any(v != 0 for v in vals):
                        has_nonzero = True
                        break
            if has_nonzero:
                break
        if has_nonzero:
            break
    if not has_nonzero:
        return

    def adjust_plot(g, this_title):
        g.set_xticklabels(rotation=90)
        g.fig.suptitle(this_title)
        g.fig.subplots_adjust(top=0.9)
        for ax in g.axes.flat:
            ax.set_ylabel('Percent supported')
            ax.set_ylim(-1, 101)

    dfs = []
    for i, mode in enumerate(available_modes):
        df = json_to_df_with_biotype(consensus_data, mode)
        if i > 0:
            df = df[mode]
        dfs.append(df)
    df = pd.concat(dfs, axis=1)
    df = pd.melt(df, value_vars=available_modes, id_vars=['genome', 'biotype'])
    af = atomic_file(tgt)
    with PdfPages(af.temp_path) as pdf:
        if len(ordered_genomes) > 1:
            g = sns.catplot(data=df, y='value', x='genome', col='variable', col_wrap=2, kind='violin', sharex=True,
                               sharey=True, row_order=ordered_genomes, cut=0)
        else:
            g = sns.catplot(data=df, y='value', x='variable', kind='violin', sharex=True,
                               sharey=True, row_order=ordered_genomes, cut=0)
        adjust_plot(g, title)
        multipage_close(pdf, tight_layout=False)
        title += ' for {}'
        for biotype in biotypes:
            this_title = title.format(biotype)
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                if len(ordered_genomes) > 1:
                    g = sns.catplot(data=biotype_df, y='value', x='genome', col='variable', col_wrap=2,
                                       kind='violin', sharex=True, sharey=True, row_order=ordered_genomes, cut=0)
                else:
                    g = sns.catplot(data=df, y='value', x='variable', kind='violin', sharex=True,
                                       sharey=True, row_order=ordered_genomes, cut=0)
                adjust_plot(g, this_title)
                multipage_close(pdf, tight_layout=False)
    af.move_to_final_destination()


def tm_para_plot(tm_data, ordered_genomes, biotypes, para_tgt, unfiltered_para_tgt):
    """transMap paralogy plots"""
    for key, tgt in [['Paralogy', para_tgt], ['UnfilteredParalogy', unfiltered_para_tgt]]:
        legend_labels = ['= 1', '= 2', '= 3', '\u2265 4']
        title_string = 'Proportion of transcripts that have multiple alignments'
        biotype_title_string = 'Proportion of {} transcripts that have multiple alignments'
        df = json_biotype_nested_counter_to_df(tm_data, key)
        # we want a dataframe where each row is the counts, in genome order
        # we construct the transpose first
        r = []
        df[key] = pd.to_numeric(df[key])
        # make sure genomes are in order
        df['genome'] = pd.Categorical(df['genome'], ordered_genomes, ordered=True)
        df = df.sort_values('genome')
        for biotype, biotype_df in df.groupby('biotype'):
            for genome, genome_df in biotype_df.groupby('genome'):
                high_para = genome_df[genome_df[key] >= 4]['count'].sum()
                counts = dict(list(zip(genome_df[key], genome_df['count'])))
                r.append([biotype, genome, counts.get(1, 0), counts.get(2, 0), counts.get(3, 0), high_para])
        df = pd.DataFrame(r, columns=['biotype', 'genome', '1', '2', '3', '\u2265 4'])
        numeric_cols = ['1', '2', '3', '\u2265 4']
        sum_df = df.groupby('genome', sort=False)[numeric_cols].sum().T

        plot_fn = generic_unstacked_barplot if len(df.columns) <= 5 else generic_stacked_barplot
        box_label = 'Number of\nalignments'
        af = atomic_file(tgt)
        with PdfPages(af.temp_path) as pdf:
            plot_fn(sum_df, pdf, title_string, legend_labels, 'Number of transcripts', ordered_genomes, box_label)
            for biotype in biotypes:
                biotype_df = biotype_filter(df, biotype)
                if biotype_df is not None:
                    biotype_df = biotype_df.drop(['genome', 'biotype'], axis=1).T
                    title_string = biotype_title_string.format(biotype)
                    plot_fn(biotype_df, pdf, title_string, legend_labels, 'Number of transcripts', ordered_genomes,
                            box_label)
        af.move_to_final_destination()


def tm_gene_family_plot(tm_data, ordered_genomes, biotypes, gene_family_tgt):
    """transMap gene family collapse plots."""
    try:
        df = json_biotype_nested_counter_to_df(tm_data, 'Gene Family Collapse')
        if df.empty:
            raise ValueError("empty")
    except (ValueError, KeyError):
        return
    df['Gene Family Collapse'] = pd.to_numeric(df['Gene Family Collapse'])
    tot_df = df[['Gene Family Collapse', 'genome', 'count']].\
        groupby(['genome', 'Gene Family Collapse']).aggregate(sum).reset_index()
    tot_df = tot_df.sort_values('Gene Family Collapse')
    af = atomic_file(gene_family_tgt)
    with PdfPages(af.temp_path) as pdf:
        g = sns.catplot(y='count', col='genome', x='Gene Family Collapse', data=tot_df, kind='bar',
                           col_order=ordered_genomes, col_wrap=4)
        g.fig.suptitle('Number of genes collapsed during gene family collapse')
        g.set_xlabels('Number of genes collapsed to one locus')
        g.set_ylabels('Number of genes')
        g.fig.subplots_adjust(top=0.9)
        multipage_close(pdf, tight_layout=False)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is None:
                continue
            biotype_df = biotype_df.sort_values('Gene Family Collapse')
            g = sns.catplot(y='count', col='genome', x='Gene Family Collapse', data=biotype_df, kind='bar',
                               col_order=[x for x in ordered_genomes if x in set(biotype_df.genome)], col_wrap=4)
            g.fig.suptitle('Number of genes collapsed during gene family collapse for {}'.format(biotype))
            g.set_xlabels('Number of genes collapsed to one locus')
            g.set_ylabels('Number of genes')
            g.fig.subplots_adjust(top=0.9)
            multipage_close(pdf, tight_layout=False)
    af.move_to_final_destination()


def missing_rate_plot(consensus_data, ordered_genomes, biotypes, missing_plot_tgt):
    """Missing genes/transcripts"""
    base_title = 'Number of missing orthologs in consensus set'
    gene_missing_df = json_biotype_counter_to_df(consensus_data, 'Gene Missing')
    transcript_missing_df = json_biotype_counter_to_df(consensus_data, 'Transcript Missing')
    if gene_missing_df.empty or transcript_missing_df.empty:
        return
    gene_missing_df.columns = ['biotype', 'Genes', 'genome']
    transcript_missing_df.columns = ['biotype', 'Transcripts', 'genome']
    df = transcript_missing_df.merge(gene_missing_df, on=['genome', 'biotype'])
    if df.empty:
        return
    df = pd.melt(df, id_vars=['biotype', 'genome'])
    ylabel = 'Number of genes or transcripts'
    af = atomic_file(missing_plot_tgt)
    with PdfPages(af.temp_path) as pdf:
        tot_df = df.groupby(['genome', 'biotype', 'variable'])['value'].sum().reset_index()
        generic_barplot(tot_df, pdf, '', ylabel, base_title, x='genome', y='value',
                        col='variable', row_order=ordered_genomes)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is None:
                continue
            biotype_df = biotype_df.groupby(['genome', 'variable'])['value'].sum().reset_index()
            title = base_title + ' for biotype {}'.format(biotype)
            generic_barplot(biotype_df, pdf, '', ylabel, title, x='genome', y='value',
                            col='variable', row_order=ordered_genomes)
    af.move_to_final_destination()


def tx_modes_plot(consensus_data, ordered_genomes, tx_mode_plot_tgt):
    if not any('Transcript Modes' in d for d in consensus_data.values()):
        return
    ordered_groups = ['transMap', 'transMap+TM', 'transMap+TMR', 'transMap+TM+TMR',
                      'TM', 'TMR', 'TM+TMR',
                      'augTM', 'augTMR', 'augPB',
                      'PB', 'strg', 'txTM', 'augMP', 'exRef', 'other']
    ordered_groups = OrderedDict([[frozenset(x.split('+')), x] for x in ordered_groups])

    def split_fn(s):
        mode = str(s['Transcript Modes'])
        mode = mode.replace('_pairwise', '')
        key = frozenset(mode.split(','))
        if key in ordered_groups:
            return ordered_groups[key]
        stripped = mode.replace('aug', '')
        return ordered_groups.get(frozenset(stripped.split(',')), 'other')

    modes_df = json_biotype_counter_to_df(consensus_data, 'Transcript Modes')
    if len(modes_df) == 0:
        return
    pivoted = modes_df.pivot(index='genome', columns='Transcript Modes')
    # Genomes whose consensus produced no transcripts (e.g. tiny test data)
    # are absent from the pivot; reindex so every requested column exists.
    pivoted = pivoted.reindex(ordered_genomes, fill_value=0)
    df = pivoted.transpose().reset_index()
    df['Modes'] = df.apply(split_fn, axis=1)
    df = df[['Modes'] + ordered_genomes]
    ordered_values = [x for x in ordered_groups.values() if x in set(df['Modes'])]
    af = atomic_file(tx_mode_plot_tgt)
    with PdfPages(af.temp_path) as pdf:
        title_string = 'Transcript modes in protein coding consensus gene set'
        ylabel = 'Number of transcripts'
        if len(ordered_genomes) > 1:
            df['Ordered Modes'] = pd.Categorical(df['Modes'], ordered_values, ordered=True)
            df = df.sort_values('Ordered Modes')
            df = df[['Ordered Modes'] + ordered_genomes].set_index('Ordered Modes')
            df = df.fillna(0)
            generic_stacked_barplot(df, pdf, title_string, df.index, ylabel, ordered_genomes, 'Transcript mode(s)',
                                    bbox_to_anchor=(1.25, 0.7))

        else:
            generic_barplot(pd.melt(df, id_vars='Modes'), pdf, 'Transcript mode(s)', ylabel, title_string, x='Modes',
                            y='value', order=ordered_values)
    af.move_to_final_destination()


def denovo_plot(consensus_data, ordered_genomes, denovo_tgt):
    try:
        df = json_biotype_nested_counter_to_df(consensus_data, 'denovo')
    except (ValueError, KeyError):
        return
    if df.empty:
        return
    df.columns = ['Result', 'Number of transcripts', 'De novo mode', 'genome']
    non_zero = df[df['Number of transcripts'] > 0]
    if non_zero.empty:
        return
    af = atomic_file(denovo_tgt)
    with PdfPages(af.temp_path) as pdf:
        has_multiple_modes = len(set(df['De novo mode'])) > 1
        if len(set(df.genome)) > 1:
            if has_multiple_modes:
                ax = sns.catplot(data=df, x='genome', y='Number of transcripts', kind='bar', col='Result',
                                    hue='De novo mode', col_wrap=2, row_order=ordered_genomes, sharex=True,
                                    sharey=False)
            else:
                ax = sns.catplot(data=df, x='genome', y='Number of transcripts', kind='bar', col='Result',
                                    col_wrap=2, row_order=ordered_genomes, sharex=True, sharey=False)
        else:
            if has_multiple_modes:
                ax = sns.catplot(data=df, x='Result', y='Number of transcripts', kind='bar', hue='De novo mode')
            else:
                ax = sns.catplot(data=df, x='Result', y='Number of transcripts', kind='bar')
        ax.set_xticklabels(rotation=90)
        ax.fig.suptitle('Incorporation of de-novo predictions')
        ax.fig.subplots_adjust(top=0.9)
        multipage_close(pdf, tight_layout=False)
    af.move_to_final_destination()


def split_genes_plot(tm_data, ordered_genomes, split_plot_tgt):
    df = json_biotype_counter_to_df(tm_data, 'Split Genes')
    if df.empty:
        return
    df.columns = ['category', 'count', 'genome']
    if df['count'].sum() == 0:
        return
    af = atomic_file(split_plot_tgt)
    with PdfPages(af.temp_path) as pdf:
        title = 'Split genes'
        if len(ordered_genomes) > 1:
            g = generic_barplot(pdf=pdf, data=df, x='genome', y='count', col='category', xlabel='', col_wrap=2,
                                sharey=False, ylabel='Number of transcripts or genes', row_order=ordered_genomes,
                                title=title)
        else:
            g = generic_barplot(pdf=pdf, data=df, x='category', y='count', ylabel='Number of transcripts or genes',
                                title=title, xlabel='Category')
    af.move_to_final_destination()


def pb_support_plot(consensus_data, ordered_genomes, pb_genomes, pb_support_tgt):
    pb_genomes = [x for x in ordered_genomes if x in pb_genomes]
    df = json_biotype_counter_to_df(consensus_data, 'IsoSeq Transcript Validation')
    if len(df) == 0:
        return
    df.columns = ['IsoSeq Transcript Validation', 'Number of transcripts', 'genome']
    af = atomic_file(pb_support_tgt)
    with PdfPages(af.temp_path) as pdf:
        ax = sns.catplot(data=df, x='genome', y='Number of transcripts', hue='IsoSeq Transcript Validation',
                            kind='bar', row_order=pb_genomes)
        ax.set_xticklabels(rotation=90)
        ax.fig.suptitle('Isoforms validated by at least one IsoSeq read')
        multipage_close(pdf, tight_layout=False)
    af.move_to_final_destination()


def completeness_plot(consensus_data, ordered_genomes, biotypes, completeness_plot_tgt, gene_biotype_map,
                      transcript_biotype_map):
    def adjust_plot(g, gene_count, tx_count):
        for ax, c in zip(*[g.axes[0], [gene_count, tx_count]]):
            _ = ax.set_ylim(0, c)
            ax.spines['top'].set_edgecolor('#e74c3c')
            ax.spines['top'].set_linewidth(2)
            ax.spines['top'].set_visible(True)
            ax.spines['top'].set_linestyle('dashed')

    df = json_grouped_biotype_nested_counter_to_df(consensus_data, 'Completeness')
    if df.empty:
        return
    af = atomic_file(completeness_plot_tgt)
    with PdfPages(af.temp_path) as pdf:
        tot_df = df.groupby(by=['genome', 'category'])['count'].sum().reset_index()
        tot_df = sort_long_df(tot_df, ordered_genomes)
        title = 'Number of comparative genes/transcripts present'
        g = generic_barplot(pdf=pdf, data=tot_df, x='genome', y='count', col='category', xlabel='',
                            sharey=False, ylabel='Number of genes/transcripts', title=title,
                            col_order=['Gene', 'Transcript'], close=False, palette=choose_palette(ordered_genomes))
        adjust_plot(g, len(gene_biotype_map), len(transcript_biotype_map))
        multipage_close(pdf, tight_layout=False)
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                biotype_df = sort_long_df(biotype_df, ordered_genomes)
                gene_biotype_count = len({i for i, b in gene_biotype_map.items() if b == biotype})
                tx_biotype_count = len({i for i, b in transcript_biotype_map.items() if b == biotype})
                title = 'Number of comparative genes/transcripts present for biotype {}'.format(biotype)
                g = generic_barplot(pdf=pdf, data=biotype_df, x='genome', y='count', col='category', xlabel='',
                                    sharey=False, ylabel='Number of genes/transcripts',
                                    title=title, col_order=['Gene', 'Transcript'], close=False,
                                    palette=choose_palette(ordered_genomes))
                adjust_plot(g, gene_biotype_count, tx_biotype_count)
                multipage_close(pdf, tight_layout=False)
    af.move_to_final_destination()


def improvement_plot(consensus_data, ordered_genomes, improvement_tgt):
    def do_kdeplot(x, y, ax, n_levels=None, bw='scott'):
        try:
            sns.kdeplot(x, y, ax=ax, cut=0, cmap='Purples_d', shade=True, shade_lowest=False, n_levels=n_levels, bw=bw,
                        rasterized=True)
        except:
            logger.warning('Unable to do a KDE fit to AUGUSTUS improvement.')
            pass

    af = atomic_file(improvement_tgt)
    with PdfPages(af.temp_path) as pdf, sns.axes_style("whitegrid"):
        for genome in ordered_genomes:
            data = pd.DataFrame(consensus_data[genome]['Evaluation Improvement']['changes'])
            unchanged = consensus_data[genome]['Evaluation Improvement']['unchanged']
            if len(data) == 0:
                continue
            data.columns = ['transMap original introns',
                            'transMap intron annotation support',
                            'transMap intron RNA support',
                            'Original introns',
                            'Intron annotation support',
                            'Intron RNA support',
                            'transMap alignment goodness',
                            'Alignment goodness']
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2)
            for ax in [ax1, ax2, ax3, ax4]: 
                ax.set_xlim(0, 100)
                ax.set_ylim(0, 100)
            
            do_kdeplot(data['transMap original introns'], data['Original introns'], ax1, n_levels=25, bw=2)
            sns.regplot(x=data['transMap original introns'], y=data['Original introns'], ax=ax1,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7, 'rasterized': True}, fit_reg=False)
            do_kdeplot(data['transMap intron annotation support'], data['Intron annotation support'], ax2,
                       n_levels=25, bw=2)
            sns.regplot(x=data['transMap intron annotation support'], y=data['Intron annotation support'], ax=ax2,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7, 'rasterized': True}, fit_reg=False)          
            do_kdeplot(data['transMap intron RNA support'], data['Intron RNA support'], ax3, n_levels=25, bw=2)
            sns.regplot(x=data['transMap intron RNA support'], y=data['Intron RNA support'], ax=ax3,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7, 'rasterized': True}, fit_reg=False)
            
            do_kdeplot(data['transMap alignment goodness'], data['Alignment goodness'], ax4, n_levels=20, bw=1)
            sns.regplot(x=data['transMap alignment goodness'], y=data['Alignment goodness'], ax=ax4,
                        color='#A9B36F', scatter_kws={"s": 3, 'alpha': 0.7, 'rasterized': True}, fit_reg=False)

            fig.suptitle('AUGUSTUS metric improvements for {:,} transcripts in {}.\n'
                         '{:,} transMap transcripts were chosen.'.format(len(data), genome, unchanged))
            
            for ax in [ax1, ax2, ax3, ax4]:
                ax.set(adjustable='box', aspect='equal')
            fig.subplots_adjust(hspace=0.3)
            multipage_close(pdf, tight_layout=False)
    af.move_to_final_destination()


def indel_plot(consensus_data, ordered_genomes, indel_plot_tgt):
    af = atomic_file(indel_plot_tgt)
    with PdfPages(af.temp_path) as pdf:
        tm_df = pd.concat([pd.DataFrame.from_dict(consensus_data[genome]['transMap Indels'], orient='index').T
                           for genome in ordered_genomes])
        try:  # this is a hack to deal with weird small input datasets
            tm_df['genome'] = ordered_genomes
        except:
            return
        tm_df['transcript set'] = ['transMap'] * len(tm_df)
        consensus_df = pd.concat([pd.DataFrame.from_dict(consensus_data[genome]['Consensus Indels'], orient='index').T
                                  for genome in ordered_genomes])
        consensus_df['genome'] = ordered_genomes
        consensus_df['transcript set'] = ['Consensus'] * len(consensus_df)
        df = pd.concat([consensus_df, tm_df], ignore_index=True)
        df = pd.melt(df, id_vars=['genome', 'transcript set'],
                     value_vars=['CodingDeletion', 'CodingInsertion', 'CodingMult3Indel'])
        df.columns = ['Genome', 'Transcript set', 'Type', 'Percent of transcripts']
        g = sns.catplot(data=df, x='Genome', y='Percent of transcripts', col='Transcript set',
                           hue='Type', kind='bar', row_order=ordered_genomes,
                           col_order=['transMap', 'Consensus'])
        g.set_xticklabels(rotation=90)
        g.fig.subplots_adjust(top=.8)
        g.fig.suptitle('Coding indels')
        multipage_close(pdf, tight_layout=False)
    af.move_to_final_destination()


###
# shared plotting functions
###


def cov_ident_plot(biotypes, ordered_genomes, mode, tgt, df, x=None, y=None, xlabel=None):
    """violin plots for coverage and identity."""
    if xlabel is None:
        xlabel = 'Percent {}'.format(mode)
    af = atomic_file(tgt)
    with PdfPages(af.temp_path) as pdf:
        title = 'Overall {}'.format(mode)
        xmin = int(min(df[mode]))
        horizontal_violin_plot(df, ordered_genomes, title, xlabel, pdf, x=x, y=y, xlim=(xmin, 100))
        for biotype in biotypes:
            biotype_df = biotype_filter(df, biotype)
            if biotype_df is not None:
                title = '{} for biotype {}'.format(mode, biotype)
                xmin = int(min(df[mode]))
                horizontal_violin_plot(biotype_df, ordered_genomes, title, xlabel, pdf, x=x, y=y, xlim=(xmin, 100))
    af.move_to_final_destination()

###
# generic plotting functions
###


def generic_barplot(data, pdf, xlabel, ylabel, title, row_order=None, x=None, y=None, hue=None, hue_order=None,
                    order=None, col=None, col_wrap=None, sharex=True, sharey=True, col_order=None, palette=None,
                    close=True):
    g = sns.catplot(data=data, x=x, y=y, hue=hue, ci=None, kind='bar', hue_order=hue_order, row_order=row_order,
                       col=col, col_wrap=col_wrap, sharex=sharex, sharey=sharey, col_order=col_order, palette=palette,
                       order=order)
    g.set_xticklabels(rotation=90)
    g.fig.suptitle(title)
    g.fig.subplots_adjust(top=.8)
    g.set_axis_labels(xlabel, ylabel)
    try:  # depending on columns, axes could be flat or not
        axes = list(itertools.chain.from_iterable(g.axes))
    except TypeError:
        axes = g.axes
    for ax in axes:
        ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=10, steps=[1, 2, 5, 10], integer=True))
        ax.margins(y=0.15)
        ax.autoscale(enable=True, axis='y', tight=False)
        ax.set_ylim(0, ax.get_ylim()[1])
    if close == True:
        multipage_close(pdf, tight_layout=False)
    return g


def horizontal_violin_plot(data, ordered_genomes, title, xlabel, pdf, hue=None, x=None, y=None, xlim=None):
    """not so generic function that specifically produces a paired boxplot/violinplot"""
    fig, ax = plt.subplots()
    sns.violinplot(data=data, x=x, y=y, hue=hue, order=ordered_genomes, palette=choose_palette(ordered_genomes),
                   saturation=BOXPLOT_SATURATION, orient='h', cut=0, scale='count', ax=ax)
    fig.suptitle(title)
    ax.set_xlabel(xlabel)
    if xlim is not None:
        ax.set_xlim(xlim)
    multipage_close(pdf, tight_layout=False)


def _generic_histogram(bars, legend_labels, title_string, pdf, ax, fig, ylabel, names, box_label, bbox_to_anchor):
    fig.legend([x[0] for x in bars[::-1]], legend_labels[::-1], bbox_to_anchor=bbox_to_anchor, frameon=True,
               title=box_label)
    ax.set_title(title_string)
    ax.set_ylabel(ylabel)
    set_ticks(names, ax)
    ax.xaxis.set_ticks(np.arange(0, len(names)) + BAR_WIDTH / 2.0)
    sns.despine(top=True, right=True)
    multipage_close(pdf)


def generic_unstacked_barplot(df, pdf, title_string, legend_labels, ylabel, names, box_label,
                              bbox_to_anchor=(1.12, 0.7)):
    fig, ax = plt.subplots()
    bars = []
    shorter_bar_width = BAR_WIDTH / len(df)
    for i, (_, d) in enumerate(df.iterrows()):
        bars.append(ax.bar(np.arange(len(df.columns)) + shorter_bar_width * i, d, shorter_bar_width,
                           color=sns.color_palette()[i], linewidth=0.0))
    _generic_histogram(bars, legend_labels, title_string, pdf, ax, fig, ylabel, names, box_label, bbox_to_anchor)


def generic_stacked_barplot(df, pdf, title_string, legend_labels, ylabel, names, box_label, bbox_to_anchor=(1.12, 0.7)):
    fig, ax = plt.subplots()
    bars = []
    cumulative = np.zeros(len(df.columns))
    color_palette = choose_palette(legend_labels)
    for i, (_, d) in enumerate(df.iterrows()):
        bars.append(ax.bar(np.arange(len(df.columns)), d, BAR_WIDTH, bottom=cumulative,
                           color=color_palette[i], linewidth=0.0))
        cumulative += d
    _generic_histogram(bars, legend_labels, title_string, pdf, ax, fig, ylabel, names, box_label, bbox_to_anchor)


###
# Shared functions
###


def json_flat_to_df(consensus_data, key):
    """converts cases where we have exactly genome:value pairs"""
    r = []
    for genome, d in consensus_data.items():
        r.append([genome, d[key]])
    return pd.DataFrame(r)


def json_to_df_with_biotype(consensus_data, key):
    """converts JSON entries with many transcripts, such as those for coverage/identity"""
    dfs = []
    for genome, d in consensus_data.items():
        if key not in d:
            continue
        for biotype, vals in d[key].items():
            df = pd.DataFrame(vals)
            if len(df) > 0:
                df.columns = [key]
                df = df.assign(biotype=[biotype] * len(df), genome=[genome] * len(df))
            dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def json_biotype_nested_counter_to_df(consensus_data, key):
    """converts the JSON entries with nested counts. Expects the first level keys to be biotypes"""
    dfs = []
    for genome, d in consensus_data.items():
        if key in d:
            for biotype, vals in d[key].items():
                df = pd.DataFrame(list(vals.items()))
                if len(df) > 0:
                    df.columns = [key, 'count']
                    df = df.assign(biotype=[biotype] * len(df), genome=[genome] * len(df))
                    dfs.append(df)
    if not dfs:
        raise ValueError(f"No data found for key '{key}'")
    return pd.concat(dfs, ignore_index=True)


def json_grouped_biotype_nested_counter_to_df(consensus_data, key):
    """converts the JSON entries with nested counts. Expects the second level keys to be biotypes"""
    dfs = []
    for genome, d in consensus_data.items():
        if key not in d:
            continue
        for group, vals in d[key].items():
            df = pd.DataFrame(list(vals.items()))
            if len(df) > 0:
                df.columns = ['biotype', 'count']
                df = df.assign(category=[group] * len(df), genome=[genome] * len(df))
            dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def json_biotype_counter_to_df(consensus_data, key):
    """converts the JSON entries with nested counts. Expects the first level keys to be biotypes"""
    dfs = []
    for genome, d in consensus_data.items():
        if key not in d:
            continue
        vals = d[key]
        df = pd.DataFrame(list(vals.items()))
        if len(df) > 0:
            df.columns = [key, 'count']
            df = df.assign(genome=[genome] * len(df))
        dfs.append(df)
    if not dfs:
        return pd.DataFrame()
    return pd.concat(dfs, ignore_index=True)


def dict_to_df_with_biotype(data, transcript_biotype_map):
    df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in data.items()]))
    try:
        df['biotype'] = [transcript_biotype_map[tx] for tx in df.index]
    except KeyError:
        # try removing names
        df['biotype'] = [transcript_biotype_map[tools.nameConversions.strip_alignment_numbers(tx)] for tx in df.index]
    return df


def biotype_filter(df, biotype):
    df = df[df.biotype == biotype]
    return df if len(df) > 0 else None


def multipage_close(pdf, tight_layout=True):
    """convenience function for closing up a pdf page"""
    if tight_layout:
        plt.tight_layout()
    pdf.savefig(bbox_inches='tight')
    plt.close('all')


def choose_palette(ordered_genomes):
    """choose palette in cases where genomes get different colors"""
    if len(ordered_genomes) <= 6:
        return sns.color_palette()
    else:
        return sns.color_palette("Set2", len(ordered_genomes))


def set_ticks(names, ax, nbins=10.0):
    ax.margins(y=0.15)
    ax.autoscale(enable=True, axis='y', tight=False)
    ax.set_ylim(0, plt.ylim()[1])
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(nbins=nbins, steps=[1, 2, 5, 10], integer=True))
    ax.xaxis.set_major_locator(matplotlib.ticker.LinearLocator(len(names)))
    ax.yaxis.set_minor_locator(matplotlib.ticker.AutoMinorLocator())
    ax.xaxis.set_ticklabels(names, rotation=90)


def sort_long_df(df, ordered_genomes):
    """sorts a long form dataframe by ordered genomes"""
    ordered_index = dict(list(zip(ordered_genomes, list(range(len(ordered_genomes))))))
    df['order'] = df['genome'].map(ordered_index)
    df = df.sort_values('order')
    return df.drop('order', axis=1)

if __name__ == "__main__":
    main()
