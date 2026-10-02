"""
Functions to interface with the sqlite databases produced by various steps of the annotation pipeline
"""
from . import transcripts

import pandas as pd
from sqlalchemy import Column, Integer, Text, Float, Boolean, func, create_engine, inspect
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker

###
# Data model
###


Base = declarative_base()


class Annotation(Base):
    """Table for the annotation table. Only exists in ref_genome"""
    __tablename__ = 'annotation'
    GeneId = Column(Text, primary_key=True)
    TranscriptId = Column(Text, primary_key=True)
    TranscriptName = Column(Text)
    GeneName = Column(Text)
    GeneBiotype = Column(Text)
    TranscriptBiotype = Column(Text)
    ExtraTags = Column(Text)


class Bed12(object):
    """General table description for storing BED12 features"""
    chromosome = Column(Text)
    start = Column(Integer)
    stop = Column(Integer)
    name = Column(Text)
    score = Column(Integer)
    strand = Column(Text)
    thickStart = Column(Integer)
    thickStop = Column(Integer)
    rgb = Column(Text)
    blockCount = Column(Integer)
    blockSizes = Column(Text)
    blockStarts = Column(Text)


class EvaluationColumns(Bed12):
    """Mixin class for all TranscriptEvaluation module tables. Represents a bed12 with a leading ID column"""
    AlignmentId = Column(Text, primary_key=True)


class MrnaTmEval(EvaluationColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from transMap"""
    __tablename__ = 'mRNA_transMap_Evaluation'


class MrnaAugTmEval(EvaluationColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusTM"""
    __tablename__ = 'mRNA_augTM_Evaluation'


class MrnaAugTmrEval(EvaluationColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusTMR"""
    __tablename__ = 'mRNA_augTMR_Evaluation'

class MrnaAugMpEval(EvaluationColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusMP"""
    __tablename__ = 'mRNA_augMP_Evaluation'

class MrnaTxTmEval(EvaluationColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from txTM (in-house transcript-level transMap; formerly the external Liftoff tool)"""
    __tablename__ = 'mRNA_txTM_Evaluation'


class CdsTmEval(EvaluationColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from transMap"""
    __tablename__ = 'CDS_transMap_Evaluation'


class CdsAugTmEval(EvaluationColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusTM"""
    __tablename__ = 'CDS_augTM_Evaluation'


class CdsAugTmrEval(EvaluationColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusTMR"""
    __tablename__ = 'CDS_augTMR_Evaluation'

class CdsAugMpEval(EvaluationColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusMP"""
    __tablename__ = 'CDS_augMP_Evaluation'

class CdsTxTmEval(EvaluationColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from txTM (in-house transcript-level transMap; formerly the external Liftoff tool)"""
    __tablename__ = 'CDS_txTM_Evaluation'


class MetricsColumns(object):
    """Mixin class for all TranscriptMetrics module tables"""
    AlignmentId = Column(Text, primary_key=True)
    classifier = Column(Text)
    value = Column(Float)


class TmEval(MetricsColumns, Base):
    """Table for evaluations from TransMapEvaluation module"""
    __tablename__ = 'TransMapEvaluation'
    TranscriptId = Column(Text, primary_key=True)
    GeneId = Column(Text, primary_key=True)


class TmFilterEval(MetricsColumns, Base):
    """Table for evaluations from FilterTransMap module. This table is stored in a stacked format for simplicity."""
    __tablename__ = 'TransMapFilterEvaluation'
    GeneId = Column(Text, primary_key=True)
    TranscriptId = Column(Text, primary_key=True)
    AlignmentId = Column(Text, primary_key=True)
    GeneAlternateContigs = Column(Text)
    GeneAlternateLoci = Column(Text)
    CollapsedGeneNames = Column(Text)
    CollapsedGeneIds = Column(Text)
    Paralogy = Column(Text)


class TmPwFilterEval(MetricsColumns, Base):
    """Table for evaluations from FilterTransMap module for pairwise (BAM-based) chains. This table is stored in a stacked format for simplicity."""
    __tablename__ = 'TransMapPairwiseFilterEvaluation'
    GeneId = Column(Text, primary_key=True)
    TranscriptId = Column(Text, primary_key=True)
    AlignmentId = Column(Text, primary_key=True)
    GeneAlternateContigs = Column(Text)
    GeneAlternateLoci = Column(Text)
    CollapsedGeneNames = Column(Text)
    CollapsedGeneIds = Column(Text)
    Paralogy = Column(Text)
    UnfilteredParalogy = Column(Text)
    PossibleSplitGeneLocations = Column(Text)


class TmMetrics(MetricsColumns, Base):
    """Table for evaluations from TransMapMetrics module"""
    __tablename__ = 'TransMapMetrics'


class MrnaTmMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from transMap"""
    __tablename__ = 'mRNA_transMap_Metrics'


class MrnaTmPairwiseMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from transMap_pairwise"""
    __tablename__ = 'mRNA_transMap_pairwise_Metrics'


class MrnaAugTmMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusTM"""
    __tablename__ = 'mRNA_augTM_Metrics'


class MrnaAugTmPairwiseMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusTM_pairwise"""
    __tablename__ = 'mRNA_augTM_pairwise_Metrics'


class MrnaAugTmrMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusTMR"""
    __tablename__ = 'mRNA_augTMR_Metrics'


class MrnaAugTmrPairwiseMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusTMR_pairwise"""
    __tablename__ = 'mRNA_augTMR_pairwise_Metrics'

class MrnaAugMpMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from AugustusMP"""
    __tablename__ = 'mRNA_augMP_Metrics'

class MrnaTxTmMetrics(MetricsColumns, Base):
    """Table for evaluations of mRNA alignments of transcripts derived from txTM (in-house transcript-level transMap; formerly the external Liftoff tool)"""
    __tablename__ = 'mRNA_txTM_Metrics'


class CdsTmMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from transMap"""
    __tablename__ = 'CDS_transMap_Metrics'


class CdsTmPairwiseMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from transMap_pairwise"""
    __tablename__ = 'CDS_transMap_pairwise_Metrics'


class CdsAugTmMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusTM"""
    __tablename__ = 'CDS_augTM_Metrics'


class CdsAugTmPairwiseMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusTM_pairwise"""
    __tablename__ = 'CDS_augTM_pairwise_Metrics'


class CdsAugTmrMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusTMR"""
    __tablename__ = 'CDS_augTMR_Metrics'


class CdsAugTmrPairwiseMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusTMR_pairwise"""
    __tablename__ = 'CDS_augTMR_pairwise_Metrics'

class CdsAugMpMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from AugustusMP"""
    __tablename__ = 'CDS_augMP_Metrics'

class CdsTxTmMetrics(MetricsColumns, Base):
    """Table for evaluations of CDS alignments of transcripts derived from txTM (in-house transcript-level transMap; formerly the external Liftoff tool)"""
    __tablename__ = 'CDS_txTM_Metrics'


class AlternativeGeneIdColumns(object):
    """mixin class for AlternativeGenes"""
    TranscriptId = Column(Text, primary_key=True)
    AssignedGeneId = Column(Text)
    AlternativeGeneIds = Column(Text)
    ResolutionMethod = Column(Text)


class AugPbAlternativeGenes(AlternativeGeneIdColumns, Base):
    """Table for recording a list of alternative parental genes for IsoSeq"""
    __tablename__ = 'AugPbAlternativeGenes'


class StrgAlternativeGenes(AlternativeGeneIdColumns, Base):
    """Table for recording a list of alternative parental genes for StringTie (strg mode)"""
    __tablename__ = 'StrgAlternativeGenes'


class ExRefAlternativeGenes(AlternativeGeneIdColumns, Base):
    """Table for recording a list of alternative parental genes for external references"""
    __tablename__ = 'ExRef_AlternativeGenes'


class IsoSeqExonStructures(Bed12, Base):
    """Table for recording all distinct exon structures present in a IsoSeq hints file"""
    __tablename__ = 'IsoSeqExonStructures'
    index = Column(Integer, primary_key=True)


###
# Wrapper functions for setting up sessions
###


# SQLite connections are thread-affined; when SQLAlchemy's engine pool is
# garbage-collected from a finalizer/atexit thread, sqlite3 raises
# "SQLite objects created in a thread can only be used in that same thread".
# We only use one session per thread at a time, so disabling the check is safe
# and silences harmless teardown noise across the pipeline.
_SQLITE_CONNECT_ARGS = {"check_same_thread": False}


def _make_engine(db_path):
    return create_engine('sqlite:///' + db_path, connect_args=_SQLITE_CONNECT_ARGS)


def start_session(db_path):
    """basic script for starting a session"""
    engine = _make_engine(db_path)
    Session = sessionmaker(bind=engine)
    return Session()


###
# Dictionary mapping tables to their respective transcript/alignment modes
###


tables = {'CDS': {'augTM': {'metrics': CdsAugTmMetrics, 'evaluation': CdsAugTmEval},
                  'augTMR': {'metrics': CdsAugTmrMetrics, 'evaluation': CdsAugTmrEval},
                  'augMP': {'metrics': CdsAugMpMetrics, 'evaluation': CdsAugMpEval},
                  'transMap': {'metrics': CdsTmMetrics, 'evaluation': CdsTmEval},
                  'transMap_pairwise': {'metrics': CdsTmPairwiseMetrics, 'evaluation': CdsTmEval},
                  'augTM_pairwise': {'metrics': CdsAugTmPairwiseMetrics, 'evaluation': CdsAugTmEval},
                  'augTMR_pairwise': {'metrics': CdsAugTmrPairwiseMetrics, 'evaluation': CdsAugTmrEval},
                  'txTM': {'metrics': CdsTxTmMetrics, 'evaluation': CdsTxTmEval}},
          'mRNA': {'augTM': {'metrics': MrnaAugTmMetrics, 'evaluation': MrnaAugTmEval},
                   'augTMR': {'metrics': MrnaAugTmrMetrics, 'evaluation': MrnaAugTmrEval},
                   'augMP': {'metrics': MrnaAugMpMetrics, 'evaluation': MrnaAugMpEval},
                   'transMap': {'metrics': MrnaTmMetrics, 'evaluation': MrnaTmEval},
                   'transMap_pairwise': {'metrics': MrnaTmPairwiseMetrics, 'evaluation': MrnaTmEval},
                   'augTM_pairwise': {'metrics': MrnaAugTmPairwiseMetrics, 'evaluation': MrnaAugTmEval},
                   'augTMR_pairwise': {'metrics': MrnaAugTmrPairwiseMetrics, 'evaluation': MrnaAugTmrEval},
                   'txTM': {'metrics': MrnaTxTmMetrics, 'evaluation': MrnaTxTmEval}},
          'alt_names': {'exRef': ExRefAlternativeGenes,
                        'augPB': AugPbAlternativeGenes,
                        'strg': StrgAlternativeGenes}}


###
# Attributes functions -- read data from the annotation table
###


def read_attrs(db_path, table=Annotation.__tablename__, index_col='TranscriptId'):
    """
    Read the attributes database file into a pandas DataFrame
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: pandas DataFrame
    """
    engine = _make_engine(db_path)
    return pd.read_sql_table(table, engine, index_col=index_col)


def get_transcript_gene_map(db_path, table=Annotation.__tablename__, index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to gene IDs.
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: GeneId}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(list(zip(df.index, df.GeneId)))


def get_gene_transcript_map(db_path, table=Annotation.__tablename__, index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to gene IDs.
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {GeneId: [tx_id1, tx_id2, etc]}
    """
    df = read_attrs(db_path, table, index_col).reset_index()
    r = {}
    for gene_id, s in df.groupby('GeneId'):
        r[gene_id] = s.TranscriptId.tolist()
    return r


def get_transcript_biotype_map(db_path, table=Annotation.__tablename__, index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping transcript IDs to their biotype
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: tx_biotype}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(list(zip(df.index, df.TranscriptBiotype)))


def get_gene_biotype_map(db_path, table=Annotation.__tablename__, index_col='TranscriptId'):
    """
    Convenience wrapper for read_attrs that returns a dictionary mapping gene IDs to their biotype
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :param index_col: column to index on. should generally be tx_id.
    :return: dictionary {tx_id: tx_biotype}
    """
    df = read_attrs(db_path, table, index_col)
    return dict(list(zip(df.GeneId, df.GeneBiotype)))


def get_transcript_biotypes(db_path, table=Annotation):
    """
    Returns a set of transcript biotypes seen in this annotation set
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :return: dictionary {tx_id: tx_biotype}
    """
    session = start_session(db_path)
    query = session.query(table.TranscriptBiotype).distinct()
    return {x[0] for x in query.all()}


def get_gene_biotypes(db_path, table=Annotation):
    """
    Returns a set of transcript biotypes seen in this annotation set
    :param db_path: path to the attributes database
    :param table: table name. should generally be annotation
    :return: dictionary {tx_id: tx_biotype}
    """
    session = start_session(db_path)
    query = session.query(table.GeneBiotype).distinct()
    return {x[0] for x in query.all()}


###
# Loading entire tables
###


def load_annotation(ref_db_path):
    """
    Load the reference annotation table
    :param ref_db_path: path to reference genome database. Must have table Annotation.__tablename__
    :return: DataFrame
    """
    engine = _make_engine(ref_db_path)
    df = pd.read_sql_table(Annotation.__tablename__, engine)
    return df


def load_alignment_evaluation(db_path):
    """
    Loads the transMap alignment evaluation table
    :param db_path: path to genome database
    :return: DataFrame
    """
    engine = _make_engine(db_path)
    df = pd.read_sql_table(TmEval.__tablename__, engine)
    df = pd.pivot_table(df, index=['TranscriptId', 'AlignmentId'], columns='classifier', values='value')
    return df.reset_index()


def load_filter_evaluation(db_path):
    """
    Loads the transMap alignment filtering evaluation table
    :param db_path: path to genome database
    :return: DataFrame
    """
    engine = _make_engine(db_path)
    return pd.read_sql_table(TmFilterEval.__tablename__, engine)


def load_pairwise_filter_evaluation(db_path):
    """
    Loads the transMap pairwise (BAM-based) alignment filtering evaluation table
    :param db_path: path to genome database
    :return: DataFrame
    """
    engine = _make_engine(db_path)
    try:
        return pd.read_sql_table(TmPwFilterEval.__tablename__, engine)
    except ValueError:
        # Table doesn't exist yet (for older databases)
        return pd.DataFrame()


def load_isoseq_txs(db_path):
    """
    Loads the table IsoSeqExonStructures, constructing actual ChromosomeInterval objects.
    :param db_path: path to genome db
    :return: list of Transcript objects
    """
    engine = _make_engine(db_path)
    
    # Check if the table exists
    inspector = inspect(engine)
    if not inspector.has_table(IsoSeqExonStructures.__tablename__):
        return []  # Return empty list if table doesn't exist
    
    df = pd.read_sql_table(IsoSeqExonStructures.__tablename__, engine, index_col='index')
    if df.empty:
        return []  # Return empty list if table is empty
    
    txs = [transcripts.Transcript(list(s)) for _, s in df.iterrows()]
    return txs


def load_evaluation(table, session):
    """
    load evaluation entries for this gene. Makes use of count() and group by to get the # of times the classifier failed
    :param table: One of the evaluation tables
    :param session: Active sqlalchemy session.
    :return: DataFrame
    """
    assert any(table == cls for cls in (MrnaAugTmrEval, MrnaAugTmEval, MrnaAugMpEval, MrnaTmEval, MrnaTxTmEval,
                                        CdsAugTmrEval, CdsAugTmEval, CdsAugMpEval, CdsTmEval, CdsTxTmEval))
    query = session.query(table.AlignmentId, table.name, func.count(table.name).label('value')). \
        group_by(table.AlignmentId, table.name)
    return pd.read_sql(query.statement, session.bind)


def load_metrics(table, session):
    """
    load metrics entries for this gene. Wrapper for generic_gene_query.
    :param table: One of the metrics tables
    :param session: Active sqlalchemy session.
    :return: DataFrame
    """
    assert any(table == cls for cls in (MrnaAugTmrMetrics, MrnaAugTmMetrics, MrnaAugMpMetrics, MrnaTmMetrics, MrnaTxTmMetrics,
                                        MrnaTmPairwiseMetrics, MrnaAugTmPairwiseMetrics, MrnaAugTmrPairwiseMetrics,
                                        CdsAugTmrMetrics, CdsAugTmMetrics, CdsAugMpMetrics, CdsTmMetrics, CdsTxTmMetrics,
                                        CdsTmPairwiseMetrics, CdsAugTmPairwiseMetrics, CdsAugTmrPairwiseMetrics))
    query = session.query(table)
    return pd.read_sql(query.statement, session.bind)


def load_alternatives(table, session):
    """
    load de novo parental assignment + alternative parents
    :param table: One of AugPbAlternativeGenes, StrgAlternativeGenes, or ExRefAlternativeGenes
    :param session: Active sqlalchemy session.
    :return: DataFrame
    """
    assert table in (AugPbAlternativeGenes, StrgAlternativeGenes, ExRefAlternativeGenes)
    query = session.query(table)
    return pd.read_sql(query.statement, session.bind)


###
# Stats functions
###

def load_luigi_stats(db_path, table):
    """
    Loads the luigi stats from the stats db
    :param db_path: path to database
    :return: DataFrame
    """
    engine = _make_engine(db_path)
    return pd.read_sql_table(table, engine)
