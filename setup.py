import glob
from setuptools import setup

setup(
    name='cat2',
    version='2.0',
    packages=['cat', 'tools'],
    python_requires='>=3.7.0',
    # When using conda (environment.yaml), install with: pip install -e . --no-deps
    # so pip does not re-resolve these packages on top of conda builds.
    install_requires=[
        'toil>=7.0.0',
        'seaborn>=0.13.0',
        'pandas>=2.0.0',
        'frozendict>=2.4.1',
        'configobj>=5.0',
        'sqlalchemy>=2.0',
        'ete3>=3.0',
        'pysam>=0.19.1',
        'numpy>=2.2.0',
        'scipy>=1.13.3',
        'bx-python>=0.13',
        'gffutils>=0.10',
        'biopython>=1.76',
        'snakemake>=8.0',
    ],
    entry_points={'console_scripts': ['cat2=cat.cli:main']},
    scripts=['programs/cat_to_ncbi_submit', 'programs/translate_gene_pred',
             'programs/validate_gff3', 'programs/cat_parse_ncbi_genbank',
             'programs/cat_parse_ncbi_refseq', 'programs/cat_parse_prokka_gff3'],
    data_files=[
        ('share/cat', (
            glob.glob('augustus_cfgs/*.cfg') +
            ['standalones/vi2-7k.kan', 'standalones/vi2-7k.kan.cali']
        )),
        ('share/cat2', ['Snakefile']),
    ],
    author='Ian Fiddes, Prajna Hebbar',
    description='Comparative Annotation Toolkit',
    url='https://github.com/ComparativeGenomicsToolkit/Comparative-Annotation-Toolkit',
    license='Apache 2.0',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'License :: OSI Approved :: Apache Software License',
        'Programming Language :: Python :: 3',
    ],
    keywords='bioinformatics comparative genomics annotation',
)
