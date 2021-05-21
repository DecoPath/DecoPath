# -*- coding: utf-8 -*-

"""Constants module."""

import os
import time

HERE = os.path.dirname(os.path.realpath(__file__))
APP_DIR = os.path.abspath(os.path.join(HERE, os.pardir))
PROJECT_DIR = os.path.abspath(os.path.join(APP_DIR, os.pardir))

# Logs path
LOGS_DIR = os.path.join(PROJECT_DIR, 'logs')

# JSON files
PATHWAY_HIERARCHY_JSON = os.path.join(APP_DIR, 'static', 'json', 'pathway_hierarchy.json')
HGNC_MAPPINGS = os.path.join(APP_DIR, 'static', 'json', 'hgnc_mappings.json')

# Decopath Ontology
DECOPATH_ONTOLOGY = os.path.join(APP_DIR, 'static', 'decopath_ontology', 'decopath_ontology.xlsx')


def make_log_dir():
    """Ensure that log directory exists."""
    os.makedirs(LOGS_DIR, exist_ok=True)


make_log_dir()

TODAY = time.strftime("%d_%m_%Y")

#: Round floats to 2 decimals
FLOAT_ROUND = 2

#: Significance value
P_VALUE = 0.05
PADJ = 0.05
SIGNIFICANCE_VALUE = 0.05

#: Number of HGNC gene symbols
GENE_UNIVERSE = 42345

#: Databases
KEGG = 'kegg'
REACTOME = 'reactome'
WIKIPATHWAYS = 'wikipathways'
PATHBANK = 'pathbank'
DECOPATH = 'decopath'
CUSTOM = 'custom'

DATABASES = {
    KEGG: 'KEGG',
    REACTOME: 'Reactome',
    WIKIPATHWAYS: 'WikiPathways',
    PATHBANK: 'PathBank',
    DECOPATH: 'DecoPath',
    CUSTOM: 'Custom'
}

#: DecoPath run options
RESULTS = 'results'
GSEA = 'gsea'
ORA = 'ora'

#: Database pairwise mappings
KEGG_REACTOME = 'kegg_reactome'
KEGG_WIKIPATHWAYS = 'kegg_wikipathways'
REACTOME_WIKIPATHWAYS = 'reactome_wikipathways'
PATHBANK_KEGG = "pathbank_kegg"
PATHBANK_REACTOME = "pathbank_reactome"
PATHBANK_WIKIPATHWAYS = "pathbank_wikipathways"

#: Maximum number of jobs a non-staff user can have
MAX_JOBS = 2

"""GMT files"""

GMT_FILES_DIR = os.path.join(APP_DIR, 'static', 'gmt_files')
GMT_FILE_EXTENSION = '.gmt'


def make_geneset_dir():
    """Ensure that geneset directory exists."""
    os.makedirs(GMT_FILES_DIR, exist_ok=True)


make_geneset_dir()

REACTOME_GMT = 'reactome.gmt'
REACTOME_MAPPINGS_FILE = 'reactome_mappings.gmt'

PATHWAY_FORTE_BASE_URL = "https://raw.githubusercontent.com/pathwayforte/pathway-forte/master/data/gmt_files"

#: GMT files URLs
DEFAULT_DATABASES_URLS = {
    KEGG: f"{PATHWAY_FORTE_BASE_URL}/kegg.gmt",
    REACTOME: f"{PATHWAY_FORTE_BASE_URL}/{REACTOME_GMT}",
    WIKIPATHWAYS: f"{PATHWAY_FORTE_BASE_URL}/wikipathways.gmt",
    PATHBANK: f"{PATHWAY_FORTE_BASE_URL}/pathbank.gmt",
    DECOPATH: f"{PATHWAY_FORTE_BASE_URL}/decopath.gmt",
}

"""PathBank downloads"""

PATHBANK_BASE_URL = "https://pathbank.org/downloads/"
PATHBANK_PATHWAY_URL = f"{PATHBANK_BASE_URL}'pathbank_all_pathways.csv.zip"
PATHBANK_PROTEIN_URL = f"{PATHBANK_BASE_URL}'pathbank_all_pathways.csv.zip'"

PATHBANK_PATHWAYS_FILE = "pathbank_pathways.csv"

"""JSON files with ID to name mappings"""

NAME_ID_MAPPINGS = os.path.join(APP_DIR, 'static', 'json', 'name_id_mappings')
PATHWAY_NAMES = os.path.join(NAME_ID_MAPPINGS, 'pathway_names.json')
DECOPATH_NAME_MAPPINGS = os.path.join(NAME_ID_MAPPINGS, 'decopath.json')

"""ComPath mappings"""

COMPATH_BASE_URL = "https://raw.githubusercontent.com/ComPath/compath-resources/master/mappings"

#: ComPath mapping URLs
COMPATH_MAPPING_URLS = {
    KEGG_REACTOME: f"{COMPATH_BASE_URL}/kegg_reactome.csv",
    KEGG_WIKIPATHWAYS: f"{COMPATH_BASE_URL}/kegg_wikipathways.csv",
    REACTOME_WIKIPATHWAYS: f"{COMPATH_BASE_URL}/wikipathways_reactome.csv",
    PATHBANK_KEGG: f"{COMPATH_BASE_URL}/pathbank_kegg.csv",
    PATHBANK_REACTOME: f"{COMPATH_BASE_URL}/pathbank_reactome.csv",
    PATHBANK_WIKIPATHWAYS: f"{COMPATH_BASE_URL}/pathbank_wikipathways.csv"
}

"""Output files with results for GSEA"""


def make_gsea_export_directories():
    """Ensure that gsea results export directories exist."""
    os.makedirs(RESULTS, exist_ok=True)


GSEA_RESULTS = os.path.join(PROJECT_DIR, RESULTS)
DECOPATH_GMT = os.path.join(GMT_FILES_DIR, 'decopath.gmt')
DECOPATH_CSV = os.path.join(PROJECT_DIR, 'viewer', 'static', 'gmt_files', 'decopath.csv')
CONCATENATE_GMT = os.path.join(PROJECT_DIR, 'viewer', 'static', 'dummy_files', 'concatenate.gmt')
TEST_GENESET = os.path.join(PROJECT_DIR, 'test_files', 'gsea/test_geneset.gmt')
TEST_FILES = os.path.join(PROJECT_DIR, 'test_files')

"""Mapping file columns and values"""

#: Equivalent mapping
IS_PART_OF = "BFO:0000050"
EQUIVALENT_TO = "skos:exactMatch"
MAPPING_TYPE = 'Mapping Type'
SOURCE_ID = 'Source ID'
TARGET_ID = 'Target ID'
SOURCE_NAME = 'Source Name'
TARGET_NAME = 'Target Name'
SOURCE_RESOURCE = 'Source Resource'
TARGET_RESOURCE = 'Target Resource'

"""Enrichment Parameters"""

#: Default upper and lower limit of genes to run GSEA
DEFAULT_MAX_GSEA = 500
DEFAULT_MIN_GSEA = 10

#: Default upper and lower limit of genes to run ORA
DEFAULT_MAX_ORA = 1000
DEFAULT_MIN_ORA = 10

#: Number of permutations to run GSEA
PERMUTATION_NUM = (
    ("100", "default"),
    ("10", "10"),
    ("100", "100"),
    ("1000", "1000"),
)

#: Type of permutations to run GSEA
PERMUTATION_TYPE = (
    ("phenotype", "default"),
    ("phenotype", "phenotype"),
    ("gene_set", "gene set"),
)

#: Pathway enrichment methods
ENRICHMENT_METHOD = (
    (GSEA, "GSEA"),
    (ORA, "ORA")
)

#: GSEA methods
METHOD = (
    ("log2_ratio_of_classes", "default"),
    ("log2_ratio_of_classes", "log2 ratio of classes"),
    ("signal_to_noise", "signal to noise"),
    ("t_test", "t test"),
    ("ratio_of_classes", "ratio of classes"),
    ("diff_of_classes", "difference of classes")
)

#: GSEA methods
METHOD_MAPPING = {
    "log2_ratio_of_classes": "log2 ratio of classes",
    "signal_to_noise": "signal to noise",
    "t_test": "t test",
    "ratio_of_classes": "ratio of classes",
    "diff_of_classes": "difference of classes"
}

#: Status of current job
STATUS_CODE_MAPPING = {
    0: '<span class="fas fa-times"  style="color: #FF0000"></span>',  # failure
    1: '<span class="fas fa-redo" style="color: #fce702"></span>',  # running
    2: '<span class="fas fa-check" style="color: #49e52d"></span>',  # success
}

"""Display results"""

#: Summary table supplemental information
RESULTS_METADATA_HEADER = [
    'Experiment',
    'Status code',
    'Error message',
    'Status',
    'Results',
    'Consensus',
    'Explore',
    'Analysis',
    'Databases',
    'Dataset',
    'Classes',
    'Class file',
    'Total samples',
    'Method',
    'Maximum genes',
    'Minimum genes',
    'Significance threshold',
    'Permutation type',
    'Number of permutations',
    'Significance threshold log2 fold changes',
    'Differential gene expression file',
]

"""Visualization"""

COLORMAP_VALUES = [
    (0.2519971417644415, 0.4987337088076726, 0.5751602783606602),
    (0.43026136111758173, 0.6200066482697917, 0.6787801878373952),
    (0.6085255804707219, 0.7412795877319109, 0.7824000973141302),
    (0.786789799823862, 0.86255252719403, 0.8860200067908652),
    (0.95, 0.95, 0.95),
    (0.954577257933482, 0.7665309859226215, 0.7803256889894359),
    (0.9197182699854205, 0.5873587656270927, 0.6117400023569117),
    (0.884859282037359, 0.40818654533156384, 0.4431543157243877),
    (0.8510408608937171, 0.23436274952246883, 0.2796010376480583),
]

"""Form formats"""

#: Separators
CSV = 'csv'
TSV = 'tsv'

"""User submitted results columns"""

#: User submitted GSEA results expected header
USER_RESULTS_GSEA_COLUMN_NAMES = {
    'es',
    'nes',
    'p_value',
    'q_value',  # adjusted_p_value? FDR?
}

#: User submitted ORA results expected header
USER_RESULTS_ORA_COLUMN_NAMES = {
    'p_value',
    'q_value',
}

"""GSEA results column from gseapy"""

#: GSEA run results header
DEFAULT_GSEA_RESULTS_COLUMN_NAMES = {
    'es', 'nes', 'pval', 'fdr', 'geneset_size', 'matched_size', 'genes', 'ledge_genes'
}

#: Fold changes header
USER_FOLD_CHANGES_COLUMN_NAME = {
    'log2fc', 'q_value'
}

GENE_SYMBOL = 'gene_symbol'

CONSENSUS_MAPPINGS = {
    0: 'Discordant',
    1: 'no-mappings',
    2: 'Concordant'
}
