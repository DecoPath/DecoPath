# -*- coding: utf-8 -*-

"""Utilities module."""

import functools
import io
import json
import logging
import operator
import random
import sqlite3
import string
import tempfile
import urllib.request
import zipfile
from collections import defaultdict
from os.path import isfile
from pathlib import Path
from typing import Any, Dict, List, Tuple, Iterable, Union

import matplotlib
import matplotlib.cm
import matplotlib.colors
import networkx as nx
import pandas as pd
import requests
from bio2bel_reactome import Manager as ReactomeManager
from django.core.files.uploadedfile import UploadedFile
from django.core.files.uploadhandler import FileUploadHandler
from networkx import DiGraph
from networkx.readwrite.json_graph import tree_data

from viewer.models import PathwayDatabase, User, Pathway
from viewer.src.constants import *

logger = logging.getLogger(__name__)


def handle_file_download(url: str, filename: str):
    """Handle file download from url, write to file and save to static directory."""
    if isfile(filename) is False:
        response = io.StringIO(urllib.request.urlopen(url).read().decode('utf-8'))
        with open(filename, 'w') as f:
            f.write(response.read())


def read_data_file(file_path: str, filename: str) -> Union[pd.DataFrame, str]:
    """Check read data file."""
    logger.info(f"Reading {file_path}")
    try:
        if file_path.endswith(CSV):
            return pd.read_csv(file_path, sep=",")

        elif file_path.endswith(TSV):
            return pd.read_csv(file_path, sep="\t")

    except IOError:
        logger.error(f"Failed to read {filename} {file_path}. File exists: {os.path.isfile(file_path)}")
        return f'There is a problem with your file {filename}. please check that it meets the criteria.'


@functools.lru_cache()
def _get_hgnc_mapping_dict(hgnc_mappings=HGNC_MAPPINGS):
    """Load HGNC name-id mappings."""
    with open(hgnc_mappings) as json_file:
        return json.load(json_file)


def concatenate_files(files_list: List, databases_list: List) -> str:
    """Concatenate GMT files in a list of files and write to a new file."""
    databases = "_".join(str(x) for x in sorted(databases_list))
    concatenated_file = databases + GMT_FILE_EXTENSION
    file_path = os.path.join(GMT_FILES_DIR, concatenated_file)

    # Check if GMT file already exists
    if not isfile(file_path):

        # Concatenate list of text files
        with open(file_path, 'w') as outfile:
            for file in files_list:
                with open(file) as infile:
                    for line in infile:
                        outfile.write(line)

    return file_path


def spliterate(lines: Iterable[str], sep='\t') -> Iterable[Tuple[str, ...]]:
    """Split each line in the iterable by the given separator."""
    for line in lines:
        yield line.strip().split(sep)


def query_pathway_db_model():
    """Query pathway database model for all databases to display as multiple choice field for user selection."""
    return tuple([
        (index, DATABASES[str(obj)])
        for index, obj in enumerate(PathwayDatabase.objects.filter(
            database_name__in=[KEGG, REACTOME, PATHBANK, WIKIPATHWAYS]
        ))
    ])


def get_database_by_id(row: pd.Series):
    """Get database by identifier suffix if default databases used for analysis."""
    if row.startswith('hsa'):
        return KEGG
    elif row.startswith('PW'):
        return PATHBANK
    elif row.startswith('R-HSA'):
        return REACTOME
    elif row.startswith('WP'):
        return WIKIPATHWAYS
    elif row.startswith('DC'):
        return DECOPATH

    return CUSTOM


def get_missing_columns(column_labels: set, df: pd.DataFrame) -> List:
    """Get missing column labels."""
    return [
        column
        for column in column_labels
        if column not in df
    ]


def get_equivalent_dicts(equivalent_rows):
    """Return equivalent pathways."""
    equivalent_pathways = defaultdict(list)
    for source_db, source_id, source_name, mapping, target_db, target_id, target_name in equivalent_rows.values:
        equivalent_pathways[(source_db, source_id)].append((target_db, target_id))
        equivalent_pathways[(target_db, target_id)].append((source_db, source_id))

    # Warning if the number of mappings is not the same between two pathways (it should always be)
    for source_db, source_id, source_name, mapping, target_db, target_id, target_name in equivalent_rows.values:
        if (
            len(equivalent_pathways[(source_db, source_id)]) != len(equivalent_pathways[(target_db, target_id)])
        ):
            raise ValueError(
                f'{(source_db, source_id)} has {len(equivalent_pathways[(source_db, source_id)])} mappings'
                f'{(target_db, target_id)} has {len(equivalent_pathways[(target_db, target_id)])} mappings'
            )

    return dict(equivalent_pathways)


def _check_duplicates(df):
    """Check quality of df."""
    hierarchy_df = df[df[MAPPING_TYPE] == IS_PART_OF]

    duplicates_hierarchy = [
        df.iloc[i]
        for i, duplicated in hierarchy_df[SOURCE_ID].duplicated().items()
        if duplicated
    ]

    if duplicates_hierarchy:
        raise ValueError(f'Duplicate hierarchy: {duplicates_hierarchy}')


def _add_gsea_cmap(nodes, node_to_score, fdr_dict, significance_value):
    """Return dictionary with the nodes and their corresponding normalized colors."""
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("DecoPath ColorMap", COLORMAP_VALUES)
    norm = matplotlib.colors.Normalize(
        min(node_to_score.items(), key=operator.itemgetter(1))[1],  # min value
        max(node_to_score.items(), key=operator.itemgetter(1))[1],  # max value
    )

    # Return the colors for each node after doing a normalization step
    color_map = {}
    for node_id in nodes:
        if node_id in node_to_score:
            if fdr_dict[node_id] < significance_value:
                color_map[node_id] = matplotlib.colors.rgb2hex(
                    cmap(norm(float(node_to_score[node_id])))
                )
            else:
                color_map[node_id] = '#b4b4b4'  # grey color for nodes that don't pass significance threshold

        else:
            color_map[node_id] = '#dadada'  # grey color for not mapped nodes

    return color_map


def _add_pvalue_cmap(nodes, node_to_score, significance_value):
    """Return dictionary with the nodes and their corresponding normalized colors according to p value."""
    cmap = matplotlib.cm.get_cmap('Reds').reversed()
    norm = matplotlib.colors.Normalize(
        -0.01,  # min value
        0.05,  # max value
    )

    # Return the colors for each node after doing a normalization step
    color_map = {}

    for node_id in nodes:
        if node_id in node_to_score and node_to_score[node_id] < significance_value:
            color_map[node_id] = matplotlib.colors.rgb2hex(
                cmap(norm(float(node_to_score[node_id])))
            )
        else:
            color_map[node_id] = '#94989c'  # grey color for not mapped nodes

    return color_map


def _label_tree(network, id_to_database, id_to_name):
    """Add attributes to the nodes in the tree."""
    # Set attributes
    nx.set_node_attributes(network, id_to_database, 'database')
    nx.set_node_attributes(
        network,
        {
            node_id: id_to_name[node_id] if node_id in id_to_name else ''  # TODO: add something or not?
            for node_id in network.nodes()
        },
        'name',
    )

    # Remove KEGG prefixes from the canonical name
    mapping = {
        node_name: node_name.replace("path:", "")
        for node_name in network
        if "path:" in node_name
    }
    nx.relabel_nodes(network, mapping, copy=False)


def parse_hierarchy_excel(url):
    """Parse hierarchy file."""
    xls = pd.ExcelFile(
        DECOPATH_ONTOLOGY,
    )

    # Initialize nested dictionary for d3, then recursively iterate through tree
    root_node = 'SuperPathway'
    tree_network = DiGraph()

    #: Dictionaries to store later the node attributes
    id_to_database = {}
    id_to_name = {}

    for sheet_name in xls.sheet_names:
        if sheet_name == 'equivalence_same_db':
            continue

        # Get the dataframe for the given sheet
        df = pd.read_excel(
            io=xls,
            sheet_name=sheet_name,
            usecols=[
                SOURCE_RESOURCE, SOURCE_ID, SOURCE_NAME,
                MAPPING_TYPE,
                TARGET_RESOURCE, TARGET_ID, TARGET_NAME,
            ],
            dtype=str,
        )

        # Get equivalent pathways
        equivalent_pathways_dict = get_equivalent_dicts(df[df[MAPPING_TYPE] == EQUIVALENT_TO])

        #: Check quality of the excel
        _check_duplicates(df)

        # Filter

        for source_db, source_id, source_name, mapping, target_db, target_id, target_name in df.values:

            # Skip equivalentTo and empty lines
            if mapping == EQUIVALENT_TO or pd.isna(mapping):
                continue

            id_to_name[target_id] = target_name
            id_to_database[target_id] = target_db

            id_to_name[source_id] = source_name
            id_to_database[source_id] = source_db

            """Logic to generate the tree"""
            # Create super pathway node
            if mapping == 'SuperPathway':
                tree_network.add_edge(root_node, target_id)

            # Create hierarchy
            elif mapping == IS_PART_OF:
                # Create the hierarchical relation

                tree_network.add_edge(target_id, source_id)

            else:
                raise ValueError(f'invalid {mapping}')

    # Add properties to the nodes in the network
    _label_tree(tree_network, id_to_database, id_to_name)

    return tree_data(tree_network, root_node), tree_network, equivalent_pathways_dict, root_node


def map_results_to_hierarchy(
    network_hierarchy,
    root_node,
    results: List[Tuple[str, str, str, str, str]],
    enrichment_method: str,
    significance_value: float,
) -> Dict[Any, Any]:
    """Map results of GSEA experiment to hierarchy.

    :param network_hierarchy: network reprsenting the hierarchy
    :param root_node: hierarchy root
    :param results: results of gsea experiment
    :param enrichment_method: enrichment method
    :return: hierarchy with mapped results
    """
    geneset_size_mapping = {}
    node_to_score = {}
    #: only used for GSEA
    fdr_dict = {}

    # Iterate through nodes in the pathway hierarchy graph and annotate them
    for db, pathway_id, score, fdr, geneset_size in results:

        if pathway_id not in network_hierarchy:
            continue

        if enrichment_method.lower() == ORA:
            node_to_score[pathway_id] = fdr

        elif enrichment_method.lower() == GSEA:
            node_to_score[pathway_id] = score
            fdr_dict[pathway_id] = fdr

        else:
            logger.warning(f'unknown enrichment method {enrichment_method}')
            raise ValueError()
        geneset_size_mapping[pathway_id] = geneset_size

    if not node_to_score:
        raise ValueError('could not map to any pathway in the hierarchy')

    elif all(i is None for i in node_to_score.values()):
        logger.warning('error with pathway scores')
        raise ValueError()

    # Add the attributes
    nx.set_node_attributes(network_hierarchy, geneset_size_mapping, 'geneset_size')
    nx.set_node_attributes(network_hierarchy, node_to_score, 'direction')

    if enrichment_method.lower() == GSEA:
        nx.set_node_attributes(network_hierarchy, fdr_dict, 'fdr')
        nx.set_node_attributes(
            network_hierarchy,
            _add_gsea_cmap(network_hierarchy.nodes(), node_to_score, fdr_dict, significance_value),
            'color',
        )
    #: ORA
    else:
        nx.set_node_attributes(
            network_hierarchy,
            _add_pvalue_cmap(network_hierarchy.nodes(), node_to_score, significance_value),
            'color',
        )

    return tree_data(network_hierarchy, root_node)


def export_geneset(geneset_dict, database, temp_file, gmt_file):
    """Export gene set to gmt file format."""
    df = pd.DataFrame.from_dict(data=geneset_dict, orient='index')

    # Add resource column
    df['Resource'] = pd.Series(database, index=df.index)
    df = df[['Resource'] + [col for col in df.columns if col != 'Resource']]

    df.to_csv(temp_file, header=False, sep='\t')

    _export_geneset_to_gmt(temp_file, gmt_file)

    # Remove intermediary file
    os.remove(temp_file)


def _export_geneset_to_gmt(geneset_file, outfile):
    """Export gene set to gmt file format."""
    with open(geneset_file, 'r') as file:
        with open(outfile, 'w') as f:
            for line in file:
                line = line.rstrip()
                f.write(line + '\n')

    f.close()


def get_reactome_mappings():
    query_reactome = Pathway.objects.filter(pathway_database='reactome')

    reactome_ids = [
        obj.pathway_id
        for obj in query_reactome
    ]

    gmt_dict = _get_gmt_dict(os.path.join(GMT_FILES_DIR, 'reactome.gmt'))

    reactome_gmt = {
        k: v
        for k, v in gmt_dict.items()
        if k in reactome_ids
    }

    export_geneset(reactome_gmt, 'reactome', 'reactome_mappings.csv', 'reactome_mappings.gmt')


def get_equivalent_pathway_dc_ids(decopath_ontology):
    """Parse DecoPath ontology file and get DC IDs for equivalent super pathways."""
    sheets_dict = pd.read_excel(
        io=decopath_ontology,
        engine='openpyxl',
        sheet_name=None,
        usecols=[
            SOURCE_RESOURCE, SOURCE_ID, SOURCE_NAME,
            MAPPING_TYPE,
            TARGET_RESOURCE, TARGET_ID, TARGET_NAME,
        ],
        dtype=str,
        index_col=None,
    )

    sheets_dict.pop('equivalence_same_db', None)

    frames = [sheet for name, sheet in sheets_dict.items()]
    df = pd.concat(frames)

    # Remove kegg prefixes
    df[SOURCE_ID] = df[SOURCE_ID].str.replace('path:', '')
    df[TARGET_ID] = df[TARGET_ID].str.replace('path:', '')

    # Get equivalent pathways
    equivalence_df = df.loc[df[MAPPING_TYPE] == EQUIVALENT_TO]
    equivalent_pathways = equivalence_df[SOURCE_ID].to_list() + equivalence_df[TARGET_ID].to_list()

    # Get DC IDs for equivalent pathway
    id_to_dc_id = pd.DataFrame.from_dict({
        source_id: {'pathway_id': source_id, 'dc_id': target_id, 'dc_name': target_name}
        for source_db, source_id, source_name, mapping_type, target_db, target_id, target_name in df.values
        if mapping_type == IS_PART_OF
        if source_id in equivalent_pathways
    }, orient='index')

    id_to_dc_id = id_to_dc_id.reset_index(drop=True)

    return id_to_dc_id, df


def get_dc_pathway_resources(decopath_ontology):
    """Parse DecoPath ontology file and get source resources for DC super pathways."""
    sheets_dict = pd.read_excel(
        io=decopath_ontology,
        engine='openpyxl',
        sheet_name=None,
        usecols=[
            SOURCE_RESOURCE, SOURCE_ID, SOURCE_NAME,
            MAPPING_TYPE,
            TARGET_RESOURCE, TARGET_ID, TARGET_NAME,
        ],
        dtype=str,
        index_col=None,
    )

    sheets_dict.pop('equivalence_same_db', None)

    frames = [sheet for name, sheet in sheets_dict.items()]
    df = pd.concat(frames)

    # Remove kegg prefixes
    df[SOURCE_ID] = df[SOURCE_ID].str.replace('path:', '')
    df[TARGET_ID] = df[TARGET_ID].str.replace('path:', '')

    # Get equivalent pathways
    equivalence_df = df.loc[df[MAPPING_TYPE] == EQUIVALENT_TO]
    equivalent_pathways = equivalence_df[SOURCE_ID].to_list() + equivalence_df[TARGET_ID].to_list()

    # Get pathway resources for DC ID
    dc_source_dict = defaultdict(set)
    dc_source_id_dict = defaultdict(set)
    for source_db, source_id, source_name, mapping_type, target_db, target_id, target_name in df.values:

        if mapping_type == IS_PART_OF:
            if source_id in equivalent_pathways:
                dc_source_dict[target_id].add(source_db)
                dc_source_id_dict[target_id].add(source_id)

    return dc_source_dict, dc_source_id_dict


def _get_gmt_dict(filename):
    """Parse gmt files and get gene sets."""
    with open(filename, 'r') as f:
        content = [line.strip().split('\t') for line in f]

    return {pathway[0]: pathway[2:] for pathway in content}


def get_decopath_genesets(decopath_ontology, gmt_dir: str):
    """Generate DecoPath gene sets with super pathways."""
    concatenated_genesets_dict = {}
    dc_mapping = defaultdict(list)

    if not os.path.isdir(gmt_dir):
        make_geneset_dir()

    super_pathway_mappings, ontology_df = get_equivalent_pathway_dc_ids(decopath_ontology)

    # Get gene sets from individual databases
    gmt_files = [os.path.join(GMT_FILES_DIR, filename) for filename in os.listdir(gmt_dir) if filename.endswith('.gmt')]

    genesets = [_get_gmt_dict(file) for file in gmt_files]

    # Concatenate gene sets from individual databases
    for geneset in genesets:
        concatenated_genesets_dict.update(geneset)

    # Get super pathway gene sets with DecoPath IDs
    for pathway_id, dc_id, dc_name in super_pathway_mappings.values:

        if pathway_id in concatenated_genesets_dict:
            dc_mapping[dc_id].append(concatenated_genesets_dict[pathway_id])

    # Return DecoPath gene sets
    return {
        pathway_id: {gene for sublist in geneset for gene in sublist}
        for pathway_id, geneset in dc_mapping.items()
    }


def export_to_json(obj, outfile):
    with open(outfile, "w") as f:
        json.dump(obj, f)


def get_decopath_id_name_mappings(outfile, decopath_ontology: str = DECOPATH_ONTOLOGY):
    sheets_dict = pd.read_excel(
        io=decopath_ontology,
        engine='openpyxl',
        sheet_name=None,
        usecols=[
            SOURCE_RESOURCE, SOURCE_ID, SOURCE_NAME,
            MAPPING_TYPE,
            TARGET_RESOURCE, TARGET_ID, TARGET_NAME,
        ],
        dtype=str,
        index_col=None,
    )

    sheets_dict.pop('equivalence_same_db', None)

    frames = [sheet for name, sheet in sheets_dict.items()]
    df = pd.concat(frames)

    decopath_mapping_dict = {}

    for source_db, source_id, source_name, mapping_type, target_db, target_id, target_name in df.values:
        if source_db == DECOPATH:
            decopath_mapping_dict[source_id] = source_name
        if target_db == DECOPATH:
            decopath_mapping_dict[target_id] = target_name

    export_to_json(decopath_mapping_dict, outfile)


def get_reactome_id_name_mappings(outfile):
    reactome_manager = ReactomeManager()

    reactome_mappings = reactome_manager.get_pathway_names_to_ids()

    reactome_mappings_invert = {v: k for k, v in reactome_mappings.items()}

    export_to_json(reactome_mappings_invert, outfile)


def get_name_id_mapping(db_file, pathway, outfile):
    """Create a database connection to the SQLite database specified by the db_file and get pathway name/ID mappings."""
    conn = sqlite3.connect(db_file)

    cur = conn.cursor()
    cur.execute("SELECT * FROM " + pathway)
    rows = cur.fetchall()
    name_id_mappings = {identifier: name for _, identifier, name in rows}

    export_to_json(name_id_mappings, outfile)


def handle_zipfile_download(path: str) -> None:
    """Download and extract zip file content."""
    r = requests.get(path)
    z = zipfile.ZipFile(io.BytesIO(r.content))
    z.extractall()


def export_pathbank_name_id_mappings(outfile, pathbank_pathways_url: str = PATHBANK_PATHWAY_URL):
    """Download PathBank pathway metadata and export as json."""

    logger.info('Downloading PathBank content.')

    # Get PathBank pathways
    handle_zipfile_download(pathbank_pathways_url)

    df = pd.read_csv(PATHBANK_PATHWAYS_FILE, sep=',')

    pathbank_mappings = pd.Series(df['Name'].values, index=df['PW ID']).to_dict()

    export_to_json(pathbank_mappings, outfile)


def get_pathway_names(outfile, mapping_dir: str = NAME_ID_MAPPINGS):
    name_id_mappings = []

    for filename in os.listdir(mapping_dir):
        if '.json' in filename:
            file_path = os.path.join(os.getcwd(), mapping_dir, filename)
            with open(file_path) as f:
                data = json.load(f)

        name_id_mappings.append(data)

    id_name_dict = {
        identifier: name
        for pathway_dict in name_id_mappings
        for identifier, name in pathway_dict.items()
    }

    with open(outfile, "w") as f:
        json.dump(id_name_dict, f)


class TransientFileUploadHandler(FileUploadHandler):
    """
    Upload handler that streams data into a temporary file.
    """
    def new_file(self, *args, **kwargs):
        """
        Create the file object to append to as data is coming in.
        """
        super().new_file(*args, **kwargs)
        self.file = TransientUploadedFile(self.file_name, self.content_type, 0, self.charset, self.content_type_extra)

    def receive_data_chunk(self, raw_data, start):
        self.file.write(raw_data)

    def file_complete(self, file_size):
        self.file.seek(0)
        self.file.size = file_size
        return self.file


class TransientUploadedFile(UploadedFile):
    """
    A file uploaded to a temporary location (i.e. stream-to-disk).
    """
    def __init__(self, name, content_type, size, charset, content_type_extra=None):
        _, ext = os.path.splitext(name)
        temp_name = ''.join(random.choices(string.ascii_letters + string.digits, k=16))
        Path(os.path.join(tempfile.gettempdir(), f'{temp_name}.upload{ext}')).touch()
        file = open(os.path.join(tempfile.gettempdir(), f'{temp_name}.upload{ext}'), mode='w+b')
        super().__init__(file, name, content_type, size, charset, content_type_extra)

    def transient_file_path(self):
        """Return the full path of this file."""
        return self.file.name

    def close(self):
        try:
            return self.file.close()
        except FileNotFoundError:
            # The file was moved or deleted before the tempfile could unlink
            # it. Still sets self.file.close_called and calls
            # self.file.file.close() before the exception.
            pass


def del_transient_files(*args):
    for file in args:
        os.remove(file)


def del_user(user_email):
    User.objects.get(email__exact=user_email).delete()
    return 1
