# -*- coding: utf-8 -*-

"""Handle results module."""
import json
import logging
from collections import defaultdict
from itertools import combinations
from typing import List, Dict, Union

import pandas as pd
from django.core.exceptions import ObjectDoesNotExist

from viewer.models import EnrichmentResult, Pathway
from viewer.src.constants import *
from viewer.src.results_utils import (
    clean_none_values,
    check_consensus,
    round_float,
    process_database_names,
    get_dc_info_gsea,
    get_consensus_gsea, get_dc_info_ora, get_database_by_pathway_id,
)
from viewer.src.utils import get_database_by_id, spliterate

logger = logging.getLogger(__name__)


def create_summary_table(current_user):
    """Generate experiment summary table for user with link to results for a given experiment."""
    if current_user.is_staff:
        query_results = EnrichmentResult.objects.all()
    else:
        query_results = EnrichmentResult.objects.filter(user=current_user)

    # Get results and experiment metadata summary table on accounts page with links to results pages
    table_body = []
    objs = []

    for index, obj in enumerate(query_results):

        enrichment_method = obj.enrichment_method

        for class_label in obj.get_phenotypes():

            if isinstance(class_label, int):
                classes = '|'.join(str(label) for label in obj.get_phenotypes())
            elif len(class_label) == 1:
                for label in class_label:
                    classes = label
            else:
                classes = '|'.join(obj.get_phenotypes())

        if obj.result_status != 2:

            if enrichment_method == ORA:
                results_link = f'<a href="/results_{ORA}/{obj.result_id}" class="button disabled">Load Results</a>'
                consensus_link = f'<a href="/consensus_{ORA}/{obj.result_id}" class="button disabled">Load Consensus Table</a>'
            else:
                results_link = f'<a href="/results_{GSEA}/{obj.result_id}" class="button disabled">Load Results</a>'
                consensus_link = f'<a href="/consensus_{GSEA}/{obj.result_id}" class="button disabled">Load Consensus Table</a>'

            circles_link = f'<a href="/viz/circles/{obj.result_id}" class="button disabled">Visualize Consensus</a>'

        else:

            if enrichment_method == ORA:
                results_link = f'<a href="/results_{ORA}/{obj.result_id}" class="button">Load Results</a>'
                consensus_link = f'<a href="/consensus_{ORA}/{obj.result_id}" class="button">Load Consensus Table</a>'
            else:
                results_link = f'<a href="/results_{GSEA}/{obj.result_id}" class="button">Load Results</a>'
                consensus_link = f'<a href="/consensus_{GSEA}/{obj.result_id}" class="button">Load Consensus Table</a>'

            circles_link = f'<a href="/viz/circles/{obj.result_id}" class="button">Visualize Consensus</a>'

        if enrichment_method != PRERANK:
            enrichment_method = obj.enrichment_method.upper()

        row = [
            index + 1,
            obj.result_status,
            obj.error_message,
            STATUS_CODE_MAPPING[obj.result_status],
            results_link,
            consensus_link,
            circles_link,
            enrichment_method,
            ' '.join([DATABASES[database] for database in sorted(obj.get_databases())]),
            obj.data_filename,
            classes,
            obj.class_filename,
            clean_none_values(obj.sample_number),
            obj.calculation_method,
            clean_none_values(obj.max_genes),
            clean_none_values(obj.min_genes),
            clean_none_values(obj.significance_threshold),
            obj.permutation_type,
            clean_none_values(obj.permutation_number),
            clean_none_values(obj.significance_threshold_fc),
            obj.fold_changes_filename,
        ]
        table_body.append(row)
        objs.append(obj)

    return table_body, objs


def get_ranking_table(results_df: pd.DataFrame):
    """Process results and display rankings."""
    # Get results table
    df = _get_results_table(results_df)

    df['Database'] = df['Database'].map(DATABASES)

    # Check if nes in columns if GSEA run
    if 'nes' in df.columns:

        df.nes = df.nes.round(2)
        df.es = df.es.round(2)

        # Check if results are uploaded by user
        if set(USER_RESULTS_GSEA_COLUMN_NAMES).issubset(set(df.columns)):

            df = df[['Identifier', 'Pathway', 'Database', 'nes', 'es', 'p_value', 'q_value']]
            df.rename(columns={
                'nes': 'NES',
                'es': 'ES',
                'p_value': 'p-value',
                'q_value': 'q-value',
            }, inplace=True)

            df['p-value'] = df['p-value'].map(round_float)
            df['q-value'] = df['q-value'].map(round_float)

            return df.to_html(
                classes='table table-hover table-bordered',
                table_id='ranking_table',
                justify='left',
                index=False,
                render_links=True,
                border=0,
            )

        # Results are from GSEA run
        else:
            df_dict = {}

            databases = sorted(list(df.Database.unique()))

            df = df[['Identifier', 'Pathway', 'Database', 'nes', 'es', 'pval', 'fdr', 'geneset_size', 'matched_size']]
            df.rename(columns={
                'pval': 'p-value',
                'fdr': 'q-value',
                'nes': 'NES',
                'es': 'ES',
                'geneset_size': 'Gene Set Size',
                'matched_size': 'Mapped Genes',
            }, inplace=True)

            df['p-value'] = df['p-value'].map(round_float)
            df['q-value'] = df['q-value'].map(round_float)

            # Get results table for each individual pathway database
            for db in databases:
                df_dict[f'{db}_ranking_table'] = df.loc[df['Database'] == db].to_html(
                    classes='table table-hover table-bordered datatable-display',
                    table_id=f'{db}',
                    justify='left',
                    index=False,
                    render_links=True,
                    border=0,
                )

            return df_dict, databases

    # ORA results
    else:
        df = df[['Identifier', 'Pathway', 'Database', 'p_value', 'q_value']]
        df.rename(columns={
            'p_value': 'p-value',
            'q_value': 'q-value',
        }, inplace=True)

    df['p-value'] = df['p-value'].map(round_float)
    df['q-value'] = df['q-value'].map(round_float)

    return df.to_html(
        classes='table table-hover table-bordered',
        table_id='ranking_table',
        justify='left',
        index=False,
        render_links=True,
        border=0,
    )


def _get_results_table(df: pd.DataFrame, json_file: str = PATHWAY_NAMES):
    """Generate generic results table."""
    with open(json_file) as f:
        pathway_id_name_dict = json.load(f)

    if 'pathway_id' in df:
        df['Pathway'] = df.pathway_id.map(pathway_id_name_dict)
        df['Database'] = df.pathway_id.map(get_database_by_id).to_list()

        df.rename(columns={'pathway_id': "Identifier"}, inplace=True)

    else:
        df['Pathway'] = df.index.map(pathway_id_name_dict)
        df['Database'] = df.index.map(get_database_by_id).to_list()

        df = df.rename_axis('Identifier').reset_index()

    return df


def get_genesets(databases: List, size=None):
    """Get genesets from gmt files for databases selected by user."""
    geneset_list = []

    # Get size of gene set for each pathway from user-selected databases
    for database in databases:
        path = os.path.join(GMT_FILES_DIR, f"{database}.gmt")

        with open(path) as file:
            # Get dictionary with pathway and corresponding gene set
            genesets_dict = {
                name: genes
                for name, _, *genes in spliterate(file)
            }

        if size is True:
            geneset_list.append({k: len(v) for k, v in genesets_dict.items()})

        else:
            geneset_list.append({k: v for k, v in genesets_dict.items()})

    geneset_dict = {
        pathway_id: geneset_size
        for geneset_sizes in geneset_list
        for pathway_id, geneset_size in geneset_sizes.items()
    }

    return geneset_dict


def get_results_circle_viz(results_df: pd.DataFrame, databases: List, enrichment_method: str):
    """Generate circles viz."""
    pathway_dict = defaultdict(dict)

    df = _get_results_table(results_df)

    geneset_size_dict = get_genesets(databases, size=True)

    if enrichment_method == ORA:
        # ensure same order in dataframe
        df = df[['Identifier', 'p_value', 'q_value', 'Pathway', 'Database']]

        for pathway_id, p_val, fdr, name, database in df.values:
            metadata = {'name': name, 'database': database, 'fdr': fdr}

            if pathway_id in geneset_size_dict:
                metadata['geneset_size'] = geneset_size_dict[pathway_id]

            pathway_dict[pathway_id] = metadata

    # method is GSEA or GSEA preranked
    else:
        if 'geneset_size' in df:
            df = df[['Identifier', 'Pathway', 'Database', 'nes', 'fdr', 'geneset_size']]

            for pathway_id, name, database, nes, fdr, geneset_size in df.values:
                metadata = {
                    'name': name,
                    'database': database,
                    'nes': nes,
                    'fdr': fdr,
                    'geneset_size': geneset_size
                }

                pathway_dict[pathway_id] = metadata

        else:
            df = df[['Identifier', 'Pathway', 'Database', 'nes', 'q_value']]

            for pathway_id, name, database, nes, fdr in df.values:
                metadata = {'name': name, 'database': database, 'nes': nes, 'fdr': fdr}

                if pathway_id in geneset_size_dict:
                    metadata['geneset_size'] = geneset_size_dict[pathway_id]

                pathway_dict[pathway_id] = metadata

    if enrichment_method == ORA:
        pathway_results = [
            (metadata['database'], pathway_id, None, metadata['fdr'], metadata['geneset_size'])
            for pathway_id, metadata in pathway_dict.items()
        ]

    else:
        pathway_results = [
            (metadata['database'], pathway_id, metadata['nes'], metadata['fdr'], metadata['geneset_size'])
            for pathway_id, metadata in pathway_dict.items()
        ]

    return pathway_results


def get_fold_changes_dict(df: pd.DataFrame, databases: List):
    fc_pathway_dict = defaultdict(dict)

    # Get genesets for user selected databases
    geneset_dict = get_genesets(databases)

    # TODO: get user defined q-value to filter out genes that do not pass significance
    df = df[df.padj < PADJ]
    fc_results_dict = df.set_index('gene_symbol').to_dict('index')

    fold_changes = dict(zip(df.gene_symbol, df.log2FoldChange))

    # Get fold changes for genes in genesets
    for pathway_id, geneset in geneset_dict.items():

        gene_fc_dict = {}

        for gene in geneset:

            for gene_symbol, metadata_fc_dict in fc_results_dict.items():

                if gene == gene_symbol:
                    gene_fc_dict[gene] = metadata_fc_dict['log2FoldChange']

        if any(gene_fc_dict):
            fc_pathway_dict[pathway_id] = gene_fc_dict

    return fc_pathway_dict, fold_changes


def process_overlap_for_venn_diagram(
    pathway_id: str,
    databases: List,
) -> Union[List[Dict], str]:
    """Calculate gene sets overlaps and process structure to render venn diagram."""
    # Creates future js array with gene sets' lengths
    overlaps_venn_diagram = []
    equivalent_pathways = set()
    equivalent_genesets = {}
    id_name_mappings = {}

    pathway_to_index = {}
    index = 0

    geneset_dict = get_genesets(databases)

    # Query pathway mappings model to filter pathways and retain only those with mappings
    if pathway_id.startswith('DC'):
        try:
            queryset = Pathway.objects.filter(decopath_id=pathway_id)
        except ValueError:
            return f'The pathway ID {pathway_id} does not have any equivalent pathways in the selected pathway' \
                   f' databases.'
    else:
        try:
            queryset = Pathway.objects.filter(pathway_id=pathway_id)
        except ValueError:
            return f'The pathway ID {pathway_id} does not have any equivalent pathways in the selected pathway' \
                   f' databases.'

    if not queryset:
        return f'The pathway ID {pathway_id} does not have any equivalent pathways in the selected pathway databases.'

    for obj in queryset:

        # Add current pathway to set of equivalent pathways if pathway has mappings
        if pathway_id == obj.pathway_id:
            equivalent_pathways.add(pathway_id)

        id_name_mappings[obj.pathway_id] = obj.pathway_name
        id_name_mappings[obj.decopath_id] = obj.decopath_name

        # Get all pathways equivalent to the queried pathway
        mapping_queryset = obj.mapping_pathway.all()

        for mapping_obj in mapping_queryset:

            # Skip pathways with mappings not in user selected databases
            if mapping_obj.pathway_database not in databases:
                continue

            # Get pathway IDs of equivalent pathways
            equivalent_pathways.add(mapping_obj.pathway_id)

            id_name_mappings[mapping_obj.pathway_id] = mapping_obj.pathway_name

    # Get genesets of equivalent pathways
    for pathway_id, gene_set in geneset_dict.items():

        if pathway_id in equivalent_pathways and not pathway_id.startswith('DC'):
            equivalent_genesets[pathway_id] = gene_set

    # Get geneset information for each equivalent pathway
    for pathway_id, gene_set in equivalent_genesets.items():
        database = get_database_by_pathway_id(pathway_id)

        overlaps_venn_diagram.append(
            {
                'sets': [index],
                'size': len(gene_set),
                'pathway_id': pathway_id,
                'label': f'{id_name_mappings[pathway_id]} ({database})',
                'gene_set': list(gene_set),
            }
        )

        pathway_to_index[pathway_id] = index

        index += 1

    # Get geneset overlap/intersection information
    for (set_1_name, set_1_values), (set_2_name, set_2_values) in combinations(equivalent_genesets.items(), r=2):
        genesets_intersection = set(set_1_values).intersection(set_2_values)

        overlaps_venn_diagram.append(
            {
                'sets': [pathway_to_index[set_1_name], pathway_to_index[set_2_name]],
                'size': len(genesets_intersection),
                'gene_set': list(genesets_intersection),
                'intersection': set_1_name + ' &#8745 ' + set_2_name
            }
        )

    return overlaps_venn_diagram


def generate_consensus_table_gsea(results_df: pd.DataFrame, significance_value):
    """Generate results table showing pathway database consensus."""
    df = _get_results_table(results_df)

    if 'pval' in df.columns:
        df.rename(columns={'pval': 'p_value'}, inplace=True)

    if 'fdr' in df.columns:
        df.rename(columns={'fdr': 'q_value'}, inplace=True)

    df = df[['Identifier', 'Pathway', 'Database', 'nes', 'es', 'p_value', 'q_value']]

    # Sort dataFrame by database column
    sorted_df = df.sort_values('Database')

    header = ['Pathway', 'Consensus']
    body = []
    metadata_ids = []
    metadata_qvals = []
    metadata_cons = []
    metadata_dc_cons = []
    cons_dict = defaultdict(list)

    header, databases, db_order_dict = process_database_names(sorted_df, header)

    # List of pathways already processed
    skip_ids = []

    pathway_data = []

    # For each pathway, check if it has an equivalent pathway in one of the user-selected databases
    for identifier, pathway, db, nes, es, pval, fdr in sorted_df.values:

        # Continue if equivalent pathway has already been processed
        if identifier in skip_ids:
            continue

        # Continue if pathway is a DC pathway
        if identifier.startswith('DC'):
            continue

        row = ['NA' for i in range(len(databases) + 2)]
        identifiers = ['NA' for i in range(len(databases))]
        qvals = ['NA' for i in range(len(databases))]

        scores = []
        fdr_list = []

        # Query pathway mappings model to filter pathways and retain only those with mappings
        queryset = Pathway.objects.filter(pathway_id=identifier)

        for obj in queryset:

            row, identifiers, qvals, dc_nes, dc_fdr = get_dc_info_gsea(
                obj,
                sorted_df,
                databases,
                db_order_dict,
                row,
                identifiers,
                qvals
            )

            # Get current pathway score, ID and FDR
            row[db_order_dict[db] + 2] = round(nes, 2)
            identifiers[db_order_dict[db]] = identifier
            qvals[db_order_dict[db]] = fdr

            # Get pathway score of pathway with mappings
            scores.append(nes)

            # Get q-value of pathway with mappings
            fdr_list.append(fdr)

            # Get all pathways equivalent to the queried pathway
            mapping_queryset = obj.mapping_pathway.all().order_by('pathway_database')

            for mapping_obj in mapping_queryset:

                skip_ids.append(mapping_obj.pathway_id)

                # Skip pathways with mappings not from databases in the results dataFrame
                if mapping_obj.pathway_database not in databases:
                    continue

                # Get info for the equivalent pathway from the results dataFrame
                for id_equiv, name_equiv, db_equiv, nes_equiv, es_equiv, pval_equiv, fdr_equiv in sorted_df.values:

                    # Skip already processed DecoPath pathways
                    if obj.decopath_id == id_equiv:
                        continue

                    # Get scores, ID and FDR for pathways with mappings
                    if mapping_obj.pathway_id == id_equiv:
                        row[db_order_dict[db_equiv] + 2] = round(nes_equiv, 2)
                        identifiers[db_order_dict[db_equiv]] = id_equiv
                        qvals[db_order_dict[db_equiv]] = fdr_equiv

                        # Get pathway score of pathway with mappings
                        scores.append(nes_equiv)

                        # Get FDR of pathway with mappings
                        fdr_list.append(fdr_equiv)

            consensus, dc_db_consensus = get_consensus_gsea(row, scores, fdr_list, significance_value, dc_fdr, dc_nes)

            cons_dict[CONSENSUS_MAPPINGS[consensus]].append(row[0])

            # Ensure minimum of 2 database columns exist for comparison
            if len(row) < 4:
                continue

            # Skip row if there are no scores for at least 2 databases
            if len(row[2:]) - row.count('NA') < 2:
                continue

            body.append(row[:])
            metadata_ids.append(identifiers[:])
            metadata_qvals.append(qvals)
            metadata_cons.append(consensus)
            metadata_dc_cons.append(dc_db_consensus)

            pathway_data.append((row[:2] + list(sum(list(zip(identifiers, row[2:], qvals)), ()))))

    row_db = []

    if DECOPATH in databases:
        databases.remove(DECOPATH)
        for db in databases:
            row_db.append((f'{db}_id', f'{db}_nes', f'{db}_qval'))
        row_db.append(('decopath_id', 'decopath_nes', 'decopath_qval'))
    else:
        for db in databases:
            row_db.append((f'{db}_id', f'{db}_nes', f'{db}_qval'))

    full_consensus_df = pd.DataFrame.from_records(pathway_data)
    full_consensus_df.columns = ['pathway', 'consensus'] + list(sum(row_db, ()))

    df_list = [i for i in range(len(full_consensus_df.columns))]

    df_to_html = full_consensus_df.to_html(
        classes='table',
        table_id='full_consensus_gsea_df',
        index=False,
        render_links=True,
    )

    return header, body, metadata_ids, metadata_qvals, metadata_cons, metadata_dc_cons, df_to_html, df_list, cons_dict


# TODO: test with custom mappings
def generate_consensus_table_ora(results_df: pd.DataFrame, significance_value):
    """Generate results table showing pathway database consensus."""
    df = _get_results_table(results_df)

    df = df[['Identifier', 'Pathway', 'Database', 'p_value', 'q_value']]

    # Sort dataFrame by database column
    sorted_df = df.sort_values('Database')

    header = ['Pathway', 'Consensus']
    body = []
    metadata_ids = []
    metadata_consensus = []
    metadata_dc_consensus = []
    consensus_dict = defaultdict(list)

    header, databases, db_order_dict = process_database_names(sorted_df, header)

    # List of pathways already processed
    skip_ids = []

    pathway_data = []

    # For each pathway, check if it has an equivalent pathway in one of the user-selected databases
    for identifier, pathway, db, pval, fdr in sorted_df.values:

        # Continue if equivalent pathway has already been processed
        if identifier in skip_ids:
            continue

        # Continue if pathway is a DC pathway
        if identifier.startswith('DC'):
            continue

        fdr_list = []

        row = ['NA' for i in range(len(databases) + 2)]
        identifiers = ['NA' for i in range(len(databases))]

        # Query pathway mappings model to filter pathways and retain only those with mappings
        queryset = Pathway.objects.filter(pathway_id=identifier)

        for obj in queryset:

            # Get DC pathway info if exists in user ORA results
            row, identifiers, dc_fdr = get_dc_info_ora(obj, sorted_df, databases, db_order_dict, row, identifiers)

            # Get current pathway FDR and ID
            row[db_order_dict[db] + 2] = fdr
            identifiers[db_order_dict[db]] = identifier

            # Get q-value of pathway with mappings
            fdr_list.append(fdr)

            # Get all pathways equivalent to the queried pathway
            mapping_queryset = obj.mapping_pathway.all().order_by('pathway_database')

            for mapping_obj in mapping_queryset:
                skip_ids.append(mapping_obj.pathway_id)

                # Skip pathways with mappings not in databases in the results dataFrame
                if mapping_obj.pathway_database not in databases:
                    continue

                # Get info for the equivalent pathway from the results dataFrame
                for identifier_equiv, name_equiv, db_equiv, pval_equiv, fdr_equiv in sorted_df.values:

                    # Skip already processed DecoPath pathways
                    if obj.decopath_id == identifier_equiv:
                        continue

                    if mapping_obj.pathway_id == identifier_equiv:
                        # Get scores and FDR for pathways with mappings
                        row[db_order_dict[db_equiv] + 2] = fdr_equiv
                        identifiers[db_order_dict[db_equiv]] = identifier_equiv

                        # Get q-value of pathway with mappings
                        fdr_list.append(fdr_equiv)

            # Check if mappings are in consensus
            sig_consensus = check_consensus(fdr_list, threshold=significance_value, max_num_of_mappings=2)

            row[1] = CONSENSUS_MAPPINGS[sig_consensus]

            consensus_dict[row[1]].append(row[0])

            # Check if mappings are in consensus with DC
            sig_consensus_dc = check_consensus(fdr_list, threshold=significance_value, max_num_of_mappings=1)

            if dc_fdr is None or sig_consensus_dc == 1:
                dc_db_consensus = 1  # no mappings
            elif dc_fdr < significance_value and all(fdr_val < significance_value for fdr_val in fdr_list):
                dc_db_consensus = 2  # consensus green
            elif dc_fdr > significance_value and all(fdr_val > significance_value for fdr_val in fdr_list):
                dc_db_consensus = 2  # consensus green
            else:
                dc_db_consensus = 0  # disagreement red

            rounded_row = []

            for elem in row:
                rounded_row.append(round_float(elem))

            # Ensure minimum of 2 database columns exist for comparison
            if len(row) < 4:
                continue

            # Skip row if there are no scores for at least 2 databases
            if len(row[2:]) - row.count('NA') < 2:
                continue

            body.append(rounded_row[:])
            metadata_ids.append(identifiers[:])
            metadata_consensus.append(sig_consensus)
            metadata_dc_consensus.append(dc_db_consensus)

            # Get pathway row with pathway name, consensus info and fdr vals
            pathway_data.append(
                (row[:2] + list(sum(list(zip(identifiers, row[2:])), ())))
            )

    row_db = []

    if DECOPATH in databases:
        databases.remove(DECOPATH)
        for db in databases:
            row_db.append((f'{db}_id', f'{db}_qval'))
        row_db.append(('decopath_id', 'decopath_qval'))
    else:
        for db in databases:
            row_db.append((f'{db}_id', f'{db}_qval'))

    full_consensus_df = pd.DataFrame.from_records(pathway_data)
    full_consensus_df.columns = ['pathway', 'consensus'] + list(sum(row_db, ()))

    df_list = [i for i in range(len(full_consensus_df.columns))]

    df_to_html = full_consensus_df.to_html(
        classes='table',
        table_id='full_ora_consensus_df',
        justify='left',
        index=False,
        render_links=True,
        border=0,
    )

    return (
        header,
        body,
        metadata_ids,
        metadata_consensus,
        metadata_dc_consensus,
        df_to_html,
        df_list,
        consensus_dict,
    )


def query_results_model(result_id, current_user):
    """Query results model."""
    try:
        result_object = EnrichmentResult.objects.get(result_id=result_id, user=current_user)
    except ObjectDoesNotExist:
        return 'Your experiment was not found.'

    if not result_object:
        return f'It appears {result_object} are empty. Please ensure the correct files were submitted.'

    # Get results data
    df = result_object.get_df()
    databases = result_object.get_databases()
    significance_value = result_object.significance_threshold
    enrichment_method = result_object.enrichment_method
    data_filename = result_object.data_filename
    symbol_to_fold_change = result_object.fold_change_results
    fc_filename = result_object.fold_changes_filename

    return df, databases, significance_value, enrichment_method, data_filename, symbol_to_fold_change, fc_filename
