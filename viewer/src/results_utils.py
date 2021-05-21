# -*- coding: utf-8 -*-

"""Handle results module."""
from typing import List, Union, Tuple, Optional, Dict

import pandas as pd

from viewer.src.constants import *


def clean_none_values(item: Union[int, None]):
    if item is None:
        return 'NA'
    return item


def check_consensus(value_list: List, threshold: float, max_num_of_mappings: int) -> int:
    """Check if results of equivalent pathways are in consensus."""
    # No mappings
    if len(value_list) < max_num_of_mappings:
        return 1

    # All databases are in consensus
    elif all(i > threshold for i in value_list) or all(
            i < threshold for i in value_list) or all(
            i == threshold for i in value_list):
        return 2

    # Databases are in disagreement
    else:
        return 0


def round_float(elem: float):
    if isinstance(elem, float) and abs(elem) < 0.01 and elem != 0:
        return "{:.2e}".format(elem)

    elif isinstance(elem, float) and abs(elem) >= 0.01:
        return round(elem, 3)

    else:
        return elem


def process_database_names(df: pd.DataFrame, header: List) -> Tuple[List, List, Dict]:
    """Get sorted databases and clean names."""
    databases = [db for db in sorted(df.Database.unique())]

    db_order_dict = {
        database: index
        for index, database in enumerate([db for db in databases if db is not DECOPATH])
    }

    if DECOPATH in databases:
        db_order_dict[DECOPATH] = len(databases) - 1

    for db in db_order_dict:
        header.append(DATABASES[db])

    return header, databases, db_order_dict


def get_consensus_gsea(
        row: List,
        scores: List,
        fdr_list: List,
        significance_value: float,
        dc_fdr: Optional[float],
        dc_nes: Optional[float]
) -> Tuple[int, int]:
    """Check if all mappings are in consensus."""
    score_consensus = check_consensus(scores, threshold=0, max_num_of_mappings=2)
    sig_consensus = check_consensus(fdr_list, threshold=significance_value, max_num_of_mappings=2)

    if sig_consensus == 1:
        row[1] = 'no-mappings'
        consensus = 1
    elif sig_consensus == 2 and score_consensus == 2:
        row[1] = 'Concordant'
        consensus = 2
    else:
        row[1] = 'Discordant'
        consensus = 0

    # Check if database mappings are in consensus with DC pathways
    dc_sig_consensus = check_consensus(fdr_list, threshold=significance_value, max_num_of_mappings=1)

    if dc_fdr is None or dc_nes is None or dc_sig_consensus == 1:
        dc_db_consensus = 1  # no mappings

    elif dc_fdr < significance_value and all(
            fdr_val < significance_value for fdr_val in fdr_list) or dc_fdr > significance_value and all(
            fdr_val > significance_value for fdr_val in fdr_list):

        if dc_nes > 0 and all(score > 0 for score in scores):
            dc_db_consensus = 2  # consensus green
        elif dc_nes < 0 and all(score < 0 for score in scores):
            dc_db_consensus = 2  # consensus green
        else:
            dc_db_consensus = 0  # disagreement red

    else:
        dc_db_consensus = 0  # disagreement red

    return consensus, dc_db_consensus


def get_dc_info_gsea(
        obj: object,
        df: pd.DataFrame,
        databases: List,
        database_order: dict,
        row: List,
        identifiers: List,
        qvals: List,
) -> Tuple[List, List, List, Optional[float], Optional[float]]:
    """Get DecoPath pathway info if mapping exists for GSEA."""
    dc_nes = None
    dc_fdr = None

    # Check if DecoPath in list of databases and move to front of table
    if DECOPATH in databases:

        # Check if pathway has a DC ID
        if obj.decopath_id:

            # Get DC pathway name
            row[0] = obj.decopath_name

            # If results dataFrame has the DC ID, get DC pathway info
            if obj.decopath_id in df.Identifier.values:

                # Get DC score
                dc_nes = df.loc[df['Identifier'] == obj.decopath_id, 'nes'].iloc[0]
                row[database_order[DECOPATH] + 2] = round(dc_nes, 2)

                # Get DC FDR
                dc_fdr = df.loc[df['Identifier'] == obj.decopath_id, 'q_value'].iloc[0]
                identifiers[database_order[DECOPATH]] = obj.decopath_id
                qvals[database_order[DECOPATH]] = dc_fdr

            else:
                row[database_order[DECOPATH] + 2] = 'NA'
                identifiers[database_order[DECOPATH]] = 'NA'
                qvals[database_order[DECOPATH]] = 'NA'

        else:
            row[0] = obj.pathway_name
            row[database_order[DECOPATH] + 2] = 'NA'
            identifiers[database_order[DECOPATH]] = 'NA'
            qvals[database_order[DECOPATH]] = 'NA'

    else:
        # Get current pathway name
        row[0] = obj.pathway_name

    return row, identifiers, qvals, dc_nes, dc_fdr


def get_dc_info_ora(
        obj: object,
        df: pd.DataFrame,
        databases: List,
        database_order: dict,
        row: List,
        identifiers: List,
) -> Tuple[List, List, Optional[float]]:
    """Get DecoPath pathway info if mapping exists for ORA."""
    dc_fdr = None

    # Check if DecoPath in list of databases and add to front of table
    if DECOPATH in databases:

        # Check if pathway has a DC ID
        if obj.decopath_id:

            # Get DC pathway name
            row[0] = obj.decopath_name

            # Get DC q-value
            if obj.decopath_id in df.Identifier.values:
                dc_fdr = df.loc[df['Identifier'] == obj.decopath_id, 'q_value'].iloc[0]
                row[database_order[DECOPATH] + 2] = dc_fdr
                identifiers[database_order[DECOPATH]] = obj.decopath_id

            else:
                row[database_order[DECOPATH] + 2] = 'NA'
                identifiers[database_order[DECOPATH]] = 'NA'

        else:
            row[0] = obj.pathway_name
            row[database_order[DECOPATH] + 2] = 'NA'
            identifiers[database_order[DECOPATH]] = 'NA'

    else:
        # Get current pathway name
        row[0] = obj.pathway_name

    return row, identifiers, dc_fdr
