# -*- coding: utf-8 -*-

"""Data processing module."""
import logging
from typing import Dict, Union, List, Optional

import pandas as pd

from viewer.src.constants import *
from viewer.src.response_handler import *
from viewer.src.utils import _get_hgnc_mapping_dict, get_missing_columns

logger = logging.getLogger(__name__)

"""Check if data files submitted by user are valid."""


def _read_text_file(file_path: str, filename, index_col: int = 0) -> Union[pd.DataFrame, str]:
    """Check read data file."""
    try:
        if file_path.endswith(CSV):
            df = pd.read_csv(file_path, sep=",", index_col=index_col)

        elif file_path.endswith(TSV):
            df = pd.read_csv(file_path, sep="\t", index_col=index_col)

        else:
            df = pd.read_csv(file_path, index_col=index_col, sep=None, engine='python')

    except pd.errors.ParserError:
        return f'There is a problem with your {filename} file. please ensure it contains the correct number of columns.'

    except IOError:
        return f'There is a problem with your {filename} file. please check that it meets the criteria.'

    return df


def _check_df_validity(filename: str, df: pd.DataFrame) -> Union[pd.DataFrame, str]:
    """Check if file is valid."""
    # Check dimensions of dataFrame and ensure it is not empty
    if df.empty is True or df.shape[1] == 0:
        return f'Your file {filename} appears to be empty. Please ensure it meets the criteria.'

    return df


"""Check if user submitted file is valid."""


def check_text_file(file_path: str, filename) -> Union[pd.DataFrame, str]:
    """Check user submitted results file."""
    check_read_file = _read_text_file(file_path, filename)

    # Check if file cannot be read and throw error
    if isinstance(check_read_file, str):
        return check_read_file

    # Check if file is empty
    result = _check_df_validity(filename, check_read_file)

    # Throw error if file is empty
    if isinstance(result, str):
        return result

    # Return results dataFrame
    return result


"""Check if user submitted class file corresponds to results file."""


def check_label_compliance(expr_df: pd.DataFrame, class_label_df: pd.DataFrame) -> Union[List, str]:
    """Ensure sample IDs in class labels dataFrame correspond to IDs in expression dataFrame."""
    class_labels = []

    if 'gene_symbol' in expr_df:
        expr_df.set_index('gene_symbol', inplace=True)

    # Ensure sample labels are ordered identically in expression and class label files
    expr_ids = expr_df.columns.to_list()
    label_ids = class_label_df.index.to_list()

    id_label_dict = dict(zip(class_label_df.index, class_label_df.class_label))

    for sample_id in expr_df:
        try:
            class_labels.append(id_label_dict[sample_id])
        except KeyError as e:
            return 'The following sample identifier could not be found in the class labels file: ' + str(e)

    # if expr_ids != label_ids:
    #     return CLS_COMPLIANCE_MSG

    return class_labels


"""Check if expression data submitted by user is valid."""


def check_hgnc_symbols(df: pd.DataFrame, hgnc_symbol_threshold=0.4) -> Union[pd.DataFrame, str]:
    """Check if expression dataFrame contains threshold number of HGNC symbols."""
    hgnc_symbols = _get_hgnc_mapping_dict()

    if 'gene_symbol' in df.columns:
        df.reset_index(inplace=True)
        df.set_index('gene_symbol', inplace=True)

    blacklist = [
        symbol
        for symbol in df.index.values if symbol not in hgnc_symbols
    ]

    if len(blacklist) / len(df.index) > hgnc_symbol_threshold:
        return HGNC_SYMBOL_CHECK_MSG

    # Remove rows with blacklisted symbols
    df.drop(index=blacklist, inplace=True)

    return df


def process_data_file(file_path: str, filename: str) -> Union[pd.DataFrame, str]:
    """Check if data file is valid."""
    # Check read data file
    check_read_res = _read_text_file(file_path, filename)

    if isinstance(check_read_res, str):
        return check_read_res

    # Check if dataFrame is empty
    check_df = _check_df_validity(filename, check_read_res)

    if isinstance(check_df, str):
        return check_df

    # Check if gene identifiers are HGNC symbols
    check_symbols = check_hgnc_symbols(check_df)

    if isinstance(check_symbols, str):
        return check_symbols

    return check_symbols


def process_data_ora(file_path: str, filename: str) -> Union[pd.DataFrame, str]:
    """Check if ORA data file is valid."""
    # Check read ORA data file
    check_read_res = _read_text_file(file_path, filename)

    if isinstance(check_read_res, str):
        return check_read_res

    if len(check_read_res.index) == 0:
        return EMPTY_DF_ERROR_MSG

    # Check if gene identifiers are HGNC symbols
    check_symbols = check_hgnc_symbols(check_read_res)

    if isinstance(check_symbols, str):
        return check_symbols

    return check_symbols


"""Check if .gmt file is valid."""


def _check_symbols_gmt(geneset_dict: Dict[str, List], hgnc_symbol_threshold=0.4) -> Union[None, str]:
    """Check if gene sets contain valid HGNC symbols."""
    hgnc_symbols = _get_hgnc_mapping_dict()

    try:

        blacklist = [
            symbol
            for geneset in geneset_dict.values()
            for symbol in geneset
            if symbol not in hgnc_symbols
        ]

        if len(blacklist) / len(geneset_dict) > hgnc_symbol_threshold:
            return HGNC_SYMBOL_GMT_MSG

    except ValueError:
        return GMT_FILE_ERROR_MSG

    return None


def parse_custom_gmt(file: str) -> Union[Dict[str, list], str]:
    """Parse custom gene sets in GMT file format."""
    # Check read GMT file
    try:
        df = pd.read_table(file, header=None, sep='\n')
        df[0].str.split('\t', expand=True)

    except IOError:
        return CUSTOM_GMT_MSG

    # Check parse gmt file
    check_genesets = parse_gmt_file(file)

    if isinstance(check_genesets, str):
        return check_genesets

    # Check if dictionary of gene sets is empty
    if not check_genesets:
        return EMPTY_GMT_MSG

    # Check if gene sets contain valid HGNC symbols
    check_symbols = _check_symbols_gmt(check_genesets)

    if isinstance(check_symbols, str):
        return HGNC_SYMBOL_GMT_MSG

    return check_genesets


def parse_gmt_file(gmt_path: str, min_size=3, max_size=3000) -> Union[Dict[str, list], str]:
    """Parse gmt file."""
    try:
        with open(gmt_path) as f:

            genesets_dict = {
                line.strip().split("\t")[0]: line.strip().split("\t")[2:]
                for line in f
            }

        return {
            k: v for k, v in genesets_dict.items() if min_size <= len(v) <= max_size}

    except IOError:
        return GMT_FILE_ERROR_MSG


"""Check if mapping file submitted by user is valid."""


# TODO: complete mapping check
def _check_mapping_df(df: pd.DataFrame, filename: str) -> Union[pd.DataFrame, str]:
    """Check if mandatory columns present in mapping file."""
    mapping_columns = {SOURCE_ID, SOURCE_NAME, SOURCE_RESOURCE, TARGET_ID, TARGET_NAME, TARGET_RESOURCE}

    for column in mapping_columns:
        if column == df.index.name:
            df.reset_index(inplace=True)

    # Check column names of dataFrame submitted by user
    missing_columns = get_missing_columns(mapping_columns, df)

    if missing_columns:
        return f'Your mapping file {filename} is missing some information. Please ensure your file contains' \
               f' the following column(s): {", ".join(missing_columns)}. See FAQs and sample files for details.'

    return df


def check_instance(result: Union[pd.DataFrame, str]) -> Optional[str]:
    """Check if instance of input is of type string."""
    if isinstance(result, str):
        response = result
        return response
