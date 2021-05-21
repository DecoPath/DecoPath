# -*- coding: utf-8 -*-

"""This module contains the functions to run Over Representation Analysis (ORA)."""

import logging
from typing import Mapping, Set, List

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

from viewer.src.constants import GENE_UNIVERSE
from viewer.src.utils import spliterate

logger = logging.getLogger(__name__)


def run_ora(gmt_path: str, set_gene_symbols: Set[str], min_size: int, max_size: int):
    """Run hyper-geometric test."""
    # Note that by default, parser filters out gene sets smaller than 3 and larger than 5000
    gene_sets = gmt_parser(gmt_path, min_size=min_size, max_size=max_size)

    df = perform_hypergeometric_test(
        genes_to_test=set_gene_symbols,
        pathway_dict=gene_sets,
    )

    logger.info(f'# of pathways enriched {len(df.index)}')

    return df


def _prepare_hypergeometric_test(
    query_gene_set: Set[str],
    pathway_gene_set: Set[str],
    gene_universe: int,
) -> np.ndarray:
    """Prepare the matrix for hypergeometric test calculations.
    :param query_gene_set: set of genes to test against pathway
    :param pathway_gene_set: pathway gene set
    :param gene_universe: total number of HGNC symbols
    :return: 2x2 matrix
    """
    # Cast lists to sets
    if not isinstance(query_gene_set, set):
        query_gene_set = set(query_gene_set)
    if not isinstance(pathway_gene_set, set):
        pathway_gene_set = set(pathway_gene_set)

    # Return matrix to test hyper-geometric test
    return np.array([
        [
            len(query_gene_set.intersection(pathway_gene_set)),
            len(query_gene_set.difference(pathway_gene_set)),
        ],
        [
            len(pathway_gene_set.difference(query_gene_set)),
            gene_universe - len(pathway_gene_set.union(query_gene_set)),
        ],
    ])


def perform_hypergeometric_test(
    genes_to_test: Set[str],
    pathway_dict: Mapping[str, List[str]],
    gene_universe: int = GENE_UNIVERSE,
) -> pd.DataFrame:
    """Perform hypergeometric tests.
    :param genes_to_test: gene set to test against pathway
    :param pathway_dict: pathway name to gene set
    :param gene_universe: number of HGNC symbols
    """
    rows = []

    for pathway_id, pathway_gene_set in pathway_dict.items():
        # Prepare the test table to conduct the fisher test
        test_table = _prepare_hypergeometric_test(genes_to_test, pathway_gene_set, gene_universe)
        # Calculate fisher test (returns tuple of odds ratio and p_value
        p_value = fisher_exact(test_table, alternative='greater')[1]
        rows.append((pathway_id, p_value))

    df = pd.DataFrame(rows, columns=['pathway_id', 'p_value'])

    correction_test = multipletests(df.p_value, method='fdr_bh')
    df['q_value'] = correction_test[1]

    return df


def gmt_parser(
    path: str,
    min_size: int,
    max_size: int,
    gene_list=None,
) -> dict:
    """Parse GMT file."""
    with open(path) as file:

        # Get dictionary with pathway and corresponding gene set
        genesets_dict = {
            name: genes
            for name, _, *genes in spliterate(file)
        }
    # Apply gene set filter
    genesets_filter = {
        key: genes
        for key, genes in genesets_dict.items()
        if min_size < len(genes) <= max_size
    }

    if gene_list is not None:
        subsets = sorted(genesets_filter.keys())
        for subset in subsets:
            tag_indicator = np.in1d(gene_list, genesets_filter.get(subset), assume_unique=True)
            tag_len = sum(tag_indicator)
            if tag_len <= min_size or tag_len >= max_size:
                del genesets_filter[subset]
            else:
                continue

    filsets_num = len(genesets_dict) - len(genesets_filter)
    logging.info(f"{filsets_num} gene sets were removed with filters: max_size={max_size} and min_size={min_size}")

    if filsets_num == len(genesets_dict):
        logging.error(
            "No gene sets passed filtering condition!! Try new parameters!\n" + "GSEApy gene names are case sensitive."
        )
        raise Exception(
            "No gene sets passed filtering condition!! Try new parameters!\n" + "GSEApy gene names are case sensitive.")
    else:
        return genesets_filter
