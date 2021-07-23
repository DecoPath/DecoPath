# -*- coding: utf-8 -*-

"""This module runs GSEA and pre-ranked GSEA."""

from typing import List, Union

import gseapy
import pandas as pd


def perform_gsea(
    data: Union[str, pd.DataFrame],
    gmt: str,
    class_vector: List,
    output_dir: str,
    min_size: int,
    max_size: int,
    permutation_type: str,
    permutation_num: int,
    method: str
):
    """Run GSEA on a given dataset and geneset."""
    return gseapy.gsea(
        data=data,
        gene_sets=gmt,
        cls=class_vector,
        min_size=min_size,
        max_size=max_size,
        permutation_type=permutation_type,  # set permutation_type to phenotype if samples >=15
        permutation_num=permutation_num,  # reduce number to speed up test
        method=method,
        outdir=output_dir,
        no_plot=True,  # Skip plotting
        processes=1,
    )


def perform_prerank(
    rnk: pd.DataFrame,
    gmt: str,
    output_dir: str,
    min_size: int,
    max_size: int,
    permutation_num: int,
):
    """Run GSEA on a pre-ranked list of genes."""
    return gseapy.prerank(
        rnk=rnk,
        gene_sets=gmt,
        min_size=min_size,
        max_size=max_size,
        permutation_num=permutation_num,
        outdir=output_dir,
        no_plot=True,  # Skip plotting
        processes=1,
    )
