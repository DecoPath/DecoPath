# -*- coding: utf-8 -*-

"""This module runs GSEA."""

from typing import List, Union

import gseapy
import pandas as pd


def perform_gsea(
    data: Union[str, pd.DataFrame],
    gmt: str,
    class_vector: List,
    output_dir: str,
    max_size: int,
    min_size: int,
    permutation_type: str,
    permutation_num: int,
    method: str
):
    """Run GSEA on a given dataset and geneset."""
    return gseapy.gsea(
        data=data,
        gene_sets=gmt,
        cls=class_vector,
        max_size=max_size,
        min_size=min_size,
        # set permutation_type to phenotype if samples >=15
        permutation_type=permutation_type,
        permutation_num=permutation_num,  # reduce number to speed up test
        method=method,
        outdir=output_dir,
        no_plot=True,  # Skip plotting
        processes=1,
        format='png',
    )
