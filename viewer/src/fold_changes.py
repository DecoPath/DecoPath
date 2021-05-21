# -*- coding: utf-8 -*-

"""Run DESeq2 through rpy2 module."""
import sys

import rpy2
from rpy2.robjects import pandas2ri, Formula, r, conversion, default_converter

pandas2ri.activate()
from rpy2.robjects.packages import importr

to_dataframe = r('function(x) data.frame(x)')

deseq = importr('DESeq2')
pandas2ri.activate()


class py_DESeq2:

    def __init__(self, count_matrix, design_matrix, design_formula, gene_column='id'):
        try:
            assert gene_column in count_matrix.columns, 'Wrong gene id column name'

        except AttributeError:
            sys.exit('Wrong Pandas dataframe?')

        self.dds = None
        self.deseq_result = None
        self.comparison = None
        self.normalized_count_matrix = None
        self.gene_column = gene_column
        self.gene_id = count_matrix[self.gene_column]

        count_matrix = count_matrix.drop(gene_column, axis=1)

        print(
            f'Number of columns in counts data {count_matrix.shape[1]} | '
            f'Number of rows in design matrix {design_matrix.shape[0]}'
        )

        # Load dataframe into R environment
        # Important: Change to r.data() if you use numpys and rpy2 latests versions
        self.count_matrix = rpy2.robjects.conversion.py2rpy(count_matrix)

        # Assign columns to NULL
        self.count_matrix.names = rpy2.rinterface.NULL

        self.count_matrix = count_matrix

        self.design_matrix = rpy2.robjects.conversion.py2rpy(design_matrix)

        self.design_formula = Formula(design_formula)

    def run_deseq(self, **kwargs):
        self.dds = deseq.DESeqDataSetFromMatrix(
            countData=self.count_matrix,
            colData=self.design_matrix,
            design=self.design_formula
        )
        self.dds = deseq.DESeq(self.dds, **kwargs)
        # Previous script had "deseq.counts" instead
        # TODO check standard practice
        self.normalized_count_matrix = deseq.counts_DESeqDataSet(self.dds, normalized=True)

    def get_deseq_result(self, **kwargs):

        self.comparison = deseq.resultsNames(self.dds)

        self.deseq_result = deseq.results(self.dds, **kwargs)
        self.deseq_result = to_dataframe(self.deseq_result)
        self.deseq_result = rpy2.robjects.conversion.py2rpy(self.deseq_result)  # back to pandas dataframe

        with conversion.localconverter(default_converter + pandas2ri.converter):
            results = rpy2.robjects.conversion.rpy2py(self.deseq_result)

        results[self.gene_column] = self.gene_id.values

        return results
