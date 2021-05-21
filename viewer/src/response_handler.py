# -*- coding: utf-8 -*-

"""Response headers and messages module."""

"""Data-related"""

EMPTY_DF_ERROR_MSG = 'Your file appears to be empty. Please ensure it meets the criteria.'

EXP_DATA_VALIDITY_MSG = 'There is a problem with your expression data file. please ensure that it meets the criteria.'

RESULTS_UPLOAD_MSG = 'Your file is missing some information. Please ensure your file contains columns named as follows:' \
                     ' pathway, es, nes, p_value and q_value. See example files for details.'

RESULTS_PATHWAY_IDS_MSG = 'Your file is missing some information. Please ensure the first column of your file' \
                          ' contains appropriate pathway identifiers. See the example files for details.'

RESULTS_CHECK_DATABASES_MSG = 'It appears the file you have submitted does not contain any pathway databases included' \
                              ' in DecoPath. Please ensure your results have pathways from at least 2 of the ' \
                              'following databases: KEGG (pathway prefix: hsa), PathBank (pathway prefix: PW), ' \
                              'Reactome (pathway prefix: R-HSA), and WikiPathways (pathway prefix: WP). See FAQs for ' \
                              'more details.'

CLS_COMPLIANCE_MSG = 'Something went wrong. Please ensure that the sample identifiers in the expression data file you' \
                     ' have uploaded correspond to the sample identifiers you have reported in the class labels file' \
                     ' and samples in both files are in the same order.'

HGNC_SYMBOL_CHECK_MSG = 'Please input a dataset where the majority of gene symbols are from HGNC.'

ORA_EMPTY_FORM_MSG = 'Looks like no files were uploaded to run ORA. Please ensure either a gene list is uploaded or' \
                     ' files for differential gene expression analysis.'

ORA_UPLOAD_MSG = 'To run ORA, please ensure only one of the following is uploaded: a file containing a list of genes,' \
                 ' files to run differential gene expression analysis or results of differential gene expression' \
                 ' analysis. If the problem still persists, try refreshing the home page.'

"""Database-related"""

DB_SELECTION_MSG = 'Please ensure at least two databases are selected to run DecoPath.'

MAPPING_SELECTION_MSG = 'Please ensure mappings between selected databases have been uploaded.'

"""Other"""

REGISTRATION_MSG = 'You have now successfully registered. Confirm you email address before logging in.'

CHECK_DETS_MSG = 'The details were not correct. Please check your details once again'

EXCEED_JOBS_MSG = 'Job quota exceeded. Please delete a job before starting another job.'

FORM_COMPLETION_MSG = 'Something went wrong. To run DecoPath, please ensure at least two databases are selected and' \
                      ' the data files you have submitted are correctly formatted. Please see example files for ' \
                      'correct formatting and acceptable file extensions.'

MAX_MIN_GENES_MSG = "The number of genes you have selected must be within a range from 3-10,000."

INVALID_RESULTS_FORM_MSG = 'Something went wrong. Please ensure the results you have submitted are correctly ' \
                           'formatted. Please see the FAQs and example files for correct formatting and acceptable' \
                           ' file extensions.'

INVALID_FOLD_CHANGES_MSG = 'Something went wrong. Please ensure the file(s) you have submitted to run or submit ' \
                           'results for differential gene expression analysis are correctly formatted. Please see the' \
                           ' FAQs and example files for correct formatting and acceptable file extensions.'

RUN_DGE_ANALYSIS_MSG = 'Please ensure both the read counts file and its corresponding class labels files are submitted.'

DGE_OPTIONS_MSG = 'Please ensure only 1 of the 2 options is selected. You can either submit a file containing ' \
                  'log 2 fold changes or the files required to perform differential gene expression analysis. ' \
                  'See the FAQs section for more details.'

DOWNLOAD_ONTOLOGY_MSG = "Looks like the DecoPath ontology has not been loaded. Please make sure it's downloaded first."

DOWNLOAD_DECOPATH_NAMES = "Looks like the DecoPath ID to name mapping file is missing. Please make sure it's " \
                          "downloaded first."

"""Custom database-related"""

CUSTOM_DB_MSG = 'Something went wrong. Please ensure each database you add is unique and differs from the default' \
                ' databases, KEGG, Reactome and WikiPathways. If you wish to upload different versions of the' \
                ' same database, please ensure each database has a unique name (e.g., KEGG_V2).'

GMT_FILE_ERROR_MSG = 'There is a problem with your GMT file. Please ensure a correctly formatted file was submitted' \
                     ' and all fields have been entered.'

CUSTOM_GMT_MSG = 'There is a problem with your file. Please ensure the file is correctly formatted with a .csv or ' \
                 '.tsv file extension.'

EMPTY_GMT_MSG = 'Your GMT file is empty. Please enter a valid GMT file.'

DB_PRELOADED_ERROR = 'It appears the database you have selected has already been loaded into DecoPath. Please' \
                     ' navigate to the Home page to to run your experiments.'

HGNC_SYMBOL_GMT_MSG = 'Please input a gmt file where the majority of gene symbols are from HGNC.'

CUSTOM_MAPPINGS_MSG = 'Something went wrong. Please ensure all fields are filled and the mapping file is correctly ' \
                      'formatted with a .csv or .tsv file extension.'

CUSTOM_MAPPING_DATABASES_MSG = 'It appears the database you added does not correspond to any mapping files. Please' \
                               ' ensure mappings are available across each database pair.'

