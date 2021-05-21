# -*- coding: utf-8 -*-

"""Form processing utils module."""

import pickle
from typing import List, Tuple, Union

import pandas as pd

from viewer.models import Pathway, PathwayDatabase
from viewer.src.constants import *
from viewer.src.data_preprocessing import (
    check_text_file,
    process_data_file,
    check_label_compliance, _read_text_file, _check_df_validity, _check_mapping_df
)
from viewer.src.response_handler import MAPPING_SELECTION_MSG, RUN_DGE_ANALYSIS_MSG, DGE_OPTIONS_MSG
from viewer.src.utils import concatenate_files, get_missing_columns


def get_genesets(databases: List, method: str) -> Union[str, Tuple[bool, List, str]]:
    """Process forms submitted by user and generate gene sets."""
    # Check if mappings are valid
    mapping_is_valid = True

    # Get list of all pathway databases in mapping table
    pathway_mapping_databases = Pathway.objects.all().values_list('pathway_database').distinct()

    # Flatten list of tuples to list of pathway databases in mapping table
    mapping_pathways = list(sum(pathway_mapping_databases, ()))
    gmt_files_path = set()

    # Ensure databases selected exist in mapping table
    for db in databases:
        if db not in mapping_pathways:
            return MAPPING_SELECTION_MSG

    databases.append(DECOPATH)

    # Get queryset of PathwayDatabase objects filtered by selected databases
    pathway_database_queryset = PathwayDatabase.objects.filter(database_name__in=databases)

    # Concatenate gene sets from selected databases
    for obj in pathway_database_queryset:
        gmt_file_obj = obj.gmt_file
        gmt_file = pickle.loads(gmt_file_obj)
        gmt_files_path.add(gmt_file)

    gmt_files_list = list(gmt_files_path)

    if method == ORA:
        # Concatenate gmt files selected by user
        gmt_files = concatenate_files(gmt_files_list, databases)

    else:
        gmt_files = [
            file.replace(REACTOME_GMT, REACTOME_MAPPINGS_FILE) if file.endswith(REACTOME_GMT) else file
            for file in gmt_files_list
        ]

    return mapping_is_valid, databases, gmt_files


def add_forms_to_formset(request, formset, extra: int):
    """Add additional forms to formset."""

    # Generate copy of request.POST dictionary
    formset_dictionary_copy = request.POST.copy()

    # Update formset management form with additional forms
    formset_dictionary_copy['form-TOTAL_FORMS'] = int(formset_dictionary_copy['form-TOTAL_FORMS']) + extra

    # Update formset sent to user with additional empty forms
    return formset(formset_dictionary_copy)


def check_mapping_validity(form) -> Union[bool, Union[pd.DataFrame, str]]:
    """Check if mapping forms submitted are valid."""
    clean_mapping_file = form.cleaned_data['mapping_file']

    # Get paths to temporary uploaded files
    mapping_path = clean_mapping_file.transient_file_path()

    # Check if file can be read
    check_read = _read_text_file(mapping_path, str(clean_mapping_file))

    if isinstance(check_read, str):
        return check_read

    # Check if file is empty
    check_df = _check_df_validity(str(clean_mapping_file), check_read)

    if isinstance(check_df, str):
        return check_df

    # Check if all necessary columns present in mappings dataFrame
    check_mapping = _check_mapping_df(check_df, str(clean_mapping_file))

    if isinstance(check_mapping, str):
        return check_mapping

    return check_mapping


def _check_fold_changes_form(
    read_counts,
    significance_val_analysis,
    clean_fold_changes,
    significance_val_results,
    **kwargs,
):
    """Check fold change forms validity."""
    if 'class_labels' in kwargs:
        run_dge = [read_counts, kwargs['class_labels'], significance_val_analysis]
    else:
        run_dge = [read_counts, significance_val_analysis]

    upload_fc = [clean_fold_changes, significance_val_results]

    # Ensure user optionally submits files to only run fold change analysis or submit fold changes results but not both
    if not all(i is None for i in run_dge) and not all(i is None for i in upload_fc):
        return DGE_OPTIONS_MSG

    # Check if fold changes have been uploaded and if the file is valid
    if clean_fold_changes:
        fold_changes = process_fold_changes_upload(clean_fold_changes)

        if isinstance(fold_changes, str):
            return fold_changes

        if significance_val_results is None:
            significance_val_results = PADJ

        return fold_changes, significance_val_results

    if 'class_labels' in kwargs:
        # If user submits read counts file, ensure class label file is also submitted (and vice versa) or throw error
        if read_counts and kwargs['class_labels'] is None or kwargs['class_labels'] and read_counts is None:
            return RUN_DGE_ANALYSIS_MSG

    if read_counts:

        # Get paths to temporary uploaded files
        data_path = read_counts.transient_file_path()

        # Check if expression data is valid
        expression_file_val = process_data_file(data_path, str(read_counts))

        if isinstance(expression_file_val, str):
            return expression_file_val

        read_counts_df = expression_file_val

        if 'class_labels' in kwargs:

            class_path = kwargs['class_labels'].transient_file_path()

            # Check if class labels file is valid
            class_file = check_text_file(class_path, str(kwargs['class_labels']))

            # Check if file cannot be read and throw error
            if isinstance(class_file, str):
                return class_file

            class_label_df = class_file

            # Get list of class labels in order corresponding to read counts file
            class_labels = check_label_compliance(read_counts_df, class_label_df)

            if isinstance(class_labels, str):
                return class_labels

        else:
            class_labels = ['NA']
            class_path = ""

        if significance_val_analysis is None:
            significance_val_analysis = PADJ

        return data_path, significance_val_analysis, class_labels, class_path

    return None


def process_fold_changes_upload(clean_fold_changes) -> Union[pd.DataFrame, str]:
    """Check if fold changes form submitted is valid."""
    fold_changes_path = clean_fold_changes.transient_file_path()

    # Check if fold changes file is valid and gene identifiers on index are HGNC symbols
    fc_file_val = process_data_file(fold_changes_path, str(clean_fold_changes))

    # Check if file cannot be read and throw error
    if isinstance(fc_file_val, str):
        return fc_file_val

    fold_changes_df = fc_file_val

    missing_columns = get_missing_columns(USER_FOLD_CHANGES_COLUMN_NAME, fold_changes_df)

    if missing_columns:
        return f'Your file containing fold changes is missing some information. Please ensure your file contains' \
               f' the following column(s): {", ".join(missing_columns)}. See FAQs and sample files for details.'

    return fold_changes_df
