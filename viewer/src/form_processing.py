# -*- coding: utf-8 -*-

"""Form processing module."""

from typing import Tuple, Union, Any

import pandas as pd

from viewer.models import EnrichmentResult
from viewer.src.constants import *
from viewer.src.data_preprocessing import (
    check_text_file,
    process_data_file,
    check_label_compliance,
    process_data_ora
)
from viewer.src.db_utils import load_results_metadata
from viewer.src.form_processing_utils import get_genesets, _check_fold_changes_form
from viewer.src.response_handler import *
from viewer.src.utils import get_database_by_id, get_missing_columns, read_data_file
from viewer.tasks import deploy_gsea, deploy_ora, deploy_deseq


def process_user_results(
    current_user,
    user_email,
    results_form,
    results_form_fc
) -> Union[Tuple[EnrichmentResult, Any], str, bool]:
    """Process user submitted results file and load enrichment result model."""
    fold_changes_df = None
    read_counts_path = None

    # Get cleaned data
    clean_results_data = results_form.cleaned_data['results_file']
    enrichment_method = results_form.cleaned_data['select_enrichment_method']
    sig_threshold_results = results_form.cleaned_data['sig_threshold_results']

    # Get optional read counts matrix
    clean_read_counts = results_form_fc.cleaned_data['read_counts_file']
    clean_class_labels = results_form_fc.cleaned_data['class_file']
    significance_val_analysis = results_form_fc.cleaned_data['significance_value_run']

    # Get optional fold changes
    clean_fold_changes = results_form_fc.cleaned_data['fold_changes_file']
    significance_val_results = results_form_fc.cleaned_data['significance_value_upload']

    # Ensure either FC are uploaded or analysis is run and in latter case ensure read counts and labels both submitted
    fc_check = _check_fold_changes_form(
        clean_read_counts,
        significance_val_analysis,
        clean_fold_changes,
        significance_val_results,
        class_labels=clean_class_labels,
    )

    # Check if fold changes or files to run DESeq2 submitted
    if fc_check is None:
        pass

    # Process fold change file upload
    elif len(fc_check) == 2:
        fold_changes_df, sig_cutoff = fc_check

    # Process files to run DGE analysis
    elif len(fc_check) == 4:
        read_counts_path, sig_cutoff, class_labels_list, class_label_path = fc_check

    # return HTTPBadRequest
    else:
        return fc_check

    # Get path to temporary uploaded results file
    results_path = clean_results_data.transient_file_path()

    # Check if results file is valid
    results_file_val = check_text_file(results_path, str(clean_results_data))

    # Check if file cannot be read and throw error
    if isinstance(results_file_val, str):
        return results_file_val

    results_df = results_file_val

    # Get databases by pathway ID prefixes
    databases = set(results_df.index.map(get_database_by_id).to_list())

    dc_databases = {KEGG, REACTOME, PATHBANK, WIKIPATHWAYS, DECOPATH}

    # Ensure at least 2 of the 5 possible pathway databases are present in user uploaded results
    if len(dc_databases.intersection(databases)) < 2:
        return RESULTS_CHECK_DATABASES_MSG

    if sig_threshold_results is None:
        sig_threshold_results = PADJ

    # Check if column headers are valid
    if enrichment_method == GSEA:

        missing_columns = get_missing_columns(USER_RESULTS_GSEA_COLUMN_NAMES, results_df)

        if missing_columns:
            return f'Your GSEA results file {str(clean_results_data)} is missing some information. Please ensure your' \
                   f' file contains the column(s): {", ".join(missing_columns)}. See FAQs and example files for more' \
                   f' details.'

    else:

        missing_columns = get_missing_columns(USER_RESULTS_ORA_COLUMN_NAMES, results_df)

        if missing_columns:
            return f'Your ORA results file {str(clean_results_data)} is missing some information. Please ensure your' \
                   f' file contains the column(s): {", ".join(missing_columns)}. See FAQs and example files for more' \
                   f' details.'

    # Load results, metadata and FCs into results model
    if isinstance(fold_changes_df, pd.DataFrame):
        job = load_results_metadata(
            results=results_df.copy(),
            current_user=current_user,
            data_filename=str(clean_results_data),
            class_filename='NA',
            sample_number=None,
            class_labels=['NA'],
            databases=list(databases),
            significance_threshold=sig_threshold_results,
            max_size=None,
            min_size=None,
            permutation_number=None,
            permutation_type='NA',
            calculation_method='NA',
            enrichment_method=enrichment_method.upper(),
            significance_threshold_fc=sig_cutoff,
            fold_change_results=fold_changes_df.copy(),
            read_counts_path=None,
            fold_changes_filename=str(clean_fold_changes),
        )

        job.result_status = 2
        job.save()

    # Load results and metadata into results model
    elif read_counts_path:
        job = load_results_metadata(
            results=results_df.copy(),
            current_user=current_user,
            data_filename=str(clean_results_data),
            class_filename=str(clean_class_labels),
            sample_number=len(class_labels_list),
            class_labels=list(set(class_labels_list)),
            databases=list(databases),
            significance_threshold=sig_threshold_results,
            max_size=None,
            min_size=None,
            permutation_number=None,
            permutation_type='NA',
            calculation_method='NA',
            enrichment_method=enrichment_method.upper(),
            significance_threshold_fc=sig_cutoff,
            fold_change_results=None,
            fold_changes_filename=str(clean_read_counts),
            read_counts_path=read_counts_path,
        )

        # Run DESeq2
        task = deploy_deseq.delay(
            read_counts_path=read_counts_path,
            design_matrix_path=class_label_path,
            design_matrix_filename=str(clean_class_labels),
            read_counts_filename=str(clean_read_counts),
            user_mail=user_email,
            job_id=job.get_job_id()
        )

        return job, task

    # Load results and metadata directly into Enrichment Results model
    else:
        job = load_results_metadata(
            results=results_df.copy(),
            current_user=current_user,
            data_filename=str(clean_results_data),
            class_filename='NA',
            sample_number=None,
            class_labels=['NA'],
            databases=list(databases),
            significance_threshold=sig_threshold_results,
            max_size=None,
            min_size=None,
            permutation_number=None,
            permutation_type='NA',
            calculation_method='NA',
            enrichment_method=enrichment_method.upper(),
            significance_threshold_fc=None,
            fold_change_results=None,
            read_counts_path=None,
            fold_changes_filename='NA',
        )

        job.result_status = 2
        job.save()


def process_files_run_ora(
    current_user,
    user_email,
    form,
    db_form,
    parameters_form,
) -> Union[Tuple[EnrichmentResult, Any], str, bool]:
    """Process user submitted files to run ORA and load enrichment results model."""
    gene_list_path = None
    fold_changes_path = None
    read_counts_path = None

    # Get cleaned database selection
    select_database = db_form.cleaned_data['select_databases']

    # Get cleaned optional parameters
    min_size = parameters_form.cleaned_data['minimum_size']
    max_size = parameters_form.cleaned_data['maximum_size']
    sig_threshold_ora = parameters_form.cleaned_data['sig_threshold_ora']

    # Get cleaned optional gene list form
    clean_gene_list = form.cleaned_data['gene_list']

    # Get cleaned optional read counts matrix form
    clean_read_counts = form.cleaned_data['read_counts_file_ora']
    clean_class_labels = form.cleaned_data['class_file_ora']
    significance_val_analysis = form.cleaned_data['significance_value_run_ora']

    # Get cleaned optional fold changes form
    clean_fold_changes = form.cleaned_data['fold_changes_file_ora']
    significance_val_results = form.cleaned_data['significance_value_upload_ora']

    # return HTTPBadRequest if no files submitted
    if not (clean_gene_list or clean_read_counts or clean_fold_changes):
        return ORA_EMPTY_FORM_MSG

    # return HTTPBadRequest if more than 1 file submitted
    elif [clean_gene_list, clean_read_counts, clean_fold_changes].count(None) < 2:
        return ORA_UPLOAD_MSG

    # Ensure either FC are uploaded or analysis is run and in latter case ensure read counts and labels both submitted
    fc_check = _check_fold_changes_form(
        clean_read_counts,
        significance_val_analysis,
        clean_fold_changes,
        significance_val_results,
        class_labels=clean_class_labels,
    )

    if fc_check is None:
        pass

    # Process fold change file upload
    elif len(fc_check) == 2:
        fold_changes_df, sig_cutoff = fc_check
        fold_changes_path = clean_fold_changes.transient_file_path()

    # Process files to run DGE analysis
    elif len(fc_check) == 4:
        read_counts_path, sig_cutoff, class_labels_list, class_label_path = fc_check

    # return HTTPBadRequest
    else:
        return fc_check

    # Get path to temporary uploaded gene list file
    if clean_gene_list:
        gene_list_path = clean_gene_list.transient_file_path()

        # Check if gene list file is valid
        ora_file_val = process_data_ora(gene_list_path, str(clean_gene_list))

        # return HTTPBadRequest
        if isinstance(ora_file_val, str):
            return ora_file_val

        set_gene_symbols = set(ora_file_val.index.values)

    # Ensure at least two databases are selected
    if len(select_database) < 2:
        return False

    databases = [db.lower() for db in select_database]

    # Concatenate genesets in GMT format and process forms
    selected_dbs_val = get_genesets(
        databases, method=ORA
    )

    if isinstance(selected_dbs_val, str):
        return selected_dbs_val

    mapping_is_valid, database_list, gmt_file_path = selected_dbs_val

    if mapping_is_valid is False:
        return MAPPING_SELECTION_MSG

    if max_size is None:
        max_size = DEFAULT_MAX_ORA
    if min_size is None:
        min_size = DEFAULT_MIN_ORA

    if sig_threshold_ora is None:
        sig_threshold_ora = PADJ

    # Run ORA on user uploaded gene list
    if gene_list_path:
        # Load enrichment results model
        job = load_results_metadata(
            results=None,
            current_user=current_user,
            data_filename=str(clean_gene_list),
            class_filename="NA",
            class_labels=['NA'],
            sample_number=None,
            databases=database_list,
            min_size=min_size,
            max_size=max_size,
            significance_threshold=sig_threshold_ora,
            permutation_type="NA",
            permutation_number=None,
            calculation_method="NA",
            enrichment_method=ORA.upper(),
            significance_threshold_fc=None,
            fold_changes_filename="NA",
            fold_change_results=None,
            read_counts_path=None,
        )

        task = deploy_ora.delay(
            gmt_file_path=gmt_file_path,
            min_size=min_size,
            max_size=max_size,
            user_mail=user_email,
            job_id=job.get_job_id(),
            set_gene_symbols=list(set_gene_symbols),
        )

    # Run ORA on DEGs from fold change results file
    elif fold_changes_path:
        # Load enrichment results model with metadata and fold changes
        job = load_results_metadata(
            results=None,
            current_user=current_user,
            data_filename=str(clean_fold_changes),
            class_filename="NA",
            class_labels=['NA'],
            sample_number=None,
            databases=database_list,
            min_size=min_size,
            max_size=max_size,
            significance_threshold=sig_threshold_ora,
            permutation_type="NA",
            permutation_number=None,
            calculation_method="NA",
            enrichment_method=ORA.upper(),
            significance_threshold_fc=sig_cutoff,
            fold_change_results=fold_changes_df.copy(),
            fold_changes_filename=str(clean_fold_changes),
            read_counts_path=None,
        )

        task = deploy_ora.delay(
            gmt_file_path=gmt_file_path,
            min_size=min_size,
            max_size=max_size,
            user_mail=user_email,
            job_id=job.get_job_id(),
            sig_threshold_fc=sig_cutoff,
            fold_changes_path=fold_changes_path,
        )

    # Run DESeq2 and ORA on significant DEGs
    else:
        # Load enrichment results model
        job = load_results_metadata(
            results=None,
            current_user=current_user,
            data_filename=str(clean_read_counts),
            class_filename=str(clean_class_labels),
            class_labels=list(set(class_labels_list)),
            sample_number=len(class_labels_list),
            databases=database_list,
            min_size=min_size,
            max_size=max_size,
            significance_threshold=sig_threshold_ora,
            permutation_type="NA",
            permutation_number=None,
            calculation_method="NA",
            enrichment_method=ORA.upper(),
            significance_threshold_fc=sig_cutoff,
            fold_change_results=None,
            fold_changes_filename=str(clean_read_counts),
            read_counts_path=read_counts_path,
        )

        task = deploy_ora.delay(
            gmt_file_path=gmt_file_path,
            min_size=min_size,
            max_size=max_size,
            user_mail=user_email,
            job_id=job.get_job_id(),
            sig_threshold_fc=sig_cutoff,
            read_counts_path=read_counts_path,
            design_matrix_path=class_label_path,
        )

    return job, task


def process_files_run_gsea(
    current_user,
    user_email,
    form,
    db_form,
    parameters_form,
    fc_form
) -> Union[Tuple[EnrichmentResult, Any], str, bool]:
    """Process user submitted files to run enrichment and load GSEA model."""
    read_counts_path = None
    fold_changes_path = None

    # Get cleaned data
    select_database = db_form.cleaned_data['select_databases']
    clean_exp_data = form.cleaned_data['data_file']
    clean_cls = form.cleaned_data['class_file']

    # Get optional parameters
    permutation_num = parameters_form.cleaned_data['permutation_num']
    permutation_type = parameters_form.cleaned_data['permutation_type']
    method = parameters_form.cleaned_data['method']
    max_size = parameters_form.cleaned_data['maximum_size']
    min_size = parameters_form.cleaned_data['minimum_size']
    sig_threshold_gsea = parameters_form.cleaned_data['sig_threshold_gsea']

    # Get optional read counts matrix
    clean_read_counts = fc_form.cleaned_data['read_counts_file_gsea']
    significance_val_analysis = fc_form.cleaned_data['significance_value_run_gsea']

    # Get optional fold changes
    clean_fold_changes = fc_form.cleaned_data['fold_changes_file_gsea']
    significance_val_results = fc_form.cleaned_data['significance_value_upload_gsea']

    # Ensure either FC are uploaded or analysis is run and in latter case ensure read counts and labels both submitted
    fc_check = _check_fold_changes_form(
        clean_read_counts,
        significance_val_analysis,
        clean_fold_changes,
        significance_val_results,
    )

    # Check if fold changes submission is valid
    if fc_check is None:
        sig_cutoff = None
    elif len(fc_check) == 2:
        fold_changes_path, sig_cutoff = fc_check
    elif len(fc_check) == 4:
        read_counts_path, sig_cutoff, class_labels_list, _ = fc_check
        read_counts_df = read_data_file(read_counts_path, str(clean_read_counts))

    else:
        return fc_check

    # Get paths to temporary uploaded files
    data_path = clean_exp_data.transient_file_path()
    class_path = clean_cls.transient_file_path()

    # Ensure at least two databases are selected
    if len(select_database) < 2:
        return False

    # Check if expression data is valid
    expression_file_val = process_data_file(data_path, str(clean_exp_data))

    if isinstance(expression_file_val, str):
        return expression_file_val

    expression_df = expression_file_val

    databases = [db.lower() for db in select_database]

    # Concatenate genesets in GMT format and process forms
    selected_dbs_val = get_genesets(
        databases, method=GSEA
    )

    if isinstance(selected_dbs_val, str):
        return selected_dbs_val

    mapping_is_valid, database_list, gmt_file_path = selected_dbs_val

    if mapping_is_valid is False:
        return MAPPING_SELECTION_MSG

    # Check if class labels file is valid
    class_file_val = check_text_file(class_path, str(clean_cls))

    # Check if file cannot be read and throw error
    if isinstance(class_file_val, str):
        return class_file_val

    # Get list of class labels in order corresponding to expression data file
    class_labels_gsea = check_label_compliance(expression_df, class_file_val)

    # Check if file cannot be read and throw error
    if isinstance(class_labels_gsea, str):
        return class_labels_gsea

    if clean_read_counts:
        # Get list of class labels in order corresponding to read counts file
        class_labels_fc = check_label_compliance(read_counts_df, class_file_val)

        if isinstance(class_labels_fc, str):
            return class_labels_fc

    calculation_method = METHOD_MAPPING[method]

    if max_size is None:
        max_size = DEFAULT_MAX_GSEA
    if min_size is None:
        min_size = DEFAULT_MIN_GSEA

    if sig_threshold_gsea is None:
        sig_threshold_gsea = PADJ

    sample_number = len(class_labels_gsea)
    unique_labels = list(set(class_labels_gsea))

    if isinstance(fold_changes_path, pd.DataFrame):
        # Load metadata to Enrichment results model
        job = load_results_metadata(
            results=None,
            current_user=current_user,
            data_filename=str(clean_exp_data),
            class_filename=str(clean_cls),
            sample_number=sample_number,
            class_labels=unique_labels,
            databases=database_list,
            max_size=max_size,
            min_size=min_size,
            significance_threshold=sig_threshold_gsea,
            permutation_number=permutation_num,
            permutation_type=permutation_type,
            calculation_method=calculation_method,
            enrichment_method=GSEA.upper(),
            significance_threshold_fc=sig_cutoff,
            fold_changes_filename=str(clean_fold_changes),
            fold_change_results=fold_changes_path.copy(),
            read_counts_path=None,
        )
        # Run GSEA and load results to enrichment results model
        task = deploy_gsea.delay(
            data_path=data_path,
            class_labels_path=class_path,
            data_filename=str(clean_exp_data),
            class_filename=str(clean_cls),
            gmt_files=gmt_file_path,
            output_dir=GSEA_RESULTS,
            max_size=max_size,
            min_size=min_size,
            permutation_type=permutation_type,
            permutation_num=permutation_num,
            method=method,
            user_mail=user_email,
            job_id=job.get_job_id(),
        )

    else:
        # Load metadata to Enrichment results model
        job = load_results_metadata(
            results=None,
            current_user=current_user,
            data_filename=str(clean_exp_data),
            class_filename=str(clean_cls),
            sample_number=sample_number,
            class_labels=unique_labels,
            databases=database_list,
            max_size=max_size,
            min_size=min_size,
            significance_threshold=sig_threshold_gsea,
            permutation_number=permutation_num,
            permutation_type=permutation_type,
            calculation_method=calculation_method,
            enrichment_method=GSEA.upper(),
            significance_threshold_fc=sig_cutoff,
            fold_change_results=None,
            read_counts_path=read_counts_path,
            fold_changes_filename=str(clean_read_counts),
        )

        # Run GSEA and load results to enrichment results model
        task = deploy_gsea.delay(
            data_path=data_path,
            class_labels_path=class_path,
            data_filename=str(clean_exp_data),
            class_filename=str(clean_cls),
            gmt_files=gmt_file_path,
            output_dir=GSEA_RESULTS,
            max_size=max_size,
            min_size=min_size,
            permutation_type=permutation_type,
            permutation_num=permutation_num,
            method=method,
            user_mail=user_email,
            job_id=job.get_job_id(),
            read_counts_path=read_counts_path,
            read_counts_filename=str(clean_read_counts),
        )

    return job, task
