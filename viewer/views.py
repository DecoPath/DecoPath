# -*- coding: utf-8 -*-

"""Views module."""

import json
import os.path
import pickle

from celery.result import AsyncResult
from django.contrib import messages
from django.contrib.auth import authenticate, login, logout, update_session_auth_hash
from django.contrib.auth.decorators import login_required
from django.contrib.auth.forms import AuthenticationForm, PasswordChangeForm
from django.core.exceptions import ObjectDoesNotExist
from django.db.models import F
from django.forms import formset_factory
from django.http import HttpResponseBadRequest
from django.shortcuts import render, redirect
from django.utils.safestring import SafeString
from pandas import DataFrame, notna

from viewer.forms import *
from viewer.models import PathwayHierarchy, User, Pathway
from viewer.src.constants import *
from viewer.src.data_preprocessing import (
    parse_custom_gmt, parse_gmt_file)
from viewer.src.db_utils import (
    load_custom_pathway_databases,
    load_custom_pathway_mappings,
)
from viewer.src.form_processing import (
    process_user_results,
    process_files_run_ora,
    process_files_run_gsea
)
from viewer.src.form_processing_utils import add_forms_to_formset, check_mapping_validity
from viewer.src.handle_results import (
    create_summary_table,
    query_results_model,
    get_ranking_table,
    generate_consensus_table_gsea,
    generate_consensus_table_ora,
    get_results_circle_viz,
    process_overlap_for_venn_diagram,
)
from viewer.src.response_handler import *
from viewer.src.utils import (map_results_to_hierarchy, _get_gmt_dict, handle_file_download, get_dc_pathway_resources,
                              del_user)

"""HTML Pages"""


def home(request):
    """Render home page."""
    context = {
        'results_form': UploadResultsForm(),
        'fold_changes_form': UploadFoldChangesForm(),
        'ora_db_form': ORADatabaseSelectionForm(),
        'ora_parameters_form': ORAParametersForm(),
        'ora_form': OraOptions(),
        'gsea_form': UploadGSEAForm(),
        'gsea_db_form': GSEADatabaseSelectionForm(),
        'gsea_parameters_form': GSEAParametersForm(),
        'fold_changes_form_gsea': UploadFoldChangesFormGSEA(),
    }
    return render(request, "viewer/home.html", context)


def login_request(request):
    """."""
    if request.method == 'POST':
        form = AuthenticationForm(request=request, data=request.POST)
        if form.is_valid():
            username = form.cleaned_data.get('username')
            password = form.cleaned_data.get('password')
            user = authenticate(username=username, password=password)
            if user is not None:
                login(request, user)
                messages.success(request, f"You are now logged in as {username}")
                return redirect('/')
            else:
                messages.error(request, "Invalid username or password.")
        else:
            messages.error(request, "Invalid username or password.")
    form = AuthenticationForm()
    return render(request=request, template_name="viewer/login.html", context={"form": form})


@login_required
def logout_request(request):
    logout(request)
    messages.success(request, "Logged out successfully!")
    return redirect("/")


def register(request):
    if request.method == 'POST':
        form = UserCreationForm(data=request.POST)
        if form.is_valid():
            form.save()
            messages.success(request, REGISTRATION_MSG)
            return redirect('/')
        else:
            messages.error(request, CHECK_DETS_MSG)
    form = UserCreationForm()
    return render(request=request, template_name='viewer/registration.html', context={'form': form})


def overview(request):
    """Render overview page."""
    return render(request, "viewer/overview.html")


def user_guide(request):
    """Render user guide page."""
    return render(request, "viewer/user_guide.html")


def faqs(request):
    """Render FAQs page."""
    return render(request, "viewer/faqs.html")


def imprint(request):
    """Render imprint page."""
    return render(request, "viewer/imprint.html")


def data_protection(request):
    """Render data protection page."""
    return render(request, "viewer/data_protection.html")


@login_required
def account(request):
    """Render accounts' page for user to manage their account."""
    details = {
        'Email Address': request.user.email,
        'Password': request.user.password,
        'Number of currently running jobs': request.user.get_num_of_jobs(),
    }

    if "pass_change" in request.POST:
        form = PasswordChangeForm(request.user, request.POST)
        if form.is_valid():
            user = form.save()
            update_session_auth_hash(request, user)  # Important!
            messages.success(request, 'Your password was successfully updated!')
            return redirect('/account')
        else:
            messages.error(request, 'Please correct the error below.')

    elif "del_acc" in request.POST:
        form = AuthenticationForm(request=request, data=request.POST)
        if form.is_valid():
            username = form.cleaned_data.get('username')
            password = form.cleaned_data.get('password')
            user = authenticate(username=username, password=password)
            current_user = request.user
            if user is not None:
                if user == current_user:
                    user_email = user.email
                    logout(request)
                    res = del_user(user_email)
                    if res == 1:
                        messages.success(request, 'Your account has been successfully deleted!')
                        return redirect('/')
                else:
                    messages.error(request, "Invalid user.")
            else:
                messages.error(request, "Invalid username or password.")
        else:
            messages.error(request, "Invalid username or password.")

    pass_form = PasswordChangeForm(request.user)
    del_acc_form = AuthenticationForm()
    return render(
        request=request,
        template_name="viewer/account.html",
        context={
            'details': details,
            'pass_form': pass_form,
            'del_acc_form': del_acc_form,
        }
    )


@login_required
def run_decopath(request):
    """Process file submission page."""
    # Get current user
    current_user = request.user
    user_email = current_user.email

    # Check if user submits pathway analysis results
    if 'submit_results' in request.POST:

        # Populate form with data from the request
        results_form = UploadResultsForm(request.POST, request.FILES)
        results_form_fc = UploadFoldChangesForm(request.POST, request.FILES)

        # If form is not valid, return HttpResponseBadRequest
        if not results_form.is_valid():
            return HttpResponseBadRequest(INVALID_RESULTS_FORM_MSG)

        # If form is not valid, return HttpResponseBadRequest
        if not results_form_fc.is_valid():
            return HttpResponseBadRequest(INVALID_FOLD_CHANGES_MSG)

        # If form is valid, process file data in request.FILES and throw HTTPBadRequest if improper submission
        results_file_val = process_user_results(current_user, user_email, results_form, results_form_fc)

        # Check if file cannot be read and throw error
        if isinstance(results_file_val, str):
            return HttpResponseBadRequest(results_file_val)

    # Check if user submits form to run ORA
    elif 'run_ora' in request.POST:

        db_form = ORADatabaseSelectionForm(request.POST, request.FILES)
        parameters_form = ORAParametersForm(request.POST, request.FILES)
        form = OraOptions(request.POST, request.FILES)  # dict

        # Check if form is valid
        if not db_form.is_valid():
            return HttpResponseBadRequest(DB_SELECTION_MSG)

        # Check if parameters form is valid
        if not parameters_form.is_valid():
            return HttpResponseBadRequest(MAX_MIN_GENES_MSG)

        # Check if form is valid
        if not form.is_valid():
            return HttpResponseBadRequest(FORM_COMPLETION_MSG)

        # If form is valid, process file data in request.FILES and throw HTTPBadRequest if improper submission
        run_ora_val = process_files_run_ora(current_user, user_email, form, db_form, parameters_form)

        if isinstance(run_ora_val, str):
            return HttpResponseBadRequest(run_ora_val)

        elif isinstance(run_ora_val, bool):
            messages.error(request, DB_SELECTION_MSG)
            return redirect("/")
        else:
            _update_job_user(request, current_user, run_ora_val[0], run_ora_val[1])

    # Check if user submits form to run GSEA
    else:
        form = UploadGSEAForm(request.POST, request.FILES)
        db_form = GSEADatabaseSelectionForm(request.POST, request.FILES)
        parameters_form = GSEAParametersForm(request.POST, request.FILES)
        fc_form = UploadFoldChangesFormGSEA(request.POST, request.FILES)

        # Check if form is valid
        if not form.is_valid():
            return HttpResponseBadRequest(FORM_COMPLETION_MSG)

        # Check if form is valid
        if not db_form.is_valid():
            return HttpResponseBadRequest(DB_SELECTION_MSG)

        # Check if form is valid
        if not parameters_form.is_valid():
            return HttpResponseBadRequest(MAX_MIN_GENES_MSG)

        # Check if form is valid
        if not fc_form.is_valid():
            return HttpResponseBadRequest(INVALID_FOLD_CHANGES_MSG)

        # If form is valid, process file data in request.FILES and throw HTTPBadRequest if improper submission
        run_gsea_val = process_files_run_gsea(current_user, user_email, form, db_form, parameters_form, fc_form)

        if isinstance(run_gsea_val, str):
            return HttpResponseBadRequest(run_gsea_val)

        elif isinstance(run_gsea_val, bool):
            messages.error(request, DB_SELECTION_MSG)
            return redirect("/")

        else:
            _update_job_user(request, current_user, run_gsea_val[0], run_gsea_val[1])

    return redirect("/experiments")


@login_required
def experiments(request):
    """Render user summary page with information on results of each experiment."""
    current_user = request.user

    if request.method == "POST":
        # Get details on user experiment to check POST Query
        summary_table_body, objs = create_summary_table(current_user)

        for idx, obj in enumerate(objs):
            if request.POST.get(f"Delete_{idx}"):
                obj.delete()
            if request.POST.get(f"Stop_{idx}"):
                task = AsyncResult(id=obj.get_task_id())
                task.revoke(terminate=True)
                parent = task.parent
                while parent is not None:
                    parent.revoke(terminate=True)
                    parent = parent.parent

                obj.result_status = 0
                obj.error_message = "Stopped by User."
                obj.save()

    summary_table_body, objs = create_summary_table(current_user)

    context = {"headers": RESULTS_METADATA_HEADER, "body": summary_table_body}

    return render(request=request, template_name="viewer/experiments.html", context=context)


# TODO: Add it to celery task
def _update_job_user(request, current_user, job, task):
    job.task_id = task.id
    job.save()

    User.objects.filter(email__exact=current_user.email).update(num_of_jobs=F('num_of_jobs') + 1)

    if current_user.get_num_of_jobs == MAX_JOBS:
        User.objects.filter(email__exact=current_user.email).update(has_exceeded_quota=True)
        messages.warning(request=request, message=EXCEED_JOBS_MSG)


"""Generate result views"""


@login_required
def results(request, result_id):
    """Render results page."""
    current_user = request.user

    # Get results dataframe and databases from results model
    query_results_val = query_results_model(result_id, current_user)

    # Check if file cannot be read and throw error
    if isinstance(query_results_val, str):
        return HttpResponseBadRequest(query_results_val)
    else:
        df, _, _, enrichment_method, data_filename, _ = query_results_val

    # Get results table to display
    data = get_ranking_table(df)

    context = {
        'ranking_table': data,
        'result_id': result_id,
        'enrichment_method': enrichment_method,
        'data_filename': data_filename
    }

    return render(request, 'viewer/results.html', context=context)


@login_required
def gsea_results(request, result_id):
    """Render results page."""
    current_user = request.user

    # Get results dataframe and databases from results model
    query_results_val = query_results_model(result_id, current_user)

    # Check if file cannot be read and throw error
    if isinstance(query_results_val, str):
        return HttpResponseBadRequest(query_results_val)
    else:
        df, _, _, enrichment_method, data_filename, _ = query_results_val

    # Get results table to display
    data_dict, databases = get_ranking_table(df)

    context = {
        'ranking_tables': data_dict,
        'databases': databases,
        'result_id': result_id,
        'enrichment_method': enrichment_method,
        'data_filename': data_filename,
    }

    return render(request, 'viewer/gsea_results.html', context=context)


@login_required
def dc_genesets(request, pathway_id):
    """Render genesets page."""

    if os.path.isfile(DECOPATH_GMT):
        geneset_dict = _get_gmt_dict(DECOPATH_GMT)
    else:
        handle_file_download(DEFAULT_DATABASES_URLS['DECOPATH'], DECOPATH_GMT)
        geneset_dict = parse_gmt_file(DECOPATH_GMT)

    # Get decopath source pathway resources
    try:
        dc_source_dict, dc_source_id_dict = get_dc_pathway_resources(DECOPATH_ONTOLOGY)

    except IOError:
        return HttpResponseBadRequest(DOWNLOAD_ONTOLOGY_MSG)

    # Get decopath super pathway name
    try:
        with open(DECOPATH_NAME_MAPPINGS) as json_file:
            dc_id_name_mappings = json.load(json_file)

    except IOError:
        return HttpResponseBadRequest(DOWNLOAD_DECOPATH_NAMES)

    resource_list = list(dc_source_dict[pathway_id])
    resource_list.sort()
    resources = '|'.join(DATABASES[database] for database in resource_list)

    source_ids_list = sorted(list(dc_source_id_dict[pathway_id]), key=str.casefold)

    pathway_name = dc_id_name_mappings[pathway_id]
    dc_geneset_size = (len(geneset_dict[pathway_id]))
    dc_geneset = json.dumps(geneset_dict[pathway_id])

    context = {
        'pathway_name': pathway_name,
        'pathway_id': pathway_id,
        'geneset': dc_geneset,
        'resources': resources,
        'source_ids': source_ids_list,
        'geneset_size': dc_geneset_size,
    }

    return render(request, "viewer/dc_genesets.html", context)


@login_required
def consensus_ora(request, result_id):
    """Render consensus page for ORA."""
    current_user = request.user

    # Get results dataframe and databases from results model
    query_results_val = query_results_model(result_id, current_user)

    # Check if file cannot be read and throw error
    if isinstance(query_results_val, str):
        return HttpResponseBadRequest(query_results_val)
    else:
        df, databases, significance_value, enrichment_method, data_filename, symbol_to_fold_change = query_results_val

    (table_header,
     table_body,
     metadata_ids,
     metadata_consensus,
     metadata_dc_consensus,
     full_consensus_df,
     df_list,
     consensus_dict,
     ) = generate_consensus_table_ora(df, significance_value)

    header = json.dumps(table_header)
    identifiers = json.dumps(metadata_ids)
    consensus = json.dumps(metadata_consensus)
    dc_consensus = json.dumps(metadata_dc_consensus)
    df_len = json.dumps(df_list)
    pie_chart_data = json.dumps([
        [k, len(v)]
        for k, v in consensus_dict.items()
    ])

    context = {
        'table_header': table_header,
        'table_body': table_body,
        'header': header,
        'metadata_ids': identifiers,
        'default_qval': PADJ,
        'metadata_consensus': consensus,
        'metadata_dc_consensus': dc_consensus,
        'enrichment_method': enrichment_method,
        'data_filename': data_filename,
        'full_consensus_df': full_consensus_df,
        'df_list': df_len,
        'pie_chart_data': pie_chart_data,
    }

    return render(request, "viewer/consensus_ora.html", context=context)


@login_required
def consensus_gsea(request, result_id):
    """Render consensus page for GSEA."""
    current_user = request.user

    # Get results dataframe and databases from results model
    query_results_val = query_results_model(result_id, current_user)

    # Check if file cannot be read and throw error
    if isinstance(query_results_val, str):
        return HttpResponseBadRequest(query_results_val)
    else:
        df, databases, significance_value, enrichment_method, data_filename, symbol_to_fold_change = query_results_val

    (
        table_header,
        table_body,
        metadata_ids,
        metadata_qvals,
        metadata_consensus,
        metadata_dc_consensus,
        full_consensus_df,
        df_list,
        consensus_dict,
    ) = generate_consensus_table_gsea(df, significance_value)

    header = json.dumps(table_header)
    identifiers = json.dumps(metadata_ids)
    qvals = json.dumps(metadata_qvals)
    consensus = json.dumps(metadata_consensus)
    dc_consensus = json.dumps(metadata_dc_consensus)
    df_len = json.dumps(df_list)
    pie_chart_data = json.dumps([
        [k, len(v)]
        for k, v in consensus_dict.items()
    ])

    context = {
        'table_header': table_header,
        'table_body': table_body,
        'header': header,
        'metadata_ids': identifiers,
        'metadata_qvals': qvals,
        'default_qval': PADJ,
        'metadata_consensus': consensus,
        'metadata_dc_consensus': dc_consensus,
        'enrichment_method': enrichment_method,
        'data_filename': data_filename,
        'full_consensus_df': full_consensus_df,
        'df_list': df_len,
        'pie_chart_data': pie_chart_data,
    }

    return render(request, "viewer/consensus_gsea.html", context=context)


@login_required
def circles_viz(request, result_id):
    """Render circles viz."""
    current_user = request.user

    # Get results dataframe and databases from results model
    query_results_val = query_results_model(result_id, current_user)

    # Check if file cannot be read and throw error
    if isinstance(query_results_val, str):
        return HttpResponseBadRequest(query_results_val)
    else:
        df, dbs, significance_value, enrichment_method, data_filename, symbol_to_fold_change = query_results_val

    # Access pathway hierarchy
    try:
        pathway_hierarchy = PathwayHierarchy.objects.get(name='default')
        pathway_network = pathway_hierarchy.get_network()
        composite_pathway = pathway_hierarchy.super_pathway
    except ObjectDoesNotExist:
        raise ValueError('Please run "python manage.py load_db" to load the database')

    pathway_results = get_results_circle_viz(df, dbs, enrichment_method)

    try:
        tree_json = map_results_to_hierarchy(
            network_hierarchy=pathway_network,
            root_node=composite_pathway,
            results=pathway_results,
            enrichment_method=enrichment_method,
            significance_value=significance_value
        )
    except ValueError:
        return HttpResponseBadRequest('Could not map to any pathways.')

    return render(
        request,
        "viewer/viz/circles.html",
        context={
            'tree_json': tree_json,
            'result_id': result_id,
            'enrichment_method': enrichment_method,
            'data_filename': data_filename,
        }
    )


@login_required
def zoom_in(request, result_id, pathway_id):
    """Render zoom-in page."""
    query_results_val = query_results_model(result_id, request.user)

    # Check if file cannot be read and throw error
    if isinstance(query_results_val, str):
        return HttpResponseBadRequest(query_results_val)
    else:
        _, databases, _, _, _, symbol_to_fold_change = query_results_val
    # Get changes
    if symbol_to_fold_change:
        symbol_to_fold_change = pickle.loads(symbol_to_fold_change)
        if isinstance(symbol_to_fold_change, DataFrame):
            if 'log2FoldChange' in symbol_to_fold_change.columns:
                symbol_to_fold_change.rename(columns={
                    'log2FoldChange': 'log2fc',
                }, inplace=True)

            symbol_to_fold_change = {
                gene_symbol: row['log2fc']
                for gene_symbol, row in symbol_to_fold_change.iterrows()
                if notna(row['log2fc'])
            }
        else:
            raise ValueError('Something went wrong.')
    else:
        symbol_to_fold_change = {}

    venn_diagram_data = process_overlap_for_venn_diagram(
        pathway_id=pathway_id,
        databases=databases,
    )

    # Return HTTP bad request if pathway ID does not exist in equivalent pathway databases
    if isinstance(venn_diagram_data, str):
        return HttpResponseBadRequest(venn_diagram_data)

    venn_diagram = venn_diagram_data

    return render(
        request,
        "viewer/viz/zoom_in.html",
        context={
            'venn_diagram_data': SafeString(venn_diagram),
            'fold_changes': symbol_to_fold_change,
        }
    )


"""Custom user gene sets and mappings"""


def custom_databases(request):
    """Render custom gene set submission page."""
    # Initialize formset with 1 empty form
    extra = 1

    # Create a custom gene set formset class
    CustomFormSet = formset_factory(UploadCustomDataForm, extra=extra)

    # Check if request method is POST
    if request.method == 'POST':

        # Check if button clicked to add more forms to formset
        if 'addfile' in request.POST and request.POST['addfile'] == 'true':
            formset = add_forms_to_formset(request, CustomFormSet, extra)

        # Handle actual form submission
        else:
            # POST; populate form with data from the request (standard submit)
            formset = CustomFormSet(request.POST, request.FILES)

            # Check if the formset is valid
            if formset.is_valid():

                for form in formset:

                    # Get cleaned data
                    database_name = form.cleaned_data['database']
                    cleaned_gmt_file = form.cleaned_data['gene_set']

                    # Get paths to temporary uploaded files
                    gmt_path = cleaned_gmt_file.transient_file_path()

                    # Check if GMT file is valid
                    check_gmt = parse_custom_gmt(gmt_path)

                    # Check if file cannot be read and throw error
                    if isinstance(check_gmt, str):
                        return HttpResponseBadRequest(check_gmt)

                    # Check if database already exists in PathwayDatabase model
                    pathway_db_queryset = PathwayDatabase.objects.filter(database_name=database_name)
                    if pathway_db_queryset:
                        return HttpResponseBadRequest(DB_PRELOADED_ERROR)

                    # Check if custom databases are concordant with custom mappings
                    pathway_queryset = Pathway.objects.all()

                    set_mapping_dbs = set(
                        obj.pathway_database
                        for obj in pathway_queryset
                    )

                    if database_name not in set_mapping_dbs:
                        return HttpResponseBadRequest(CUSTOM_MAPPING_DATABASES_MSG)

                    # Load custom pathway database gene sets
                    load_custom_pathway_databases(database_name, check_gmt, gmt_path)

                # Redirect to custom mappings submission page
                return redirect("/")

            else:
                return HttpResponseBadRequest(GMT_FILE_ERROR_MSG)

    # """GET; first time you access the query."""
    else:
        # Generate blank submission form
        formset = CustomFormSet()

    return render(request, "viewer/custom_databases.html", {'formset': formset})


def mappings(request):
    """Render custom mappings submission page."""
    # Initialize formset with 1 empty form
    extra = 1

    # Create a mapping formset class
    MappingFormSet = formset_factory(UploadMappingsForm, extra=extra)

    # Check if request method is POST
    if request.method == 'POST':

        # Check if button clicked to add more forms to formset
        if 'addfile' in request.POST and request.POST['addfile'] == 'true':
            formset = add_forms_to_formset(request, MappingFormSet, extra)

        # Handle actual form submission
        else:
            # POST; populate form with data from the request (standard submit)
            formset = MappingFormSet(request.POST, request.FILES)

            # Check if the formset is valid
            if formset.is_valid():

                for form in formset:

                    result = check_mapping_validity(form)

                    # Check if mapping file is valid before loading into database
                    if isinstance(result, str):
                        return HttpResponseBadRequest(result)

                    load_custom_pathway_mappings(result)

                # Redirect to custom database submission page
                return redirect("custom_databases")

            else:
                return HttpResponseBadRequest(CUSTOM_MAPPINGS_MSG)

    # """GET; first time you access the query."""
    else:
        # Generate blank submission form
        formset = MappingFormSet()

    return render(request, "viewer/mappings.html", {'formset': formset})
