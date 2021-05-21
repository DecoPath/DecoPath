# -*- coding: utf-8 -*-

"""Utilities to load the DecoPath models."""

import json
import pickle
from typing import List, Optional

import pandas as pd

from viewer.models import PathwayDatabase, EnrichmentResult, Pathway, PathwayHierarchy
from viewer.src.constants import *
from viewer.src.data_preprocessing import parse_gmt_file
from viewer.src.utils import (
    handle_file_download,
    parse_hierarchy_excel,
    get_equivalent_pathway_dc_ids,
    get_decopath_genesets,
    export_geneset
)

"""Populate pathway database table with gene sets and gmt files per database."""


def load_standard_pathway_databases(database_to_gmt: dict = DEFAULT_DATABASES_URLS):
    """Load the pathway database model with the standard pathway databases."""
    for database_name, url in database_to_gmt.items():
        filename = database_name + GMT_FILE_EXTENSION
        file_path = os.path.join(GMT_FILES_DIR, filename)

        # Download gmt from url and save to directory
        handle_file_download(url, file_path)

        # Get dictionary of pathways and corresponding gene sets
        geneset_dict = parse_gmt_file(file_path)

        # Load Pathway Database model with databases, gene set dictionaries and gmt files from ComPath
        pathway_database_object, created = PathwayDatabase.objects.get_or_create(
            database_name=database_name,
            gene_set=pickle.dumps(geneset_dict),
            gmt_file=pickle.dumps(file_path),
        )

        pathway_database_object.save()


def load_decopath_pathway_databases(
    decopath_ontology: str = DECOPATH_ONTOLOGY,
    gmt_dir: str = GMT_FILES_DIR,
    database: str = DECOPATH,
    outfile: str = DECOPATH_CSV,
    gmt_file: str = DECOPATH_GMT
):
    """Load Pathway Database model with databases, gene set dictionaries and gmt files from the DecoPath ontology."""
    decopath_genesets = get_decopath_genesets(decopath_ontology, gmt_dir)
    export_geneset(decopath_genesets, database, outfile, gmt_file)

    pathway_database_object, created = PathwayDatabase.objects.get_or_create(
        database_name=database,
        gene_set=pickle.dumps(decopath_genesets),
        gmt_file=pickle.dumps(gmt_file),
    )
    pathway_database_object.save()


"""Populate pathway mapping table."""


def load_pathway_mappings(df: pd.DataFrame):
    """Load pathway model with equivalent pathway mappings and DC IDs."""
    source_dict = {}
    target_dict = {}
    decopath_dict = {}

    # Convert dataFrame to a list where each row is dictionary
    for mapping_dict in df.to_dict('records'):
        source_dict['resource'] = mapping_dict[SOURCE_RESOURCE]
        source_dict['id'] = mapping_dict[SOURCE_ID].replace('path:', '')
        source_dict['name'] = mapping_dict[SOURCE_NAME]

        decopath_dict[source_dict['id']] = {
            'decopath_id': mapping_dict['dc_id'],
            'decopath_name': mapping_dict['dc_name']
        }

        target_dict['resource'] = mapping_dict[TARGET_RESOURCE]
        target_dict['id'] = mapping_dict[TARGET_ID].replace('path:', '')
        target_dict['name'] = mapping_dict[TARGET_NAME]

        pathway_source_object, created = Pathway.objects.get_or_create(
            **{
                'pathway_id': source_dict['id'],
                'pathway_name': source_dict['name'],
                'pathway_database': source_dict['resource'],
                'decopath_id': mapping_dict['dc_id'],
                'decopath_name': mapping_dict['dc_name']
            }
        )

        pathway_target_object, created = Pathway.objects.get_or_create(
            **{
                'pathway_id': target_dict['id'],
                'pathway_name': target_dict['name'],
                'pathway_database': target_dict['resource'],
                'decopath_id': mapping_dict['dc_id'],
                'decopath_name': mapping_dict['dc_name']
            }
        )

        pathway_target_object.mapping_pathway.add(pathway_source_object)
        pathway_target_object.save()
        pathway_source_object.mapping_pathway.add(pathway_target_object)
        pathway_source_object.save()


def load_standard_pathway_mappings(decopath_ontology: str = DECOPATH_ONTOLOGY):
    """Get pairwise pathway mappings from DecoPath ontology and load pathway model."""
    # Get DC IDs for equivalent pathways
    id_to_dc_id, df = get_equivalent_pathway_dc_ids(decopath_ontology)

    # Get equivalent pathways
    equivalence_df = df.loc[df[MAPPING_TYPE] == EQUIVALENT_TO]

    # Merge equivalent mappings dataFrame with DC IDs dataFrame
    merged_df = pd.merge(equivalence_df, id_to_dc_id, left_on=SOURCE_ID, right_on='pathway_id')

    # Load pathway model with equivalent mappings
    load_pathway_mappings(merged_df)


"""Handle Enrichment Results model loading."""


def load_results_metadata(
    results: Optional[pd.DataFrame],
    current_user,
    data_filename: str,
    class_filename: str,
    class_labels: List,
    databases: List,
    min_size: int,
    max_size: int,
    significance_threshold: Optional[float],
    permutation_type: str,
    calculation_method: str,
    sample_number: Optional[int],
    permutation_number: Optional[int],
    enrichment_method: str,
    significance_threshold_fc: Optional[float],
    fold_changes_filename: str,
    fold_change_results: Optional[pd.DataFrame],
    read_counts_path: Optional[str],
) -> EnrichmentResult:
    """Get metadata from files submitted by user and load enrichment results model."""
    # Load Enrichment Results model with results and fold changes uploaded by user
    if isinstance(fold_change_results, pd.DataFrame) and isinstance(results, pd.DataFrame):
        enrichment_results_object, created = EnrichmentResult.objects.get_or_create(
            result=pickle.dumps(results),
            user=current_user,
            data_filename=data_filename,
            class_filename=class_filename,
            sample_number=sample_number,
            phenotype_classes=json.dumps(class_labels),
            databases=json.dumps(databases),
            min_genes=min_size,
            max_genes=max_size,
            significance_threshold=significance_threshold,
            permutation_number=permutation_number,
            permutation_type=permutation_type,
            calculation_method=calculation_method,
            enrichment_method=enrichment_method,
            significance_threshold_fc=significance_threshold_fc,
            fold_change_results=pickle.dumps(fold_change_results),
            fold_changes_filename=fold_changes_filename,
        )
    # Load Enrichment Results model with results and metadata to run DEG uploaded by user
    elif read_counts_path and isinstance(results, pd.DataFrame):
        enrichment_results_object, created = EnrichmentResult.objects.get_or_create(
            result=pickle.dumps(results),
            user=current_user,
            data_filename=data_filename,
            class_filename=class_filename,
            sample_number=sample_number,
            phenotype_classes=json.dumps(class_labels),
            databases=json.dumps(databases),
            min_genes=min_size,
            max_genes=max_size,
            significance_threshold=significance_threshold,
            permutation_number=permutation_number,
            permutation_type=permutation_type,
            calculation_method=calculation_method,
            enrichment_method=enrichment_method,
            significance_threshold_fc=significance_threshold_fc,
            fold_changes_filename=fold_changes_filename,
        )

    elif isinstance(results, pd.DataFrame) and not (read_counts_path or fold_change_results):
        enrichment_results_object, created = EnrichmentResult.objects.get_or_create(
            result=pickle.dumps(results),
            user=current_user,
            phenotype_classes=json.dumps(class_labels),
            data_filename=data_filename,
            class_filename=class_filename,
            sample_number=sample_number,
            databases=json.dumps(databases),
            max_genes=max_size,
            min_genes=min_size,
            significance_threshold=significance_threshold,
            permutation_number=permutation_number,
            permutation_type=permutation_type,
            calculation_method=calculation_method,
            enrichment_method=enrichment_method,
            significance_threshold_fc=significance_threshold_fc,
            fold_changes_filename=fold_changes_filename,
        )

    elif isinstance(fold_change_results, pd.DataFrame) and results is None:
        enrichment_results_object, created = EnrichmentResult.objects.get_or_create(
            user=current_user,
            phenotype_classes=json.dumps(class_labels),
            data_filename=data_filename,
            class_filename=class_filename,
            sample_number=sample_number,
            databases=json.dumps(databases),
            max_genes=max_size,
            min_genes=min_size,
            significance_threshold=significance_threshold,
            permutation_number=permutation_number,
            permutation_type=permutation_type,
            calculation_method=calculation_method,
            enrichment_method=enrichment_method,
            significance_threshold_fc=significance_threshold_fc,
            fold_change_results=pickle.dumps(fold_change_results),
            fold_changes_filename=fold_changes_filename,
        )

    elif read_counts_path and results is None:
        enrichment_results_object, created = EnrichmentResult.objects.get_or_create(
            user=current_user,
            phenotype_classes=json.dumps(class_labels),
            data_filename=data_filename,
            class_filename=class_filename,
            sample_number=sample_number,
            databases=json.dumps(databases),
            max_genes=max_size,
            min_genes=min_size,
            significance_threshold=significance_threshold,
            permutation_number=permutation_number,
            permutation_type=permutation_type,
            calculation_method=calculation_method,
            enrichment_method=enrichment_method,
            significance_threshold_fc=significance_threshold_fc,
            fold_changes_filename=fold_changes_filename,
        )

    # Load Enrichment Result model with metadata
    else:
        enrichment_results_object, created = EnrichmentResult.objects.get_or_create(
            user=current_user,
            phenotype_classes=json.dumps(class_labels),
            data_filename=data_filename,
            class_filename=class_filename,
            sample_number=sample_number,
            databases=json.dumps(databases),
            significance_threshold=significance_threshold,
            max_genes=max_size,
            min_genes=min_size,
            permutation_number=permutation_number,
            permutation_type=permutation_type,
            calculation_method=calculation_method,
            enrichment_method=enrichment_method,
        )

    enrichment_results_object.save()
    return enrichment_results_object


def get_results_per_user(current_user):
    """Get user specific results."""
    objects = EnrichmentResult.objects.filter(user=current_user)

    return [
        experiment.get_df()
        for experiment in objects
    ]


"""Load Pathway hierarchy model."""


def add_tree_to_database():
    """Insert tree-hierarchy to the database."""
    tree, tree_network, equivalent_pathways, super_pathway = parse_hierarchy_excel(DECOPATH_ONTOLOGY)

    PathwayHierarchy.objects.get_or_create(
        name='default',
        json_tree=pickle.dumps(tree),
        network=pickle.dumps(tree_network),
        equivalent_pathways=pickle.dumps(equivalent_pathways),
        super_pathway=super_pathway,
    )


"""Load custom databases and mappings."""


def load_custom_pathway_databases(database_name: str, geneset_dict: dict, file_path: str):
    """Load the pathway database model with custom, user-defined pathway databases."""

    pathway_database_object, created = PathwayDatabase.objects.get_or_create(
        database_name=database_name,
        gene_set=pickle.dumps(geneset_dict),
        gmt_file=pickle.dumps(file_path),
    )

    pathway_database_object.save()


def load_custom_pathway_mappings(df: pd.DataFrame):
    """Load pathway model with equivalent pathway mappings and DC IDs."""
    source_dict = {}
    target_dict = {}

    # Convert dataFrame to a list where each row is dictionary
    for mapping_dict in df.to_dict('records'):
        source_dict['resource'] = mapping_dict[SOURCE_RESOURCE]
        source_dict['id'] = mapping_dict[SOURCE_ID].replace('path:', '')
        source_dict['name'] = mapping_dict[SOURCE_NAME]

        target_dict['resource'] = mapping_dict[TARGET_RESOURCE]
        target_dict['id'] = mapping_dict[TARGET_ID].replace('path:', '')
        target_dict['name'] = mapping_dict[TARGET_NAME]

        pathway_source_object, created = Pathway.objects.get_or_create(
            **{
                'pathway_id': source_dict['id'],
                'pathway_name': source_dict['name'],
                'pathway_database': source_dict['resource'],
            }
        )

        pathway_target_object, created = Pathway.objects.get_or_create(
            **{
                'pathway_id': target_dict['id'],
                'pathway_name': target_dict['name'],
                'pathway_database': target_dict['resource'],
            }
        )

        pathway_target_object.mapping_pathway.add(pathway_source_object)
        pathway_target_object.save()
        pathway_source_object.mapping_pathway.add(pathway_target_object)
        pathway_source_object.save()


"""Delete existing database."""


def erase_db():
    """Drops all objects from the database."""
    # TODO: Add new models here. E.g: model_name.objects.all().delete()
    PathwayDatabase.objects.all().delete()
    EnrichmentResult.objects.all().delete()
    Pathway.objects.all().delete()
    PathwayHierarchy.objects.all().delete()
