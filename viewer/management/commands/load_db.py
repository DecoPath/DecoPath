# -*- coding: utf-8 -*-

"""Command to load the database with app data models"""

import logging
import os

from django.core.management.base import BaseCommand

from viewer.models import PathwayDatabase, User
from viewer.src.constants import LOGS_DIR
from viewer.src.db_utils import (
    load_standard_pathway_databases,
    load_standard_pathway_mappings,
    erase_db,
    add_tree_to_database,
    load_decopath_pathway_databases,
)

# TODO: Check if database is already populated.
class Command(BaseCommand):
    help = 'Load DecoPath database'

    def add_arguments(self, parser):
        parser.add_argument(
            '--keep_db',
            help='Do not drop DB',
            type=bool
        )

    def handle(self, *args, **options):
        """--keep_db=True avoid dropping the db"""
        # Create logger
        logger = logging.getLogger('decopath.parser')
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(os.path.join(LOGS_DIR, 'decopath.log'), mode='a')
        fh.setLevel(logging.DEBUG)

        date_fmt = '%m/%d/%Y %I:%M:%S %p'

        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', date_fmt)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        logger.warning('Deleting existing database')
        erase_db()
        logger.warning('Finished deleting existing database')

        logger.info('Populating PathwayDatabase table')
        load_standard_pathway_databases()
        logger.info('Finished populating PathwayDatabase table')

        logger.info('Populating PathwayDatabase table with DecoPath super pathways')
        load_decopath_pathway_databases()
        logger.info('Finished populating PathwayDatabase table with DecoPath super pathways')

        logger.info('Populating Pathway table with ComPath mappings')
        load_standard_pathway_mappings()
        logger.info('Finished populating Pathway table with ComPath mappings')

        logger.info(f'{len(PathwayDatabase.objects.all())} databases have been loaded')

        # TODO: load gene set size and load MPath equivalent representations with decopath ids
        add_tree_to_database()
        logger.info('Databases have been loaded')
