# -*- coding: utf-8 -*-

"""Command to initialise users in the database with app user model."""

import logging
import os

from django.core.management.base import BaseCommand
from django.conf import settings

from viewer.models import User
from viewer.src.constants import LOGS_DIR


class Command(BaseCommand):
    help = 'Initialise admin and test user in DecoPath.'

    def add_arguments(self, parser):
        parser.add_argument(
            '--superuser', default=None, type=str, help='Email for the superuser.',
        )

    def handle(self, *args, **options):
        # Create logger
        logger = logging.getLogger('decopath.parser')
        logger.setLevel(logging.DEBUG)
        fh = logging.FileHandler(os.path.join(LOGS_DIR, 'decopath.log'), mode='a')
        fh.setLevel(logging.DEBUG)

        date_fmt = '%m/%d/%Y %I:%M:%S %p'

        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s', date_fmt)
        fh.setFormatter(formatter)
        logger.addHandler(fh)

        # Add Admin
        admin_mail = options.get('superuser')
        if admin_mail:
            if not User.objects.filter(email__exact=admin_mail).exists():
                admin = User.objects.create_superuser(admin_mail, password='admin')
                admin.save()

        # Add default user
        user_mail = getattr(settings, "USER_EMAIL", 'testuser@user.com')
        if not User.objects.filter(email__exact=user_mail).exists():
            user = User.objects.create_user(user_mail, password='password')
            user.save()
