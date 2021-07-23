# -*- coding: utf-8 -*-

"""Models module."""

import json
import pickle

from django.contrib.auth.base_user import BaseUserManager
from django.contrib.auth.models import AbstractBaseUser
from django.db import models
from django.utils import timezone


class PathwayDatabase(models.Model):
    """Pathway Database class storing GMT files."""
    database_name = models.CharField(max_length=360)
    version = models.DateTimeField(default=timezone.now)
    gene_set = models.BinaryField()
    gmt_file = models.BinaryField(null=True, blank=True)

    def __str__(self):
        return self.database_name

    def attrs(self):
        """Dictionary with instance attributes."""
        for attr, value in self.__dict__.items():
            yield attr, value

    @property
    def number_of_pathways(self):
        """Return the number of gene sets in gmt file."""
        return len(pickle.loads(self.gene_set))


class Pathway(models.Model):
    """Pathway class storing pathway mappings between different databases."""
    pathway_id = models.CharField(max_length=360)
    pathway_name = models.CharField(max_length=360)
    pathway_database = models.CharField(max_length=360)
    decopath_id = models.CharField(max_length=360, default="NA")
    decopath_name = models.CharField(max_length=360, default="NA")

    # Create a recursive relationship
    mapping_pathway = models.ManyToManyField('self', blank=True, symmetrical=False, related_name='equivalent_pathway')

    def __str__(self):
        return f'{self.pathway_id} {self.pathway_name} {self.pathway_database}'

    def attrs(self):
        """Dictionary with instance attributes."""
        for attr, value in self.__dict__.items():
            yield attr, value

    def number_of_equivalent_pathways(self):
        """Return the number of equivalent pathways."""
        return len(self.mapping_pathway.all())

    def get_equivalent_pathways(self):
        """Admin method."""
        return " |||| ".join(
            [
                mapping.__str__()
                for mapping in self.mapping_pathway.all()
            ]
        )


class UserManager(BaseUserManager):
    def create_user(self, email, password=None):
        """Creates and saves a User with the given email and password."""
        if not email:
            raise ValueError('Users must have an email address')

        user = self.model(
            email=self.normalize_email(email),
        )

        user.set_password(password)
        user.save(using=self._db)
        return user

    def create_superuser(self, email, password=None):
        """Creates and saves a User with the given email and password."""
        user = self.create_user(
            email,
            password=password,
        )
        user.is_admin = True
        user.save(using=self._db)
        return user


class User(AbstractBaseUser):
    """User class."""

    email = models.EmailField(
        verbose_name='email address',
        max_length=255,
        unique=True,
        db_index=True,
    )
    num_of_jobs = models.IntegerField(default=0)
    has_exceeded_quota = models.BooleanField(default=False)
    is_active = models.BooleanField(default=True)
    is_admin = models.BooleanField(default=False)
    verified = models.BooleanField(default=False)

    objects = UserManager()

    USERNAME_FIELD = 'email'

    def __str__(self):
        return f"{self.email} ({'Staff' if self.is_staff else 'Normal User'})"

    def __unicode__(self):
        return f"{self.email}"

    def get_email(self):
        """Return the email of a given user"""
        return f"{self.email}"

    def get_num_of_jobs(self):
        """Return the number of jobs of a given user"""
        return self.num_of_jobs

    def is_quota_left(self):
        return not self.has_exceeded_quota

    def is_verified(self):
        """Returns a boolean value to indicate if the user's email verified"""
        return self.verified

    def has_perm(self, perm, obj=None):
        """Does the user have a specific permission?"""
        return True

    def has_module_perms(self, app_label):
        """Does the user have permissions to view the app `app_label`?"""
        return True

    @property
    def is_staff(self):
        """Is the user a member of staff?"""
        return self.is_admin


class EnrichmentResult(models.Model):
    """Results class storing enrichment results matrix."""
    result_id = models.AutoField(primary_key=True)  # (Job ID)
    user = models.ForeignKey(User, on_delete=models.CASCADE)
    date = models.DateTimeField(default=timezone.now)
    result = models.BinaryField(null=True, blank=False)
    result_status = models.IntegerField(default=1)  # 0 means failed, 1 means processing and 2 is success
    error_message = models.CharField(max_length=1000, default="NA")
    enrichment_method = models.CharField(max_length=360, default="NA")  # GSEA, GSEA Pre-Ranked, ORA
    phenotype_classes = models.JSONField(default=list)
    data_filename = models.CharField(max_length=360, default="NA")
    class_filename = models.CharField(max_length=360, default="NA")
    sample_number = models.IntegerField(null=True)
    databases = models.JSONField(default=list)
    significance_threshold = models.FloatField(null=True)
    min_genes = models.IntegerField(null=True)
    max_genes = models.IntegerField(null=True)
    permutation_number = models.IntegerField(null=True)
    permutation_type = models.CharField(max_length=360, default="NA")
    calculation_method = models.CharField(max_length=360, default="NA")
    significance_threshold_fc = models.FloatField(null=True)
    fold_change_results = models.BinaryField(null=True, blank=False)
    fold_changes_filename = models.CharField(max_length=360, default="NA")
    task_id = models.CharField(max_length=450, null=True, blank=False)  # Celery task ID

    def __str__(self):
        return f'results from {self.user} on {self.date}'

    def attrs(self):
        """Dictionary with instance attributes."""
        for attr, value in self.__dict__.items():
            yield attr, value

    def get_df(self):
        return pickle.loads(self.result)

    def get_job_id(self):
        return f'{self.result_id}'

    def get_task_id(self):
        return f'{self.task_id}'

    def get_databases(self):
        return json.loads(self.databases)

    def get_phenotypes(self):
        return json.loads(self.phenotype_classes)

    def get_fold_changes(self):
        df = pickle.loads(self.fold_change_results)

        if not df:
            return None

        if 'log2fc' in df.columns and 'gene_symbol' in df.columns:
            return {
                row['gene_symbol']: df['log2fc']
                for index, row in df.iterrows()
            }
        else:
            raise ValueError('Your file is missing some information. Please ensure the following columns exist: '
                             'gene_symbol, log2fc and q-value. See the FAQs for more details.')


class PathwayHierarchy(models.Model):
    """Pathway hierarchy to be rendered."""
    name = models.CharField(max_length=120)
    json_tree = models.BinaryField(null=True, blank=False)
    network = models.BinaryField(null=True, blank=False)
    equivalent_pathways = models.BinaryField(null=True, blank=False)
    super_pathway = models.CharField(max_length=120)

    def get_json_tree(self):
        return pickle.loads(self.json_tree)

    def get_network(self):
        return pickle.loads(self.network)

    def get_equivalent_pathways(self):
        return pickle.loads(self.equivalent_pathways)
