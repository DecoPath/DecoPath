# -*- coding: utf-8 -*-

"""Forms module."""

from django import forms
from django.contrib.auth.admin import UserAdmin as BaseUserAdmin
from django.contrib.auth.forms import ReadOnlyPasswordHashField
from django.core import validators

from viewer.models import PathwayDatabase, User
from viewer.src.constants import (
    PERMUTATION_NUM,
    PERMUTATION_TYPE,
    METHOD,
    ENRICHMENT_METHOD,
    CSV,
    TSV,
    DECOPATH, DATABASES, REACTOME,
)

"""User uploaded results forms."""


class UploadResultsForm(forms.Form):
    """Upload results file form class"""
    select_enrichment_method = forms.ChoiceField(
        label='Enrichment method',
        choices=ENRICHMENT_METHOD,
    )

    results_file = forms.FileField(
        label='Upload results',
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Supported file formats:</strong> *.csv or *.tsv."
    )

    sig_threshold_results = forms.FloatField(
        label='Significance threshold (adjusted <i>p</i>-value)',
        required=False,
        help_text="Significance threshold to filter enrichment results. <strong>Default:</strong> 0.05"
    )


class UploadFoldChangesForm(forms.Form):
    """Upload run DGE analysis or submit fold changes form class."""
    run_analysis = forms.BooleanField(
        label='Run differential gene expression analysis',
        required=False,
    )

    read_counts_file = forms.FileField(
        label='Upload read counts',
        required=False,
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Note:</strong> this field is optional. If you would like to perform differential gene "
                  "expression analysis, please ensure a class labels file is uploaded below with samples in the same "
                  "order as the read counts file. The following file formats are supported: *.csv or *.tsv."
    )

    class_file = forms.FileField(
        label='Upload class labels',
        required=False,
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Note:</strong> this field is optional. Please ensure you have also uploaded a read counts"
                  " matrix. The following file formats are supported: *.csv or *.tsv."
    )

    significance_value_run = forms.FloatField(
        label='Significance threshold (adjusted <i>p</i>-value)',
        required=False,
        help_text="Significance threshold to filter differentially expressed genes. <strong>Default:</strong> 0.05"
    )

    upload_fold_changes = forms.BooleanField(
        label='Upload differential gene expression analysis results',
        required=False,
    )

    fold_changes_file = forms.FileField(
        label='Upload fold changes',
        required=False,
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Note:</strong> this field is optional. If you are uploading GSEA results, please ensure the "
                  "fold changes correspond to the same samples from your GSEA run. The following file formats are "
                  "supported: *.csv or *.tsv."
    )

    significance_value_upload = forms.FloatField(
        label='Significance threshold (adjusted <i>p</i>-value)',
        required=False,
        help_text="Significance threshold to filter differentially expressed genes. <strong>Default:</strong> 0.05"
    )


"""Run ORA forms."""


class ORADatabaseSelectionForm(forms.Form):
    """Select databases to run ORA form class."""
    queryset = PathwayDatabase.objects.exclude(database_name=DECOPATH)

    choices = tuple([
        (DATABASES[obj.database_name], DATABASES[obj.database_name])
        for obj in queryset
    ])

    select_databases = forms.MultipleChoiceField(
        choices=choices,
        widget=forms.CheckboxSelectMultiple,
        error_messages={
            'required': 'Please select at least two databases.'
        },
        required=True,
        help_text="<strong>Note:</strong> Please select a minimum of two databases."
    )


class ORAParametersForm(forms.Form):
    """ORA optional parameters form class."""
    minimum_size = forms.IntegerField(
        min_value=3,
        max_value=5000,
        label='Minimum number of genes in gene set',
        required=False,
        help_text="<strong>Default:</strong> 10"
    )

    maximum_size = forms.IntegerField(
        min_value=3,
        max_value=10000,
        label='Maximum number of genes in gene set',
        required=False,
        help_text="<strong>Default:</strong> 1000"
    )

    sig_threshold_ora = forms.FloatField(
        label='Significance threshold (adjusted <i>p</i>-value)',
        required=False,
        help_text="Significance threshold to filter ORA results. <strong>Default:</strong> 0.05"
    )


class OraOptions(forms.Form):
    """Upload ORA form class."""
    upload_gene_list = forms.BooleanField(
        label='Upload gene list and run ORA',
        required=False,
    )

    gene_list = forms.FileField(
        label='Upload gene list',
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        required=False,
        help_text="<strong>Supported file formats:</strong> *.csv or *.tsv."
    )

    run_dge_ora_genes = forms.BooleanField(
        label='Run ORA on differentially expressed genes',
        required=False,
    )

    run_analysis_results_ora = forms.BooleanField(
        label='Run differential gene expression analysis',
        required=False,
    )

    read_counts_file_ora = forms.FileField(
        label='Upload read counts',
        required=False,
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Note:</strong> this field is optional. If you would like to perform differential gene "
                  "expression analysis, please ensure a class labels file is uploaded below with samples in the same "
                  "order as the read counts file. The following file formats are supported: *.csv or *.tsv."
    )

    class_file_ora = forms.FileField(
        label='Upload class labels',
        required=False,
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Note:</strong> this field is optional. Please ensure you have also uploaded a read counts"
                  " matrix. The following file formats are supported: *.csv or *.tsv."
    )

    significance_value_run_ora = forms.FloatField(
        label='Significance threshold (adjusted <i>p</i>-value)',
        required=False,
        help_text="Significance threshold to filter differentially expressed genes. <strong>Default:</strong> 0.05"
    )

    upload_fold_changes_ora = forms.BooleanField(
        label='Upload differential gene expression analysis results',
        required=False,
    )

    fold_changes_file_ora = forms.FileField(
        label='Upload fold changes',
        required=False,
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Note:</strong> this field is optional. The following file formats are "
                  "supported: *.csv or *.tsv."
    )

    significance_value_upload_ora = forms.FloatField(
        label='Significance threshold (adjusted <i>p</i>-value)',
        required=False,
        help_text="Significance threshold to filter differentially expressed genes. <strong>Default:</strong> 0.05"
    )


"""Run GSEA forms."""


class GSEADatabaseSelectionForm(forms.Form):
    """Select databases to run GSEA form class."""
    queryset = PathwayDatabase.objects.exclude(database_name=DECOPATH)

    choices = tuple([
        (f'{DATABASES[obj.database_name]}', f'{DATABASES[obj.database_name]}**')
        if obj.database_name == REACTOME
        else (DATABASES[obj.database_name], DATABASES[obj.database_name])
        for obj in queryset
    ])

    select_databases = forms.MultipleChoiceField(
        choices=choices,
        widget=forms.CheckboxSelectMultiple,
        error_messages={
            'required': 'Please select at least two databases.'
        },
        required=True,
        help_text="<strong>*Note:</strong> Please select a minimum of two databases. <br><strong>**</strong> Due to the"
                  " high computational cost of running GSEA on the Reactome database, only mapped pathways are included"
                  " in the analysis. See the <a href='https://decopath.scai.fraunhofer.de/user_guide'>User Guide</a> "
                  "for alternatives to run GSEA on the full database."
    )


class UploadGSEAForm(forms.Form):
    """Upload GSEA forms class."""
    data_file = forms.FileField(
        label='Upload expression dataset',
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Supported file formats:</strong> *.csv or *.tsv."
    )

    class_file = forms.FileField(
        label='Upload class labels',
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Supported file formats:</strong> *.csv or *.tsv."
    )


class GSEAParametersForm(forms.Form):
    """GSEA optional parameters form class."""
    permutation_num = forms.ChoiceField(
        label='Number of permutations',
        choices=PERMUTATION_NUM,
        required=False,
        help_text="<strong>Default:</strong> 100"
    )

    permutation_type = forms.ChoiceField(
        label='Permutation type',
        choices=PERMUTATION_TYPE,
        required=False,
        help_text="<strong>Default:</strong> Phenotype"
    )

    method = forms.ChoiceField(
        label='Method for calculation of correlation or ranking',
        choices=METHOD,
        required=False,
        help_text="<strong>Default:</strong> Log2 ratio of classes"
    )

    minimum_size = forms.IntegerField(
        min_value=10,
        max_value=3000,
        label='Minimum number of genes in gene set',
        required=False,
        help_text="<strong>Default:</strong> 10"
    )

    maximum_size = forms.IntegerField(
        min_value=10,
        max_value=3000,
        label='Maximum number of genes in gene set',
        required=False,
        help_text="<strong>Default:</strong> 500"
    )

    sig_threshold_gsea = forms.FloatField(
        label='Significance threshold (adjusted <i>p</i>-value)',
        required=False,
        help_text="Significance threshold to filter GSEA results. <strong>Default:</strong> 0.05"
    )


class UploadFoldChangesFormGSEA(forms.Form):
    """Upload run DGE analysis or submit fold changes form class."""
    run_analysis_results_gsea = forms.BooleanField(
        label='Run differential gene expression analysis',
        required=False,
    )

    read_counts_file_gsea = forms.FileField(
        label='Upload read counts',
        required=False,
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Note:</strong> this field is optional. If you would like to perform differential gene "
                  "expression analysis and are submitting GSEA results, please ensure the counts matrix corresponds "
                  "to the class labels file above. The following file formats are supported: *.csv or *.tsv."
    )

    significance_value_run_gsea = forms.FloatField(
        label='Significance threshold (adjusted <i>p</i>-value)',
        required=False,
        help_text="Significance threshold to filter differentially expressed genes. <strong>Default:</strong> 0.05"
    )

    upload_fold_changes_gsea = forms.BooleanField(
        label='Upload differential gene expression analysis results',
        required=False,
    )

    fold_changes_file_gsea = forms.FileField(
        label='Upload fold changes',
        required=False,
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Note:</strong> this field is optional. If you are uploading GSEA results, please ensure the "
                  "fold changes correspond to the same samples from your GSEA run. The following file formats are "
                  "supported: *.csv or *.tsv."
    )

    significance_value_upload_gsea = forms.FloatField(
        label='Significance threshold (adjusted <i>p</i>-value)',
        required=False,
        help_text="Significance threshold to filter differentially expressed genes. <strong>Default:</strong> 0.05"
    )


"""Custom form submission."""


class UploadCustomDataForm(forms.Form):
    """Upload custom data files form class."""
    database = forms.CharField(label='Enter database name', validators=[validators.validate_slug])
    gene_set = forms.FileField(
        label='Gene set file',
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Supported file formats:</strong> *.csv or *.tsv."
    )


class UploadMappingsForm(forms.Form):
    """Upload custom mappings form class."""
    mapping_file = forms.FileField(
        validators=[validators.FileExtensionValidator([CSV, TSV])],
        help_text="<strong>Supported file formats:</strong> *.csv or *.tsv."
    )


class UserCreationForm(forms.ModelForm):
    """A form for creating new users. Includes all the required fields, plus a repeated password."""
    password1 = forms.CharField(label='Password', widget=forms.PasswordInput)
    password2 = forms.CharField(label='Password confirmation', widget=forms.PasswordInput)

    class Meta:
        model = User
        fields = ('email',)

    def clean_password2(self):
        # Check that the two password entries match
        password1 = self.cleaned_data.get("password1")
        password2 = self.cleaned_data.get("password2")
        if password1 and password2 and password1 != password2:
            raise forms.ValidationError("Passwords don't match")
        return password2

    def save(self, commit=True):
        # Save the provided password in hashed format
        user = super().save(commit=False)
        user.set_password(self.cleaned_data["password1"])
        if commit:
            user.save()
        return user


class UserChangeForm(forms.ModelForm):
    """A form for updating users. Includes all the fields on
    the user, but replaces the password field with admin's
    password hash display field.
    """
    password = ReadOnlyPasswordHashField()

    class Meta:
        model = User
        fields = ('email', 'password', 'is_active', 'is_admin')

    def clean_password(self):
        # Regardless of what the user provides, return the initial value.
        # This is done here, rather than on the field, because the
        # field does not have access to the initial value
        return self.initial["password"]


# TODO: Move to admin.py
class UserAdmin(BaseUserAdmin):
    """The forms to add and change user instances."""
    form = UserChangeForm
    add_form = UserCreationForm

    # The fields to be used in displaying the User model.
    # These override the definitions on the base UserAdmin
    # that reference specific fields on auth.User.
    list_display = ('email', 'is_admin')
    list_filter = ('is_admin',)
    fieldsets = (
        (None, {'fields': ('email', 'password')}),
        ('Permissions', {'fields': ('is_admin',)}),
    )
    # add_fieldsets is not a standard ModelAdmin attribute. UserAdmin
    # overrides get_fieldsets to use this attribute when creating a user.
    add_fieldsets = (
        (None, {
            'classes': ('wide',),
            'fields': ('email', 'password1', 'password2'),
        }),
    )
    search_fields = ('email',)
    ordering = ('email',)
    filter_horizontal = ()
