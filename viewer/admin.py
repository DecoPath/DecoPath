# -*- coding: utf-8 -*-

"""Admin."""

from django.contrib import admin
from django.contrib.auth.models import Group

from viewer.forms import UserAdmin
from viewer.models import PathwayDatabase, EnrichmentResult, Pathway, User, PathwayHierarchy


class PathwayDatabaseAdmin(admin.ModelAdmin):
    list_display = (
        'database_name', 'version', 'number_of_pathways',
    )


class GseaResultsAdmin(admin.ModelAdmin):
    list_display = (
        'user', 'date', 'result_id',
    )


class PathwayAdmin(admin.ModelAdmin):
    list_display = (
        'pathway_id', 'pathway_name', 'pathway_database', 'get_equivalent_pathways', 'decopath_id', 'decopath_name',
        'number_of_equivalent_pathways',
    )
    search_fields = (
        'pathway_id', 'pathway_name', 'pathway_database', 'decopath_id', 'decopath_name',
    )


class PathwayHierarchyAdmin(admin.ModelAdmin):
    list_display = (
        'name',
    )


admin.site.register(PathwayDatabase, PathwayDatabaseAdmin)
admin.site.register(EnrichmentResult, GseaResultsAdmin)
admin.site.register(Pathway, PathwayAdmin)
admin.site.register(User, UserAdmin)
admin.site.register(PathwayHierarchy, PathwayHierarchyAdmin)
admin.site.unregister(Group)
