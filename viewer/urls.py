# -*- coding: utf-8 -*-

"""Urls for the app."""

from django.urls import path

from viewer.views import *

urlpatterns = [
    # Static pages
    # Docs
    path('user_guide', user_guide, name='user_guide'),
    path('overview', overview, name='overview'),
    path('faqs', faqs, name='faqs'),
    path('imprint', imprint, name='imprint'),
    path('data_protection', data_protection, name='data_protection'),

    # Account pages
    path('login', login_request, name='login'),
    path('logout', logout_request, name='logout'),
    path('register', register, name='register'),
    path('account', account, name='account'),

    # Experiment pages
    path('run_decopath', run_decopath, name='run_decopath'),
    path('experiments', experiments, name='experiments'),
    path('results/<int:result_id>', results, name='results'),
    path('gsea_results/<int:result_id>', gsea_results, name='gsea_results'),
    path('consensus_gsea/<int:result_id>', consensus_gsea, name='consensus_gsea'),
    path('consensus_ora/<int:result_id>', consensus_ora, name='consensus_ora'),
    path('dc_genesets/<str:pathway_id>', dc_genesets, name='dc_genesets'),

    # Visualization views
    path('viz/circles/<int:result_id>', circles_viz, name='circlesviz'),
    path('viz/zoom_in/<int:result_id>/<str:pathway_id>', zoom_in, name='zoom_in'),

    # Custom pages
    path('custom_databases', custom_databases, name='custom_databases'),
    path('mappings', mappings, name='mappings'),
]
