# ACPMB007/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.acpmb007_page, name='ACPMB007_page'),
    path('api/gene_list/', views.get_gene_list),
    path('api/plot_gene_image/', views.plot_gene_image),
]
