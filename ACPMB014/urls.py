# ACPMB014/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.acpmb014_page, name='ACPMB014_page'),
    path('api/gene_list/', views.get_gene_list),
    path('api/plot_gene_image/', views.plot_gene_image),
]
