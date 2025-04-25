# ACPMB014/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB014'
urlpatterns = [
    path('', views.acpmb014_page, name='ACPMB014_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]