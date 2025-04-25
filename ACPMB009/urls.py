# ACPMB009/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB009'
urlpatterns = [
    path('', views.acpmb009_page, name='ACPMB009_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]