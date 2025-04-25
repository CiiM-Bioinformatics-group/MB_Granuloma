# ACPMB010/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB010'
urlpatterns = [
    path('', views.acpmb010_page, name='ACPMB010_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
