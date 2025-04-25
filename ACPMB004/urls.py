# ACPMB004/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB004'

urlpatterns = [
    path('', views.acpmb004_page, name='ACPMB004_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
