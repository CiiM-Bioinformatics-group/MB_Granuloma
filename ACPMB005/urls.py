# ACPMB005/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB005'
urlpatterns = [
    path('', views.acpmb005_page, name='ACPMB005_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
