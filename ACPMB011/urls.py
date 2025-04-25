# ACPMB011/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB011'
urlpatterns = [
    path('', views.acpmb011_page, name='ACPMB011_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]