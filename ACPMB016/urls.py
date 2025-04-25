# ACPMB016/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB016'
urlpatterns = [
    path('', views.acpmb016_page, name='ACPMB016_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]