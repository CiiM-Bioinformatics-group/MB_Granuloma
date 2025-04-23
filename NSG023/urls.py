# NSG023/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.nsg023_page, name='NSG023_page'),
    path('api/gene_list/', views.get_gene_list),
    path('api/plot_gene_image/', views.plot_gene_image),
]
