# NSG020_2/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.nsg020_2_page, name='NSG020_2_page'),
    path('api/gene_list/', views.get_gene_list),
    path('api/plot_gene_image/', views.plot_gene_image),
]
