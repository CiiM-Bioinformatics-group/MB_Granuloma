# NSG014/urls.py
from django.urls import path
from . import views
app_name = 'NSG014'
urlpatterns = [
    path('', views.nsg014_page, name='NSG014_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
