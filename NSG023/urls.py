# NSG023/urls.py
from django.urls import path
from . import views
app_name = 'NSG023'
urlpatterns = [
    path('', views.nsg023_page, name='NSG023_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
