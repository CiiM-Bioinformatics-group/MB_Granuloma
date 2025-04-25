# NSG019/urls.py
from django.urls import path
from . import views
app_name = 'NSG019'
urlpatterns = [
    path('', views.nsg019_page, name='NSG019_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
