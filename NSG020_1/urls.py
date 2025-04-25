# NSG020_1/urls.py
from django.urls import path
from . import views
app_name = 'NSG020_1'
urlpatterns = [
    path('', views.nsg020_1_page, name='NSG020_1_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
