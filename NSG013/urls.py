# NSG013/urls.py
from django.urls import path
from . import views
app_name = 'NSG013'
urlpatterns = [
    path('', views.nsg013_page, name='NSG013_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
