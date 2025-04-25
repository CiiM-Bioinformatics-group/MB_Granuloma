# NSG005/urls.py
from django.urls import path
from . import views
app_name = 'NSG005'
urlpatterns = [
    path('', views.nsg005_page, name='NSG005_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]

