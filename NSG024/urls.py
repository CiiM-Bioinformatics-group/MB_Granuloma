# NSG024/urls.py
from django.urls import path
from . import views
app_name = 'NSG024'
urlpatterns = [
    path('', views.nsg024_page, name='NSG024_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
