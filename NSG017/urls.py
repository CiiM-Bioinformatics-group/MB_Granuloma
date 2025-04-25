# NSG017/urls.py
from django.urls import path
from . import views
app_name = 'NSG017'
urlpatterns = [
    path('', views.nsg017_page, name='NSG017_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
