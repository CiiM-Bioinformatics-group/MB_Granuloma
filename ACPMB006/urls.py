# ACPMB006/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB006'
urlpatterns = [
    path('', views.acpmb006_page, name='ACPMB006_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
