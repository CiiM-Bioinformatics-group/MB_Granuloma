# dataview8/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.dataview8_page, name='dataview8_page'),
    path('api/gene_list/', views.get_gene_list),
    path('api/plot_gene_image/', views.plot_gene_image),
]
