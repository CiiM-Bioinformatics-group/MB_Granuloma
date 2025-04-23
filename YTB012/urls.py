# YTB012/urls.py
from django.urls import path
from . import views

urlpatterns = [
    path('', views.ytb012_page, name='YTB012_page'),
    path('api/gene_list/', views.get_gene_list),
    path('api/plot_gene_image/', views.plot_gene_image),
]
