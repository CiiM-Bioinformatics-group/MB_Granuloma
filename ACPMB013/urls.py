# ACPMB013/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB013'
urlpatterns = [
    path('', views.acpmb013_page, name='ACPMB013_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]