# ACPMB002/urls.py
from django.urls import path
from . import views
app_name = 'ACPMB002'
urlpatterns = [
    path('', views.acpmb002_page, name='ACPMB002_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
