# ACPMB003/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB003'

urlpatterns = [
    path('', views.acpmb003_page, name='ACPMB003_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
