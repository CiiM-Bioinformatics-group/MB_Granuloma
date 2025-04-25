# ACPMB008/urls.py
from django.urls import path
from . import views

app_name = 'ACPMB008'
urlpatterns = [
    path('', views.acpmb008_page, name='ACPMB008_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
