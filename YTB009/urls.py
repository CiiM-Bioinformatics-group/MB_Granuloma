# YTB009/urls.py
from django.urls import path
from . import views
app_name = 'YTB009'
urlpatterns = [
    path('', views.ytb009_page, name='YTB009_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
