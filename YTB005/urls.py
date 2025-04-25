# YTB005/urls.py
from django.urls import path
from . import views
app_name = 'YTB005'
urlpatterns = [
    path('', views.ytb005_page, name='YTB005_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
