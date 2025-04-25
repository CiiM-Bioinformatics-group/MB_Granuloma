# YTB010/urls.py
from django.urls import path
from . import views
app_name = 'YTB010'
urlpatterns = [
    path('', views.ytb010_page, name='YTB010_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
