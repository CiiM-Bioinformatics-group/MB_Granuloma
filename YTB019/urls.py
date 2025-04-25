# YTB019/urls.py
from django.urls import path
from . import views
app_name = 'YTB019'
urlpatterns = [
    path('', views.ytb019_page, name='YTB019_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
