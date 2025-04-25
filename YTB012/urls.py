# YTB012/urls.py
from django.urls import path
from . import views
app_name = 'YTB012'
urlpatterns = [
    path('', views.ytb012_page, name='YTB012_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
