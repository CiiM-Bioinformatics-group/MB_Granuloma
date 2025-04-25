# YTB001/urls.py
from django.urls import path
from . import views
app_name = 'YTB001'
urlpatterns = [
    path('', views.ytb001_page, name='YTB001_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
