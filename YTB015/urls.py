# YTB015/urls.py
from django.urls import path
from . import views
app_name = 'YTB015'
urlpatterns = [
    path('', views.ytb015_page, name='YTB015_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
