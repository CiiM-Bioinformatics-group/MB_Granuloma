# YTB008/urls.py
from django.urls import path
from . import views
app_name = 'YTB008'
urlpatterns = [
    path('', views.ytb008_page, name='YTB008_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
