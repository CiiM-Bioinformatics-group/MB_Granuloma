# YTB006/urls.py
from django.urls import path
from . import views
app_name = 'YTB006'
urlpatterns = [
    path('', views.ytb006_page, name='YTB006_page'),
    path('gene_list/', views.get_gene_list),
    path('plot_gene_image/', views.plot_gene_image),
]
