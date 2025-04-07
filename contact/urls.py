# myapp/urls.py



# contact/urls.py

from django.urls import path
from . import views
from .views import ContactCreate  # 确保类名与 views.py 中的定义一致

urlpatterns = [
    path('', views.home, name='home'),  # 主页路径
    path('contact/', ContactCreate.as_view(), name='contact-create'),
    path('about/', views.about, name='about'),  # 配置 URL 路径和视图
    path('contacts/', views.about, name='contacts'),
    path('mapping/', views.about, name='mapping'),
    path('search/', views.about, name='search'),
    path('datasets/', views.about, name='datasets'),
    path('help/', views.about, name='help'),
 
]


