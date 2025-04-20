"""
URL configuration for myproject project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""


from django.contrib import admin
from django.urls import path, include
from django.views.generic import TemplateView
from django.conf import settings
from django.conf.urls.static import static
from dataview2.views import dataview2_page
from dataview3.views import dataview3_page
from dataview4.views import dataview4_page
from dataview5.views import dataview5_page
from dataview6.views import dataview6_page
from dataview7.views import dataview7_page
from dataview8.views import dataview8_page
from dataview9.views import dataview9_page
from dataview10.views import dataview10_page

urlpatterns = [
    path('admin/', admin.site.urls),
    path('api/', include('contact.urls')),  # API 相关的路径
    path('', TemplateView.as_view(template_name="home.html"), name="home"),
    path('datasets/', TemplateView.as_view(template_name="datasets.html"), name="datasets"),
    path('about/', TemplateView.as_view(template_name="about.html"), name="about"),
    path('contacts/', TemplateView.as_view(template_name="contacts.html"), name="contacts"),
    path('lungs/', TemplateView.as_view(template_name="lungs.html"), name="lungs"),
    path('lymph/', TemplateView.as_view(template_name="lymph.html"), name="lymph"),
    path('sample/', TemplateView.as_view(template_name="sample.html"), name="sample"),

    path('dataview1/', TemplateView.as_view(template_name="dataview1.html"), name="dataview1"),
    path('dataview1/api/', include('dataview1.urls')),  
    

    #add sample pages
    path('dataview2/', dataview2_page, name='dataview2'),
    path('dataview3/', dataview2_page, name='dataview3'),
    path('dataview4/', dataview2_page, name='dataview4'),
    path('dataview5/', dataview2_page, name='dataview5'),
    path('dataview6/', dataview2_page, name='dataview6'),
    path('dataview7/', dataview2_page, name='dataview7'),
    path('dataview8/', dataview2_page, name='dataview8'),
    path('dataview9/', dataview2_page, name='dataview9'),
    path('dataview10/', dataview2_page, name='dataview10'),

    #add more apps.
    path("dataview2/api/", include("dataview2.urls")),
    path("dataview3/api/", include("dataview3.urls")),
    path("dataview4/api/", include("dataview3.urls")),
    path("dataview5/api/", include("dataview3.urls")),
    path("dataview6/api/", include("dataview3.urls")),
    path("dataview7/api/", include("dataview3.urls")),
    path("dataview8/api/", include("dataview3.urls")),
    path("dataview9/api/", include("dataview3.urls")),
    path("dataview10/api/", include("dataview10.urls")),
]


# **添加静态文件路径**
if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATICFILES_DIRS[0])

##notes
#dataviewX/ | open page | views.py → dataviewX_page()
#dataviewX/api/gene_list/ | JS fetch get  gene list | urls.py → get_gene_list
#dataviewX/api/plot_gene_image/ | JS fetch plot | urls.py → plot_gene_image
