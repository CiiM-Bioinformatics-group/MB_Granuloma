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
from django.urls import path, include


from ACPMB002.views import acpmb002_page
from ACPMB003.views import acpmb003_page
from ACPMB004.views import acpmb004_page
from ACPMB005.views import acpmb005_page
from ACPMB006.views import acpmb006_page
from ACPMB007.views import acpmb007_page
from ACPMB008.views import acpmb008_page
from ACPMB009.views import acpmb009_page
from ACPMB010.views import acpmb010_page
from ACPMB011.views import acpmb011_page
from ACPMB013.views import acpmb013_page
from ACPMB014.views import acpmb014_page
from ACPMB016.views import acpmb016_page
from NSG005.views import nsg005_page
from NSG013.views import nsg013_page
from NSG014.views import nsg014_page
from NSG017.views import nsg017_page
from NSG019.views import nsg019_page
from NSG020_1.views import nsg020_1_page
from NSG020_2.views import nsg020_2_page
from NSG020_3.views import nsg020_3_page
from NSG023.views import nsg023_page
from NSG024.views import nsg024_page
from YTB001.views import ytb001_page
from YTB005.views import ytb005_page
from YTB006.views import ytb006_page
from YTB008.views import ytb008_page
from YTB009.views import ytb009_page
from YTB010.views import ytb010_page
from YTB012.views import ytb012_page
from YTB015.views import ytb015_page
from YTB019.views import ytb019_page


urlpatterns = [
    path('admin/', admin.site.urls),
    path('api/', include('contact.urls')),  # API 相关的路径
    path('', TemplateView.as_view(template_name="home.html"), name="home"),
    path('datasets/', TemplateView.as_view(template_name="datasets.html"), name="datasets"),
    path('about/', TemplateView.as_view(template_name="about.html"), name="about"),
    path('contacts/', TemplateView.as_view(template_name="contacts.html"), name="contacts"),
    path('lungs/', TemplateView.as_view(template_name="lungs.html"), name="lungs"),
    path('lymph/', TemplateView.as_view(template_name="lymph.html"), name="lymph"),
    path('MTB/', TemplateView.as_view(template_name="MTB.html"), name="MTB"),
    path('NTM/', TemplateView.as_view(template_name="NTM.html"), name="NTM"),
    path('sample/', TemplateView.as_view(template_name="sample.html"), name="sample"),

    path('dataview1/', TemplateView.as_view(template_name="dataview1.html"), name="dataview1"),
    path('dataview1/api/', include('dataview1.urls')),  

    

    path('ACPMB002/', acpmb002_page, name='ACPMB002'),
    path('ACPMB002/api/', include('ACPMB002.urls')),
    path('ACPMB003/', acpmb003_page, name='ACPMB003'),
    path('ACPMB003/api/', include('ACPMB003.urls')),
    path('ACPMB004/', acpmb004_page, name='ACPMB004'),
    path('ACPMB004/api/', include('ACPMB004.urls')),
    path('ACPMB005/', acpmb005_page, name='ACPMB005'),
    path('ACPMB005/api/', include('ACPMB005.urls')),
    path('ACPMB006/', acpmb006_page, name='ACPMB006'),
    path('ACPMB006/api/', include('ACPMB006.urls')),
    path('ACPMB007/', acpmb007_page, name='ACPMB007'),
    path('ACPMB007/api/', include('ACPMB007.urls')),
    path('ACPMB008/', acpmb008_page, name='ACPMB008'),
    path('ACPMB008/api/', include('ACPMB008.urls')),
    path('ACPMB009/', acpmb009_page, name='ACPMB009'),
    path('ACPMB009/api/', include('ACPMB009.urls')),
    path('ACPMB010/', acpmb010_page, name='ACPMB010'),
    path('ACPMB010/api/', include('ACPMB010.urls')),
    path('ACPMB011/', acpmb011_page, name='ACPMB011'),
    path('ACPMB011/api/', include('ACPMB011.urls')),
    path('ACPMB013/', acpmb013_page, name='ACPMB013'),
    path('ACPMB013/api/', include('ACPMB013.urls')),
    path('ACPMB014/', acpmb014_page, name='ACPMB014'),
    path('ACPMB014/api/', include('ACPMB014.urls')),
    path('ACPMB016/', acpmb016_page, name='ACPMB016'),
    path('ACPMB016/api/', include('ACPMB016.urls')),
    path('NSG005/', nsg005_page, name='NSG005'),
    path('NSG005/api/', include('NSG005.urls')),
    path('NSG013/', nsg013_page, name='NSG013'),
    path('NSG013/api/', include('NSG013.urls')),
    path('NSG014/', nsg014_page, name='NSG014'),
    path('NSG014/api/', include('NSG014.urls')),
    path('NSG017/', nsg017_page, name='NSG017'),
    path('NSG017/api/', include('NSG017.urls')),
    path('NSG019/', nsg019_page, name='NSG019'),
    path('NSG019/api/', include('NSG019.urls')),
    path('NSG020_1/', nsg020_1_page, name='NSG020_1'),
    path('NSG020_1/api/', include('NSG020_1.urls')),
    path('NSG020_2/', nsg020_2_page, name='NSG020_2'),
    path('NSG020_2/api/', include('NSG020_2.urls')),
    path('NSG020_3/', nsg020_3_page, name='NSG020_3'),
    path('NSG020_3/api/', include('NSG020_3.urls')),
    path('NSG023/', nsg023_page, name='NSG023'),
    path('NSG023/api/', include('NSG023.urls')),
    path('NSG024/', nsg024_page, name='NSG024'),
    path('NSG024/api/', include('NSG024.urls')),
    path('YTB001/', ytb001_page, name='YTB001'),
    path('YTB001/api/', include('YTB001.urls')),
    path('YTB005/', ytb005_page, name='YTB005'),
    path('YTB005/api/', include('YTB005.urls')),
    path('YTB006/', ytb006_page, name='YTB006'),
    path('YTB006/api/', include('YTB006.urls')),
    path('YTB008/', ytb008_page, name='YTB008'),
    path('YTB008/api/', include('YTB008.urls')),
    path('YTB009/', ytb009_page, name='YTB009'),
    path('YTB009/api/', include('YTB009.urls')),
    path('YTB010/', ytb010_page, name='YTB010'),
    path('YTB010/api/', include('YTB010.urls')),
    path('YTB012/', ytb012_page, name='YTB012'),
    path('YTB012/api/', include('YTB012.urls')),
    path('YTB015/', ytb015_page, name='YTB015'),
    path('YTB015/api/', include('YTB015.urls')),
    path('YTB019/', ytb019_page, name='YTB019'),
    path('YTB019/api/', include('YTB019.urls')),
]





# **添加静态文件路径**
if settings.DEBUG:
    urlpatterns += static(settings.STATIC_URL, document_root=settings.STATICFILES_DIRS[0])

##notes
#dataviewX/ | open page | views.py → dataviewX_page()
#dataviewX/api/gene_list/ | JS fetch get  gene list | urls.py → get_gene_list
#dataviewX/api/plot_gene_image/ | JS fetch plot | urls.py → plot_gene_image
