

from django.urls import path
from .views import dataview1_page, get_gene_list, plot_gene_image

urlpatterns = [
    path("page/", dataview1_page, name="dataview1_page"),  # 渲染 HTML 页面
    path("gene_list/", get_gene_list, name="gene_list"),  # ✅ 修正路径
    # path("gene_expression/", get_gene_expression, name="gene_expression"),
    path("plot_gene_image/", plot_gene_image, name="plot_gene_image"),
    
]
