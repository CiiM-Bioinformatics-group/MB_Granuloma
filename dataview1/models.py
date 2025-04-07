## 这一部分可以看一下数据怎么储存的
from django.db import models

class SpatialDataset(models.Model):
    """
    存储空间转录组数据的数据库模型
    """
    name = models.CharField(max_length=255)  # 数据集名称
    description = models.TextField()  # 数据集描述
    species = models.CharField(max_length=100, default="Unknown")  # 物种
    method = models.CharField(max_length=100, default="Visium")  # 使用的方法
    date_published = models.DateField(null=True, blank=True)  # 发表日期
    source_url = models.URLField(blank=True, null=True)  # 数据来源链接
    spatial_data = models.JSONField()  # 存储空间坐标数据 (JSON 格式)
    expression_data = models.JSONField()  # 存储基因表达数据 (JSON 格式)

    def __str__(self):
        return self.name
