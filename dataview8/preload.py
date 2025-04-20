# dataview8/preload.py
import scanpy as sc
import os
from django.conf import settings

file_path = os.path.join(settings.MY_ANNDATA_DIR, 'ACPMB008.h5ad')  # 请根据实际修改文件名
adata = sc.read_h5ad(file_path)
sampleID = "ACPMB8"
