# NSG019/preload.py
import scanpy as sc
import os
from django.conf import settings

file_path = os.path.join(settings.MY_ANNDATA_DIR, 'NSG019.h5ad')
adata = sc.read_h5ad(file_path)
sampleID = "H2311362"
