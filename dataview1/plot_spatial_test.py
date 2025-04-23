import scanpy as sc
import matplotlib.pyplot as plt


plt.rcParams.update({
    "text.color": "white",
    "axes.labelcolor": "white",
    "xtick.color": "white",
    "ytick.color": "white"
})

# 设置 scanpy 参数
sc.settings.set_figure_params(dpi=120, frameon=False, figsize=(6, 6), facecolor="black")

# 读取数据
sampleID = "ACPMB3"
adata = sc.read_h5ad("/Volumes/POMMESFRITE/MTBdata/Splited_data_for_webtool/ACPMB003.h5ad")


sc.pl.spatial(adata, color=['Spots'], library_id=sampleID, alpha=0.75, alpha_img=0.3)