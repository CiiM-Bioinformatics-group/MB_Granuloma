from django.shortcuts import render
import json
import scanpy as sc
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
import os

#  ban GUI backend，solve macOS crash problem
import matplotlib
matplotlib.use("Agg")  # use backend without GUI
import matplotlib.pyplot as plt

# read HDF5 (AnnData) from usb
# adata = sc.read_h5ad(f"/Volumes/POMMESFRITE/MTBdata/Splited_data_for_webtool/ACPMB001.h5ad")

# from django.conf import settings

# # read anndata
# static_dir = settings.STATICFILES_DIRS[0]
# file_path = os.path.join(static_dir, "MTBdata/Splited_data_for_webtool/ACPMB001.h5ad")
# adata = sc.read_h5ad(file_path)

# read anndata
from django.conf import settings
file_path = os.path.join(settings.MY_ANNDATA_DIR, 'ACPMB001.h5ad')
adata = sc.read_h5ad(file_path)


sampleID = "ACPMB1"

@csrf_exempt
def get_gene_list(request):
    """
    get Visium datas gene list
    """
    genes = adata.var_names.tolist()
    return JsonResponse({"genes": genes})


plt.rcParams.update({
    "text.color": "white",
    "axes.labelcolor": "white",
    "xtick.color": "white",
    "ytick.color": "white"
})

#setting image style
sc.settings.set_figure_params(dpi=120, frameon=False, figsize=(4, 4), facecolor="black")

@csrf_exempt
def plot_gene_image(request):
    if request.method == "POST":
        try:
            data = json.loads(request.body.decode("utf-8"))
            gene = data.get("gene")

            if not gene:
                return JsonResponse({"error": "Missing gene name"}, status=400)

            if gene not in adata.var_names:
                return JsonResponse({"error": "Gene not found"}, status=404)

            # save fig path
            output_dir = os.path.join("frontend", "static", "generated")
            os.makedirs(output_dir, exist_ok=True)
            img_path = os.path.join(output_dir, f"{gene}.png")

            # generate right fig：back matplotlib Figure and save it self
            plt.ioff()
            fig = sc.pl.spatial(
                adata,
                color=[gene],
                library_id=sampleID,
                show=False,
                alpha=0.75, 
                alpha_img=0.3,
                cmap="turbo",
                return_fig=True  
            )
            
            fig.savefig(img_path, dpi=150)  # save
            plt.close(fig)

            return JsonResponse({"image_url": f"/static/generated/{gene}.png"})
            # return JsonResponse({"image_url": f"/static/generated/spots.png"})

        except Exception as e:
            return JsonResponse({"error": str(e)}, status=500)

    return JsonResponse({"error": "Invalid request method"}, status=405)


def dataview1_page(request):
    """
    渲染前端页面 dataview1.html
    """
    return render(request, "dataview1.html")
