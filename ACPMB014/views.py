# ACPMB014/views.py
import os, json
import scanpy as sc
import matplotlib.pyplot as plt
from django.conf import settings
from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt

def load_adata(pubid):
    file_path = os.path.join(settings.MY_ANNDATA_DIR, f"{pubid}.h5ad")
    return sc.read_h5ad(file_path)

def acpmb014_page(request):
    dataset = {
        'title': 'TB-Lung-Visum-ACPMB014',
        'species': 'Human',
        'tissue_type': 'Lung',
        'granuloma_description': 'Medium Granuloma + B-Cluster + Satellites',
        'age': 76,
        'sex': 'M',
        'bacteria': 'M. kansasii',
        'slide_number': 'V53F21-109',
        'sample_id': 'ACPMB14',
        'pubid': 'ACPMB014',
        'method': 'Visium',
        'default_img': f"images/spots/ACPMB14.jpeg"
        
    }
    return render(request, 'dataview/view.html', {'dataset': dataset})

@csrf_exempt
def get_gene_list(request):
    adata = load_adata("ACPMB014")
    genes = adata.var_names.tolist()
    return JsonResponse({"genes": genes})

@csrf_exempt
def plot_gene_image(request):
    if request.method == "POST":
        try:
            data = json.loads(request.body.decode("utf-8"))
            gene = data.get("gene")
            adata = load_adata("ACPMB014")

            if not gene or gene not in adata.var_names:
                return JsonResponse({"error": "Invalid gene"}, status=400)

            output_dir = os.path.join("frontend", "static", "generated", "ACPMB14")
            os.makedirs(output_dir, exist_ok=True)
            img_path = os.path.join(output_dir, f"{gene}.png")

            sc.settings.set_figure_params(dpi=120, frameon=False, figsize=(4, 4), facecolor="black")
            plt.ioff()

            figs = sc.pl.spatial(
                adata, color=[gene], library_id="ACPMB14",
                show=False, alpha=0.75, alpha_img=0.3, cmap="turbo", return_fig=True
            )
            
            # 判断返回的是 list 还是单个 Figure

            if isinstance(figs, list):
                fig = figs[0]
            else:
                fig = figs

            # 保存图像
            fig.savefig(img_path, dpi=150)
            plt.close(fig)

            return JsonResponse({"image_url": f"/static/generated/ACPMB14/{gene}.png"})
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=500)
    return JsonResponse({"error": "Invalid request method"}, status=405)
