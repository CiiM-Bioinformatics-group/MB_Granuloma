# NSG017/views.py
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

def nsg017_page(request):
    dataset = {
        'title': 'TB-Lymphnode-Visum-NSG017',
        'species': 'Human',
        'tissue_type': 'Lymphnode',
        'granuloma_description': 'Large Granuloma',
        'age': 64,
        'sex': 'M',
        'bacteria': 'M. tuberculosis',
        'slide_number': 'V53F21-121',
        'sample_id': 'H2213396',
        'pubid': 'NSG017',
        'method': 'Visium',
        'default_img': f"/static/images/spots/H2213396.jpeg"
    }
    return render(request, 'dataview/view.html', {'dataset': dataset})

@csrf_exempt
def get_gene_list(request):
    adata = load_adata("NSG017")
    genes = adata.var_names.tolist()
    return JsonResponse({"genes": genes})

@csrf_exempt
def plot_gene_image(request):
    if request.method == "POST":
        try:
            data = json.loads(request.body.decode("utf-8"))
            gene = data.get("gene")
            adata = load_adata("NSG017")

            if not gene or gene not in adata.var_names:
                return JsonResponse({"error": "Invalid gene"}, status=400)

            output_dir = os.path.join("frontend", "static", "generated", "H2213396")
            os.makedirs(output_dir, exist_ok=True)
            img_path = os.path.join(output_dir, f"{gene}.png")

            sc.settings.set_figure_params(dpi=120, frameon=False, figsize=(4, 4), facecolor="black")
            plt.ioff()
            fig = sc.pl.spatial(
                adata, color=[gene], library_id="H2213396",
                show=False, alpha=0.75, alpha_img=0.3, cmap="turbo", return_fig=True
            )
            fig.savefig(img_path, dpi=150)
            plt.close(fig)

            return JsonResponse({"image_url": f"/static/generated/H2213396/{gene}.png"})
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=500)
    return JsonResponse({"error": "Invalid request method"}, status=405)
