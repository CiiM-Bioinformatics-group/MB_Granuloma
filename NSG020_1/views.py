# NSG020_1/views.py
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

def nsg020_1_page(request):
    dataset = {
        'title': 'TB-Lung-Visum-NSG020_1',
        'species': 'Human',
        'tissue_type': 'Lung',
        'granuloma_description': 'Large Granuloma + B-cell cluster',
        'age': 71,
        'sex': 'F',
        'bacteria': 'M. avium',
        'slide_number': 'V53F21-111',
        'sample_id': 'H23_25536',
        'pubid': 'NSG020_1',
        'method': 'Visium',
        'default_img': f"images/spots/H23_25536.jpeg"
        
    }
    return render(request, 'dataview/view.html', {'dataset': dataset})

@csrf_exempt
def get_gene_list(request):
    adata = load_adata("NSG020_1")
    genes = adata.var_names.tolist()
    return JsonResponse({"genes": genes})

@csrf_exempt
def plot_gene_image(request):
    if request.method == "POST":
        try:
            data = json.loads(request.body.decode("utf-8"))
            gene = data.get("gene")
            adata = load_adata("NSG020_1")

            if not gene or gene not in adata.var_names:
                return JsonResponse({"error": "Invalid gene"}, status=400)

            output_dir = os.path.join("frontend", "static", "generated", "H23_25536")
            os.makedirs(output_dir, exist_ok=True)
            img_path = os.path.join(output_dir, f"{gene}.png")

            sc.settings.set_figure_params(dpi=120, frameon=False, figsize=(4, 4), facecolor="black")
            plt.ioff()

            figs = sc.pl.spatial(
                adata, color=[gene], library_id="H23_25536",
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

            return JsonResponse({"image_url": f"/static/generated/H23_25536/{gene}.png"})
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=500)
    return JsonResponse({"error": "Invalid request method"}, status=405)
