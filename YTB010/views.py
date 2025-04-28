# YTB010/views.py
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

def ytb010_page(request):
    dataset = {
        'title': 'TB-Lung-Visum-YTB010',
        'species': 'Human',
        'tissue_type': 'Lung',
        'granuloma_description': 'Smaller no-border Granuloma, fibrotic',
        'age': 63,
        'sex': 'M',
        'bacteria': 'M. tuberculosis complex',
        'slide_number': 'V53F21-104',
        'sample_id': 'YTB010',
        'pubid': 'YTB010',
        'method': 'Visium',
        'default_img': f"images/spots/YTB010.jpeg"
        
    }
    return render(request, 'dataview/view.html', {'dataset': dataset})

@csrf_exempt
def get_gene_list(request):
    adata = load_adata("YTB010")
    genes = adata.var_names.tolist()
    return JsonResponse({"genes": genes})

@csrf_exempt
def plot_gene_image(request):
    if request.method == "POST":
        try:
            data = json.loads(request.body.decode("utf-8"))
            gene = data.get("gene")
            adata = load_adata("YTB010")

            if not gene or gene not in adata.var_names:
                return JsonResponse({"error": "Invalid gene"}, status=400)

            output_dir = os.path.join("frontend", "static", "generated", "YTB010")
            os.makedirs(output_dir, exist_ok=True)
            img_path = os.path.join(output_dir, f"{gene}.png")

            if not os.path.exists(img_path):
                plt.ioff()
                ax = sc.pl.spatial(
                    adata,
                    color=[gene],
                    library_id='YTB010',
                    show=False,
                    alpha=0.75, 
                    alpha_img=0.3,
                    cmap="turbo",
                    return_fig=True  
                )
                fig = ax[0].get_figure() 
                fig.savefig(img_path, dpi=150)
                plt.close(fig)

            return JsonResponse({"image_url": f"/mb-granuloma/static/generated/YTB010/{gene}.png"})
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=500)
    return JsonResponse({"error": "Invalid request method"}, status=405)
