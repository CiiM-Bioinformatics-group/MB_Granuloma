# ACPMB005/views.py
import os, json
import scanpy as sc
import matplotlib.pyplot as plt
from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt
from .preload import adata, sampleID

def acpmb005_page(request):
    dataset = {
        'title': 'TB-Lung-Visum-ACPMB005',
        'species': 'Human',
        'tissue_type': 'Lung',
        'granuloma_description': 'Large Granuloma + B-cell cluster',
        'age': 76,
        'sex': 'M',
        'bacteria': 'M. avium',
        'slide_number': 'V53F21-076',
        'sample_id': sampleID,
        'pubid': 'ACPMB005',
        'method': 'Visium',
        'default_img': f'/static/generated/{sampleID}/defaultspot_plot.png',
    }
    return render(request, 'dataview/view.html', {'dataset': dataset})

@csrf_exempt
def get_gene_list(request):
    genes = adata.var_names.tolist()
    return JsonResponse({"genes": genes})

@csrf_exempt
def plot_gene_image(request):
    if request.method == "POST":
        try:
            data = json.loads(request.body.decode("utf-8"))
            gene = data.get("gene")
            if not gene or gene not in adata.var_names:
                return JsonResponse({"error": "Invalid gene"}, status=400)

            output_dir = os.path.join("frontend", "static", "generated", sampleID)
            os.makedirs(output_dir, exist_ok=True)
            img_path = os.path.join(output_dir, f"{gene}.png")

            sc.settings.set_figure_params(dpi=120, frameon=False, figsize=(4, 4), facecolor="black")
            plt.ioff()
            fig = sc.pl.spatial(
                adata, color=[gene], library_id=sampleID,
                show=False, alpha=0.75, alpha_img=0.3, cmap="turbo", return_fig=True
            )
            fig.savefig(img_path, dpi=150)
            plt.close(fig)

            return JsonResponse({"image_url": f"/static/generated/{sampleID}/{gene}.png"})
        except Exception as e:
            return JsonResponse({"error": str(e)}, status=500)
    return JsonResponse({"error": "Invalid request method"}, status=405)
