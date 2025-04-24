# import os


# apps = [
#     'ACPMB002', 'ACPMB003', 'ACPMB004', 'ACPMB005', 'ACPMB006', 'ACPMB007',
#     'ACPMB008', 'ACPMB009', 'ACPMB010', 'ACPMB011', 'ACPMB013', 'ACPMB014', 'ACPMB016',
#     'NSG005', 'NSG013', 'NSG014', 'NSG017', 'NSG019', 'NSG020_1', 'NSG020_2', 'NSG020_3',
#     'NSG023', 'NSG024', 'YTB001', 'YTB005', 'YTB006', 'YTB008', 'YTB009', 'YTB010',
#     'YTB012', 'YTB015', 'YTB019'
# ]

# for app in apps:
#     view_file = os.path.join(app, "views.py")
#     if not os.path.exists(view_file):
#         print(f"❌ Not found: {view_file}")
#         continue

#     sampleID = app.replace("ACPMB", "ACPMB").replace("NSG", "NSG").replace("YTB", "YTB")
#     title = f"TB-Lung-Visum-{app}"
#     image_path = f"/static/images/spots/{sampleID}.jpeg"

#     content = f"""\
# import os, json
# import scanpy as sc
# import matplotlib.pyplot as plt
# from django.conf import settings
# from django.shortcuts import render
# from django.http import JsonResponse
# from django.views.decorators.csrf import csrf_exempt

# def load_adata(pubid):
#     file_path = os.path.join(settings.MY_ANNDATA_DIR, f"{{pubid}}.h5ad")
#     return sc.read_h5ad(file_path)

# def {app.lower()}_page(request):
#     dataset = {{
#         'title': '{title}',
#         'species': 'Human',
#         'tissue_type': 'Lung',
#         'granuloma_description': 'See metadata',
#         'age': 0,
#         'sex': '-',
#         'bacteria': '-',
#         'slide_number': '-',
#         'sample_id': '{sampleID}',
#         'method': 'Visium',
#         'default_img': '{image_path}',
#     }}
#     return render(request, 'dataview/view.html', {{'dataset': dataset}})

# @csrf_exempt
# def get_gene_list(request):
#     adata = load_adata("{app}")
#     genes = adata.var_names.tolist()
#     return JsonResponse({{"genes": genes}})

# @csrf_exempt
# def plot_gene_image(request):
#     if request.method == "POST":
#         try:
#             data = json.loads(request.body.decode("utf-8"))
#             gene = data.get("gene")
#             sampleID = "{sampleID}"
#             adata = load_adata("{app}")

#             if not gene or gene not in adata.var_names:
#                 return JsonResponse({{"error": "Invalid gene"}}, status=400)

#             output_dir = os.path.join("frontend", "static", "generated", sampleID)
#             os.makedirs(output_dir, exist_ok=True)
#             img_path = os.path.join(output_dir, f"{{gene}}.png")

#             sc.settings.set_figure_params(dpi=120, frameon=False, figsize=(4, 4), facecolor="black")
#             plt.ioff()
#             fig = sc.pl.spatial(
#                 adata, color=[gene], library_id=sampleID,
#                 show=False, alpha=0.75, alpha_img=0.3, cmap="turbo", return_fig=True
#             )
#             fig.savefig(img_path, dpi=150)
#             plt.close(fig)

#             return JsonResponse({{"image_url": f"/static/generated/{{sampleID}}/{{gene}}.png"}})
#         except Exception as e:
#             return JsonResponse({{"error": str(e)}}, status=500)
#     return JsonResponse({{"error": "Invalid request method"}}, status=405)
# """
#     with open(view_file, "w") as f:
#         f.write(content)
#     print(f"✅ Updated: {view_file}")

import pandas as pd
import os

# Load the re-uploaded TSV file
tsv_path = "/Volumes/POMMESFRITE/MTBdata/ForPublication.tsv"
df = pd.read_csv(tsv_path, sep="\t")

# Assume app directories are in the same location as this script is executed
output_base = "/Users/boop/Desktop/MB_Granuloma"
os.makedirs(output_base, exist_ok=True)
updated_apps = []

for _, row in df.iterrows():
    pubid = row["PubID"]
    sampleID = row["sampleID"]
    app_path = os.path.join(output_base, pubid)
    os.makedirs(app_path, exist_ok=True)

    views_code = f"""\
# {pubid}/views.py
import os, json
import scanpy as sc
import matplotlib.pyplot as plt
from django.conf import settings
from django.shortcuts import render
from django.http import JsonResponse
from django.views.decorators.csrf import csrf_exempt

def load_adata(pubid):
    file_path = os.path.join(settings.MY_ANNDATA_DIR, f"{{pubid}}.h5ad")
    return sc.read_h5ad(file_path)

def {pubid.lower()}_page(request):
    dataset = {{
        'title': 'TB-{row["Tissue_Type"]}-Visum-{pubid}',
        'species': 'Human',
        'tissue_type': '{row["Tissue_Type"]}',
        'granuloma_description': '{row["Granuloma description"]}',
        'age': {int(row["Age"])},
        'sex': '{row["Sex"]}',
        'bacteria': '{row["bacteria"]}',
        'slide_number': '{row["SlideNumber"]}',
        'sample_id': '{sampleID}',
        'pubid': '{pubid}',
        'method': 'Visium',
        'default_img': f"/static/images/spots/{sampleID}.jpeg"
    }}
    return render(request, 'dataview/view.html', {{'dataset': dataset}})

@csrf_exempt
def get_gene_list(request):
    adata = load_adata("{pubid}")
    genes = adata.var_names.tolist()
    return JsonResponse({{"genes": genes}})

@csrf_exempt
def plot_gene_image(request):
    if request.method == "POST":
        try:
            data = json.loads(request.body.decode("utf-8"))
            gene = data.get("gene")
            adata = load_adata("{pubid}")

            if not gene or gene not in adata.var_names:
                return JsonResponse({{"error": "Invalid gene"}}, status=400)

            output_dir = os.path.join("frontend", "static", "generated", "{sampleID}")
            os.makedirs(output_dir, exist_ok=True)
            img_path = os.path.join(output_dir, f"{{gene}}.png")

            sc.settings.set_figure_params(dpi=120, frameon=False, figsize=(4, 4), facecolor="black")
            plt.ioff()
            fig = sc.pl.spatial(
                adata, color=[gene], library_id="{sampleID}",
                show=False, alpha=0.75, alpha_img=0.3, cmap="turbo", return_fig=True
            )
            fig.savefig(img_path, dpi=150)
            plt.close(fig)

            return JsonResponse({{"image_url": f"/static/generated/{sampleID}/{{gene}}.png"}})
        except Exception as e:
            return JsonResponse({{"error": str(e)}}, status=500)
    return JsonResponse({{"error": "Invalid request method"}}, status=405)
"""
    with open(os.path.join(app_path, "views.py"), "w") as f:
        f.write(views_code)
    updated_apps.append(pubid)

updated_apps
