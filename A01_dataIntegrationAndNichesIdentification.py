#!/usr/bin/env python
# coding: utf-8

import rapids_singlecell as rsc
import gc
import torch
import random
random.seed(1024)
import decoupler
print(decoupler.__version__)

import scanpy as sc
import cupy as cp
import time
import warnings
import anndata as an
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
warnings.filterwarnings("ignore")
print(torch.multiprocessing.cpu_count())
sc._settings.ScanpyConfig.n_jobs = -1
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(
    managed_memory=True,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)
cp.cuda.set_allocator(rmm_cupy_allocator)
sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(4, 4), facecolor="white")

# in this cell data from different samples would be combined
# Paths to your Visium datasets
df_sample = pd.read_csv("/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/sampleIndex.tsv", sep="\t")
df_sample = df_sample[df_sample.batch!="skin"]

# Paths to your Visium datasets
sample_id = df_sample[(df_sample.batch!="MTB_batch0" ) & 
                        (df_sample.quality!="bad") & 
                        (df_sample["Tissue Type"]!="Pericardium") & 
                        (df_sample["Tissue Type"]!="unclear") & 
                        (df_sample["sampleID"]!="ACPMB12") &
                        (df_sample["sampleID"]!="YTB002") & 
                        (df_sample["sampleID"]!="YTB013")].sampleID
sample_id = sample_id.to_list()
print(sample_id, "\n", len(sample_id))
data_path = "/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Outputs_fromSpaceranger/"

# Load each sample
adata_list = [sc.read_visium(data_path+sample+"/outs/") for sample in sample_id]
print(len(adata_list))
# Preprocess each sample
for adata in adata_list:
    adata.var_names_make_unique()


# Save the spatial metadata before combining
spatial_uns_metadata = {
    sample_id[i] : adata_list[i].uns['spatial'] for i in range(len(sample_id))
}

# Concatenate the datasets
adatas = adata_list[0].concatenate(adata_list[1:], batch_key='batch', batch_categories=[batch for batch in sample_id])

adatas.uns['spatial'] = {}
for sample in adatas.obs['batch'].unique():
    adatas.uns['spatial'][sample] = spatial_uns_metadata[sample][sample]

# filteration
sc.pp.filter_cells(adatas, min_counts=2000)
sc.pp.filter_cells(adatas, max_counts=100000)
sc.pp.filter_cells(adatas, min_genes=50)
sc.pp.filter_genes(adatas, min_cells=200)


# QC steps
adatas.var["mt"] = adatas.var_names.str.startswith("MT-")
adatas.var["RIBO"] = adatas.var_names.str.startswith("RPS")
sc.pp.calculate_qc_metrics(adatas, qc_vars=["mt","RIBO"], inplace=True)

fig, ax = plt.subplots(1,2, figsize=(7, 3))
sc.pl.scatter(adatas, "total_counts", "pct_counts_mt", ax=ax[0], show=False)
sc.pl.scatter(adatas, "total_counts", "n_genes_by_counts", ax=ax[1], show=False)
plt.tight_layout()
plt.show()


print(adatas.X.min(), adatas.X.max())
adatas.layers["counts"] = adatas.X.copy()
gc.collect()


adatas = adatas[adatas.obs["pct_counts_mt"] < 20]
# rsc.get.anndata_to_GPU(adatas, convert_all=True)
sc.pp.normalize_total(adatas, target_sum=1e4)
sc.pp.log1p(adatas)

adatas.raw = adatas

print(adatas.X.min(), adatas.X.max())
sc.pp.highly_variable_genes(adatas, n_top_genes=5000, batch_key="batch")
sc.pl.highly_variable_genes(adatas)

sc.pp.pca(adatas, n_comps=100)
rsc.get.anndata_to_GPU(adatas)

sc.settings.set_figure_params(dpi=100, frameon=False, figsize=(12, 4), facecolor="white")

rsc.pp.harmony_integrate(adatas, key="batch")

sc.pl.pca_variance_ratio(adatas, n_pcs=200, log=True)

rsc.pp.neighbors(adatas, n_neighbors=15, n_pcs=100, algorithm="brute")
rsc.tl.umap(adatas)
gc.collect()
torch.cuda.empty_cache()

# clustering
rsc.tl.louvain(adatas, resolution=1)
rsc.tl.leiden(adatas, resolution=1)
sc.settings.set_figure_params(dpi=180, frameon=False, figsize=(8, 7), facecolor="white")
sc.pl.umap(adatas, color=["louvain", "leiden", "batch"], legend_loc="on data")

gc.collect()
torch.cuda.empty_cache()
gc.collect()
rsc.pp.neighbors(adatas, use_rep='X_pca_harmony', n_neighbors=12, n_pcs=75, algorithm="brute")
rsc.tl.umap(adatas)
gc.collect()
torch.cuda.empty_cache()

# clustering
rsc.tl.leiden(adatas, resolution=0.6)
sc.settings.set_figure_params(dpi=180, frameon=False, figsize=(8, 7), facecolor="white")
sc.pl.umap(adatas, color=["leiden", "batch"], legend_loc="on data", palette="tab20")

sc.settings.set_figure_params(dpi=180, frameon=False, figsize=(8, 7), facecolor="white")
sc.pl.umap(adatas, color=["leiden", "SampleID"], legend_loc="on data", palette="tab20", size=2.4)
sc.pl.umap(adatas, color=["leiden", "SampleID"], legend_loc="on data", size=2.4)

adatas.obs["SampleID"] = adatas.obs['batch'].replace(dict(zip(df_sample['sampleID'], df_sample['PubID'])))

rsc.get.anndata_to_CPU(adatas, convert_all=True)
gc.collect()
torch.cuda.empty_cache()

sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(5, 5), facecolor="white")
sc.tl.embedding_density(adatas, groupby="SampleID")
sc.pl.embedding_density(adatas, groupby="SampleID", ncols=8)
plt.show()

# Assuming adatas is your integrated AnnData object with Leiden clusters
# Generate a color palette
n_clusters = adatas.obs['leiden'].unique()
palette = sns.color_palette("tab20", len(n_clusters))

# Create a dictionary mapping each cluster to a color
cluster_colors = {n_clusters[i]: palette[i] for i in range(len(n_clusters))}

# Assign colors to leiden clusters
adatas.uns['Leiden_colors'] = [cluster_colors[i] for i in n_clusters]

rsc.get.anndata_to_GPU(adatas, convert_all=True)
#rsc.tl.diffmap(adatas, n_comps=20)
rsc.tl.diffmap(adatas, n_comps=60)
sc.pl.diffmap(adatas, color="leiden", size=10, alpha=0.7, palette=cluster_colors)
rsc.tl.draw_graph(adatas)
rsc.get.anndata_to_CPU(adatas, convert_all=True)
gc.collect()
torch.cuda.empty_cache()

sc.pl.draw_graph(adatas, color="leiden", palette=cluster_colors)

fig, ax = plt.subplots(1,2, figsize=(14, 5))
sc.pl.umap(adatas, color=['leiden'], ax=ax[0], show=False, size=2, palette="tab20", alpha=0.7)
sc.pl.umap(adatas, color=['SampleID'], ax=ax[1], show=False, size=2, palette="tab20", alpha=0.7)
# Adjust layout to prevent overlap
plt.tight_layout()
plt.show()

adatas.write_h5ad('/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Anndata/AllDataHarmolized.h5ad')
gc.collect()


import decoupler as dc
import liana as li


# check markers gene of different cell types
df_markerGenes = pd.read_csv("/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Analysis/MTB_Granuloma_Marker_Genes_from_Leo.tsv", sep="\t")
markers = {key: [x for x in list(df_markerGenes[key]) if str(x) != 'nan'][:] for key in df_markerGenes.columns if key!="AT1/AT2"}

# Query Omnipath and get PanglaoDB
markers_plasmaCell = dc.get_resource('PanglaoDB')
markers_plasmaCell.organ.value_counts()
markers_plasmaCell = markers_plasmaCell[markers_plasmaCell['human'] & (markers_plasmaCell['human_sensitivity'] > 0.4) & 
                    (markers_plasmaCell["organ"].isin(["Immune system", "Lungs", "Smooth muscle"]))]

dc.run_ora(
    mat=adatas,
    net=markers_plasmaCell,
    source='cell_type',
    target='genesymbol',
    min_n=3,
    verbose=True,
    use_raw=True
)


acts = dc.get_acts(adatas, obsm_key='ora_estimate')

# We need to remove inf and set them to the maximum value observed for pvals=0
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e
df = dc.rank_sources_groups(acts, groupby='leiden', reference='rest', method='t-test_overestim_var')


n_ctypes = 3
ctypes_dict = df.groupby('group').head(n_ctypes).groupby('group')['names'].apply(lambda x: list(x)).to_dict()


sc.pl.matrixplot(acts, ctypes_dict, 'leiden', dendrogram=True, standard_scale='var',
                 colorbar_title='Z-scaled scores', cmap='RdBu_r')



for i in markers_plasmaCell[markers_plasmaCell.cell_type=="Plasma cells"].genesymbol:
    print(f'\'{i}\',', end="")


# In[46]:


markers = {'AT1': ['AGER',
  'ITLN2',
  'CAV1',
  'RTKN2',
  'EMP2',
  'SCEL',
  'CEACAM6',
  'KRT7',
  'HOPX',
  'PDPN',
  'SEMA3B',
  'SPOCK2',
  'CLDN18'],
 'AT2': ['SFTPB',
  'SFTPA1',
  'SFTPD',
  'PLA2G1B',
  'HHIP',
  'LAMP3',
  'PGC',
  'ABCA3',
  'WIF1',
  'ORM1',
  'TTN',
  'LRRK2',
  'FGG',
  'WFDC12'],
 'B': ['CD79A',
  'CD79B',
  'MS4A1',
  'FCRL4',
  'FCRLA',
  'BLK',
  'TCL1A',
  'VPREB3',
  'TNFRSF13C',
  'CD19',
  'TNFRSF13B',
  'FAM30A',
  'CR2'],
 'Basal': ['DLK2',
  'NPPC',
  'S100A2',
  'KRT5',
  'KRT6A',
  'KRT6C',
  'KRT13',
  'KRT14',
  'KRT15',
  'KRT16',
  'CALML3',
  'ADH7',
  'MMP13',
  'CLCA2',
  'SERPINB4',
  'MMP3',
  'DCAF12L2',
  'SPRR2D',
  'LGALS7'],
 'DC': ['CD207',
  'CD1A',
  'CD1C',
  'CD1E',
  'CLEC4F',
  'CLEC4C',
  'CLEC9A',
  'CLEC10A',
  'IL12B',
  'IL22RA2',
  'FCER1A',
  'PLD4',
  'KCNK10',
  'CCL22',
  'MMP12',
  'SHD',
  'S100B',
  'CCL17',
  'FCGR2B'],
 'general EC': ['PECAM1',
  'CLDN5',
  'CLEC14A',
  'PTPRB',
  'RAMP2',
  'RAMP3',
  'VWF',
  'TMEM100',
  'CALCRL',
  'LDB2',
  'A2M',
  'JAM2',
  'CYYR1',
  'CDH5',
  'SOX18',
  'S1PR1',
  'LRRC32',
  'FAM107A',
  'TIE1',
  'EGFL7',
  'F2RL3',
  'ROBO4',
  'FOXF1',
  'EMCN',
  'SLCO2A1',
  'ESAM',
  'GNG11',
  'ESM1',
  'TIMP3',
  'SPARC',
  'SPARCL1',
  'IGFBP7',
  'TSPAN7',
  'AQP1',
  'HYAL2',
  'KDR',
  'CRIP2',
  'NRN1'],
 'EC arterial': ['DKK2',
  'GJA5',
  'BMX',
  'SOX17',
  'SEMA3G',
  'CXCL12',
  'ENPP2',
  'SSTR1',
  'HEY1',
  'AIF1L',
  'EFNB2',
  'NPR3'],
 'EC lymphatic': ['CCL21',
  'MMRN1',
  'GRIN2B',
  'PROX1',
  'PHF24',
  'TM4SF18',
  'MYCN',
  'PPFIBP1',
  'KRTAP4-11',
  'TSPEAR',
  'TFF3'],
 'EC capillary': ['CA4', 'CD300LG', 'FCN3', 'SLC6A4', 'PRX', 'GPIHBP1'],
 'EC venous': ['ACKR1',
  'PLVAP',
  'SELE',
  'CCL14',
  'MEOX1',
  'ZNF385D',
  'DIPK2B',
  'ADGRL4',
  'FAM167B',
  'CRHBP',
  'EMCN',
  'TM4SF18',
  'VCAM1',
  'APLNR',
  'COL15A1',
  'CYTL1',
  'ABCC8'],
 'General Fibro': ['DCN',
  'COL1A2',
  'LUM',
  'DPT',
  'ADH1B',
  'PRELP',
  'COL6A2',
  'ABCA8',
  'LAMA2',
  'COL3A1',
  'MFAP4',
  'COL6A5',
  'SERPINF1',
  'CERCAM',
  'FBN1',
  'PDGFRA',
  'PDGFRB'],
 'Adventitial Fibro': ['SCARA5',
  'MYOC',
  'HTRA3',
  'PODN',
  'CD248',
  'PI16',
  'SFRP1',
  'SFRP2',
  'SFRP4',
  'FBLN2',
  'CXCL14',
  'MFAP5',
  'THBS2',
  'CCDC80'],
 'Alveolar Fibro': ['SCN7A',
  'NPNT',
  'GPC3',
  'ITGA8',
  'VEGFD',
  'GRIA1',
  'MUSK',
  'FGFR4',
  'CAMK2N1',
  'GDF10',
  'ADAMTS8',
  'ANGPT1',
  'BMP5',
  'IGFN1'],
 'Myofibroblast': ['CCL11',
  'ITGBL1',
  'CTHRC1',
  'COMP',
  'POSTN',
  'CDH11',
  'ITGBL1',
  'TNC',
  'COL10A1',
  'FNDC1',
  'MFAP2'],
 'Innate lymphoid cell NK': ['NKG7',
  'GNLY',
  'PRF1',
  'KLRD1',
  'FGFBP2',
  'CST7',
  'KLRF1',
  'TRDC',
  'SH2D1B',
  'CD7',
  'KLRB1',
  'CD160',
  'NCR3',
  'KIR2DL4',
  'TMIGD2',
  'TXK',
  'GZMA',
  'GZMB',
  'GZMH',
  'GZMM',
  'FCRL6',
  'CCL4',
  'CCL5',
  'CTSW',
  'CD247',
  'IFNG'],
 'Secretory': ['BPIFB1',
  'MUC4',
  'MUC5AC',
  'MUC5B',
  'MUC6',
  'MUC7',
  'MUC16',
  'MSMB',
  'WFDC2',
  'CYP2F1',
  'LCN2',
  'FAM3D',
  'TCN1',
  'TSPAN8',
  'CP'],
 'general Macrophage': ['AGRP',
  'MARCO',
  'C1QB',
  'APOC1',
  'TREM1',
  'CYBB',
  'CD68',
  'MRC1',
  'ITGAM',
  'MERTK',
  'MSR1',
  'FCGR1A',
  'LYZ'],
 'Alveolar Macrophage': ['RBP4',
  'CXCL5',
  'FABP4',
  'FABP3',
  'MME',
  'PPARG',
  'INHBA',
  'GLDN',
  'STAC',
  'TRHDE'],
 'Monocyte': ['ITGAX',
  'FCGR3A',
  'APOBEC3A',
  'CCR2',
  'FCAR',
  'S100A8',
  'S100A12',
  'FFAR3',
  'CD300E',
  'PADI4',
  'RNASE2',
  'LILRA5',
  'LILRB2',
  'PROK2',
  'CFP',
  'VCAN',
  'S100A9',
  'CD14'],
 'Mast cell': ['CPA3',
  'MS4A2',
  'SLC18A2',
  'TPSG1',
  'GATA1',
  'GCSAML',
  'HDC',
  'CDK15',
  'MS4A3',
  'CALB2'],
 'Pericyte': ['COX4I2',
  'FAM162B',
  'AGTR1',
  'KCNK3',
  'TBX5',
  'HIGD1B',
  'LAMC3',
  'CDH19',
  'MATN3',
  'KCNK17',
  'RGN',
  'LRCOL1',
  'NDUFA4L2',
  'PDGFRB',
  'FHL5',
  'LGI4',
  'BGN',
  'GJA4',
  'F10',
  'CACNA1H'],
 'SMC': ['DES',
  'CNN1',
  'MYH11',
  'ACTG2',
  'ACTA2',
  'RERGL',
  'MYLK',
  'ABRA',
  'TPM2',
  'MYL9',
  'PPP1R14A'],
 'Mesothelial': ['SLN',
  'PRG4',
  'BNC1',
  'GALNT9',
  'CPA4',
  'ITLN1',
  'RSPO1',
  'BCHE',
  'TFPI2',
  'HP',
  'ANXA8',
  'HAS1',
  'PLA2G2A',
  'CALB2',
  'CPXM1'],
 'T-Cells': ['IL17A',
  'CD2',
  'CD3D',
  'CD3E',
  'CD8A',
  'CD40LG',
  'CXCR6',
  'CCR4',
  'CCR8',
  'TRAC',
  'TRAT1',
  'CD28',
  'CTLA4'], 
 'Plasma Cells': ['CD79A','SDC1','CD27','JCHAIN','MZB1','TNFRSF17','EAF2','IGKC']        }



# Define var_group_labels that match `groupby` categories
sc.settings.set_figure_params(dpi=120, frameon=False, figsize=(32, 18), facecolor="white")
var_group_labels = adatas.obs['leiden'].unique().tolist()  # Or any other logic to get matching labels
sc.pl.heatmap(adatas, markers, groupby="leiden", cmap="turbo", dendrogram=True, var_group_labels=var_group_labels, use_raw=True)
#sc.pl.dotplot(adatas, marker_genes, groupby="Region", standard_scale="var", cmap="Blues", dendrogram=True)


# In[199]:


adatas.obs['Spots'] = adatas.obs['leiden'].replace({'0':'Layer1 (Macro & DC)',
                                                    '1':'Hub0',
                                                    '2':'Layer0 (Macro & DC & NK & T)',
                                                    '3':'Layer2 (Plasma cells & Plasmacytoid DC)',
                                                    '4': "Necrotizing", 
                                                    '5': "SMC & Myo", 
                                                    '6':'Hub0',
                                                    '7':'T & B',
                                                    '8':'rare0',
                                                    '9':'Secretion & AT',
                                                    '10':'Hub0',
                                                    '11':'AT & Plasma Cells',
                                                    '12':'rare1',
                                                    '13':'rare0',
                                                    '14':'rare1',
                                                    '16':'rare1',
                                                    '17':'T & Macro',
                                                    '18':'rare2',
                                                    '19':'rare0',
                                                    '20': "SMC & Myo", 
                                                    '21':'rare0',
                                                    '22':'rare0',
                                                 })




sc.pl.umap(adatas, color='leiden',  show=False, palette=cluster_colors, size=18, title=f'Umap {sampleID}', alpha=0.75, legend_fontsize=10, legend_loc=False)

gc.collect()
sc.settings.set_figure_params(dpi=120, frameon=True, figsize=(10,10), facecolor="white")
sample = "H2315532"
sample_data = adatas[adatas.obs['batch'] == sample]
sampleID = df_sample[df_sample.sampleID==sample]['PubID'].values[0]

# Plot spatial distribution of Leiden clusters
fig, ax = plt.subplots(1,2, figsize=(8, 3.6))
sc.pl.umap(sample_data, color='Spots', ax=ax[0], show=False, palette="tab20", size=18, title=f'Umap {sampleID}', alpha=0.75, legend_fontsize=8, legend_loc=False)
sc.pl.spatial(sample_data, color='Spots', library_id=sample, title=f'Spatial Niches {sampleID}', alpha_img=0.2, spot_size=22, ax=ax[1], show=False, palette="tab20", alpha=0.75, legend_fontsize=8, components='all')
plt.tight_layout()


for i in adatas.obs["Spots"].values.unique():
    print(i)


adatas.obs['Spots'] = adatas.obs['Spots'].replace({'Layer1 (Plasma Cells & Macro & T & Fibro)': 'Middle Niche',
                                                   'Layer0 (Macro & Mono & NK & T)': 'Inner Niche',
                                                    'Layer2 (Fibro & Plasma Cells)':'Outer Niche'
                                                 })




adatas = adatas[adatas.obs["Spots"].isin(["rare0", "rare1", "rare2"]) == False]
for i in adatas.obs["Spots"].values.unique():
    print(i)


n_clusters = adatas.obs['Spots'].unique()
print(len(n_clusters))
palette = sns.color_palette("tab20", 10)

# Create a dictionary mapping each cluster to a color
Spot_colors = {n_clusters[i]: palette[i] for i in range(len(n_clusters))}
adatas.uns['Spot_colors'] = [Spot_colors[i] for i in n_clusters]  




sc.settings.set_figure_params(dpi=120, frameon=True, figsize=(12,12), facecolor="white")
for sample in adatas.obs['batch'].unique():
    print(sample)
    gc.collect()
    sample_data = adatas[adatas.obs['batch'] == sample]

    # Plot spatial distribution of Leiden clusters
    fig, ax = plt.subplots(1,1, figsize=(6, 3.5))
    # sc.pl.umap(sample_data, color='Spots', ax=ax[0], show=False, 
    #            palette=Spot_colors, size=18, title=f'Umap {sample}', 
    #            alpha=0.75, legend_fontsize=8, legend_loc=False)
    sc.pl.spatial(sample_data, color='Spots', library_id=sample, 
                  title=f'Spatial Niches', 
                  alpha_img=0.2, spot_size=22, ax=ax, show=False, 
                  palette=Spot_colors, alpha=0.75, 
                  legend_fontsize=8, components='all')
    plt.tight_layout()
    plt.savefig(f'Figures/{sample}.jpeg', dpi=300, bbox_inches='tight')
    plt.show()


adatas.write_h5ad('/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Anndata/AllDataHarmolized_01.h5ad')

# Define var_group_labels that match `groupby` categories
var_group_labels = adatas.obs['Spots'].unique().tolist()  # Or any other logic to get matching labels

sc.tl.dendrogram(adatas, groupby="Spots")
sc.pl.dotplot(adatas, markers, groupby="Spots", standard_scale="var", cmap="PiYG_r", 
              dendrogram=True, var_group_labels=var_group_labels)
sc.pl.dotplot(adatas, marker_genes, groupby="Region", standard_scale="var", cmap="Blues", dendrogram=True)

sc.tl.rank_genes_groups(adatas, "Spots", method="t-test", n_jobs=40)
sc.pl.rank_genes_groups(adatas, n_genes=25, sharey=False)


# Extract the results from the rank_genes_groups
region_names = adatas.uns['rank_genes_groups']['names']
logfoldchanges = adatas.uns['rank_genes_groups']['logfoldchanges']
pvals = adatas.uns['rank_genes_groups']['pvals']
pvals_adj = adatas.uns['rank_genes_groups']['pvals_adj']
# Convert the results into DataFrames
region_names_df = pd.DataFrame(region_names, columns=adatas.uns['rank_genes_groups']['names'].dtype.names)
logfoldchanges_df = pd.DataFrame(logfoldchanges, columns=adatas.uns['rank_genes_groups']['logfoldchanges'].dtype.names)
pvals_df = pd.DataFrame(pvals, columns=adatas.uns['rank_genes_groups']['pvals'].dtype.names)
pvals_adj_df = pd.DataFrame(pvals_adj, columns=adatas.uns['rank_genes_groups']['pvals_adj'].dtype.names)


# The regions (groups) are in the columns of the DataFrames
regions = region_names_df.columns
# Flatten the data and combine
results_df = pd.DataFrame({
    'Region': list(regions) * (region_names_df.shape[0]),
    'Gene': region_names_df.values.flatten(),
    'LogFoldChange': logfoldchanges_df.values.flatten(),
    'PValue': pvals_df.values.flatten(),
    'AdjustedPValue': pvals_adj_df.values.flatten(),
})


# Save the results to a CSV file
results_df.to_csv("ranked_genes_by_region.csv", index=False)

sc.pl.rank_genes_groups_heatmap(adatas, n_genes=5, standard_scale='var', figsize=(16, 8), cmap="PiYG_r", show_gene_labels=True)

# Outside Group Expression: To ensure specificity, you want to limit how many cells outside this group express the same gene.
# Threshold Setting: By setting min_in_group_fraction to a value, you are specifying the minimum fraction of cells in the group that must express the gene for it to be considered as a potential marker.
# Filters out genes based on log fold change and fraction of genes expressing the gene within and outside the groupby categories.
sc.tl.filter_rank_genes_groups(adatas,
                               min_in_group_fraction=0.5,
                               max_out_group_fraction=0.30,
                               min_fold_change=0.49, #1.4 times higher
                               compare_abs=False) #log2 fold change


sc.pl.rank_genes_groups_stacked_violin(adatas, n_genes=5, key='rank_genes_groups_filtered', dendrogram=False, standard_scale='var')



sc.pl.rank_genes_groups_heatmap(adatas, n_genes=5, standard_scale='var', 
                                figsize=(16, 12), cmap="PiYG_r", dendrogram=False,
                                show_gene_labels=True, key='rank_genes_groups_filtered')








