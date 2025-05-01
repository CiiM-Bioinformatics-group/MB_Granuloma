#!/usr/bin/env python
# coding: utf-8
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import gc
import anndata as ad
import torch
import scipy
import pandas as pd
sc._settings.ScanpyConfig.n_jobs = -1
import cell2location
from matplotlib import rcParams
import subprocess as sb
from cell2location.utils.filtering import filter_genes
from cell2location.models import RegressionModel
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text for PDFs

def read_and_qc(sample_name, file, path):
    """
    Read one Visium file and add minimum metadata and QC metrics to adata.obs
    NOTE: var_names is ENSEMBL ID as it should be, you can always plot with sc.pl.scatter(gene_symbols='SYMBOL')
    """    
    adata = sc.read_visium(path,
                           count_file='filtered_feature_bc_matrix.h5',
                           load_images=True)
    # print(adata.uns["spatial"])
    adata.obs['sample'] = sample_name
    adata.var['SYMBOL'] = adata.var_names
    adata.var.rename(columns={'gene_ids': 'ENSEMBL'}, inplace=True)
    adata.var_names = adata.var['ENSEMBL']
    adata.var.drop(columns='ENSEMBL', inplace=True)
    
    # just in case there are non-unique ENSEMBL IDs
    adata.var_names_make_unique()

    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    adata.var['mt'] = [gene.startswith('mt-') for gene in adata.var['SYMBOL']]
    adata.obs['mt_frac'] = adata[:, adata.var['mt'].tolist()].X.sum(1).A.squeeze()/adata.obs['total_counts']
    
    # add sample name to obs names
    adata.obs["sample"] = [str(i) for i in adata.obs['sample']]
    adata.obs_names = 's' + adata.obs["sample"] \
                          + '_' + adata.obs_names
    adata.obs.index.name = 'spot_id'
    
    file = list(adata.uns['spatial'].keys())[0]
    adata.uns['spatial'][sample_name] = adata.uns['spatial'][file].copy()
    #del adata.uns['spatial'][file]
    print(adata.uns['spatial'].keys())
    print(file)
    
    return adata

samples = pd.read_csv("/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/sampleIndex.txt", sep=",")
samples = samples[(samples.quality!="bad") & (samples.batch!="skin")]


samples = samples.sampleID.values
print(torch.cuda.is_available())

for sample in samples:
    results_folder = '/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/cell2location/'
    # create paths and names to results folders for reference regression and cell2location models
    run_name = f'{results_folder}/{sample}_cell2location_map'
    sb.run(f'mkdir {run_name}', shell=True)
    path = "/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Outputs_fromSpaceranger/"+sample+"/outs/"
    adata_vis = read_and_qc(sample, sample, path)


    adata_ref = sc.read_h5ad("/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Analysis/results_Cell2Loc/Granuloma_analysis/reference_signatures/sc.h5ad")
    mod = cell2location.models.RegressionModel.load("/vol/projects/BIIM/Spatial_Transcriptomics_Projects/MTB_Granuloma/Analysis/results_Cell2Loc/Granuloma_analysis/reference_signatures/", 
                                                adata_ref)


    # export estimated expression in each cluster
    if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
        inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    else:
        inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                        for i in adata_ref.uns['mod']['factor_names']]].copy()
    inf_aver.columns = adata_ref.uns['mod']['factor_names']

    # find shared genes and subset both anndata and reference signatures
    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

    # create and train the model
    mod = cell2location.models.Cell2location(
        adata_vis, cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=20,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20
    )
    mod.view_anndata_setup()


    mod.train(max_epochs=30000,
              # train using full data (batch_size=None)
              batch_size=None,
              # use all data points in training because
              # we need to estimate cell abundance at all locations
              train_size=1,
              #use_gpu=True,
             )

    # In this section, we export the estimated cell abundance (summary of the posterior distribution).
    adata_vis = mod.export_posterior(
        adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
    )

    # Save model
    mod.save(f"{run_name}", overwrite=True)

    # mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

    # Save anndata object with results
    adata_file = f"{run_name}/sp.h5ad"
    adata_vis.write(adata_file)
    gc.collect()
    torch.cuda.empty_cache()
    gc.collect()

    # adata_file = f"{run_name}/sp.h5ad"
    # adata_vis = sc.read_h5ad(adata_file)
    # mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

