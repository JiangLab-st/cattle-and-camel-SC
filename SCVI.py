import sys
import warnings

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import scvi


adata = sc.read_h5ad("all5.h5ad")
print(adata)
adata.obs.tissue.value_counts()
query = np.array([s in ["caFFG", "caFAG", "caSG"] for s in adata.obs.tissue])
query
adata_ref = adata[~query].copy()
adata_query = adata[query].copy()
adata_ref
adata_query
sc.pp.highly_variable_genes(adata_ref, n_top_genes=2000, batch_key="tissue", subset=True)
adata_query = adata_query[:, adata_ref.var_names].copy()
scvi.model.SCVI.setup_anndata(adata_ref, batch_key="tissue", layer="counts")
arches_params = dict(use_layer_norm="both",use_batch_norm="none",encode_covariates=True,dropout_rate=0.2,n_layers=2)
vae_ref = scvi.model.SCVI(adata_ref, **arches_params)
#vae_ref = scvi.model.SCVI(adata_ref)
vae_ref.train()
adata_ref.obsm["X_scVI"] = vae_ref.get_latent_representation()
sc.pp.neighbors(adata_ref, use_rep="X_scVI")
sc.tl.leiden(adata_ref)
sc.tl.umap(adata_ref)
sc.pl.umap(adata_ref,color=["tissue", "celltypeB"], frameon=False, ncols=1)
plt.savefig("ref.pdf", bbox_inches = 'tight')
dir_path = "rumen_model/"
vae_ref.save(dir_path, overwrite=True)
#scvi.model.SCVI.prepare_query_anndata(adata_query, dir_path)
scvi.model.SCVI.prepare_query_anndata(adata_query, vae_ref)
#qury
vae_q = scvi.model.SCVI.load_query_data( adata_query, dir_path)
vae_q = scvi.model.SCVI.load_query_data(adata_query,vae_ref)
vae_q.train(max_epochs=200, plan_kwargs=dict(weight_decay=0.0))
adata_query.obsm["X_scVI"] = vae_q.get_latent_representation()
sc.pp.neighbors(adata_query, use_rep="X_scVI")
sc.tl.leiden(adata_query)
sc.tl.umap(adata_query)
sc.pl.umap(adata_query, color=["tissue", "celltypeB"], frameon=False, ncols=1)
plt.savefig("query.pdf", bbox_inches = 'tight')
#full
adata_full = adata_query.concatenate(adata_ref)
adata_full.obsm["X_scVI"] = vae_q.get_latent_representation(adata_full)
sc.pp.neighbors(adata_full, use_rep="X_scVI")
sc.tl.leiden(adata_full)
sc.tl.umap(adata_full)
sc.pl.umap(adata_full, color=["tissue", "celltypeB"], frameon=False, ncols=1)
plt.savefig("full.pdf", bbox_inches = 'tight')
