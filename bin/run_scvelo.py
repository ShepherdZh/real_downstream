#!/usr/bin/env python3
"""
scVelo RNA velocity analysis wrapper
Runs dynamical modeling on pre-processed single-cell data
"""
import os
print("Current working directory:", os.getcwd())
os.environ["NUMBA_CACHE_DIR"] = "./numba_cache"

import argparse
import warnings
warnings.filterwarnings('ignore')

import scanpy as sc
import scvelo as scv
import anndata as ad
import pandas as pd
from scipy.io import mmread
from pathlib import Path

# 先直接写绝对路径，还要改，base_dir和samplesheet都在之前的nextflow run的时候声明过，以后多个pipeline拼在一起，这些应该都是能走python的args或者是nextflow的params的
base_dir = Path("/opt/NAS/RD/zym/nf/results_lamanno/kallisto")
samplesheet = pd.read_csv('/opt/NAS/RD/zym/nf/samplesheet.csv')

def load_velocity_adata(sample_name):
    dir_path = base_dir / f"{sample_name}.count" / "counts_filtered"
    
    spliced = mmread(dir_path / "spliced.mtx").tocsr()
    unspliced = mmread(dir_path / "unspliced.mtx").tocsr()
    
    barcodes = pd.read_csv(dir_path / "spliced.barcodes.txt", header=None)[0].values
    genes = pd.read_csv(dir_path / "spliced.genes.txt", header=None)[0].values
    
    adata = ad.AnnData(
        X=spliced.copy(),               # 可以放 spliced（兼容旧习惯），或干脆 np.zeros_like(spliced) 让 .X 为空
        layers={'spliced': spliced, 'unspliced': unspliced},
        obs=pd.DataFrame(index=barcodes),
        var=pd.DataFrame(index=genes)
    )
    adata.obs['sample'] = sample_name
    return adata

sample_names = samplesheet['sample'].tolist()
adata_list = [load_velocity_adata(sample) for sample in sample_names]
adata = scv.utils.merge(*adata_list)

def main():
    parser = argparse.ArgumentParser(description='Run scVelo RNA velocity analysis')
    # input不写的，要自己做，给的没有spliced/unspliced这些
    parser.add_argument('--output', required=True, help='Output h5ad file with velocity results')
    parser.add_argument('--mode', default='dynamical', choices=['deterministic', 'dynamical', 'stochastic'],
                        help='Velocity estimation mode (default: dynamical)')
    parser.add_argument('--min_shared_counts', type=int, default=20,
                        help='Minimum shared counts for gene filtering')
    args = parser.parse_args()
    
    # Preprocessing for velocity
    print("Preprocessing for velocity analysis...")
    scv.pp.filter_and_normalize(adata, min_shared_counts=args.min_shared_counts)
    scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
    
    # Estimate velocity
    print(f"Computing RNA velocity using {args.mode} model...")
    if args.mode == 'dynamical':
        scv.tl.recover_dynamics(adata)
        scv.tl.velocity(adata, mode='dynamical')
    else:
        scv.tl.velocity(adata, mode=args.mode)
    
    # Compute velocity graph
    print("Computing velocity graph...")
    scv.tl.velocity_graph(adata)
    
    # Velocity embedding (project to UMAP if available)
    if 'X_umap' in adata.obsm:
        print("Computing velocity embedding on UMAP...")
        scv.tl.velocity_embedding(adata, basis='umap')
    
    # Compute additional metrics
    print("Computing velocity confidence and pseudotime...")
    scv.tl.velocity_confidence(adata)
    if args.mode == 'dynamical':
        scv.tl.latent_time(adata)
    scv.tl.velocity_pseudotime(adata)
    
    # 然后 embedding、umap 等（假设你已有 embedding）
    scv.pl.velocity_embedding_stream(adata, basis='umap')#, color='batch'
    
    # Save results
    print(f"Saving results to {args.output}...")
    adata.write_h5ad(args.output, compression='gzip')
    print("Done!")


if __name__ == '__main__':
    main()
