import scanpy as sc
import scanpy.external as sce
import numpy as np
import pandas as pd
import os
import sys

wd='D:/Desktop'

os.chdir(wd)
sc.settings.verbosity = 3
sc.logging.print_header()
sc.settings.set_figure_params(dpi=600,facecolor='white',dpi_save=600)

MT_gene_file = open('D:/Desktop/MT_genes','r')
MT_genes = [i.replace('\n','') for i in MT_gene_file]

dat = sc.read_10x_mtx("filtered_feature_bc_matrix",var_names='gene_symbols',cache=True)
dat.var_names_make_unique()
dat

all_genes = dat.var['gene_ids']
dat.var['MT'] = [(True if i in MT_genes else False) for i in all_genes]
sc.pp.calculate_qc_metrics(dat, qc_vars=['MT'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(dat, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4,multi_panel=True,save='_summary.pdf')
sc.pl.scatter(dat, x='total_counts',y='pct_counts_mt',save='_totalCount_mt.pdf')
sc.pl.scatter(dat, x='total_counts',y='n_genes_by_counts',save='_totalCount_genes.pdf')


sc.pp.filter_cells(dat,min_genes=300)
sc.pp.filter_genes(dat,min_cells=3)
dat = dat[dat.obs.n_genes_by_counts < 5000,:]
dat = dat[dat.obs.total_counts<20000,:]
dat = dat[dat.obs.pct_counts_mt<15,:]
sc.external.pp.scrublet(dat,expected_doublet_rate=0.05)
dat = dat[dat.obs.predicted_doublet==False,:]
dat = dat[dat.obs.doublet_score<0.25,:]

# output the summary statistics of corresponding sample:
mtPct = round(np.mean(dat.obs.pct_counts_mt),2)
meanGenes = round(np.mean(dat.obs.n_genes),2)
meanUMI = round(np.mean(dat.obs.total_counts),2)

outfile = open('hippocampus_summary.txt','w')
outfile.write('hippocampus' + ' ' + str(dat.n_obs) + ' ' + str(dat.n_vars) + ' ')
outfile.write(str(mtPct) + ' ' + str(meanGenes) + ' ' + str(meanUMI) + ' ')
outfile.close()

# visualize the distribution of cell types
sc.pl.violin(dat, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True,save='_summary_clean.pdf')
sc.pl.scatter(dat, x='total_counts',y='pct_counts_mt',save='_totalCount_mt_clean.pdf')
sc.pl.scatter(dat, x='total_counts',y='n_genes_by_counts',save='_totalCount_genes_clean.pdf')
sc.pp.normalize_total(dat,target_sum = 1e4) # normalization
sc.pp.log1p(dat) # log transform

# identify highly variable genes
sc.pp.highly_variable_genes(dat,min_mean=0.0125,max_mean=3,min_disp=0.5)

sc.tl.pca(dat,svd_solver='arpack')
sce.pp.harmony_integrate(dat,'SampleID')
sc.pp.neighbors(dat,n_neighbors=10,n_pcs=30,use_rep='X_pca_harmony')
sc.tl.umap(dat)

for res in [1.4,1.5,1.6,1.7,1.8,1.9,2.0]:
    sc.tl.leiden(
        dat, key_added=f"leiden_res_{res:4.2f}", resolution=res
    )

sc.pl.umap(
    dat,
    color=["leiden_res_0.40", "leiden_res_0.50", "leiden_res_0.60","leiden_res_0.70","leiden_res_0.80","leiden_res_0.90","leiden_res_1.00"],
    legend_loc="on data", legend_fontsize=8,
    save='_leiden_mouse_hippo.pdf'
)

marker_genes_all = {
    "Oligodendrocyte progenitor cells": ["PDGFRA"],
    "Oligodendrocytes": ["MBP"],
    "Committed oligodendrocyte precursors": ["FYN","GPR17"],
    "Astrocytes/Bergmann glia": ["AQP4"],
    "Microglia": ["CALCR"],
    "Olfactory ensheathing cells": ["NPY"],
    "Ependymal cells": ["CFAP126"],
    "T cells": ["PTPRC"],
    "Perivascular macrophages": ["CD163"],
    "Vascular and leptomeningeal cells": ["COL1A2","COL1A1"],
    "Vascular smooth muscle cells": ["MYH11"],
    "Pericytes": ["RGS5"],
    "Vascular endothelial cells": ["CLDN5"],
    "Hormone producing cells/Folliculostellate cells": ["GH1","PRL"],
    "Neurons": ["SNAP25"],
}

cluster_mapping = {
      '50': 'Oligodendrocytes',
      '43': 'Oligodendrocytes',
      '16': 'Oligodendrocytes',
      '51': 'Oligodendrocytes',
      '17': 'Oligodendrocytes',
      '54': 'Oligodendrocytes',
      '31': 'Oligodendrocytes',
       '0': 'Oligodendrocytes',
       '2': 'Oligodendrocytes',
       '7': 'Oligodendrocytes',
      '49': 'Oligodendrocytes',
      '11': 'Oligodendrocyte progenitor cells',
      '44': 'Oligodendrocytes',
      '12': 'Oligodendrocytes',
      '25': 'Oligodendrocytes',
      '33': 'Microglia',
      '18': 'Microglia',
       '8': 'Microglia',
       '3': 'Microglia',
      '40': 'Microglia',
      '21': 'Perivascular macrophages',
      '56': 'Pituitary cells',
      '55': 'T cells',
      '35': 'T cells',
      '66': 'Oligodendrocytes',
      '30': 'Olfactory ensheathing cells',
      '65': 'Bergmann glia',
      '15': 'Bergmann glia',
      '59': 'Bergmann glia',
      '57': 'Oligodendrocytes',
      '34': 'Astrocytes',
      '64': 'Astrocytes',
      '36': 'Astrocytes',
       '1': 'Astrocytes',
      '23': 'Astrocytes',
      '58': 'Astrocytes',
      '60': 'Astrocytes',
      '45': 'Pituitary cells',
      '39': 'Vascular endothelial cells',
      '46': 'Vascular endothelial cells',
      '24': 'Vascular endothelial cells',
       '4': 'Vascular endothelial cells',
       '5': 'Vascular endothelial cells',
      '26': 'Vascular and leptomeningeal cells',
      '47': 'Vascular and leptomeningeal cells',
      '20': 'Pericytes',
      '38': 'Pericytes',
      '63': 'Pericytes',
      '13': 'Vascular smooth muscle cells',
      '29': 'Pericytes',
      '42': 'Pericytes',
      '61': 'Pituitary cells',
      '10': 'Neurons',
      '48': 'Oligodendrocyte progenitor cells',
      '14': 'Neurons',
      '28': 'Neurons',
      '41': 'Neurons',
      '62': 'Neurons',
      '22': 'Neurons',
       '6': 'Neurons',
      '27': 'Neurons',
      '53': 'Ependymal cells',
      '37': 'Pituitary cells',
       '9': 'Pituitary cells',
      '32': 'Pituitary cells',
      '19': 'Pituitary cells',
      '52': 'Pituitary cells'
}
dat.obs['celltype'] = dat.obs['leiden_res_2.00'].map(cluster_mapping).astype('category')

marker_genes_mnn = ["CALCR","CD163","PTPRC","CLDN5","RGS5","COL1A2","PDGFRA","MBP","NPY","AQP4","GH1","PRL",
                "CADM1","TOP2A","CFAP126"]

sc.pl.dotplot(dat,marker_genes_mnn, groupby="celltype", standard_scale="var",dendrogram=True,save='_mnn_celltypemarker.pdf')