import os
import time
import datetime
import glob
import scanpy as sc
import numpy as np

tic = time.time()

################################################################################
dfp = '/gpfs/loomis/scratch60/xiting_yan/xy48/Pierce_CPB_rawdata/*/filtered*/'
out_fname = '/home/ngr4/scratch60/rpczi/rpczi_cpb.h5ad'
# modify temp_key, line 20
pp_recipe = False
################################################################################

print('Loading data...')
data_files = glob.glob(dfp)
if len(data_files) == 0:
    print('Path to data not found.\n... exiting.')
    exit()

# first load
print
adatas = {}
for i, file in enumerate(data_files):
    tic_sub = time.time()
################################################################################
    temp_key = file.split('/filtered')[0].split('/')[-1].split('-1_HHT')[0]
################################################################################
    if i==0:
        adata = sc.read_10x_mtx(file)
        batch_key = temp_key
        adata.var_names_make_unique()
    else:
        adatas[temp_key] = sc.read_10x_mtx(file)
        adatas[temp_key].var_names_make_unique()
    print('    loaded {} in {:.2f}-s \\ {:.0f}-s elapsed for {}/{} files'.format(temp_key, time.time() - tic_sub, time.time() - tic, i+1, len(data_files)))

adata = adata.concatenate(*adatas.values(), batch_categories=[batch_key]+list(adatas.keys()))
del adatas

if True:
    # save
    adata.write(out_fname)

# filter cells/genes, transform
sc.pp.calculate_qc_metrics(adata,inplace=True)
mito_genes = adata.var_names.str.startswith('MT-')
adata.obs['pmito'] = np.sum(adata[:, mito_genes].X, axis=1) / np.sum(adata.X, axis=1)
print('Ncells=%d have >10percent mt expression' % np.sum(adata.obs['pmito']>0.1))
print('Ncells=%d have <200 genes expressed' % np.sum(adata.obs['n_genes_by_counts']<200))
if pp_recipe:
    # actually filter
    adata = adata[adata.obs['pmito'] >= 0.2, :]
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3) # filtering cells gets rid of some genes of interest
    adata.raw = adata
    sc.pp.normalize_total(adata)
    sc.pp.sqrt(adata)

if True:
    # save
    adata.write(out_fname)

'''
# check data quality
varsoi = ['n_genes_by_counts', 'total_counts', 'pmito', 'pct_counts_in_top_200_genes']
p = sns.pairplot(adata.obs, hue='batch',
                 vars=varsoi,
                 plot_kws={'alpha':0.4, 'linewidth':0, 's':3},
                 corner=True)
p.savefig(os.path.join(pfp, 'pre-pp_qc.png'), dpi=600)

dt = adata.obs.loc[:, ['batch']+varsoi].groupby('batch').describe()
dt.loc[:, [(i, j) for i,j in dt.columns if j in ['count', 'mean', 'std']]]

# visualize data
sc.tl.pca(adata)
k = 30
sc.external.pp.bbknn(adata, n_pcs=50, neighbors_within_batch=int(k/len(adata.obs['batch'].unique())))
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.external.tl.phate(adata, gamma=0)
sc.pl.umap(adata, colors=['batch', 'leiden'])
sc.external.pl.phate(adata, colors=['batch', 'leiden'])
'''

print('... finished at {} after {}-min'.format(datetime.datetime.now().strftime('%y%m%d%H%M'), (time.time() - tic)/60))
