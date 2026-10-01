# NAME: scRNA_quiz2.py

import anndata as ad
import pooch
import scanpy as sc

sc.settings.set_figure_params(dpi=50, facecolor="white")

EXAMPLE_DATA = pooch.create(
    path=pooch.os_cache("scverse_tutorials"),
    base_url="doi:10.6084/m9.figshare.22716739.v1/",
)
EXAMPLE_DATA.load_registry_from_doi()

samples = {
    "s1d1": "s1d1_filtered_feature_bc_matrix.h5",
    "s1d3": "s1d3_filtered_feature_bc_matrix.h5",
}
adatas = {}

for sample_id, filename in samples.items():
    path = EXAMPLE_DATA.fetch(filename)
    sample_adata = sc.read_10x_h5(path)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()


# Quality Control
# mitochondrial genes
print("Quality Control")
adata.var["mt"] = adata.var_names.str.startswith("MT-")  # "MT-" for human, "Mt-" for mouse
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True)

# sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True)

# sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# Filter cells and genes
print("Filter Cells and Genes")
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

# Doublet detection
print("Doublet Detection")
sc.pp.scrublet(adata, batch_key="sample")

# Normalize
# Saving count data
print("Normalize")
adata.layers["counts"] = adata.X.copy()

# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data:
sc.pp.log1p(adata)

# Feature selection
print("Feature Selection")
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
#sc.pl.highly_variable_genes(adata)

# Dimension Reduction
print("Dimension Reduction")
sc.tl.pca(adata)
#sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# Visualization
print("Visualization")
sc.pp.neighbors(adata)
sc.tl.umap(adata)
# sc.pl.umap(adata, color="sample")

# Clustering
print("Clustering")
#sc.tl.leiden(adata, resolution = 0.2, partition_type="MutableVertexPartition", flavor="igraph")
sc.tl.leiden(adata, resolution = 0.8, n_iterations=2, flavor="igraph")
sc.pl.umap(adata, color=["leiden"])

print("")
print(">>>>> DONE <<<<<")













