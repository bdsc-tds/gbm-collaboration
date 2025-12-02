"""Script with utility functions to use with inferCNVpy."""

import re
import pandas as pd
import infercnvpy as cnv
import scanpy as sc
import anndata as ad
from loguru import logger
import matplotlib.pyplot as plt
from gbmcollaboration import PATH_GTF


class InferCNVAnalyzer:
    def __init__(self, adata: ad.AnnData):
        """Initializes the InferCNVAnalyzer class to use infercnvpy functions.

        Parameters:
        -----------
        adata: AnnData object
            Single cell object containing all samples and annotation columns
            for defining the reference cell population. adata.X must contain unnormalized counts
        """
        self.adata = adata
        self.counter = 0
        self.leiden_resolution = 1.0
        self.list_reference_cat: list = []
        self.tumor_clones = pd.DataFrame()
        self.bins_df = pd.DataFrame()

        # add genomic position from GTF file to the AnnData object
        logger.info("Adding genomic position from GTF file to the AnnData object.")
        cnv.io.genomic_position_from_gtf(PATH_GTF, self.adata)
        logger.info("Genomic position added successfully. Ready for next move.")

        # save raw counts
        logger.info("Saving raw counts from adata.X in adata.layers['raw_counts']")
        self.adata.layers['raw_counts'] = adata.X.copy()

        # log_normalize data
        self.log_normalize()
        logger.info("InferCNVAnalyzer initialized successfully.")

    def add_reference_column(
        self,
        cell_type_annotation_key,
        threshold_n_genes=400,
        threshold_n_counts=10000,
    ):
        """Add the reference column to the object.

        The reference column will be used in the reference_key argument of the inferCNVpy function.

        Parameters
        ----------
        self.adata : anndata object
            Annotated data object.
        cell_type_annotation_key : str
            Key in the obs to use as the cell type annotation for mapping.
        threshold_n_genes : int, optional
            Threshold to consider a cell as low_genes, by default 400.
        threshold_n_counts : int, optional
            Threshold to consider a cell as high_counts, by default 10000.

        Returns
        -------
        anndata object
            Annotated data object with the reference column.
        """
        logger.info(
            f"Using key {cell_type_annotation_key} to define the reference cell types."
        )
        # create a "reference" column in the obs to filter the healthy reference cells
        self.adata.obs["reference"] = self.adata.obs[cell_type_annotation_key].copy()

        # calculate n_genes and n_counts to further remove cells from the reference
        self.adata.obs["n_genes"] = (self.adata.layers['raw_counts'] > 0).sum(axis=1)
        self.adata.obs["n_counts"] = self.adata.layers['raw_counts'].sum(axis=1)

        self.adata.obs["reference"] = self.adata.obs["reference"].astype(str)

        # Use .loc to safely assign values based on conditions
        self.adata.obs.loc[
            self.adata.obs["n_genes"] < threshold_n_genes, "reference"
        ] = "low_genes"

        self.adata.obs.loc[
            self.adata.obs["n_counts"] > threshold_n_counts, "reference"
        ] = "high_counts"

        logger.info(
            f"Reference column added with following categories: "
            f"{self.adata.obs['reference'].unique().tolist()}"
        )

    def log_normalize(self):
        """Create a LogNormalize layer in the AnnData object.

        This layer is used for CNV calculations.
        """
        if "log1p" in self.adata.uns.keys():
            logger.warning(
                "adata.uns['log1p'] exists. adata.X is already Log Normalized. Skipping creation."
            )
            return

        logger.info("Log Normalizing adata.X.")
        sc.pp.normalize_total(
            self.adata, target_sum=1e4
        )
        sc.pp.log1p(self.adata)
        logger.info("Log Normalizing adata.X done.")

    def run_infercnv(self, list_reference_cat: list, window_size: int = 200, step: int = 10):
        """
        Run the infercnv algorithm to calculate CNV scores.

        Parameters:
        -----------
        list_reference_cat: list
            List of cell types to be used as reference for CNV calculation.
            Example for celltype level 2:
            ["MoMac", "T_NK", "Granulocyte", "Neuron", "Oligodendrocyte", "Endothelial", "Smooth_muscle"]
        window_size: int
            Size of the sliding window for CNV calculation.
        step: int
            Step size for the sliding window.
        """

        if "log1p" not in self.adata.uns:
            raise ValueError("adata.uns['log1p'] does not exist. adata.X is not Log Normalized.")

        # check if the reference column is present, other raise an error
        if "reference" not in self.adata.obs.columns:
            raise ValueError(
                "Reference column not found in obs. Please run "
                "add_reference_column() first."
            )
        self.list_reference_cat = list_reference_cat
        self.window_size = window_size
        self.step = step
        # launch inferCNVpy
        logger.info(
            f"Running inferCNV with reference categories: "
            f"{self.list_reference_cat}, "
            f"window size: {window_size}, step: {step}"
        )
        cnv.tl.infercnv(
            self.adata,
            reference_key="reference",
            reference_cat=self.list_reference_cat,
            window_size=self.window_size,
            step=self.step,
        )
        logger.info(
            "inferCNV run completed. Next steps: clustering and plotting."
        )

    def run_leiden_clustering(
        self,
        cor_method: str = "pearson",
        linkage_method: str = "complete",
        optimal_ordering: bool = True,
        resolution: float = 1.0
    ):
        """Run Leiden clustering on the CNV inference.

        Parameters:
        -----------
        resolution: float
            Resolution parameter for the Leiden algorithm. Higher values lead
            to more clusters.
        cor_method: str
            Correlation method to use for clustering. Default is "pearson".
        linkage_method: str
            Linkage method to use for clustering. Default is "complete".
        optimal_ordering: bool
            If True, the dendrogram will be ordered optimally. Default is True.
        """
        if self.counter == 0:  # meaning, never run before
            logger.info("Running CNV preprocessing steps.")
            cnv.tl.pca(self.adata)
            cnv.pp.neighbors(self.adata)

        # check if the leiden clustering has already been run
        if "cnv_leiden" in self.adata.obs.columns:
            logger.info("Leiden clustering already run. Checking resolution.")
            # check if the resolution is the same as the one already used
            if self.leiden_resolution == resolution:
                logger.info(
                    f"Leiden clustering already run with resolution "
                    f"{self.leiden_resolution}. "
                    "Skipping."
                )
            else:
                # if the resolution is different, we need to run leiden clustering again
                logger.info(
                    f"Changing Leiden resolution from {self.leiden_resolution} to {resolution}."
                )
                # if resolution is changed, we need to run leiden clustering again
                self.leiden_resolution = resolution
                logger.info(
                    f"Running leiden clustering with resolution {self.leiden_resolution}."
                )
                cnv.tl.leiden(self.adata, resolution=self.leiden_resolution,flavor="igraph", n_iterations=2)
        elif "cnv_leiden" not in self.adata.obs.columns:
            logger.info("Leiden clustering not run yet.")
            self.leiden_resolution = resolution
            cnv.tl.leiden(self.adata, resolution=self.leiden_resolution,flavor="igraph", n_iterations=2)

        self.counter += 1

        sc.tl.dendrogram(
            self.adata,
            groupby="cnv_leiden",
            use_rep="X_cnv_pca",
            cor_method=cor_method,
            linkage_method=linkage_method,
            optimal_ordering=optimal_ordering
        )

    def run_cnv_score(self):
        """Calculate the CNV score for each cluster.

        Clusters with a high score are likely affected by copy number abberations.
        The score is currently simply defined as the mean of the absolute values
        of the CNV scores in each cluster.
        """
        logger.info("Calculating CNV score.")
        cnv.tl.cnv_score(self.adata)
        logger.info("CNV score calculated successfully.")

    def composite_plots(
        self,
        save_path: str = None
    ):
        """Plot several visualizations.

        Parameters:
        -----------
        run_tl_umap: bool, optional
            If True, run the UMAP calculation on CNV profiles.
            Can be skipped if already run following leiden clustering.
        save_path: str, optional
            Path to save the figure. If None, the figure will not be saved.
        """
        if "cnv_leiden" not in self.adata.obs.columns:
            raise ValueError(
                "Leiden clustering not run yet. Please run "
                "run_leiden_clustering() first."
            )

        if 'X_umap' not in self.adata.obsm:
            logger.info("X_umap not found in adata.obsm. Running UMAP calculations.")
            sc.pp.pca(self.adata)
            sc.pp.neighbors(self.adata)
            sc.tl.umap(self.adata)

        if 'X_cnv_umap' not in self.adata.obsm:
            logger.info("X_cnv_umap not found in adata.obsm. Running UMAP calculations.")
            cnv.tl.umap(self.adata)

        if "cnv_score" not in self.adata.obs.columns:
            raise ValueError(
                "CNV score not calculated yet. Please run "
                "compute_cnv_score() first."
            )

        logger.info("Creating composite viz...")
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 12))
        # Plot normal umap colored by the cnv score or leiden clusters
        sc.pl.umap(
            self.adata,
            color="cnv_score",
            legend_loc="on data",
            legend_fontoutline=2,
            ax=ax1,
            show=False,
        )
        sc.pl.umap(self.adata, color="cnv_leiden", ax=ax2, show=False)

        # Plot UMAPs into specific axes
        cnv.pl.umap(self.adata, color="cnv_score", ax=ax3, show=False)

        cnv.pl.umap(
            self.adata,
            color="cnv_leiden",
            legend_loc="on data",
            legend_fontoutline=2,
            ax=ax4,
            show=False,
        )

        # Set the title for each subplot
        ax1.set_title("Gene expression UMAP colored by CNV score")
        ax2.set_title("Gene expression UMAP colored by CNV leiden clusters")
        ax3.set_title("CNV profiles UMAP colored by CNV scores")
        ax4.set_title("CNV profiles UMAP colored by CNV leiden clusters")

        # Reduce the title size
        ax1.title.set_size(16)
        ax2.title.set_size(16)
        ax3.title.set_size(16)
        ax4.title.set_size(16)

        # modify axes labels for ax3 and ax4
        ax3.set_xlabel("UMAP 1", fontsize=14)
        ax3.set_ylabel("UMAP 2", fontsize=14)
        ax4.set_xlabel("UMAP 1", fontsize=14)
        ax4.set_ylabel("UMAP 2", fontsize=14)

        # Adjust the layout
        fig.tight_layout()

        # Save the figure if a path is provided
        if save_path:
            fig.savefig(save_path, dpi=300)
            logger.info(f"Figure saved to {save_path}")

    def plot_gene_exp_umaps(
        self,
        cell_type_lvl: str,
        save_path: str,
        reset_umap: bool = True
    ):
        """Create UMAPs of gene expression colored by
        either cell types or leiden clusters.

        Parameters:
        -----------
        cell_type_lvl: str
            cell type level to display
        reset_umap: bool
            Recompute pca, neighbors and umap if needed, for example
            if two single cell datasets were merged.
            Default True.
        """
        if reset_umap:
            # Reset the pca, neighbors if wanted
            sc.tl.pca(self.adata)
            sc.pp.neighbors(self.adata)
            sc.tl.umap(self.adata)

        # Create umap composite figure
        fig, ((ax1, ax2)) = plt.subplots(1, 2, figsize=(11, 11))

        # Plot umaps colored by celltypes or leiden clusters
        sc.pl.umap(self.adata, color=cell_type_lvl, ax=ax1, show=False)
        sc.pl.umap(self.adata, color="cnv_leiden", ax=ax2, show=False)

        # Set the figure title
        fig.suptitle("Gene expression UMAPs", fontsize=20)

        # Adjust the layout
        fig.tight_layout()

        # Save the figure if a path is provided
        if save_path:
            fig.savefig(save_path, dpi=300)
            logger.info(f"Figure saved to {save_path}")

    def tumor_clone_annotation(self, save_path: str = None):
        """Return the cluster number and CNV score."""
        if "cnv_score" not in self.adata.obs.columns:
            raise ValueError(
                "CNV score not calculated yet. Please run "
                "compute_cnv_score() first."
            )
        if "cnv_leiden" not in self.adata.obs.columns:
            raise ValueError(
                "Leiden clustering not run yet. Please run "
                "run_leiden_clustering() first."
            )

        logger.info(
            "Annotating tumor clones based on "
            "leiden clustering and CNV scores."
        )
        # Create a new column in obs with the tumor clone annotation
        self.tumor_clones = pd.DataFrame(
            self.adata.obs[["cnv_leiden", "cnv_score"]].copy()
        )

        logger.info("Tumor clone annotation completed.")

        if save_path:
            self.tumor_clones.to_csv(save_path, index=True)
            logger.info(f"Tumor clone annotation saved to {save_path}")

        return self.tumor_clones

    def get_bin_positions(
        self,
        save_path: str = None
    ):
        """
        Function that mimic infercnvpy creation of the genomic window
        (or bins).


        Returns:
        --------
        pd.DataFrame:
            Returns a dataframe with one row per bin (starting at bin 0),
            with columns:
            - bin:         integer bin index (global, across chromosomes)
            - chromosome:  chromosome name
            - start:       genomic start coordinate of that bin
            - end:         genomic end coordinate of that bin

        Mirrors infercnvpy’s behavior:
        • For chromosomes with >= window_size genes, slide a window of size
            window_size in steps of `step` (mode='valid').
        • For chromosomes with < window_size genes, emit exactly one window
            spanning from first to last gene.
        """
        var = self.adata.var.copy()
        # replicate infercnvpy’s “natural” sort:
        def alphanum_key(s):
            return [int(t) if t.isdigit() else t.lower()
                    for t in re.split(r"([0-9]+)", s)]
        chromosomes = sorted(
            [c for c in var["chromosome"].unique()
            if isinstance(c, str) and c.startswith("chr") and c != "chrM"],
            key=alphanum_key,
        )

        bins = []
        bin_id = 0
        for chrom in chromosomes:
            df_chr = var[var["chromosome"] == chrom].sort_values("start")
            n_genes = df_chr.shape[0]

            if n_genes < self.window_size:
                # single “short” window
                bins.append({
                    "bin":        bin_id,
                    "chromosome": chrom,
                    "start":      int(df_chr["start"].iat[0]),
                    "end":        int(df_chr["end"].iat[-1]),
                })
                bin_id += 1
            else:
                # full sliding windows in steps of `step`
                n_valid = n_genes - self.window_size + 1
                for i in range(0, n_valid, self.step):
                    bins.append({
                        "bin":        bin_id,
                        "chromosome": chrom,
                        "start":      int(df_chr["start"].iat[i]),
                        "end":        int(df_chr["end"].iat[i + self.window_size - 1]),
                    })
                    bin_id += 1

        self.bins_df = pd.DataFrame(bins)

        logger.info("Genomic windows recomputed successfully.")

        if save_path:
            self.bins_df.to_csv(save_path, index=True)
            logger.info(f"Bins dataframe saved to {save_path}")

        return self.bins_df

    def run_ith_scores(self, groupby, label_key: str = None, label_to_score: str = None, save_path: str = None):
        """ Intratumoral heterogeneity scores based on CNAs and gene expressions """
        if label_key is not None:
            logger.info(f"Subsetting to {label_to_score} cells based on {label_key} cells for ITH scores")
            adata_score = self.adata[self.adata.obs[label_key] == label_to_score].copy()
        else:
            adata_score = self.adata

        logger.info("Computing ITH scores based on CNAs")
        cnv.tl.ithcna(adata_score, groupby, use_rep='X_cnv', key_added='ithcna', inplace=True)
        logger.info("Computing ITH scores based on GEX")
        cnv.tl.ithcna(adata_score, groupby, key_added='ithgex', inplace=True)

        if save_path is not None:
            adata_score.obs[['ithcna','ithgex']].to_csv(save_path)
