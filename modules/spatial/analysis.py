"""
Spatial Omics Module for OmniOmicsAI
==================================
Analysis for spatial transcriptomics and proteomics.
"""

import os
import logging
from typing import Dict, Any, List, Optional, Tuple, Union
from pathlib import Path

import pandas as pd
import numpy as np
from scipy import sparse
from scipy.spatial.distance import cdist
import json

logger = logging.getLogger(__name__)


class SpatialData:
    """
    Container for spatial transcriptomics/proteomics data.
    """
    
    # Supported platforms
    PLATFORMS = [
        "visium", "visium_hd", "xenium", "merfish", 
        "cosmx", "slideseq", "stereoseq"
    ]
    
    def __init__(
        self,
        counts: Optional[pd.DataFrame] = None,
        positions: Optional[pd.DataFrame] = None,
        metadata: Optional[Dict] = None,
        platform: str = "visium"
    ):
        """
        Initialize spatial data container.
        
        Args:
            counts: Spot/cell by gene expression matrix
            positions: Spatial coordinates (x, y) for each spot
            metadata: Additional metadata
            platform: Spatial platform (visium, xenium, etc.)
        """
        self.counts = counts
        self.positions = positions
        self.metadata = metadata or {}
        self.platform = platform
        
        # Normalized data
        self.normalized = None
        
        # Clustering results
        self.clusters = None
        
        # Cell type annotations
        self.cell_types = None
    
    @classmethod
    def from_10x_visium(
        cls,
        counts_path: str,
        positions_path: Optional[str] = None,
        metadata_path: Optional[str] = None
    ) -> "SpatialData":
        """
        Load 10x Visium data.
        
        Args:
            counts_path: Path to counts matrix (CSV or mtx directory)
            positions_path: Path to spatial coordinates
            metadata_path: Path to metadata
            
        Returns:
            SpatialData object
        """
        logger.info(f"Loading Visium data from {counts_path}")
        
        # Load counts
        if counts_path.endswith(".csv"):
            counts = pd.read_csv(counts_path, index_col=0)
        elif counts_path.endswith(".mtx") or os.path.isdir(counts_path):
            # Would need to parse matrix market format
            logger.warning("MTX format not fully implemented")
            counts = None
        else:
            counts = pd.read_csv(counts_path, sep="\t", index_col=0)
        
        # Load positions if provided
        positions = None
        if positions_path and os.path.exists(positions_path):
            positions = pd.read_csv(positions_path, index_col=0)
        
        # Load metadata
        metadata = {}
        if metadata_path and os.path.exists(metadata_path):
            metadata = json.load(open(metadata_path))
        
        return cls(
            counts=counts,
            positions=positions,
            metadata=metadata,
            platform="visium"
        )
    
    @classmethod
    def from_xenium(
        cls,
        cell_by_gene_path: str,
        cell_locations_path: str
    ) -> "SpatialData":
        """
        Load Xenium data.
        
        Args:
            cell_by_gene_path: Cell by gene expression matrix
            cell_locations_path: Cell locations (X, Y, cell_id)
            
        Returns:
            SpatialData object
        """
        logger.info("Loading Xenium data")
        
        # Load expression matrix
        if cell_by_gene_path.endswith(".csv"):
            counts = pd.read_csv(cell_by_gene_path, index_col=0)
        else:
            counts = pd.read_csv(cell_by_gene_path, sep="\t", index_col=0)
        
        # Load cell locations
        positions = pd.read_csv(cell_locations_path)
        
        return cls(
            counts=counts,
            positions=positions,
            platform="xenium"
        )
    
    def normalize(
        self,
        method: str = "log_norm"
    ) -> pd.DataFrame:
        """
        Normalize spatial expression data.
        
        Args:
            method: Normalization method
            
        Returns:
            Normalized expression matrix
        """
        if self.counts is None:
            raise ValueError("No counts data")
        
        if method == "log_norm":
            # Log normalization (common for spatial)
            # Library size normalization first
            lib_sizes = self.counts.sum(axis=0)
            lib_sizes = lib_sizes.replace(0, 1)  # Avoid division by zero
            
            # Counts per million
            cpm = self.counts / lib_sizes * 1e6
            
            # Log transform
            normalized = np.log1p(cpm)
            
        elif method == "scran":
            # Placeholder for scran-like pooling
            logger.warning("Scran normalization not implemented, using log_norm")
            normalized = np.log1p(self.counts)
            
        else:
            normalized = self.counts
        
        self.normalized = normalized
        logger.info(f"Applied {method} normalization")
        
        return normalized
    
    def find_neighbors(
        self,
        radius: Optional[float] = None,
        n_neighbors: int = 10
    ) -> Dict[str, np.ndarray]:
        """
        Find spatial neighbors for each spot/cell.
        
        Args:
            radius: Maximum distance for neighbors (if None, use kNN)
            n_neighbors: Number of neighbors for kNN
            
        Returns:
            Dictionary with neighbor indices and distances
        """
        if self.positions is None:
            raise ValueError("No spatial coordinates")
        
        # Get coordinates
        coords = self.positions[["x", "y"]].values
        
        if radius is not None:
            # Distance-based neighbors
            distances = cdist(coords, coords)
            neighbors = np.where(distances <= radius)
            
            # Build sparse neighbor structure
            neighbor_dict = {}
            for i in range(len(coords)):
                idx = np.where(neighbors[0] == i)[0]
                neighbor_dict[i] = neighbors[1][idx]
            
        else:
            # kNN-based neighbors
            distances = cdist(coords, coords)
            
            # Get k nearest (excluding self)
            sorted_idx = np.argsort(distances, axis=1)
            neighbor_indices = sorted_idx[:, 1:n_neighbors+1]
            neighbor_distances = np.take_along_axis(
                distances, sorted_idx[:, 1:n_neighbors+1], axis=1
            )
            
            neighbor_dict = {
                i: neighbor_indices[i] 
                for i in range(len(coords))
            }
        
        logger.info(f"Found neighbors for {len(coords)} spots")
        
        return {
            "indices": neighbor_dict,
            "distances": neighbor_distances if radius is None else None
        }
    
    def spatial_clustering(
        self,
        n_clusters: int = 10,
        method: str = "leiden"
    ) -> np.ndarray:
        """
        Perform spatial clustering.
        
        Args:
            n_clusters: Number of clusters
            method: Clustering method
            
        Returns:
            Cluster labels
        """
        if self.normalized is None:
            self.normalize()
        
        # Simple k-means for now (would integrate with scanpy/leiden)
        from sklearn.cluster import KMeans
        
        # Use PCA for dimensionality reduction
        from sklearn.decomposition import PCA
        
        # PCA
        n_components = min(50, self.normalized.shape[1] - 1)
        pca = PCA(n_components=n_components)
        X_pca = pca.fit_transform(self.normalized.T)
        
        # K-means clustering
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        clusters = kmeans.fit_predict(X_pca)
        
        self.clusters = clusters
        
        logger.info(f"Spatial clustering complete: {n_clusters} clusters")
        
        return clusters
    
    def calculate_moran_i(
        self,
        gene: str,
        n_permutations: int = 1000
    ) -> Dict[str, float]:
        """
        Calculate Moran's I for spatial autocorrelation.
        
        Args:
            gene: Gene name to analyze
            n_permutations: Number of permutations
            
        Returns:
            Dictionary with Moran's I and p-value
        """
        if gene not in self.counts.index:
            raise ValueError(f"Gene {gene} not found")
        
        if self.positions is None:
            raise ValueError("No spatial coordinates")
        
        # Get expression values
        expr = self.counts.loc[gene].values
        
        # Get coordinates
        coords = self.positions[["x", "y"]].values
        
        # Calculate spatial weights (inverse distance)
        distances = cdist(coords, coords)
        np.fill_diagonal(distances, np.inf)
        
        # Weight matrix (inverse distance)
        weights = 1 / distances
        np.fill_diagonal(weights, 0)
        weights = weights / weights.sum(axis=1, keepdims=True)
        
        # Moran's I
        mean_expr = expr.mean()
        z = expr - mean_expr
        
        numerator = np.sum(weights * np.outer(z, z))
        denominator = np.sum(z ** 2)
        
        n = len(expr)
        moran_i = (n / weights.sum()) * (numerator / denominator)
        
        # Permutation test
        perm_morans = []
        for _ in range(n_permutations):
            perm_expr = np.random.permutation(expr)
            z_perm = perm_expr - mean_expr
            perm_num = np.sum(weights * np.outer(z_perm, z_perm))
            perm_morans.append((n / weights.sum()) * (perm_num / denominator))
        
        p_value = np.mean(np.array(perm_morans) >= moran_i)
        
        return {
            "moran_i": moran_i,
            "p_value": p_value,
            "gene": gene,
            "mean_expression": mean_expr
        }
    
    def spot_deconvolution(
        self,
        referenceProfiles: pd.DataFrame,
        method: str = "nnls"
    ) -> pd.DataFrame:
        """
        Deconvolute spot expression into cell type proportions.
        
        Args:
            referenceProfiles: Reference cell type profiles (genes x cell types)
            method: Deconvolution method
            
        Returns:
            Cell type proportions per spot
        """
        if self.normalized is None:
            self.normalize()
        
        # Find common genes
        common_genes = list(
            set(self.normalized.index) & set(referenceProfiles.index)
        )
        
        if len(common_genes) < 10:
            logger.warning("Too few common genes for deconvolution")
            return None
        
        # Subset to common genes
        spot_expr = self.normalized.loc[common_genes]
        ref = referenceProfiles.loc[common_genes]
        
        if method == "nnls":
            from scipy.optimize import nnls
            
            proportions = []
            for i in range(spot_expr.shape[1]):
                spot = spot_expr.iloc[:, i].values
                try:
                    prop, _ = nnls(ref.values, spot)
                    prop = prop / prop.sum()  # Normalize
                except:
                    prop = np.zeros(ref.shape[1])
                proportions.append(prop)
            
            result = pd.DataFrame(
                proportions,
                index=spot_expr.columns,
                columns=ref.columns
            )
        
        else:
            # Placeholder for other methods
            logger.warning(f"Method {method} not implemented")
            result = None
        
        logger.info(f"Deconvolution complete for {result.shape[0]} spots")
        
        return result


class SpatialAnalysis:
    """
    Complete spatial analysis pipeline.
    """
    
    def __init__(self, data: SpatialData, llm_analyzer: Optional[Any] = None):
        self.data = data
        self.llm_analyzer = llm_analyzer
        self.results = {}
    
    def run_complete_analysis(
        self,
        n_clusters: int = 10,
        genes_of_interest: Optional[List[str]] = None
    ) -> Dict[str, Any]:
        """
        Run complete spatial analysis.
        
        Args:
            n_clusters: Number of spatial clusters
            genes_of_interest: Genes for spatial autocorrelation
            
        Returns:
            Complete analysis results
        """
        logger.info("Running complete spatial analysis")
        
        results = {}
        
        # 1. Normalization
        self.data.normalize()
        results["normalization"] = "log_norm"
        
        # 2. Spatial clustering
        clusters = self.data.spatial_clustering(n_clusters=n_clusters)
        results["clusters"] = clusters.tolist()
        results["n_clusters"] = n_clusters
        
        # 3. Find neighbors
        neighbors = self.data.find_neighbors(n_neighbors=10)
        results["neighbors"] = "computed"
        
        # 4. Spatial autocorrelation for genes of interest
        if genes_of_interest:
            moran_results = []
            for gene in genes_of_interest[:10]:  # Limit for performance
                try:
                    moran = self.data.calculate_moran_i(gene)
                    moran_results.append(moran)
                except Exception as e:
                    logger.warning(f"Could not calculate Moran's I for {gene}: {e}")

            results["moran_i"] = moran_results
            results["moran_i"] = moran_results
        
        # 5. LLM analysis
        if self.llm_analyzer:
            llm_results = self.llm_analyzer.analyze_spatial(
                clusters=clusters,
                n_clusters=n_clusters,
                moran_results=results.get("moran_i", [])
            )
            results["llm_analysis"] = llm_results
        
        self.results = results
        
        logger.info("Spatial analysis complete")
        
        return results
    
    def cell_cell_communication(
        self,
        ligand_receptor_db: Optional[pd.DataFrame] = None
    ) -> pd.DataFrame:
        """
        Analyze cell-cell communication.
        
        Args:
            ligand_receptor_db: L-R database (optional)
            
        Returns:
            Communication scores
        """
        if self.data.normalized is None:
            self.data.normalize()
        
        # Placeholder - would integrate with CellChat/LIANA
        logger.info("Cell-cell communication analysis")
        
        # Simple correlation-based approach
        if self.data.clusters is not None:
            # Aggregate by cluster
            cluster_expr = pd.DataFrame([
                self.data.normalized.iloc[:, self.data.clusters == i].mean(axis=1)
                for i in range(self.data.clusters.max() + 1)
            ]).T
            
            # Calculate inter-cluster correlations (simplified)
            corr = cluster_expr.corr()
            
            return corr
        
        return None


# Spatial visualization helpers
def create_spatial_plot(
    data: SpatialData,
    gene: str,
    style: str = "heatmap"
) -> Dict[str, Any]:
    """
    Create spatial expression plot data.
    
    Args:
        data: SpatialData object
        gene: Gene to visualize
        style: Visualization style
        
    Returns:
        Plot data dictionary
    """
    if gene not in data.counts.index:
        raise ValueError(f"Gene {gene} not found")
    
    if data.positions is None:
        raise ValueError("No spatial coordinates")
    
    # Get expression values
    expr = data.counts.loc[gene].values
    
    # Get coordinates
    x = data.positions["x"].values
    y = data.positions["y"].values
    
    return {
        "x": x.tolist(),
        "y": y.tolist(),
        "expression": expr.tolist(),
        "gene": gene,
        "style": style
    }


# Example usage
if __name__ == "__main__":
    print("Spatial omics module loaded.")
