"""
Multi-Omics Integration Module for OmniOmicsAI
=============================================
Integrates proteomics, transcriptomics, and other omics data.
"""

import os
import logging
from typing import Dict, Any, List, Optional, Tuple
import json

import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from scipy.cluster.hierarchy import linkage, dendrogram
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

logger = logging.getLogger(__name__)


class MultiOmicsIntegrator:
    """
    Multi-omics data integration and analysis.
    """

    def __init__(
        self,
        omics_data: Dict[str, pd.DataFrame],
        sample_mapping: Optional[Dict[str, str]] = None,
    ):
        """
        Initialize multi-omics integrator.

        Args:
            Dictionary with omics_data: omics type as key and DataFrame as value
            sample_mapping: Optional mapping of sample names across omics
        """
        self.omics_data = omics_data
        self.sample_mapping = sample_mapping or {}
        self.integrated_data = None
        self.correlations = {}
        self.latent_factors = {}

    @classmethod
    def from_files(
        cls, files: Dict[str, str], sample_mapping: Optional[Dict[str, str]] = None
    ) -> "MultiOmicsIntegrator":
        """
        Load multi-omics data from files.

        Args:
            files: Dictionary mapping omics type to file paths
            sample_mapping: Optional sample name mapping

        Returns:
            MultiOmicsIntegrator instance
        """
        omics_data = {}

        for omics_type, filepath in files.items():
            logger.info(f"Loading {omics_type} from {filepath}")

            if filepath.endswith(".csv"):
                df = pd.read_csv(filepath)
            elif filepath.endswith(".tsv") or filepath.endswith(".txt"):
                df = pd.read_csv(filepath, sep="\t")
            elif filepath.endswith(".parquet"):
                df = pd.read_parquet(filepath)
            else:
                logger.warning(f"Unknown file format for {omics_type}")
                continue

            omics_data[omics_type] = df

        return cls(omics_data, sample_mapping)

    def harmonize_samples(
        self, method: str = "intersection"
    ) -> Dict[str, pd.DataFrame]:
        """
        Harmonize samples across omics datasets.

        Args:
            method: Harmonization method ('intersection', 'union', 'mapping')

        Returns:
            Dictionary of harmonized DataFrames
        """
        if method == "intersection":
            # Find common samples
            sample_sets = []

            for omics_type, df in self.omics_data.items():
                # Assume samples are columns or there's a sample column
                if "Sample" in df.columns:
                    samples = set(df["Sample"])
                else:
                    samples = set(df.columns)
                sample_sets.append(samples)

            common_samples = set.intersection(*sample_sets) if sample_sets else set()

            harmonized = {}
            for omics_type, df in self.omics_data.items():
                if "Sample" in df.columns:
                    harmonized[omics_type] = df[
                        df["Sample"].isin(common_samples)
                    ].set_index("Sample")
                else:
                    harmonized[omics_type] = df[list(common_samples)]

            logger.info(f"Harmonized to {len(common_samples)} common samples")
            return harmonized

        elif method == "mapping":
            # Use provided sample mapping
            harmonized = {}
            for omics_type, df in self.omics_data.items():
                # Apply sample mapping
                if self.sample_mapping:
                    df = df.rename(columns=self.sample_mapping)
                harmonized[omics_type] = df

            return harmonized

        return self.omics_data

    def calculate_correlations(
        self, omics1: str, omics2: str, method: str = "pearson"
    ) -> pd.DataFrame:
        """
        Calculate cross-omics correlations.

        Args:
            omics1: First omics type
            omics2: Second omics type
            method: Correlation method ('pearson' or 'spearman')

        Returns:
            DataFrame with correlations
        """
        if omics1 not in self.omics_data or omics2 not in self.omics_data:
            raise ValueError(f"Unknown omics type: {omics1} or {omics2}")

        df1 = self.omics_data[omics1]
        df2 = self.omics_data[omics2]

        # Get common features (if gene names in both)
        common_features = []
        if "Gene" in df1.columns and "Gene" in df2.columns:
            common_features = list(set(df1["Gene"]) & set(df2["Gene"]))

        if not common_features:
            # Use index intersection
            common_features = list(df1.index & df2.index)

        # Calculate correlations
        correlations = []

        for feature in common_features[:1000]:  # Limit to 1000 for performance
            try:
                vals1 = df1.loc[feature].values.astype(float)
                vals2 = df2.loc[feature].values.astype(float)

                # Remove NaN
                mask = ~(np.isnan(vals1) | np.isnan(vals2))

                if sum(mask) < 3:
                    continue

                if method == "pearson":
                    corr, pval = pearsonr(vals1[mask], vals2[mask])
                else:
                    corr, pval = spearmanr(vals1[mask], vals2[mask])

                correlations.append(
                    {
                        "Feature": feature,
                        "Correlation": corr,
                        "PValue": pval,
                        "N": sum(mask),
                    }
                )
            except Exception as e:
                logger.debug(f"Skipping {feature}: {e}")

        result = pd.DataFrame(correlations)

        if len(result) > 0:
            # FDR correction
            from statsmodels.stats.multitest import multipletests

            reject, adj_pval, _, _ = multipletests(result["PValue"].values, method="bh")
            result["AdjPValue"] = adj_pval

        logger.info(
            f"Calculated {len(result)} correlations between {omics1} and {omics2}"
        )

        return result

    def perform_pca(
        self, omics_types: Optional[List[str]] = None, n_components: int = 10
    ) -> Dict[str, Any]:
        """
        Perform PCA on combined or individual omics data.

        Args:
            omics_types: List of omics types to include (None = all)
            n_components: Number of principal components

        Returns:
            Dictionary with PCA results
        """
        if omics_types is None:
            omics_types = list(self.omics_data.keys())

        # Combine data matrices
        combined_dfs = []

        for omics_type in omics_types:
            df = self.omics_data[omics_type]

            # Extract numeric matrix
            if "Gene" in df.columns:
                numeric_df = df.set_index("Gene")
            else:
                numeric_df = df

            numeric_df = numeric_df.select_dtypes(include=[np.number])

            # Standardize
            scaler = StandardScaler()
            scaled = scaler.fit_transform(numeric_df.T).T

            combined_dfs.append(
                pd.DataFrame(
                    scaled,
                    index=numeric_df.index,
                    columns=[f"{omics_type}_{c}" for c in numeric_df.columns],
                )
            )

        # Concatenate
        combined = pd.concat(combined_dfs, axis=1)

        # Fill NaN
        combined = combined.fillna(0)

        # Perform PCA
        pca = PCA(n_components=min(n_components, combined.shape[1] - 1))
        pcs = pca.fit_transform(combined)

        results = {
            "pcs": pcs,
            "explained_variance": pca.explained_variance_ratio_,
            "cumulative_variance": np.cumsum(pca.explained_variance_ratio_),
            "samples": combined.columns.tolist(),
            "features": combined.index.tolist(),
        }

        logger.info(
            f"PCA complete: {results['cumulative_variance'][-1]:.1%} "
            f"variance explained by {n_components} components"
        )

        return results

    def cluster_samples(
        self, method: str = "hierarchical", metric: str = "euclidean"
    ) -> Dict[str, Any]:
        """
        Cluster samples based on multi-omics data.

        Args:
            method: Clustering method ('hierarchical', 'kmeans')
            metric: Distance metric

        Returns:
            Clustering results
        """
        # Combine all omics
        all_dfs = []

        for omics_type, df in self.omics_data.items():
            if "Gene" in df.columns:
                numeric_df = df.set_index("Gene")
            else:
                numeric_df = df

            numeric_df = numeric_df.select_dtypes(include=[np.number])
            all_dfs.append(numeric_df)

        combined = pd.concat(all_dfs, axis=1).fillna(0)

        # Standardize
        scaler = StandardScaler()
        scaled = scaler.fit_transform(combined.T).T

        if method == "hierarchical":
            linkage_matrix = linkage(scaled.T, method="ward", metric=metric)

            return {
                "linkage": linkage_matrix,
                "samples": combined.columns.tolist(),
                "method": "hierarchical",
            }

        return {"status": "placeholder"}

    def integrate_mofa(self, n_factors: int = 10) -> Dict[str, Any]:
        """
        Perform MOFA-style integration.

        Note: This is a simplified version. Use MOFA2 R package for production.

        Args:
            n_factors: Number of latent factors

        Returns:
            Integration results
        """
        logger.info("Running MOFA-style integration")

        # Combine all omics
        all_dfs = []

        for omics_type, df in self.omics_data.items():
            if "Gene" in df.columns:
                numeric_df = df.set_index("Gene")
            else:
                numeric_df = df

            numeric_df = numeric_df.select_dtypes(include=[np.number])

            # Take top variable features
            variances = numeric_df.var(axis=1).nlargest(500)
            numeric_df = numeric_df.loc[variances.index]

            all_dfs.append(numeric_df)

        combined = pd.concat(all_dfs, axis=1).fillna(0)

        # PCA for each omics
        omics_factors = {}

        for i, (omics_type, df) in enumerate(self.omics_data.items()):
            if "Gene" in df.columns:
                numeric_df = df.set_index("Gene")
            else:
                numeric_df = df

            numeric_df = numeric_df.select_dtypes(include=[np.number])

            pca = PCA(n_components=min(n_factors, numeric_df.shape[1] - 1))
            pcs = pca.fit_transform(numeric_df.T)

            omics_factors[omics_type] = {
                "pcs": pcs,
                "variance": pca.explained_variance_ratio_,
            }

        # Joint factors using all omics combined
        joint_pca = PCA(n_components=n_factors)
        joint_pcs = joint_pca.fit_transform(combined.T)

        return {
            "omics_factors": omics_factors,
            "joint_factors": joint_pcs,
            "joint_variance": joint_pca.explained_variance_ratio_,
            "samples": combined.columns.tolist(),
        }


class CrossOmicsAnalyzer:
    """
    Analyzer for cross-omics insights and interpretations.
    """

    def __init__(self, integrator: MultiOmicsIntegrator):
        self.integrator = integrator

    def find_coordinated_features(
        self, corr_threshold: float = 0.7
    ) -> List[Dict[str, Any]]:
        """
        Find features with high cross-omics correlation.

        Args:
            corr_threshold: Minimum correlation magnitude

        Returns:
            List of coordinated feature pairs
        """
        coordinated = []

        omics_types = list(self.integrator.omics_data.keys())

        for i in range(len(omics_types)):
            for j in range(i + 1, len(omics_types)):
                omics1, omics2 = omics_types[i], omics_types[j]

                correlations = self.integrator.calculate_correlations(
                    omics1, omics2, method="pearson"
                )

                if len(correlations) > 0:
                    high_corr = correlations[
                        abs(correlations["Correlation"]) > corr_threshold
                    ]

                    for _, row in high_corr.iterrows():
                        coordinated.append(
                            {
                                "feature": row["Feature"],
                                "omics1": omics1,
                                "omics2": omics2,
                                "correlation": row["Correlation"],
                                "pvalue": row["PValue"],
                            }
                        )

        return coordinated

    def generate_insights(self) -> Dict[str, Any]:
        """
        Generate comprehensive cross-omics insights.

        Returns:
            Dictionary with insights
        """
        insights = {
            "datasets": list(self.integrator.omics_data.keys()),
            "n_samples": {},
            "n_features": {},
            "coordinated_features": self.find_coordinated_features(),
        }

        for omics_type, df in self.integrator.omics_data.items():
            if "Gene" in df.columns:
                insights["n_features"][omics_type] = len(df)
                insights["n_samples"][omics_type] = len(df.columns)
            else:
                insights["n_features"][omics_type] = len(df.index)
                insights["n_samples"][omics_type] = len(df.columns)

        return insights


# Example usage
if __name__ == "__main__":
    print("Multi-omics integration module loaded.")
