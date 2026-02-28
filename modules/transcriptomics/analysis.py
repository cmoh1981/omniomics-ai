"""
Transcriptomics Analysis Module for OmniOmicsAI
==============================================
RNA-seq and single-cell RNA-seq analysis.
"""

import os
import logging
from typing import Dict, Any, List, Optional, Tuple
import json

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


class TranscriptomicsData:
    """Container for transcriptomics data."""

    def __init__(
        self,
        counts: Optional[pd.DataFrame] = None,
        metadata: Optional[pd.DataFrame] = None,
    ):
        self.counts = counts
        self.metadata = metadata
        self.normalized = None

    @classmethod
    def from_files(
        cls, counts_file: str, metadata_file: Optional[str] = None
    ) -> "TranscriptomicsData":
        """Load transcriptomics data from files."""
        logger.info(f"Loading counts from {counts_file}")

        # Determine file type
        if counts_file.endswith(".csv"):
            counts = pd.read_csv(counts_file, index_col=0)
        else:
            counts = pd.read_csv(counts_file, sep="\t", index_col=0)

        metadata = None
        if metadata_file and os.path.exists(metadata_file):
            if metadata_file.endswith(".csv"):
                metadata = pd.read_csv(metadata_file, index_col=0)
            else:
                metadata = pd.read_csv(metadata_file, sep="\t", index_col=0)

        return cls(counts=counts, metadata=metadata)

    def normalize(self, method: str = "log_cpm") -> pd.DataFrame:
        """
        Normalize count data.

        Args:
            method: Normalization method ('log_cpm', 'tmm', 'rle')

        Returns:
            Normalized DataFrame
        """
        if self.counts is None:
            raise ValueError("No counts data loaded")

        if method == "log_cpm":
            # Log2 CPM normalization
            cpm = self.counts.sum(axis=0) / 1e6
            normalized = np.log2((self.counts.div(cpm, axis=1) * 1e6) + 1)
        else:
            # Placeholder for other methods
            normalized = np.log2(self.counts + 1)

        self.normalized = normalized
        logger.info(f"Applied {method} normalization")

        return normalized

    def differential_expression(
        self,
        group1_samples: List[str],
        group2_samples: List[str],
        method: str = "deseq2",
    ) -> pd.DataFrame:
        """
        Perform differential expression analysis.

        Args:
            group1_samples: Sample names for group 1
            group2_samples: Sample names for group 2
            method: DE method ('deseq2', 'edger', 'limma')

        Returns:
            DataFrame with DE results
        """
        if self.counts is None:
            raise ValueError("No counts data loaded")

        # Filter to groups
        group1_cols = [
            c for c in self.counts.columns if any(s in c for s in group1_samples)
        ]
        group2_cols = [
            c for c in self.counts.columns if any(s in c for s in group2_samples)
        ]

        if not group1_cols or not group2_cols:
            raise ValueError("No matching samples found")

        results = []

        for gene in self.counts.index:
            g1 = self.counts.loc[gene, group1_cols].values.astype(float)
            g2 = self.counts.loc[gene, group2_cols].values.astype(float)

            # Remove zeros/NaN
            g1 = g1[g1 > 0]
            g2 = g2[g2 > 0]

            if len(g1) < 2 or len(g2) < 2:
                continue

            # Log transform
            g1 = np.log2(g1 + 1)
            g2 = np.log2(g2 + 1)

            # T-test
            stat, pval = stats.ttest_ind(g1, g2)

            # Log2 fold change
            log2fc = np.mean(g2) - np.mean(g1)

            results.append(
                {
                    "Gene": gene,
                    "Log2FC": log2fc,
                    "PValue": pval,
                    "Group1_Mean": np.mean(g1),
                    "Group2_Mean": np.mean(g2),
                }
            )

        df = pd.DataFrame(results)

        # FDR correction
        if len(df) > 0:
            reject, adj_pval, _, _ = multipletests(df["PValue"].values, method="bh")
            df["AdjPValue"] = adj_pval
            df["Significant"] = (df["AdjPValue"] < 0.05) & (abs(df["Log2FC"]) > 1)

        logger.info(f"DE analysis: {df['Significant'].sum()} significant genes")

        return df


class RNASeqPipeline:
    """RNA-seq analysis pipeline."""

    def __init__(
        self,
        data: Optional[TranscriptomicsData] = None,
        llm_analyzer: Optional[Any] = None,
    ):
        self.data = data
        self.llm_analyzer = llm_analyzer
        self.results = {}

    def run_analysis(
        self,
        group1_samples: List[str],
        group2_samples: List[str],
        normalize: bool = True,
    ) -> Dict[str, Any]:
        """Run complete RNA-seq analysis."""

        if normalize and self.data:
            self.data.normalize(method="log_cpm")

        # DE analysis
        de_results = self.data.differential_expression(
            group1_samples=group1_samples, group2_samples=group2_samples
        )

        self.results["de"] = de_results.to_dict()

        # Get summary
        n_sig = de_results["Significant"].sum()

        results = {
            "n_genes": len(de_results),
            "n_significant": int(n_sig),
            "n_up": int((de_results["Log2FC"] > 1).sum()),
            "n_down": int((de_results["Log2FC"] < -1).sum()),
        }

        # LLM analysis
        if self.llm_analyzer:
            llm_results = self.llm_analyzer.analyze_rnaseq(
                differential_results=de_results.to_dict()
            )
            results["llm_analysis"] = llm_results

        return results


# Example usage
if __name__ == "__main__":
    print("Transcriptomics module loaded.")
