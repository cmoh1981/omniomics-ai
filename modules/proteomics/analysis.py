"""
Proteomics Analysis Module for OmniOmicsAI
===========================================
Comprehensive proteomics data processing and analysis.
"""

import os
import logging
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple, Union
import json

import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

logger = logging.getLogger(__name__)


class ProteomicsData:
    """
    Container for proteomics data with methods for processing and analysis.
    """

    def __init__(
        self,
        protein_groups: Optional[pd.DataFrame] = None,
        evidence: Optional[pd.DataFrame] = None,
        metadata: Optional[Dict] = None,
    ):
        """
        Initialize proteomics data container.

        Args:
            protein_groups: DataFrame with protein quantification data
            evidence: DataFrame with peptide-level evidence
            metadata: Dictionary with experiment metadata
        """
        self.protein_groups = protein_groups
        self.evidence = evidence
        self.metadata = metadata or {}
        self.intensity_columns = []
        self.sample_groups = {}

    @classmethod
    def from_maxquant(
        cls,
        protein_groups_path: str,
        evidence_path: Optional[str] = None,
        remove_contaminants: bool = True,
    ) -> "ProteomicsData":
        """
        Load data from MaxQuant output files.

        Args:
            protein_groups_path: Path to proteinGroups.txt
            evidence_path: Optional path to evidence.txt
            remove_contaminants: Whether to filter out contaminants

        Returns:
            ProteomicsData object
        """
        logger.info(f"Loading MaxQuant data from {protein_groups_path}")

        # Load protein groups
        protein_groups = pd.read_csv(protein_groups_path, sep="\t", low_memory=False)

        # Filter contaminants
        if remove_contaminants:
            original_count = len(protein_groups)
            # Check if columns exist before filtering
            filter_cols = ["Potential contaminant", "Reverse", "Only identified by site"]
            for col in filter_cols:
                if col in protein_groups.columns:
                    protein_groups = protein_groups[protein_groups[col] != "+"]
            logger.info(f"Filtered {original_count - len(protein_groups)} contaminants")

        # Load evidence if provided
        evidence = None
        if evidence_path and os.path.exists(evidence_path):
            evidence = pd.read_csv(evidence_path, sep="\t", low_memory=False)

        return cls(protein_groups=protein_groups, evidence=evidence)

    @classmethod
    def from_diann(cls, report_path: str) -> "ProteomicsData":
        """
        Load data from DIA-NN/Spectronaut output.

        Args:
            report_path: Path to report.tsv

        Returns:
            ProteomicsData object
        """
        logger.info(f"Loading DIA-NN data from {report_path}")

        # Load report
        report = pd.read_csv(report_path, sep="\t")

        # Pivot to protein-level matrix
        protein_groups = report.pivot_table(
            index="Protein.Group", columns="Run", values="PG.MaxLFQ", aggfunc="first"
        ).reset_index()

        protein_groups.columns.name = None

        return cls(protein_groups=protein_groups)

    def identify_intensity_columns(self) -> List[str]:
        """
        Identify intensity quantification columns.

        Returns:
            List of intensity column names
        """
        if self.protein_groups is None:
            return []

        # Common intensity column patterns
        patterns = ["Intensity ", "LFQ intensity ", "iBAQ ", "Max intensity "]

        self.intensity_columns = []
        for col in self.protein_groups.columns:
            for pattern in patterns:
                if pattern in col:
                    self.intensity_columns.append(col)
                    break

        logger.info(f"Identified {len(self.intensity_columns)} intensity columns")
        return self.intensity_columns

    def set_sample_groups(self, groups: Dict[str, str]):
        """
        Define experimental groups for samples.

        Args:
            groups: Dictionary mapping sample name to group
        """
        self.sample_groups = groups
        logger.info(f"Set {len(groups)} sample groups")

    def normalize(self, method: str = "median") -> pd.DataFrame:
        """
        Normalize protein intensities.

        Args:
            method: Normalization method ('median', 'mean', 'vsn', 'quantile')

        Returns:
            Normalized DataFrame
        """
        if not self.intensity_columns:
            self.identify_intensity_columns()

        if not self.intensity_columns:
            raise ValueError("No intensity columns found")

        data = self.protein_groups[self.intensity_columns].copy()

        if method == "median":
            # Median normalization
            medians = data.median()
            global_median = medians.median()
            normalization_factors = global_median / medians
            normalized = data * normalization_factors

        elif method == "mean":
            # Mean normalization
            means = data.mean()
            global_mean = means.mean()
            normalization_factors = global_mean / means
            normalized = data * normalization_factors

        elif method == "vsn":
            # Variance stabilizing normalization (placeholder)
            logger.warning("VSN not implemented, using log2 transformation")
            normalized = np.log2(data + 1)

        elif method == "quantile":
            # Quantile normalization
            normalized = self._quantile_normalize(data)

        else:
            raise ValueError(f"Unknown normalization method: {method}")

        # Replace in protein_groups
        result = self.protein_groups.copy()
        result[self.intensity_columns] = normalized

        logger.info(f"Applied {method} normalization")
        return result

    def _quantile_normalize(self, df: pd.DataFrame) -> pd.DataFrame:
        """Perform quantile normalization."""
        from scipy.stats import rankdata

        # Get ranks
        rank_mean = df.rank(method="average").mean(axis=1)

        # Sort by mean and assign average values
        sorted_df = df.apply(lambda x: rankdata(x), axis=0, result_type="expand")
        column_means = sorted_df.mean(axis=1)

        # Apply
        result = df.copy()
        for col in df.columns:
            result[col] = df[col].rank(method="average").map(column_means)

        return result

    def filter_by_coverage(
        self, min_proteins: int = 1000, min_unique_peptides: int = 1
    ) -> pd.DataFrame:
        """
        Filter proteins by identification quality.

        Args:
            min_proteins: Minimum number of proteins to keep
            min_unique_peptides: Minimum unique peptides

        Returns:
            Filtered DataFrame
        """
        if self.protein_groups is None:
            return None

        # Get peptide count column
        peptide_col = [c for c in self.protein_groups.columns if "Unique peptides" in c]

        if peptide_col:
            filtered = self.protein_groups[
                self.protein_groups[peptide_col[0]] >= min_unique_peptides
            ]
        else:
            filtered = self.protein_groups

        logger.info(
            f"Filtered from {len(self.protein_groups)} to {len(filtered)} proteins"
        )

        return filtered

    def differential_expression(
        self,
        group1: str,
        group2: str,
        group1_samples: List[str],
        group2_samples: List[str],
        method: str = "ttest",
        fdr_method: str = "fdr_bh",  # Benjamini-Hochberg for statsmodels
    ) -> pd.DataFrame:
        """
        Perform differential expression analysis.

        Args:
            group1: Name for first group
            group2: Name for second group
            group1_samples: List of sample names in group 1
            group2_samples: List of sample names in group 2
            method: Statistical test ('ttest', 'mannwhitney')
            fdr_method: FDR correction method

        Returns:
            DataFrame with DE results
        """
        if not self.intensity_columns:
            self.identify_intensity_columns()

        # Filter to groups that exist in data
        group1_cols = [
            c for c in self.intensity_columns if any(s in c for s in group1_samples)
        ]
        group2_cols = [
            c for c in self.intensity_columns if any(s in c for s in group2_samples)
        ]

        if not group1_cols or not group2_cols:
            raise ValueError("No matching samples found for groups")

        results = []

        for idx, row in self.protein_groups.iterrows():
            protein = row.get("Protein IDs", idx)

            g1_values = row[group1_cols].values.astype(float)
            g2_values = row[group2_cols].values.astype(float)

            # Remove NaN
            g1 = g1_values[~np.isnan(g1_values)]
            g2 = g2_values[~np.isnan(g2_values)]

            if len(g1) < 2 or len(g2) < 2:
                continue

            # Statistical test
            if method == "ttest":
                stat, pval = stats.ttest_ind(g1, g2)
            else:
                stat, pval = stats.mannwhitneyu(g1, g2, alternative="two-sided")

            # Log2 fold change
            log2fc = np.mean(g2) - np.mean(g1)

            results.append(
                {
                    "Protein": protein,
                    "Gene": row.get("Gene names", ""),
                    "Log2FC": log2fc,
                    "PValue": pval,
                    "Group1_Mean": np.mean(g1),
                    "Group2_Mean": np.mean(g2),
                    "Group1_N": len(g1),
                    "Group2_N": len(g2),
                }
            )

        df = pd.DataFrame(results)

        # FDR correction
        if len(df) > 0 and "PValue" in df.columns:
            reject, adj_pval, _, _ = multipletests(
                df["PValue"].values, method=fdr_method
            )
            df["AdjPValue"] = adj_pval
            df["Significant"] = (df["AdjPValue"] < 0.05) & (abs(df["Log2FC"]) > 1)

        logger.info(
            f"DE analysis: {df['Significant'].sum()} significant proteins "
            f"out of {len(df)}"
        )

        return df

    def get_top_proteins(
        self, de_results: pd.DataFrame, n: int = 50, direction: str = "both"
    ) -> Dict[str, List[str]]:
        """
        Get top differentially expressed proteins.

        Args:
            de_results: DataFrame from differential_expression()
            n: Number of proteins to return
            direction: 'up', 'down', or 'both'

        Returns:
            Dictionary with 'up' and 'down' protein lists
        """
        if de_results is None or len(de_results) == 0:
            return {"up": [], "down": []}

        sig = de_results[de_results.get("Significant", False)]

        if direction in ["up", "both"]:
            up = sig[sig["Log2FC"] > 0].nlargest(n, "Log2FC")
            top_up = up["Protein"].tolist()
        else:
            top_up = []

        if direction in ["down", "both"]:
            down = sig[sig["Log2FC"] < 0].nsmallest(n, "Log2FC")
            top_down = down["Protein"].tolist()
        else:
            top_down = []

        return {"up": top_up, "down": top_down}

    def pathway_enrichment(
        self, proteins: List[str], database: str = "kegg", organism: str = "hsa"
    ) -> Dict[str, Any]:
        """
        Perform pathway enrichment analysis.

        Args:
            proteins: List of protein/gene names
            database: Database to use ('kegg', 'go', 'reactome')
            organism: Organism code ('hsa', 'mmu', etc.)

        Returns:
            Dictionary with enrichment results
        """
        # This is a placeholder - would integrate with gseapy or enrichr
        logger.info(f"Running pathway enrichment for {len(proteins)} proteins")

        return {
            "status": "placeholder",
            "message": "Integrate with gseapy for actual enrichment",
            "n_proteins": len(proteins),
        }


class ProteomicsPipeline:
    """
    Complete proteomics analysis pipeline.
    """

    def __init__(
        self, data: Optional[ProteomicsData] = None, llm_analyzer: Optional[Any] = None
    ):
        """
        Initialize proteomics pipeline.

        Args:
            data: ProteomicsData object
            llm_analyzer: LLM analyzer for intelligent insights
        """
        self.data = data
        self.llm_analyzer = llm_analyzer
        self.results = {}

    def run_complete_analysis(
        self,
        group1_samples: List[str],
        group2_samples: List[str],
        group1_name: str = "Control",
        group2_name: str = "Treatment",
        normalize: bool = True,
        min_peptides: int = 1,
    ) -> Dict[str, Any]:
        """
        Run complete proteomics analysis workflow.

        Args:
            group1_samples: Sample names for first group
            group2_samples: Sample names for second group
            group1_name: Name for first group
            group2_name: Name for second group
            normalize: Whether to normalize data
            min_peptides: Minimum unique peptides

        Returns:
            Complete analysis results
        """
        logger.info("Starting complete proteomics analysis")

        results = {
            "groups": {
                "group1": {"name": group1_name, "samples": group1_samples},
                "group2": {"name": group2_name, "samples": group2_samples},
            }
        }

        # Step 1: Normalization
        if normalize and self.data:
            self.data.normalize(method="median")
            results["normalization"] = "median"

        # Step 2: Filtering
        if self.data:
            self.data.filter_by_coverage(min_unique_peptides=min_peptides)
            results["filtering"] = f"min_peptides={min_peptides}"

        # Step 3: Differential expression
        de_results = self.data.differential_expression(
            group1=group1_name,
            group2=group2_name,
            group1_samples=group1_samples,
            group2_samples=group2_samples,
        )

        results["differential_expression"] = de_results.to_dict()
        results["n_significant"] = int(de_results["Significant"].sum())

        # Step 4: Get top proteins
        top_proteins = self.data.get_top_proteins(de_results)
        results["top_proteins"] = top_proteins

        # Step 5: LLM analysis (if available)
        if self.llm_analyzer:
            llm_results = self.llm_analyzer.analyze_proteomics(
                differential_results=de_results.to_dict()
            )
            results["llm_analysis"] = llm_results

        self.results = results

        logger.info(
            f"Analysis complete: {results.get('n_significant', 0)} significant proteins"
        )

        return results

    def save_results(self, output_dir: str):
        """
        Save analysis results to files.

        Args:
            output_dir: Directory to save results
        """
        os.makedirs(output_dir, exist_ok=True)

        # Save DE results
        if "differential_expression" in self.results:
            de_df = pd.DataFrame(self.results["differential_expression"])
            de_df.to_csv(
                os.path.join(output_dir, "differential_expression.csv"), index=False
            )

        # Save summary
        summary = {
            k: v for k, v in self.results.items() if k != "differential_expression"
        }

        with open(os.path.join(output_dir, "analysis_summary.json"), "w") as f:
            json.dump(summary, f, indent=2, default=str)

        logger.info(f"Results saved to {output_dir}")


# Example usage
if __name__ == "__main__":
    # This would be used with actual data
    print("Proteomics module loaded. Use with MaxQuant or DIA-NN output.")
