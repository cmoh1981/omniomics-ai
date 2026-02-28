"""
Workflow Orchestration for OmniOmics AI
Simple Python-based workflow runner for common analysis pipelines
"""

from pathlib import Path
from typing import Dict, List, Optional, Any
import json
import logging
from dataclasses import dataclass
from enum import Enum

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class PipelineType(Enum):
    """Available pipeline types"""
    PROTEOMICS = "proteomics"
    TRANSCRIPTOMICS = "transcriptomics"
    SPATIAL = "spatial"
    MULTIOMICS = "multiomics"
    FULL = "full"


@dataclass
class PipelineConfig:
    """Configuration for a pipeline run"""
    input_dir: Path
    output_dir: Path
    design_file: Optional[Path] = None
    normalization: str = "median"
    llm_enabled: bool = True
    parallel: bool = True


class WorkflowOrchestrator:
    """Main workflow orchestration class"""
    
    def __init__(self, config: PipelineConfig):
        self.config = config
        self.output_dir = config.output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
    def run_proteomics_pipeline(
        self,
        input_file: Path,
        design_file: Optional[Path] = None
    ) -> Dict[str, Any]:
        """Run complete proteomics analysis pipeline"""
        from omniomics_ai.modules.proteomics.analysis import (
            ProteomicsData, ProteomicsPipeline
        )
        
        logger.info(f"Loading proteomics data from {input_file}")
        data = ProteomicsData.from_maxquant(str(input_file))
        
        if design_file:
            logger.info(f"Loading design from {design_file}")
            data.load_design(str(design_file))
        
        # Run pipeline
        pipeline = ProteomicsPipeline(data)
        
        logger.info("Running differential expression analysis")
        results = pipeline.run_differential_analysis(
            normalization=self.config.normalization
        )
        
        # Save results
        output_file = self.output_dir / "proteomics_results.json"
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2, default=str)
        
        logger.info(f"Results saved to {output_file}")
        return results
    
    def run_transcriptomics_pipeline(
        self,
        count_file: Path,
        design_file: Optional[Path] = None
    ) -> Dict[str, Any]:
        """Run complete transcriptomics analysis pipeline"""
        from omniomics_ai.modules.transcriptomics.analysis import (
            TranscriptomicsData, RNASeqPipeline
        )
        
        logger.info(f"Loading transcriptomics data from {count_file}")
        data = TranscriptomicsData.from_counts(str(count_file))
        
        if design_file:
           Loading design from {design_file}")
            data.load_design(str(design_file logger.info(f"))
        
        # Run pipeline
        pipeline = RNASeqPipeline(data)
        results = pipeline.run_differential_expression()
        
        # Save results
        output_file = self.output_dir / "transcriptomics_results.json"
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2, default=str)
        
        return results
    
    def run_spatial_pipeline(
        self,
        expression_file: Path,
        coordinates_file: Path,
        platform: str = "visium"
    ) -> Dict[str, Any]:
        """Run spatial omics analysis pipeline"""
        from omniomics_ai.modules.spatial.analysis import SpatialData, SpatialAnalysis
        
        logger.info(f"Loading spatial data from {expression_file}")
        data = SpatialData.from_files(
            expression_file=str(expression_file),
            coordinates_file=str(coordinates_file),
            platform=platform
        )
        
        # Run analysis
        analysis = SpatialAnalysis(data)
        results = analysis.run_full_analysis()
        
        # Save results
        output_file = self.output_dir / "spatial_results.json"
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2, default=str)
        
        return results
    
    def run_multiomics_pipeline(
        self,
        proteomics_file: Path,
        transcriptomics_file: Path,
        method: str = "pca"
    ) -> Dict[str, Any]:
        """Run multi-omics integration pipeline"""
        from omniomics_ai.modules.integration.multiomics import MultiOmicsIntegrator
        
        logger.info("Loading omics data for integration")
        
        # Load proteomics
        with open(proteomics_file) as f:
            proteomics_data = json.load(f)
        
        # Load transcriptomics
        with open(transcriptomics_file) as f:
            transcriptomics_data = json.load(f)
        
        # Integrate
        integrator = MultiOmicsIntegrator()
        integrator.add_omics("proteomics", proteomics_data)
        integrator.add_omics("transcriptomics", transcriptomics_data)
        
        results = integrator.integrate(method=method)
        
        # Save results
        output_file = self.output_dir / "multiomics_results.json"
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2, default=str)
        
        return results
    
    def run_full_pipeline(
        self,
        proteomics_file: Path,
        transcriptomics_file: Path,
        design_file: Optional[Path] = None
    ) -> Dict[str, Any]:
        """Run complete multi-omics pipeline with all analyses"""
        logger.info("Starting full pipeline")
        
        results = {}
        
        # 1. Proteomics
        logger.info("=" * 50)
        logger.info("Step 1: Proteomics Analysis")
        logger.info("=" * 50)
        proteomics_results = self.run_proteomics_pipeline(
            proteomics_file, design_file
        )
        results["proteomics"] = proteomics_results
        
        # 2. Transcriptomics
        logger.info("=" * 50)
        logger.info("Step 2: Transcriptomics Analysis")
        logger.info("=" * 50)
        transcriptomics_results = self.run_transcriptomics_pipeline(
            transcriptomics_file, design_file
        )
        results["transcriptomics"] = transcriptomics_results
        
        # 3. Multi-omics Integration
        logger.info("=" * 50)
        logger.info("Step 3: Multi-omics Integration")
        logger.info("=" * 50)
        proteomics_file = self.output_dir / "proteomics_results.json"
        transcriptomics_file = self.output_dir / "transcriptomics_results.json"
        
        multiomics_results = self.run_multiomics_pipeline(
            proteomics_file,
            transcriptomics_file,
            method="correlation"
        )
        results["multiomics"] = multiomics_results
        
        # 4. Summary
        logger.info("=" * 50)
        logger.info("Pipeline Complete!")
        logger.info("=" * 50)
        
        summary_file = self.output_dir / "pipeline_summary.json"
        with open(summary_file, "w") as f:
            json.dump(results, f, indent=2, default=str)
        
        logger.info(f"Summary saved to {summary_file}")
        
        return results


def run_workflow(
    pipeline_type: PipelineType,
    config: PipelineConfig,
    **kwargs
) -> Dict[str, Any]:
    """Main entry point for running workflows"""
    
    orchestrator = WorkflowOrchestrator(config)
    
    if pipeline_type == PipelineType.PROTEOMICS:
        return orchestrator.run_proteomics_pipeline(
            kwargs.get("input_file"),
            kwargs.get("design_file")
        )
    elif pipeline_type == PipelineType.TRANSCRIPTOMICS:
        return orchestrator.run_transcriptomics_pipeline(
            kwargs.get("input_file"),
            kwargs.get("design_file")
        )
    elif pipeline_type == PipelineType.SPATIAL:
        return orchestrator.run_spatial_pipeline(
            kwargs.get("expression_file"),
            kwargs.get("coordinates_file"),
            kwargs.get("platform", "visium")
        )
    elif pipeline_type == PipelineType.MULTIOMICS:
        return orchestrator.run_multiomics_pipeline(
            kwargs.get("proteomics_file"),
            kwargs.get("transcriptomics_file"),
            kwargs.get("method", "pca")
        )
    elif pipeline_type == PipelineType.FULL:
        return orchestrator.run_full_pipeline(
            kwargs.get("proteomics_file"),
            kwargs.get("transcriptomics_file"),
            kwargs.get("design_file")
        )
    else:
        raise ValueError(f"Unknown pipeline type: {pipeline_type}")


if __name__ == "__main__":
    # Example usage
    config = PipelineConfig(
        input_dir=Path("data"),
        output_dir=Path("results"),
        llm_enabled=True
    )
    
    # Run full pipeline
    results = run_workflow(
        PipelineType.FULL,
        config,
        proteomics_file=Path("data/proteomics.txt"),
        transcriptomics_file=Path("data/transcriptomics.csv"),
        design_file=Path("data/design.csv")
    )
