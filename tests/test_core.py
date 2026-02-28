"""
Test Suite for OmniOmics AI
"""

import pytest
import pandas as pd
import numpy as np
from pathlib import Path
import tempfile
import json


class TestProteomicsModule:
    """Tests for proteomics analysis module"""

    @pytest.fixture
    def sample_proteomics_data(self):
        """Create sample proteomics data for testing"""
        np.random.seed(42)
        n_proteins = 100
        n_samples = 6

        # Create expression data
        data = pd.DataFrame(
            np.random.lognormal(mean=5, sigma=1, size=(n_proteins, n_samples)),
            index=[f"PROT_{i:04d}" for i in range(n_proteins)],
            columns=[f"Sample_{i}" for i in range(n_samples)],
        )

        # Add protein names
        data.index.name = "Protein"

        # Add additional columns like in MaxQuant
        data["Protein names"] = [f"Protein {i}" for i in range(n_proteins)]
        data["Gene names"] = [f"Gene {i}" for i in range(n_proteins)]

        return data

    def test_proteomics_data_creation(self, sample_proteomics_data):
        """Test creating ProteomicsData object"""
        from omniomics_ai.modules.proteomics.analysis import ProteomicsData

        # Save to temp file
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            sample_proteomics_data.to_csv(f, sep="\t")
            temp_file = f.name

        # Load data
        data = ProteomicsData.from_maxquant(temp_file)

        assert data is not None
        assert data.expression is not None
        assert len(data.expression) > 0

        # Cleanup
        Path(temp_file).unlink()

    def test_normalization(self, sample_proteomics_data):
        """Test normalization methods"""
        from omniomics_ai.modules.proteomics.analysis import ProteomicsData

        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            sample_proteomics_data.to_csv(f, sep="\t")
            temp_file = f.name

        data = ProteomicsData.from_maxquant(temp_file)

        # Test median normalization
        normalized = data.normalize(method="median")
        assert normalized is not None

        # Test mean normalization
        normalized = data.normalize(method="mean")
        assert normalized is not None

        # Cleanup
        Path(temp_file).unlink()


class TestTranscriptomicsModule:
    """Tests for transcriptomics analysis module"""

    @pytest.fixture
    def sample_transcriptomics_data(self):
        """Create sample transcriptomics data for testing"""
        np.random.seed(42)
        n_genes = 100
        n_samples = 6

        data = pd.DataFrame(
            np.random.lognormal(mean=5, sigma=1, size=(n_genes, n_samples)),
            index=[f"Gene_{i:04d}" for i in range(n_genes)],
            columns=[f"Sample_{i}" for i in range(n_samples)],
        )

        data.index.name = "Gene"
        return data

    def test_transcriptomics_data_creation(self, sample_transcriptomics_data):
        """Test creating TranscriptomicsData object"""
        from omniomics_ai.modules.transcriptomics.analysis import TranscriptomicsData

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            sample_transcriptomics_data.to_csv(f)
            temp_file = f.name

        data = TranscriptomicsData.from_counts(temp_file)

        assert data is not None
        assert data.counts is not None

        # Cleanup
        Path(temp_file).unlink()

    def test_log_cpm_normalization(self, sample_transcriptomics_data):
        """Test log-CPM normalization"""
        from omniomics_ai.modules.transcriptomics.analysis import TranscriptomicsData

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            sample_transcriptomics_data.to_csv(f)
            temp_file = f.name

        data = TranscriptomicsData.from_counts(temp_file)
        normalized = data.normalize()

        assert normalized is not None
        # Log-transformed data should have mostly negative values
        assert normalized.min().min() < 10  # CPM values should be reasonable

        # Cleanup
        Path(temp_file).unlink()


class TestSpatialModule:
    """Tests for spatial omics analysis module"""

    @pytest.fixture
    def sample_spatial_data(self):
        """Create sample spatial data for testing"""
        np.random.seed(42)
        n_spots = 100
        n_genes = 50

        # Create expression matrix
        expression = pd.DataFrame(
            np.random.lognormal(mean=3, sigma=1, size=(n_spots, n_genes)),
            index=[f"spot_{i:03d}" for i in range(n_spots)],
            columns=[f"Gene_{i:03d}" for i in range(n_genes)],
        )

        # Create spatial coordinates
        coordinates = pd.DataFrame(
            {
                "x": np.random.uniform(0, 10, n_spots),
                "y": np.random.uniform(0, 10, n_spots),
            },
            index=expression.index,
        )

        return expression, coordinates

    def test_spatial_data_creation(self, sample_spatial_data):
        """Test creating SpatialData object"""
        from omniomics_ai.modules.spatial.analysis import SpatialData

        expression, coordinates = sample_spatial_data

        data = SpatialData(
            expression=expression, coordinates=coordinates, platform="visium"
        )

        assert data is not None
        assert data.expression is not None
        assert data.coordinates is not None

    def test_spatial_normalization(self, sample_spatial_data):
        """Test spatial data normalization"""
        from omniomics_ai.modules.spatial.analysis import SpatialData

        expression, coordinates = sample_spatial_data
        data = SpatialData(expression=expression, coordinates=coordinates)

        normalized = data.normalize()

        assert normalized is not None


class TestMultiOmicsIntegration:
    """Tests for multi-omics integration"""

    def test_integrator_creation(self):
        """Test creating MultiOmicsIntegrator"""
        from omniomics_ai.modules.integration.multiomics import MultiOmicsIntegrator

        integrator = MultiOmicsIntegrator()

        assert integrator is not None
        assert hasattr(integrator, "omics_data")

    def test_add_omics(self):
        """Test adding omics data"""
        from omniomics_ai.modules.integration.multiomics import MultiOmicsIntegrator

        integrator = MultiOmicsIntegrator()

        # Add proteomics data
        proteomics_data = {
            "expression": pd.DataFrame(np.random.rand(10, 5)),
            "sample_info": pd.DataFrame({"condition": ["A", "B"] * 3}),
        }

        integrator.add_omics("proteomics", proteomics_data)

        assert "proteomics" in integrator.omics_data

    def test_pca_integration(self):
        """Test PCA integration method"""
        from omniomics_ai.modules.integration.multiomics import MultiOmicsIntegrator

        integrator = MultiOmicsIntegrator()

        # Add sample data
        data1 = pd.DataFrame(
            np.random.rand(20, 10), columns=[f"feat_{i}" for i in range(10)]
        )
        data2 = pd.DataFrame(
            np.random.rand(20, 8), columns=[f"feat_{i}" for i in range(8)]
        )

        integrator.add_omics("omics1", data1)
        integrator.add_omics("omics2", data2)

        result = integrator.integrate(method="pca")

        assert result is not None


class TestLLMOrchestrator:
    """Tests for LLM orchestration"""

    def test_omnillm_creation(self):
        """Test creating OmniLLM"""
        from omniomics_ai.modules.llm.orchestrator import OmniLLM

        llm = OmniLLM()

        assert llm is not None

    def test_provider_initialization(self):
        """Test LLM provider initialization"""
        from omniomics_ai.modules.llm.orchestrator import GeminiProvider, KiloProvider

        # Test with mock API keys (should not fail)
        gemini = GeminiProvider(api_key="test_key")
        kilo = KiloProvider(api_key="test_key")

        assert gemini is not None
        assert kilo is not None


class TestCLI:
    """Tests for CLI interface"""

    def test_cli_import(self):
        """Test CLI can be imported"""
        from omniomics_ai import cli

        assert cli is not None

    def test_cli_help(self):
        """Test CLI help command works"""
        from typer.testing import CliRunner
        from omniomics_ai.cli import app

        runner = CliRunner()
        result = runner.invoke(app, ["--help"])

        assert result.exit_code == 0
        assert "omniomics" in result.stdout.lower()


class TestAPIServer:
    """Tests for API server"""

    def test_api_import(self):
        """Test API can be imported"""
        from omniomics_ai.api import main

        assert main is not None

    def test_api_app_creation(self):
        """Test FastAPI app creation"""
        from omniomics_ai.api.main import app

        assert app is not None
        assert app.title == "OmniOmics AI API"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
