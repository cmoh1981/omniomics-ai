"""
CLI Interface for OmniOmics AI
Provides command-line interface for running analyses and managing the platform.
"""

import typer
from typing import Optional
from pathlib import Path
import json

app = typer.Typer(
    name="omniomics", help="OmniOmics AI - LLM-powered Multi-Omics Analysis Platform"
)

# Sub-apps for different commands
server_app = typer.Typer(help="Server management commands")
analysis_app = typer.Typer(help="Run analysis commands")
app.add_typer(server_app, name="server")
app.add_typer(analysis_app, name="analyze")


@server_app.command("dashboard")
def start_dashboard(
    port: int = typer.Option(8501, help="Port for Streamlit dashboard"),
    host: str = typer.Option("localhost", help="Host for dashboard"),
):
    """Start the Streamlit dashboard"""
    import subprocess
    import sys

    typer.echo(f"Starting Streamlit dashboard on {host}:{port}...")
    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        "omniomics_ai/ui/app.py",
        "--server.port",
        str(port),
        "--server.address",
        host,
    ]
    subprocess.run(cmd)


@server_app.command("api")
def start_api(
    host: str = typer.Option("0.0.0.0", help="API server host"),
    port: int = typer.Option(8000, help="API server port"),
    reload: bool = typer.Option(False, help="Enable auto-reload"),
):
    """Start the FastAPI server"""
    import subprocess
    import sys

    typer.echo(f"Starting API server on {host}:{port}...")
    cmd = [
        sys.executable,
        "-m",
        "uvicorn",
        "omniomics_ai.api.main:app",
        "--host",
        host,
        "--port",
        str(port),
    ]
    if reload:
        cmd.append("--reload")
    subprocess.run(cmd)


@analysis_app.command("proteomics")
def run_proteomics(
    input_file: Path = typer.Option(..., help="Input protein groups file (MaxQuant)"),
    design_file: Optional[Path] = typer.Option(None, help="Experimental design CSV"),
    output_dir: Path = typer.Option(Path("results"), help="Output directory"),
    normalize: str = typer.Option("median", help="Normalization method"),
):
    """Run proteomics differential expression analysis"""
    from omniomics_ai.modules.proteomics.analysis import (
        ProteomicsData,
        ProteomicsPipeline,
    )

    typer.echo(f"Loading proteomics data from {input_file}...")
    data = ProteomicsData.from_maxquant(str(input_file))

    if design_file:
        typer.echo(f"Loading experimental design from {design_file}...")
        data.load_design(str(design_file))

    typer.echo(f"Running proteomics pipeline with {normalize} normalization...")
    pipeline = ProteomicsPipeline(data)
    results = pipeline.run_differential_analysis(normalization=normalize)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "proteomics_results.json"

    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, default=str)

    typer.echo(f"Results saved to {output_file}")


@analysis_app.command("transcriptomics")
def run_transcriptomics(
    count_file: Path = typer.Option(..., help="Input count matrix"),
    design_file: Optional[Path] = typer.Option(None, help="Experimental design CSV"),
    output_dir: Path = typer.Option(Path("results"), help="Output directory"),
):
    """Run transcriptomics differential expression analysis"""
    from omniomics_ai.modules.transcriptomics.analysis import TranscriptomicsData

    typer.echo(f"Loading transcriptomics data from {count_file}...")
    data = TranscriptomicsData.from_counts(str(count_file))

    if design_file:
        typer.echo(f"Loading experimental design from {design_file}...")
        data.load_design(str(design_file))

    results = data.differential_expression()

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "transcriptomics_results.json"

    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, default=str)

    typer.echo(f"Results saved to {output_file}")


@analysis_app.command("multiomics")
def run_multiomics(
    proteomics_file: Path = typer.Option(..., help="Proteomics results JSON"),
    transcriptomics_file: Path = typer.Option(..., help="Transcriptomics results JSON"),
    output_dir: Path = typer.Option(Path("results"), help="Output directory"),
    method: str = typer.Option(
        "pca", help="Integration method (pca, correlation, mofa)"
    ),
):
    """Run multi-omics integration analysis"""
    from omniomics_ai.modules.integration.multiomics import MultiOmicsIntegrator

    typer.echo("Loading omics data...")
    with open(proteomics_file) as f:
        proteomics_data = json.load(f)
    with open(transcriptomics_file) as f:
        transcriptomics_data = json.load(f)

    typer.echo(f"Running multi-omics integration with {method}...")
    integrator = MultiOmicsIntegrator()
    integrator.add_omics("proteomics", proteomics_data)
    integrator.add_omics("transcriptomics", transcriptomics_data)

    results = integrator.integrate(method=method)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "multiomics_results.json"

    with open(output_file, "w") as f:
        json.dump(results, f, indent=2, default=str)

    typer.echo(f"Results saved to {output_file}")


@analysis_app.command("llm")
def run_llm_analysis(
    data_file: Path = typer.Option(..., help="Analysis results JSON"),
    query: str = typer.Option(..., help="Natural language query"),
    output_dir: Path = typer.Option(Path("results"), help="Output directory"),
):
    """Run LLM-powered analysis on omics results"""
    from omniomics_ai.modules.llm.orchestrator import OmniLLM

    typer.echo("Loading data...")
    with open(data_file) as f:
        data = json.load(f)

    typer.echo("Initializing LLM analyzer...")
    llm = OmniLLM()

    typer.echo(f"Running query: {query}")
    result = llm.analyze_proteomics(data, query)

    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "llm_analysis.txt"

    with open(output_file, "w") as f:
        f.write(result)

    typer.echo(f"Analysis saved to {output_file}")
    typer.echo(f"\nResult:\n{result}")


@app.command("init")
def init_project(
    project_name: str = typer.Argument(..., help="Name of the project"),
):
    """Initialize a new OmniOmics project with required directory structure"""
    from omniomics_ai import __version__

    typer.echo(f"Initializing OmniOmics AI project: {project_name}")
    typer.echo(f"Version: {__version__}")

    # Create directory structure
    dirs = [
        Path(project_name) / "data",
        Path(project_name) / "data/proteomics",
        Path(project_name) / "data/transcriptomics",
        Path(project_name) / "data/spatial",
        Path(project_name) / "results",
        Path(project_name) / "config",
    ]

    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
        typer.echo(f"  Created: {d}")

    typer.echo(f"\nProject '{project_name}' initialized successfully!")


@app.command("version")
def show_version():
    """Show version information"""
    from omniomics_ai import __version__

    typer.echo(f"OmniOmics AI v{__version__}")


if __name__ == "__main__":
    app()
