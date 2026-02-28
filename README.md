# 🧬 OmniOmics AI

LLM-Powered Multiomics Analysis Platform for proteomics, transcriptomics, and spatial omics integration.

## Features

### 🔬 Proteomics
- MaxQuant & DIA-NN data import
- Differential expression analysis
- LLM-powered result interpretation
- Pathway enrichment

### 🧬 Multi-omics Integration
- Cross-omics correlation analysis
- MOFA-style integration
- PCA & clustering
- LLM-guided insights

### 🤖 LLM Integration
- Gemini API support
- Kilo API support
- Natural language queries
- Automated report generation

## Installation

```bash
# Clone repository
git clone https://github.com/yourusername/omniomics_ai.git
cd omniomics_ai

# Create virtual environment
python -m venv venv
source venv/bin/activate  # Linux/Mac
# or
venv\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt
```

## Configuration

Set API keys as environment variables:

```bash
export GEMINI_API_KEY="your_gemini_api_key"
export KILO_API_KEY="your_kilo_api_key"
```

Or configure in `config/config.yaml`.

## Quick Start

### Run Dashboard

```bash
streamlit run omniomics_ai/ui/app.py
```

### Python API

```python
from omniomics_ai.modules.llm.orchestrator import OmniLLM
from omniomics_ai.modules.proteomics.analysis import ProteomicsPipeline

# Initialize LLM
llm = OmniLLM(provider="gemini")

# Analyze proteomics data
pipeline = ProteomicsPipeline(llm_analyzer=llm)
results = pipeline.run_complete_analysis(
    group1_samples=["control_1", "control_2"],
    group2_samples=["treatment_1", "treatment_2"]
)
```

## Project Structure

```
omniomics_ai/
├── config/              # Configuration files
├── modules/
│   ├── llm/          # LLM orchestration
│   ├── proteomics/   # Proteomics analysis
│   └── integration/  # Multi-omics integration
├── workflows/         # Analysis workflows
├── ui/               # Streamlit dashboard
├── api/              # REST API
└── tests/            # Unit tests
```

## Supported Data Formats

### Proteomics
- MaxQuant: `proteinGroups.txt`, `evidence.txt`
- DIA-NN: `report.tsv`
- Spectronaut: `report.tsv`

### Transcriptomics
- FeatureCounts
- Salmon/Kallisto
- HTSeq

### Multi-omics
- CSV, TSV, Parquet

## LLM Capabilities

| Task | Description |
|------|-------------|
| Data Query | Natural language questions about your data |
| Method Selection | Statistical approach recommendations |
| Result Interpretation | Explain analysis outputs |
| Literature Search | Find relevant papers |
| Hypothesis Generation | Propose testable hypotheses |
| Report Writing | Generate analysis summaries |

## Roadmap

- [x] Core LLM orchestration
- [x] Proteomics analysis module
- [x] Multi-omics integration
- [x] Streamlit dashboard
- [ ] Transcriptomics module
- [ ] Spatial omics module
- [ ] REST API
- [ ] Cloud deployment

## License

MIT License

## Citation

If you use OmniOmics AI in your research, please cite:

```
OmniOmics AI: LLM-Powered Multiomics Analysis Platform
```

## Contact

- GitHub Issues: https://github.com/yourusername/omniomics_ai/issues
- Email: contact@omniomics.ai
