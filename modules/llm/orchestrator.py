"""
LLM Orchestration Layer for OmniOmicsAI
=======================================
Unified interface for Gemini, Kilo, and Claude APIs for multiomics analysis.
"""

import os
import json
import logging
from typing import Optional, Dict, Any, List
from abc import ABC, abstractmethod
import requests

logger = logging.getLogger(__name__)


class LLMProvider(ABC):
    """Abstract base class for LLM providers."""

    @abstractmethod
    def generate(self, prompt: str, **kwargs) -> str:
        """Generate text from prompt."""
        pass

    @abstractmethod
    def analyze_data(self, data_summary: str, question: str) -> str:
        """Analyze omics data with natural language."""
        pass


class GeminiProvider(LLMProvider):
    """Google Gemini API provider."""

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or os.environ.get("GEMINI_API_KEY")
        self.base_url = "https://generativelanguage.googleapis.com/v1beta"
        self.model = "gemini-2.0-flash"

    def generate(self, prompt: str, **kwargs) -> str:
        """Generate text using Gemini API."""
        if not self.api_key:
            raise ValueError("GEMINI_API_KEY not provided")

        url = f"{self.base_url}/models/{self.model}:generateContent"

        headers = {"Content-Type": "application/json"}

        payload = {
            "contents": [{"parts": [{"text": prompt}]}],
            "generationConfig": {
                "temperature": kwargs.get("temperature", 0.7),
                "maxOutputTokens": kwargs.get("max_tokens", 2048),
                "topP": kwargs.get("top_p", 0.95),
            },
        }

        try:
            response = requests.post(
                url,
                headers=headers,
                params={"key": self.api_key},
                json=payload,
                timeout=120,
            )
            response.raise_for_status()
            result = response.json()

            return result["candidates"][0]["content"]["parts"][0]["text"]

        except requests.exceptions.RequestException as e:
            logger.error(f"Gemini API error: {e}")
            raise

    def analyze_data(self, data_summary: str, question: str) -> str:
        """Analyze omics data with context."""
        prompt = f"""You are an expert bioinformatics analyst specializing in proteomics, 
transcriptomics, and multi-omics integration.

DATA SUMMARY:
{data_summary}

QUESTION:
{question}

Provide a detailed, scientifically accurate analysis. Include:
1. Key findings from the data
2. Statistical interpretation
3. Biological significance
4. Recommendations for further analysis

Be specific and reference the data provided."""

        return self.generate(prompt, temperature=0.3)


class KiloProvider(LLMProvider):
    """Kilo API provider."""

    def __init__(self, api_key: Optional[str] = None):
        self.api_key = api_key or os.environ.get("KILO_API_KEY")
        self.base_url = "https://api.kilo.ai/api/gateway"
        self.model = "anthropic/claude-sonnet-4.5"

    def generate(self, prompt: str, **kwargs) -> str:
        """Generate text using Kilo API."""
        if not self.api_key:
            raise ValueError("KILO_API_KEY not provided")

        url = f"{self.base_url}/chat/completions"

        headers = {
            "Content-Type": "application/json",
            "Authorization": f"Bearer {self.api_key}",
        }

        payload = {
            "model": self.model,
            "messages": [{"role": "user", "content": prompt}],
            "temperature": kwargs.get("temperature", 0.7),
            "max_tokens": kwargs.get("max_tokens", 2048),
        }

        try:
            response = requests.post(url, headers=headers, json=payload, timeout=120)
            response.raise_for_status()
            result = response.json()

            return result["choices"][0]["message"]["content"]

        except requests.exceptions.RequestException as e:
            logger.error(f"Kilo API error: {e}")
            raise

    def analyze_data(self, data_summary: str, question: str) -> str:
        """Analyze omics data with context."""
        prompt = f"""You are an expert bioinformatics analyst specializing in proteomics, 
transcriptomics, and multi-omics integration.

DATA SUMMARY:
{data_summary}

QUESTION:
{question}

Provide a detailed, scientifically accurate analysis. Include:
1. Key findings from the data
2. Statistical interpretation
3. Biological significance
4. Recommendations for further analysis

Be specific and reference the data provided."""

        return self.generate(prompt, temperature=0.3)


class OmniLLM:
    """
    Unified LLM orchestrator for OmniOmicsAI.
    Provides a single interface for all LLM operations.
    """

    PROVIDERS = {
        "gemini": GeminiProvider,
        "kilo": KiloProvider,
    }

    def __init__(self, provider: str = "gemini", api_key: Optional[str] = None):
        """
        Initialize LLM orchestrator.

        Args:
            provider: LLM provider name ('gemini', 'kilo')
            api_key: API key for the provider
        """
        self.provider_name = provider.lower()

        if self.provider_name not in self.PROVIDERS:
            raise ValueError(
                f"Unknown provider: {provider}. "
                f"Available: {list(self.PROVIDERS.keys())}"
            )

        self.provider = self.PROVIDERS[self.provider_name](api_key)
        logger.info(f"Initialized {provider} LLM provider")

    def generate(self, prompt: str, **kwargs) -> str:
        """Generate text from prompt."""
        return self.provider.generate(prompt, **kwargs)

    def analyze_proteomics(
        self,
        differential_results: Dict[str, Any],
        volcano_plot_data: Optional[str] = None,
    ) -> Dict[str, Any]:
        """
        Analyze proteomics differential expression results.

        Args:
            differential_results: DataFrame or dict of DE results
            volcano_plot_data: Optional path or data for volcano plot

        Returns:
            Dictionary with analysis and recommendations
        """
        # Create summary
        n_up = sum(1 for r in differential_results.get("log2FC", []) if r > 1)
        n_down = sum(1 for r in differential_results.get("log2FC", []) if r < -1)
        n_sig = sum(1 for r in differential_results.get("adj_pvalue", []) if r < 0.05)

        summary = f"""
PROTEOMICS DIFFERENTIAL EXPRESSION ANALYSIS SUMMARY:
- Total proteins analyzed: {differential_results.get("n_proteins", "N/A")}
- Upregulated (log2FC > 1, p < 0.05): {n_up}
- Downregulated (log2FC < -1, p < 0.05): {n_down}
- Significant proteins: {n_sig}
- Top upregulated proteins: {", ".join(differential_results.get("top_up", [])[:5])}
- Top downregulated proteins: {", ".join(differential_results.get("top_down", [])[:5])}
"""

        question = """What are the key biological insights from this proteomics data?
Which pathways might be affected?
What follow-up experiments would you recommend?"""

        analysis = self.provider.analyze_data(summary, question)

        return {
            "summary": summary,
            "analysis": analysis,
            "statistics": {
                "n_upregulated": n_up,
                "n_downregulated": n_down,
                "n_significant": n_sig,
            },
        }

    def analyze_multiomics(
        self,
        omics_data: Dict[str, Dict[str, Any]],
        correlation_matrix: Optional[Dict] = None,
    ) -> Dict[str, Any]:
        """
        Analyze multi-omics integration results.

        Args:
            omics_data: Dictionary with omics type as key and data as value
            correlation_matrix: Optional cross-omics correlation data

        Returns:
            Dictionary with integrated analysis
        """
        summary_parts = []

        for omics_type, data in omics_data.items():
            summary_parts.append(
                f"\n{omics_type.upper()}:\n"
                f"- Features: {data.get('n_features', 'N/A')}\n"
                f"- Samples: {data.get('n_samples', 'N/A')}\n"
                f"- Variance explained: {data.get('var_explained', 'N/A')}%"
            )

        summary = f"MULTI-OMICS INTEGRATION ANALYSIS:{''.join(summary_parts)}"

        if correlation_matrix:
            summary += (
                f"\n\nCross-omics correlations detected: {len(correlation_matrix)}"
            )

        question = """What are the key insights from this multi-omics integration?
Which molecular pathways show cross-omics coordination?
What are the potential diagnostic or therapeutic implications?"""

        analysis = self.provider.analyze_data(summary, question)

        return {
            "summary": summary,
            "analysis": analysis,
            "recommendations": self._generate_recommendations(analysis),
        }

    def _generate_recommendations(self, analysis: str) -> List[str]:
        """Generate actionable recommendations from analysis."""
        prompt = f"""Based on this analysis, provide 5 specific, actionable 
recommendations for follow-up research or clinical application:

{analysis}

Format as a numbered list."""

        try:
            recommendations = self.generate(prompt, temperature=0.5)
            return [r.strip() for r in recommendations.split("\n") if r.strip()]
        except:
            return [
                "Review results with domain expert",
                "Validate findings with orthogonal methods",
            ]

    def generate_report(
        self, analysis_type: str, data: Dict[str, Any], format: str = "markdown"
    ) -> str:
        """
        Generate a complete analysis report.

        Args:
            analysis_type: Type of analysis (proteomics, transcriptomics, etc.)
            data: Analysis results and data
            format: Output format (markdown, html)

        Returns:
            Formatted report string
        """
        prompt = f"""Generate a professional analysis report for {analysis_type}.

DATA:
{json.dumps(data, indent=2)}

Format as a {format} report with:
1. Executive Summary
2. Methods
3. Results
4. Conclusions
5. References (where appropriate)

Use scientific writing style."""

        return self.generate(prompt, temperature=0.3)


def create_llm_provider(provider: str = "gemini", **kwargs) -> LLMProvider:
    """Factory function to create LLM provider."""
    if provider.lower() not in OmniLLM.PROVIDERS:
        raise ValueError(f"Unknown provider: {provider}")

    return OmniLLM.PROVIDERS[provider.lower()](**kwargs)


# Example usage
if __name__ == "__main__":
    # Initialize with Gemini
    llm = OmniLLM(provider="gemini", api_key=os.environ.get("GEMINI_API_KEY"))

    # Example analysis
    sample_data = {
        "n_proteins": 5000,
        "log2FC": [2.5, 1.8, -1.2, -2.1, 0.5, 3.2],
        "adj_pvalue": [0.001, 0.01, 0.05, 0.001, 0.3, 0.002],
        "top_up": ["TP53", "EGFR", "KRAS"],
        "top_down": ["PTEN", "RB1", "CDKN2A"],
    }

    result = llm.analyze_proteomics(sample_data)
    print(result["analysis"])
