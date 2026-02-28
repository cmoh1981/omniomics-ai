"""
OmniOmicsAI - Main Entry Point
================================
LLM-Powered Multiomics Analysis Platform
"""

import os
import sys
import logging
from pathlib import Path

# Add project root to path
project_root = Path(__file__).parent
sys.path.insert(0, str(project_root))

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

from omniomics_ai.modules.llm.orchestrator import OmniLLM
from omniomics_ai.modules.proteomics.analysis import ProteomicsPipeline
from omniomics_ai.modules.integration.multiomics import MultiOmicsIntegrator


def main():
    """Main entry point for OmniOmicsAI."""

    print("""
    ╔═══════════════════════════════════════════════════════════╗
    ║                                                           ║
    ║   🧬  OmniOmics AI  -  LLM-Powered Multiomics          ║
    ║                                                           ║
    ║   Version: 0.1.0                                         ║
    ║   Description: Comprehensive proteomics & multiomics   ║
    ║                analysis platform with LLM integration      ║
    ║                                                           ║
    ╚═══════════════════════════════════════════════════════════╝
    """)

    # Check for API keys
    gemini_key = os.environ.get("GEMINI_API_KEY")
    kilo_key = os.environ.get("KILO_API_KEY")

    if gemini_key:
        print("✓ Gemini API configured")
    else:
        print("⚠ GEMINI_API_KEY not set (optional)")

    if kilo_key:
        print("✓ Kilo API configured")
    else:
        print("⚠ KILO_API_KEY not set (optional)")

    print("\n" + "=" * 60)
    print("Available Commands:")
    print("=" * 60)
    print("1. dashboard  - Launch Streamlit dashboard")
    print("2. analyze   - Run proteomics analysis")
    print("3. integrate - Run multi-omics integration")
    print("4. llm       - Interactive LLM analysis")
    print("\n")

    # Import CLI
    try:
        import typer
        from omniomics_ai.cli import app

        app()
    except ImportError:
        # Fallback to simple CLI
        cmd = input("Enter command (dashboard/analyze/integrate/llm): ").strip().lower()

        if cmd == "dashboard":
            print("\nStarting Streamlit dashboard...")
            os.system("streamlit run omniomics_ai/ui/app.py")

        elif cmd == "analyze":
            print("\nRunning proteomics analysis...")
            # Placeholder for actual analysis
            print("Configure data path and run from dashboard")

        elif cmd == "integrate":
            print("\nRunning multi-omics integration...")
            print("Configure data files and run from dashboard")

        elif cmd == "llm":
            print("\nStarting interactive LLM...")

            llm = OmniLLM(provider="gemini")

            while True:
                query = input("\nYou: ").strip()

                if query.lower() in ["exit", "quit", "q"]:
                    break

                if query:
                    try:
                        response = llm.generate(query)
                        print(f"\nLLM: {response}")
                    except Exception as e:
                        print(f"Error: {e}")

        else:
            print(f"Unknown command: {cmd}")


if __name__ == "__main__":
    main()
