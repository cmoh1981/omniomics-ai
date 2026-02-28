"""
Pytest configuration and fixtures
"""

import pytest
import sys
from pathlib import Path

# Add parent directory to path for imports
omniomics_path = Path(__file__).parent.parent
sys.path.insert(0, str(omniomics_path))
