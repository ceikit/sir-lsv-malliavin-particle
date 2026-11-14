"""Puts the repository root on the import path, so that the scripts in this folder
can be run directly (python checks/<script>.py) without installing the sirlsv package."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
