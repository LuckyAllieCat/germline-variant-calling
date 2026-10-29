import sys
from pathlib import Path

# Make scripts/ importable from the tests
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
