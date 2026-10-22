"""Run a portable config with the exact source package for its experiment."""
from pathlib import Path
import argparse
import runpy
import os
import sys

root = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("profile", choices=["core", "fingerprint", "both-cold", "growth", "decoder"])
parser.add_argument("--cost-only", action="store_true", help="Load inputs and count fits; do not fit")
args = parser.parse_args()
source_dir = {"growth": "growth_src", "decoder": "decoder_src"}.get(args.profile, "src")
sys.path.insert(0, str(root / source_dir))
if args.profile == "decoder":
    # Frozen decoder source resolves its frozen-results path from the working directory.
    os.chdir(root)
sys.argv = ["profile_information.run", "--config", str(root / "configs" / (args.profile + ".json"))]
if args.cost_only:
    sys.argv.append("--cost-only")
module = "profile_information.decoder_sensitivity" if args.profile == "decoder" else "profile_information.run"
runpy.run_module(module, run_name="__main__")
