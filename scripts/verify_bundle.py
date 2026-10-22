"""Verify prepared-input and frozen-source hashes without fitting or network use."""
from pathlib import Path
import hashlib
import json

root = Path(__file__).resolve().parents[1]
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
inputs = json.loads((root / "provenance/prepared_inputs.json").read_text())["data_sha256"]
for name, expected in inputs.items():
    path = root / "data" / name
    if sha(path) != expected:
        raise SystemExit("Prepared input mismatch: " + name)
modules = json.loads((root / "provenance/source_code.json").read_text())["modules"]
for name, record in modules.items():
    if sha(root / name) != record["sha256"]:
        raise SystemExit("Frozen source mismatch: " + name)
stages = root / "provenance/replay_source_packages.json"
if stages.exists():
    for records in json.loads(stages.read_text())["packages"].values():
        for item in records:
            if sha(root / item["portable_path"]) != item["sha256"]:
                raise SystemExit("Frozen replay module mismatch: " + item["portable_path"])
print("Verified", len(inputs), "prepared input files and", len(modules), "frozen source modules; no fitting.")
