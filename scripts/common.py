"""Release-local configuration and labels. No parent-workspace imports."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
CONFIG = json.loads((ROOT / 'config.json').read_text())
CONTRASTS = [c['id'] for c in CONFIG['contrasts']]
LABELS = {c['id']: c['label'] for c in CONFIG['contrasts']}
TICKS = ['Whole-gene\nKO', 'Knock-in\nKO', 'Degron\n(matched\ncontrol)',
         'Degron\n(parental\ncontrol)', 'CRISPRi']
DE2_LABEL = 'Anterior-enriched\ndefinitive endoderm (DE2)'

def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def verify_inputs():
    records = json.loads((ROOT / 'provenance/source_manifest.json').read_text())
    for record in records:
        path = ROOT / record['release_path']
        if sha256(path) != record['sha256']:
            raise ValueError(f'Frozen input checksum mismatch: {path}')
    return len(records)

def load_sets():
    return {fields[0]: set(fields[2:]) for line in
            (ROOT / 'inputs/reference/cell_type_signatures.gmt').read_text().splitlines()
            if (fields := line.split('\t'))}
