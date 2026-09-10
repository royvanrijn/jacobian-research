#!/usr/bin/env python3
"""Immutable class-1 intake; no searches, candidate admission or V3 dispatch."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / 'artifacts/generated-results/elliptic-curves'
DEST = ART / 'x1092_class1_arithmetic_gate_v1'
PREFIX = 'x1092_class1_realization_'
PINNED_PARENT = '7c6ee40c46f5a1f3d1fc464b5685a0e4c77ed1f3347b67865d7c9f93b2acd862'
FILES = {'parent.json': 'compact_parent', 'realization-manifest.json': 'manifest',
         'section-certificate.json': 'sections', 'marking-certificate.json': 'marking'}

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def write_once(path, data):
    if path.exists():
        if path.read_bytes() != data:
            raise ValueError('immutable file differs: '+str(path))
    else:
        with path.open('xb') as f:
            f.write(data)

def frozen_manifest():
    parent = ART / (PREFIX+'compact_parent_v1.json')
    if sha(parent) != PINNED_PARENT:
        raise ValueError('certified class1 parent changed')
    realization = json.loads((ART/(PREFIX+'manifest_v1.json')).read_text())
    for name, digest in {**realization['proof_packets'], **realization['software']}.items():
        if sha(ROOT/name) != digest:
            raise ValueError('realization proof binding changed: '+name)
    return {'schema':'x1092.class1.arithmetic-freeze.v1',
      'status':'FROZEN_CERTIFIED_PARENT_GENERIC_INPUTS_ONLY',
      'parent_sha256':PINNED_PARENT,
      'files':{name:sha(ART/(PREFIX+stem+'_v1.json')) for name,stem in FILES.items()},
      'allowed_class_indices':[1], 'parameter_panel':'DISABLED',
      'next_gate':'AGENT1_PROSPECTIVE_STRICT_CLASS_NOVELTY_CHECKER',
      'required_order':['prospective_class_novelty','cover_construction',
                        'rational_solubility','independent_quotient_direction','V3'],
      'V3':'DISABLED_UNTIL_EXACT_RATIONAL_POINT_AND_QUOTIENT_CERTIFICATES',
      'generic_input_scope':'Exact cubic algebra, inherited generic section representatives and norm-square identities; no parameter selection, points on specialized fibres or arithmetic-space completeness claims.',
      'freeze_script_sha256':sha(Path(__file__))}

def freeze():
    manifest = frozen_manifest()
    DEST.mkdir(exist_ok=True)
    for name,stem in FILES.items():
        write_once(DEST/name,(ART/(PREFIX+stem+'_v1.json')).read_bytes())
    write_once(DEST/'freeze.json',(json.dumps(manifest,indent=2,sort_keys=True)+'\n').encode())
    return manifest

def check():
    manifest=json.loads((DEST/'freeze.json').read_text())
    if manifest != frozen_manifest():
        raise ValueError('freeze differs from certified source')
    for name,digest in manifest['files'].items():
        if sha(DEST/name)!=digest:
            raise ValueError('frozen input changed: '+name)
    return manifest

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true')
    a=p.parse_args();result=check() if a.check else freeze()
    print(result['status'],result['parent_sha256'])
