#!/usr/bin/env sage-python
"""Use the unchanged two-implementation conductor checker on the NFS supplement."""
import argparse
from importlib.machinery import SourceFileLoader
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
source = ROOT/'elliptic-curves/cas/certify_submitted627_630_conductors.sage'
loader = SourceFileLoader('submitted_conductor_proof',str(source))
spec = importlib.util.spec_from_loader(loader.name,loader)
proof = importlib.util.module_from_spec(spec)
loader.exec_module(proof)
proof.WORK = ROOT/'artifacts/local/elliptic-curves/submitted628-nfs-v1'
proof.OUT = ROOT/'artifacts/generated-results/elliptic-curves/submitted628_conductor_v2'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--check',action='store_true')
    args = parser.parse_args()
    protocol = json.loads((proof.WORK/'protocol.json').read_text())
    for path,digest in protocol['certificate_sources'].items():
        if proof.cert.hashed(ROOT/path) != digest:
            raise ArithmeticError('NFS supplement source changed')
    proof.prove(protocol['roster'][0],args.check)
