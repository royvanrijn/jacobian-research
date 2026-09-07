#!/usr/bin/env sage-python
"""Portable exact coefficients of the verified global ordinate branches."""
import sys,argparse
from pathlib import Path
from sage.all import QQ,PolynomialRing
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ART=ROOT/'artifacts/generated-results/elliptic-curves';SOURCE=ART/'mestre_component_label_audit_v1.json';OUT=ART/'mestre_component_coherent_sections_v1.json'

def expected():
    r=cert.read(SOURCE);assert r['status']=='PASS' and r['coherent_stacked_rank']==11
    K=PolynomialRing(QQ,'u').fraction_field()
    def encode(s):
        f=K(s);return {'numerator':list(map(str,f.numerator().list())),'denominator':list(map(str,f.denominator().list()))}
    return {'schema':'elliptic-curves.mestre-coherent-sections.v1','source_sha256':cert.hashed(SOURCE),'exporter_sha256':cert.hashed(Path(__file__).resolve()),'fixed_source_root_order':r['fixed_source_root_order'],'roots':[encode(s) for s in r['normalized_root_functions']],
        'extra_ordinate_coefficients':[[encode(s) for s in y] for y in r['global_extra_ordinate_coefficients']],
        'generic_rank_lower_bound_certified':11,'boundary':'Fixed labelled roots and global rational ordinate branches. Fourteen supplied Jacobian points are not fourteen independent directions. Certify each specialized subgroup separately.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();r=expected()
    if a.check:assert r==cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve coherent section export')
        checkpoint(OUT,r)
    print('PASS portable coherent sections')
