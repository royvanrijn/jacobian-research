#!/usr/bin/env sage-python
"""Bounded numerical proposals; rational group identities are checked separately."""
import sys,math
from pathlib import Path
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from sage.all import matrix,vector,RealField
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_endpoint_section_relations as control
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
geometry=SourceFileLoader('endpoint_relation_geometry',str(CAS/'prospective_half_lattice_v2.sage')).load_module()

def main():
    p=control.protocol();out=control.D/'proposals.json'
    if out.exists():raise FileExistsError('preserve finite relation proposals')
    data={'status':'RUNNING','rows':[]};checkpoint(out,data);real=RealField(p['precision_bits'])
    for row in cert.read(control.INPUT)['rows']:
        if row['status']!='CERTIFIED_SPECIALIZED_SUBGROUP':continue
        model=tuple(map(cert.F,row['curve']));basis=[tuple(map(cert.F,P)) for P in row['points']];generic=[tuple(map(cert.F,P)) for P in row['generic_points']];r=len(basis)
        gram,asym=geometry.canonical_height_gram(model,[*basis,*generic]);g=matrix(real,[[str(gram[i][j]) for j in range(r)] for i in range(r)]);relations=[]
        for j,Q in enumerate(generic):
            coordinates=g.solve_right(vector(real,[str(gram[i][r+j]) for i in range(r)]));fractions=[Fraction(str(v)).limit_denominator(p['maximum_denominator']) for v in coordinates];d=math.lcm(*(v.denominator for v in fractions));word=[int(v*d) for v in fractions]
            if d>p['maximum_denominator'] or any(abs(c)>p['maximum_relation_coefficient'] for c in word):relations.append({'section_index':j,'status':'UNKNOWN','reason':'bounded rational relation proposal unavailable'});continue
            relations.append({'section_index':j,'status':'PROPOSED_INTEGER_RELATION','denominator':d,'coefficients':word})
        data['rows'].append({'family':row['family'],'endpoint':row['endpoint'],'metric_gram':[[str(v) for v in row] for row in gram],'maximum_asymmetry':str(asym),'relations':relations});checkpoint(out,data);print('PROPOSED ENDPOINT',row['family'],row['endpoint'],flush=True)
    data['status']='COMPLETE_DECLARED_PROPOSALS';checkpoint(out,data)
if __name__=='__main__':main()
