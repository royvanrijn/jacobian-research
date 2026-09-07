#!/usr/bin/env sage-python
"""Propose numerical span coordinates, then prove every relation exactly."""
import json,sys
from pathlib import Path
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from sage.all import QQ,ZZ,RealField,EllipticCurve,matrix,lcm
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
D=ROOT/'artifacts/local/elliptic-curves/kihara-positive-seed-spans-v2'
def main():
    protocol=cert.read(D/'protocol.json')
    for p,h in protocol['sources'].items():assert cert.hashed(ROOT/p)==h
    geo=SourceFileLoader('positive_span_geometry',str(CAS/'prospective_half_lattice_v3.sage')).load_module();out=[]
    for i in protocol['seed_indices']:
        seed=cert.read(ROOT/'artifacts/local/elliptic-curves/kihara-positive-fibre-intake-v1'/('seed'+str(i)+'.json'));audit=cert.read(ROOT/'artifacts/local/elliptic-curves/kihara-positive-seed-gaps-v1'/('seed'+str(i)+'.json'));indices=audit['audits'][1]['independent_indices']
        model=tuple(map(cert.F,seed['curve']));cloud=[tuple(map(cert.F,p)) for p in seed['points']];E=EllipticCurve(QQ,list(map(QQ,model)));points=[E([QQ(c) for c in p]) for p in cloud];basis=[points[j] for j in indices]
        heights,asym=geo.canonical_height_gram(model,cloud);H=matrix(RealField(350),heights);B=H.matrix_from_rows_and_columns(indices,indices);relations=[]
        torsion=[E(0)]+[E([QQ(c) for c in p]) for p in audit['two_torsion_points']]
        for j,P in enumerate(points):
            numerical=B.solve_right(H.matrix_from_rows_and_columns(indices,[j]).column(0));coeff=[QQ(Fraction(str(c)).limit_denominator(128)) for c in numerical];den=lcm(c.denominator() for c in coeff);word=[ZZ(c*den) for c in coeff]
            residual=den*P-sum((c*Q for c,Q in zip(word,basis)),E(0));assert residual in torsion
            relations.append({'point_index':j,'multiplier':int(den),'basis_word':list(map(int,word)),'torsion_remainder':list(map(str,residual))})
        row={'id':seed['id'],'status':'PASS','span_rank':len(indices),'basis_indices':indices,'relations':relations,'maximum_numerical_asymmetry':str(asym)}
        with (D/('seed'+str(i)+'.json')).open('x') as f:json.dump(row,f,indent=2);f.write('\n')
        out.append(row);print(seed['id'],'EXACT SPAN',len(indices),'relations',len(relations),flush=True)
    with (D/'result.json').open('x') as f:json.dump({'status':'PASS','rows':out,'sources':protocol['sources'],'scope':'Every input point lies in the rational span of the independently certified odd-prime basis, modulo exact2-torsion. Numerical heights only propose rational words; exact elliptic group equations prove all relations. No point search or complete curve rank assertion.'},f,indent=2);f.write('\n')
if __name__=='__main__':main()
