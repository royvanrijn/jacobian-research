#!/usr/bin/env sage-python
"""One bounded height proposal for the unresolved nineteen-point anchor span."""
import sys,math
from pathlib import Path
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from sage.all import matrix,vector,RealField,QQ,EllipticCurve
ROOT=Path(__file__).resolve().parents[2];CAS=ROOT/'elliptic-curves/cas';sys.path.insert(0,str(CAS))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
geometry=SourceFileLoader('carrier_anchor_relation_geometry',str(CAS/'prospective_half_lattice_v2.sage')).load_module()
ART=ROOT/'artifacts/generated-results/elliptic-curves';D=ROOT/'artifacts/local/elliptic-curves/carrier-anchor-relation-v1'
INPUT=ART/'soluble_pair_carrier_height_v1.json'
if __name__=='__main__':
    if (D/'protocol.json').exists():raise FileExistsError('preserve single bounded anchor relation proposal')
    paths=[Path(__file__).resolve(),INPUT,CAS/'prospective_half_lattice_v2.sage']
    p={'sources':{str(q.relative_to(ROOT)):cert.hashed(q) for q in paths},'word':[1,1],
        'maximum_denominator':64,'maximum_coefficient':64,'precision_bits':384,'wall_seconds':120,
        'scope':'The two supplied points plus17 generic sections have lower bound18 modulo2,3,5 at the existing anchor. Propose one rational expression for the second supplied point in the certified18-point basis. Exact rational group equality must verify it. No parameter expansion, point search, whole-curve rank or automatic retry.'}
    checkpoint(D/'protocol.json',p)
    row=next(r for r in cert.read(INPUT)['rows'] if r['word']==p['word']);points=row['independent_points'];target=row['supplied_points'][1]
    if points!=row['generic_points']+[row['supplied_points'][0]] or len(points)!=18:raise ArithmeticError('fixed reference basis differs')
    model=tuple(map(cert.F,row['curve']));allpoints=[tuple(map(cert.F,P)) for P in points+[target]]
    gram,asym=geometry.canonical_height_gram(model,allpoints);real=RealField(p['precision_bits'])
    g=matrix(real,[[str(gram[i][j]) for j in range(18)] for i in range(18)])
    coordinates=g.solve_right(vector(real,[str(gram[i][18]) for i in range(18)]))
    fractions=[Fraction(str(v)).limit_denominator(p['maximum_denominator']) for v in coordinates]
    denominator=math.lcm(*(v.denominator for v in fractions));word=[int(v*denominator) for v in fractions]
    result={'status':'UNKNOWN','protocol_sha256':cert.hashed(D/'protocol.json'),'curve':row['curve'],
        'parameter':row['compact_parameter'],'basis':points,'target':target,
        'denominator':denominator,'coefficients':word,'rank_certificate':row['rank_certificate'],
        'numerical_coordinates':list(map(str,coordinates)),'maximum_gram_asymmetry':str(asym)}
    if denominator<=p['maximum_denominator'] and max(map(abs,word))<=p['maximum_coefficient']:
        E=EllipticCurve(QQ,row['curve']);basis=[E([QQ(v) for v in P]) for P in points];Q=E([QQ(v) for v in target])
        if denominator*Q==sum((c*P for c,P in zip(word,basis)),E(0)):result['status']='EXACT_GROUP_RELATION_PENDING_INDEPENDENT_REPLAY'
    checkpoint(D/'proposal.json',result);print('CARRIER ANCHOR RELATION',result['status'],denominator,word,flush=True)
