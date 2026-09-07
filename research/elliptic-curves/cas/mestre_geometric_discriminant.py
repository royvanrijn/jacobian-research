#!/usr/bin/env python3
"""Exact determinant accounting for the one missing geometric divisor line.

The proof uses the geometrically rational elliptic quotient, height scaling
under a degree-two base change, and integral involution eigenspace glue.
No lattice/frame enumeration or new section search is performed.
"""
import argparse
from fractions import Fraction as Q
from pathlib import Path
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
SOURCE=ART/'mestre_parent_picard_and_saturation_v1.json'
GEOMETRY=ART/'mestre_parent_fibre_geometry_v1.json'
OUT=ART/'mestre_geometric_discriminant_v1.json'

def compute():
    source=cert.read(SOURCE);geometry=cert.read(GEOMETRY)
    if source['status']!='PASS' or geometry['status']!='PASS':raise ArithmeticError('exact rank and geometry inputs required')
    rows=[]
    for r,g in zip(source['rows'],geometry['rows']):
        if r['outer_u']!=g['outer_u'] or r['geometric_NS_rank']!=19 or r['arithmetic_NS_rank']!=18 or r['certified_index']!=1:raise ArithmeticError('rank19/18 and rational saturation required')
        if g['quotient_geometric_mw_rank']!=6 or g['quotient_discriminant_profile']!=[[1,2],[8,1]] or r['arithmetic_invariant_rank']!=5:raise ArithmeticError('rank6 rational quotient and rank5 rational invariant subgroup required')
        # The rational quotient has NS determinant1, root determinant4,
        # no torsion since h>=2-(1/2+1/2)=1, hence MW determinant1/4.
        quotient_full_det=Q(1,4)
        quotient_fixed_det=Q(r['invariant_height_determinant'])/2**5
        if quotient_fixed_det!=Q(1,2):raise ArithmeticError('degree-two invariant height scaling differs')
        possibilities=[{'index':i,'anti_height':str(quotient_full_det*i*i/quotient_fixed_det)} for i in (1,2)]
        eligible=[a for a in possibilities if Q(a['anti_height'])>=1]
        if eligible!=[{'index':2,'anti_height':'2'}]:raise ArithmeticError('height lower bound did not determine anti line')
        # h(P)=2+2(P.O)-contr(P)=2 forces P.O=0 and contr=0,
        # since quotient correction is at most1. Thus P is narrow.
        quotient_anti_height=Q(2);pullback_anti_height=2*quotient_anti_height
        if pullback_anti_height!=4:raise ArithmeticError('pulled-back anti class norm differs')
        # The resulting anti NS class has square -4. It is primitive in
        # an even integral K3 lattice: division by2 would have square -1.
        # The quotient's index2 glue pulls back nontrivially, so the full
        # K3 NS eigenspace index is2, not1. Rank-one anti space bounds it by2.
        ns_index=2;geometric_det=Q(r['full_arithmetic_NS_absolute_determinant'])*pullback_anti_height/ns_index**2
        if geometric_det!=468:raise ArithmeticError('geometric discriminant differs')
        rows.append({'outer_u':r['outer_u'],'quotient_full_MW_determinant':str(quotient_full_det),
          'quotient_rational_MW_determinant':str(quotient_fixed_det),'quotient_anti_height_possibilities':possibilities,
          'quotient_anti_height':2,'quotient_MW_Galois_eigenspace_index':2,
          'primitive_K3_anti_NS_square':-4,'K3_NS_Galois_eigenspace_index':ns_index,
          'geometric_NS_rank':19,'geometric_NS_discriminant':-468,
          'arithmetic_NS_rank':18,'arithmetic_NS_discriminant':-468,
          'different_geometric_NS_from_determinant948':True,
          'all_Q_fibrations_MW_upper_bound':16,
          'full_geometric_NS_Gram':'UNKNOWN','transcendental_Gram':'UNKNOWN',
          'different_NS_arithmetic_MW17_foundry_admission':False})
    paths=[Path(__file__).resolve(),SOURCE,GEOMETRY]
    return {'schema':'elliptic-curves.mestre-geometric-discriminant.v1','status':'PASS','rows':rows,
      'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
      'proof':'The quotient geometric MW lattice has rank6 and determinant1/4; its Galois-fixed rank5 lattice has determinant1/2. An involution with rank-one anti space has eigenspace index1 or2. The quotient height lower bound1 forces anti height2 and index2. This generator is narrow. Its degree-two pullback has height4 and integral primitive anti NS class of square-4. The quotient index2 glue survives pullback, and the K3 anti space is still rank one, forcing NS eigenspace index2. Thus geometric NS determinant is468*4/4=468. The rational fixed lattice has rank18, so these determinant468 surfaces cannot supply an arithmetic MW17 fibration. Their full Gram and a useful new fibration remain unconstructed.',
      'scope':'A theorem-derived exact discriminant, not a full marked NS isometry or the arithmetic MW17 foundry milestone. Six Q-distinct parents now provably have geometric NS discriminant different from948. No new elliptic point, fibration or rootless frame is constructed.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--check',action='store_true');a=parser.parse_args();result=compute()
    if a.check:assert result==cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve geometric determinant proof')
        checkpoint(OUT,result)
    print('PASS6 GEOMETRIC NS DISC -468; ARITHMETIC RANK18; MW CAP16')
