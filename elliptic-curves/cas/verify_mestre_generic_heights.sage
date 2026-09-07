#!/usr/bin/env sage-python
"""Independent height replay by quadrupling into the narrow subgroup."""
import argparse,sys
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'elliptic-curves/cas'))
import certify_compact_r17_candidates as cert
from research_runtime.store import checkpoint
D=ROOT/'artifacts/local/elliptic-curves/mestre-generic-height-v1'
OUT=ROOT/'artifacts/generated-results/elliptic-curves/mestre_generic_height_independent_v1.json'

def compute():
    rows=[];paths=[Path(__file__).resolve(),D/'protocol.json',D/'ledger.json']
    if cert.read(D/'ledger.json')['status']!='PASS':raise ArithmeticError('all six original computations required')
    for u in cert.read(D/'protocol.json')['outer_parameters']:
        path=D/('u'+u)/'height.json';d=cert.read(path);paths.append(path)
        R=PolynomialRing(QQ,'T');T=R.gen();K=R.fraction_field();A=R(d['curve_A_coefficients']);B=R(d['curve_B_coefficients']);E=EllipticCurve(K,[A,B])
        P=[E([K(x),K(y)]) for x,y in d['covariant_points']]
        # Exponent4 kills all component groups (I2,I2,I4). The Shioda
        # height of4P is4+2(4P.O); no local component selector is reused.
        def h(point):
            if point.is_zero():return QQ(0)
            Q=4*point
            if Q.is_zero():return QQ(0)
            x=Q[0];n,den=x.numerator(),x.denominator()
            finite=QQ(den.degree())/2;infinite=max(QQ(0),QQ(n.degree()-den.degree()-4)/2)
            if finite.denominator()!=1 or infinite.denominator()!=1:raise ArithmeticError('integral zero-section intersection required')
            return (4+2*(finite+infinite))/16
        H=[h(point) for point in P];G=matrix(QQ,14,14)
        for i in range(14):
            G[i,i]=H[i]
            for j in range(i):G[i,j]=G[j,i]=(h(P[i]+P[j])-H[i]-H[j])/2
        if G!=matrix(QQ,d['covariant_height_gram']):raise ArithmeticError('quadrupling-degree height Gram differs')
        for word in d['generic_covariant_relations']:
            if sum((ZZ(a)*point for a,point in zip(word,P)),E(0))!=E(0):raise ArithmeticError('generic relation replay failed')
        C=matrix(QQ,d['divisor_change_matrix']);H=C*G*C.transpose();basis=H.matrix_from_rows_and_columns(d['seed_indices'],d['seed_indices'])
        if H!=matrix(QQ,d['divisor_height_gram']) or basis!=matrix(QQ,d['seed_height_gram']) or not basis.is_positive_definite() or str(16*basis.det())!=d['known_rank18_divisor_lattice_absolute_determinant']:raise ArithmeticError('divisor lattice determinant differs')
        rows.append({'outer_u':u,'quadrupled_height_checks':119,'generic_relations':3,'section_rank':int(G.rank()),'height_determinant':str(basis.det()),'rank18_divisor_sublattice_absolute_determinant':str(16*basis.det())})
        print('INDEPENDENT u'+u,'119 HEIGHTS; DET',16*basis.det(),flush=True)
    return {'schema':'elliptic-curves.mestre-generic-height-independent.v1','status':'PASS','rows':rows,
      'sources':{str(p.relative_to(ROOT)):cert.hashed(p) for p in paths},
      'scope':'Point membership and exact rational-function group law on saved section equations. All119 heights per parent are recomputed from degrees after multiplication by4 clears the I2/I2/I4 component groups. No producer local component selector is imported. Three generic group relations per parent and every Gram/determinant replay. No full NS or MW saturation or ambient-rank upper bound.'}

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');a=p.parse_args();result=compute()
    if a.check:assert result==cert.read(OUT)
    else:
        if OUT.exists():raise FileExistsError('preserve independent height proof')
        checkpoint(OUT,result)
