#!/usr/bin/env python3
"""Bounded generic height and Kummer saturation certificate at the retained prime."""
import argparse
from pathlib import Path
import retrospective as r
import residual_quartic_incidence as source

PROTOCOL=Path(__file__).with_name('RESIDUAL_QUARTIC_AT_131_PROTOCOL.json')
MOMENTS=r.OUT/'rank_jump_native_twist_frobenius_v1.json'
REPLAY=r.OUT/'rank_jump_native_twist_frobenius_verification_v2.json'
OUTPUT=r.OUT/'rank_jump_residual_quartic_at_131_v1.json'


def compute():
    from sage.all import QQ, GF, PolynomialRing, matrix
    data=r.read(source.INPUT); p=131
    def model(field):
        R=PolynomialRing(field,'t')
        convert=lambda cs:R([field(QQ(c)) for c in cs])
        A,B=convert(data['A']),convert(data['B'])
        sections=[(convert(s['x']),convert(s['y'])) for s in data['sections']]
        D=-4*A**3-27*B**2
        assert D.degree()==24 and D.is_squarefree() and D.gcd(A)==1
        assert 4*A[8]**3+27*B[12]**2 != 0
        for x,y in sections:assert y*y==x**3+A*x+B and x.degree()==4 and y.degree()==6
        gram=matrix(QQ,17,17); intersections=[]
        for i,(x,y) in enumerate(sections):
            gram[i,i]=4
            for j in range(i):
                xx,yy=sections[j];dx=x-xx;dy=y-yy
                assert dx and dy
                finite=int(dx.gcd(dy).degree())
                infinity=int(min(4-dx.degree(),6-dy.degree()))
                gram[i,j]=gram[j,i]=2-finite-infinity
                intersections.append([i,j,finite,infinity])
        assert gram.is_positive_definite()
        return A,B,D,sections,gram,intersections
    aq,bq,dq,sq,hq,iq=model(QQ)
    A,B,D,sections,hp,ip=model(GF(p))
    assert hq==hp and hq.det()==948
    rows=[];skipped=[];irreducible=[]
    R=PolynomialRing(GF(p),'X');X=R.gen()
    for t in range(p):
        if D(t)==0 or any(y(t)==0 for x,y in sections):
            skipped.append(t);continue
        f=X**3+A(t)*X+B(t)
        roots=f.roots(multiplicities=False)
        if not roots:irreducible.append(t)
        for theta in sorted(roots):
            values=[x(t)-theta for x,y in sections]
            assert all(values)
            bits=[int(not v.is_square()) for v in values]
            rows.append({'t':t,'theta':int(theta),'values':[int(v) for v in values], 'mask':r.pack(bits)})
    masks=[z['mask'] for z in rows];rank=r.rank(masks)
    witness=[];basis=[]
    for i,m in enumerate(masks):
        if r.rank(basis+[m])>len(basis):basis.append(m);witness.append(i)
    ledger=r.read(MOMENTS);replay=r.read(REPLAY)
    assert ledger['status']==replay['status']=='PASS'
    assert replay['bindings'][str(MOMENTS.relative_to(r.ROOT))]==r.digest(MOMENTS.read_bytes())
    traces=[-sum(row[3]*row[4] for row in ledger['fibre_trace_ledger'] if row[0]==n) for n in (1,2)]
    assert traces==[1884,319520] and traces==ledger['rows'][0]['Frobenius_traces']
    residual=[traces[n-1]-17*p**n for n in (1,2)]
    signs=[e for e in (-1,1) if (residual[0]-e*p)**2==residual[1]+p*p]
    assert signs==[-1] and residual==[-343,27783]
    a=residual[0]+p
    factor=QQ(2)*(2-QQ(a)/p)
    assert factor==QQ(948)/p and p*factor==hq.det()
    status='PASS' if rank==17 and irreducible else 'UNKNOWN'
    return {'schema':'rank-jump.residual-quartic-at-131.v1','status':status,'prime':p,
        'generic_sections':17,'gram':[[int(v) for v in row] for row in hq.rows()],
        'gram_determinant':948,'intersections_Q':iq,'intersections_F131':ip,
        'kummer_characters':rows,'kummer_character_rank':rank,'independent_character_row_indices':witness,
        'skipped_base_values':skipped,'irreducible_cubic_base_values':irreducible,
        'Frobenius_traces':traces,'residual_traces':residual,
        'residual_characteristic_polynomial_factors':[[p,1],[p*p,-a,1]],
        'normalized_Artin_Tate_factor':str(factor),'p_times_factor':int(p*factor),
        'finite_field_generic_rank':17,'two_saturation':'PROVED' if rank==17 else 'UNKNOWN',
        'bindings':source.bindings([Path(__file__),PROTOCOL,source.INPUT,MOMENTS,REPLAY,Path(r.__file__)]),
        'boundary':r.read(PROTOCOL)['boundary']}


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['build','check']);args=parser.parse_args()
    result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print(result['status'],'Kummer rank',result['kummer_character_rank'],'determinant',result['gram_determinant'])
