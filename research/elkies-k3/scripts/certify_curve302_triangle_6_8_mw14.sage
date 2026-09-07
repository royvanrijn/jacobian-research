#!/usr/bin/env sage-python
"""Certify one degree-three A2+A1/MW14 pencil over QQ on the 11952 K3.

One fixed divisor, norm-two frame enumeration (cap 10000), 120 seconds.
Produces the rational pencil and abstract full MW lattice, not Weierstrass
section coordinates. Checks good reduction against each available modular
probe at the frozen primes 1009, 1013, 1021.
"""
import argparse
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import runpy
import signal
from sage.all import EllipticCurve,GF,PolynomialRing,QQ,ZZ,block_diagonal_matrix,matrix,pari,vector

ROOT=Path(__file__).resolve().parents[2]
PROBE=ROOT/'elkies-k3/scripts/probe_curve302_triangle_6_8_mw14.sage'
SOURCE=ROOT/'artifacts/generated-results/elkies-k3-r17-norm12-orbit11952-direct-fibration-v1.json'
OUT=ROOT/'artifacts/generated-results/elkies-k3-curve302-triangle-6-8-mw14-qq-v1.json'


def rows(A):return [list(map(str,r)) for r in A.rows()]
def rat(a):return {'numerator_coefficients_low_to_high':list(map(str,a.numerator().list())), 'denominator_coefficients_low_to_high':list(map(str,a.denominator().list()))}


def exact_roots(H):
    """Rational LDL enumeration, independent of PARI's short-vector routine."""
    n=H.nrows();L=matrix.identity(QQ,n);diag=[]
    for i in range(n):
        diag.append(QQ(H[i,i])-sum(L[i,k]**2*diag[k] for k in range(i)))
        assert diag[-1]>0
        for j in range(i+1,n):
            L[j,i]=(H[j,i]-sum(L[j,k]*L[i,k]*diag[k] for k in range(i)))/diag[i]
    assert L*matrix.diagonal(QQ,diag)*L.transpose()==H
    values=[0]*n;answer=[];nodes=0
    def visit(i,remaining):
        nonlocal nodes
        nodes+=1;assert nodes<=1000000
        if i<0:
            if remaining==0:answer.append(tuple(values))
            return
        center=-sum((L[j,i]*values[j] for j in range(i+1,n)),QQ(0))
        radius=isqrt(int((remaining/diag[i]).floor()))+1
        for v in range(int(center.floor())-radius,int(center.ceil())+radius+1):
            cost=diag[i]*(v-center)**2
            if cost<=remaining:
                values[i]=v;visit(i-1,remaining-cost)
    visit(n-1,QQ(2))
    return set(answer),nodes


def build():
    d=json.loads(SOURCE.read_text());G=matrix(ZZ,d['sections']['height_gram']);assert G.det()==948
    L=block_diagonal_matrix(matrix(ZZ,[[-2,1],[1,0]]),-G)
    e=matrix.identity(ZZ,17)
    def section(w):return vector(ZZ,[1,(w*G*w)//2]+list(w))
    O=section(vector(ZZ,17));P=section(e[5]);Q=section(e[7]);Z=section(e[6]);F=O+P+Q
    assert [O*L*P,O*L*Q,P*L*Q]==[1,1,1]
    assert F*L*F==0 and Z*L*Z==-2 and F*L*Z==1
    C=matrix(ZZ,[F*L,Z*L]).right_kernel_matrix();H=-C*L*C.transpose()
    U=matrix(ZZ,pari(H).qflllgram());assert abs(U.det())==1
    C=U.transpose()*C;H=-C*L*C.transpose();assert H.det()==948 and H.is_positive_definite()
    data=pari(H).qfminim(2,10000,2);V=matrix(ZZ,data[2]);assert int(data[0])==2*V.ncols()==8
    roots=[v for w in V.columns() for v in [w,-w]];assert all(v*H*v==2 for v in roots)
    exact,enumeration_nodes=exact_roots(H);assert exact=={tuple(v) for v in roots}
    pos=[v for v in roots if next(c for c in v if c)!=0 and next(c for c in v if c)>0]
    pset={tuple(v) for v in pos}
    simple=[v for v in pos if not any(tuple(v-w) in pset for w in pos)]
    B=matrix(ZZ,simple);RG=B*H*B.transpose();assert B.rank()==3 and RG.det()==6
    D,A,V=B.smith_form();assert list(D.diagonal())==[1,1,1]
    complement=V.inverse()[3:,:]
    projection=complement-complement*H*B.transpose()*RG.inverse()*B
    MW=projection*H*projection.transpose();assert MW.nrows()==14 and MW.det()==158
    mod=runpy.run_path(str(PROBE));K,f,coeff,zparam,c,meet=mod['make_pencil'](QQ)
    # Exact cancellation at the only three exceptional finite old fibres.
    R=PolynomialRing(QQ,'u');u=R.gen()
    def decode(v):return K(R(list(map(QQ,v['numerator_coefficients_low_to_high'])))/R(list(map(QQ,v['denominator_coefficients_low_to_high']))))
    E=EllipticCurve(K,[R(list(map(QQ,d['weierstrass_model'][name+'_coefficients_low_to_high']))) for name in ['A','B']])
    p,q=[E(decode(d['sections']['records'][i]['X']),decode(d['sections']['records'][i]['Y'])) for i in [5,7]]
    rp,rq,lam=meet
    assert p[0](lam)==q[0](lam) and p[1](lam)==q[1](lam)
    for point,r in [(p,rp),(q,rq)]:
        v=u-r
        assert point[0].valuation(v)==-2 and point[1].valuation(v)==-3
        a=-point[1]/point[0];principal=((v*a)(r))/v+(v*a).derivative()(r)
        assert (a-principal).valuation(v)>=1
        # At infinity these sections stay in the minimal old chart.
        assert point[0].numerator().degree()-point[0].denominator().degree()<=4
        assert point[1].numerator().degree()-point[1].denominator().degree()<=6
    assert c[0]+c[1] # There is a pole at O, as well as P and Q.
    witnesses=[]
    for prime in [1009,1013,1021]:
        path=OUT.with_name(f'elkies-k3-curve302-triangle-6-8-mw14-mod{prime}-v1.json')
        if not path.exists():continue
        witness=json.loads(path.read_text())
        assert witness['script_sha256']==sha256(PROBE.read_bytes()).hexdigest()
        field=GF(prime);Kp,fp,cp,zp,cc,mp=mod['make_pencil'](field)
        def reduce(a):return Kp(PolynomialRing(field,'u')([field(v) for v in a.numerator().list()])/PolynomialRing(field,'u')([field(v) for v in a.denominator().list()]))
        assert all(reduce(a)==b for a,b in zip(c,cc)) and reduce(zparam)==zp
        # Compare the actual rational polynomial before content normalization.
        assert all(all(reduce(f[i][j])==fp[i][j] for j in range(3)) for i in range(4))
        reduced=[[[int(t) for t in a.numerator().list()] for a in p.list()] for p in cp]
        assert reduced==witness['trigonal_coefficients_x_s_u_low_to_high']
        witnesses.append({'path':str(path.relative_to(ROOT)),'sha256':sha256(path.read_bytes()).hexdigest(),
                          'prime':prime,'j_degree':witness['j_degree'],'finite_target_roots':witness['finite_target_roots'],
                          'infinity_possible':witness['infinity_possible']})
    exclusions=[w['prime'] for w in witnesses if w['j_degree']==24 and not w['finite_target_roots'] and not w['infinity_possible']]
    return {'schema':'curve302.triangle-6-8-mw14-qq.v1','status':'PASS_QQ_MW14_AND_TARGET_EXCLUSION' if exclusions else 'PASS_QQ_MW14_TARGET_UNRESOLVED',
        'input_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),PROBE,SOURCE]},
        'source_authority':'EC-K3-R17-NORM12-11952-DIRECT-Q80-EQUATION',
        'limits':{'divisors':1,'root_norm':2,'root_storage_cap':10000,'seconds':120},
        'old_NS_basis':'O,F,v1,...,v17; Gram block diagonal ([[-2,1],[1,0]],-G)',
        'fibre':list(map(str,F)),'zero':list(map(str,Z)),'old_degree':3,
        'frame_basis':rows(C),'frame_gram':rows(H),'signed_root_count':8,'exact_LDL_enumeration_nodes':enumeration_nodes,'positive_root_vectors':rows(matrix(ZZ,pos)),
        'simple_root_vectors':rows(B),'simple_root_gram':rows(RG),'root_type':'A2+A1',
        'generic_arithmetic_MW_rank':14,'generic_torsion_order':1,'abstract_MW_gram':rows(MW),
        'abstract_MW_determinant':'158','MW_quotient_representatives_in_old_NS':rows(complement*C),
        'pencil_description':'s=c1*(y+yP6)/(x-xP6)+c2*(y+yP8)/(x-xP8)+c0',
        'pencil_coefficients':list(map(rat,c)),'zero_parameter':rat(zparam),'intersection_parameters':list(map(str,meet)),
        'trigonal_coefficients_x_s_u_low_to_high':[[list(map(str,a.numerator().list())) for a in p.list()] for p in coeff],
        'modular_witnesses':witnesses,'target_exclusion_primes':exclusions,
        'boundary':'This is an explicit rational genus-one pencil with a section and a complete abstract MW lattice. No characteristic-zero Weierstrass equation or fourteen explicit Weierstrass section coordinates are supplied. The root census is checked by both PARI and exact rational LDL enumeration; the remaining pencil and modular proof share implementation. A target exclusion, when present, applies only to this fibration. No302 parent, new surface isomorphism class, or universal MW14 exclusion is asserted.'}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--check',action='store_true');args=parser.parse_args();signal.alarm(120)
    result=build()
    if args.check:assert result==json.loads(OUT.read_text())
    else:OUT.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print(result['status'])
