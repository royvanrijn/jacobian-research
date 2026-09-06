#!/usr/bin/env python3
"""Explicit genus-one 2-covers and symplectic pairs for the minus reference."""
import argparse
from pathlib import Path
import retrospective as r
import fixed_cubic_minus_reference as source
import bounded_gain_reference_completion as completed

OUTPUT=r.OUT/'rank_jump_minus_reference_carriers_v1.json'


def compute():
    from sage.all import QQ,PolynomialRing,matrix,GF,EllipticCurve
    data=r.read(source.OUTPUT);loc=r.read(completed.OUTPUT)['stages']['local'];M=data['twist_CT_matrix']
    def pairing(a,b):return sum(((a>>i)&1)*M[i][j]*((b>>j)&1) for i in range(6) for j in range(6))%2
    remaining=[1<<i for i in range(6)];pairs=[]
    while remaining:
        a=remaining[0];b=next(v for v in remaining[1:] if pairing(a,v));pairs.append([a,b])
        remaining=[v^(a if pairing(v,b) else 0)^(b if pairing(v,a) else 0) for v in remaining if v not in (a,b)]
    change=[v for pair in pairs for v in pair];T=matrix(GF(2),[[v>>i&1 for i in range(6)] for v in change])
    J=[[int(i//2==j//2 and i!=j) for j in range(6)] for i in range(6)]
    assert T.rank()==6 and T*matrix(GF(2),M)*T.transpose()==matrix(GF(2),J)
    R=PolynomialRing(QQ,'z');f=R(data['arms'][1]['cubic_ascending']);z=R.gen()
    P=PolynomialRing(QQ,names=('u0','u1','u2','v'));u0,u1,u2,v=P.gens()
    Q=PolynomialRing(P,'zz');zz=Q.gen();h=Q(f.list());u=u0+u1*zz+u2*zz**2
    L=PolynomialRing(QQ,'t');t=L.gen();rows=[]
    E=EllipticCurve(list(map(QQ,data['arms'][1]['model'])))
    for i,rec in enumerate(loc['class_records']):
        beta=[QQ(c)*(-1)**j for j,c in enumerate(rec['beta_ascending'])]
        product=(Q(beta)*u*u)%h;quadrics=[product[1]+v*v,product[2]];matrices=[]
        for q in quadrics:
            H=matrix(QQ,4)
            for j,x in enumerate(P.gens()):
                for k,y in enumerate(P.gens()):H[j,k]=q.monomial_coefficient(x*y)/(1 if j==k else 2)
            assert sum(H[j,k]*P.gen(j)*P.gen(k) for j in range(4) for k in range(4))==q
            matrices.append(H)
        pencil=(matrices[0].change_ring(L)+t*matrices[1].change_ring(L)).det()
        assert pencil.degree()==3 and pencil.discriminant()!=0
        aa,bb,cc,dd,ee=0,pencil[3],pencil[2],pencil[1],pencil[0]
        I=12*aa*ee-3*bb*dd+cc**2
        JJ=72*aa*cc*ee+9*bb*cc*dd-27*aa*dd**2-27*bb**2*ee-2*cc**3
        assert 6912*I**3/(4*I**3-JJ**2)==E.j_invariant()
        rows.append({'index':i,'generic_mask':rec['generic_mask'],'beta_minus_ascending':list(map(str,beta)),
            'quadric_matrices':[[list(map(str,row)) for row in H.rows()] for H in matrices],
            'pencil_ascending':list(map(str,pencil.list())),'pencil_discriminant':str(pencil.discriminant()),
            'genus':1,'projective_degree':4,'locally_soluble':True,'rationally_soluble':False})
    files=(Path(__file__),source.OUTPUT,completed.OUTPUT)
    return {'schema':'rank-jump.minus-reference-carriers.v1','status':'PASS','rows':rows,
        'symplectic_pairs_in_strict_basis':pairs,'standard_CT_matrix':J,
        'simultaneous_carrier':'Product of the six genus-one torsors over Q, a torsor under the sixth power of the twist elliptic curve.',
        'minimal_abelian_torsor_carrier_dimension':6,'common_smooth_curve_carrier_genus_lower_bound':6,
        'minimality_scope':'Carriers with Q-morphisms to all six torsors. Abelian-torsor minimality follows from six independent H1(Q,E)[2] classes and End_Q(E)=Z; see proof note. Exact minimal curve genus is not computed.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();data=compute()
    if args.mode=='build':r.write_new(OUTPUT,data)
    else:assert data==r.read(OUTPUT)
    print('PASS six genus-one carriers',data['symplectic_pairs_in_strict_basis'],flush=True)
