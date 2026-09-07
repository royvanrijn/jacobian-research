#!/usr/bin/env python3
"""Independent sign-module and local-power verification; no point search."""
import argparse
from collections import deque
from itertools import product
from pathlib import Path
import retrospective as r
import four_division_class_creation as run

OUTPUT=r.OUT/'rank_jump_four_division_class_creation_verification_v1.json'


def module_check(saved):
    # Enumerate all subspaces as sets, independently of binary row reduction.
    zero=(0,0,0,0);space=list(product((0,1),repeat=4))
    plus=lambda a,b:tuple(x^y for x,y in zip(a,b))
    def S(v):z,a,b,c=v;return z,z^a,c,b
    def T(v):z,a,b,c=v;return z,c,z^a,z^b
    def phi(v):z,a,b,c=v;return z^a^b,a^c,z^b^c
    subs={frozenset([zero])};pending=[frozenset([zero])]
    while pending:
        h=pending.pop()
        for v in space:
            k=h|frozenset(plus(w,v) for w in h)
            if k not in subs:subs.add(k);pending.append(k)
    assert len(subs)==67
    stable=[h for h in subs if {S(v) for v in h}==set(h) and {T(v) for v in h}==set(h)]
    rows=[]
    for h in stable:
        fixed={v for v in h if T(v)==v};standard={plus(v,T(v)) for v in h}
        assert len(h)==len(fixed)*len(standard)
        assert all(phi(v)==(0,0,0) for v in fixed)
        assert len(standard) in (1,4)
        assert len({phi(v) for v in standard})==len(standard)
        rows.append({'subspace':sorted(r.pack(v) for v in h),'standard_dimension':0 if len(standard)==1 else 2})
    assert saved['sign_action_transposition']==[r.pack(S(tuple(v>>i&1 for i in range(4)))) for v in range(16)]
    assert saved['sign_action_three_cycle']==[r.pack(T(tuple(v>>i&1 for i in range(4)))) for v in range(16)]
    assert saved['derivative_map']==[r.pack(phi(tuple(v>>i&1 for i in range(4)))) for v in range(16)]
    return {'all_subspaces':67,'stable_subspaces':sorted(rows,key=lambda x:x['subspace']),
        'conclusion':'Every stable subgroup has zero or one standard summand, detected exactly by the derivative map.'}


def group_check(row):
    from sage.all import GF, matrix
    n=row['modulus'];gens=list(map(tuple,row['generators']));I=(1,0,0,1)
    def mul(a,b):return ((a[0]*b[0]+a[1]*b[2])%n,(a[0]*b[1]+a[1]*b[3])%n,
                         (a[2]*b[0]+a[3]*b[2])%n,(a[2]*b[1]+a[3]*b[3])%n)
    # Sparse equation set followed by Sage GF(2) kernel, instead of packed-rank elimination.
    expressions={I:(0,0)};queue=deque([I]);equations=set()
    while queue:
        a=queue.popleft()
        for j,b in enumerate(gens):
            c=mul(a,b);u=expressions[a]
            v=(u[0]^((a[0]%2)<<(2*j))^((a[1]%2)<<(2*j+1)),
               u[1]^((a[2]%2)<<(2*j))^((a[3]%2)<<(2*j+1)))
            if c in expressions:equations.update(x^y for x,y in zip(expressions[c],v))
            else:expressions[c]=v;queue.append(c)
    assert len(expressions)==row['group_order']==6*(n//2)**4
    M=matrix(GF(2),[[v>>j&1 for j in range(2*len(gens))] for v in equations])
    K=M.right_kernel();assignments=sorted(r.pack(v) for v in K)
    assert assignments==row['cocycle_generator_assignments'] and K.dimension()==3
    cob=[]
    for x,y in product((0,1),repeat=2):
        vals=[]
        for a,b,c,d in gens:vals.extend(((a*x+b*y+x)%2,(c*x+d*y+y)%2))
        cob.append(r.pack(vals))
    assert len(set(cob))==4 and set(cob)<=set(assignments)
    return {'modulus':n,'order':len(expressions),'Z1_dimension':3,'B1_dimension':2,'H1_dimension':1}


def verify():
    from sage.all import QQ, ZZ, PolynomialRing, EllipticCurve, lcm, GF, pari
    pari.allocatemem(64000000,256000000,silent=True)
    out=r.read(run.OUTPUT);assert out['status']=='PASS'
    for name,sha in out['bindings'].items():assert r.digest((r.ROOT/name).read_bytes())==sha,name
    raw={x['token']:x for x in r.read(run.panel.INPUT)['cases']}
    gal={x['token']:x['galois'] for x in r.read(run.octics.OUTPUT)['rows']}
    R=PolynomialRing(QQ,'z');records=[]
    for row in out['rows']:
        token=row['token'];model=list(map(QQ,raw[token]['model']));E=EllipticCurve(model)
        A=-E.c4()/48;B=-E.c6()/864;scale=lcm(A.denominator(),B.denominator())
        f=R([B*scale**6,A*scale**4,0,1]);assert list(map(str,f.list()))==row['cubic_ascending']
        delta=f.discriminant();assert str(delta)==row['discriminant'] and not ZZ(delta).is_square()
        p=gal[token]['irreducibility_prime'];assert ZZ(p).is_prime(proof=True)
        fp=f.change_ring(GF(p));assert fp.degree()==3 and fp.is_irreducible()
        beta=R(list(map(QQ,row['derivative_ascending'])))
        assert beta==-delta*f.derivative() and f.resultant(beta)==delta**4
        place=row['witness']['place'];record={'token':token,'irreducibility_prime':p,'witness_place':place}
        if place=='infinity':
            assert delta>0 and row['witness']['derivative_signs']==[1,0,1]
            assert [1,0,1] not in row['witness']['full_point_image']
        else:
            # Independent PARI local-power oracle, not LocalSquareclasses signatures.
            nf=pari.nfinit([pari(f),[place]]);primes=pari.idealprimedec(nf,place)
            theta=pari.Mod('z',pari(f));b=pari.Mod(pari(beta),pari(f))
            def local_square(a):
                # PARI's odd-prime local-power path can reject a nonintegral
                # basis vector. Clear its denominator by a rational square.
                coordinates=pari.nfalgtobasis(nf,a)
                denominator=lcm([QQ(c).denominator() for c in coordinates])
                integral=a*pari(denominator)**2
                return all(pari.nfislocalpower(nf,P,integral,2) for P in primes)
            representatives=[pari.Mod(1,pari(f))];indices=[]
            d=len(primes)-1+int(place==2)
            for i,(x,y) in enumerate(raw[token]['generic_sections']):
                x,y=QQ(x),QQ(y);assert E.is_on_curve(x,y)
                shifted=(x+E.b2()/12)*scale**2;gamma=pari(shifted)-theta
                if not any(local_square(gamma/a) for a in representatives):
                    representatives+= [gamma*a for a in representatives];indices.append(i)
                if len(representatives)==2**d:break
            assert len(representatives)==2**d
            assert all(not local_square(b/a) for a in representatives)
            record.update(point_dimension=d,generic_basis_indices=indices,
                          excluded_local_point_classes=len(representatives))
        assert row['four_division_class_dimension']==1
        assert row['four_division_Selmer_dimension']==row['four_division_strict_dimension']==0
        records.append(record);print(token,'PASS',place,flush=True)
    module=module_check(out['module']);groups=[group_check(x) for x in out['abstract_full_group_controls']]
    return {'schema':'rank-jump.four-division-class-creation-verification.v1','status':'PASS',
        'module':module,'abstract_full_group_controls':groups,'rows':records,
        'method':'All 67 sign-module subspaces; independent Weierstrass invariants, finite cubic irreducibility, exact norms and PARI local-power tests against complete generic local images; Sage GF(2) cohomology kernels.',
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),run.OUTPUT,run.panel.INPUT,run.octics.OUTPUT,Path(r.__file__)]}}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=verify()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert result==r.read(OUTPUT)
    print('PASS independent sixteen-fibre exclusion and sign-module proof')
