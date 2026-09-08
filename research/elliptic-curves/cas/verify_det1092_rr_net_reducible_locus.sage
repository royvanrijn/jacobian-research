#!/usr/bin/env sage-python
"""Independent exact ball enumeration and divisor-word classification.

No import from the original CVP implementation. Enumeration uses a distinct
column permutation, outward integer bounds and a fixed closed radius10.
"""
import hashlib,json,math
from pathlib import Path
from sage.all import QQ,ZZ,PolynomialRing,EllipticCurve,matrix,vector
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
CERT=ART/'det1092_rr_net_reducible_locus_v2.json';OUT=ART/'det1092_rr_net_reducible_locus_replay_v2.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def verify():
    d=json.loads(CERT.read_text())
    for p,h in d['inputs'].items():assert sha(ROOT/p)==h
    pp=ART/'curve302_recovered_mw17_parent_v1.json';np=ART/'det1092_first_centre_rr_net_v1.json'
    op=ART/'det1092_historical_unlock_obstruction_v1.json'
    parent=json.loads(pp.read_text());net=json.loads(np.read_text());obs=json.loads(op.read_text())
    G=matrix(ZZ,parent['generic_height_gram']);w=vector(ZZ,d['centre_word']);assert w==-vector(ZZ,net['trace_word'])
    U=matrix(ZZ,d['LLL_columns']);assert abs(U.det())==1
    U.swap_columns(0,1);gram=U.transpose()*G*U;parity=vector(ZZ,U.inverse()*w)
    n=17;lower=matrix.identity(QQ,n);diagonal=[]
    for j in range(n):
        dj=gram[j,j]-sum(lower[j,k]**2*diagonal[k] for k in range(j));assert dj>0
        diagonal.append(dj)
        for i in range(j+1,n):lower[i,j]=(gram[i,j]-sum(lower[i,k]*lower[j,k]*diagonal[k] for k in range(j)))/dj
    assert lower*matrix.diagonal(QQ,diagonal)*lower.transpose()==gram
    words=[];z=[ZZ(0)]*n;nodes=0
    def visit(i,remaining):
        nonlocal nodes
        nodes+=1
        if nodes>200000:raise RuntimeError('independent exact node cap reached')
        if i<0:
            v=U*vector(ZZ,z);assert v*G*v==10;words.append(tuple(v));return
        shift=QQ(sum(lower[j,i]*z[j] for j in range(i+1,n)));radius2=QQ(remaining/diagonal[i])
        bound=ZZ(math.isqrt(int(radius2.numerator()//radius2.denominator())))+1
        centre=ZZ((-shift).floor());lo=centre-bound-1;hi=centre+bound+1
        for value in range(int(lo),int(hi)+1):
            if (value-parity[i])%2:continue
            cost=diagonal[i]*(value+shift)**2
            if cost>remaining:continue
            z[i]=ZZ(value);visit(i-1,remaining-cost)
    visit(16,QQ(10));words=set(words)
    assert words=={tuple(row) for row in d['exact_enumeration']['short_representatives']}
    pairs=set()
    for v in words:
        assert all((a-b)%2==0 for a,b in zip(w,v))
        x=vector(ZZ,[(a+b)//2 for a,b in zip(w,v)]);y=w-x
        assert x*G*x+y*G*y==10
        pairs.add(tuple(sorted([tuple(x),tuple(y)])))
    assert pairs=={tuple(tuple(v) for v in row['words']) for row in d['section_pairs']}
    assert len(pairs)==d['section_pair_count']==22 and len(words)==44
    for row in d['section_pairs']:
        assert row['heights']==[int(vector(ZZ,v)*G*vector(ZZ,v)) for v in row['words']]
    R=PolynomialRing(QQ,'t');K=R.fraction_field();t=R.gen()
    def dec(v):return K(R(v['numerator']))/R(v['denominator'])
    E=EllipticCurve(K,[dec(v) for v in parent['a_invariants']]);base=[E([dec(v) for v in row]) for row in parent['basis_weierstrass_coordinates']]
    T=-sum((v*P for v,P in zip(w,base)),E(0));A=list(map(R,net['A']));B=list(map(R,net['B']))
    canonical=d['canonical_O_plus_C_member'];V=list(map(R,canonical['line_coefficients']));cs=list(map(QQ,canonical['coefficients_in_A_tA_B']))
    assert V[2]==0 and V[1] and V[0]+V[1]*T[0]==0
    assert all(V[i]==cs[0]*A[i]+cs[1]*t*A[i]+cs[2]*B[i] for i in range(3))
    # Retain the old bisection obstruction, not the overbroad v1 statement.
    q0=QQ(obs['zero_q']);assert q0 and not q0.is_square()
    return {'classification':'verified application and new deduction','status':'PASS_INDEPENDENT_COMPLETE_FIXED_NET_REDUCIBLE_LOCUS',
        'short_representatives':len(words),'section_pairs':len(pairs),'independent_nodes':nodes,
        'complete_radius':10,'independent_column_permutation':[1,0,*range(2,17)],
        'canonical_O_plus_C_equation_verified':True,'old_bisection_nonsplit_at302':True,
        'claim_boundary':'22 known-section pairs and the old-bisection-plus-fibre pencil exhaust reducible divisors in this fixed net. Vertical components are not point constructions. No exclusion of useful old-bisection specializations on other fibres or of irreducible net members.',
        'v1_claim_boundary_rejected':True,'point_searches':0,
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [CERT,pp,np,op]},'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    d=verify();text=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if OUT.exists():assert OUT.read_text()==text
    else:OUT.write_text(text)
    print(d['status'],'nodes',d['independent_nodes'],flush=True)
