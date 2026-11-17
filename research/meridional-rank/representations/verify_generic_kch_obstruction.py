#!/usr/bin/env python3
"""Signed generic-mu rank-one obstruction for13a373; standard library only.

Replays exact integer Laurent identities, literal oriented diagram/seed
binding, a rank-three quadratic-field cord witness, and a genuine generic
T(4,5) irreducible rank-four control. --census additionally re-reads the
pinned independently packaged named diagram. No finite sampling is used.
"""
import argparse
from copy import deepcopy
from fractions import Fraction
from itertools import permutations
import hashlib
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
INPUT=HERE/'generic_kch_obstruction_inputs.json'


class Ring:
    """Sparse integer Laurent polynomials; all negative exponents explicit."""
    def __init__(self,names):
        self.names=names.split();self.zero=(0,)*len(self.names)
    def c(self,n):return {self.zero:n} if n else {}
    def v(self,name,power=1):
        e=list(self.zero);e[self.names.index(name)]=power
        return {tuple(e):1}
    def add(self,*ps):
        out={}
        for p in ps:
            for e,a in p.items():
                out[e]=out.get(e,0)+a
                if not out[e]:del out[e]
        return out
    def neg(self,p):return {e:-a for e,a in p.items()}
    def sub(self,p,q):return self.add(p,self.neg(q))
    def mul(self,*ps):
        out=self.c(1)
        for p in ps:
            tmp={}
            for e,a in out.items():
                for f,b in p.items():
                    g=tuple(i+j for i,j in zip(e,f))
                    tmp[g]=tmp.get(g,0)+a*b
                    if not tmp[g]:del tmp[g]
            out=tmp
        return out
    def power(self,p,n):
        assert n>=0
        out=self.c(1)
        while n:
            if n&1:out=self.mul(out,p)
            p=self.mul(p,p);n//=2
        return out
    def substitute(self,p,images):
        out={}
        for e,a in p.items():
            term=self.c(a)
            for name,n in zip(self.names,e):
                if name in images:
                    assert n>=0
                    factor=self.power(images[name],n)
                else:factor=self.v(name,n)
                term=self.mul(term,factor)
            out=self.add(out,term)
        return out
    def dump(self,p):return [[list(e),a] for e,a in sorted(p.items())]


def basis(R,n):return [[R.c(int(i==j)) for j in range(n)] for i in range(n)]
def vadd(R,u,v):return [R.add(a,b) for a,b in zip(u,v)]
def vscale(R,a,u):return [R.mul(a,b) for b in u]
def vsub(R,u,v):return vadd(R,u,vscale(R,R.c(-1),v))
def dot(R,B,f,u):
    return R.add(*(R.mul(f[i],B[i][j],u[j]) for i in range(len(f)) for j in range(len(u))))


def propagate(R,B,roots,forms,step):
    o,a,b,exponent=step
    u,f=roots[a],forms[a];v,g=roots[o],forms[o]
    gu=dot(R,B,g,u);fv=dot(R,B,f,v)
    if exponent==1:
        roots[b]=vsub(R,u,vscale(R,gu,v))
        forms[b]=vadd(R,f,vscale(R,R.mul(fv,R.v('mu',-1)),g))
    else:
        assert exponent==-1
        roots[b]=vadd(R,u,vscale(R,R.mul(gu,R.v('mu',-1)),v))
        forms[b]=vsub(R,f,vscale(R,fv,g))
    assert dot(R,B,forms[b],roots[b])==R.sub(R.c(1),R.v('mu'))


def raw_first_loop(data):
    R=Ring('mu r R x X a b')
    mu,r,rr,x,X,a,b=[R.v(n) for n in R.names]
    c,add,sub,mul=R.c,R.add,R.sub,R.mul;mi=R.v('mu',-1);d=sub(c(1),mu)
    B=[[d,r,x],[rr,d,a],[X,b,d]]
    e=basis(R,3);roots={k:e[i] for i,k in enumerate([1,2,6])}
    forms=deepcopy(roots)
    assert data['signed_trace'][:4]==[[2,6,7,1],[7,2,3,1],[3,1,0,1],[1,3,4,-1]]
    for step in data['signed_trace'][:3]:propagate(R,B,roots,forms,step)
    k=add(c(1),mul(a,b,mi))
    ss=add(mul(k,rr),mul(a,X));tt=sub(mul(k,r),mul(b,x,mi))
    assert roots[3]==[{},k,R.neg(mul(b,mi))] and forms[3]==[{},k,a]
    assert roots[0]==vsub(R,e[0],vscale(R,ss,roots[3]))
    assert forms[0]==vadd(R,e[0],vscale(R,mul(tt,mi),forms[3]))
    H=add(d,mul(ss,tt,mi))
    assert dot(R,B,forms[0],e[0])==H
    closing_root=vadd(R,e[0],vscale(R,mul(H,mi),roots[0]))
    assert closing_root==[add(c(1),mul(H,mi)),R.neg(mul(H,ss,k,mi)),mul(H,ss,b,mi,mi)]
    a_u0=sub(d,mul(ss,tt))
    closing_form=vsub(R,e[0],vscale(R,a_u0,forms[0]))
    assert closing_form==vsub(R,e[0],vscale(R,a_u0,vadd(R,e[0],vscale(R,mul(tt,mi),forms[3]))))
    # The written coefficient proof forces a=b=0, r*R=-mu. Replay the
    # admissible normal form R=1,r=-mu and both pairs on0/4/6/7.
    normalized=[[d,R.neg(mu),x],[c(1),d,{}],[X,{},d]]
    roots={k:e[i] for i,k in enumerate([1,2,6])};forms=deepcopy(roots)
    for step in data['signed_trace'][:4]:propagate(R,normalized,roots,forms,step)
    U=vsub(R,e[0],e[1])
    assert roots[0]==forms[0]==U
    assert roots[4]==forms[4]==vscale(R,c(-1),U)
    assert roots[7]==forms[7]==e[2]
    return {'signed_first_clasp_steps':3,'signed_normalized_steps':4,
            'constraints':['A12=0','A21=0','A01*A10=-mu'],
            'normalized_A10':1,'normalized_A01':'-mu',
            'seed_root_and_covector_basis_changes_determinant':1}


def paired_core(data):
    R=Ring('mu x X y Y z Z');mu,x,X,y,Y,z,Z=[R.v(n) for n in R.names]
    c,add,sub,mul,neg=R.c,R.add,R.sub,R.mul,R.neg
    mi=R.v('mu',-1);d=sub(c(1),mu)
    B=[[d,x,y],[X,d,z],[Y,Z,d]];U,V,W=basis(R,3)
    roots={0:U,4:vscale(R,c(-1),U),6:V,7:V,8:W};forms=deepcopy(roots)
    assert data['signed_trace'][4:]==[[8,4,5,-1],[5,8,9,-1],[6,9,10,1],[10,0,12,1],[4,10,11,-1]]
    for step in data['signed_trace'][4:]:propagate(R,B,roots,forms,step)
    h=add(c(1),mul(y,Y,mi));v=add(mul(y,X),mul(h,z));C=sub(mul(h,Z),mul(Y,x,mi))
    w=sub(mul(y,add(d,h)),mul(x,v));Q=add(mul(Y,sub(h,mul(d,mi))),mul(C,X,mi))
    expected_u={5:[c(-1),{},neg(mul(Y,mi))],9:[y,{},h],10:[y,neg(v),h],
                11:[add(y,mul(w,mi)),neg(v),h],12:[sub(c(1),mul(Q,y)),mul(Q,v),neg(mul(Q,h))]}
    expected_f={5:[c(-1),{},y],9:[neg(mul(Y,mi)),{},h],10:[neg(mul(Y,mi)),mul(C,mi),h],
                11:[sub(neg(mul(Y,mi)),Q),mul(C,mi),h],
                12:[sub(c(1),mul(w,Y,mi,mi)),mul(w,C,mi,mi),mul(w,h,mi)]}
    assert all(roots[i]==u for i,u in expected_u.items())
    assert all(forms[i]==f for i,f in expected_f.items())
    return {'generic_paired_core_steps':5,'independent_directed_pairings':6,
            'paired_formula_rows_checked':10,'all_propagated_norms':'1-mu'}


def unit_certificate(data):
    R=Ring('mu h y Y v C w Q p q r s alpha beta')
    M,H,y,Y,v,C,w,Q,p,q,r,s,alpha,beta=[R.v(n) for n in R.names]
    mi,ai,bi=R.v('mu',-1),R.v('alpha',-1),R.v('beta',-1)
    c,add,sub,mul,neg=R.c,R.add,R.sub,R.mul,R.neg
    # First closing equation: compare in independent U,V,W and F,G,J.
    E=basis(R,3);U,V,W=E
    u11=[add(y,mul(w,mi)),neg(v),H]
    f11=[sub(neg(mul(Y,mi)),Q),mul(C,mi),H]
    root11=vsub(R,vsub(R,V,vscale(R,p,u11)),vscale(R,beta,W))
    form11=vsub(R,vadd(R,V,vscale(R,mul(q,mi),f11)),vscale(R,bi,W))
    assert root11==[neg(mul(p,add(y,mul(w,mi)))),add(c(1),mul(p,v)),sub(neg(mul(p,H)),beta)]
    assert form11==[mul(q,mi,sub(neg(mul(Y,mi)),Q)),add(c(1),mul(q,C,mi,mi)),sub(mul(q,H,mi),bi)]
    # Nonzero p,q follow from their V/G coefficients. Thus w=-mu*y,
    # Q=-Y/mu. Verify the exact cleaned pair before using its coefficients.
    u5=[c(-1),{},neg(mul(Y,mi))];f5=[c(-1),{},y]
    u12=[sub(c(1),mul(Q,y)),mul(Q,v),neg(mul(Q,H))]
    f12=[sub(c(1),mul(w,Y,mi,mi)),mul(w,C,mi,mi),mul(w,H,mi)]
    images={'w':neg(mul(M,y)),'Q':neg(mul(Y,mi)),'h':add(c(1),mul(y,Y,mi))}
    diff_u=vadd(R,vadd(R,u12,vscale(R,H,u5)),vscale(R,mul(Y,v,mi),V))
    diff_f=vadd(R,vadd(R,f12,vscale(R,H,f5)),vscale(R,mul(y,C,mi),V))
    assert all(not R.substitute(a,images) for a in diff_u+diff_f)
    clean_u12=vsub(R,vscale(R,neg(H),u5),vscale(R,mul(Y,v,mi),V))
    clean_f12=vsub(R,vscale(R,neg(H),f5),vscale(R,mul(y,C,mi),V))
    root12=vsub(R,vadd(R,u5,vscale(R,mul(r,mi),clean_u12)),vscale(R,alpha,V))
    form12=vsub(R,vsub(R,f5,vscale(R,s,clean_f12)),vscale(R,ai,V))
    assert root12==vadd(R,vscale(R,sub(c(1),mul(r,H,mi)),u5),
                       vscale(R,sub(neg(mul(r,Y,v,mi,mi)),alpha),V))
    assert form12==vadd(R,vscale(R,add(c(1),mul(s,H)),f5),
                       vscale(R,sub(mul(s,y,C,mi),ai),V))
    F=[sub(sub(H,c(1)),mul(y,Y,mi)),
       mul(p,add(y,mul(w,mi))),add(c(1),mul(p,v)),sub(neg(mul(p,H)),beta),
       mul(q,mi,sub(neg(mul(Y,mi)),Q)),add(c(1),mul(q,C,mi,mi)),sub(mul(q,H,mi),bi),
       sub(c(1),mul(r,H,mi)),sub(neg(mul(r,Y,v,mi,mi)),alpha),
       add(c(1),mul(s,H)),sub(mul(s,y,C,mi),ai)]
    def generator(i):return F[i],[c(int(i==j)) for j in range(len(F))]
    def scale(a,e):return mul(a,e[0]),[mul(a,b) for b in e[1]]
    def total(*es):return add(*(e[0] for e in es)),[add(*(e[1][j] for e in es)) for j in range(len(F))]
    f=list(map(generator,range(len(F))))
    ev=total(scale(neg(mul(bi,H)),f[2]),scale(neg(mul(bi,v)),f[3]))
    assert ev[0]==sub(v,mul(bi,H))
    ec=total(scale(neg(mul(M,beta,H)),f[5]),scale(mul(beta,C),f[6]))
    assert ec[0]==sub(neg(C),mul(M,beta,H))
    ea=total(f[8],scale(neg(mul(Y,bi,mi)),f[7]),scale(mul(r,Y,mi,mi),ev))
    assert ea[0]==sub(neg(alpha),mul(Y,bi,mi))
    eai=total(f[10],scale(mul(y,beta),f[9]),scale(mul(s,y,mi),ec))
    assert eai[0]==sub(mul(y,beta),ai)
    ey=total(scale(neg(ai),ea),scale(mul(Y,bi,mi),eai))
    assert ey[0]==add(c(1),mul(y,Y,mi))
    eh=total(f[0],ey);assert eh[0]==H
    one=total(f[7],scale(mul(r,mi),eh));assert one[0]==c(1)
    out={'variables':R.names,'inverted_variables':['mu','alpha','beta'],
         'equations':[R.dump(a) for a in F],'multipliers':[R.dump(a) for a in one[1]],'identity_constant':1}
    assert out==data['unit_identity']
    decode=lambda p:{tuple(e):a for e,a in p}
    assert add(*(mul(decode(a),decode(b)) for a,b in zip(out['multipliers'],out['equations'])))==c(1)
    mod2=lambda p:{e:a%2 for e,a in p.items() if a%2}
    assert mod2(add(*(mul(mod2(a),mod2(b)) for a,b in zip(one[1],F))))==c(1)
    bad=deepcopy(one[1]);bad[0][R.zero]=bad[0].get(R.zero,0)+1
    assert add(*(mul(a,b) for a,b in zip(bad,F)))!=c(1)
    return {'integer_Laurent_identity':1,'equations':len(F),
            'multiplier_monomials':sum(map(len,one[1])),'identity_mod2':1,
            'closing_coefficient_binding_checked':True,'corrupted_multiplier_rejected':True}


def diagram(pd):
    ends={};adj={e:set() for row in pd for e in row}
    for i,row in enumerate(pd):
        assert len(row)==4 and len(set(row))==4
        for j,e in enumerate(row):ends.setdefault(e,[]).append(4*i+j)
        a,b,c,d=row;adj[b].add(d);adj[d].add(b)
    assert len(ends)==2*len(pd) and all(len(v)==2 for v in ends.values())
    alpha={a:b for pair in ends.values() for a,b in (pair,pair[::-1])}
    unused=set(alpha);faces=0
    while unused:
        start=x=min(unused)
        while True:
            unused.remove(x);y=alpha[x];x=4*(y//4)+(y+1)%4
            if x==start:break
        faces+=1
    assert faces==len(pd)+2
    # Coherent orientation from original dart(0,0); it crosses each crossing
    # twice and assigns exactly one incoming corner on each strand.
    incoming={};x=0;seen=set()
    while x not in seen:
        seen.add(x);i,j=divmod(x,4);incoming.setdefault(i,[]).append(j)
        x=alpha[4*i+(j+2)%4]
    assert x==0 and len(seen)==2*len(pd) and all(len(v)==2 for v in incoming.values())
    components=[];unused=set(adj)
    while unused:
        cc=set();queue=[min(unused)]
        for e in queue:
            if e in cc:continue
            cc.add(e);queue.extend(adj[e]-cc)
        components.append(sorted(cc));unused-=cc
    components.sort();assert len(components)==len(pd)
    arc={e:i for i,cc in enumerate(components) for e in cc}
    signed=[]
    for i,row in enumerate(pd):
        ui=next(j for j in incoming[i] if j%2==0)
        oi=next(j for j in incoming[i] if j%2==1)
        sign=1 if (oi-ui)%4==3 else -1
        signed.append({'crossing':i,'over':arc[row[1]],'input':arc[row[ui]],
                       'output':arc[row[(ui+2)%4]],'exponent':-sign,'crossing_sign':sign})
    return signed,alpha


def census_bijection(pd,census_pd,cert):
    _,a=diagram(pd);_,b=diagram(census_pd)
    f=cert['dart_images'];orientation=cert['orientation'];swap=cert['over_under_swap']
    assert orientation in (-1,1) and swap in (0,1)
    assert len(pd)==len(census_pd) and sorted(f)==list(range(4*len(pd)))
    for x,y in enumerate(f):
        assert f[a[x]]==b[y]
        assert (x%2)^(y%2)==swap
        assert f[4*(x//4)+(x+1)%4]==4*(y//4)+(y+orientation)%4


def case_binding(data):
    case=data['case'];assert case['name']=='13a_373'
    source=json.loads((HERE/'seven_sheet_cover_inputs.json').read_text())
    old=next(c for c in source['cases'] if c['name']==case['name'])
    assert all(case[k]==old[k] for k in case)
    signed,_=diagram(case['pd']);assert signed==data['signed_crossings']
    census_bijection(case['pd'],case['census_pd'],case['census_map_witness'])
    known=set(case['seed_arcs']);assert case['seed_arcs']==[1,2,6,8]
    used=set()
    for o,a,b,e in data['signed_trace']:
        assert o in known and a in known and b not in known and e in (-1,1)
        hits=[r for r in signed if r['over']==o and {r['input'],r['output']}=={a,b}]
        assert len(hits)==1
        r=hits[0];assert e==(r['exponent'] if a==r['input'] else -r['exponent'])
        known.add(b);used.add(r['crossing'])
    assert known==set(range(13)) and len(used)==9
    left=[[r['over'],r['input'],r['output'],r['exponent']] for r in signed if r['crossing'] not in used]
    assert left==data['remaining_signed_crossings']
    assert [0,1,2,-1] in left and [11,7,8,1] in left and [12,5,6,-1] in left
    return {'literal_PD_crossings':13,'signed_crossings_checked':13,
            'four_seed_signed_trace_steps':9,'remaining_signed_crossings':4,
            'exact_named_planar_map_checked':True}


# Exact Q(phi), phi^2=phi+1, only needed for the separate rank3 attainment.
Q0=(Fraction(0),Fraction(0));Q1=(Fraction(1),Fraction(0))
def qa(a,b):return a[0]+b[0],a[1]+b[1]
def qn(a):return -a[0],-a[1]
def qm(a,b):return a[0]*b[0]+a[1]*b[1],a[0]*b[1]+a[1]*b[0]+a[1]*b[1]
def qi(a):
    norm=a[0]*a[0]+a[0]*a[1]-a[1]*a[1];assert norm
    return (a[0]+a[1])/norm,-a[1]/norm
def qsum(xs):
    out=Q0
    for a in xs:out=qa(out,a)
    return out
def qrank(matrix):
    A=deepcopy(matrix);row=0
    for j in range(len(A[0])):
        k=next((k for k in range(row,len(A)) if A[k][j]!=Q0),None)
        if k is None:continue
        A[row],A[k]=A[k],A[row];v=qi(A[row][j]);A[row]=[qm(v,a) for a in A[row]]
        for k in range(len(A)):
            if k!=row and A[k][j]!=Q0:
                v=A[k][j];A[k]=[qa(a,qn(qm(v,b))) for a,b in zip(A[k],A[row])]
        row+=1
        if row==len(A):break
    return row


def rank_three_witness(data):
    record=data['rank_three_witness'];assert record['name']=='13a_373'
    original=json.loads((HERE/'reflective_augmentation_inputs.json').read_text())
    assert record==next(c for c in original['cases'] if c['name']==record['name'])
    vectors=[[tuple(map(Fraction,a)) for a in v] for v in record['normal_vectors']]
    assert len(vectors)==13 and all(len(v)==4 for v in vectors)
    half=(Fraction(1,2),Fraction(0))
    C=[[qm(half,qsum(qm(a,b) for a,b in zip(u,v))) for v in vectors] for u in vectors]
    assert all(C[i][i]==(Fraction(2),Fraction(0)) for i in range(13))
    crossings,_=diagram(data['case']['pd'])
    for r in crossings:
        o,a,b=r['over'],r['input'],r['output']
        for l in range(13):assert qa(C[l][a],C[l][b])==qm(C[l][o],C[o][a])
    assert qrank(vectors)==qrank(C)==3
    return {'rank_three_cord_witness_field':'Q(phi), phi^2-phi-1=0',
            'normal_vector_rank':3,'cord_matrix_rank':3,'cord_equations_checked':13*13,
            'used_only_for_augmentation_lower_bound':True}


def matrix_product(R,A,B):
    return [[R.add(*(R.mul(A[i][k],B[k][j]) for k in range(len(B))))
             for j in range(len(B[0]))] for i in range(len(A))]
def matrix_power(R,A,n):
    out=basis(R,len(A))
    for _ in range(n):out=matrix_product(R,out,A)
    return out
def determinant(R,A):
    out=R.c(0)
    for p in permutations(range(len(A))):
        sign=(-1)**sum(p[i]>p[j] for i in range(len(A)) for j in range(i+1,len(A)))
        out=R.add(out,R.mul(R.c(sign),*(A[i][p[i]] for i in range(len(A)))))
    return out


def torus_control(data):
    R=Ring('mu');mu=R.v('mu');d=R.sub(R.c(1),mu);I=basis(R,4)
    A=[[d if i==j else R.neg(mu) if j==i+1 else R.c(1) if i==j+1 else {}
        for j in range(4)] for i in range(4)]
    matrices=[]
    for i in range(4):
        T=deepcopy(I);T[i]=[R.sub(T[i][j],A[i][j]) for j in range(4)];matrices.append(T)
    Y=deepcopy(I)
    for T in matrices:Y=matrix_product(R,Y,T)
    X=matrix_product(R,matrices[0],Y)
    scalar=[[R.neg(R.v('mu',5)) if i==j else {} for j in range(4)] for i in range(4)]
    assert matrix_power(R,X,4)==matrix_power(R,Y,5)==scalar
    assert all(matrix_product(R,Y,matrices[i])==matrix_product(R,matrices[i+1],Y) for i in range(3))
    det=R.add(R.c(1),R.neg(mu),R.v('mu',2),R.neg(R.v('mu',3)),R.v('mu',4))
    assert determinant(R,A)==det and det
    # These nonzero adjacent pairings prove irreducibility: any invariant
    # subspace contains some root, then all four roots. The image contains
    # every T_i because Y transports them and T0=XY^-1.
    assert all(A[i+1][i]==R.c(1) and A[i][i+1]==R.neg(mu) for i in range(3))
    output={'field':'Q(mu)','knot':'T(4,5)','meridian_word':'X Y^-1',
            'torus_relation':'X^4=Y^5=-mu^5 I','pairing_determinant':R.dump(det),
            'dimension':4,'irreducibility_checked_by_nonzero_adjacent_pairings':True}
    assert output==data['positive_control']
    # Invalid torus relation must be rejected rather than silently accepted.
    bad=deepcopy(X);bad[0][0]=R.add(bad[0][0],R.c(1))
    assert matrix_power(R,bad,4)!=matrix_power(R,Y,5)
    return output|{'corrupted_torus_matrix_rejected':True}


def torus_PD_control(data):
    record=data['positive_PD_control'];assert record['name']=='T(4,5)'
    R=Ring('mu');mu=R.v('mu');mi=R.v('mu',-1);d=R.sub(R.c(1),mu);I=basis(R,4)
    decode=lambda p:{(e,):a for e,a in p}
    roots=[[decode(p) for p in u] for u in record['arc_roots']]
    forms=[[decode(p) for p in f] for f in record['arc_covectors']]
    assert len(roots)==len(forms)==15
    matrices=[];inverses=[]
    for u,f in zip(roots,forms):
        assert any(u) and any(f) and dot(R,I,f,u)==d
        matrices.append([[R.sub(I[i][j],R.mul(u[i],f[j])) for j in range(4)] for i in range(4)])
        inverses.append([[R.add(I[i][j],R.mul(mi,u[i],f[j])) for j in range(4)] for i in range(4)])
        assert matrix_product(R,matrices[-1],inverses[-1])==I
    signed,_=diagram(record['pd']);assert len(signed)==15
    for crossing in signed:
        o,a,b,e=[crossing[k] for k in ['over','input','output','exponent']]
        first,last=(matrices[o],inverses[o]) if e==1 else (inverses[o],matrices[o])
        assert matrix_product(R,matrix_product(R,first,matrices[a]),last)==matrices[b]
    # The literal PD is exactly the closure of (sigma1 sigma2 sigma3)^5.
    # Check the dart-map binding, and check all incoming/outgoing wire labels
    # through the explicitly retained braid, including the final closure.
    census_bijection(record['pd'],record['braid_pd'],
                     {'dart_images':record['dart_map_to_braid'],'orientation':1,'over_under_swap':0})
    assert record['braid_word']==[1,2,3]*5 and len(record['braid_pd'])==15
    top=[None]*4;current=[None]*4
    for generator,row in zip(record['braid_word'],record['braid_pd']):
        i=generator-1
        for wire,edge in [(i,row[0]),(i+1,row[3])]:
            if current[wire] is None:top[wire]=current[wire]=edge
            assert current[wire]==edge
        current[i],current[i+1]=row[1],row[2]
    assert current==top and None not in top
    A=[[d if i==j else R.neg(mu) if j==i+1 else R.c(1) if i==j+1 else {}
        for j in range(4)] for i in range(4)]
    chain=[]
    for i in range(4):
        T=deepcopy(I);T[i]=[R.sub(T[i][j],A[i][j]) for j in range(4)];chain.append(T)
    assert all(T in matrices for T in chain)
    seeds=record['seed_arcs'];assert seeds==[0,1,5,9]
    root_det=determinant(R,[roots[i] for i in seeds])
    form_det=determinant(R,[forms[i] for i in seeds])
    assert root_det==decode(record['seed_root_determinant'])==R.c(1)
    assert form_det==decode(record['seed_covector_determinant'])==determinant(R,A)
    assert R.mul(root_det,form_det)==decode(record['seed_pairing_determinant']) and form_det
    # A deliberately reversed individual sign must fail the same literal
    # crossing equations; a genuine noninvolutory generic control detects it.
    rejected=False
    for crossing in signed:
        o,a,b=[crossing[k] for k in ['over','input','output']]
        if matrix_product(R,matrix_product(R,matrices[o],matrices[a]),inverses[o])!=matrices[b]:
            rejected=True;break
    assert rejected
    return {'generic_torus_literal_crossings_checked':15,'torus_braid_word_binding_checked':True,
            'torus_seed_root_determinant':1,'torus_seed_covector_determinant':'1-mu+mu^2-mu^3+mu^4',
            'torus_image_contains_irreducible_four_root_chain':True,
            'reversed_generic_control_crossing_sign_rejected':True}


def trefoil_sum_control(data):
    assert data['connected_sum_scope']['common_parameter']==-1
    # Rational2D S3 representation; exact noncommuting involutions, with
    # displacement roots and covectors both spanning, is irreducible.
    R=Ring('mu');I=basis(R,2);c=R.c
    A=[[c(-1),c(1)],[{},c(1)]];B=[[c(1),{}],[c(1),c(-1)]]
    assert matrix_product(R,A,A)==matrix_product(R,B,B)==I
    assert matrix_product(R,matrix_product(R,A,B),A)==matrix_product(R,matrix_product(R,B,A),B)
    assert matrix_product(R,A,B)!=matrix_product(R,B,A)
    # Roots e0,e1, covectors(2,-1),(-1,2) have determinant3 !=0 overQ;
    # both cross-pairings are nonzero, so the root propagation proof applies.
    assert determinant(R,[[c(2),c(-1)],[c(-1),c(2)]])==c(3)
    return {'trefoil_rank2_minus_one_control':True,
            'sum_rank_statement':'ar_C(13a373 #^s trefoil)=s+3; independently mr=b=s+4',
            'sum_augmentation_import':'COR17 Theorem1.9; same-parameter lower construction and unconditional upper additivity',
            'sum_meridional_import':'INVOLUTION_COVER_SUM_FORMULAS.tex cor:hundredone'}


def named_census(data):
    import importlib.metadata as metadata
    import snappy_15_knots
    from spherogram import Link
    source=data['census_provenance'];case=data['case']
    assert metadata.version('snappy_15_knots')==source['version']
    assert metadata.version('spherogram')==source['spherogram_version']
    root=Path(snappy_15_knots.__file__).parent
    for pin in source['source_files']:
        raw=(root/pin['path']).read_bytes()
        assert len(raw)==pin['bytes'] and hashlib.sha256(raw).hexdigest()==pin['sha256']
    link=Link('K'+case['name'].replace('_',''))
    assert [list(row) for row in link.PD_code()]==case['census_pd']
    assert json.loads(json.dumps(link.DT_code()))==case['census_DT']
    control=data['positive_PD_control']
    assert [list(row) for row in Link('T(4,5)').PD_code()]==control['pd']
    assert [list(row) for row in Link(braid_closure=control['braid_word']).PD_code()]==control['braid_pd']
    return 1


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--census',action='store_true');args=parser.parse_args()
    data=json.loads(INPUT.read_text());assert data['names']==['13a_373']
    for path,sha in data['source_sha256'].items():
        assert hashlib.sha256((HERE/path).read_bytes()).hexdigest()==sha,path
    results={}
    for fn in (case_binding,raw_first_loop,paired_core,unit_certificate,rank_three_witness,torus_control,torus_PD_control,trefoil_sum_control):
        results.update(fn(data))
    # Source identity corruption and signed crossing corruption are both
    # rejected by the exact diagram binding, independently of the algebra.
    bad=deepcopy(data);bad['signed_trace'][0][3]*=-1
    try:case_binding(bad)
    except AssertionError:pass
    else:raise AssertionError('corrupted crossing sign accepted')
    bad=deepcopy(data);bad['case']['census_map_witness']['dart_images'][0]=bad['case']['census_map_witness']['dart_images'][1]
    try:case_binding(bad)
    except AssertionError:pass
    else:raise AssertionError('corrupted named map accepted')
    count=named_census(data) if args.census else None
    files=[Path(__file__),INPUT,HERE/'GENERIC_KCH_OBSTRUCTION.tex']
    receipt={'status':'PASS','names':data['names'],'checks':results,
             'external_named_census_reread_count':count,
             'corruptions_rejected':['unit_multiplier','crossing_sign','named_planar_map','torus_control_matrix','generic_positive_control_crossing_sign'],
             'independent_meridional_rank_dependency':data['independent_meridional_rank_dependency'],
             'input_sha256':{p.relative_to(HERE).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}|data['source_sha256'],
             'representation_scope':data['representation_scope'],'augmentation_scope':data['augmentation_scope'],
             'assurance':'Exact internal algebra/diagram replay. The written irreducibility and nonzero-cancellation arguments, primary augmentation correspondence, and independent meridional lower certificate remain distinct proof dependencies; no worldwide-priority or formal-verification claim.'}
    receipt_path=HERE/'generic_kch_obstruction_verification.json'
    if receipt_path.is_file():
        previous=json.loads(receipt_path.read_text())
        if 'internal_pdf_reviews' in previous:
            # Reviews retain their own source/PDF hashes. A later source edit
            # does not turn an earlier edition's review into current assurance.
            receipt['internal_pdf_reviews']=previous['internal_pdf_reviews']
    receipt_path.write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','names':data['names'],'integer_Laurent_identity':1,
                      'generic_torus_rank4_control':True,'rank_three_cord_attainment':True,
                      'optional_census_rereads':count},indent=2))


if __name__=='__main__':main()
