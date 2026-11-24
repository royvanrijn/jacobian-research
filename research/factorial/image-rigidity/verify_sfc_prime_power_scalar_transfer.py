#!/usr/bin/env python3
"""Exact finite controls for the uniform prime-power proof, not its quantifiers."""
if not __debug__:
    raise RuntimeError("This verifier requires assertions; do not run Python with -O.")

from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from math import comb, factorial, prod
from pathlib import Path
import json
import sys


EXPECTED_SOURCE_SHA256 = {
    "SFC_PRIME_POWER_SCALAR_TRANSFER.md": "9fa457f230397c670c7b4152c0510df93da24df03fc08112f1c5bf5c73bc0349",
    "ARITHMETIC_FRONTIER_TERMINAL_ORDERS.tex": "9721978bcb795af9169c695c5fb8348974069c5148e5e0caa8bd3c4d65650e07",
    "LOWER_TWO_FOURTH_POWER_GAP_SFC.tex": "e238569c21a0ccc7fb507d910608a734461fefaab4d96d5df0124afe73407631",
    "DIAGONAL_PRIME_GAP_SFC.tex": "db034b5fb29ece9c3a0cb8cbce4c82f692b933e3679b27ee823db5ca5e7f98eb"
}


def root_dir():
    for parent in Path(__file__).resolve().parents:
        if (parent / "research/factorial/image-rigidity").is_dir():
            return parent
    raise RuntimeError("Cannot locate the research worktree.")


def verify_sources():
    proofdir = root_dir()/"research/factorial/image-rigidity"
    found = {}
    if len(EXPECTED_SOURCE_SHA256) != 4:
        raise RuntimeError("Expected the own-note pin and exactly three retained-source pins.")
    for name, expected in EXPECTED_SOURCE_SHA256.items():
        actual = sha256((proofdir/name).read_bytes()).hexdigest()
        if actual != expected:
            raise RuntimeError("Source SHA256 mismatch: "+name)
        found[name] = actual
    return found


def trim(a):
    a = list(a)
    while len(a)>1 and not a[-1]: a.pop()
    return a


def add(a,b):
    c=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a): c[i]+=x
    for i,x in enumerate(b): c[i]+=x
    return trim(c)


def mul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]+=x*y
    return trim(c)


def ev(a,x):
    c=F(0)
    for y in reversed(a): c=c*x+y
    return c


def vp(x,p):
    if not x: return 10**9
    x=F(x); a,b=abs(x.numerator),x.denominator; v=0
    while a%p==0: a//=p; v+=1
    while b%p==0: b//=p; v-=1
    return v


def residue(x,modulus):
    x=F(x)
    return x.numerator*pow(x.denominator,-1,modulus)%modulus


def rising(x,n):
    return prod((F(x)+j for j in range(n)), start=F(1))


def series_mul(a,b,D):
    return (mul(a,b)+[F(0)]*(D+1))[:D+1]


def series_div(a,b,D):
    out=[F(0)]*(D+1)
    for k in range(D+1):
        out[k]=(F(a[k] if k<len(a) else 0)-sum(b[j]*out[k-j]
                    for j in range(1,min(k,len(b)-1)+1)))/b[0]
    return out


@lru_cache(None)
def ah_inverse(p,D):
    e=[F(1)]
    for m in range(1,D+1):
        powers=[]; k=1
        while k<=m: powers.append(k); k*=p
        e.append(-sum(e[m-k] for k in powers)/m)
    assert all(vp(x,p)>=0 for x in e)
    return tuple(e)


def leaves(p,s,c,D,seed=None):
    q=p**s; initial_degree=q*(D+1)-1
    h=[F(1)]+[F(0)]*initial_degree if seed is None else (
        list(map(F,seed))+[F(0)]*(initial_degree+1))[:initial_degree+1]
    states=[dict(c=F(c),F={1:F(1)},h=h,j=0,e=0,scale=F(1))]
    for level in range(s):
        nxt=[]
        for state in states:
            degree=len(state['h'])-1; child=(degree-(p-1))//p
            b=F(state['F'].get(p,0)-1)
            assert vp(b,p)==0
            g=series_mul(ah_inverse(p,degree),state['h'],degree)
            newF={1:F(1)}; power=p
            while power<=child:
                newF[power]=(state['F'].get(p*power,0)-1)*F(p**(power-1),1)/b**power
                assert vp(newF[power],p)==power-1
                power*=p
            for j in range(p):
                e=residue(state['c']-j,p)
                hnew=[g[p*k+j]*(p/b)**k for k in range(child+1)]
                if seed is None:
                    assert vp(hnew[0],p)==0
                    assert all(vp(hnew[k],p)>=k for k in range(1,child+1))
                nxt.append(dict(c=(state['c']-j-e)/p,F=newF,h=hnew,
                    j=state['j']+p**level*j,e=state['e']+p**level*e,
                    scale=b*state['scale']**p/p))
        states=nxt
    b_s=(q-1)//(p-1)
    scales={a['scale'] for a in states}
    assert len(scales)==1
    nu=next(iter(scales))*p**b_s
    assert vp(nu,p)==0
    for state in states:
        assert q*state['c']+state['e']==c-state['j']
        assert state['e']==residue(c-state['j'],q)
    return sorted(states,key=lambda z:z['j']),nu


def apply_operator(h,a,Fpart,step,D):
    out=[(a-step*k)*h[k] for k in range(D+1)]
    for power,coef in Fpart.items():
        for k in range(power,D+1): out[k]+=step*coef*h[k-power]
    return out


def poly_operator(poly,h,a,Fpart,step,D):
    out=[F(0)]*(D+1); v=h[:]
    for x in poly:
        out=[a0+x*b0 for a0,b0 in zip(out,v)]
        v=apply_operator(v,a,Fpart,step,D)
    return out


def ordinary_operator(poly,c):
    D=len(poly)-1
    return trim(poly_operator(poly,[F(1)]+[F(0)]*D,F(c),{1:F(1)},1,D))


def upper_functional(P,n,gamma,q,h):
    c=gamma+q-1; C=gamma+q*n+q
    remainder=list(map(F,P))+[F(0)]*(q-len(P)); ds=[F(0)]*(q-1)
    for j in range(q-2,-1,-1):
        ell=[F(q*n),F(1)]
        for r in range(j): ell=mul(ell,[c-r,-1])
        d=remainder[j+1]*(-1)**j; ds[j]=d
        for k,x in enumerate(ell): remainder[k]-=d*x
    assert all(not x for x in remainder[1:])
    return sum(ds[j]*rising(C+q-j,h-q+j)
               for j in range(max(0,q-h),q-1))


def interpolation_control(p,s,n,gamma,D=2):
    gamma=F(gamma); assert gamma>0 and vp(gamma,p)>=0
    q=p**s; c=gamma+q-1; b_s=(q-1)//(p-1); B=b_s-s
    all_leaves,nu=leaves(p,s,c,D)
    ls=[state for state in all_leaves if state['e']]
    nodes=[c-state['j'] for state in ls]
    Q0=[F(1)]
    for a in nodes: Q0=mul(Q0,[-a,1])
    Gs=[]
    for i in range(q-1):
        g=[F(1)]
        for j,a in enumerate(nodes):
            if i!=j:g=mul(g,[-a,1])
        Gs.append(g)
    U=[[[F(0)]*(q-1) for _ in range(q-1)] for _ in range(D+1)]
    forcing=[[F(0)]*(q-1) for _ in range(D+1)]
    for i,state in enumerate(ls):
        a=nodes[i]; H=[x/state['h'][0] for x in state['h']]
        for ell in range(n):
            H=[x/(a+q*ell) for x in apply_operator(H,a+q*ell,state['F'],q,D)]
        assert H[0]==1
        assert all(vp(x,p)>=1 for x in H[1:])
        w=ev(Gs[i],a)
        for j,g in enumerate(Gs):
            row=series_div(poly_operator(g,H,a,state['F'],q,D),H,D)
            for k in range(D+1): U[k][i][j]=row[k]*nu**k/w
            assert U[0][i][j]==int(i==j)
            assert all(vp(U[k][i][j],p)>=1 for k in range(1,D+1))
        row=series_div(poly_operator(Q0,H,a,state['F'],q,D),H,D)
        for k in range(D+1): forcing[k][i]=row[k]*nu**k/w
        assert forcing[0][i]==0
        assert all(vp(forcing[k][i],p)>=s for k in range(1,D+1))
    bs=[[F(0)]*(q-1) for _ in range(D+1)]
    for k in range(1,D+1):
        for i in range(q-1):
            bs[k][i]=-forcing[k][i]-sum(U[r][i][j]*bs[k-r][j]
                            for r in range(1,k+1) for j in range(q-1))
            assert vp(bs[k][i],p)>=s
            assert residue(bs[k][i]/q,p)==(residue(-nu,p) if k==1 else 0)
    assert vp(ev(Q0,-q*n),p)==B
    assert all(vp(ev(g,-q*n),p)>=B-s+1 for g in Gs)
    for k in range(1,D+1):
        correction=sum(bs[k][i]*ev(Gs[i],-q*n) for i in range(q-1))
        assert vp(correction,p)>=B+1
    assert upper_functional(Q0,n,gamma,q,1)==0
    assert upper_functional(Q0,n,gamma,q,2)==(-1)**q
    assert all(upper_functional(g,n,gamma,q,2)==0 for g in Gs)
    if q>=3:
        jstar=residue(c,q)
        assert upper_functional(Q0,n,gamma,q,3)==(-1)**q*(2*(gamma+q*n)+4+jstar)
        assert all(upper_functional(g,n,gamma,q,3)==-(-1)**q for g in Gs)
    # Split the actual original polynomial output independently of the
    # composed leaf operator evaluation, at all q leaves.
    Pi=[F(1)]
    for ell in range(n): Pi=mul(Pi,[q*ell,1])
    probe=mul(Pi,[2,-3,1]); original=ordinary_operator(probe,c)
    assert len(original)-1<=q*(D+1)-1
    output,_=leaves(p,s,c,D,seed=original)
    for lhs,seed in zip(output,all_leaves):
        rhs=poly_operator(probe,seed['h'],c-seed['j'],seed['F'],q,D)
        assert lhs['h']==rhs
    return dict(p=p,s=s,q=q,n=n,gamma=str(gamma),nu=str(nu),
                checked_outer_degree=D,actual_split_leaves=q,
                restricted_seed_coefficients_checked=D+1,
                first_variation_coefficients_checked=D)


def sharp_controls():
    controls=[]
    for p,s in [(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(5,1)]:
        q=p**s;b=(q-1)//(p-1);B=b-s
        chi=[F(1)]
        for j in range(1,q+1): chi=mul(chi,[-j,1])
        u=(-1)**q*factorial(q);chi[0]-=u
        assert chi[0]==0
        Q=chi[1:]
        normalized=[x*p**((s-1)*j)/p**B for j,x in enumerate(Q)]
        assert all(vp(x,p)>=0 for x in normalized)
        assert vp(normalized[0],p)==0 and vp(Q[-1],p)==0
        assert all(vp(ev(Q,q*a),p)==B for a in range(-3,4))
        actual=ordinary_operator([0]+Q,q)
        assert actual==[-u]+[F(0)]*(q-1)+[F(1)]
        controls.append(dict(p=p,s=s,q=q,B=B,u=u))
    return controls


def difficult_actual_case():
    p,q,s,n,gamma,h=5,5,1,1,F(11),5
    c=gamma+q-1;jstar=residue(c,q)
    nodes=[c-j for j in range(q) if j!=jstar]
    Q0=[F(1)]
    for a in nodes:Q0=mul(Q0,[-a,1])
    Ts=[]
    for j in range(q-1):
        G=[F(1)]
        for k,a in enumerate(nodes):
            if k!=j:G=mul(G,[-a,1])
        Ts.append(upper_functional(G,n,gamma,q,h))
    T0=upper_functional(Q0,n,gamma,q,h)
    C=gamma+q*n+q;kappa=F(q*(n+1),1)/rising(C,h)
    assert T0==-32400 and Ts==[1146,1184,1224,1266]
    assert sum(Ts)==4820
    assert vp(kappa,5)==-1 and vp(T0,5)==2
    assert min(vp(x,5) for x in Ts)==0 and vp(sum(Ts),5)==1
    root=-1/rising(gamma,q)
    upper=sum(F(comb(2,k))*rising(gamma+h,q*k)*root**k for k in range(3))
    assert upper>0
    return dict(p=5,s=1,q=5,n=1,gamma='11',h=5,T0=int(T0),
                Tj=[int(x) for x in Ts],sum_Tj=int(sum(Ts)),
                kappa=str(kappa),kappa_valuation=-1,
                exact_lower_root=str(root),upper_at_lower_root=str(upper),
                meaning='actual nondegeneracy failure; exact pair still coprime')


def arithmetic_controls():
    count=0; period_count=0
    for p,s in [(2,1),(2,2),(2,3),(2,4),(3,1),(3,2),(5,1),(5,2)]:
        q=p**s
        for a in range(1,q):
            w=vp(a,p);A=a//p**w;K=2*s-vp(factorial(a),p)
            for N in range(2,22):
                kappa=F(q,(q+a)*prod((q+a)*N-j for j in range(1,a)))
                rhs=s-vp(factorial(a),p)-vp(comb((A+p**(s-w))*N-1,A-1),p)
                assert vp(kappa,p)==rhs
                if A==1:assert rhs==s-(a-1)//(p-1)
                if 1<A<p and K>=1:
                    expected=any(((A+p**(s-w))*N-j)%p**K==0 for j in range(1,A))
                    assert (vp(kappa,p)<=-s)==expected
                count+=1
            if A>1 and K>0:
                m=A-1;L=0
                while p**(L+1)<=m:L+=1
                modulus=p**(K+L);factor=A+p**(s-w)
                for N in range(2,7):
                    assert (comb(factor*(N+modulus)-1,m)-comb(factor*N-1,m))%p**K==0
                    if K>L:
                        balls=[]
                        for j in range(1,m+1):
                            delta=vp(m*comb(m-1,j-1),p)
                            balls.append((factor*N-j)%p**(K+delta)==0)
                        assert sum(balls)<=1
                        assert (vp(comb(factor*N-1,m),p)>=K)==any(balls)
                    period_count+=1
    return dict(ordinary_valuation_identity_controls=count,
                finite_precision_and_disjoint_ball_controls=period_count)


def gram_controls():
    count=0
    for beta in [F(1,3),F(1,2),F(1),F(5,3),F(7,2)]:
        for i in range(7):
            for j in range(7):
                rhs=sum(comb(i,k)*rising(beta+k,i-k)*factorial(k)*rising(beta,k)
                        *comb(j,k)*rising(beta+k,j-k) for k in range(min(i,j)+1))
                assert rhs==rising(beta,i+j)
                count+=1
        assert all(factorial(k)*rising(beta,k)>0 for k in range(7))
    return count


def main():
    pins=verify_sources()
    params=[(2,1,1,F(1,3)),(3,1,2,F(1,2)),(2,2,3,F(5,3)),
            (5,1,1,F(11)),(2,3,1,F(1,3)),(3,2,2,F(1,2)),
            (2,4,1,F(3,5)),(3,3,1,F(2,5))]
    records=[interpolation_control(*row) for row in params]
    # A formal unit constant alone does not imply a closed-disk unit.
    assert ev([1,-1],1)==0
    report=dict(status='PASS',python=sys.version.split()[0],source_sha256=pins,
        scope='finite exact controls; all-parameter conclusions rely on the written proofs',
        actual_leaf_controls=records,sharp_stage_controls=sharp_controls(),
        actual_deep_response_countercontrol=difficult_actual_case(),
        arithmetic_controls=arithmetic_controls(),gamma_gram_entry_controls=gram_controls(),
        formal_unit_closed_disk_countercontrol={'series':'1-v','constant':1,'value_at_1':0},
        all_parameter_quantifiers_certified_by_finite_tests=False)
    path=root_dir()/"research/factorial/image-rigidity/sfc_prime_power_scalar_transfer_verification.json"
    path.write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
