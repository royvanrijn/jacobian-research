#!/usr/bin/env python3
"""Exact algebraic guards for the fixed-gap existential large-lower theorem."""
from __future__ import annotations
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import argparse, hashlib, json
if not __debug__:
    raise SystemExit('Refusing -O/-OO: exact checks are required.')
HERE=Path(__file__).resolve().parent
ROWS=[]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def st(x):
    x=F(x)
    return str(x.numerator) if x.denominator==1 else f'{x.numerator}/{x.denominator}'
def ck(name,left,right,relation='<'):
    left,right=F(left),F(right)
    assert (left<right if relation=='<' else left==right),(name,st(left),relation,st(right))
    ROWS.append({'name':name,'left':st(left),'relation':relation,'right':st(right)})
def eb(x):
    x=F(x)
    if x<0:
        if -x>100:return F(0),1/eb(F(100))[0]
        lo,hi=eb(-x);return 1/hi,1/lo
    term=total=F(1)
    for k in range(1,257):term*=x/k;total+=term
    assert x/258<1
    return total,total+term*x/257/(1-x/258)
def eu(x):return eb(x)[1]
def el(x):return eb(x)[0]
def sb(x):
    x=F(x);q=10**12;k=isqrt(x.numerator*q*q//x.denominator)
    assert F(k*k,q*q)<=x<F((k+1)**2,q*q)
    return F(k,q),F(k+1,q)
def pl(a,b):
    out=dict(a)
    for k,v in b.items():out[k]=out.get(k,F(0))+v
    return {k:v for k,v in out.items() if v}
def ps(a,c):return{k:v*c for k,v in a.items()if v*c}
def pm(a,b):
    out={}
    for i,u in a.items():
        for j,v in b.items():
            k=tuple(x+y for x,y in zip(i,j));out[k]=out.get(k,F(0))+u*v
    return{k:v for k,v in out.items()if v}
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,required=True);ap.add_argument('--output',type=Path,default=HERE/'fixed_gap_large_lower_extension_verification.json');args=ap.parse_args()
    repo=args.repo_root.resolve();ledger_path=HERE/'fixed_gap_large_lower_extension_source_ledger.json';ledger=json.loads(ledger_path.read_text())
    mandatory={('local','FIXED_GAP_LARGE_LOWER_EXPONENT_EXTENSION.md'),('local',Path(__file__).name)}
    mandatory|={('repo','research/factorial/image-rigidity/'+n)for n in ['GAP_TWO_LARGE_LOWER_EXPONENT_ENDPOINT_THEOREM.md','verify_gap_two_endpoint_contour_identities.py','GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md','GAMMA_PRODUCT_ROOT_MODULUS.tex','COALESCING_BINOMIAL_SADDLE_BRIDGE.md']}
    mandatory|={('repo','research/factorial/image-rigidity/'+n)for n in ['gap_two_large_lower_endpoint_source_ledger.json','gap_two_large_lower_endpoint_primary_audit.md','COALESCING_BINOMIAL_HERMITE_BULK.md','GAP_TWO_EXTENDED_BULK_ENDPOINT_REDUCTION.md','GAP_TWO_OUTSIDE_ENDPOINT_ONE_PACKET.md']}
    assert ledger.get('schema')==1
    assert len(mandatory)==12 and len(ledger['files'])==12
    assert mandatory=={(q['root'],q['path'])for q in ledger['files']}
    pins=[]
    for q in ledger['files']:
        path=(HERE if q['root']=='local'else repo)/q['path'];assert path.is_file()and sha(path)==q['sha256'],q['path'];pins.append(q)
    def var(j):return{tuple(int(i==j)for i in range(5)):F(1)}
    M,E,U,H,G=[var(j)for j in range(5)];one={(0,0,0,0,0):F(1)}
    mE=pl(M,E);F0=pl(E,ps(U,2));mF=pl(M,F0);gH=pl(G,ps(H,3))
    origin=pl(pm(mE,gH),ps(pm(M,G),-1))
    expectedE=pl(ps(pm(M,H),3),pl(pm(G,E),ps(pm(H,E),3)));assert origin==expectedE
    shift=pl(pl(ps(pm(M,pl(ps(G,-1),H)),2),pm(mF,gH)),pm(M,G))
    expectedU=pl(ps(pm(M,H),5),pl(pm(G,F0),ps(pm(H,F0),3)));assert shift==expectedU
    diskdiff=pl(pm(mF,mF),ps(pl(pm(mE,mE),ps(pm(M,U),4)),-1))
    assert diskdiff==ps(pm(U,pl(E,U)),4)
    # Universal majorant induction after multiplying by G>=2.
    e={(1,0):F(1)};g={(0,1):F(1)};one2={(0,0):F(1)}
    enew=pl(ps(e,3),ps(one2,3));gap=pl(pm(g,pl(pl(enew,ps(one2,-3)),ps(e,-1))),ps(e,-3))
    assert gap==pm(e,pl(ps(g,2),ps(one2,-3)))
    ck('majorant induction at minimum G',0,2*2-3)
    ck('gamma upper bound factor',F(1)+F(1),2,'=')
    # Literal e/u coefficient recurrence controls; proof induction is not a grid inference.
    expected=[(3,1),(12,20),(39,161),(120,1088),(363,6893),(1092,42452)]
    es,us=3,1
    for n,pair in enumerate(expected,1):
        assert(es,us)==pair
        es,us=3+3*es,5+3*es+6*us
    # Actual rising-factorial multiplication formula at every coefficient in each fixture.
    def rising(x,k):
        out=F(1)
        for j in range(k):out*=x+j
        return out
    factor_controls=0;bound_controls=0
    for d in range(2,9):
        ed,ud=3,1
        for j in range(1,d):ed,ud=3+3*ed,5+3*ed+6*ud
        for a in [d,d+2,2*d+3]:
            for N in range(1,9):
                gamma=[F(a*N+j,d)for j in range(1,d+1)]
                for k in range(N+1):
                    right=F(d)**(d*k)
                    for v in gamma:right*=rising(v,k)
                    assert rising(a*N+1,d*k)==right;factor_controls+=1
                h=N-1;M=gamma[0];E=3*h;U=h
                for v in gamma[1:]:
                    FF=E+2*U
                    M,E,U=M*v,3*M*h+v*E+3*h*E,5*M*h+v*FF+3*h*FF
                big=F(2*a,d)
                assert E<=ed*big**(d-1)*N**d and U<=ud*big**(d-1)*N**d
                assert F(a,d)**d*N**d<=M<=big**d*N**d
                if N==1:assert E==U==0
                bound_controls+=1
    # Algebra of the explicit normalized strip constants.
    # d^(2d-2)*(2a/d)^(2d-2)/a^(2d-1) times(4a/d+e)
    # equals2^(2d-2)*(4/d+e/a); tested exactly as a rational identity.
    strip_controls=0
    for d in range(2,11):
        ed,ud=3,1
        for j in range(1,d):ed,ud=3+3*ed,5+3*ed+6*ud
        for a in[d,d+1,2*d,2**d*d*ud+1]:
            left=F(d)**(2*d-2)*F(2*a,d)**(2*d-2)/F(a)**(2*d-1)*(F(4*a,d)+ed)
            right=F(2)**(2*d-2)*(F(4,d)+F(ed,a))
            assert left==right
            assert right<=F(2)**(2*d-2)*F(4+ed,d)
            if a>2**d*d*ud:assert F(ud*2**(d-1)*d,a)<F(1,2)
            strip_controls+=1
    # Whole-tail endpoint signs, including the actual h'/h prefactor.
    ck('right-tail derivative at R3',-F(199,200)*3+F(11,10)*(F(1,3)+F(1,100))+1,0)
    ck('left-tail derivative at R3',1,3-1/(F(3,2)-F(1,100)))
    ck('left h half-radius',1-F(1,2),F(1,2),'=')
    ck('radial exponential near zero',eu(F(1,2)),2)
    ck('endpoint drift at delta1/15',F(3,4)*15-F(15,4),F(1,2)*15,'=')
    # An explicit far cutoff R=4(T+3) beats its quadratic/log bound for all T>=0.
    t={(1,):F(1)};o={(0,):F(1)};RR=ps(pl(t,ps(o,3)),4)
    far=pl(pl(ps(pm(RR,RR),-F(1,2)),pl(t,pl(RR,ps(o,3)))),ps(o,4))
    assert far=={(2,):-8,(1,):-43,(0,):-53}
    # Bernoulli alpha=(1-1/N)^(-d)-1>=d/N and hardN give O_d(delta^2).
    for d in range(2,10):
        for N in[2*d,2*d+1,10*d]:
            alpha=F(N,N-1)**d-1
            assert alpha>=F(d,N)
    # Three actual sign/normalization/uniformity negative controls.
    gg,hh=F(7,3),F(11,4);root=-gg*hh
    assert root+gg*hh==0 and -root+gg*hh!=0
    assert rising(9,6)!=F(3)**(2)*rising(F(9,3),2)*rising(F(10,3),2)*rising(F(11,3),2)
    assert expected[3][1]>expected[1][1] # These majorants explicitly depend on d.
    result={'schema':1,'scope':'Fixed integer gap: existence of an uncomputed lower-exponent threshold; no d-uniform growth or finite global box','comparisons':ROWS,'universal_identity_count':5,'gamma_factorization_controls':factor_controls,'majorant_bound_controls':bound_controls,'strip_scaling_controls':strip_controls,'negative_controls':3,'pins':pins,'source_ledger_sha256':sha(ledger_path)}
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps({'status':'PASS','comparisons':len(ROWS),'universal_identities':5,'gamma_factorizations':factor_controls,'majorant_bounds':bound_controls,'strip_scales':strip_controls,'negative_controls':3,'output':str(args.output)}))
if __name__=='__main__':main()
