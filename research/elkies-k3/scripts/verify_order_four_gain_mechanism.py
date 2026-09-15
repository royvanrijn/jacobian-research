#!/usr/bin/env python3
"""Exact order-four control and the fixed-field gate for the retained Inose involution."""
import argparse
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import resource
import time
import sympy as S

ROOT=Path(__file__).resolve().parents[2]
SOURCE='artifacts/generated-results/elkies-k3-det622-inose-source-v1/certificate.json'
t,w,s,x,y,k,a,b=S.symbols('t w s x y k a b')
def zero(f):assert S.cancel(f)==0

def verify():
 start=time.monotonic()
 A=-t**6+5*t**4-t**2;B=t**12+t**10+3*t**7-3*t**5+t**2+1;d=t*t+1
 x1=t;g1=t**5+1;x2=-t**3;g2=1-t**5
 zero(d*g1*g1-x1**3-A*x1-B);zero(d*g2*g2-x2**3-A*x2-B)
 zero(t**8*A.subs(t,-1/t)-A);zero(t**12*B.subs(t,-1/t)-B)
 zero(t**4*x1.subs(t,-1/t)-x2);zero(-t**5*g1.subs(t,-1/t)-g2)
 # beta(t,w)=(-1/t,w/t) preserves the conic and squares to its deck involution.
 zero((w/t)**2-((-1/t)**2+1)-(w*w-d)/t**2)
 zero(-1/(-1/t)-t);zero((w/t)/(-1/t)+w)
 tt=(s*s-1)/(2*s);ww=(s*s+1)/(2*s)
 zero(ww*ww-tt*tt-1);zero((ww-1)/tt-(s-1)/(s+1))
 D=S.Poly(4*A**3+27*B**2,t);assert D.degree()==24
 assert S.gcd(D,S.Poly(d,t)).degree()==0
 zero(x1-x2-t*d);zero(g1-g2-2*t**5)
 p=13;dp=S.Poly(D,t,modulus=p);assert dp.degree()==24 and S.gcd(dp,dp.diff()).degree()==0
 counts=[]
 for z in range(p+1):
  aa=int(A.subs(t,z))%p if z<p else 0;bb=int(B.subs(t,z))%p if z<p else 1
  legendre=p+1+sum(0 if (v:=(xx**3+aa*xx+bb)%p)==0 else 1 if pow(v,(p-1)//2,p)==1 else -1 for xx in range(p))
  explicit=1+sum((yy*yy-xx**3-aa*xx-bb)%p==0 for xx in range(p) for yy in range(p))
  assert legendre==explicit;counts.append(explicit)
 surface=sum(counts);trace=surface-1-p*p;ns_upper=(trace+22*p)//(2*p);mw_upper=ns_upper-2
 assert (surface,trace,ns_upper,mw_upper)==(200,30,12,10)
 # General Inose reversal: k is the last coefficient of t^5(t²+b*t+k).
 AA=a*t**4;BB=t**5*(t*t+b*t+k)
 zero(AA.subs(t,k/t)-k**4/t**8*AA);zero(BB.subs(t,k/t)-k**6/t**12*BB)
 # Its negative-y lift eta squares to identity and has symplectic differential factor1.
 zero((-k/t**2)*(k**2/t**4)/(-k**3/t**6)-1)
 data=json.loads((ROOT/SOURCE).read_text());I=S.Rational(data['I']);J=S.Rational(data['J']);kk=I**2*J**3
 assert I<0 and J<0 and kk<0
 q=Fraction(str(-J));lo=isqrt(q.numerator*q.denominator)
 # Here -J is integral, making the strict nonsquare witness particularly short.
 assert q.denominator==1 and lo*lo<q.numerator<(lo+1)**2
 residual=S.Poly(4*(-3*I*J)**3*t*t+27*(t*t-2*I*J**2*t+kk)**2,t)
 assert S.gcd(residual,S.Poly(t*t-kk,t)).degree()==0
 return {'schema':'order-four-gain-mechanism-v1','script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'control':{'A':str(A),'B':str(B),'cover':'w^2=t^2+1','P':['t','w*(t^5+1)'],'Q':['-t^3','w*(1-t^5)'],
                    'point_and_action_identities':True,'genus':0,'height_gram':[[8,0],[0,8]],
                    'height_boundary':'Written pole/intersection proof; polynomial support checks retained.',
                    'good_prime':p,'fibre_counts':counts,'surface_count':surface,'H2_trace':trace,'rational_NS_upper':ns_upper,'arithmetic_MW_upper':mw_upper},
         'inose':{'source':SOURCE,'source_sha256':sha256((ROOT/SOURCE).read_bytes()).hexdigest(),'I':str(I),'J':str(J),'k':str(kk),
                  'negative_J_floor_sqrt':str(lo),'negative_J_is_square':False,'fixed_fibres_smooth':True,
                  'fixed_point_field':'Q(sqrt(J)), imaginary quadratic and not Q(i)',
                  'scope':'Canonical Inose reversal and its Q-conjugates, on every equivariant rational fibration and every constant twist; other involutions are not excluded.'},
         'limits':{'cpu_seconds':20,'memory_bytes':1024**3},'elapsed_seconds':time.monotonic()-start,'software':{'sympy':S.__version__},
         'boundary':'Conditional even-rank/orthogonal-orbit theorem and explicit low-rank K3 control. No arithmetic MW17 positive endpoint. Geometry, height and Weil bounds remain written proof dependencies.'}

if __name__=='__main__':
 resource.setrlimit(resource.RLIMIT_CPU,(20,20));resource.setrlimit(resource.RLIMIT_AS,(1024**3,1024**3))
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--record',type=Path,required=True);args=p.parse_args();result=verify()
 args.record.parent.mkdir(parents=True,exist_ok=True)
 with args.record.open('x') as f:json.dump(result,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps({'control_parent_rank_upper':10,'control_gram':[[8,0],[0,8]],'inose_canonical_involution_passes':False,'elapsed_seconds':result['elapsed_seconds']}))
