#!/usr/bin/env python3
"""Complete independent actual-moment terminal-six Bezout replay."""
if not __debug__:raise SystemExit('Exact assertions require unoptimized Python')
from pathlib import Path
from fractions import Fraction as F
from math import comb,prod,isqrt,factorial
from hashlib import sha256
import json,struct,time
import flint
from flint import nmod_poly

N=6;P=1000033
assert all(P%j for j in range(2,isqrt(P)+1))
e=0
while 2**e<4*N*(N-1):e+=1
D=2*(N-1)*e
assert (e,D)==(7,70)

def tables(d,last):
    assert last<P
    iv=[0]*(last+1);iv[1]=1
    for j in range(2,last+1):iv[j]=-(P//j)*iv[P%j]%P
    rise=[0]*(last+1);q=prod(range(1,d+1))%P;rise[1]=q
    for j in range(1,last-d+1):
        q=q*(j+d)%P*iv[j]%P;rise[j+1]=q
    return iv,rise

def coeff(a,d,m,t):
    iv,rise=t;v=1;c=[1]
    for k in range(m):
        v=v*(m-k)%P*iv[k+1]%P*rise[a*m+1+d*k]%P
        c.append(v)
    return c

start=time.monotonic();controls=0;stream=sha256();rows=[];total=0;maxarg=0
for d in range(2,D):
    lasta=d*(2*N-3)-1;limit=(lasta+d)*N;t=tables(d,limit)
    for a in range(1,lasta+1):
        c5=coeff(a,d,5,t);c6=coeff(a,d,6,t)
        assert len(c5)==6 and len(c6)==7 and c5[0]==c6[0]==1 and c5[-1] and c6[-1]
        f5=nmod_poly(c5,P);f6=nmod_poly(c6,P)
        assert f5.degree()==5 and f6.degree()==6
        g,s,q=f6.xgcd(f5)
        assert g==nmod_poly([1],P) and s*f6+q*f5==g
        stream.update(struct.pack('<IIII',a,d,P,0));total+=1;maxarg=max(maxarg,limit)
        if a in {1,lasta}:
            for m,c in [(5,c5),(6,c6)]:
                literal=[comb(m,k)*prod(range(a*m+1,a*m+d*k+1))%P for k in range(m+1)]
                assert c==literal;controls+=1
    rows.append({'gap':d,'first_lower_exponent':1,'last_lower_exponent':lasta,'cases':lasta})
assert total==21658==9*(69*70//2-1)-68 and len(rows)==68 and maxarg==4134

# Literal ordinary moments at both actual physical shapes, not an auxiliary
# discriminant or a substitution of one fixed Gamma shape into both orders.
for a,d in [(1,2),(2,3)]:
    t=tables(d,(a+d)*6)
    for m in [5,6]:
        c=coeff(a,d,m,t);A=F(2,3);B=F(-3,5)
        raw=sum(F(comb(m,k))*A**(m-k)*B**k*factorial(a*m+d*k)for k in range(m+1))
        # Recompute the normalized polynomial with literal integer coefficients.
        exact=[comb(m,k)*prod(range(a*m+1,a*m+d*k+1))for k in range(m+1)]
        assert raw==factorial(a*m)*A**m*sum(F(v)*(B/A)**k for k,v in enumerate(exact))
f=nmod_poly([1,1],P);g,s,t=f.xgcd(f*f)
assert g.degree()==1 and s*f+t*f*f==g
assert prod(range(3,5))!=prod(range(2,4))

out={'status':'PASS','terminal':N,'binary_exponent':e,'gap_cutoff':D,
     'prime':P,'tested_cases':total,'full_polynomial_Bezout_identities':total,
     'all_full_degrees':[5,6],'maximum_factorial_argument':maxarg,
     'coefficient_construction':'c0=1; c(k+1)=c(k)*(m-k)/(k+1)*(a*m+1+d*k)_(d), inverse-table and sliding rising products; no factorial tables',
     'actual_physical_shapes':['5*a+1','6*a+1'],
     'literal_integer_coefficient_controls':controls,'gaps':rows,
     'case_stream_sha256':stream.hexdigest(),'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
     'versions':{'python_flint':flint.__version__,'FLINT':flint.__FLINT_VERSION__},
     'seconds':time.monotonic()-start,
     'scope':'Completed finite actual terminal-six triangle; universal all-exponent terminal-six conclusion additionally uses retained high-gap/large-lower-exponent complements. No all-terminal conclusion.'}
Path('/tmp/independent_terminal_six_triangle_replay.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items()if k!='gaps'},indent=2))
