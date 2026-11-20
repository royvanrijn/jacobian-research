#!/usr/bin/env python3
"""Independent complete fixed-shape lower-one bridge, with checked Bezout identities.

Construct coefficients by a binomial/rising-block recurrence, not factorial
tables. The tested polynomial T_N(t)=sum binom(N,k)(N)_(dk)t^k is the
reciprocal of A_N(z;N). We explicitly preserve both degrees and its nonzero
constant coefficient. A Bezout identity for T_N and T_N' over a prime field
certifies squarefreeness over QQ. This is the FIXED-SHAPE auxiliary, not
the actual upper moment of shape N+1. The latter satisfies
N*A_N(z;N+1)=N*(d+1)*A_N(z;N)-d*z*A_N'(z;N), while
A_N'(z;N)=N*A_(N-1)(z;N). Thus auxiliary squarefreeness excludes a
common root of the two actual ordinary factorial moments.
"""
if not __debug__:
    raise SystemExit("Unoptimized Python required for exact checks")

from pathlib import Path
from math import isqrt, factorial, comb
from hashlib import sha256
from concurrent.futures import ProcessPoolExecutor, as_completed
import argparse, json, struct, time
import flint
from flint import nmod_poly

PRIMES=(1000033,1000037,1000039,1000003)

def prime(p):
    return p>1 and all(p%i for i in range(2,isqrt(p)+1))

def recurrence_tables(d,limit,p):
    assert 1<=d<=limit<p and prime(p)
    inv=[0]*(limit+1);inv[1]=1
    for j in range(2,limit+1):
        inv[j]=(-(p//j)*inv[p%j])%p
    rising=[0]*(limit+1)
    q=1
    for j in range(1,d+1):q=q*j%p
    rising[1]=q
    # rising[t]=t(t+1)...(t+d-1), for all required starting points.
    for t in range(1,limit-d+1):
        q=q*(t+d)%p*inv[t]%p;rising[t+1]=q
    return inv,rising

def coefficients(N,d,p,tables):
    inv,rising=tables
    c=[1];v=1
    for k in range(N):
        v=v*(N-k)%p*inv[k+1]%p*rising[N+d*k]%p
        c.append(v)
    return c

def worker(args):
    d,last=args;start=time.monotonic();limit=(d+1)*last-1
    cache={};counts={};fallbacks=[];stream=sha256()
    for N in range(2,last+1):
        attempts=[]
        for p in PRIMES:
            if p not in cache:cache[p]=recurrence_tables(d,limit,p)
            c=coefficients(N,d,p,cache[p]);assert len(c)==N+1 and c[0]==1 and c[-1]
            f=nmod_poly(c,p)
            # Independent literal differentiation, not f.derivative().
            df=nmod_poly([k*c[k]%p for k in range(1,N+1)],p)
            assert f.degree()==N and df.degree()==N-1
            g,s,t=f.xgcd(df)
            # Every full-size xgcd result is independently verified as a
            # polynomial identity, not accepted solely from the gcd call.
            assert s*f+t*df==g
            attempts.append([p,g.degree()])
            if g.degree()==0:
                assert g==nmod_poly([1],p)
                counts[str(p)]=counts.get(str(p),0)+1
                stream.update(struct.pack('<III',d,N,p))
                if len(attempts)>1:fallbacks.append({'terminal':N,'attempts':attempts})
                break
        else:raise AssertionError(('Unresolved required case',d,N,attempts))
    return {'gap':d,'first_terminal':2,'last_terminal':last,'tested':last-1,
        'verified_full_Bezout_identities':last-1+sum(len(v['attempts'])-1 for v in fallbacks),
        'successful_primes':counts,'fallbacks':fallbacks,
        'full_degree':last,'case_stream_sha256':stream.hexdigest(),
        'seconds':time.monotonic()-start}

def controls():
    count=0
    for d in (1,2,6,10,43):
        p=PRIMES[0];tt=recurrence_tables(d,10*(d+1),p)
        for N in range(2,10):
            c=coefficients(N,d,p,tt)
            assert c==[comb(N,k)*factorial(N+d*k-1)//factorial(N-1)%p for k in range(N+1)]
            # Literal actual-upper shape N+1 and the fixed-shape identity.
            upper=[comb(N,k)*factorial(N+d*k)//factorial(N)%p for k in range(N+1)]
            assert [N*v%p for v in upper]==[(N+d*k)*v%p for k,v in enumerate(c)]
            count+=1
    p=PRIMES[0];f=nmod_poly([1,2,1],p);g,s,t=f.xgcd(nmod_poly([2,2],p))
    assert g.degree()==1 and s*f+t*nmod_poly([2,2],p)==g
    assert s*f+(t+1)*nmod_poly([2,2],p)!=g
    return count

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--jobs',type=int,default=4)
    ap.add_argument('--output',type=Path,default=Path('/tmp/independent_lower_one_practical_finite_box.json'))
    ap.add_argument('--benchmark-last',type=int)
    a=ap.parse_args();start=time.monotonic();controls_count=controls()
    expected=[(d,max(5000,128*(d+1))-1)for d in range(10,44)]
    assert len(expected)==34 and sum(last-1 for d,last in expected)==171812
    domain=expected if a.benchmark_last is None else [(d,min(last,a.benchmark_last))for d,last in expected]
    rows=[]
    with ProcessPoolExecutor(max_workers=a.jobs)as pool:
        futures={pool.submit(worker,c):c for c in domain}
        for f in as_completed(futures):
            r=f.result();rows.append(r);print(json.dumps(r),flush=True)
            a.output.with_suffix('.running.json').write_text(json.dumps({'status':'RUNNING','completed_gaps':sorted(rows,key=lambda r:r['gap'])},indent=2)+'\n')
    rows.sort(key=lambda r:r['gap']);total=sum(r['tested']for r in rows)
    full=a.benchmark_last is None
    if full:
        assert [(r['gap'],r['last_terminal'])for r in rows]==expected
        assert total==171812
    out={'status':'PASS','complete_box':full,'tested_cases':total,
        'verified_full_Bezout_identities':sum(r['verified_full_Bezout_identities']for r in rows),
        'coefficient_construction':'c0=1; c(k+1)=c(k)*(N-k)/(k+1)*(N+d*k)_(d), using independently generated modular integer inverses and sliding rising products',
        'squarefreeness_certificate':'Full-degree modular Bezout identity s*T_N+t*T_N_prime=1 checked for every required case',
        'actual_physical_shape':'Upper actual shape N+1, previous actual shape N; fixed-shape auxiliary N links them by the exact derivative identity',
        'control_identities':controls_count,'flint_version':flint.__version__,
        'gaps':rows,'seconds':time.monotonic()-start,
        'checker_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    a.output.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'PASS','complete_box':full,'tested_cases':total,'verified_full_Bezout_identities':out['verified_full_Bezout_identities'],'seconds':out['seconds']}),flush=True)

if __name__=='__main__':main()
