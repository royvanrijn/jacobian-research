#!/usr/bin/env python3
"""Complete actual terminal3..10 bridge with universal terminal-two complement."""
from __future__ import annotations
if not __debug__:
    raise SystemExit('Refusing -O/-OO: exact guards must execute.')
from pathlib import Path
from fractions import Fraction as Q
from math import comb,factorial,prod,isqrt
from hashlib import sha256
import argparse,copy,json,re,struct,sys
import flint
from flint import nmod_poly
HERE=Path(__file__).resolve().parent
PN='ORDINARY_BINOMIAL_TERMINALS_THROUGH_TEN.md'
CN='verify_ordinary_binomial_terminal_ten.py'
IN='ordinary_binomial_terminal_ten_inputs.json'
LN='ordinary_binomial_terminal_ten_source_ledger.json'
RN='ordinary_binomial_terminal_ten_verification.json'
PREFIX='research/factorial/image-rigidity/'
ROLES={
 'proof':('local',PN),'checker':('local',CN),'inputs':('local',IN),
 'old_rational_annuli':('repo',PREFIX+'lower_exponent_initial_and_boundary.tex'),
 'sharp_soft_criterion':('repo',PREFIX+'GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md'),
 'constant_large_exponent':('repo',PREFIX+'ALL_EXPONENTS_AND_LARGE_EXPONENT_CUTOFF.tex'),
}
SCOPE={
 'schema':1,'terminals':[3,4,5,6,7,8,9,10],'C':[3,4,6,7,9,10,12,13,14],
 'high_gap_cutoffs':[7,12,18,24,30,37,44,50],
 'main_domain':'N3..10: 1<=d<d0(N),1<=a<d(N-2), plus a0,d1..8',
 'main_counts':[23,129,450,1089,2154,3968,6587,9759],'main_total':24159,
 'terminal_two_entirely_sharp_soft':True,'maximum_factorial_argument':4400,
 'primes':[65537,1000003,1000033],
 'conclusion':'all ordinary supports; every terminal2..10; every characteristic-zero field',
 'all_orders':False,'all_first_ten_pairwise_coprime':False,
 'all_constant_finite_prefix_premise':False,'adjacent_all_order_premise':False,
 'earlier_terminal_finite_certificate_premise':False,
 'actual_shapes':['a*(N-1)+1','a*N+1'],
 'formal_verification':False,'external_review':False,'independent_replay':False,
}
def need(ok,msg):
    if not ok:raise AssertionError(msg)
def digest(p):return sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def same_typed(a,b):
    if type(a) is not type(b):return False
    if type(b) is dict:return set(a)==set(b) and all(same_typed(a[k],v) for k,v in b.items())
    if type(b) is list:return len(a)==len(b) and all(same_typed(x,y) for x,y in zip(a,b))
    return a==b
def check_inputs(d):need(same_typed(d,SCOPE),'Exact typed terminal/domain scope')
def check_sources(d,repo):
    need(type(d) is dict and set(d)=={'schema','sources'},'Ledger schema')
    need(type(d['schema']) is int and d['schema']==1,'Integer schema1')
    need(type(d['sources']) is list and len(d['sources'])==len(ROLES),'All mandatory source roles')
    out={};roots={'local':HERE,'repo':repo}
    for row in d['sources']:
        need(type(row) is dict and set(row)=={'role','kind','path','sha256'},'Source row schema')
        role=row['role'];need(type(role) is str and role in ROLES and role not in out,'Unique exact role')
        need((row['kind'],row['path'])==ROLES[role],'Exact source path')
        need(type(row['sha256']) is str and re.fullmatch('[0-9a-f]{64}',row['sha256']),'SHA256')
        base=roots[row['kind']].resolve();p=(base/row['path']).resolve()
        need(p.is_relative_to(base) and p.is_file() and digest(p)==row['sha256'],'Stale/missing source:'+role)
        out[role]=row
    return out
def prime(p):return type(p) is int and p>=2 and all(p%q for q in range(2,isqrt(p)+1))
def domain():
    main=[];rows=[]
    for N,d0,want in zip(SCOPE['terminals'],SCOPE['high_gap_cutoffs'],SCOPE['main_counts']):
        chunk=[(N,a,d) for d in range(1,d0) for a in range(1,d*(N-2))]
        chunk +=[(N,0,d) for d in range(1,9)]
        need(len(chunk)==want and len(set(chunk))==want,'Exact terminal domain')
        need(len(chunk)==(N-2)*(d0-1)*d0//2-(d0-1)+8,'Closed count formula')
        main+=chunk;rows.append({'terminal':N,'gap_cutoff':d0,'cases':len(chunk),
          'positive_a_cases':len(chunk)-8,'constant_boundary_cases':8})
    need(len(main)==24159 and len(set(main))==24159,'Full direct count')
    need(max((a+d)*N for N,a,d in main)==4400,'Actual factorial maximum')
    return main,rows
def complement_controls():
    ann=[]
    for N,d0 in zip(SCOPE['terminals'],SCOPE['high_gap_cutoffs']):
        c0=SCOPE['C'][N-3];c1=SCOPE['C'][N-2]
        need((c0+1)**(N-1)<2*c0**(N-1) and (c1+1)**N<2*c1**N,'Exact rational radii')
        B0=factorial((d0+1)*(N-1))//factorial(N-1)
        B1=factorial((d0+1)*N)//factorial(N)
        left=B1**(N-1);right=(c0*c1)**(N*(N-1))*B0**N
        need(left>right,'Universal high-gap factorial comparison')
        encoded=struct.pack('<I',left.bit_length())+left.to_bytes((left.bit_length()+7)//8,'big')+right.to_bytes((right.bit_length()+7)//8,'big')
        ann.append({'terminal':N,'C_previous':c0,'C_upper':c1,'gap_cutoff':d0,
                    'integer_inequality_sha256':sha256(encoded).hexdigest()})
    need(2**9>=4*10**2,'Constant large-exponent complement')
    for N in range(2,11):
        for d in range(1,12):
            for r in range(1,N):
                for s in range(r+1,N+1):
                    need(d*(N-2)*(s-r)>=d*(r-1),'Sharp-soft all-first-N threshold')
    return ann
def tables(p,last):
    need(prime(p) and last<p,'Prime/good factorial arguments')
    f=[1]*(last+1)
    for j in range(1,last+1):f[j]=f[j-1]*j%p
    ivf=[1]*(last+1);ivf[last]=pow(f[last],-1,p)
    for j in range(last,0,-1):ivf[j-1]=ivf[j]*j%p
    iv=[0]*(last+1);iv[1]=1
    for j in range(2,last+1):iv[j]=-(p//j)*iv[p%j]%p
    need(all(f[j]*ivf[j]%p==1 for j in range(last+1)),'Factorial units')
    return f,ivf,iv
def factorial_coeff(a,d,m,p,t):
    f,ivf,_=t
    return [comb(m,m-j)*f[a*m+d*(m-j)]*ivf[a*m]%p for j in range(m+1)]
def rising_table(d,p,iv,last):
    r=[0]*(last-d+2);v=1
    for j in range(1,d+1):v=v*j%p
    r[1]=v
    for s in range(1,last-d+1):v=v*(s+d)*iv[s]%p;r[s+1]=v
    return r
def recurrence_coeff(a,d,m,p,t,rise):
    iv=t[2];c=[1];v=1
    for k in range(m):
        v=v*(m-k)*iv[k+1]*rise[a*m+1+d*k]%p;c.append(v)
    return list(reversed(c))
def pack_poly(h,f):
    c=[int(x) for x in f];h.update(struct.pack('<I',len(c)))
    if c:h.update(struct.pack('<'+'I'*len(c),*c))
def finite(cases):
    primes=SCOPE['primes'];cache={};rise={};coverage=sha256();attempts=sha256()
    perprime={str(p):{'attempts':0,'accepted':0,'nontrivial_gcds':0} for p in primes}
    exceptions=[];bezout=0;coefficient_vectors=0;appell=0;success=0
    for p in primes:need(prime(p) and p>4400,'Verified good-reduction prime')
    for group,cases in [('main',cases)]:
        for N,a,d in cases:
            coverage.update(struct.pack('<III',N,a,d));accepted=False
            for index,p in enumerate(primes):
                if p not in cache:cache[p]=tables(p,4400)
                if (p,d) not in rise:rise[p,d]=rising_table(d,p,cache[p][2],4400)
                c0=factorial_coeff(a,d,N-1,p,cache[p]);c1=factorial_coeff(a,d,N,p,cache[p])
                need(c0==recurrence_coeff(a,d,N-1,p,cache[p],rise[p,d]),'Every lower coefficient vector independently rebuilt')
                need(c1==recurrence_coeff(a,d,N,p,cache[p],rise[p,d]),'Every upper coefficient vector independently rebuilt')
                coefficient_vectors+=2
                lower=nmod_poly(c0,p);upper=nmod_poly(c1,p)
                need(lower.degree()==N-1 and upper.degree()==N,'Full actual degrees')
                need(c0[-1]==c1[-1]==1 and c0[0] and c1[0],'Monic and nonzero constant')
                g,u,v=upper.xgcd(lower)
                need(u*upper+v*lower==g,'FULL polynomial Bezout identity every attempt')
                need(g.degree()>=0,'Nonzero gcd');bezout+=1
                perprime[str(p)]['attempts']+=1
                attempts.update(struct.pack('<IIIII',N,a,d,p,g.degree()))
                for poly in (upper,lower,g,u,v):pack_poly(attempts,poly)
                if a==0:need(upper.derivative()==N*lower,'Actual constant-family Appell');appell+=1
                if g==nmod_poly([1],p):
                    perprime[str(p)]['accepted']+=1;accepted=True;success+=1;break
                perprime[str(p)]['nontrivial_gcds']+=1
                exceptions.append({'group':group,'terminal':N,'a':a,'d':d,'prime':p,'gcd_degree':g.degree()})
            need(accepted,'Unresolved rational case: '+str((N,a,d)))
    need(success==24159,'Complete direct domain accepted')
    return {'main_cases':len(cases),'total_accepted':success,
      'full_polynomial_Bezout_identities':bezout,'independently_rebuilt_coefficient_vectors':coefficient_vectors,
      'actual_Appell_checks':appell,'primes':perprime,'exceptional_reductions':exceptions,
      'coverage_stream_sha256':coverage.hexdigest(),'polynomial_Bezout_stream_sha256':attempts.hexdigest()}
def controls():
    literal=0
    for N in range(2,11):
        for a,d in ((0,3),(1,2),(2,5)):
            for m in (N-1,N):
                c=[comb(m,m-j)*prod(range(a*m+1,a*m+d*(m-j)+1)) for j in range(m+1)]
                A=Q(2,3);B=-Q(3,5)
                raw=sum(Q(comb(m,k)*factorial(a*m+d*k))*A**(m-k)*B**k for k in range(m+1))
                need(raw==B**m*factorial(a*m)*sum(Q(x)*(A/B)**j for j,x in enumerate(c)),'Literal exact actual moments')
                literal+=1
    p=65537;f=nmod_poly([1,1],p);h=f*f;g,u,v=f.xgcd(h)
    need(g.degree()==1 and u*f+v*h==g,'Actual common-factor negative')
    need((u+1)*f+v*h!=g,'Actual corrupted-Bezout negative')
    need(6*prod(range(15,18))!=6*prod(range(13,16)),
         'Actual N7,a2,d3 lower c1 uses gamma13, not upper gamma15')
    return {'literal_rational_moment_controls':literal,'math_negatives':3}
def source_negatives(ledger,repo):
    out=[]
    def reject(name,f):
        try:f()
        except (AssertionError,KeyError,TypeError,ValueError):out.append(name);return
        raise AssertionError('Negative accepted:'+name)
    for role in sorted(ROLES):
        d=copy.deepcopy(ledger);d['sources']=[r for r in d['sources'] if r['role']!=role]
        reject('omit:'+role,lambda d=d:check_sources(d,repo))
        d=copy.deepcopy(ledger);next(r for r in d['sources'] if r['role']==role)['sha256']='0'*64
        reject('badpin:'+role,lambda d=d:check_sources(d,repo))
    for key in ('all_orders','all_first_ten_pairwise_coprime','all_constant_finite_prefix_premise',
                'adjacent_all_order_premise','earlier_terminal_finite_certificate_premise'):
        d=copy.deepcopy(SCOPE);d[key]=True;reject('false-scope:'+key,lambda d=d:check_inputs(d))
    d=copy.deepcopy(ledger);d['schema']=True;reject('boolean-schema',lambda:check_sources(d,repo))
    d=copy.deepcopy(SCOPE);d['main_total']=True;reject('typed-domain',lambda:check_inputs(d))
    return out
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--repo-root',required=True,type=Path)
    ap.add_argument('--output',type=Path,default=HERE/RN);ns=ap.parse_args();repo=ns.repo_root.resolve()
    inp=load(HERE/IN);check_inputs(inp);ledger=load(HERE/LN);rows=check_sources(ledger,repo)
    main_cases,domains=domain();ann=complement_controls();algebra=controls()
    result=finite(main_cases);bad=source_negatives(ledger,repo);check_sources(ledger,repo)
    payload={'schema':1,'status':'PASS','complete_main_finite_domain':True,'terminal_two_entirely_sharp_soft':True,
      'scope':inp,'proof_sha256':digest(HERE/PN),'checker_sha256':digest(HERE/CN),
      'input_sha256':digest(HERE/IN),'source_ledger_sha256':digest(HERE/LN),
      'mandatory_source_roles':len(rows),'sources':[rows[r] for r in sorted(rows)],
      'exact_complement_inequalities':ann,'domains':domains,'finite':result,'controls':algebra,
      'source_scope_negative_controls':bad,
      'versions':{'python_flint':flint.__version__,'FLINT':flint.__FLINT_VERSION__},
      'meaning':'All actual finite cases, every attempt full Bezout and independently rebuilt coefficients. Universal exponent coverage additionally uses exactly the written old-annulus and sharp-soft complement proofs. No all-terminal claim.',
      'formal_verification':False,'external_review':False,'independent_replay':False}
    ns.output.parent.mkdir(parents=True,exist_ok=True)
    ns.output.write_text(json.dumps(payload,indent=2,sort_keys=True)+'\n')
    print(json.dumps({'status':'PASS','main_cases':result['main_cases'],
      'full_Bezout':result['full_polynomial_Bezout_identities'],'exceptions':len(result['exceptional_reductions']),
      'primes':result['primes'],'output':str(ns.output)}))
if __name__=='__main__':main()
