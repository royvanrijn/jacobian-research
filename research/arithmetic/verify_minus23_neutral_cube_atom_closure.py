#!/usr/bin/env python3
"""Exact finite bridge and rational constants for the nonprincipal cube atoms.

The certificate checks literal witnesses, not the GP search algorithm. The
infinite conclusion also uses the written retained Hecke/local proof.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import gcd, isqrt, prod
from pathlib import Path
import sys

if not __debug__:
    raise SystemExit('Run without -O: exact assertion checks must remain enabled.')

HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--repo-root',type=Path)
parser.add_argument('--output',type=Path,default=HERE/'minus23_neutral_cube_atom_closure_verification.json')
args=parser.parse_args()
root=args.repo_root.resolve() if args.repo_root else next((p for p in [HERE,*HERE.parents] if(p/'research/MATH_STATUS.json').is_file()),None)
if root is None:raise SystemExit('Use --repo-root outside the repository.')

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def primes(bound):
    sieve=bytearray(b'\1')*(bound+1)
    sieve[:2]=b'\0\0'
    for p in range(2,isqrt(bound)+1):
        if sieve[p]:
            sieve[p*p:bound+1:p]=b'\0'*((bound-p*p)//p+1)
    return [p for p in range(2,bound+1)if sieve[p]]

def chi(p):
    if p==2:return 1
    if p==23:return 0
    return 1 if pow((-23)%p,(p-1)//2,p)==1 else -1

def Q(x,y):return x*x+x*y+6*y*y

def primitive_flags(bound):
    flags=set()
    ymax=isqrt(4*bound//23)
    for y in range(-ymax,ymax+1):
        smax=isqrt(4*bound-23*y*y)
        for s in range(-smax,smax+1):
            if (s-y)%2:continue
            x=(s-y)//2
            if gcd(x,y)==1:
                n=Q(x,y)
                if n>0:flags.add(n)
    return flags

packet_path=HERE/'minus23_neutral_cube_atom_witnesses.json'
packet=json.loads(packet_path.read_text())
for rel,pin in packet['retained_sources'].items():
    assert digest(root/rel)==pin,('stale retained source',rel)
bound=2**14
assert packet['prime_bound_exclusive']==bound
all_primes=primes(bound-1)
principal=primitive_flags(bound-1)
expected={p for p in all_primes if chi(p)==1 and p not in principal}
assert len(expected)==633
target_domain={(p,c)for p in expected for c in (5,10,15)}
rows=packet['witnesses']
assert len(rows)==len(target_domain)==1899
seen=set()
for p,c,x,y,u,v in rows:
    assert (p,c)not in seen
    seen.add((p,c))
    assert (p,c)in target_domain
    assert gcd(x,y)==gcd(u,v)==1
    a,b=Q(x,y),Q(u,v)
    assert min(a,b)>0 and a%23 and b%23
    assert a+b==c*p**3
assert seen==target_domain

# Negative controls: every domain row matters and imprimitive norms are refused.
assert {(p,c)for p,c,*_ in rows[:-1]}!=target_domain
p,c,x,y,u,v=rows[0]
assert Q(x,y)+Q(u+1,v)!=c*p**3
assert Q(2,0)==4 and gcd(2,0)!=1

# Complete divisor-envelope maximization, no floating-point inequalities.
euler_primes=primes(31)
exponents={p:max(range(33),key=lambda a:F((a+1)**5,p**a))for p in euler_primes}
assert list(exponents.values())==[6,4,2,2,1,1,1,1,1,1,1]
assert F(34,33)**5<2
divisor_product=prod(F((a+1)**5,p**a)for p,a in exponents.items())
assert divisor_product==F(16780370116608000,955049953)
assert divisor_product<29**5

# Exact short lattice table, with complete square-completion coordinate bounds.
short=[]
for n in range(1,25):
    total=0
    ymax=isqrt(4*n//23)
    for y in range(-ymax,ymax+1):
        smax=isqrt(4*n-23*y*y)
        for s in range(-smax,smax+1):
            if (s-y)%2==0:total+=1
    short.append(total-1)
    assert total-1<=2*n
assert short==[2,2,2,4,4,8,8,12,14,14,14,18,18,18,18,20,20,24,24,24,24,24,26,30]
assert F(19,4)**2<23
assert F(54,19)/F(90,133)==F(21,5)
assert F(21,5)**2<25

# Actual 2/3 baseline constants; exact exclusion of 5 saves its 1/10 mass.
baselines=[F(72,55),F(108,55),F(96,55)]
selected=[F(5,9),F(25,54),F(1,2)]
main=[b*(r-F(7,50))for b,r in zip(baselines,selected)]
assert main==[F(748,1375),F(872,1375),F(864,1375)]
assert all(v>F(1,2)for v in main)
assert F(6,25)-F(1,10)==F(7,50)

# All elementary onset comparisons are exact integer/rational statements.
assert 7*43>=15*20
assert 87**20<=2**129
assert 3*(2**15-1)>F(11,4)**11
assert F(3)*(F(1)-F(1,2**15))>=F(29,10)
cube_margin=F(87,2)-F(256,11)-F(640,33)-F(3584,2**15)
assert cube_margin==F(139,192)>0
assert 5*2**42>2**43

# e<87/32 from its factorial series and a geometric tail.
factorial=1
esum=F(0)
for k in range(7):
    if k:factorial*=k
    esum+=F(1,factorial)
assert esum+F(1,4410)<F(87,32)
assert F(87,32)**11<4*(2**14-1)
assert F(4)*(F(1)-F(1,2**14))>=F(39,10)
principal_margin=F(32)-F(512,33)-F(18560,1287)-F(9,128)
assert principal_margin==F(328385,164736)>0

ledger_path=HERE/'minus23_neutral_cube_atom_source_ledger.json'
ledger=json.loads(ledger_path.read_text())
primary_pins=[]
for entry in ledger['retained_imports']:
    assert digest(root/entry['path'])==entry['sha256'],entry['path']
for entry in ledger['primary_sources']:
    assert digest(HERE/entry['local_path'])==entry['sha256'],entry['local_path']
    primary_pins.append({'path':entry['local_path'],'sha256':entry['sha256'],'role':entry['role']})
for name,pin in ledger['packet_sha256'].items():
    assert digest(HERE/name)==pin,name

proof=HERE/'MINUS23_NEUTRAL_CUBE_ATOM_CLOSURE.md'
payload={
    'status':'PASS',
    'scope':'Complete literal finite nonprincipal-cube atom bridge plus exact constants/onset arithmetic. The infinite theorem uses the written retained Hecke and local operator inputs; no principal-atom completeness is claimed.',
    'prime_bound_exclusive':bound,
    'complete_nonprincipal_split_prime_count':len(expected),
    'literal_primitive_23_unit_decompositions':len(rows),
    'cube_target_onset_power_of_two':43,
    'principal_target_sufficient_onset_power_of_two':40,
    'cube_scaled_positive_surplus':str(cube_margin),
    'principal_scaled_positive_surplus':str(principal_margin),
    'euler_maximizing_exponents':{str(p):a for p,a in exponents.items()},
    'divisor_constant':29,
    'nonzero_binary_lattice_count_short_table':short,
    'main_term_constants':[str(v)for v in main],
    'negative_controls':['omitted_domain_row','changed_literal_coordinate','ordinary_imprimitive_four'],
    'source_sha256':{'proof':digest(proof),'checker':digest(Path(__file__).resolve()),'literal_packet':digest(packet_path)},
    'retained_source_sha256':packet['retained_sources'],
    'source_ledger_sha256':digest(ledger_path),
    'verified_primary_source_pins':primary_pins,
    'external_or_formal_review':False,
}
args.output.write_text(json.dumps(payload,indent=2)+'\n')
print(json.dumps(payload,indent=2))
