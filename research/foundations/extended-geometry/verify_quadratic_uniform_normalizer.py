#!/usr/bin/env python3
"""Bounded exact replay; the written proof supplies the universal theorem."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
from math import comb, factorial
from pathlib import Path, PurePosixPath
import json
import sys

HERE=Path(__file__).resolve().parent
PROOF='QUADRATIC_COMPLEMENT_UNIFORM_NORMALIZER_COUNTEREXAMPLES.md'
INPUT='quadratic_uniform_normalizer_inputs.json'
LEDGER='quadratic_uniform_normalizer_source_ledger.json'
PINS={
 PROOF:'69582d9ee0446768701cb3fb0375ffefd8e8631e1382c2923bfcbfa28e0ae82c',
 'sources/long-spherical-e9edac2.pdf':'b2bf25fbc2d6211333d4ff3e1f6de0b6562990e869c6078a2e219ef746f1026b',
 'sources/milne-algebraic-groups-2017.pdf':'503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67',
 'REDUCTIVE_HOMOGENEOUS_AUTOMORPHISM_INVARIANT_TRACE.md':'b4b223a00cd9ed00a3f3957749b488dc0daea32f6e3491cb047d6d6dcac0aba5'}
ROLES={PROOF:'written_universal_proof',Path(__file__).name:'bounded_exact_checker',INPUT:'typed_replay_scope',
 'sources/long-spherical-e9edac2.pdf':'classical_method_context_not_logical_premise',
 'sources/milne-algebraic-groups-2017.pdf':'same_field_complete_reducibility',
 'REDUCTIVE_HOMOGENEOUS_AUTOMORPHISM_INVARIANT_TRACE.md':'narrow_normalized_trace_alternate_route'}
SCOPE={'characteristic':0,'dimensions':'all_N_ge_3','forms':'all_nondegenerate_over_given_field',
 'extensions':'all_characteristic_zero_extensions','requires_rational_point':False,'rank_A':2,'rank_B':'N-2',
 'phase':'(q+1)*(q-(q+1)^2*B)^2/q^3','marker':'q^(N-2)','nonpower_phase':'P^(N-2)*(1+P)',
 'absolute_nonpower_up_to_units':True,'ternary_direct_phase':True,'global_failure_descent_claimed':False,
 'sharp_pure_power_iff':'a>=N-2','eventual_integer_marker_iff':'ell>=N-2','all_orders_power_marker_iff':'N-2<=ell<=a',
 'exact_mixed_formula_domain':'k>=1,ell>=-k','long_is_logical_premise':False,
 'circuit_coefficients':[2,2,2], 'probe_dimensions':[3,4,5,6,8,9], 'probe_orders':[1,2,3,4,5,6,7,8,9,10],
 'assurance_flags':{'external_review':False,'formal_verification':False,'independent_replay':False,'worldwide_priority':False}}

def require(ok,msg):
 if not ok: raise ValueError(msg)
def unique(pairs):
 d={}
 for k,v in pairs:
  require(k not in d,'duplicate JSON key: '+k);d[k]=v
 return d
def readjson(path):
 return json.loads(path.read_text(encoding='utf8'),object_pairs_hook=unique,parse_constant=lambda x: (_ for _ in ()).throw(ValueError('nonfinite JSON')))
def typed_equal(a,b):
 if type(a) is not type(b): return False
 if isinstance(a,dict): return a.keys()==b.keys() and all(typed_equal(a[k],b[k]) for k in a)
 if isinstance(a,list): return len(a)==len(b) and all(typed_equal(x,y) for x,y in zip(a,b))
 return a==b
def digest(path): return sha256(path.read_bytes()).hexdigest()
def local(name):
 pp=PurePosixPath(name)
 require(type(name) is str and name==str(pp) and not pp.is_absolute() and '..' not in pp.parts and '\\' not in name,'nonlocal source path')
 path=HERE.joinpath(*pp.parts)
 require(path.resolve().is_relative_to(HERE),'source escapes package')
 require(path.is_file(),'missing source: '+name)
 return path

def add(a,b):
 r=a.copy()
 for e,c in b.items():
  r[e]=r.get(e,F(0))+c
  if not r[e]: del r[e]
 return r
def scale(a,c): return {e:v*c for e,v in a.items() if v*c}
def mul(a,b):
 r={}
 for (u,s),x in a.items():
  for (v,t),y in b.items():
   e=(u+v,s+t);r[e]=r.get(e,F(0))+x*y
 return {e:c for e,c in r.items() if c}
def mon(u,s=0,c=1): return {(u,s):F(c)}
def power(a,n):
 r=mon(0)
 for _ in range(n): r=mul(r,a)
 return r
def trace(a,N,ell=0):
 M=N-2
 return sum((c*F(M,M+2*s) for (u,s),c in a.items() if u+ell==0),F(0))
def jvalue(M,k):
 return sum((F(M*(-1)**v*comb(2*k,v),M+2*v) for v in range(2*k+1)),F(0))
def gbinom(n,r):
 if r<0:return F(0)
 v=F(1)
 for i in range(r):v*=F(n-i,i+1)
 return v
def endpoint(N,k,ell):
 require(ell>=-k,'endpoint formula outside its proved domain')
 return jvalue(N-2,k)*gbinom(k-(N-2),k-ell)

def frac(x): return str(x.numerator)+'/'+str(x.denominator)

def replay(scope):
 one,U,S=mon(0),mon(1),mon(0,1)
 T=add(one,U)
 V=mul(mon(-1),power(add(one,scale(S,-1)),2))
 R=mul(mul(S,add(U,scale(one,2))),add(scale(one,2),scale(mul(S,add(add(power(U,2),scale(U,2)),scale(one,2))),-1)))
 circuit=mul(T,add(V,scale(R,-1)))
 square=mul(mul(mon(-1),T),power(add(one,scale(mul(power(T,2),S),-1)),2))
 require(circuit==square,'regular circuit identity')
 radial={(u-s,s):v for (u,s),v in circuit.items()}
 literal=mul(mul(mon(-3),T),power(add(U,scale(mul(power(T,2),mon(0,1)),-1)),2))
 require(radial==literal,'radial q,B identity')
 powers=[one]
 maxk=max(max(scope['probe_orders']),2*(max(scope['probe_dimensions'])-2)+2)
 for k in range(1,maxk+1):powers.append(mul(powers[-1],circuit))
 rows=[];comparisons=0
 for N in scope['probe_dimensions']:
  M=N-2
  for k in scope['probe_orders']:
   for ell in sorted(set((-k,1-k,0,1,M,k,k+1))):
    require(trace(powers[k],N,ell)==endpoint(N,k,ell),'direct expansion differs from endpoint');comparisons+=1
   require((trace(powers[k],N)==0)==(k>=M),'sharp pure-power finite model')
   if k>=M:
    require(trace(powers[k],N,M)==jvalue(M,k)>0,'fixed marker')
    for ell in range(M+3):require((trace(powers[k],N,ell)!=0)==(M<=ell<=k),'sharp allorders-marker finite model')
   rows.append({'N':N,'k':k,'pure_moment':frac(trace(powers[k],N)),'marked_moment':frac(trace(powers[k],N,M))})
  for m in (1,2):
   expanded={}
   for r in range(m+1):expanded=add(expanded,scale(powers[M*m+r],comb(m,r)))
   expected=sum((comb(m,r)*jvalue(M,M*m+r) for r in range(m+1)),F(0))
   require(trace(expanded,N)==0 and trace(expanded,N,M)==expected and expected>0,'nonpower phase moments')
 # At B=0, P=(q+1)/q, so the minus-one level and transverse derivative are rational identities.
 require({(u,0):v for (u,s),v in literal.items() if s==0}==add(one,mon(-1)),'simple-point radial restriction')
 q=F(-1,2)
 require(1+1/q==-1 and -1/q**2==-4,'simple point value/derivative')
 require(sum((v*(-1)**u for (u,s),v in square.items() if s==0),F(0))==0,'ternary q+1 factor')
 # A wrong phase is detected mathematically, and an invalid endpoint extrapolation is excluded.
 require(trace(add(circuit,one),3)==1,'bad coefficient negative control')
 require(trace(powers[1],3,-2)==F(4,3) and trace(powers[1],3,-2)!=jvalue(1,1)*gbinom(0,3),'negative-marker domain control')
 return {'direct_endpoint_comparisons':comparisons,'direct_nonpower_phase_models':2*len(scope['probe_dimensions']),
 'regular_factor_and_simple_point_identities':True,'bad_phase_negative_control':True,
 'outside_endpoint_domain_counterexample':'N3,k1,ell-2:4/3','rows':rows}

def main():
 require(__debug__,'optimized Python is unsupported')
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root',type=Path,default=HERE);ap.add_argument('--output',type=Path,default=HERE/'quadratic_uniform_normalizer_verification.json');args=ap.parse_args()
 require(args.repo_root.is_dir(),'repo-root must exist')
 scope=readjson(local(INPUT));require(typed_equal(scope,SCOPE),'typed replay scope drift')
 ledger=readjson(local(LEDGER));require(set(ledger)=={'schema','entries'} and type(ledger['schema']) is int and ledger['schema']==1,'ledger schema')
 require(type(ledger['entries']) is list and len(ledger['entries'])==len(ROLES),'ledger rows')
 seen=set();authenticated={}
 for row in ledger['entries']:
  require(type(row) is dict and set(row)=={'path','sha256','role'},'ledger entry schema')
  name=row['path'];require(type(name) is str and name in ROLES and name not in seen,'unexpected/duplicate source role')
  require(row['role']==ROLES[name],'source role drift');path=local(name)
  require(type(row['sha256']) is str and digest(path)==row['sha256'],'source body drift: '+name)
  if name in PINS:require(row['sha256']==PINS[name],'pinned source drift: '+name)
  authenticated[name]=row['sha256'];seen.add(name)
 require(seen==set(ROLES),'source inventory')
 out=args.output.resolve()
 for name in list(ROLES)+[LEDGER]:
  src=local(name)
  require(out!=src.resolve() and not (out.exists() and out.samefile(src)),'output would overwrite protected source')
 require(out.parent.is_dir(),'output parent must exist')
 result={'schema':1,'status':'PASS','scope':scope,'authenticated_sources':authenticated,'bounded_models':replay(scope),
 'universal_proof':PROOF,'computation_is_not_universal_proof':True,'assurance_flags':SCOPE['assurance_flags']}
 out.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n',encoding='utf8')
 print('PASS: bounded exact replay; universal claims are in '+PROOF)
if __name__=='__main__':
 try:main()
 except (ValueError,OSError,TypeError,KeyError,json.JSONDecodeError) as exc:
  print('FAIL: '+str(exc),file=sys.stderr);sys.exit(1)
