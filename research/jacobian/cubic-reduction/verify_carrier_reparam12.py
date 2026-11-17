#!/usr/bin/env python3
"""Independent Fraction audit of the47-coefficient cubic carrier realization.
No SymPy, FLINT or Singular is imported. Completeness uses the stated elementary
algebraic elimination; it is not proof-assistant formalization.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement
from collections import Counter
import hashlib,json,time
from verify import add,mul,variable,derivative,decode
from verify_power_single_compress import scale
from weighted_verify_utils import expression,substitution,coefficient_buckets,equal_generators,normal
from verify_weighted_unmarked7 import laurent
HERE=Path(__file__).resolve().parent
ROOT=HERE/'carrier-reparam-wide'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def const(c,n):return {(0,)*n:Q(c)} if c else {}
def main():
 started=time.monotonic();source=HERE/'cubic12-integer-witness.json';old=load(source);base_validation=load(HERE/'integer-validation12.json');assert base_validation['status']=='pass' and base_validation['witness_sha256']==sha(source) and base_validation['determinant_one'] and base_validation['identity_derivative'];F=[decode(p,12) for p in old['coefficients']];assert sum(sum(e)>=2 for f in F for e in f)==36
 # All-degree one-axis nonlinear source shear obstruction. Subsequent changes
 # in target coordinate i do not affect the listed off-axis witness row.
 witnesses={0:5,1:4,2:3,3:1,4:2,5:0,6:0,7:1,8:2,9:1,10:2,11:1};records=[]
 assert F[5]==expression('x6+x1**2',['x'+str(i) for i in range(1,13)]) and F[4]==expression('x5+x2**3',['x'+str(i) for i in range(1,13)])
 for i,row in witnesses.items():
  assert i!=row
  if i>=2:
   D=derivative(F[row],i);assert max(sum(e) for e in D)==2 and not derivative(D,i);assert not any(e[i]>1 for e in F[row])
  records.append(dict(source_coordinate=i+1,unmodified_output=row+1,degree_after_nonlinear_shear='2*d' if i==0 else ('3*d' if i==1 else 'd+2')))
 ip=ROOT/'input.json';d=load(ip);assert d['source_sha256']==sha(HERE/'carrier_reparam_wide.py') and d['input_sha256']==sha(source);Z=d['variables'];assert Z==['z'+str(i) for i in range(1,48)];xn=['x'+str(i) for i in range(1,13)];names=xn+Z;N=len(names);atoms=d['atoms'];weights=[1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0];original={(i,e):c for i in range(3,12) for e,c in F[i].items() if sum(e)>=2};expected=set(original)
 for i in range(3,12):
  for a in range(4):
   for b in range(4):
    for c in range(4):
     if 2<=a+b+c<=3 and a-b-2*c==weights[i]:expected.add((i,(a,b,c)+(0,)*9))
  for j,k in combinations_with_replacement(range(i),2):
   if weights[j]+weights[k]==weights[i]:expected.add((i,tuple(int(a==j)+int(a==k) for a in range(12))))
 assert {(a['row'],tuple(a['exponent'])) for a in atoms}==expected and len(atoms)==len(expected)==47
 new=[variable(i,N) for i in range(12)]
 for j,a in enumerate(atoms):
  assert a['parameter']==Z[j];e=tuple(a['exponent']);c=Q(a['base_coefficient']);assert c==original.get((a['row'],e),Q(1));assert not any(e[a['row']:]);term={e+tuple(int(k==j) for k in range(47)):c};new[a['row']]=add(new[a['row']],term)
 E=[expression(text,names) for text in d['source_substitution']];assert len(E)==12 and E[:3]==new[:3]
 for i in range(3,12):assert derivative(E[i],i)==const(1,N) and not any(any(e[i+1:12]) for e in E[i]);assert not any(sum(e[:12])<2 for e in add(E[i],variable(i,N),Q(-1)))
 G=[substitution(f,E,N) for f in F];assert G[3:]==new[3:];assert [expression(text,names) for text in d['carrier_outputs']]==new[3:]
 high=[];low=[]
 for i,f in enumerate(G):
  for e,p in coefficient_buckets(f,12).items():
   if sum(e)>3:high.append(p)
   elif sum(e)>=2:low.append((i,e,p))
 gen=[{tuple(e):Q(c) for e,c in terms} for terms in d['high_equations']];equal_generators(high,gen);claimed={(a['row'],tuple(a['exponent'])):expression(a['coefficient'],Z) for a in d['low_coefficients']};assert {(i,e):p for i,e,p in low}==claimed;assert len(gen)==93 and len(low)==96
 cp=ROOT/'classification-certificate.json';c=load(cp);assert c['source_sha256']==sha(HERE/'carrier_reparam_lift.py') and c['input_sha256']==sha(ip) and c['singular_source_sha256']==sha(ROOT/'classification-lift.sing') and c['variables']==Z
 active={1,2,3,4,5,6,7,8,9,10,11,12,13,17,18,22,38,42};texts=['z'+str(j) for j in range(1,48) if j not in active]+['z3-1','z6-1','z1-z2','z7-z9','z8-z9','z10-z12','z2*z9-1','z5*z12-z1','z4-z9*z12','z11-z4-3*z22','z5*z38-z17','z2*z17-2*z13','z2*z18-z5*z22','z5*z42+z22'];targets=[expression(text,Z) for text in texts];assert all(all(k>=0 for k in e) for poly in targets+gen for e in poly);assert c['target_texts']==texts and targets==[{tuple(e):Q(v) for e,v in terms} for terms in c['targets']];assert len(targets)==43;results=[{} for _ in targets];terms=0;maxdegree=0;seen=set()
 for entry in c['entries']:
  t,g=entry['target'],entry['generator'];assert (t,g) not in seen;seen.add((t,g));m={tuple(e):Q(v) for e,v in entry['multiplier']};assert all(len(e)==47 and all(k>=0 for k in e) for e in m);results[t]=add(results[t],mul(m,gen[g]));terms+=len(m);maxdegree=max(maxdegree,max(sum(e) for e in m))
 assert results==targets and terms==123;print('REPLAYED 43 MEMBERSHIPS',terms,'TERMS',flush=True)
 fp=ROOT/'family.json';family=load(fp);assert family['source_sha256']==sha(HERE/'carrier_reparam_family.py') and family['input_sha256']==sha(ip);pn=['p','q','s','t'];assert family['variables']==pn and family['nonzero']==['p','q'];values=[laurent(family['parameterization'][name],pn) for name in Z]
 explicit={1:'1/p',2:'1/p',3:'1',4:'p*q',5:'1/(p*q)',6:'1',7:'p',8:'p',9:'p',10:'q',11:'p*q+3*t',12:'q',13:'s/(2*p)',17:'s',18:'t/q',22:'t',38:'p*q*s',42:'-p*q*t'};assert values==[laurent(explicit.get(j,'0'),pn) for j in range(1,48)]
 assert not any(substitution(poly,values,4) for poly in targets+gen)
 reduced={(i,e):q for i,e,p in low if (q:=substitution(p,values,4))};assert reduced=={(a['row'],tuple(a['exponent'])):laurent(a['coefficient'],pn) for a in family['low_coefficients']};assert len(reduced)==45
 # p,q are nonzero. Classify coefficient factors after removing their Laurent
 # monomial units, retaining the exact s,t zero conditions.
 factors=Counter()
 for f in reduced.values():
  shift=[-min(e[i] for e in f) for i in range(2)]+[0,0];q={tuple(a+b for a,b in zip(e,shift)):v for e,v in f.items()};key=normal(q)
  if key==normal(expression('1',pn)):factors['unit']+=1
  elif key==normal(expression('s',pn)):factors['s']+=1
  elif key==normal(expression('t',pn)):factors['t']+=1
  elif key==normal(expression('3*p+2*s',pn)):factors['3p+2s']+=1
  elif key==normal(expression('p*q+3*t',pn)):factors['pq+3t']+=1
  else:raise AssertionError(q)
 assert factors==Counter({'unit':33,'s':4,'t':5,'3p+2s':1,'pq+3t':2})
 # At s=t=0 the entire map is the diagonal conjugate of KCUB3, not a new
 # stable-equivalence construction. Check every complete coefficient.
 D=[laurent(text,pn) for text in ['1','1','1','p','1','1/(p*q)','p*q','1','1/p','1/p','1/(p*q)','1/(p*q)']];Dinv=[laurent(text,pn) for text in ['1','1','1','1/p','1','p*q','1/(p*q)','1','p','p','p*q','p*q']];at0=[variable(0,4),variable(1,4),{},{}];base={(i,e):{ex:c0 for ex,c0 in co.items() if ex[2:]==(0,0)} for (i,e),co in reduced.items()};base={k:v for k,v in base.items() if v};diag={}
 for i,f in enumerate(F):
  for e,c0 in f.items():
   if sum(e)<2:continue
   co=scale(Dinv[i],c0)
   for j,k in enumerate(e):
    for a in range(k):co=mul(co,D[j])
   diag[i,e]=co
 assert base==diag and len(base)==36
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),laurent_utility_sha256=sha(HERE/'verify_weighted_unmarked7.py'),base_witness_sha256=sha(source),base_validation_sha256=sha(HERE/'integer-validation12.json'),input_sha256=sha(ip),certificate_sha256=sha(cp),family_sha256=sha(fp),parameters=47,high_equations=93,low_slots_before_reduction=96,exact_memberships=43,multiplier_monomials=terms,max_multiplier_degree=maxdegree,all_memberships_multiplied=True,all_source_substitutions_triangular=True,entire_carrier_composition_reconstructed=True,four_parameter_family_sufficient=True,family_nonzero_parameters=['p','q'],family_free_parameters=['s','t'],generic_nonlinear_total=45,coefficient_factor_counts=dict(factors),exact_support_formula='36+4*[s!=0]+5*[t!=0]-[s=-3*p/2]-2*[t=-p*q/3]',minimum_nonlinear_total=36,minimizers='s=t=0,p*q!=0',minimum_family_diagonal_conjugacy_checked=True,one_axis_arbitrary_nonlinear_source_obstructions=records,all_degree_one_axis_obstruction=True,ordinary_degree=3,dimension=12,scope='Fixed KCUB3 map with all47 specified weighted triangular carrier-output coefficients varied; C_old(E)=C_new and G=K E. All cubic members are the written four-parameter family. Its minimum displayed nonlinear count is36 and minimizers are diagonal conjugates of KCUB3. Separately any nonlinear one-axis elementary source shear (with any subsequent correction of that same target coordinate) leaves some unchanged output of degree>3. No multistage arbitrary transformation, affine-orbit optimum, other carrier graph/order/support or global simplicity/dimension optimum is asserted.',written_elimination_formalized=False,seconds=round(time.monotonic()-started,3))
 (ROOT/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
