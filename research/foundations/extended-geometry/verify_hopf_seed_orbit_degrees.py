#!/usr/bin/env python3
"""Exact polynomial/source support for three full Hopf automorphism-orbit spectra."""
if not __debug__:raise SystemExit('optimization unsupported: exact controls and source guards are required')
from pathlib import Path
from fractions import Fraction as F
from math import factorial,comb,prod
import argparse,copy,hashlib,json
HERE=Path(__file__).resolve().parent
PN='HOPF_KNOWN_SEED_AUTOMORPHISM_ORBIT_DEGREES.md';CN='verify_hopf_seed_orbit_degrees.py'
IN='hopf_seed_orbit_degree_inputs.json';LN='hopf_seed_orbit_degree_source_ledger.json'
PROOF='5056557fcfccc64d642b1baf6ad62e053c0a323de21d99abf22a69e26da6f01d'
PDF='quadric_orbit_sources/blanc_bot_poloni2023.pdf'
PRIMARY='91b7e296d7c04c3bcbcf550b19f5748d709295c479337819065aa94ff076438e'
BOUNDARY={'all_characteristic_zero_fields_and_extensions':True,'full_automorphism_orbit_degrees':True,
 'quadratic_case_primary_applies':True,'quartic_exact_even_degrees_at_least_four':True,
 'degree_six_exact_multiples_three_at_least_six':True,'degree_eight_exact_odd_times_first_shear_spectrum':True,
 'degree_eight_minimum_eight_first_prime_eleven':True,'degree_eight_prime_classes_two_mod_three_or_eight_mod_thirteen':True,
 'generic_H_polynomials_and_constant_T_absorption':True,'degree_growth_only_for_reduced_words':True,
 'pointwise_height_degree_lemma_requires_matching_top_height':True,'PGL2_frames_require_nonzero_determinant':True,
 'trace_or_purity_theorem_logical_premise':False,'old_seed_unsafe_theorem_logical_premise':False,
 'all_unsafe_phases_are_in_these_orbits':False,'global_degree_five_seven_exclusion':False,
 'unconditional_original_Image_counterexample_conclusion':False,'degree_eight_orbit_contains_prime_thirteen_or_nineteen':False,
 'quartic_or_degree_six_orbit_has_prime_degree':False,'every_nonlinear_word_increases_degree':False,
 'positive_characteristic_scope':False,'prime_density_claim':False,'formal_verification':False,
 'independent_replay':False,'external_review':False}
PRIMARY_METADATA={'title':'Real forms of some Gizatullin surfaces and Koras-Russell threefolds',
 'authors':['Jeremy Blanc','Anna Bot','Pierre-Marie Poloni'],'year':2023,'journal':'Publicacions Matematiques',
 'volume':67,'pages':'851-890','doi':'10.5565/PUBLMAT6722314','url':'https://algebra.dmi.unibas.ch/bot/publmat.pdf',
 'locator':'Theorems3.5-3.6 printed860-862; Lemma3.9 printed863-864',
 'role':'mathematical_quadratic_automorphism_generation_only'}
FRAMES=[
 {'name':'native','matrix':[['1','0','0'],['0','1','0'],['0','0','1']], 'heights':[2,3,3], 'last_top':['-1','-6','-6']},
 {'name':'opposite','matrix':[['0','1','0'],['1','0','0'],['0','0','1']], 'heights':[6,9,13], 'last_top':['-90',str(F(factorial(9),8)),str(F(factorial(13),32))]},
 {'name':'generic','matrix':[['1','-1/2','-1'],['0','1','0'],['0','1','1']], 'heights':[8,12,16], 'last_top':['-315',str(F(factorial(12),512)),str(F(factorial(16),8192))]},
 {'name':'native_nontrivial_T','matrix':[['2','0','0'],['-9/4','1/2','-3/2'],['3','0','1']], 'heights':[2,3,3], 'last_top':['-4','-48','-192']},
 {'name':'opposite_nontrivial_T','matrix':[['-9/2','4','-6'],['1/4','0','0'],['-3/2','0','-1']], 'heights':[6,9,13], 'last_top':['-1440',str(-factorial(9)*8),str(-factorial(13)*32)]},
 {'name':'generic_nontrivial_T','matrix':[['4/3','-2/3','4/3'],['-1/6','4/3','-2/3'],['2/3','-4/3','5/3']], 'heights':[8,12,16], 'last_top':[str(F(-factorial(8)*2,3**4)),str(F(-factorial(12)*8,3**6)),str(F(-factorial(16)*8,3**8))]}]
EXPECTED={'schema':1,'logical_boundary':BOUNDARY,'primary_metadata':PRIMARY_METADATA,
 'universal_frame_contexts':{'native':'b=0 and a*d!=0','opposite':'d=0 and b*c!=0','generic':'b*d*(a*d-b*c)!=0'},
 'seed_degrees':[4,6,8],'linear_frames':FRAMES,
 'H_coefficients':[['1','1'],['-2','3','-1'],['2','-1','0','2']],
 'realization_k_controls':list(range(1,13)),'odd_multiplier_controls':[1,3,5,7,9,11,13],
 'prime_controls':[p for p in range(2,151) if all(p%j for j in range(2,p))]}
def need(v,m):
 if not v:raise ValueError(m)
def strict(a,b):
 return type(a) is type(b) and (set(a)==set(b) and all(strict(a[k],b[k]) for k in a) if isinstance(a,dict)
   else len(a)==len(b) and all(strict(x,y) for x,y in zip(a,b)) if isinstance(a,list) else a==b)
def h(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def validate_input(d):need(strict(d,EXPECTED),'typed immutable inputs/scope mismatch')
def source_check(d):
 need(type(d) is dict and set(d)=={'schema','logical_boundary','sources'},'ledger fields')
 need(type(d['schema']) is int and d['schema']==1,'integer schema')
 need(strict(d['logical_boundary'],BOUNDARY),'typed logical boundary')
 roles={'written_proof':(PN,'mathematical_proof','Theorem and Sections1-8'),
  'exact_checker':(CN,'exact_identity_and_source_support','Universal frame/top identities and finite degree controls'),
  'typed_inputs':(IN,'exact_control_data','Immutable typed frames, H polynomials and scope'),
  'quadratic_automorphism_primary':(PDF,'mathematical_quadratic_automorphism_generation_only',PRIMARY_METADATA['locator'])}
 need(type(d['sources']) is list and len(d['sources'])==4,'exact four source roles')
 out={}
 for row in d['sources']:
  need(type(row) is dict and set(row)=={'role','path','relation','locator','sha256'},'source record fields')
  need(all(type(v) is str for v in row.values()),'source string types');role=row['role']
  need(role in roles and role not in out,'source role');name,relation,locator=roles[role]
  need(row['path']==name and row['relation']==relation and row['locator']==locator,'source path/relation/locator')
  need(len(row['sha256'])==64 and all(c in '0123456789abcdef' for c in row['sha256']),'sha format')
  need((HERE/name).is_file() and h(HERE/name)==row['sha256'],'source hash '+role)
  if role=='written_proof':need(row['sha256']==PROOF,'reviewed proof hardpin')
  if role=='quadratic_automorphism_primary':need(row['sha256']==PRIMARY,'actual quadratic primary hardpin')
  out[role]=dict(row)
 need(set(out)==set(roles),'mandatory source closure');return out

# Sparse Fraction polynomials; tuple coordinates are x,y,t, or a,b,c,d,x,y,t.
def clean(p):return {e:F(c) for e,c in p.items() if c}
def con(v,n=3):return {(0,)*n:F(v)} if v else {}
def var(i,n=3):
 e=[0]*n;e[i]=1;return {tuple(e):F(1)}
def add(*ps):
 out={}
 for p in ps:
  for e,c in p.items():out[e]=out.get(e,F(0))+c
 return clean(out)
def sc(p,c):return clean({e:v*c for e,v in p.items()})
def mul(a,b):
 out={}
 for e,c in a.items():
  for f,d in b.items():
   g=tuple(x+y for x,y in zip(e,f));out[g]=out.get(g,F(0))+c*d
 return clean(out)
def pw(p,k,n=3):
 out=con(1,n)
 for _ in range(k):out=mul(out,p)
 return out
def product(ps,n=3):
 out=con(1,n)
 for p in ps:out=mul(out,p)
 return out
def deriv(p,i):
 out={}
 for e,c in p.items():
  if e[i]:
   z=list(e);z[i]-=1;out[tuple(z)]=c*e[i]
 return clean(out)
def jop(p,n=3):
 return add(mul(var(n-3,n),deriv(p,n-1)),sc(mul(var(n-1,n),deriv(p,n-2)),-1))
def jpower(p,k,n=3):
 for _ in range(k):p=jop(p,n)
 return p
def sub(p,images,n=3):
 out={}
 for e,c in p.items():out=add(out,sc(product([pw(images[i],k,n) for i,k in enumerate(e)],n),c))
 return out
def nf(p):
 out={}
 for (a,b,c),v in p.items():
  k=min(a,b)
  for r in range(k+1):
   e=(a-k,b-k,c+2*r);out[e]=out.get(e,F(0))+v*F((-1)**r*comb(k,r),2**k)
 return clean(out)
def degree(p):return max(map(sum,p),default=-1)
def top(p):
 D=degree(p);return {e:v for e,v in p.items() if sum(e)==D}
def height(p):
 p=nf(p);k=0
 while True:
  z=nf(jop(p))
  if not z:return k,p
  p=z;k+=1;need(k<=100,'unexpected height')
def images(H):
 x,y,t=[var(i) for i in range(3)]
 return [x,add(y,sc(mul(H,t),-1),sc(mul(x,pw(H,2)),F(-1,2))),add(t,mul(x,H))]
def beta(p,H):return nf(sub(p,images(H)))
def univariate(coeffs):return clean({(j,0,0):F(c) for j,c in enumerate(coeffs)})
def seeds():
 x,y,t=[var(i) for i in range(3)];u=add(con(1),x)
 f4=add(mul(u,y),sc(product([u,add(con(2),x),pw(t,2)]),F(-1,2)))
 f6=add(sc(mul(y,t),2),sc(t,2),mul(x,t),sc(mul(univariate([4,6,4,1]),pw(t,3)),-1))
 f8=add(sc(mul(y,t),8),sc(t,8),sc(mul(x,t),4),sc(mul(univariate([20,41,44,26,8,1]),pw(t,3)),-1))
 need(not nf(add(f4,sc(mul(pw(u,3),y),-1),sc(mul(u,add(con(2),x)),F(1,2)))),'f4 regular identity')
 need(not nf(add(mul(x,f6),sc(mul(t,add(pw(u,2),sc(mul(pw(u,4),pw(t,2)),-1))),-1))),'f6 regular identity')
 need(not nf(add(mul(x,f8),sc(mul(t,add(sc(pw(u,2),4),sc(product([pw(u,4),pw(add(con(2),x),2),pw(t,2)]),-1))),-1))),'f8 regular identity')
 return {4:nf(f4),6:nf(f6),8:nf(f8)}

def universal_controls():
 a,b,c,d,x,y,t=[var(i,7) for i in range(7)];Delta=add(mul(a,d),sc(mul(b,c),-1))
 L=add(mul(pw(a,2,7),x),sc(mul(pw(b,2,7),y),-2),sc(product([a,b,t],7),2))
 M=add(sc(mul(pw(c,2,7),x),F(-1,2)),mul(pw(d,2,7),y),sc(product([c,d,t],7),-1))
 N=add(product([a,c,x],7),sc(product([b,d,y],7),-2),mul(add(mul(a,d),mul(b,c)),t))
 ids=[]
 def identity(name,p):need(not p,'universal identity '+name);ids.append(name)
 identity('PGL2 quadratic-frame preservation with denominator cleared',add(pw(N,2,7),sc(mul(L,M),2),sc(mul(pw(Delta,2,7),add(pw(t,2,7),sc(mul(x,y),2))),-1)))
 def zero_parameter(p,i):return {e:v for e,v in p.items() if e[i]==0}
 for D in [4,6,8]:
  leading=mul(pw(L,3,7),M) if D==4 else sc(mul(pw(L,D-3,7),pw(N,3,7)),-1)
  for name,index,nu in [('native',1,2 if D==4 else 3),('opposite',3,6 if D==4 else 2*D-3),('generic',None,2*D)]:
   p=zero_parameter(leading,index) if index is not None else leading
   if D==4:
    expected=sc(product([pw(a,6,7),pw(d,2,7),pw(x,4,7)],7),-1) if name=='native' else sc(product([pw(b,6,7),pw(c,2,7),pw(x,4,7)],7),-360) if name=='opposite' else sc(product([pw(b,6,7),pw(d,2,7),pw(x,4,7)],7),-20160)
   else:
    expected=sc(product([pw(a,2*D-3,7),pw(d,3,7),pw(x,D,7)],7),-6) if name=='native' else sc(product([pw(b,2*D-3,7),pw(c,3,7),pw(x,D,7)],7),-factorial(nu)) if name=='opposite' else sc(product([pw(b,2*D-3,7),pw(d,3,7),pw(x,D,7)],7),-factorial(nu))
   last=jpower(p,nu,7);identity(str(D)+' '+name+' universal last top derivative',add(last,sc(expected,-1)))
   identity(str(D)+' '+name+' universal next derivative zero',jop(last,7))
 return ids

def exact_controls(inp):
 fs=seeds();x,y,t=[var(i) for i in range(3)];q=add(pw(t,2),sc(mul(x,y),2),con(-1));rows=[];coordinate=0
 for frame in inp['linear_frames']:
  ims=[add(*[sc(var(j),F(v)) for j,v in enumerate(row)]) for row in frame['matrix']]
  need(not nf(sub(q,ims)),'literal orthogonal frame')
  for i,D in enumerate(inp['seed_degrees']):
   g=nf(sub(fs[D],ims));nu,last=height(g);need(nu==frame['heights'][i] and degree(g)==D and degree(last)==D,'actual frame heights/c degree')
   need(top(last)=={(D,0,0):F(frame['last_top'][i])},'actual frame leading coefficient')
   for coeffs in inp['H_coefficients']:
    H=univariate(coeffs);k=len(coeffs)-1;actual=beta(g,H);expected=D+k*nu
    need(degree(actual)==expected,'first nonmonomial shear degree')
    lead=F(frame['last_top'][i])*F(coeffs[-1])**nu/factorial(nu)
    need(top(actual)=={(expected,0,0):lead},'unique pure-x leading coefficient')
    rows.append({'seed_degree':D,'frame':frame['name'],'height':nu,'H':coeffs,'degree':expected,'leading_coefficient':str(lead)})
   if frame['name']=='native':
    H=x;left=beta(g,H);right={}
    for j in range(nu+1):right=add(right,sc(mul(pw(H,j),jpower(g,j)),F(1,factorial(j))))
    need(left==nf(right),'independent direct-coordinate/exponential formula');coordinate+=1
 # Polynomial H additivity, constant-T absorption, and nonreduced inverse cancellation.
 absorption=0
 for coeffs in inp['H_coefficients']:
  H=univariate(coeffs);rest=add(H,con(-F(coeffs[0])))
  for v in [x,y,t]:
   need(beta(beta(v,rest),con(F(coeffs[0])))==beta(v,H),'constant T absorption')
   need(beta(beta(v,H),sc(H,-1))==v,'inverse shear word cancellation');absorption+=2
 for D,f in fs.items():
  g=beta(f,x);sw=nf(sub(g,[y,x,t]));nu,last=height(sw);need(nu==2*degree(g),'later non-T frame height')
  out=beta(sw,x);need(degree(out)==3*degree(g),'later reduced odd multiplier')
  need(beta(beta(beta(f,x),con(1)),sc(x,-1))==beta(f,con(1)),'interior T word has no forced growth')
  need(degree(beta(f,con(1)))==D,'T preserves degree');absorption+=2
 # An unrelated quartic can have odd transformed degree: actual/top heights must match.
 unrelated=add(mul(pw(x,3),y),pw(t,3));need(degree(unrelated)==4 and degree(beta(unrelated,pw(x,2)))==9,'no parity theorem for arbitrary quartic')
 arithmetic=0
 for k in inp['realization_k_controls']:
  need((4+2*k)%2==0 and 4+2*k>=6,'quartic forward/converse arithmetic')
  need((6+3*k)%3==0 and 6+3*k>=9,'six forward/converse arithmetic');arithmetic+=2
  for nu in [3,13,16]:
   for s in inp['odd_multiplier_controls']:
    D=s*(8+nu*k);need(D>=11 and (s==1 or D%s==0),'eight degree bound')
    need(s==1 or 2*((s-1)//2)+1==s,'odd multiplier converse');arithmetic+=1
 primes=0
 for p in inp['prime_controls']:
  actual=p==8 or any(p>=8+nu and (p-8)%nu==0 for nu in [3,13,16])
  want=p>=11 and (p%3==2 or p%13==8);need(actual==want,'exact prime-class finite control');primes+=1
 need(not any((p-8)>0 and (p-8)%nu==0 for p in [13,19] for nu in [3,13,16]),'13/19 excluded')
 need(jop(t)==x and not {e:int(v)%2 for e,v in jop(pw(t,2)).items() if int(v)%2},'actual characteristic2 product-height obstruction')
 return {'seed_polynomial_identities':3,'full_frame_height_cases':18,'nonmonomial_first_shears':rows,
  'coordinate_exponential_comparisons':coordinate,'T_absorption_and_inverse_word_controls':absorption,
  'later_full_reduced_word_controls':3,'degree_bound_and_converse_arithmetic_controls':arithmetic,
  'prime_class_controls':primes,'mathematical_scope_countercontrols':5}

def reject(fn,name):
 try:fn()
 except (ValueError,TypeError,KeyError):return
 raise ValueError('negative accepted '+name)
def negatives(inp,ledger):
 n=0
 for row in ledger['sources']:
  z=copy.deepcopy(ledger);z['sources']=[v for v in z['sources'] if v['role']!=row['role']];reject(lambda:source_check(z),'omission');n+=1
  z=copy.deepcopy(ledger);next(v for v in z['sources'] if v['role']==row['role'])['sha256']='0'*64;reject(lambda:source_check(z),'badpin');n+=1
 for key in BOUNDARY:
  z=copy.deepcopy(inp);z['logical_boundary'][key]=not z['logical_boundary'][key];reject(lambda:validate_input(z),'false scope '+key);n+=1
 z=copy.deepcopy(inp);z['schema']=True;reject(lambda:validate_input(z),'bool schema');n+=1
 z=copy.deepcopy(ledger);z['logical_boundary']['full_automorphism_orbit_degrees']=1;reject(lambda:source_check(z),'int bool scope');n+=1
 z=copy.deepcopy(ledger);z['sources'][0]['path']='../'+PN;reject(lambda:source_check(z),'path forgery');n+=1
 z=copy.deepcopy(inp);z['primary_metadata']['role']='classical_context_only';reject(lambda:validate_input(z),'generation source improperly demoted');n+=1
 return n
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--repo-root');ap.add_argument('--output');args=ap.parse_args()
 inp=load(HERE/IN);validate_input(inp);ledger=load(HERE/LN);src=source_check(ledger)
 ids=universal_controls();controls=exact_controls(inp);count=negatives(inp,ledger)
 result={'schema':1,'status':'PASS_HOPF_SEED_ORBIT_EXACT_DEGREE_AND_SOURCE_SUPPORT','logical_boundary':BOUNDARY,
  'sources':src,'source_ledger_sha256':h(HERE/LN),'universal_expansion_identities':ids,'finite_exact_controls':controls,
  'source_scope_negative_controls':count,'primary_metadata':PRIMARY_METADATA,
  'assurance':{'formal_verification':False,'independent_replay':False,'external_review':False}}
 text=json.dumps(result,indent=2,sort_keys=True)+'\n'
 if args.output:Path(args.output).write_text(text,encoding='utf-8',newline='\n')
 else:print(text,end='')
if __name__=='__main__':main()
