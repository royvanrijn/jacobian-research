#!/usr/bin/env python3
"""Independent exact classification of all135 weighted cubic triangular carriers.
No SymPy/Singular import. Polynomial input powers are guarded; Laurent parameters
are used only as substituted values and multiplicative coefficient units.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import combinations_with_replacement,permutations
from collections import Counter
import hashlib,json,time
from verify import add,mul,variable,derivative,decode
from verify_power_single_compress import scale
from weighted_verify_utils import expression,substitution,coefficient_buckets,equal_generators,normal,power
from verify_weighted_unmarked7 import laurent
HERE=Path(__file__).resolve().parent;ROOT=HERE/'carrier-reparam-full'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def const(c,n):return {(0,)*n:Q(c)} if c else {}
def sub(poly,values,n):
 assert all(all(k>=0 for k in e) for e in poly)
 return substitution(poly,values,n)
def det(M,n):
 answer={};k=len(M)
 for perm in permutations(range(k)):
  term=const((-1)**sum(perm[i]>perm[j] for i in range(k) for j in range(i+1,k)),n)
  for i,j in enumerate(perm):term=mul(term,M[i][j])
  answer=add(answer,term)
 return answer

def main():
 started=time.monotonic();source=HERE/'cubic12-integer-witness.json';F=[decode(p,12) for p in load(source)['coefficients']];base=load(HERE/'integer-validation12.json');assert base['status']=='pass' and base['witness_sha256']==sha(source) and base['determinant_one'] and base['identity_derivative'];ip=ROOT/'input.json';d=load(ip);assert d['source_sha256']==sha(HERE/'carrier_reparam_full.py') and d['input_sha256']==sha(source);Z=d['variables'];assert Z==['z'+str(i) for i in range(1,136)];xn=['x'+str(i) for i in range(1,13)];names=xn+Z;N=len(names);atoms=d['atoms'];weights=[1,-1,-2,0,-3,2,-1,-2,-2,-1,-1,0];original={(i,e):c for i in range(3,12) for e,c in F[i].items() if sum(e)>=2};expected=set()
 for i in range(3,12):
  for degree in [2,3]:
   for js in combinations_with_replacement(range(i),degree):
    if sum(weights[j] for j in js)==weights[i]:expected.add((i,tuple(js.count(a) for a in range(12))))
 assert len(expected)==len(atoms)==135 and {(a['row'],tuple(a['exponent'])) for a in atoms}==expected
 new=[variable(i,N) for i in range(12)]
 for j,a in enumerate(atoms):
  assert a['parameter']==Z[j];e=tuple(a['exponent']);c=Q(a['base_coefficient']);assert c==original.get((a['row'],e),Q(1));new[a['row']]=add(new[a['row']],{e+tuple(int(k==j) for k in range(135)):c})
 E=[expression(text,names) for text in d['source_substitution']];assert len(E)==12 and E[:3]==new[:3]
 for i in range(3,12):assert derivative(E[i],i)==const(1,N) and not any(any(e[i+1:12]) for e in E[i]);assert not any(sum(e[:12])<2 for e in add(E[i],variable(i,N),Q(-1)))
 G=[sub(f,E,N) for f in F];assert G[3:]==new[3:] and [expression(text,names) for text in d['carrier_outputs']]==new[3:];high=[];low=[]
 for i,f in enumerate(G):
  for e,p in coefficient_buckets(f,12).items():
   if sum(e)>3:high.append(p)
   elif sum(e)>=2:low.append((i,e,p))
 gen=[{tuple(e):Q(c) for e,c in terms} for terms in d['high_equations']];equal_generators(high,gen);assert {(i,e):p for i,e,p in low}=={(a['row'],tuple(a['exponent'])):expression(a['coefficient'],Z) for a in d['low_coefficients']};assert len(gen)==540 and len(low)==184;print('RECONSTRUCTED135/540/184',flush=True)
 cp=ROOT/'classification-certificate.json';c=load(cp);assert c['source_sha256']==sha(HERE/'carrier_full_lift.py') and c['input_sha256']==sha(ip) and c['singular_source_sha256']==sha(ROOT/'classification-lift.sing') and c['variables']==Z
 active={1,2,3,4,5,6,7,8,9,10,11,12,13,17,18,19,22,26,27,28,30,34,35,38,39,42,45,79,80,83,96,117,118};texts=['z'+str(j) for j in range(1,136) if j not in active]+['z3-1','z1-z2','z7-z9','z8-z9','z10-z12','z2*z9-1','z5*z12-z1','z4-z9*z12','z11-z4-3*z22','z5*z42+z22','z2*z17-2*z13','z2*z18-z5*z22','z5*z38-3*z5*z96-z17','z12*z30-z96','z12*z34-z117','z12*z35-z118','z2*z80-z27','z9*z28-z118','z9*z26-z79','z12*z79-z12*z117+3*z28','z9*z19+3*z34','z6*z9-z9-z19','z5*z39+z96','z5*z45-z80','z83-z118'];targets=[expression(text,Z) for text in texts];assert len(targets)==127 and c['target_texts']==texts and targets==[{tuple(e):Q(v) for e,v in terms} for terms in c['targets']];results=[{} for _ in targets];terms=0;degree=0;seen=set()
 for entry in c['entries']:
  t,g=entry['target'],entry['generator'];assert (t,g) not in seen;seen.add((t,g));m={tuple(e):Q(v) for e,v in entry['multiplier']};assert all(len(e)==135 and all(k>=0 for k in e) for e in m);results[t]=add(results[t],mul(m,gen[g]));terms+=len(m);degree=max(degree,max(sum(e) for e in m))
 assert results==targets and terms==410;print('REPLAYED127 MEMBERSHIPS',terms,flush=True)
 fp=ROOT/'family.json';family=load(fp);assert family['source_sha256']==sha(HERE/'carrier_full_family.py') and family['input_sha256']==sha(ip);pn=['p','q','s','t','u','v','w','z'];assert family['variables']==pn and family['nonzero']==['p','q'];explicit={1:'1/p',2:'1/p',3:'1',4:'p*q',5:'1/(p*q)',6:'1-3*u/(p**2*q)',7:'p',8:'p',9:'p',10:'q',11:'p*q+3*t',12:'q',13:'s/(2*p)',17:'s',18:'t/q',19:'-3*u/(p*q)',22:'t',26:'(u-3*v/(p*q))/p',27:'w/p',28:'v/p',30:'z/q',34:'u/q',35:'v/q',38:'p*q*s+3*z',39:'-p*q*z',42:'-p*q*t',45:'p*q*w',79:'u-3*v/(p*q)',80:'w',83:'v',96:'z',117:'u',118:'v'};values=[laurent(explicit.get(j,'0'),pn) for j in range(1,136)];assert values==[laurent(family['parameterization'][name],pn) for name in Z];assert not any(sub(f,values,8) for f in targets+gen)
 reduced={(i,e):q for i,e,p in low if (q:=sub(p,values,8))};assert len(reduced)==64 and reduced=={(a['row'],tuple(a['exponent'])):laurent(a['coefficient'],pn) for a in family['low_coefficients']}
 # Normalize by p,q without ever passing Laurent input to polynomial-only sub.
 normalized={}
 for key,f in reduced.items():
  out={}
  for e,c0 in f.items():out[e[2:]]=out.get(e[2:],Q(0))+c0
  normalized[key]={e:c0 for e,c0 in out.items() if c0}
 nn=pn[2:];assert all(all(k>=0 for k in e) and sum(e)<=1 for f in normalized.values() for e in f)
 D=[laurent(text,pn) for text in ['1','1','1','p','1','1/(p*q)','p*q','1','1/p','1/p','1/(p*q)','1/(p*q)']];Dinv=[laurent(text,pn) for text in ['1','1','1','1/p','1','p*q','1/(p*q)','1','p','p','p*q','p*q']];invariants=[laurent(text,pn) for text in ['s/p','t/(p*q)','u/(p**2*q)','v/(p**3*q**2)','w/p**2','z/(p**2*q)']]
 for (i,e),f in normalized.items():
  scale0=Dinv[i]
  for j,k in enumerate(e):scale0=mul(scale0,power(D[j],k,8))
  assert reduced[i,e]==mul(scale0,sub(f,invariants,8))
 factors=Counter(normal(f) for f in normalized.values());expected_texts={'1':31,'u-1':1,'t-u/3+v+1/3':1,'t+w':1,'t+v':1,'u':4,'v-w':1,'s+3/2':1,'s+3*z':2,'z':5,'s':2,'u-1/3':1,'u-3*v':2,'w':3,'v':4,'t':3,'t+1/3':1};assert factors==Counter({normal(expression(text,nn)):count for text,count in expected_texts.items()});at0={key:f.get((0,)*6,Q(0)) for key,f in normalized.items()};at0={key:value for key,value in at0.items() if value};assert at0=={(i,e):value for i,f in enumerate(F) for e,value in f.items() if sum(e)>=2} and len(at0)==36
 # Full nonlinear-output rank: coefficient minors within disjoint weight blocks.
 def mon(text):return next(iter(expression(text,xn)))
 def block(rows,words):return [[reduced.get((i,mon(text)),{}) for text in words] for i in rows]
 rank=[]
 for rows,words in [([0],['x1**2*x2']),([4],['x2**3']),([5],['x1**2']),([3,11],['x1**2*x3','x1*x10']),([1,6,9,10],['x4*x10','x1*x3','x1*x8','x1*x9'])]:
  determinant=det(block(rows,words),8);assert len(determinant)==1 and all(e[j]==0 for e in determinant for j in range(2,8));rank.append(dict(rows=[i+1 for i in rows],monomials=words,determinant=[[list(e),str(c0)] for e,c0 in determinant.items()]))
 quad_words=['x2**2','x2*x7','x7**2'];quad=[reduced[7,mon(text)] for text in quad_words];comb=add(add(quad[0],mul(laurent('3/(p*q)',pn),quad[1])),mul(laurent('9/(p**2*q**2)',pn),quad[2]));assert comb==const(1,8)
 for word,co in zip(quad_words,quad):assert det(block([2,7,8],['x1*x4*x5','x1*x5',word]),8)==scale(mul(laurent('p**2',pn),co),Q(-16))
 # The standard homogeneous companion uses n+rank(C)+1 variables.
 # Five disjoint-weight cubic output pivots are constant in normalized variables.
 def nblock(rows,words):return [[normalized.get((i,mon(text)),{}) for text in words] for i in rows]
 cubic_rows=[0,1,2,3,4];cubic_words=['x1**2*x2','x1*x4*x8','x1*x4*x5','x1**2*x3','x2**3'];assert det(nblock(cubic_rows,cubic_words),6)==const(-288,6)
 assert det(nblock([2,8],['x1*x4*x5','x1*x2*x3']),6)==expression('8*s',nn)
 assert det(nblock([3,11],['x1**2*x3','x1*x2*x4']),6)==expression('u',nn)
 assert det(nblock([3,11],['x1**2*x3','x1*x4*x7']),6)==expression('v',nn)
 assert det(nblock([1,9],['x1*x4*x8','x1*x3*x4']),6)==expression('12*w',nn)
 assert det(nblock([1,10],['x1*x4*x8','x1*x3*x4']),6)==expression('12*z',nn)
 # Select only t as a polynomial parameter; no Laurent polynomial is input.
 only_t=[{},variable(1,6),{},{},{},{}];r5={(i,e):sub(f,only_t,6) for (i,e),f in normalized.items() if sum(e)==3};r5={key:value for key,value in r5.items() if value}
 assert all(i in cubic_rows+[11] for i,e in r5)
 for e in {e for i,e in r5 if i in [3,11]}:assert r5.get((11,e),{})==mul(variable(1,6),r5.get((3,e),{}))
 # For p=q1 on the rank5 locus, G_t=S(t) K S(t)^-1,
 # S=I+t E_(10,7)+t E_(12,4). Verify the full map, not only cubic parts.
 width=20;tp=variable(15,width);source_values=[variable(i,width) for i in range(12)];source_values[9]=add(source_values[9],mul(tp,source_values[6]),Q(-1));source_values[11]=add(source_values[11],mul(tp,source_values[3]),Q(-1));conjugate=[sub(f,source_values,width) for f in F];conjugate[9]=add(conjugate[9],mul(tp,conjugate[6]));conjugate[11]=add(conjugate[11],mul(tp,conjugate[3]));expected_map=[variable(i,width) for i in range(12)];tvals=[{},tp,{},{},{},{}]
 for (i,e),f in normalized.items():expected_map[i]=add(expected_map[i],mul({e+(0,)*8:Q(1)},sub(f,tvals,width)))
 assert conjugate==expected_map
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),laurent_utility_sha256=sha(HERE/'verify_weighted_unmarked7.py'),base_witness_sha256=sha(source),base_validation_sha256=sha(HERE/'integer-validation12.json'),input_sha256=sha(ip),certificate_sha256=sha(cp),family_sha256=sha(fp),parameters=135,high_equations=540,low_slots_before_reduction=184,exact_memberships=127,multiplier_monomials=terms,max_multiplier_degree=degree,all_memberships_multiplied=True,complete_weighted_degree23_triangular_atoms=True,all_source_substitutions_triangular=True,entire_composition_reconstructed=True,complete_eight_parameter_family_sufficiency_checked=True,nonzero_family_parameters=['p','q'],free_family_parameters=['s','t','u','v','w','z'],generic_nonlinear_total=64,diagonal_normalization_of_entire_family_checked=True,normalized_affine_coefficient_groups=expected_texts,minimum_nonlinear_total=36,minimum_equality='s=t=u=v=w=z=0,p*q!=0',minimum_family_is_diagonal_KCUB3=True,next_possible_support=39,rank_blocks=rank,remaining_rank_block_rows=[3,8,9],remaining_rank_block_three_minors_checked=True,remaining_rank_block_coefficient_bezout_checked=True,nonlinear_output_rank=12,no_affine_linear_output_component=True,polynomial_substitution_domain_guarded=True,minimum_cubic_output_rank=5,minimum_standard_homogeneous_companion_dimension=18,cubic_rank5_locus='s=u=v=w=z=0,p*q!=0,t arbitrary',rank5_full_linear_conjugacy_checked=True,written_elimination_and_counting_formalized=False,dimension=12,ordinary_degree=3,seconds=round(time.monotonic()-started,3),scope='All135 weighted degree2/3 carrier monomials triangular in the fixed original order, for G=KCUB3 composed with C_old^-1 C_new fixing the three base variables. Complete eight-parameter family, sharp36 displayed terms, equality only the diagonal KCUB3 family and no37/38 terms. Full nonlinear-output rank12 excludes an affine-linear output component in any affine frame, hence affine-linear component deletion within this class. The cubic-output rank is at least5, with rank5 exactly s=u=v=w=z=0; the whole rank5 locus is linearly conjugate to KCUB3, so the standard n+rank(C)+1 homogeneous companion has minimum18 and supplies no new tensor there up to linear equivalence. No nonlinear graph deletion, other order/base/frame/carrier realization, other map or global simplicity/dimension minimum follows. Base determinant/collision is inherited from pinned KCUB3, not re-proved.')
 (ROOT/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
