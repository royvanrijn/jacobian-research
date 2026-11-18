#!/usr/bin/env python3
"""Independently verify a group-algebra Fox-ideal Bezout identity in char2."""
from pathlib import Path
from collections import Counter
from copy import deepcopy
import hashlib,json,time
HERE=Path(__file__).resolve().parent
PRESENTATION=HERE/'uniform_fox_ideal_inputs.json'
CERTIFICATE=HERE/'uniform_fox_a8_char2_certificate.json'
RECEIPT=HERE/'uniform_fox_ideal_verification.json'
ONE=tuple(range(8))
def mul(a,b):return tuple(a[b[i]]for i in range(8))
def red(w):
 out=[]
 for a in w:
  if out and out[-1]==-a:out.pop()
  else:out.append(a)
 return out
def inv(w):return [-a for a in reversed(w)]
def odd(p):return sum(p[i]>p[j]for i in range(8)for j in range(i+1,8))%2
def reconstruct():
 d=json.loads(PRESENTATION.read_text())
 for q in d['source_inputs']:
  assert hashlib.sha256((HERE/q['path']).read_bytes()).hexdigest()==q['sha256'],q['path']
 case=json.loads((HERE/d['source_inputs'][0]['path']).read_text())['case']
 signed=json.loads((HERE/d['source_inputs'][1]['path']).read_text())['case']
 assert case['literal_braid_word']==signed['literal_braid_word']==d['literal_braid_word']
 assert case['pd']==signed['pd'] and len(d['literal_braid_word'])==17
 fam=d['family'];assert (fam['tail_generator'],fam['tail_initial'],fam['tail_step'])==(2,3,6)
 assert fam['prefix']+[2]*3==d['literal_braid_word']
 top=[tuple(case['arc_permutations'][i])for i in case['seed_arcs']]
 assert [list(p)for p in top]==d['ordered_A8_top']
 assert all(mul(p,p)==ONE and not odd(p)for p in top)
 f=[[1],[2],[3],[4]];labels=list(case['seed_arcs']);trace=[]
 for e,row in zip(d['literal_braid_word'],signed['canonical_signed_relations']):
  o,a,b,sign=row;j=abs(e)-1;old=labels[:];u,v=f[j:j+2]
  assert(sign==-1)==(e>0)
  if e>0:
   assert labels[j:j+2]==[a,o];labels[j:j+2]=[o,b];f[j:j+2]=[v,red(inv(v)+u+v)]
  else:
   assert labels[j:j+2]==[o,a];labels[j:j+2]=[b,o];f[j:j+2]=[red(u+v+inv(u)),u]
  trace.append({'braid_letter':e,'signed_PD_relation':row,'incoming_arc_labels':old,'outgoing_arc_labels':labels[:]})
 assert labels==case['seed_arcs'] and trace==d['signed_PD_arc_tracking']
 assert f==d['bottom_free_words'] and red([x for w in f for x in w])==[1,2,3,4]
 rels=[red(w+[-i-1])for i,w in enumerate(f)]
 assert rels==d['all_four_closure_relators'] and d['deficiency_one_relator_indices']==[0,1,2]
 allD=[]
 for rel in rels:
  D=[Counter()for _ in top];cur=ONE
  for e in rel:
   i=abs(e)-1
   if e>0:D[i][cur]+=1;cur=mul(cur,top[i])
   else:cur=mul(cur,top[i]);D[i][cur]-=1
  assert cur==ONE
  identity=Counter()
  for p,F in zip(top,D):
   for g,c in F.items():identity[mul(g,p)]+=c;identity[g]-=c
  assert not any(identity.values()),'integer Fox fundamental identity'
  allD.extend([{g:c for g,c in F.items()if c}for F in D])
 G={ONE};queue=[ONE]
 for g in queue:
  for h in top:
   z=mul(g,h)
   if z not in G:G.add(z);queue.append(z)
 assert len(G)==20160
 # A single finite return is the inductive step, not a sample of family values.
 tail=top[:];tail_trace=[]
 for _ in range(6):
  u,v=tail[1:3];tail[1:3]=[v,mul(mul(v,u),v)]
  tail_trace.append([list(g)for g in tail])
 assert tail==top
 # Every terminal sigma2 changes only coordinates1,2. Coordinate0 (and R0)
 # remains exactly fixed for any number6n of such steps by induction.
 assert 0 not in [1,2]
 return d,allD[:12],G,tail_trace

def convolution(certificate,derivatives,group,right_side=False):
 assert certificate['field_prime']==2 and certificate['group_order']==20160
 indices=certificate['derivative_pairs'];coeffs=certificate['left_coefficients_permutation_tuples']
 assert len(indices)==len(coeffs)==2
 result=set();counts=[]
 for pair,tuples in zip(indices,coeffs):
  j,i=pair;assert 0<=j<3 and 0<=i<4
  left=[tuple(g)for g in tuples]
  assert len(left)==len(set(left)) and all(sorted(g)==list(range(8)) and not odd(g) and g in group for g in left)
  right=[g for g,c in derivatives[4*j+i].items()if c%2]
  counts.append((len(left),len(right)))
  for a in left:
   for b in right:
    z=mul(b,a)if right_side else mul(a,b)
    if z in result:result.remove(z)
    else:result.add(z)
 assert result=={ONE},('Not the unit',len(result))
 return counts


def padd(a,b):
 out=dict(a)
 for e,c in b.items():out[e]=out.get(e,0)+c
 return {e:c for e,c in out.items()if c}
def pmul(a,b):
 out={}
 for e,c in a.items():
  for f,k in b.items():out[e+f]=out.get(e+f,0)+c*k
 return {e:c for e,c in out.items()if c}
def scalar_identity(d):
 H=[{}for _ in range(4)];e=0
 for x in d['all_four_closure_relators'][0]:
  i=abs(x)-1
  if x>0:H[i]=padd(H[i],{e:1});e+=1
  else:e-=1;H[i]=padd(H[i],{e:-1})
 assert e==0
 assert H==[{int(e):c for e,c in q.items()}for q in d['R0_Fox_entries_ascending_laurent']]
 u=d['unit_identity'];assert u['zero_based_relator_index']==0
 assert u['zero_based_derivative_generator_indices']==[1,2]and u['normalizing_t_powers']==[4,3]
 f={e+4:c for e,c in H[1].items()};g={e+3:c for e,c in H[2].items()}
 expected=[{e:c for e,c in enumerate(q)if c}for q in u['normalized_polynomials_ascending']]
 assert [f,g]==expected
 A,B=[{e:c for e,c in enumerate(q)if c}for q in u['integer_Bezout_polynomials_ascending']]
 assert padd(pmul(A,f),pmul(B,g))=={0:1}
 return {'integer_Laurent_identity':'t^4*A*H01+t^3*B*H02=1','unit_pair_zero_based':[1,2],'integer_coefficients_only':True,'all_commutative_scalar_models_excluded':True}
def free_action(braid):
 f=[[1],[2],[3],[4]]
 for e in braid:
  j=abs(e)-1;u,v=f[j:j+2]
  f[j:j+2]=[v,red(inv(v)+u+v)]if e>0 else[red(u+v+inv(u)),u]
 return f

def run_verification():
 start=time.monotonic();d,D,G,tail_trace=reconstruct()
 c=json.loads(CERTIFICATE.read_text());assert c['presentation_sha256']==hashlib.sha256(PRESENTATION.read_bytes()).hexdigest()
 assert c['derivative_pairs']==[[0,0],[0,1]]
 u=c['uniform_family'];assert(u['tail_generator'],u['tail_crossings_per_n'],u['fixed_relator_index'])==(2,6,0)
 counts=convolution(c,D,G);assert [a for a,b in counts]==c['coefficient_supports']==[9792,369]
 scalar=scalar_identity(d)
 # Finite controls of the written local substitution induction.
 prefix=d['family']['prefix'];fixed=free_action(prefix)[0]
 for n in [0,1,2]:assert free_action(prefix+[2]*(3+6*n))[0]==fixed
 assert red(fixed+[-1])==d['all_four_closure_relators'][0]
 rejected=[]
 for name,mut in [('remove_one_left_coefficient',lambda x:x['left_coefficients_permutation_tuples'][0].pop()),('omit_second_term',lambda x:x['left_coefficients_permutation_tuples'].__setitem__(1,[])),('wrong_derivative',lambda x:x['derivative_pairs'][1].__setitem__(1,2))]:
  bad=deepcopy(c);mut(bad)
  try:convolution(bad,D,G)
  except(AssertionError,ValueError,IndexError):rejected.append(name)
  else:raise AssertionError('Accepted corrupt certificate '+name)
 try:convolution(c,D,G,right_side=True)
 except AssertionError:rejected.append('wrong_multiplication_side')
 else:raise AssertionError('Wrong multiplication side accepted')
 for name,mut in [('scalar_Bezout_coefficient',lambda x:x['unit_identity']['integer_Bezout_polynomials_ascending'][0].__setitem__(0,17)),('scalar_wrong_unit_pair',lambda x:x['unit_identity'].__setitem__('zero_based_derivative_generator_indices',[0,1]))]:
  bad=deepcopy(d);mut(bad)
  try:scalar_identity(bad)
  except(AssertionError,ValueError,IndexError):rejected.append(name)
  else:raise AssertionError('Accepted corrupt scalar certificate '+name)
 files=[PRESENTATION.name,CERTIFICATE.name,Path(__file__).name,'UNIFORM_FOX_IDEAL_MODEL_BARRIERS.tex']+[q['path']for q in d['source_inputs']]
 edition={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()for name in files}
 return {'status':'PASS exact uniform scalar and A8 Fox-ideal unit certificates','edition_sha256':edition,'scalar_verification':scalar,'A8_verification':{'field_prime':2,'group_order':len(G),'identity':'X*D00+Y*D01=1','left_coefficient_supports':[a for a,b in counts],'derivative_supports_mod2':[b for a,b in counts],'all_four_relators_evaluate_to_identity':True,'all_integer_Fox_fundamental_identities_verified':True,'all_signed_PD_crossings_reconstructed':True,'permutation_tuples_are_literal_not_indices':True,'six_step_return_trace':tail_trace,'regular_rank_or_simple_module_classification_needed':False},'all_n_uniformity':'Written induction: terminal sigma2 changes only slots1,2, so R0 is fixed; one six-step permutation return iterates to a valid fixed quotient for all n>=0. Representative word controls are not the universal proof.','rejected_controls':rejected,'ordinary_rank_status':'UNRESOLVED','scope_limits':['No ordinary group rank upper improvement or rank-four lower proof.','The full noncommutative ideal in Z[G_n] is not shown to contain1.','Only this presentation and fixed A8 quotient kernel; other matrix models remain open.','A8 unit extends to 2-nilpotent rings and local rings with fixed-A8 characteristic-two residue action by the separate written lifting proof.'],'assurance':{'whole_proof_independent_replay':False,'formal_verification':False,'external_specialist_review':False},'elapsed_seconds':round(time.monotonic()-start,3)}
def main():
 out=run_verification();old=json.loads(RECEIPT.read_text())if RECEIPT.exists()else{}
 for key in ['internal_pdf_reviews','root_review']:
  if key in old:out[key]=old[key]
 RECEIPT.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
if __name__=='__main__':main()
