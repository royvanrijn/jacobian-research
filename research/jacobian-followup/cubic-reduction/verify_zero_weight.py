#!/usr/bin/env python3
"""Independent exact source/principal-part audit for the all-degree zero-weight proof."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,permutations
import json,hashlib,time
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3,evaluate
HERE=Path(__file__).resolve().parent;OUT=HERE/'zero-weight';ROOT=HERE.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def total(*ps):
 out={}
 for p in ps:out=add(out,p)
 return out

def con(c,n):return {(0,)*n:Q(c)} if c else {}
def power(p,k,n):
 out=con(1,n)
 for _ in range(k):out=mul(out,p)
 return out

def cof(p,i,k):return {e[:i]+e[i+1:]:c for e,c in p.items() if e[i]==k}
def kill(p,indices):return {e:c for e,c in p.items() if all(e[j]==0 for j in indices)}
def determinant(A,n):
 out={};size=len(A)
 for perm in permutations(range(size)):
  term=con((-1)**sum(perm[i]>perm[j] for i in range(size) for j in range(i+1,size)),n)
  for i in range(size):term=mul(term,A[i][perm[i]])
  out=add(out,term)
 return out

def main():
 started=time.monotonic();cp=OUT/'checkpoint.json';check=json.loads(cp.read_text());assert check['status']=='complete' and check['source_sha256']==sha(HERE/'zero_weight.py');pairs=[(1,1),(1,2),(1,3),(1,4),(1,5),(1,6),(2,3),(2,5),(3,4)];assert check['pairs']==list(map(list,pairs));hashes={};atoms=0;removed=0;reduced=0;faces=0;ngens=0;controls=[]
 for rec,pair in zip(check['models'],pairs):
  pp,qq=pair;ip=OUT/rec['input'];d=json.loads(ip.read_text());assert rec['pair']==list(pair) and d['pair']==list(pair) and rec['input_sha256']==sha(ip) and d['source_sha256']==sha(HERE/'zero_weight.py');basis=[]
  for i in range(3):
   for degree in range(2,8):
    for a in range(degree+1):
     for b in range(degree-a+1):
      e=(a,b,degree-a-b)
      if pp*e[1]-qq*e[2]==[0,pp,-qq][i]:basis.append((i,e))
  basis.sort();N=len(basis);assert d['parameter_names']==['c'+str(j) for j in range(N)] and d['basis']==[{'output':i,'exponent':list(e)} for i,e in basis];ns=3+N;source=[variable(j,ns) for j in range(3)];F=list(source);ni=2+N;x,t=variable(0,ni),variable(1,ni);forms=[x,con(1,ni),con(1,ni)];module=[];zero=[]
  for j,(i,e) in enumerate(basis):
   ey,ez=e[1]-int(i==1),e[2]-int(i==2);assert ey%qq==0 and ez%pp==0 and ey//qq==ez//pp;k=ey//qq;module.append({'output':i,'x_degree':e[0],'t_degree':k,'parameter':j});F[i]=add(F[i],mul({e+(0,)*N:Q(1)},variable(3+j,ns)));forms[i]=add(forms[i],mul({(e[0],k)+(0,)*N:Q(1)},variable(2+j,ni)))
   if k==0:zero.append(j)
  assert module==d['module_terms'] and zero==d['pure_x_parameters'] and len(zero)==18;assert F==[decode(z,ns) for z in d['source_outputs']] and forms==[decode(z,ni) for z in d['forms']];J=determinant3([[derivative(z,j) for j in range(3)] for z in F],ns);D={}
  for e,c in J.items():assert e[1]%qq==0 and e[2]%pp==0 and e[1]//qq==e[2]//pp;D[(e[0],e[1]//qq)+e[3:]]=c
  A,B,C=forms;formula=add(mul(derivative(A,0),total(mul(B,C),mul(t,total(scale(mul(C,derivative(B,1)),qq),scale(mul(B,derivative(C,1)),pp))))),mul(derivative(A,1),mul(t,total(scale(mul(C,derivative(B,0)),qq),scale(mul(B,derivative(C,0)),pp)))),Q(-1));assert D==formula==decode(d['quotient_jacobian'],ni);I=add(D,con(1,ni),Q(-1));eqs=[{e[2:]:c for e,c in I.items() if e[:2]==mon} for mon in sorted({e[:2] for e in I},reverse=True)];assert eqs==[decode(z,N) for z in d['native_generators']];ngens+=len(eqs)
  # t0 gives exactly the complete pure-x unit product.
  A0,B0,C0=[cof(z,1,0) for z in forms];assert cof(D,1,0)==mul(derivative(A0,0),mul(B0,C0));RR=[kill(z,[2+j for j in zero]) for z in forms];assert RR==[decode(z,ni) for z in d['boundary_reduced_forms']]
  slopes=sorted({Q(m['x_degree']-int(m['output']==0),m['t_degree']) for m in module if m['t_degree'] and m['x_degree']-int(m['output']==0)>0});assert len(slopes)==len(d['positive_principal_faces'])
  for face,slope in zip(d['positive_principal_faces'],slopes):
   r,s=slope.numerator,slope.denominator;omega=[qq*s,-r,0];v=[r,qq*s,pp*s];assert face['slope']==str(slope) and face['auxiliary_weights']==omega and face['ray_exponent']==v;kept=[m['parameter'] for m in module if m['t_degree'] and Q(m['x_degree']-int(m['output']==0),m['t_degree'])<=slope];top=[m['parameter'] for m in module if m['t_degree'] and Q(m['x_degree']-int(m['output']==0),m['t_degree'])==slope];assert kept==face['retained_parameters'] and top==face['top_parameters'] and top;excluded=[j for j in range(N) if j not in kept];FR=[kill(z,[3+j for j in excluded]) for z in F];initial=[]
   for i,f in enumerate(FR):
    assert all(sum(a*b for a,b in zip(e[:3],omega))<=omega[i] for e in f);initial.append({e:c for e,c in f.items() if sum(a*b for a,b in zip(e[:3],omega))==omega[i]})
   actual_top=determinant3([[derivative(z,j) for j in range(3)] for z in initial],ns);specialized=kill(J,[3+j for j in excluded]);assert actual_top=={e:c for e,c in specialized.items() if sum(a*b for a,b in zip(e[:3],omega))==0};nf=1+N;Z=variable(0,nf);fi=[]
   for i,f in enumerate(initial):
    out={}
    for e,c in f.items():
     ev=[e[j]-int(j==i) for j in range(3)];assert min(ev)>=0 and ev[0]%r==0;k=ev[0]//r;assert ev==[k*a for a in v];out[(k,)+e[3:]]=c
    assert max(e[0] for e in out)<=2;fi.append(out)
   rankone=mul(fi[0],mul(fi[1],fi[2]))
   for i in range(3):
    term=scale(mul(Z,derivative(fi[i],0)),v[i])
    for j in range(3):
     if j!=i:term=mul(term,fi[j])
    rankone=add(rankone,term)
   projected={}
   for e,c in actual_top.items():assert e[0]%r==0 and list(e[:3])==[e[0]//r*a for a in v];projected[(e[0]//r,)+e[3:]]=c
   assert projected==rankone;faces+=1
  # Once the universal leading-face proof removes every positive slope, verify the final univariate factor.
  positive=[m['parameter'] for m in module if m['t_degree'] and m['x_degree']-int(m['output']==0)>0];AA,BB,CC=[kill(z,[2+j for j in zero+positive]) for z in forms];fx=derivative(AA,0);assert not derivative(fx,0) and not derivative(BB,0) and not derivative(CC,0);tail=add(mul(BB,CC),mul(t,total(scale(mul(CC,derivative(BB,1)),qq),scale(mul(BB,derivative(CC,1)),pp))));assert kill(D,[2+j for j in zero+positive])==mul(fx,tail)
  # Exact final generic source family and inverse, with every allowed degree-seven coefficient free.
  n=3+7//(pp+qq);x,y,z=[variable(j,n) for j in range(3)];T=mul(power(y,qq,n),power(z,pp,n));shear=total(*(mul(variable(3+j-1,n),power(T,j,n)) for j in range(1,1+7//(pp+qq))));G=[add(x,shear),y,z];inverse=[add(x,shear,Q(-1)),y,z];params=[variable(j,n) for j in range(3,n)];assert [compose(a,inverse+params,n) for a in G]==[x,y,z] and [compose(a,G+params,n) for a in inverse]==[x,y,z] and determinant3([[derivative(a,j) for j in range(3)] for a in G],n)==con(1,n);assert all(pp*e[1]-qq*e[2]==w for a,w in zip(G,[0,pp,-qq]) for e in a);controls.append({'pair':list(pair),'allowed_shear_terms':7//(pp+qq),'both_generic_inverse_compositions':True});hashes[ip.name]=sha(ip);atoms+=N;removed+=18;reduced+=N-18;print('MODEL PASS',pp,qq,N,'atoms',len(slopes),'principal faces',flush=True)
 assert (atoms,removed,reduced,faces)==(240,162,78,12)
 # Complete degree triples for the rank-one principal maps, with arbitrary lower coefficients.
 leading=0
 for degrees in product(range(3),repeat=3):
  n=10;Z=variable(0,n);vv=[variable(j,n) for j in range(1,4)];coeffs=[[variable(4+2*i,n),variable(5+2*i,n)] for i in range(3)];fi=[total(con(1,n),*(mul(coeffs[i][j-1],power(Z,j,n)) for j in range(1,d+1))) for i,d in enumerate(degrees)];expr=mul(fi[0],mul(fi[1],fi[2]))
  for i in range(3):
   term=mul(Z,mul(vv[i],derivative(fi[i],0)))
   for j in range(3):
    if j!=i:term=mul(term,fi[j])
   expr=add(expr,term)
  target=con(1,9)
  for i,d in enumerate(degrees):target=add(target,scale(variable(i,9),d))
  for i,d in enumerate(degrees):
   if d:target=mul(target,variable(3+2*i+d-1,9))
  assert cof(expr,0,sum(degrees))==target;leading+=1
 # All degree<=3 pairs for the final B,C bracket, with symbolic original positive weights.
 tail_cases=0
 for db,dc in product(range(4),repeat=2):
  n=9;T,p,q=[variable(j,n) for j in range(3)];B=total(con(1,n),*(mul(variable(2+j,n),power(T,j,n)) for j in range(1,db+1)));C=total(con(1,n),*(mul(variable(5+j,n),power(T,j,n)) for j in range(1,dc+1)));expr=total(mul(B,C),mul(T,total(mul(q,mul(C,derivative(B,0))),mul(p,mul(B,derivative(C,0))))));target=total(con(1,8),scale(variable(1,8),db),scale(variable(0,8),dc));target=mul(target,con(1,8) if db==0 else variable(1+db,8));target=mul(target,con(1,8) if dc==0 else variable(4+dc,8));assert cof(expr,0,db+dc)==target;tail_cases+=1
 # Generic same-sign equal-weight block determinant and adjugate identities.
 n=7;a,b,c,d,e,f,g=[variable(j,n) for j in range(n)];delta=add(mul(c,g),mul(d,f),Q(-1));assert determinant3([[a,{},{}],[b,c,d],[e,f,g]],n)==mul(a,delta);assert add(mul(c,g),mul(d,f),Q(-1))==delta and add(mul(c,scale(d,-1)),mul(d,c))=={} and add(mul(f,g),mul(g,scale(f,-1)))=={} and add(mul(g,c),mul(f,d),Q(-1))==delta
 # Independent four-dimensional identity stabilization of the existing13-term control.
 n=4;x,a,b,c=[variable(j,n) for j in range(n)];one=con(1,n);v=add(one,scale(mul(a,b),2));v2=power(v,2,n);v3=power(v,3,n);F=[x,total(a,scale(mul(power(a,2,n),b),-3),scale(mul(power(a,3,n),c),-1)),add(scale(mul(v2,mul(a,c)),3),mul(total(scale(v2,9),scale(v,-6),scale(one,-2)),b)),add(mul(v3,c),mul(add(scale(v2,6),scale(v,2)),power(b,2,n)))];assert determinant([[derivative(z,j) for j in range(n)] for z in F],n)==one;assert [{e:co for e,co in z.items() if sum(e)<=1} for z in F]==[x,a,b,c];assert [sum(sum(e)>=2 for e in z) for z in F]==[0,2,5,6] and max(sum(e) for z in F for e in z)==7;row_degrees=[max(sum(e) for e in z) for z in F];assert row_degrees==[1,4,6,7];old_control=HERE/'raw13-specialization/input.json';raw=json.loads(old_control.read_text());assert {(i-1,e[1:]):co for i,z in enumerate(F) for e,co in z.items() if i>0 and sum(e)>=2}=={(z['row'],tuple(z['exponent'])):Q(z['coefficient']) for z in raw['atoms']};weights=[0,1,-1,-2];assert all(sum(u*v for u,v in zip(e,weights))==w for z,w in zip(F,weights) for e in z);points=[[Q(0),Q(0),Q(1),Q(-8)],[Q(0),Q(-1),Q(1),Q(4)]];assert [[evaluate(z,p) for z in F] for p in points]==[[Q(0),Q(0),Q(1),Q(0)]]*2
 parent=HERE/'weighted-312-412/validation.json';P=json.loads(parent.read_text());assert P['status']=='pass' and len(P['remaining_diagonal_weight_improvement_types'])==10;remaining=[z for z in P['remaining_diagonal_weight_improvement_types'] if z['mixed_positive_negative_magnitudes']==[1,1,1]];assert len(remaining)==1;ledger=json.loads((ROOT/'MATH_STATUS.json').read_text());entries={e['id']:e for e in ledger['entries']};pins={}
 for key in ['KCUB38','KCUB23','KCUB36']:
  entry=entries[key];assert entry['state']=='proved';pins[key]={'scope_fingerprint':hashlib.sha256((entry['state']+'\0'+entry['title']+'\0'+entry['scope']).encode()).hexdigest()[:16],'proof_sha256':sha(ROOT/entry['canonical_source'])};proof=HERE/'ZERO_WEIGHT_CLASSIFICATION.md';receipt={'status':'pass','claim':'KCUB39','checker_sha256':sha(Path(__file__)),'producer_sha256':sha(HERE/'zero_weight.py'),'arithmetic_sha256':sha(HERE/'verify.py'),'checkpoint_sha256':sha(cp),'model_hashes':hashes,'complete_degree7_models':9,'native_atoms':atoms,'native_coefficient_generators':ngens,'boundary_removed_pure_x_atoms':removed,'boundary_reduced_atoms':reduced,'principal_face_specializations':faces,'rank_one_leading_degree_cases':leading,'tail_leading_degree_cases':tail_cases,'generic_shear_inverse_controls':controls,'same_sign_matrix_unit_identities_checked':True,'all_degree_one_zero_theorem_has_written_filtration_proof':True,'minimum_collision_dimension_exactly_one_zero_form':4,'dimension_lower_bound_is_written_all_degree_one_zero_and_dimension2_unit_proof':True,'existing_four_dimensional_upper_control':{'dimension':4,'weights':weights,'degree':7,'nonlinear_counts':[0,2,5,6],'full_determinant_one':True,'row_degrees':row_degrees,'old_raw13_input_sha256':sha(old_control),'points':[[str(z) for z in p] for p in points],'common_image':['0','0','1','0'],'new_witness':False},'parent_frontier_receipt_sha256':sha(parent),'imports':pins,'remaining_diagonal_weight_raw12_type':remaining,'strict_flag_below764_has_no_nontrivial_linear_Gm_symmetry':True,'strict_corollary_uses_imported_KCUB23_and_known112_flag':True,'old111_timeouts_untouched':True,'no_new_solver_query':True,'goal_complete':False,'global_positive_bounds_and_unrestricted_profile_counts_unchanged':True,'trust_boundary':'Independent complete degree7 sources,principal Jacobian identities,leading coefficients and upper-control determinant/collision checks; universal all-degree normalized-support maximum,filtered determinant,field-unit and same-sign proofs are written,not formalized or externally reviewed. KCUB38/KCUB36 imported for the full frontier and known112 classification;KCUB23 imported only for its strict flag branch,not full111.','seconds':time.monotonic()-started}
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (OUT/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS all-degree one-zero three-variable closure:240 native atoms,12 principal faces;exactly-one-zero dimension minimum4 with checked stabilized control;degree7 weighted raw<=12 frontier10->1(full111);global goal open',flush=True)
if __name__=='__main__':main()
