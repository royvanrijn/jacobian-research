#!/usr/bin/env python3
"""Independent degree9 base-field minor and full symbolic source-shear obstructions."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from math import isqrt
import json,hashlib,time
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3,evaluate
from verify_raw13_target_field import madd,ms,mm,eye,inv_adj
B=Path(__file__).resolve().parent;O=B/'positive-source-quadratic';ROOT=B.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def con(c,n):return {(0,)*n:Q(c)} if c else {}
def total(*ps):
 z={}
 for p in ps:z=add(z,p)
 return z

def pw(p,k,n):
 z=con(1,n)
 for _ in range(k):z=mul(z,p)
 return z

def hom(p,d):return {e:c for e,c in p.items() if sum(e[:3])==d}
def main():
 start=time.monotonic();raw=B/'raw13-specialization/input.json';F=[variable(i,3) for i in range(3)]
 for atom in json.loads(raw.read_text())['atoms']:F[atom['row']]=add(F[atom['row']],{tuple(atom['exponent']):Q(atom['coefficient'])})
 assert determinant3([[derivative(f,j) for j in range(3)] for f in F],3)==con(1,3)
 Q8=pw(F[0],2,3);S8=total(pw(F[1],3,3),scale(mul(F[0],pw(F[2],2,3)),27),scale(mul(F[1],F[2]),-9));T8=add(pw(F[1],2,3),scale(mul(F[0],S8),-4));known=[con(1,3)]+F+[Q8,S8,T8];assert all(max(map(sum,f))<=8 for f in known)
 # New complete degree9 base-field intersection:165 old-source functions do not exhaust V9;
 # all220 columns are independently reconstructed and a213 minor is checked.
 base=O/'base-field9';cp9=base/'checkpoint.json';C9=json.loads(cp9.read_text());ip9=base/'input.json';D9=json.loads(ip9.read_text());assert C9['status']=='complete' and C9['source_sha256']==D9['source_sha256']==sha(B/'raw13_target_field9.py') and C9['input_sha256']==sha(ip9) and C9['source_input_sha256']==D9['source_input_sha256']==sha(raw)
 p=D9['modulus'];assert p==1000000007 and all(p%i for i in range(2,isqrt(p)+1));basis=sorted([(i,j,k) for i in range(10) for j in range(10) for k in range(10) if i+j+k<=9],key=lambda e:(sum(e),-e[0],-e[1]));assert len(basis)==220 and list(map(list,basis))==D9['source_basis'];positions={e:i for i,e in enumerate(basis)};rows=[];edets=[]
 for ss,tt,uu in D9['target_points']:
  assert ss%p;I=eye();W=[[0,0,2*ss*ss*uu%p],[1,0,-ss*tt%p],[0,1,1]];W2=mm(W,W,p);W3=mm(W2,W,p);assert madd(madd(W3,ms(W2,-1,p),p),madd(ms(W,ss*tt,p),ms(I,-2*ss*ss*uu,p),p),p)==[[0]*3 for _ in range(3)]
  E=madd(madd(ms(W2,3,p),ms(W,-2,p),p),ms(I,ss*tt,p),p);Ei,ed=inv_adj(E,p);edets.append(ed);is_=pow(ss,p-2,p);i2=pow(2,p-2,p);A=ms(Ei,ss,p);BB=ms(madd(W,ms(E,-1,p),p),i2*is_%p,p);E2=mm(E,E,p);E3=mm(E2,E,p);CC=ms(madd(madd(ms(E2,5,p),ms(mm(W,E,p),-3,p),p),ms(E3,-2,p),p),i2*is_*is_%p,p);powers=[]
  for val in [A,BB,CC]:
   stack=[I]
   for i in range(9):stack.append(mm(stack[-1],val,p))
   powers.append(stack)
  values=[mm(mm(powers[0][i],powers[1][j],p),powers[2][k],p) for i,j,k in basis]
  for ki,z in enumerate(known):
   v=[[0]*3 for _ in range(3)]
   for e,c in z.items():v=madd(v,ms(values[positions[e]],int(c),p),p)
   assert v==ms(I,v[0][0],p)
   if 1<=ki<=3:assert v==ms(I,[ss,tt,uu][ki-1],p)
  rows.extend([[z[j][0] for z in values] for j in [1,2]])
 assert len(D9['target_points'])==107 and len(rows)==214 and len(D9['selected_rows'])==len(D9['minor_columns'])==len(set(D9['minor_columns']))==213
 selected=[]
 for rec in D9['selected_rows']:
  assert rec['component'] in [1,2] and rec['row_index']==2*rec['point']+rec['component']-1;selected.append([rows[rec['row_index']][j] for j in D9['minor_columns']])
 M=[v[:] for v in selected];det=1
 for i in range(213):
  j=next((j for j in range(i,213) if M[i][j]),None);assert j is not None
  if j!=i:
   for row in M:row[i],row[j]=row[j],row[i]
   det=-det
  pivot=M[i][i];det=det*pivot%p;iv=pow(pivot,p-2,p)
  for j in range(i+1,213):
   fac=M[i][j]*iv%p
   if fac:
    for k in range(i,213):M[k][j]=(M[k][j]-fac*M[k][i])%p
 assert det
 # All known7 functions are independent:four independent constant/linear jets plus3 distinct degree8 heads.
 assert [{e:c for e,c in f.items() if sum(e)<=1} for f in F]==[variable(i,3) for i in range(3)]
 heads=[hom(f,8) for f in known[-3:]];assert [list(h) for h in heads]==[[(6,0,2)],[(3,4,1)],[(4,2,2)]]
 # Independent full3-parameter source composition and polynomial inverse.
 cp=O/'checkpoint.json';C=json.loads(cp.read_text());ip=O/'input.json';D=json.loads(ip.read_text());assert C['status']=='complete' and C['source_sha256']==D['source_sha256']==sha(B/'positive_source_quadratic.py') and C['input_sha256']==sha(ip) and D['raw_input_sha256']==C['raw_input_sha256']==sha(raw)
 n=6;a,b,c,n1,n2,n3=[variable(i,n) for i in range(n)];q=total(mul(n1,pw(b,2,n)),mul(n2,mul(b,c)),mul(n3,pw(c,2,n)));E=[add(a,q),b,c];Ei=[add(a,q,Q(-1)),b,c];nu=[n1,n2,n3];G=[compose(f,E,n) for f in F]
 assert q==decode(D['q'],n) and E==[decode(v,n) for v in D['source_map']] and Ei==[decode(v,n) for v in D['source_inverse']] and G==[decode(v,n) for v in D['composed_map']]
 assert [compose(f,Ei+nu,n) for f in E]==[a,b,c] and [compose(f,E+nu,n) for f in Ei]==[a,b,c] and determinant3([[derivative(f,j) for j in range(3)] for f in E],n)==con(1,n)
 assert [hom(f,d) for f,d in zip(G,[7,9,10])]==[decode(v,n) for v in D['composed_top_forms']]
 assert [hom(f,d) for f,d in zip(G,[7,9,10])]==[scale(mul(pw(q,3,n),c),-1),scale(mul(pw(q,3,n),mul(pw(b,2,n),c)),12),scale(mul(pw(q,3,n),mul(pw(b,3,n),c)),8)]
 weights=[1,-1,-2,3,4,5];assert all(sum(a*b for a,b in zip(e,weights))==w for f,w in zip(G,[1,-1,-2]) for e in f)
 lifted=[compose(f,[a,b,c],n) for f in known[-3:]];obstructions=[hom(mul(q,derivative(f,0)),9) for f in lifted];d2=derivative(derivative(compose(F[2],[a,b,c],n),0),0);d3=derivative(derivative(derivative(compose(F[1],[a,b,c],n),0),0),0);obstructions+=[hom(scale(mul(pw(q,2,n),d2),Q(1,2)),9),hom(scale(mul(pw(q,3,n),d3),Q(1,6)),9)]
 expected=[scale(mul(mul(pw(a,5,n),pw(c,2,n)),q),6),scale(mul(mul(pw(a,2,n),mul(pw(b,4,n),c)),q),-108),scale(mul(mul(pw(a,3,n),mul(pw(b,2,n),pw(c,2,n))),q),144),scale(mul(mul(a,mul(pw(b,3,n),c)),pw(q,2,n)),24),scale(mul(mul(pw(b,2,n),c),pw(q,3,n)),12)]
 assert obstructions==expected==[decode(v,n) for v in D['degree9_obstruction_heads']] and [sorted({e[0] for e in f}) for f in obstructions]==[[5],[2],[3],[1],[0]]
 assert known[-3:]==[decode(v,3) for v in D['base_nonlinear_functions']]
 # Finite bookkeeping controls for the written induction:all possible Taylor orders,not sampled coefficients.
 cases=[]
 for d in [3,4,5]:
  for order in range(3*d+1):
   choices=[(m,es,order-sum(es)) for m in range(4) for es in product(range(d,6),repeat=m) if sum(es)<=order]
   assert all(m<=order//d for m,es,i in choices)
   if order<2*d:assert all(m<=1 for m,es,i in choices)
   if order<3*d:assert all(m<=2 for m,es,i in choices)
   if order==3*d:assert [(es,i) for m,es,i in choices if m==3]==[((d,d,d),0)]
  cases.append({'leading_order':d,'last_order_used':3*d})
 assert [v['order'] for v in D['leading_parameter_strata']]==[3,4,5]
 # The zero quadratic is the known base map;retain its different low-degree kernels.
 assert [compose(f,[variable(j,3) for j in range(3)]+[{}, {}, {}],3) for f in G]==F
 # Degree10 endpoint is a control,not a global improvement.
 Ec=[add(variable(0,3),pw(variable(1,3),2,3)),variable(1,3),variable(2,3)];control=[compose(f,Ec,3) for f in F];degrees=[max(map(sum,f)) for f in control];counts=[sum(sum(e)>=2 for e in f) for f in control];assert degrees==[7,9,10] and counts==[8,14,15] and determinant3([[derivative(f,j) for j in range(3)] for f in control],3)==con(1,3)
 points=[[-1,1,-8],[-2,1,4]];assert [[evaluate(f,pt) for f in control] for pt in points]==[[Q(0),Q(1),Q(0)]]*2
 # All-parameter collision transport uses the polynomial inverse;check E(E^-1 p)=p in parameter polynomials.
 pv=[scale(add(add(variable(0,3),scale(variable(1,3),-8)),scale(variable(2,3),64)),-1),con(1,3),con(-8,3)];qv=[add(con(-1,3),total(variable(0,3),scale(variable(1,3),4),scale(variable(2,3),16)),Q(-1)),con(1,3),con(4,3)]
 for pt,old in [(pv,[0,1,-8]),(qv,[-1,1,4])]:assert [compose(f,pt+[variable(i,3) for i in range(3)],3) for f in E]==[con(x,3) for x in old]
 ledger=json.loads((ROOT/'MATH_STATUS.json').read_text());entries={e['id']:e for e in ledger['entries']};pins={}
 for k in ['KCUB42','KCUB43']:
  e=entries[k];assert e['state']=='proved';pins[k]={'scope_fingerprint':hashlib.sha256((e['state']+'\0'+e['title']+'\0'+e['scope']).encode()).hexdigest()[:16],'proof_sha256':sha(ROOT/e['canonical_source'])}
 receipt={'status':'pass','claim':'KCUB44','checker_sha256':sha(Path(__file__)),'producer_sha256':C['source_sha256'],'degree9_producer_sha256':C9['source_sha256'],'arithmetic_sha256':sha(B/'verify.py'),'companion_helper_sha256':sha(B/'verify_raw13_target_field.py'),'input_sha256':sha(ip),'checkpoint_sha256':sha(cp),'degree9_input_sha256':sha(ip9),'degree9_checkpoint_sha256':sha(cp9),'discovery_sha256':sha(O/'discovery.json'),'degree9_base_source_columns':220,'degree9_base_field_kernel_dimension':7,'degree9_base_field_kernel_basis':['1','H1','H2','H3','H1^2','S','T'],'degree9_independent_minor_size':213,'exact_rational_specialization_points':107,'prime':p,'prime_checked':True,'minor_determinant_mod_prime':det,'minor_sha256':hashlib.sha256(json.dumps(selected,separators=(',',':')).encode()).hexdigest(),'all_specialization_E_determinants_nonzero':bool(all(edets)),'all_source_quadratic_parameters_retained':True,'exact_source_inverse_and_conjugacy_weights_checked':True,'five_Taylor_obstruction_heads_checked':True,'obstruction_a_exponents':[5,2,3,1,0],'formal_leading_order_cases':cases,'universal_formal_induction_and_weighted_initial_space_proof_written_not_formalized':True,'every_nonzero_quadratic_low_source8_field_intersection_dimension':2,'intersection_basis':['1','H1 composed E_q'],'source_degree7_intersection_also_dimension2':True,'target_degree_cap':None,'target_denominator_cap':None,'no_Keller_rational_target_repair_through_source_degree8':True,'extends_to_lower_affine_terms_and_affine_frames_after_the_source_shear':True,'does_not_classify_nonlinear_source_shears_in_arbitrary_intermediate_frames':True,'zero_quadratic_base_kernel_dimensions_degree7_8_9':[4,7,7],'transformed_source_degree9_not_classified':True,'endpoint_degree_interval_for_nonzero_source_quadratic':[9,10],'matching_degree10_control':{'source_shear':'a->a+b^2','row_degrees':degrees,'nonlinear_counts':counts,'full_determinant_one':True,'points':points,'image':[0,1,0],'global_improvement':False},'all_parameter_collision_transport_checked':True,'pilot_parameter_triples_are_discovery_only':True,'imports':pins,'goal_complete':False,'global_positive_bounds_and_unrestricted_profiles_unchanged':True,'trust_boundary':'Exact full source family,inverse,weight conjugacy,five Taylor heads,all finite induction index cases and one213-square rational rank minor are independently checked. Coefficientwise target-field descent,formal induction and finite weighted-filtration dimension argument are written,not formal or externally reviewed. The universal parameter theorem is not inferred from fixed-coefficient pilots. The original base field degree9 lemma is not transferred to the transformed field at degree9;only source cap8 is closed. KCUB42 supplies the generic cubic field presentation;KCUB43 preserves the zero-shear/global affine baseline boundaries.','seconds':time.monotonic()-start}
 proof=B/'POSITIVE_SOURCE_QUADRATIC_FIELD_OBSTRUCTION.md'
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (O/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS all positive-axis quadratics:base-field9 rank213/220 det',det,';five exact Taylor heads,all weighted branches;nonzero low-field8 dimension2,no Keller repair<=8;endpoint[9,10];global goal open',flush=True)
if __name__=='__main__':main()
