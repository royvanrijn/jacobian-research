#!/usr/bin/env python3
"""Independent rational reconstruction for weights(1,-1,-r),r>=4,degree<=7."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json,hashlib
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3
HERE=Path(__file__).resolve().parent;OUT=HERE/'exceptional-weight1'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def total(*ps):
 o={}
 for p in ps:o=add(o,p)
 return o

def power(p,k,n):
 o={(0,)*n:Q(1)}
 for _ in range(k):o=mul(o,p)
 return o

def const(c,n):return {(0,)*n:Q(c)} if c else {}
def bracket(f,g):return add(mul(derivative(f,0),derivative(g,1)),mul(derivative(f,1),derivative(g,0)),Q(-1))
def dictionary(f,g,h,r):return total(mul(g,bracket(f,h)),mul(r,mul(h,bracket(g,f))),mul(f,bracket(g,h)))
def coefficient(p,i,k):
 return {e[:i]+e[i+1:]:c for e,c in p.items() if e[i]==k}
def degree(p,i):return max((e[i] for e in p),default=-1)
def main():
 cp=OUT/'checkpoint.json';check=json.loads(cp.read_text());assert check['status']=='complete' and check['source_sha256']==sha(HERE/'exceptional_weight1.py');assert [x['r'] for x in check['models']]==[4,5,6,7]
 atom_total=0;generators_total=0;model_hashes={};shape_checks=[]
 for meta in check['models']:
  r=meta['r'];mp=OUT/meta['input'];d=json.loads(mp.read_text());assert meta['input_sha256']==sha(mp) and d['source_sha256']==sha(HERE/'exceptional_weight1.py');es=[[e for e in product(range(8),repeat=3) if 2<=sum(e)<=7 and e[0]-e[1]-r*e[2]==t] for t in [1,-1,-r]];assert d['basis']==[[list(e) for e in row] for row in es];names=['c'+str(j+1) for j in range(sum(map(len,es)))];assert d['parameter_names']==names;n=3+len(names);x=[variable(j,n) for j in range(3)];F=list(x);pos=0
  for j,row in enumerate(es):
   for e in row:F[j]=add(F[j],mul({e+(0,)*len(names):Q(1)},variable(3+pos,n)));pos+=1
  assert F==[decode(p,n) for p in d['source_outputs']];ni=2+len(names);o=const(1,ni);u,w=variable(0,ni),variable(1,ni);sub=[o,u,w]+[variable(2+j,ni) for j in range(len(names))];f,g,h=[compose(z,sub,ni) for z in F];assert [f,g,h]==[decode(p,ni) for p in d['quotient_outputs']];D=dictionary(f,g,h,const(r,ni));assert D==decode(d['quotient_jacobian'],ni);lift=[mul(x[0],x[1]),mul(power(x[0],r,n),x[2])]+[variable(3+j,n) for j in range(len(names))];actual=determinant3([[derivative(p,j) for j in range(3)] for p in F],n);assert compose(D,lift,n)==actual
  eq=[]
  for mon in sorted({e[:2] for e in add(D,o,Q(-1))},reverse=True):eq.append({e[2:]:c for e,c in add(D,o,Q(-1)).items() if e[:2]==mon})
  assert eq==[decode(p,len(names)) for p in d['native_generators']]
  assert degree(f,1)==int(r<6) and degree(g,1)<=1 and degree(h,1)<=2
  assert degree(coefficient(g,1,1),0)<=int(r<6) and degree(coefficient(g,1,0),0)<=4
  assert degree(coefficient(h,1,1),0)<=3 and degree(coefficient(h,1,2),0)<=0
  E=coefficient(h,1,0);assert all(r<=e[0]<=r+int(r<6) for e in E)
  shape_checks.append({'r':r,'slots':list(map(len,es)),'complete_dictionary':True});atom_total+=pos;generators_total+=len(eq);model_hashes[mp.name]=sha(mp);print('MODEL PASS',r,pos,len(eq),flush=True)
 assert (atom_total,generators_total)==(52,72)
 up=OUT/'structural-identities.json';data=json.loads(up.read_text());assert data['status']=='complete-primary-identities' and data['source_sha256']==sha(HERE/'exceptional_weight1_structure.py')
 d=data['univariate'];names=['r','mu','p','z','a','b','b2','b3','b4','C','d1','d2','d3']+['e'+str(j) for j in range(1,8)];assert d['parameters']==names;n=2+len(names);u,w=variable(0,n),variable(1,n);P={z:variable(2+j,n) for j,z in enumerate(names)};one=const(1,n);phi=total(one,mul(P['mu'],u),mul(P['p'],power(u,2,n)),mul(P['z'],power(u,3,n)));A=add(P['a'],mul(P['b'],u));B=total(u,*(mul(P['b'+str(j)],power(u,j,n)) for j in range(2,5)));C=P['C'];D=total(one,*(mul(P['d'+str(j)],power(u,j,n)) for j in range(1,4)));E=total(*(mul(P['e'+str(j)],power(u,j,n)) for j in range(1,8)));g=add(mul(A,w),B);h=total(mul(C,power(w,2,n)),mul(D,w),E);assert [phi,g,h]==[decode(p,n) for p in d['forms']];J=dictionary(phi,g,h,P['r']);assert J==decode(d['dictionary'],n)
 dp=lambda p:derivative(p,0)
 wanted=[add(mul(dp(phi),add(mul(B,D),mul(P['r'],mul(A,E)),Q(-1))),mul(phi,add(mul(dp(B),D),mul(A,dp(E)),Q(-1)))),add(mul(dp(phi),add(mul(add(one,P['r'],Q(-1)),mul(A,D)),scale(mul(B,C),2))),mul(phi,total(mul(dp(A),D),scale(mul(dp(B),C),2),scale(mul(A,dp(D)),-1)))),add(mul(add(scale(one,2),P['r'],Q(-1)),mul(dp(phi),mul(A,C))),scale(mul(phi,mul(dp(A),C)),2))]
 assert [coefficient(J,1,j) for j in range(3)]==[{e[:1]+e[2:]:c for e,c in p.items()} for p in wanted]==[decode(p,n-1) for p in d['w_coefficients']]
 # Constant first quotient specializes to the elementary plane identities.
 S=[variable(j,n) for j in range(n)]
 for z in ['mu','p','z']:S[2+names.index(z)]={}
 pure=compose(J,S,n);expected=total(scale(mul(dp(A),mul(C,power(w,2,n))),2),mul(total(mul(dp(A),D),scale(mul(dp(B),C),2),scale(mul(A,dp(D)),-1)),w),mul(dp(B),D),scale(mul(A,dp(E)),-1));assert pure==expected
 # Graph coefficient equations, retaining symbolic r and every degree-six boundary tail.
 d=data['graph'];names=['r','A','B','alpha','beta','T']+['h'+str(j) for j in range(7)]+['k'+str(j) for j in range(4)];assert d['parameters']==names;n=2+len(names);u,v=variable(0,n),variable(1,n);P={z:variable(2+j,n) for j,z in enumerate(names)};one=const(1,n);G0=add(mul(P['A'],u),P['B']);G1=add(mul(P['alpha'],u),P['beta']);H0=total(*(mul(P['h'+str(j)],power(u,j,n)) for j in range(7)));H1=total(*(mul(P['k'+str(j)],power(u,j,n)) for j in range(4)));G=add(G0,mul(v,G1));H=total(H0,mul(v,H1),mul(P['T'],power(v,2,n)));graph=add(mul(add(G,mul(v,derivative(G,1))),dp(H)),mul(dp(G),add(mul(P['r'],H),mul(v,derivative(H,1)))),Q(-1));assert [G0,G1,H0,H1,G,H]==[decode(p,n) for p in d['forms']] and graph==decode(d['dictionary'],n)
 cofs=[add(mul(G0,dp(H0)),mul(P['r'],mul(P['A'],H0)),Q(-1)),total(scale(mul(G1,dp(H0)),2),mul(G0,dp(H1)),scale(mul(P['r'],mul(P['alpha'],H0)),-1),scale(mul(add(P['r'],one),mul(P['A'],H1)),-1)),total(scale(mul(G1,dp(H1)),2),scale(mul(add(P['r'],one),mul(P['alpha'],H1)),-1),scale(mul(add(P['r'],scale(one,2)),mul(P['A'],P['T'])),-1)),scale(mul(add(P['r'],scale(one,2)),mul(P['alpha'],P['T'])),-1)]
 assert [coefficient(graph,1,j) for j in range(4)]==[{e[:1]+e[2:]:c for e,c in p.items()} for p in cofs]==[decode(p,n-1) for p in d['v_coefficients']]
 # Affine-G vector field and every monomial eigenvalue needed through total degree six.
 n=5;u,v,a,b,r=[variable(j,n) for j in range(n)];one=const(1,n);q=add(mul(a,u),mul(b,add(v,one,Q(-1))));Xu=add(q,mul(v,b));Xv=scale(mul(a,v),-1)
 act=lambda p:add(mul(Xu,derivative(p,0)),mul(Xv,derivative(p,1)))
 assert act(q)==mul(a,q) and act(v)==scale(mul(a,v),-1)
 eigen=0
 for i in range(7):
  for j in range(7-i):
   mon=mul(power(q,i,n),power(v,j,n));assert add(act(mon),mul(a,mul(r,mon)),Q(-1))==mul(a,mul(add(const(i-j,n),r,Q(-1)),mon));eigen+=1
 # Complete normalized r5 residual model and denominator-free triangular elimination.
 d=data['r5'];assert d['parameters']==['alpha','P','R','T'];n=5;u,alpha,P,R,T=[variable(j,n) for j in range(n)];one=const(1,n);A=scale(add(one,alpha),Q(-1,2));mu=scale(add(scale(one,3),alpha,Q(-1)),Q(1,2));N=total(one,mul(mu,u),mul(P,power(u,2,n)),mul(R,power(u,3,n)));G0=add(mul(A,u),one,Q(-1));G1=add(one,mul(alpha,u));H1=scale(power(G1,3,n),-1);temp=mul(N,power(G1,3,n));H0=add({e:c for e,c in temp.items() if e[0]<5},mul(T,power(u,5,n)));assert [N,G0,G1,H0,H1]==[decode(p,n) for p in d['forms']];E0=total(mul(G0,dp(H0)),scale(mul(A,H0),-5),scale(one,-1));E1=total(scale(mul(G1,dp(H0)),2),mul(G0,dp(H1)),scale(mul(alpha,H0),-5),scale(mul(A,H1),-6));assert E0==decode(d['E0'],n) and E1==decode(d['E1'],n)
 # g's source linear condition and constant boundary equation force exactly these A,mu.
 assert total(A,alpha,mu,scale(one,-1))=={} and total(mu,scale(alpha,3),scale(A,5),one)=={}
 sub=[variable(j,4) for j in range(4)]
 assert d['triangular_elimination_indices']==[1,2,4]
 for j,(name,k,lead) in enumerate(zip(['P','R','T'],[1,2,4],[-2,-3,-5]),1):
  eq=compose(coefficient(E0,0,k),sub,4);assert derivative(eq,j)==const(lead,4)
  sol=decode(d['substitutions'][name],1);lift={(e[0],0,0,0):c for e,c in sol.items()};sub[j]=lift;assert compose(coefficient(E0,0,k),sub,4)=={}
 residual=[compose(coefficient(E1,0,k),[variable(0,1)]+[decode(d['substitutions'][z],1) for z in ['P','R','T']],1) for k in [1,2]];assert d['residual_indices']==[1,2] and residual==[decode(p,1) for p in d['residuals']];multipliers=[decode(p,1) for p in d['unit_multipliers']];assert total(*(mul(a,b) for a,b in zip(residual,multipliers)))==const(1,1);unit_terms=sum(map(len,multipliers));assert unit_terms==d['multiplier_monomials']==5
 # Constant-first quadratic tails: coefficients imply 2t^2=0 or5t^3=0.
 d=data['constant_first_quadratic_exclusions'];n=4;u,t,p,q=[variable(j,n) for j in range(n)];B4=add(u,mul(p,power(u,2,n)));B5=add(B4,mul(q,power(u,3,n)));eq4=[coefficient(add(add(B4,mul(t,power(B4,2,n))),u,Q(-1)),0,k) for k in [2,3]];eq5=[coefficient(add(add(B5,mul(t,power(B5,2,n))),u,Q(-1)),0,k) for k in [2,3,4]];assert eq4==[decode(p,3) for p in d['r4_equations']] and eq5==[decode(p,3) for p in d['r5_equations']];t,p,q=[variable(j,3) for j in range(3)];assert add(scale(mul(t,eq4[0]),2),eq4[1],Q(-1))==decode(d['r4_target'],3)==scale(power(t,2,3),2);assert total(eq5[2],scale(mul(t,eq5[1]),-2),scale(mul(t,mul(add(p,scale(t,5),Q(-1)),eq5[0])),-1))==decode(d['r5_target'],3)==scale(power(t,3,3),5)
 # Complete two-chart normal forms and generic two-sided polynomial inverses.
 controls=[]
 for r in [4,5,6,7]:
  for chart in ['A=0','A free']:
   names=['A','e0','e1'];n=6;x=[variable(j,n) for j in range(3)];A,e0,e1=[variable(j,n) for j in range(3,6)]
   if chart=='A=0':A={}
   elif r==4:e1={}
   else:e0={};e1={}
   if r>=6:e1={}
   shear=add(mul(e0,power(x[1],r,n)),mul(e1,mul(x[0],power(x[1],r+1,n))));third=add(x[2],shear);second=add(x[1],mul(A,mul(power(x[0],r-1,n),third)));F=[x[0],second,third];bi=add(x[1],mul(A,mul(power(x[0],r-1,n),x[2])),Q(-1));ci=add(x[2],add(mul(e0,power(bi,r,n)),mul(e1,mul(x[0],power(bi,r+1,n)))),Q(-1));inv=[x[0],bi,ci];extra=[variable(j,n) for j in range(3,6)];assert [compose(p,inv+extra,n) for p in F]==x and [compose(p,F+extra,n) for p in inv]==x;assert determinant3([[derivative(p,j) for j in range(3)] for p in F],n)==const(1,n);deg=max(sum(e[:3]) for p in F for e in p);assert deg<=7;assert all(e[0]-e[1]-r*e[2]==target for p,target in zip(F,[1,-1,-r]) for e in p);controls.append({'r':r,'chart':chart,'degree':deg,'generic_both_inverse_compositions':True});print('CHART PASS',r,chart,'degree',deg,flush=True)
 proof=HERE/'EXCEPTIONAL_UNIT_WEIGHT_CLASSIFICATION.md';receipt={'status':'pass','claim':'KCUB34','checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'model_producer_sha256':sha(HERE/'exceptional_weight1.py'),'structure_producer_sha256':sha(HERE/'exceptional_weight1_structure.py'),'checkpoint_sha256':sha(cp),'structural_identities_sha256':sha(up),'model_hashes':model_hashes,'full_models':shape_checks,'complete_source_atoms':atom_total,'native_coefficient_generators':generators_total,'univariate_w_coefficient_identities':3,'graph_v_coefficient_identities':4,'affine_vectorfield_eigenvalue_checks':eigen,'r5_denominator_free_eliminations':3,'r5_univariate_unit_multiplier_monomials':unit_terms,'quadratic_tail_exclusion_identities':2,'generic_classified_inverse_charts':controls,'all_r_at_least8_covered_by_written_degree_bound':True,'all_r_at_least4_classification':'composition of two elementary shears with exact source-degree product restrictions','no_external_plane_theorem_imported':True,'prior_r2_r3_not_recomputed':True,'positive_bounds_and_full_profile_counts_unchanged':True,'goal_complete':False,'trust_boundary':'Exact finite source and polynomial reconstruction plus written all-weight module, Wronskian/field-unit, rational-derivative and degree arguments; not formal or externally reviewed.'}
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (OUT/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS complete(1,-1,-r),r>=4,degree<=7 classification;52 atoms,72 generators,5-term unit,8 generic inverse charts; global goal open',flush=True)
if __name__=='__main__':main()
