#!/usr/bin/env python3
"""Independent rational replay of both complete twelve-slot frontier closures."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json,hashlib
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3
HERE=Path(__file__).resolve().parent;OUT=HERE/'weighted-124-125';ROOT=HERE.parents[1]
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
def deg(p,i):return max((e[i] for e in p),default=-1)
def bracket(f,g):return add(mul(derivative(f,0),derivative(g,1)),mul(derivative(f,1),derivative(g,0)),Q(-1))
def dictionary(f,g,h,r):return total(mul(f,bracket(g,h)),scale(mul(g,bracket(f,h)),2),mul(r,mul(h,bracket(g,f))))
def main():
 cp=OUT/'checkpoint.json';check=json.loads(cp.read_text());assert check['status']=='complete' and check['source_sha256']==sha(HERE/'weighted124125.py');model_hashes={};atoms=0;native=0
 for meta in check['models']:
  r=meta['r'];assert r in [4,5];ip=OUT/meta['input'];d=json.loads(ip.read_text());assert d['source_sha256']==sha(HERE/'weighted124125.py') and meta['input_sha256']==sha(ip);names=['mu','P','nu','A2','A3','L','B']+(['E2','E3'] if r==4 else ['E3','E4'])+['H1','H2','T'];assert d['parameter_names']==names;n=3+len(names);x=[variable(i,n) for i in range(3)];P={z:variable(3+j,n) for j,z in enumerate(names)};supports=[{(1,0):'mu',(2,0):'P',(0,1):'nu'},{(2,0):'A2',(3,0):'A3',(0,1):'L',(1,1):'B'},{(r-2,0):'E'+str(r-2),(r-1,0):'E'+str(r-1),(1,1):'H1',(2,1):'H2',(0,2):'T'}];basis=[];F=list(x);seen=[]
  for i in range(3):
   for degree in range(2,8):
    for a in range(degree+1):
     for b in range(degree-a+1):
      e=(a,b,degree-a-b)
      if e[0]-2*e[1]-r*e[2]!=[1,-2,-r][i]:continue
      name=supports[i][e[1:]];basis.append((i,e));seen.append(name);F[i]=add(F[i],mul({e+(0,)*len(names):Q(1)},P[name]))
  assert len(basis)==len(seen)==12 and set(seen)==set(names);basis.sort();assert d['basis']==[{'output':i,'exponent':list(e)} for i,e in basis];assert F==[decode(z,n) for z in d['source_outputs']];ni=2+len(names);u,w=variable(0,ni),variable(1,ni);sub=[con(1,ni),u,w]+[variable(2+j,ni) for j in range(len(names))];f,g,h=[compose(z,sub,ni) for z in F];assert [f,g,h]==[decode(z,ni) for z in d['quotient_outputs']];D=dictionary(f,g,h,con(r,ni));assert D==decode(d['quotient_jacobian'],ni);lift=[mul(power(x[0],2,n),x[1]),mul(power(x[0],r,n),x[2])]+[variable(3+j,n) for j in range(len(names))];assert compose(D,lift,n)==determinant3([[derivative(z,j) for j in range(3)] for z in F],n)
  eqs=[{e[2:]:c for e,c in add(D,con(1,ni),Q(-1)).items() if e[:2]==m} for m in sorted({e[:2] for e in add(D,con(1,ni),Q(-1))},reverse=True)];assert eqs==[decode(z,len(names)) for z in d['native_generators']];atoms+=12;native+=len(eqs);parent=HERE/'weight-resonance7/mixed-cases'/('p1q2r'+str(r)+'.json');old=json.loads(parent.read_text());assert d['parent_timeout_model_sha256']==sha(parent) and old['basis']==d['basis']
  # Full graph transformation at nu=-1, before any inferred boundary reductions.
  S=[variable(j,ni) for j in range(ni)];S[2+names.index('nu')]=con(-1,ni);ff,gg,hh=[compose(z,S,ni) for z in [f,g,h]];N=add(ff,w);GS=[u,add(N,w,Q(-1))]+[variable(2+j,ni) for j in range(len(names))];GG,HH=[compose(z,GS,ni) for z in [gg,hh]];assert compose(ff,GS,ni)==w;graph=add(mul(add(scale(GG,2),mul(w,derivative(GG,1))),derivative(HH,0)),mul(derivative(GG,0),add(scale(HH,r),mul(w,derivative(HH,1)))),Q(-1));assert graph==compose(compose(D,S,ni),GS,ni);assert deg(GG,1)<=1 and deg(HH,1)<=2 and deg(cof(GG,1,0),0)<=3 and deg(cof(HH,1,0),0)<=4 and deg(cof(GG,1,1),0)<=1 and deg(cof(HH,1,1),0)<=2 and deg(cof(HH,1,2),0)<=0;model_hashes[ip.name]=sha(ip);print('MODEL PASS',r,12,len(eqs),flush=True)
 assert (atoms,native)==(24,30)
 sp=OUT/'structural-identities.json';struct=json.loads(sp.read_text());assert struct['source_sha256']==sha(HERE/'weighted124125.py') and check['structural_identities_sha256']==sha(sp)
 # Universal q=2 coefficient identities and the weighted target-square shear.
 d=struct['univariate'];names=['r','mu','P','A0','A1','B2','B3','C','D1','D2','E2','E3','E4','rho'];assert d['parameters']==names;n=2+len(names);u,w=variable(0,n),variable(1,n);P={z:variable(2+j,n) for j,z in enumerate(names)};one=con(1,n);phi=total(one,mul(P['mu'],u),mul(P['P'],power(u,2,n)));A=add(P['A0'],mul(P['A1'],u));B=total(u,mul(P['B2'],power(u,2,n)),mul(P['B3'],power(u,3,n)));DD=total(one,mul(P['D1'],u),mul(P['D2'],power(u,2,n)));E=total(*(mul(P['E'+str(j)],power(u,j,n)) for j in [2,3,4]));C=P['C'];g=add(mul(A,w),B);h=total(mul(C,power(w,2,n)),mul(DD,w),E);assert [phi,A,B,DD,E,g,h]==[decode(z,n) for z in d['forms']];J=dictionary(phi,g,h,P['r']);assert J==decode(d['jacobian'],n);dp=lambda z:derivative(z,0);cofs=[add(mul(dp(phi),add(scale(mul(B,DD),2),mul(P['r'],mul(A,E)),Q(-1))),mul(phi,add(mul(dp(B),DD),mul(A,dp(E)),Q(-1)))),add(mul(dp(phi),add(mul(add(scale(one,2),P['r'],Q(-1)),mul(A,DD)),scale(mul(B,C),4))),mul(phi,total(mul(dp(A),DD),scale(mul(dp(B),C),2),scale(mul(A,dp(DD)),-1)))),add(mul(add(scale(one,4),P['r'],Q(-1)),mul(dp(phi),mul(A,C))),scale(mul(phi,mul(dp(A),C)),2))];assert [cof(J,1,j) for j in range(3)]==[{e[:1]+e[2:]:c for e,c in z.items()} for z in cofs]==[decode(z,n-1) for z in d['w_coefficients']];difference=add(dictionary(phi,g,add(h,mul(P['rho'],power(g,2,n)),Q(-1)),P['r']),J,Q(-1));assert difference==decode(d['target_shear_difference'],n)==mul(add(P['r'],scale(one,4),Q(-1)),mul(P['rho'],mul(power(g,2,n),bracket(phi,g))))
 # For C=0 the polynomial numerator of d(A/(phi^(r-2)*D)) is phi^(r-3)*J_w.
 for rr in [4,5]:
  S=[variable(j,n) for j in range(n)];S[2+names.index('r')]=con(rr,n);S[2+names.index('C')]={};W=compose(cofs[1],S,n);numer=add(mul(dp(A),mul(power(phi,rr-2,n),DD)),mul(A,dp(mul(power(phi,rr-2,n),DD))),Q(-1));assert numer==mul(power(phi,rr-3,n),W)
 # The degree-raising r4 shear needs no bound on its new D: check the exact primitive numerator.
 S=[variable(j,n) for j in range(n)];S[2+names.index('r')]=con(4,n);S[2+names.index('A1')]={};S[2+names.index('C')]=mul(power(P['A0'],2,n),P['rho']);gs,hs=[compose(z,S,n) for z in [g,h]];hs=add(hs,mul(P['rho'],power(gs,2,n)),Q(-1));assert deg(cof(hs,1,2),0)==-1;Ds=cof(hs,1,1);W=cof(dictionary(phi,gs,hs,con(4,n)),1,1);ph={e[:1]+e[2:]:c for e,c in phi.items()};a0={e[:1]+e[2:]:c for e,c in P['A0'].items()};assert mul(ph,W)==scale(mul(a0,derivative(mul(power(ph,2,n-1),Ds),0)),-1)
 # Complete Wronskian degree pairs, plus the exceptional r4 power subtraction.
 allowed={}
 for rr in [4,5]:
  retained=[]
  for m in range(4):
   for nn in range(5):
    if 2*nn!=rr*m:
     if m+nn==1:retained.append([m,nn])
    elif m>0:
     if rr==4 and m==2 and nn==4:continue # eliminated by the independently checked power subtraction below
     retained.append([m,nn])
  allowed[rr]=retained
 assert allowed=={4:[[0,1],[1,0],[1,2]],5:[[0,1],[1,0]]}
 n=6;u,a,b,c,k,t=[variable(j,n) for j in range(n)];G=total(mul(a,power(u,2,n)),mul(b,u),c);R=total(con(1,n),u,power(u,2,n),power(u,3,n));H=add(mul(k,power(G,2,n)),R);wr=lambda H:add(scale(mul(G,derivative(H,0)),2),scale(mul(derivative(G,0),H),4),Q(-1));assert wr(H)==wr(R)
 for degree in range(4):
  RR=mul(t,power(u,degree,n));W=wr(RR);assert deg(W,0)==degree+1 and cof(W,0,degree+1)==scale(mul(variable(0,n-1),variable(4,n-1)),2*(degree-4))
 # Universal graph coefficients after the proved degree reductions.
 d=struct['graph'];names=['r','A','B','alpha','beta','T']+['h'+str(j) for j in range(3)]+['k'+str(j) for j in range(3)];assert d['parameters']==names;n=2+len(names);u,v=variable(0,n),variable(1,n);P={z:variable(2+j,n) for j,z in enumerate(names)};one=con(1,n);G0=add(mul(P['A'],u),P['B']);G1=add(mul(P['alpha'],u),P['beta']);H0=total(*(mul(P['h'+str(j)],power(u,j,n)) for j in range(3)));H1=total(*(mul(P['k'+str(j)],power(u,j,n)) for j in range(3)));G=add(G0,mul(v,G1));H=total(H0,mul(v,H1),mul(P['T'],power(v,2,n)));J=add(mul(add(scale(G,2),mul(v,derivative(G,1))),dp(H)),mul(dp(G),add(mul(P['r'],H),mul(v,derivative(H,1)))),Q(-1));assert [G0,G1,H0,H1,G,H]==[decode(z,n) for z in d['forms']] and J==decode(d['jacobian'],n);expected=[add(scale(mul(G0,dp(H0)),2),mul(P['r'],mul(P['A'],H0)),Q(-1)),total(scale(mul(G1,dp(H0)),3),scale(mul(G0,dp(H1)),2),scale(mul(P['r'],mul(P['alpha'],H0)),-1),scale(mul(add(P['r'],one),mul(P['A'],H1)),-1)),total(scale(mul(G1,dp(H1)),3),scale(mul(add(P['r'],one),mul(P['alpha'],H1)),-1),scale(mul(add(P['r'],scale(one,2)),mul(P['A'],P['T'])),-1)),scale(mul(add(P['r'],scale(one,2)),mul(P['alpha'],P['T'])),-1)];assert [cof(J,1,j) for j in range(4)]==[{e[:1]+e[2:]:c for e,c in z.items()} for z in expected]==[decode(z,n-1) for z in d['v_coefficients']]
 # Affine-G vector field: q has eigenvalue2A, v has eigenvalue-A.
 n=5;u,v,a,b,r=[variable(j,n) for j in range(n)];one=con(1,n);q=add(mul(a,u),mul(b,add(v,one,Q(-1))));act=lambda z:add(mul(add(scale(q,2),mul(b,v)),derivative(z,0)),mul(mul(a,v),derivative(z,1)),Q(-1));assert act(q)==scale(mul(a,q),2) and act(v)==scale(mul(a,v),-1);eigen=0
 for i in range(5):
  for j in range(5-i):
   mon=mul(power(q,i,n),power(v,j,n));assert add(act(mon),mul(a,mul(r,mon)),Q(-1))==mul(a,mul(add(con(2*i-j,n),r,Q(-1)),mon));eigen+=1
 # Complete normalized r5 one-parameter contradiction with a four-term unit.
 d=struct['r5_residual'];n=2;u,A=variable(0,n),variable(1,n);one=con(1,n);ee=scale(add(one,scale(A,5)),Q(-1,2));alpha=scale(add(one,A),Q(-3,2));mu=scale(add(scale(one,5),A),Q(1,2));assert ee=={(0,)+e:c for e,c in decode(d['e'],1).items()} and alpha=={(0,)+e:c for e,c in decode(d['alpha'],1).items()} and mu=={(0,)+e:c for e,c in decode(d['mu'],1).items()};assert total(scale(ee,-2),scale(A,-5),scale(one,-1))=={} and total(A,alpha,mu,scale(one,-1))=={} and total(ee,scale(alpha,-2),scale(mu,-1))=={};assert determinant3([[con(-2,1),{},{}],[{},con(1,1),con(1,1)],[con(1,1),con(-2,1),con(-1,1)]],1)==con(-2,1);G0=add(mul(A,u),one,Q(-1));G1=add(one,mul(alpha,u));H0=add(one,mul(ee,u));H1=scale(power(G1,2,n),-1);assert [G0,G1,H0,H1]==[decode(z,n) for z in d['forms']];E0=total(scale(mul(G0,dp(H0)),2),scale(mul(A,H0),-5),scale(one,-1));E1=total(scale(mul(G1,dp(H0)),3),scale(mul(G0,dp(H1)),2),scale(mul(alpha,H0),-5),scale(mul(A,H1),-6));assert E0==decode(d['E0'],n) and E1==decode(d['E1'],n);residuals=[cof(z,0,1) for z in [E0,E1]];assert residuals==[decode(z,1) for z in d['residuals']];multipliers=[decode(z,1) for z in d['multipliers']];assert total(*(mul(x,y) for x,y in zip(residuals,multipliers)))==con(1,1);assert sum(map(len,multipliers))==d['multiplier_terms']==4
 # All four normal-form families, with five complete generic parameter charts.
 controls=[]
 for r,family,chart in [(4,'I','free'),(4,'II','free'),(5,'I','L=0'),(5,'I','theta=0'),(5,'II','free')]:
  n=6;x=[variable(j,n) for j in range(3)];L,M,N=[variable(j,n) for j in range(3,6)];a,b,c=x
  if family=='I':
   if chart=='L=0':L={}
   if chart=='theta=0':N={}
   shear=add(mul(M,mul(power(a,r-4,n),power(b,r-2,n))),mul(N,mul(power(a,r-2,n),power(b,r-1,n))));third=add(c,shear);second=add(b,mul(L,mul(power(a,r-2,n),third)));F=[a,second,third];bi=add(b,mul(L,mul(power(a,r-2,n),c)),Q(-1));ci=add(c,add(mul(M,mul(power(a,r-4,n),power(bi,r-2,n))),mul(N,mul(power(a,r-2,n),power(bi,r-1,n)))),Q(-1));inv=[a,bi,ci]
  elif r==4:
   second=add(b,mul(L,mul(power(a,2,n),c)));F=[a,second,add(c,mul(M,power(second,2,n)))];ci=add(c,mul(M,power(b,2,n)),Q(-1));inv=[a,add(b,mul(L,mul(power(a,2,n),ci)),Q(-1)),ci]
  else:
   quad=add(mul(M,power(b,2,n)),mul(a,c),Q(-1));second=add(b,mul(L,mul(power(a,2,n),quad)),Q(-1));third=total(c,scale(mul(M,mul(L,mul(a,mul(b,quad)))),-2),mul(M,mul(power(L,2,n),mul(power(a,3,n),power(quad,2,n)))));F=[a,second,third];inv=[a,add(b,mul(L,mul(power(a,2,n),quad))),total(c,scale(mul(M,mul(L,mul(a,mul(b,quad)))),2),mul(M,mul(power(L,2,n),mul(power(a,3,n),power(quad,2,n)))))];assert compose(quad,F+[variable(j,n) for j in range(3,6)],n)==quad
  params=[variable(j,n) for j in range(3,6)];assert [compose(z,inv+params,n) for z in F]==x and [compose(z,F+params,n) for z in inv]==x and determinant3([[derivative(z,j) for j in range(3)] for z in F],n)==con(1,n);degree=max(sum(e[:3]) for z in F for e in z);assert degree<=7 and all(e[0]-2*e[1]-r*e[2]==target for z,target in zip(F,[1,-2,-r]) for e in z);controls.append({'r':r,'family':family,'chart':chart,'degree':degree,'both_generic_inverse_compositions':True});print('INVERSE PASS',r,family,chart,degree,flush=True)
 oldcp=HERE/'weight-resonance7/mixed-lift-checkpoint.json';old=json.loads(oldcp.read_text());timeouts=[z for z in old['cases'] if z['status']=='timeout'];assert len(timeouts)==4 and {tuple(z['magnitudes']) for z in timeouts}=={(1,2,4),(1,2,5),(3,1,2),(4,1,2)};frontier=HERE/'weight-resonance7/validation.json';oldfront=json.loads(frontier.read_text());assert oldfront['status']=='pass' and len(oldfront['remaining_weight_types'])==14;remaining=[z for z in oldfront['remaining_weight_types'] if z['mixed_positive_negative_magnitudes'] not in [[1,2,4],[1,2,5]]];assert len(remaining)==12
 proof=HERE/'WEIGHTED_124_125_EXCLUSION.md';receipt={'status':'pass','claim':'KCUB37','checker_sha256':sha(Path(__file__)),'producer_sha256':sha(HERE/'weighted124125.py'),'arithmetic_sha256':sha(HERE/'verify.py'),'checkpoint_sha256':sha(cp),'structural_identities_sha256':sha(sp),'model_hashes':model_hashes,'source_atoms':atoms,'native_coefficient_generators':native,'complete_graph_transformations':2,'univariate_coefficients':3,'weighted_square_shear_identity':True,'power_relation_numerators':2,'integer_boundary_power_subtraction_leading_cases':4,'graph_coefficients':4,'affine_vectorfield_eigenvalues':eigen,'constant_pivot_determinant':-2,'one_variable_unit_multiplier_terms':4,'generic_inverse_charts':controls,'quadratic_invariant_flow_checked':True,'old_four_timeout_checkpoint_sha256':sha(oldcp),'old_timeouts_preserved':True,'frontier_KCUB36_receipt_sha256':sha(frontier),'remaining_diagonal_weight_improvement_types':remaining,'no_external_theorem_imported_for_two_weight_closures':True,'KCUB36_import_only_for_frontier_reduction':True,'global_profile_counts_and_positive_bounds_unchanged':True,'goal_complete':False,'trust_boundary':'Independent rational full-source and identity/inverse reconstruction plus written field/unit, degree, complete graph and univariate case proofs. Not formalized or externally reviewed.'}
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (OUT/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS complete124/125 closure:24 atoms,30 native equations,4-term graph unit,5 generic inverse charts; diagonal-weight raw<=12 frontier14->12; global goal open',flush=True)
if __name__=='__main__':main()
