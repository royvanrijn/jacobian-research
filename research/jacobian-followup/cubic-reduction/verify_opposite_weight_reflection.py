#!/usr/bin/env python3
"""Independent polynomial checks for the complete primitive opposite-weight class."""
from pathlib import Path
from fractions import Fraction as Q
from math import gcd
from itertools import product
import json,hashlib
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def total(*ps):
 out={}
 for p in ps:out=add(out,p)
 return out

def power(p,k,n):
 out={(0,)*n:Q(1)}
 for j in range(k):out=mul(out,p)
 return out

def bracket(f,g):return add(mul(derivative(f,0),derivative(g,1)),mul(derivative(f,1),derivative(g,0)),Q(-1))
def dictionary(f,g,h,k,u,v,m,q):
 R=add(mul(u,h),mul(v,k));return add(mul(g,total(mul(f,derivative(R,0)),mul(R,derivative(f,0)),mul(q,mul(v,bracket(R,f))))),mul(m,mul(v,add(mul(f,bracket(R,g)),mul(R,bracket(f,g))))))
def module(m,q,support,ksupport):
 names=[p+str(i) for p in ['f','g','h'] for i in range(len(support))]+['k'+str(i) for i in range(len(ksupport))];n=3+len(names);X=[variable(i,n) for i in range(3)];P={name:variable(3+i,n) for i,name in enumerate(names)};F=list(X);seen=[]
 # Full independent source-monomial enumeration, rather than the producer's invariant expansion.
 for out,wgt in enumerate([m,-q,-m]):
  for e in product(range(8),repeat=3):
   if not 2<=sum(e)<=7 or m*(e[0]-e[2])-q*e[1]!=wgt:continue
   if out==0:assert e[1]%m==0;key=(e[2],e[1]//m);name='f'+str(support.index(key))
   elif out==1:assert (e[1]-1)%m==0;key=(e[2],(e[1]-1)//m);name='g'+str(support.index(key))
   elif e[2]:assert e[1]%m==0;key=(e[2]-1,e[1]//m);name='h'+str(support.index(key))
   else:assert e[1]%m==0;key=e[1]//m-1;name='k'+str(ksupport.index(key))
   assert name in P;F[out]=add(F[out],mul({e+(0,)*len(names):Q(1)},P[name]));seen.append(name)
 assert len(seen)==len(names) and set(seen)==set(names)
 return names,F

def invariant_forms(support,ksupport,names):
 n=2+len(names);u,v=variable(0,n),variable(1,n);P={name:variable(2+i,n) for i,name in enumerate(names)};one={(0,)*n:Q(1)};forms=[total(one,*(mul(P[p+str(j)],mul(power(u,i,n),power(v,k,n))) for j,(i,k) in enumerate(support))) for p in ['f','g','h']];K=total(*(mul(P['k'+str(i)],power(v,j,n)) for i,j in enumerate(ksupport)));return u,v,P,forms,K

def main():
 folder=HERE/'opposite-weight-reflection';cp=folder/'checkpoint.json';check=json.loads(cp.read_text());assert check['status']=='complete' and check['source_sha256']==sha(HERE/'opposite_weight_reflection.py');pairs=[(m,s-m) for s in range(3,9) for m in range(2,s) if gcd(m,s-m)==1];assert len(pairs)==14 and len(check['models'])==14
 # Formal weight parameters in the invariant dictionary, not a finite weight-value probe.
 up=folder/'universal-identities.json';universal=json.loads(up.read_text());assert universal['status']=='complete-primary-universal-identities' and universal['source_sha256']==sha(HERE/'opposite_weight_structure.py');names=['m','q','A','B','C','D','E','G','H','I','K','L','M'];assert universal['parameter_names']==names;n=2+len(names);u,v=variable(0,n),variable(1,n);P={name:variable(2+i,n) for i,name in enumerate(names)};one={(0,)*n:Q(1)}
 f=total(one,mul(P['A'],v),mul(P['B'],mul(u,v)),mul(P['C'],power(v,2,n)));g=total(one,mul(P['D'],v),mul(P['E'],mul(u,v)),mul(P['G'],power(v,2,n)));h=total(one,mul(P['H'],v),mul(P['I'],mul(u,v)),mul(P['K'],power(v,2,n)));k=add(P['L'],mul(P['M'],v));J=dictionary(f,g,h,k,u,v,P['m'],P['q']);assert J==decode(universal['dictionary'],n)
 pp={name:variable(i,len(names)) for i,name in enumerate(names)};o={(0,)*len(names):Q(1)};cc=lambda a:scale(o,Q(a));m,q=pp['m'],pp['q'];B,E,I=pp['B'],pp['E'],pp['I'];expected=[total(mul(add(q,cc(2)),B),mul(add(m,cc(1)),E),scale(I,2)),total(mul(total(m,q,cc(2)),mul(B,E)),mul(add(q,cc(3)),mul(B,I)),mul(add(m,cc(2)),mul(E,I))),mul(total(m,q,cc(3)),mul(B,mul(E,I)))];diagonal=[]
 for j in [1,2,3]:diagonal.append({e[2:]:c for e,c in J.items() if e[:2]==(j,j)})
 assert diagonal==expected==[decode(p,len(names)) for p in universal['diagonal_coefficients']]
 assert [a['mask'] for a in universal['support_masks']]==list(range(1,8))
 for rec in universal['support_masks']:
  mask=rec['mask'];active=[int(bool(mask&(1<<j))) for j in range(3)];degree=sum(active);S=[variable(j,len(names)) for j in range(len(names))]
  for name,on in zip(['B','E','I'],active):
   if not on:S[names.index(name)]={}
  coef=total(cc(1+active[0]+active[2]),scale(q,active[0]),scale(m,active[1]));wanted=coef
  for p,on in zip([B,E,I],active):
   if on:wanted=mul(wanted,p)
  assert compose(diagonal[degree-1],S,len(names))==wanted and rec['degree']==degree and decode(rec['coefficient'],len(names))==coef
  # The verified scalar has positive integer coefficients and a positive constant; nonzero for all allowed weights.
  assert all(c>0 for c in coef.values()) and coef.get((0,)*len(names),0)>=1
 S=[variable(j,n) for j in range(n)]
 for name in ['B','E','I']:S[2+names.index(name)]={}
 pure=compose(J,S,n);ff,gg,hh=[compose(p,S,n) for p in [f,g,h]];rhs=mul(hh,add(mul(gg,add(ff,mul(P['q'],mul(v,derivative(ff,1))))),mul(P['m'],mul(v,mul(ff,derivative(gg,1))))));assert pure==rhs and pure=={(0,)+e:c for e,c in decode(universal['pure_v_factorization'],n-1).items()}
 assert len(universal['leading_cases'])==9
 for rec in universal['leading_cases']:
  df,dg=rec['degrees'];F=[one,add(one,mul(P['A'],v)),ff][df];G=[one,add(one,mul(P['D'],v)),gg][dg];expr=add(mul(G,add(F,mul(P['q'],mul(v,derivative(F,1))))),mul(P['m'],mul(v,mul(F,derivative(G,1)))));coef={e[2:]:c for e,c in expr.items() if e[:2]==(0,df+dg)};leading_f=[o,pp['A'],pp['C']][df];leading_g=[o,pp['D'],pp['G']][dg];wanted=mul(total(o,scale(q,df),scale(m,dg)),mul(leading_f,leading_g));assert coef==wanted==decode(rec['coefficient'],len(names))
 model_hashes={};full_atoms=0;reduced_atoms=0;controls=[]
 for meta,(mw,qw) in zip(check['models'],pairs):
  name='m%dq%d'%(mw,qw);assert meta['name']==name;mp=folder/meta['input'];model=json.loads(mp.read_text());assert model['source_sha256']==sha(HERE/'opposite_weight_reflection.py') and meta['input_sha256']==sha(mp) and model['weights']==[mw,-qw,-mw];sm=mw+qw;support=[(i,j) for j in range(3) for i in range(4) if 0<2*i+sm*j<=6];ksupport=[j for j in range(3) if sm-1+sm*j<=7];assert model['support']==list(map(list,support)) and model['k_support']==ksupport
  fnames,F=module(mw,qw,support,ksupport);nf=3+len(fnames);assert model['parameter_names']==fnames and F==[decode(p,nf) for p in model['outputs']];full_atoms+=len(fnames)
  # Direct b=0 source Jacobian, with no primary determinant expression trusted.
  Jb=[[{e:c for e,c in derivative(p,j).items() if e[1]==0} for j in range(3)] for p in F];boundary=determinant3(Jb,nf);assert boundary=={(e[0],0)+e[1:]:c for e,c in decode(model['boundary_determinant'],nf-1).items()}
  nu=1+len(fnames);U=variable(0,nu);po={name:variable(1+i,nu) for i,name in enumerate(fnames)};ou={(0,)*nu:Q(1)};parts=[total(ou,*(mul(po[p+str(j)],power(U,i,nu)) for j,(i,kp) in enumerate(support) if kp==0)) for p in ['f','g','h']];bp=mul(parts[1],derivative(mul(U,mul(parts[0],parts[2])),0));assert bp==decode(model['boundary_product'],nu);xx=[variable(j,nf) for j in range(3)];lift=[mul(xx[0],xx[2])]+[variable(3+j,nf) for j in range(len(fnames))];assert compose(bp,lift,nf)==boundary
  reduced=[e for e in support if e[1]];assert model['reduced_support']==list(map(list,reduced))
  # Rebuild the boundary-reduced maps from invariant forms (the full source enumeration above proves support coverage).
  rnames=[p+str(j) for p in ['f','g','h'] for j in range(len(reduced))]+['k'+str(j) for j in range(len(ksupport))];ui,vi,rp,rforms,rk=invariant_forms(reduced,ksupport,rnames);nr=3+len(rnames);xr=[variable(j,nr) for j in range(3)];us=mul(xr[0],xr[2]);vs=mul(power(xr[0],qw,nr),power(xr[1],mw,nr));isub=[us,vs]+[variable(3+j,nr) for j in range(len(rnames))];RF=[mul(xr[0],compose(rforms[0],isub,nr)),mul(xr[1],compose(rforms[1],isub,nr)),add(mul(xr[2],compose(rforms[2],isub,nr)),mul(mul(power(xr[0],qw-1,nr),power(xr[1],mw,nr)),compose(rk,isub,nr)))];assert rnames==model['reduced_parameter_names'] and RF==[decode(p,nr) for p in model['reduced_outputs']];reduced_atoms+=len(rnames)
  ni=2+len(rnames);oi={(0,)*ni:Q(1)};D=dictionary(*rforms,rk,ui,vi,scale(oi,mw),scale(oi,qw));assert D==decode(model['invariant_dictionary'],ni);actual=determinant3([[derivative(p,j) for j in range(3)] for p in RF],nr);assert actual==compose(D,isub,nr)
  # Every reduced source model is an exact specialization of the universal eleven-coefficient form.
  rcoeff={}
  for prefix,letters in [('f',['A','B','C']),('g',['D','E','G']),('h',['H','I','K'])]:
   for letter,key in zip(letters,[(0,1),(1,1),(0,2)]):rcoeff[letter]=rp[prefix+str(reduced.index(key))] if key in reduced else {}
  rcoeff['L']=rp['k'+str(ksupport.index(0))] if 0 in ksupport else {};rcoeff['M']=rp['k'+str(ksupport.index(1))] if 1 in ksupport else {};proto_sub=[ui,vi,scale(oi,mw),scale(oi,qw)]+[rcoeff[z] for z in names[2:]];assert compose(J,proto_sub,ni)==D
  # All classified shear coefficients are allowed; choose every permitted one=1 and check the actual inverse.
  x=[variable(j,3) for j in range(3)];vv=mul(power(x[0],qw,3),power(x[1],mw,3));shear=mul(mul(power(x[0],qw-1,3),power(x[1],mw,3)),total(*(power(vv,j,3) for j in ksupport)));Hmap=[x[0],x[1],add(x[2],shear)];inverse=[x[0],x[1],add(x[2],shear,Q(-1))];assert [compose(p,inverse,3) for p in Hmap]==x and [compose(p,Hmap,3) for p in inverse]==x and determinant3([[derivative(p,j) for j in range(3)] for p in Hmap],3)=={(0,)*3:Q(1)};degree=max(sum(e) for p in Hmap for e in p);assert degree<=7;controls.append({'weights':[mw,-qw,-mw],'degree':degree,'shear_terms':len(ksupport),'both_inverse_compositions':True});model_hashes[name]=sha(mp);print('MODEL PASS',name,len(fnames),'->',len(rnames),'inverse degree',degree,flush=True)
 assert (full_atoms,reduced_atoms)==(169,43)
 receipt={'status':'pass','claim':'KCUB33','checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'model_producer_sha256':sha(HERE/'opposite_weight_reflection.py'),'universal_producer_sha256':sha(HERE/'opposite_weight_structure.py'),'universal_identities_sha256':sha(up),'checkpoint_sha256':sha(cp),'model_hashes':model_hashes,'primitive_small_weight_pairs':pairs,'full_nonlinear_atoms':full_atoms,'boundary_reduced_atoms':reduced_atoms,'complete_boundary_and_invariant_dictionaries':14,'universal_diagonal_identities':3,'all_nonzero_support_masks':7,'univariate_degree_cases':9,'classified_triangular_inverse_controls':controls,'large_primary_unit_certificates_not_used':True,'no_plane_or_external_theorem_imported':True,'normalized_classification':'(a,b,c+a^(q-1)*b^m*k(a^q*b^m)), with gcd(m,q)=1,m>=2,q>=1 and every retained monomial of degree<=7','m_equals_one_not_covered':True,'large_sum_at_least9_covered_by_written_module_degree_bound':True,'global_profile_counts_unchanged':[62,12],'positive_counterexample_bounds_unchanged':True,'goal_complete':False,'trust_boundary':'Exact finite polynomial reconstruction plus written primitive-module, univariate-unit and positive-integer leading-coefficient arguments. No formal or external-review claim.'}
 proof=HERE/'OPPOSITE_WEIGHT_CLASSIFICATION.md'
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (folder/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS complete primitive(m,-q,-m),m>=2 throughdegree7 classification: all triangular;169 atoms,14 boundary dictionaries,3 universal coefficients,7 support masks,9 degree cases; global goal open',flush=True)
if __name__=='__main__':main()
