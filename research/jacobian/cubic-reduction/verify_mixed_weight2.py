#!/usr/bin/env python3
"""Independent rational replay of the complete weight(2,-1,-r),r>=4 class."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
import json,hashlib
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3
HERE=Path(__file__).resolve().parent;OUT=HERE/'mixed-weight2'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def total(*ps):
 o={}
 for p in ps:o=add(o,p)
 return o

def power(p,k,n):
 o={(0,)*n:Q(1)}
 for _ in range(k):o=mul(o,p)
 return o

def con(c,n):return {(0,)*n:Q(c)} if c else {}
def cof(p,i,k):return {e[:i]+e[i+1:]:c for e,c in p.items() if e[i]==k}
def bracket(f,g):return add(mul(derivative(f,0),derivative(g,1)),mul(derivative(f,1),derivative(g,0)),Q(-1))
def dictionary(f,g,h,r):return total(mul(f,bracket(g,h)),scale(mul(g,bracket(f,h)),Q(1,2)),scale(mul(r,mul(h,bracket(g,f))),Q(1,2)))
def main():
 cp=OUT/'checkpoint.json';check=json.loads(cp.read_text());assert check['status']=='complete' and check['source_sha256']==sha(HERE/'mixed_weight2.py') and check['r_values']==list(range(4,14));models={};hashes={};atoms=0;native=0;boundaries=0
 for meta in check['models']:
  r=meta['r'];mp=OUT/meta['input'];d=json.loads(mp.read_text());assert d['source_sha256']==sha(HERE/'mixed_weight2.py') and meta['input_sha256']==sha(mp);basis=[[e for e in product(range(8),repeat=3) if 2<=sum(e)<=7 and 2*e[0]-e[1]-r*e[2]==t] for t in [2,-1,-r]];assert d['basis']==[[list(e) for e in row] for row in basis];names=['c'+str(j+1) for j in range(sum(map(len,basis)))];assert d['parameter_names']==names;n=3+len(names);x=[variable(j,n) for j in range(3)];F=list(x);pos=0
  for j,row in enumerate(basis):
   for e in row:F[j]=add(F[j],mul({e+(0,)*len(names):Q(1)},variable(3+pos,n)));pos+=1
  assert F==[decode(z,n) for z in d['source_outputs']];ni=2+len(names);one=con(1,ni);u,v=variable(0,ni),variable(1,ni);sub=[one,u,v]+[variable(2+j,ni) for j in range(len(names))];forms=[compose(z,sub,ni) for z in F];assert forms==[decode(z,ni) for z in d['double_cover_outputs']];D=dictionary(*forms,con(r,ni));assert D==decode(d['quotient_jacobian'],ni)
  # Recover a exponent from invariant weight zero, verifying descent of the half-cover.
  sourceJ={}
  for e,c in D.items():
   assert (e[0]+r*e[1])%2==0;ee=((e[0]+r*e[1])//2,e[0],e[1])+e[2:];sourceJ[ee]=c
  actual=determinant3([[derivative(z,j) for j in range(3)] for z in F],n);assert sourceJ==actual
  equations=[{e[2:]:c for e,c in add(D,one,Q(-1)).items() if e[:2]==mon} for mon in sorted({e[:2] for e in add(D,one,Q(-1))},reverse=True)];assert equations==[decode(z,len(names)) for z in d['native_generators']];atoms+=pos;native+=len(equations);models[r]=(d,basis,F,forms,equations);hashes[mp.name]=sha(mp)
  if r%2==0:
   t=r//2;nf=1+len(names);f0=cof(forms[0],0,0);g0=cof(forms[1],0,1);hh=cof(forms[2],0,0);assert all(e[0]>=1 for e in hh);h0={(e[0]-1,)+e[1:]:c for e,c in hh.items()};V=variable(0,nf);bp=mul(g0,total(mul(h0,add(f0,scale(mul(V,derivative(f0,0)),t))),mul(V,mul(f0,derivative(h0,0)))));lift=[mul(power(x[0],t,n),x[2])]+[variable(3+j,n) for j in range(len(names))];assert compose(bp,lift,n)=={e:c for e,c in actual.items() if e[1]==0};boundaries+=1
  if r%2:
   fs={(0,0),(2,0),(4,0)}|({(1,1)} if r<=7 else set());gs={(1,0),(3,0),(5,0),(0,1)}|({(2,1)} if r<=7 else set());hs={(0,1),(2,1),(4,1)}|({(1,2),(r,0)} if r<=7 else set());assert [{e[:2] for e in z} for z in forms]==[fs,gs,hs];assert all((e[0]+e[1])%2==parity for z,parity in zip(forms,[0,1,1]) for e in z)
  print('MODEL PASS',r,pos,len(equations),flush=True)
 assert (atoms,native,boundaries)==(94,142,5)
 up=OUT/'structural-identities.json';struct=json.loads(up.read_text());assert struct['status']=='complete-primary-structural-identities' and struct['source_sha256']==sha(HERE/'mixed_weight2_structure.py')
 # Reconstruct the even r4 boundary-reduced universal source determinant from scratch.
 d=struct['even_r4'];names=['A','B','C','D','E','G','H','I','K','L','M'];assert d['parameters']==names;n=2+len(names);U,V=variable(0,n),variable(1,n);P={z:variable(2+j,n) for j,z in enumerate(names)};one=con(1,n);f=total(one,mul(P['A'],U),mul(P['B'],power(U,2,n)),mul(P['C'],mul(U,V)));g=total(one,mul(P['D'],U),mul(P['E'],power(U,2,n)),mul(P['G'],mul(U,V)));h=total(one,mul(P['H'],U),mul(P['I'],power(U,2,n)),mul(P['K'],mul(U,V)));k=add(P['L'],mul(P['M'],U));assert [f,g,h,k]==[decode(z,n) for z in d['forms']];ns=3+len(names);x=[variable(j,ns) for j in range(3)];sub=[mul(x[0],power(x[1],2,ns)),mul(power(x[0],2,ns),x[2])]+[variable(3+j,ns) for j in range(len(names))];F=[mul(x[0],compose(f,sub,ns)),mul(x[1],compose(g,sub,ns)),add(mul(x[2],compose(h,sub,ns)),mul(power(x[1],4,ns),compose(k,sub,ns)))];assert F==[decode(z,ns) for z in d['source_outputs']];Jsource=determinant3([[derivative(z,j) for j in range(3)] for z in F],ns);J={}
 for e,c in Jsource.items():assert e[1]%2==0 and e[0]==e[1]//2+2*e[2];J[(e[1]//2,e[2])+e[3:]]=c
 assert J==decode(d['invariant_jacobian'],n);nn=len(names);PP={z:variable(j,nn) for j,z in enumerate(names)};diag=[{e[2:]:c for e,c in J.items() if e[:2]==(j,j)} for j in [1,2,3]];wanted=[total(scale(PP['C'],4),scale(PP['G'],3),scale(PP['K'],2)),total(scale(mul(PP['C'],PP['G']),6),scale(mul(PP['C'],PP['K']),5),scale(mul(PP['G'],PP['K']),4)),scale(mul(PP['C'],mul(PP['G'],PP['K'])),7)];assert diag==wanted==[decode(z,nn) for z in d['diagonal_coefficients']]
 for mask in range(1,8):
  bits=[int(bool(mask&(1<<j))) for j in range(3)];subp=[variable(j,nn) for j in range(nn)];target=con(1+3*bits[0]+2*bits[1]+bits[2],nn)
  for z,on in zip(['C','G','K'],bits):
   if on:target=mul(target,PP[z])
   else:subp[names.index(z)]={}
  assert compose(diag[sum(bits)-1],subp,nn)==target
 S=[variable(j,n) for j in range(n)]
 for z in ['C','G','K']:S[2+names.index(z)]={}
 pure=compose(J,S,n);ff,gg,hh=[compose(z,S,n) for z in [f,g,h]];rhs=mul(hh,add(mul(gg,add(ff,mul(U,derivative(ff,0)))),scale(mul(U,mul(ff,derivative(gg,0))),2)));assert pure==rhs=={e[:1]+(0,)+e[1:]:c for e,c in decode(d['pure_U_factorization'],n-1).items()}
 # Every complete even model reduces exactly to this universal support or its univariate form.
 for r in range(4,14,2):
  md,basis,full,_,_=models[r];n0=3+len(md['parameter_names']);subp=[variable(j,n0) for j in range(n0)];pos=0;coeffs={}
  for row,es in enumerate(basis):
   for e in es:
    coeffs[(row,e)]=variable(3+pos,n0)
    if (row==0 and e[1]==0) or (row==1 and e[1]==1) or (row==2 and e[1]==0):subp[3+pos]={}
    pos+=1
  RF=[compose(z,subp,n0) for z in full];x=[variable(j,n0) for j in range(3)];us=mul(x[0],power(x[1],2,n0));vs=mul(power(x[0],r//2,n0),x[2]);spec={}
  def take(row,e):return coeffs.get((row,e),{})
  for prefix,row,shift in [('f',0,(1,0,0)),('g',1,(0,1,0)),('h',2,(0,0,1))]:
   spec[prefix]=total(con(1,n0),*(mul(take(row,(shift[0]+i+(r//2)*j,shift[1]+2*i,shift[2]+j)),mul(power(us,i,n0),power(vs,j,n0))) for i,j in [(1,0),(2,0),(1,1)] if 3*i+(r//2+1)*j<=6))
  shear=total(*(mul(take(2,(j,r+2*j,0)),mul(power(x[1],r,n0),power(us,j,n0))) for j in [0,1] if r+3*j<=7));expected=[mul(x[0],spec['f']),mul(x[1],spec['g']),add(mul(x[2],spec['h']),shear)];assert RF==expected
 # Positive leading scalars for the even boundary and the residual univariate determinant.
 n=4;z,t,L,R=[variable(j,n) for j in range(n)];one=con(1,n);leading=0
 for df in range(3):
  for dh in range(3):
   F=one if df==0 else add(one,mul(L,power(z,df,n)));H=one if dh==0 else add(one,mul(R,power(z,dh,n)));exp=add(mul(H,add(F,mul(t,mul(z,derivative(F,0))))),mul(z,mul(F,derivative(H,0))));lead=cof(exp,0,df+dh);scalar=add(con(1+dh,3),scale(variable(0,3),df));lcf=con(1,3) if df==0 else variable(1,3);lch=con(1,3) if dh==0 else variable(2,3);assert lead==mul(scalar,mul(lcf,lch));leading+=1
 for df in range(3):
  for dg in range(3):
   F=one if df==0 else add(one,mul(L,power(z,df,n)));G=one if dg==0 else add(one,mul(R,power(z,dg,n)));exp=add(mul(G,add(F,mul(z,derivative(F,0)))),scale(mul(z,mul(F,derivative(G,0))),2));lead=cof(exp,0,df+dg);assert lead==scale(mul(con(1,3) if df==0 else variable(1,3),con(1,3) if dg==0 else variable(2,3)),1+df+2*dg);leading+=1
 # Exact original odd nu1 unit identities, with every native generator independently specialized.
 lc=OUT/'odd-lift-checkpoint.json';lifts=json.loads(lc.read_text());assert lifts['status']=='complete' and lifts['source_sha256']==sha(HERE/'certify_mixed_weight2.py');unit_terms=0;unit_gens=0;unit_hashes={}
 for rec in lifts['cases']:
  r=rec['r'];assert rec['status']=='exact-primary-unit';d,basis,_,_,eqs=models[r];nu_index=basis[0].index(((r+3)//2,1,1));np=len(d['parameter_names']);S=[]
  for j in range(np):S.append(con(1,np-1) if j==nu_index else variable(j-int(j>nu_index),np-1))
  gens=[compose(z,S,np-1) for z in eqs];gens=[z for z in gens if z];ip=OUT/rec['input'];inp=json.loads(ip.read_text());assert inp['full_model_sha256']==hashes['r'+str(r)+'.json'] and inp['normalized_parameter']==d['parameter_names'][nu_index];assert inp['parameter_names']==[z for j,z in enumerate(d['parameter_names']) if j!=nu_index];assert gens==[decode(z,np-1) for z in inp['generators']];certp=OUT/('r'+str(r)+'-nu1-unit.json');cert=json.loads(certp.read_text());assert sha(ip)==rec['input_sha256']==cert['input_sha256'] and sha(certp)==rec['certificate_sha256'] and cert['source_sha256']==sha(HERE/'certify_mixed_weight2.py');script=OUT/('r'+str(r)+'-nu1.sing');assert cert['singular_source_sha256']==rec['script_sha256']==sha(script);multipliers=[decode(z,np-1) for z in cert['multipliers']];assert len(multipliers)==len(gens);assert total(*(mul(a,b) for a,b in zip(multipliers,gens)))==con(1,np-1);terms=sum(map(len,multipliers));assert terms==rec['multiplier_terms']==9 and max(sum(e) for z in multipliers for e in z)==rec['max_multiplier_degree']==1;unit_terms+=terms;unit_gens+=len(gens);unit_hashes[certp.name]=sha(certp);print('UNIT PASS',r,terms,flush=True)
 assert (unit_terms,unit_gens)==(18,34)
 # Odd univariate first quotient: formal-r coefficient identities, including C=N*u.
 d=struct['odd_univariate'];names=['r','mu','P','A0','A1','B3','B5','N','H2','H4','M5','M7'];assert d['parameters']==names;n=2+len(names);u,v=variable(0,n),variable(1,n);P={z:variable(2+j,n) for j,z in enumerate(names)};one=con(1,n);phi=total(one,mul(P['mu'],power(u,2,n)),mul(P['P'],power(u,4,n)));A=add(P['A0'],mul(P['A1'],power(u,2,n)));B=total(u,mul(P['B3'],power(u,3,n)),mul(P['B5'],power(u,5,n)));C=mul(P['N'],u);D=total(one,mul(P['H2'],power(u,2,n)),mul(P['H4'],power(u,4,n)));E=add(mul(P['M5'],power(u,5,n)),mul(P['M7'],power(u,7,n)));g=add(mul(A,v),B);h=total(mul(C,power(v,2,n)),mul(D,v),E);assert [phi,A,B,C,D,E,g,h]==[decode(z,n) for z in d['forms']];J=dictionary(phi,g,h,P['r']);assert J==decode(d['univariate_dictionary'],n);dp=lambda z:derivative(z,0)
 expected=[add(scale(mul(dp(phi),add(mul(B,D),mul(P['r'],mul(A,E)),Q(-1))),Q(1,2)),mul(phi,add(mul(dp(B),D),mul(A,dp(E)),Q(-1)))),add(scale(mul(dp(phi),add(mul(add(one,P['r'],Q(-1)),mul(A,D)),scale(mul(B,C),2))),Q(1,2)),mul(phi,total(mul(dp(A),D),scale(mul(dp(B),C),2),scale(mul(A,dp(D)),-1)))),add(scale(mul(add(scale(one,2),P['r'],Q(-1)),mul(dp(phi),mul(A,C))),Q(1,2)),mul(phi,add(scale(mul(dp(A),C),2),mul(A,dp(C)),Q(-1))))];assert [cof(J,1,j) for j in range(3)]==[{e[:1]+e[2:]:c for e,c in z.items()} for z in expected]==[decode(z,n-1) for z in d['v_coefficients']];S=[variable(j,n) for j in range(n)]
 for z in ['mu','P']:S[2+names.index(z)]={}
 top=cof(compose(J,S,n),1,2);target=mul(P['N'],add(scale(mul(P['A1'],power(u,2,n)),3),P['A0'],Q(-1)));assert top=={e[:1]+e[2:]:c for e,c in target.items()}==decode(d['constant_first_top'],n-1)
 # Complete generic classified source forms and two-sided inverses, including LM=0 at r7.
 controls=[]
 for r in range(4,14):
  charts=['single'] if r!=7 else ['L=0','M=0']
  for chart in charts:
   n=6;x=[variable(j,n) for j in range(3)];L,M,N=[variable(j,n) for j in range(3,6)]
   if r%2==0:
    if r!=4:N={}
    if r>=8:L={}
    shear=mul(power(x[1],r,n),add(L,mul(N,mul(x[0],power(x[1],2,n)))));F=[x[0],x[1],add(x[2],shear)];inv=[x[0],x[1],add(x[2],shear,Q(-1))]
   else:
    if r>=9 or chart=='M=0':M={}
    if chart=='L=0':L={}
    third=add(x[2],mul(M,power(x[1],r,n)));second=add(x[1],mul(L,mul(power(x[0],(r-1)//2,n),third)));F=[x[0],second,third];bi=add(x[1],mul(L,mul(power(x[0],(r-1)//2,n),x[2])),Q(-1));inv=[x[0],bi,add(x[2],mul(M,power(bi,r,n)),Q(-1))]
   params=[variable(j,n) for j in range(3,6)];assert [compose(z,inv+params,n) for z in F]==x and [compose(z,F+params,n) for z in inv]==x and determinant3([[derivative(z,j) for j in range(3)] for z in F],n)==con(1,n);deg=max(sum(e[:3]) for z in F for e in z);assert deg<=7 and all(2*e[0]-e[1]-r*e[2]==target for z,target in zip(F,[2,-1,-r]) for e in z);controls.append({'r':r,'chart':chart,'degree':deg,'generic_both_inverse_compositions':True})
 assert len(controls)==11
 proof=HERE/'MIXED_WEIGHT2_CLASSIFICATION.md';receipt={'status':'pass','claim':'KCUB35','checker_sha256':sha(Path(__file__)),'arithmetic_sha256':sha(HERE/'verify.py'),'model_producer_sha256':sha(HERE/'mixed_weight2.py'),'structural_producer_sha256':sha(HERE/'mixed_weight2_structure.py'),'unit_producer_sha256':sha(HERE/'certify_mixed_weight2.py'),'checkpoint_sha256':sha(cp),'structural_identities_sha256':sha(up),'unit_checkpoint_sha256':sha(lc),'model_hashes':hashes,'unit_hashes':unit_hashes,'complete_models':10,'source_atoms':atoms,'native_coefficient_generators':native,'even_boundary_factorizations':boundaries,'even_diagonal_coefficients':3,'even_nonzero_support_masks':7,'even_leading_cases':leading,'odd_nonzero_nu_units':2,'unit_generators':unit_gens,'unit_multiplier_monomials':unit_terms,'max_unit_multiplier_degree':1,'odd_univariate_coefficients':3,'odd_constant_first_top_identity':True,'generic_inverse_charts':controls,'all_r_at_least14_identity_by_written_degree_bound':True,'no_external_theorem_imported':True,'old_timeouts_untouched':True,'positive_bounds_and_full_profile_counts_unchanged':True,'goal_complete':False,'trust_boundary':'Independent finite rational identities and actual source determinants plus written all-r module, polynomial-unit, rational derivative, root and field-support arguments; not formalized or externally reviewed.'}
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (OUT/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS all weights(2,-1,-r),r>=4 throughdegree7:94 atoms,142 native generators,18 unit terms,11 generic inverse charts; global goal open',flush=True)
if __name__=='__main__':main()
