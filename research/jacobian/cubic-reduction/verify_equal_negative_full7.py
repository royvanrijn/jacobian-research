#!/usr/bin/env python3
"""Independent stdlib/Fraction audit of the full111 structural degree-seven proof."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product
from math import gcd
import json,hashlib,time
from verify import decode,variable,add,scale,mul,compose,derivative,determinant3,evaluate
B=Path(__file__).resolve().parent;O=B/'equal-negative-full7';ROOT=B.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def total(*ps):
 z={}
 for p in ps:z=add(z,p)
 return z

def con(c,n):return {(0,)*n:Q(c)} if c else {}
def pw(p,k,n):
 z=con(1,n)
 for _ in range(k):z=mul(z,p)
 return z

def br(a,b):return add(mul(derivative(a,0),derivative(b,1)),mul(derivative(a,1),derivative(b,0)),Q(-1))
def jac(f,a,b):return total(mul(f,br(a,b)),mul(a,br(f,b)),mul(b,br(a,f)))
def coeff(p,i,j):return {e[2:]:c for e,c in p.items() if e[:2]==(i,j)}
def read(name):
 d=json.loads((O/(name+'.json')).read_text());n=len(d['variables']);return d,n,{k:decode(v,n) for k,v in d['polynomials'].items()},{k:variable(i,n) for i,k in enumerate(d['variables'])}
def lift(p,shift):
 z={}
 for e,c in p.items():
  k=e[0]+e[1]+shift;assert k>=0;z[(k,e[0],e[1])+e[2:]]=c
 return z

def main():
 started=time.monotonic();cp=O/'checkpoint.json';C=json.loads(cp.read_text());assert C['status']=='complete' and C['source_sha256']==sha(B/'equal_negative_full7.py') and C['new_solver_queries']==0
 for name,h in C['files'].items():assert h==sha(O/name)
 baseline=json.loads((O/'prior-evidence.json').read_text())
 for rel,h in baseline['immutable_old_KCUB23_files'].items():assert sha(ROOT/rel)==h,rel
 d,n,P,V=read('whole');Y,Z=V['Y'],V['Z'];f=con(1,n);A=Y;BB=Z;native=[];i=0
 for deg in range(1,4):
  for j in range(deg,-1,-1):
   f=add(f,mul(V['p'+str(i)],mul(pw(Y,j,n),pw(Z,deg-j,n))));native.append({'row':0,'exponent':[deg+1,j,deg-j],'parameter':'p'+str(i)});i+=1
 for row in [1,2]:
  for deg in range(2,5):
   for j in range(deg,-1,-1):
    term=mul(V['p'+str(i)],mul(pw(Y,j,n),pw(Z,deg-j,n)))
    if row==1:A=add(A,term)
    else:BB=add(BB,term)
    native.append({'row':row,'exponent':[deg-1,j,deg-j],'parameter':'p'+str(i)});i+=1
 assert i==33 and native==d['native_atoms']
 complete={(row,(a,b,deg-a-b)) for row,w in enumerate([1,-1,-1]) for deg in range(2,8) for a in range(deg+1) for b in range(deg-a+1) if a-b-(deg-a-b)==w}
 assert complete=={(z['row'],tuple(z['exponent'])) for z in native}
 assert [f,A,BB]==[P['f'],P['A'],P['B']] and jac(f,A,BB)==P['J']
 F=[lift(f,1),lift(A,-1),lift(BB,-1)];assert determinant3([[derivative(z,j) for j in range(3)] for z in F],n+1)==lift(P['J'],0)
 print('PASS complete33-atom source determinant and invariant dictionary',flush=True)
 d,n,P,V=read('outer');Y,Z=V['Y'],V['Z'];mons=[Y,Z,pw(Y,2,n),mul(Y,Z),pw(Z,2,n),pw(Y,3,n),mul(pw(Y,2,n),Z),mul(Y,pw(Z,2,n)),pw(Z,3,n),pw(Y,4,n)]
 A=total(*(mul(V['a'+str(j)],m) for j,m in enumerate(mons)));BB=total(*(mul(V['b'+str(j)],m) for j,m in enumerate(mons[:5]+[mons[5]])));f=total(con(1,n),pw(Y,3,n),*(mul(V['l'+str(j)],m) for j,m in enumerate(mons[:5])))
 assert [f,A,BB,jac(f,A,BB)]==[P[k] for k in ['f','A','B','J']]
 assert coeff(P['J'],0,0)==add(mul(variable(5,n-2),variable(16,n-2)),mul(variable(6,n-2),variable(15,n-2)),Q(-1))
 assert all(gcd(7,e+3)==1 for e in [1,2,3])
 # The normal spaces have strict transverse coefficient bounds for the vertical root argument.
 assert max(e[0] for e in A if e[1]>=2)<=1 and max(e[0] for e in BB if e[1]>=2)==0
 structural={}
 for name in ['graph','hyperbola1','hyperbola2']:
  d,n,P,V=read(name)
  if name=='graph':
   y,W,a,u,v,c,w,beta=[V[k] for k in d['variables']]
   a0=add(mul(add(a,mul(u,W)),y),mul(v,add(W,con(-1,n))));b0=add(mul(c,y),mul(w,add(W,con(-1,n))))
   J=mul(beta,total(mul(W,br(a0,b0)),mul(b0,derivative(a0,0)),scale(mul(a0,derivative(b0,0)),-1)))
   assert [a0,b0,J]==[P[k] for k in ['A','B','J']]
   assert total(coeff(J,0,0),scale(coeff(J,0,1),Q(1,2)),scale(coeff(J,0,2),Q(1,4)))=={}
  else:
   T,U,D,a,b,c,u,v,w=[V[k] for k in ['T','U','D','a','b','c','u','v','w']];ff=add(mul(T,U),D)
   if name=='hyperbola1':a0=total(mul(a,T),mul(b,U),mul(ff,add(mul(u,T),v)));b0=total(mul(c,T),mul(V['d'],U),mul(w,ff))
   else:a0=total(mul(a,pw(T,2,n)),mul(b,pw(U,2,n)),mul(ff,add(mul(u,T),v)));b0=add(mul(c,pw(T,2,n)),mul(w,ff))
   J=jac(ff,a0,b0);assert [ff,a0,b0,J]==[P[k] for k in ['f','A','B','J']]
   if name=='hyperbola1':
    dd=variable(0,n-2);assert total(coeff(J,0,0),scale(mul(dd,coeff(J,1,1)),Q(-1,3)),scale(mul(pw(dd,2,n-2),coeff(J,2,2)),Q(1,6)))=={}
   else:assert coeff(J,0,0)=={}
  structural[name]={'source_sha256':sha(O/(name+'.json')),'all_coefficient_identities_exact':True,'determinant_monomials':len(J),'nonzero_constant_jacobian_impossible':True}
  print('PASS small constant-Jacobian obstruction',name,flush=True)
 # Complete triangular quotient normalization, without a beta-nonzero assumption.
 d,n,P,V=read('hyperbola-coordinate');assert P['f']==add(mul(P['T'],P['U']),P['D']) and br(P['T'],P['U'])==con(1,n)
 assert coeff(P['f'],0,0)==con(1,n-2)
 # Laurent image of every basis monomial in the post-leading spaces.
 n=4;T,D,alpha,beta=[variable(j,n) for j in range(n)];yy=add(T,beta,Q(-1));qq=total(pw(yy,2,n),scale(mul(beta,yy),-1),alpha,pw(beta,2,n));uu=scale(mul(D,{(-1,0,0,0):Q(1)}),-1);zz=add(uu,qq,Q(-1));am=[yy,zz,pw(yy,2,n),mul(yy,zz),pw(zz,2,n),pw(yy,3,n),mul(pw(yy,2,n),zz),mul(yy,pw(zz,2,n)),pw(zz,3,n),pw(yy,4,n)];bm=am[:5]+[am[5]]
 assert min(e[0] for z in am for e in z)==-3 and min(e[0] for z in bm for e in z)==-2
 assert [j for j,z in enumerate(am) if max(e[0] for e in z)==6]==[8]
 assert all(min(e[0] for e in z)>=-2 for j,z in enumerate(am) if j!=8)
 assert [j for j,z in enumerate(bm) if max(e[0] for e in z)==4]==[4]
 assert all(min(e[0] for e in z)>=-1 for j,z in enumerate(bm) if j!=4)
 # Euler-Wronskian differentiation identity, used in the written all-Laurent-degree lemma.
 n=6;a,b,da,db,dda,ddb=[variable(j,n) for j in range(n)];W=add(mul(a,db),mul(b,da),Q(-1));DW=total(*(mul(derivative(W,j),z) for j,z in enumerate([da,db,dda,ddb,{},{}])));assert DW==add(mul(a,ddb),mul(b,dda),Q(-1))
 laurent_controls=0
 for k in range(1,9):
  n=5;T,a,b,c,d=[variable(j,n) for j in range(n)];plus={(k,0,0,0,0):Q(1)};minus={(-k,0,0,0,0):Q(1)};a0=add(mul(a,plus),mul(b,minus));b0=add(mul(c,plus),mul(d,minus));Dop=lambda p:mul(T,derivative(p,0));assert add(mul(a0,Dop(b0)),mul(b0,Dop(a0)),Q(-1))==scale(add(mul(b,c),mul(a,d),Q(-1)),2*k);laurent_controls+=1
 # Formal one-place-at-infinity controls, all symbolic, no numerical points.
 n=4;t,p,q,r=[variable(j,n) for j in range(n)];X=total(mul(p,pw(t,2,n)),mul(q,pw(t,4,n)),mul(r,pw(t,6,n)))
 def trunc(z):return {e:c for e,c in z.items() if e[0]<=8}
 def binomial_series(a):
  z=con(1,n);c=Q(1);power=con(1,n)
  for j in range(1,5):c*=Q(a-j+1,j);power=trunc(mul(power,X));z=add(z,scale(power,c))
  return z
 S=binomial_series(Q(1,2));Si=binomial_series(Q(-1,2));assert trunc(mul(S,S))==add(con(1,n),X) and trunc(mul(S,Si))==con(1,n)
 assert {e:c for e,c in Si.items() if e[0]==0}==con(1,n)  # -dY/(2U)=S^-1 dt, order zero
 n=5;y,u,p,q,r=[variable(j,n) for j in range(n)];cubic=total(pw(y,3,n),mul(p,pw(y,2,n)),mul(q,y),r);normal_controls=0
 for i in range(5):
  for j in range(5-i):
   normal=mul(pw(y,i,n),pw(cubic,j//2,n))
   if j%2:normal=mul(normal,u)
   assert max(2*e[0]+3*e[1] for e in normal)==2*i+3*j;normal_controls+=1
 pole_controls=0
 for m in [0]+list(range(2,13)):
  for k in [0]+list(range(2,m)):
   if k>=m:continue
   n=3;t,a,b=[variable(j,n) for j in range(n)];A=mul(a,{(-m,0,0):Q(1)});BB=mul(b,{(-k,0,0):Q(1)});ww=add(mul(A,derivative(BB,0)),mul(BB,derivative(A,0)),Q(-1));assert ww==scale(mul(mul(a,b),{(-m-k-1,0,0):Q(1)}),m-k);assert -m-k-1<0;pole_controls+=1
 # Clear-denominator completion of the square, with every cubic coefficient free.
 n=9;y,z,g,h,k,p0,p1,p2,p3=[variable(j,n) for j in range(n)];linear=add(mul(h,y),k);pc=total(p0,mul(p1,y),mul(p2,pw(y,2,n)),mul(p3,pw(y,3,n)));ff=total(mul(g,pw(z,2,n)),mul(linear,z),pc)
 assert scale(mul(g,ff),4)==total(pw(add(scale(mul(g,z),2),linear),2,n),scale(mul(g,pc),4),scale(pw(linear,2,n),-1))
 # Generic affine-in-Z determinant from independently prescribed derivative jets.
 n=11;z,f,fp,a,ap,b,bp,c,ccp,d,dp=[variable(j,n) for j in range(n)];aa=add(mul(a,z),b);bb=add(mul(c,z),d);ay=add(mul(ap,z),bp);by=add(mul(ccp,z),dp)
 actual=total(mul(f,add(mul(ay,c),mul(a,by),Q(-1))),mul(fp,add(mul(aa,c),mul(bb,a),Q(-1))))
 jz=mul(f,add(mul(ap,c),mul(a,ccp),Q(-1)));j0=total(mul(f,add(mul(bp,c),mul(a,dp),Q(-1))),mul(fp,add(mul(b,c),mul(a,d),Q(-1))))
 assert actual==add(mul(z,jz),j0)
 assert add(j0,mul(a,add(mul(f,dp),mul(fp,d))))==mul(c,add(mul(f,bp),mul(fp,b)))
 # Matching old three-variable upper, reconstructed from retained coefficients.
 n=3;raw=B/'raw13-specialization/input.json';D=json.loads(raw.read_text());F=[variable(j,n) for j in range(n)]
 for atom in D['atoms']:F[atom['row']]=add(F[atom['row']],{tuple(atom['exponent']):Q(atom['coefficient'])})
 assert determinant3([[derivative(z,j) for j in range(n)] for z in F],n)==con(1,n)
 assert [sum(sum(e)>=2 for e in z) for z in F]==[2,5,6] and max(sum(e) for z in F for e in z)==7
 assert all(sum(a*b for a,b in zip(e,[1,-1,-2]))==w for z,w in zip(F,[1,-1,-2]) for e in z)
 points=[[0,1,-8],[-1,1,4]];assert [[evaluate(z,p) for z in F] for p in points]==[[Q(0),Q(1),Q(0)]]*2
 ledger=json.loads((ROOT/'MATH_STATUS.json').read_text());E={e['id']:e for e in ledger['entries']};pins={}
 for k in ['KCUB23','KCUB39','KCUB38','KCUB36']:
  e=E[k];assert e['state']=='proved';pins[k]={'scope_fingerprint':hashlib.sha256((e['state']+'\0'+e['title']+'\0'+e['scope']).encode()).hexdigest()[:16],'proof_sha256':sha(ROOT/e['canonical_source'])}
 receipt={'status':'pass','claim':'KCUB40','checker_sha256':sha(Path(__file__)),'producer_sha256':sha(B/'equal_negative_full7.py'),'arithmetic_sha256':sha(B/'verify.py'),'checkpoint_sha256':sha(cp),'model_hashes':C['files'],'complete_native_atoms':33,'post_leading_parameters':21,'full_source_dictionary_checked':True,'prime_degree_pairs':[[7,e+3] for e in [1,2,3]],'small_coefficient_obstructions':structural,'laurent_euler_jet_identity_checked':True,'generic_opposite_power_controls':laurent_controls,'elliptic_normal_monomial_controls':normal_controls,'elliptic_distinct_pole_controls':pole_controls,'elliptic_differential_has_order_zero':True,'vertical_affine_determinant_identity_checked':True,'all_degree_curve_and_laurent_lemmas_are_written_not_formalized':True,'old_evidence_sha256':sha(O/'prior-evidence.json'),'old40second_timeouts_preserved':True,'large_old_primary_certificate_unreplayed_and_unused':True,'new_solver_queries':0,'imports':pins,'remaining_degree7_diagonal_weight_raw12_types':0,'only_degree7_diagonal_weight_collision_type':[1,-1,-2],'restricted_minimum_degree':7,'restricted_minimum_nonlinear_occurrences_through_degree7':13,'restricted_minimum_collision_dimension_all_degrees':3,'matching_old_upper_control':{'input_sha256':sha(raw),'degree':7,'counts':[2,5,6],'full_determinant_one':True,'collision_points':points,'common_image':[0,1,0],'new_witness':False},'goal_complete':False,'no_global_positive_bound_or_whole_profile_change':True,'trust_boundary':'Universal leading-form coverage,pole and Laurent coprimality/Euler-support proofs,k<=2 exhaustion,lift divisibility and vertical root/unit argument are written,not formalized or externally reviewed. Independent exact finite algebraic premises and old upper are checked. KCUB23 imported for lower/strict branches and constant first quotient;KCUB39/38/36 for complete weighted frontier and rank-one low-dimensional exclusion. Raw count requires displayed diagonal source coordinates;linear diagonalization need not preserve it.','seconds':time.monotonic()-started}
 proof=B/'EQUAL_NEGATIVE_FULL_DEGREE7.md'
 if proof.exists():receipt['proof_sha256']=sha(proof)
 (O/'validation.json').write_text(json.dumps(receipt,indent=2)+'\n');print('PASS full111 degree7 structural closure;33 native atoms,three coefficient identities,all exact curve controls;full displayed diagonal-weight raw13 sharp;global goal open',flush=True)
if __name__=='__main__':main()
