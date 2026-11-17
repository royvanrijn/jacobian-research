#!/usr/bin/env python3
"""Independent sparse Fraction verification of the fixed-z coupled companion pencil.
No producer polynomial, rank or cube routines and no CAS/SMT are imported.
The base determinant is a pinned import; companion block elimination is written proof.
"""
from pathlib import Path
from fractions import Fraction as Q
from itertools import permutations
import hashlib,json,time
from verify import add,mul,variable,derivative,decode,evaluate
from verify_power_single_compress import scale
from weighted_verify_utils import expression,substitution,power,rank,coefficient_buckets
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def const(v,n):return {(0,)*n:Q(v)} if v else {}
def sub(poly,values,n):
 assert all(all(k>=0 for k in e) for e in poly)
 return substitution(poly,values,n)
def dec(p):return {tuple(e):Q(c) for e,c in p}
def determinant(M):
 n=len(M);out={}
 for perm in permutations(range(n)):
  term={(0,):Q((-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n)))}
  for i,j in enumerate(perm):term=mul(term,M[i][j])
  out=add(out,term)
 return out

def rank_replay(folder,H0,H1):
 certpath=folder/'derivative-rank-certificate.json';d=load(certpath);assert d['source_sha256']==sha(HERE/'companion_derivative_rank.py') and d['input_sha256']==sha(folder/'input.json') and d['status']=='exact-unit-rank-certificate' and d['excluded_parameter_values']==['0'];D=[(i,j,derivative(H0[i],j),derivative(H1[i],j)) for i in range(19) for j in range(19)];mons=sorted(set().union(*(set(p) for i,j,a,b in D for p in [a,b])));assert [list(e) for e in mons]==d['monomials'];rows=[];labels=[]
 for i,j,a,b in D:
  row={}
  for k,e in enumerate(mons):
   p={}
   if a.get(e):p[(0,)]=a[e]
   if b.get(e):p[(1,)]=b[e]
   if p:row[k]=p
  if row:rows.append(row);labels.append([i,j])
 assert labels==d['row_labels'];numeric=[[row.get(k,{}).get((0,),Q(0)) for k in range(len(mons))] for row in rows];zero_rank=rank(numeric,len(mons));assert zero_rank==d['rank_at_zero_primary']==37;used=set()
 for step,pivot in enumerate(d['pivots']):
  i,k=pivot['row'],pivot['column'];assert step<=i<len(rows) and k not in used;rows[step],rows[i]=rows[i],rows[step];p=rows[step][k];assert p=={(pivot['power'],):Q(pivot['coefficient'])};inv={(-pivot['power'],):1/Q(pivot['coefficient'])};rows[step]={col:mul(poly,inv) for col,poly in rows[step].items()};assert rows[step][k]==const(1,1)
  for i in range(step+1,len(rows)):
   factor=rows[i].get(k)
   if factor:
    for col,poly in rows[step].items():
     value=add(rows[i].get(col,{}),mul(factor,poly),Q(-1))
     if value:rows[i][col]=value
     else:rows[i].pop(col,None)
  used.add(k)
 assert len(d['pivots'])==d['rank']==39 and all(not row for row in rows[39:]);return dict(certificate_sha256=sha(certpath),rank_nonzero_parameter=39,rank_at_zero=zero_rank,unit_pivots=39,full_residual_zero=True)

def cube_replay(path,inp,H0,H1):
 c=load(path);sigma=Q(c['parameter_value']);forms=c['forms'];B=[[Q(v) for v in row] for row in c['B']];m=len(forms);assert c['input_sha256']==sha(inp) and c['size']==m==101 and len(B)==19 and all(len(row)==m for row in B);assert forms[c['padding_index']]==[0]*18+[1] and B[18]==[Q(i==c['padding_index']) for i in range(m)];result=[{} for _ in range(19)]
 for j,v in enumerate(forms):
  assert len(v)==19 and any(v);ell={}
  for i,value in enumerate(v):
   if value:ell=add(ell,scale(variable(i,19),Q(value)))
  cube=power(ell,3,19)
  for i in range(19):
   if B[i][j]:result[i]=add(result[i],scale(cube,B[i][j]))
 target=[add(a,b,sigma) for a,b in zip(H0,H1)];assert result==target;assert rank(B,m)==19;return c,sigma,target

def uniform_replay(folder,H0,H1):
 d=load(folder/'input.json');b=load(folder/'basis.json');assert b['input_sha256']==sha(folder/'input.json') and [dec(p) for p in d['H0']]==H0 and [dec(p) for p in d['H1']]==H1;result=[[{} for _ in range(19)] for j in range(2)];count=0
 for form_index,row in zip(b['basis'],b['B']):
  if not row:continue
  count+=1;ell={}
  for j,v in enumerate(d['forms'][form_index]):
   if v:ell=add(ell,scale(variable(j,19),Q(v)))
  cube=power(ell,3,19)
  for key,value in row.items():
   col=int(key);part,j=divmod(col,19);result[part][j]=add(result[part][j],scale(cube,Q(value)))
 assert result==[H0,H1];return dict(input_sha256=sha(folder/'input.json'),basis_sha256=sha(folder/'basis.json'),uniform_cube_count=count)

def main():
 started=time.monotonic();folder=HERE/'coupled-rankone19';ip=folder/'input.json';d=load(ip);base=HERE/'cubic12-integer-witness.json';bd=load(base);K=[decode(row,12) for row in bd['coefficients']];bv=load(HERE/'integer-validation12.json');assert bv['status']=='pass' and bv['determinant_one'] and bv['witness_sha256']==sha(base);assert d['base_sha256']==sha(base) and d['source_sha256']==sha(HERE/'prepare_coupled_rankone19.py') and d['fixed_z']=='-4/3' and not d['excluded_parameter_values'];names=['x'+str(i) for i in range(1,20)]+['w'];X=[variable(i,20) for i in range(20)];a,b,c,R=X[:4];w=X[19];z=Q(-4,3);E=list(X[:12]);E[7]=add(E[7],mul(w,mul(c,R)));E[11]=add(E[11],mul(w,power(R,2,20)));E[8]=add(E[8],scale(mul(c,R),z));E[10]=add(E[10],scale(mul(R,add(scale(b,3),X[6],Q(-1))),z));assert E==[expression(t,names) for t in d['source_shear']];G=[sub(f,E,20) for f in K];assert G==[expression(t,names) for t in d['base_family']];Q2=[{e:v for e,v in f.items() if sum(e[:19])==2} for f in G];C3=[{e:v for e,v in f.items() if sum(e[:19])==3} for f in G];assert all(G[i]==add(add(X[i],Q2[i]),C3[i]) for i in range(12));Cstar=C3[:5]+[mul(a,mul(c,R))];assert Cstar==[expression(t,names) for t in d['cubic_star']];iota=[[expression(t,names) for t in row] for row in d['embedding']];expected=[[const(int(i==j and j<5),20) for j in range(6)] for i in range(12)];expected[9][5]=w;expected[10][5]=const(z,20);assert iota==expected
 for i in range(12):
  v={}
  for j in range(6):v=add(v,mul(iota[i][j],Cstar[j]))
  assert v==C3[i]
 H=[mul(X[18],f) for f in Q2]
 for i in range(12):
  for j in range(6):H[i]=add(H[i],mul(power(X[18],2,20),mul(iota[i][j],X[12+j])),Q(-1))
 H+=Cstar+[power(X[18],3,20)];H0=[dec(p) for p in d['H0']];H1=[dec(p) for p in d['H1']];combined=[add({e+(0,):v for e,v in p.items()},{e+(1,):v for e,v in q.items()}) for p,q in zip(H0,H1)];assert H==combined and all(sum(e)==3 for f in H0+H1 for e in f)
 # Constant minor gives rank6 for every w, including zero; iota injective since z!=0.
 words=['x1**2*x2','x1*x4*x8','x1*x4*x5','x1**2*x3','x2**3','x1*x3*x4'];mons=[next(iter(expression(t,names[:19]))) for t in words];M=[[{(e[19],):v for e,v in f.items() if e[:19]==mon} for mon in mons] for f in Cstar];minor=determinant(M);assert minor=={(0,):Q(-288)},minor
 rankproof=rank_replay(folder,H0,H1);candidate,w0,target=cube_replay(folder/'best-identity.json',ip,H0,H1);assert candidate['source_sha256']==sha(HERE/'search_coupled_rankone19.py');uniform=uniform_replay(folder,H0,H1)
 points=[[Q(v) for v in p] for p in bd['points']];oldimages=[[evaluate(f,p) for f in K] for p in points];assert oldimages[0]==oldimages[1] and points[0]!=points[1]
 for p in points:p[7]-=w0*p[2]*p[3];p[11]-=w0*p[3]**2;p[8]-=z*p[2]*p[3];p[10]-=z*p[3]*(3*p[1]-p[6])
 assert [[evaluate(f,p+[Q(0)]*7+[w0]) for f in G] for p in points]==oldimages;lifted=[p+[-evaluate(f,p+[Q(0)]*7+[w0]) for f in Cstar]+[Q(1)] for p in points];companion=[add(variable(i,19),f) for i,f in enumerate(target)];companion[18]=add(companion[18],power(variable(18,19),3,19),Q(-1));images=[[evaluate(f,p) for f in companion] for p in lifted];assert lifted[0]!=lifted[1] and images[0]==images[1]==oldimages[0]+[Q(0)]*6+[Q(1)]
 cp=load(folder/'search-checkpoint.json');assert cp['source_sha256']==sha(HERE/'search_coupled_rankone19.py') and cp['input_sha256']==sha(ip) and cp['best_identity_sha256']==sha(folder/'best-identity.json') and cp['best']==101 and cp['candidate_forms']==478
 # Reconstruct grouped polynomial expressions from every recorded rank-one factor.
 grouped=[{} for _ in range(19)]
 for term in d['grouped_terms']:
  f=const(1,20)
  for values in term['factors']:
   ell={}
   for j,v in enumerate(values):ell=add(ell,scale(X[j],Q(v)))
   f=mul(f,ell)
  co=add(const(Q(term['coefficient0']),20),scale(w,Q(term['coefficient1'])));i=term['row'];grouped[i]=add(grouped[i],mul(co,f))
 assert grouped==H
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),input_sha256=sha(ip),base_witness_sha256=sha(base),base_validation_sha256=sha(HERE/'integer-validation12.json'),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),fixed_z=str(z),all_w_valid=True,full_source_composition_checked=True,cubic_degree_checked=True,cubic_rank_all_parameters=6,constant_minor=str(minor[(0,)]),companion_dimension=19,embedding_checked=True,full_companion_coefficients_checked=True,grouped_rankone_coefficients_checked=True,derivative_rank=rankproof,shared_cube_lower_bound_all_parameters=37,shared_cube_lower_bound_nonzero_w=39,cube_identity_size=101,cube_identity_parameter=str(w0),cube_identity_sha256=sha(folder/'best-identity.json'),all_cube_coefficients_checked=True,coefficient_matrix_rank=19,pure_tau_padding_checked=True,uniform_identity=uniform,base_collision_points=[[str(v) for v in p] for p in points],companion_collision_points=[[str(v) for v in p] for p in lifted],companion_common_image=[str(v) for v in images[0]],full_companion_collision_checked=True,search_steps=cp['steps'],search_seconds=cp['seconds'],global_bound_improved=False,global_cubic_linear_upper_bound=89,determinant_basis='Pinned KCUB3 determinant, four commuting unit triangular shears, and written standard companion Schur factorization. No19-variable determinant expansion.',scope='Fixed-z=-4/3 coupled companion pencil, every w including0. Exact derivative ranks39 offzero and37 atzero; verified101-cube identity atw0 and saved uniform identity. Bounded478-form basis search has no optimality consequence. No global dimension/support lower bound or improved witness.',seconds=round(time.monotonic()-started,3));(folder/'validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print('VERIFIED coupled: derivative39/37; cube101 atw',w0,'uniform',uniform['uniform_cube_count'],'seconds',receipt['seconds'])
if __name__=='__main__':main()
