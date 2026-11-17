#!/usr/bin/env python3
"""Independent sparse Fraction verification of two changed companion pencils.
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
 assert labels==d['row_labels'];numeric=[[row.get(k,{}).get((0,),Q(0)) for k in range(len(mons))] for row in rows];zero_rank=rank(numeric,len(mons));assert zero_rank==d['rank_at_zero_primary']==35;used=set()
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
 assert len(d['pivots'])==d['rank']==37 and all(not row for row in rows[37:]);return dict(certificate_sha256=sha(certpath),rank_nonzero_parameter=37,rank_at_zero=zero_rank,unit_pivots=37,full_residual_zero=True)

def cube_replay(path,inp,H0,H1):
 c=load(path);sigma=Q(c['sigma']);assert sigma;forms=c['forms'];B=[[Q(v) for v in row] for row in c['B']];m=len(forms);assert c['input_sha256']==sha(inp) and c['size']==m==98 and len(B)==19 and all(len(row)==m for row in B);assert forms[c['padding_index']]==[0]*18+[1] and B[18]==[Q(i==c['padding_index']) for i in range(m)];result=[{} for _ in range(19)]
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
 started=time.monotonic();base=HERE/'cubic12-integer-witness.json';bd=load(base);K=[decode(row,12) for row in bd['coefficients']];validation=load(HERE/'integer-validation12.json');assert validation['status']=='pass' and validation['witness_sha256']==sha(base) and validation['determinant_one'];names=['x'+str(i) for i in range(1,20)]+['sigma'];results=[]
 for tag,axis,producer,search in [('companion-pencil19','s','prepare_companion_pencil19.py','search_companion_pencil19.py'),('companion-w19','w','prepare_companion_w19.py','search_companion_w19.py')]:
  folder=HERE/tag;ip=folder/'input.json';d=load(ip);assert d['source_sha256']==sha(HERE/producer) and d['base_sha256']==sha(base);E=[variable(i,20) for i in range(12)];sig=variable(19,20)
  if axis=='s':E[4]=add(E[4],scale(mul(sig,mul(E[1],E[2])),Q(1,2)));E[10]=add(E[10],mul(sig,mul(E[1],E[3])));active=[0,1,2,3,4,8]
  else:E[7]=add(E[7],mul(sig,mul(E[2],E[3])));E[11]=add(E[11],mul(sig,power(E[3],2,20)));active=[0,1,2,3,4,9]
  assert E==[expression(text,names) for text in d['source_shear']];G=[sub(f,E,20) for f in K];assert G==[expression(text,names) for text in d['base_family']];Q2=[{e:c for e,c in f.items() if sum(e[:19])==2} for f in G];C3=[{e:c for e,c in f.items() if sum(e[:19])==3} for f in G];assert all(G[i]==add(add(variable(i,20),Q2[i]),C3[i]) for i in range(12)) and all(not f for i,f in enumerate(C3) if i not in active);H=[mul(variable(18,20),f) for f in Q2]
  for j,i in enumerate(active):H[i]=add(H[i],mul(power(variable(18,20),2,20),variable(12+j,20)),Q(-1))
  H += [C3[i] for i in active]+[power(variable(18,20),3,20)];H0=[dec(p) for p in d['H0']];H1=[dec(p) for p in d['H1']];expected=[{e+(0,):c for e,c in p.items()} for p in H0]
  for i,p in enumerate(H1):expected[i]=add(expected[i],{e+(1,):c for e,c in p.items()})
  assert H==expected and all(sum(e)==3 for p in H0+H1 for e in p);rankproof=rank_replay(folder,H0,H1);candidate,sv,target=cube_replay(folder/'best-identity.json',ip,H0,H1);assert candidate['source_sha256']==sha(HERE/search)
  # Six independent cubic output rows for every sigma!=0.
  words=['x1**2*x2','x1*x4*x8','x1*x4*x5','x1**2*x3','x2**3','x1*x2*x3' if axis=='s' else 'x1*x3*x4'];mons=[next(iter(expression(text,names[:19]))) for text in words];matrix=[[{(e[19],):c for e,c in C3[i].items() if e[:19]==mon} for mon in mons] for i in active];assert determinant(matrix)=={(1,):Q(-288)}
  points=[[Q(v) for v in point] for point in bd['points']];oldimages=[[evaluate(f,p) for f in K] for p in points];assert oldimages[0]==oldimages[1] and points[0]!=points[1]
  for p in points:
   if axis=='s':p[4]-=sv*p[1]*p[2]/2;p[10]-=sv*p[1]*p[3]
   else:p[7]-=sv*p[2]*p[3];p[11]-=sv*p[3]**2
  newimages=[[evaluate(f,p+[Q(0)]*7+[sv]) for f in G] for p in points];assert newimages==oldimages
  lifted=[]
  for p in points:
   cs=[evaluate(C3[i],p+[Q(0)]*7+[sv]) for i in active];lifted.append(p+[-c for c in cs]+[Q(1)])
  companion=[add(variable(i,19),poly) for i,poly in enumerate(target)];companion[18]=add(companion[18],power(variable(18,19),3,19),Q(-1));images=[[evaluate(f,p) for f in companion] for p in lifted];assert lifted[0]!=lifted[1] and images[0]==images[1]==oldimages[0]+[Q(0)]*6+[Q(1)]
  uniformfolder=HERE/'companion-pencil19-mixed' if axis=='s' else folder;uniform=uniform_replay(uniformfolder,H0,H1);assert uniform['uniform_cube_count']==(100 if axis=='s' else 101)
  results.append(dict(axis=axis,input_sha256=sha(ip),source_shear_two_commuting_triangular_quadratics=True,full_base_composition_checked=True,companion_reconstructed=True,cubic_output_rank_nonzero_parameter=6,derivative_rank=rankproof,shared_cube_lower_bound_nonzero_parameter=37,cube_identity_sha256=sha(folder/'best-identity.json'),cube_identity_size=98,identity_parameter=str(sv),cube_identity_all_coefficients_checked=True,coefficient_matrix_rank=19,padding_normalization_checked=True,uniform_identity=uniform,base_collision_points=[[str(v) for v in p] for p in points],companion_collision_points=[[str(v) for v in p] for p in lifted],companion_common_image=[str(v) for v in images[0]],homogeneous_companion_collision_checked=True))
  print('VERIFIED',axis,'rank37','cube98 at',sv,'uniform',uniform['uniform_cube_count'],flush=True)
 receipt=dict(status='pass',source_sha256=sha(Path(__file__)),arithmetic_sha256=sha(HERE/'verify.py'),utility_sha256=sha(HERE/'weighted_verify_utils.py'),base_witness_sha256=sha(base),base_validation_sha256=sha(HERE/'integer-validation12.json'),pencils=results,homogeneous_companion_dimension=19,nonzero_cubic_rank=6,shared_cube_lower_bound=37,best_verified_cube_identity=98,all_polynomial_and_rank_identities_replayed=True,polynomial_substitution_domain_guarded=True,determinant_basis='Pinned KCUB3 determinant, unit triangular source shears, and written standard companion Schur factorization for all lambda.',full_19_variable_determinant_expanded=False,global_bound_improved=False,global_cubic_linear_upper_bound=89,scope='Two specified changed rank6 companion pencils. Every nonzero member has derivative-space rank37, hence needs at least37 cubes. Each pencil has a verified98-cube identity at its recorded parameter and a uniform100/101-cube identity. These are conditional tensor bounds, not a new global size record or a proof that the finite pools are optimal. The base/source and full19-variable companion collisions are checked; no dominated98-variable cubic-linear matrix is advertised as an improvement.',seconds=round(time.monotonic()-started,3));(HERE/'companion-pencils19-validation.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt),flush=True)
if __name__=='__main__':main()
