"""Fail-closed verification of F7 cubic-cube predecessor identities and a full coefficient-space determinant.

Does not run Julia, F4, CRT or rational reconstruction. Uses finite modular arithmetic only for the conclusive full coefficient-space determinant.
A partial reconstruction is deliberately insufficient.
"""
from pathlib import Path
import collections,gzip,hashlib,json,platform,time,math
import sympy as S
import ijson,flint
from flint import fmpq,fmpz,fmpq_mpoly_ctx,fmpz_mpoly_ctx,nmod_mat
if not __debug__:raise RuntimeError('Run without -O')
HERE=Path(__file__).resolve().parent
START=time.monotonic();checks=[]
def yes(name,condition):
 assert condition,name
 checks.append(name)
def zero(name,expr):yes(name,S.cancel(expr)==0)

t,Y1,Y2,Y3,Y4=S.symbols('t Y1 Y2 Y3 Y4')
u,v,w,z,h,k=S.symbols('u v w z h k')
V=[u,v,w,z,h]
b=Y4-3*Y2*t*t+2*Y1*t**3-t**4
a=Y3-2*Y2*t+Y1*t*t-t**3-t*t*b
R=S.expand(a*a+a+t*b+t**3+Y1*t-Y2)
P=S.Poly(S.expand((t**3+u*t*t+v*t+w)**3*(t**3+z*t*t+h*t+k)),t)
y1=-P.nth(11)/4
y2=(P.nth(10)-4*y1*y1)/6
kfit=-12*y1*y2-2-P.nth(9).subs(k,0)
P=S.Poly(S.expand(P.as_expr().subs(k,kfit)),t)
y4=(6*y1+9*y2*y2-P.nth(8))/2
y3=(P.nth(6)-14*y1*y2+6*y2*y4-2)/2
res=S.Poly(S.expand(R.subs({Y1:y1,Y2:y2,Y3:y3,Y4:y4})-P.as_expr()),t)
for j in [12,11,10,9,8,6]:zero('eliminated coefficient '+str(j),res.nth(j))
source=json.loads((HERE/'cubic-six-input.json').read_text())
yes('source variable order',source['variables']==list(map(str,V)))
yes('six specified necessary coefficients',[row['coefficient'] for row in source['equations']]==[7,3,5,4,1,2])
for row in source['equations']:
 ratio=S.cancel(res.nth(row['coefficient'])/S.sympify(row['polynomial']))
 yes('independently rederived coefficient '+str(row['coefficient']),ratio.is_Rational and ratio!=0)

# Any claimed proof must include the terminal identity. The partial file is not read.
certificate=HERE/'cubic-six-rank-dag.json.gz'
yes('complete certificate is present',certificate.is_file())
def certificate_field(key):
 with gzip.open(certificate,'rb') as fd:return next(ijson.items(fd,key))
cert={key:certificate_field(key) for key in ['variables','source','terminal']}
def certificate_nodes():
 with gzip.open(certificate,'rb') as fd:yield from ijson.items(fd,'nodes.item')
def certificate_hash():
 with certificate.open('rb') as fd:return hashlib.file_digest(fd,'sha256').hexdigest()
yes('certificate variable and generator source match',cert['variables']==source['variables'] and cert['source']=='cubic-six-input.json')
ctx=fmpz_mpoly_ctx.get(source['variables']);originals=[]
for row in source['equations']:
 originals.append(ctx.from_dict({ex:fmpz(int(c)) for ex,c in S.Poly(S.sympify(row['polynomial']),*V).terms()}))

def integral_polynomial(rows):
 rational={};denominator=fmpz(1)
 for term in rows:
  ex=tuple(term['exponents'])
  assert len(ex)==5 and all(type(x) is int and x>=0 for x in ex) and ex not in rational
  c=fmpq(term['coefficient']);assert c!=0
  rational[ex]=c;denominator=denominator.lcm(c.q)
 assert denominator>0
 integer={}
 for ex,c in rational.items():
  assert denominator%c.q==0
  integer[ex]=c.p*(denominator//c.q)
 return ctx.from_dict(integer),denominator

polys={};denominators={};local_identities=0;source_nodes=0;relation_terms=0;polynomial_terms=0
for row in certificate_nodes():
 i=row['id'];assert type(i) is int and i>0 and i not in polys
 poly,denominator=integral_polynomial(row['polynomial']);polynomial_terms+=len(row['polynomial'])
 if 'source_generator' in row:
  index=row['source_generator'];assert type(index) is int and 0<=index<len(originals)
  scalar=fmpq(row['source_scalar']);assert scalar!=0 and 'relations' not in row
  scaled=scalar*denominator
  yes('integer-cleared source identity '+str(i),scaled.q*poly==scaled.p*originals[index]);source_nodes+=1
 else:
  scaled_relations=[];common=denominator;seen=set()
  for rel in row['relations']:
   parent=rel['parent'];assert type(parent) is int and 0<parent<i and parent in polys
   ex=tuple(rel['exponents']);assert len(ex)==5 and all(type(x) is int and x>=0 for x in ex)
   assert (parent,ex) not in seen;seen.add((parent,ex))
   c=fmpq(rel['coefficient']);assert c!=0
   scaled=c/denominators[parent];common=common.lcm(scaled.q)
   scaled_relations.append((parent,ex,scaled));relation_terms+=1
  assert common>0 and common%denominator==0
  total=-(common//denominator)*poly
  grouped=collections.defaultdict(dict)
  for parent,ex,c in scaled_relations:
   assert common%c.q==0
   grouped[parent][ex]=c.p*(common//c.q)
  for parent,multiplier in grouped.items():total+=ctx.from_dict(multiplier)*polys[parent]
  yes('literal integer-cleared local identity '+str(i),total==0);local_identities+=1
 polys[i]=poly;denominators[i]=denominator
 if len(polys)%50==0:
  progress={'stage':'integer-cleared predecessor replay','nodes_checked':len(polys),'seconds':time.monotonic()-START,'global_exclusion':'PENDING'}
  (HERE/'rank-integer-replay-progress.json').write_text(json.dumps(progress)+'\n')
  print(json.dumps(progress),flush=True)
# This is a full coefficient-space determinant, not a modular unit-ideal report.
terminal=cert['terminal']
yes('terminal is a full integer coefficient matrix',terminal['kind']=='full_integer_coefficient_matrix')
prime=terminal['prime']
yes('fixed modulus',prime==1000000007)
yes('modulus primality by complete trial division',all(prime%d for d in range(2,math.isqrt(prime)+1)))
integer_parents={};columns=[];monomials=set()
for column in terminal['columns']:
 parent=column['parent'];ex=tuple(column['exponents'])
 assert type(parent) is int and parent in polys
 assert len(ex)==5 and all(type(e) is int and e>=0 for e in ex)
 if parent not in integer_parents:integer_parents[parent]=polys[parent].to_dict()
 terms={tuple(int(a)+b for a,b in zip(e,ex)):c for e,c in integer_parents[parent].items()}
 columns.append(terms);monomials.update(terms)
# Include every actual characteristic-zero coefficient row, even if its entries
# all reduce to zero modulo the chosen prime. No row may be silently discarded.
monomials=sorted(monomials)
yes('entire coefficient support is square',len(monomials)==len(columns)>0)
yes('retained shape agrees with all actual rows',[len(monomials),len(columns)]==terminal['expected_shape']==[402,402])
yes('constant monomial belongs to the coefficient space',(0,0,0,0,0) in monomials)
indices={e:i for i,e in enumerate(monomials)}
matrix=nmod_mat(len(monomials),len(columns),prime)
for j,col in enumerate(columns):
 for e,c in col.items():matrix[indices[e],j]=int(c%prime)
determinant=int(matrix.det())
yes('nonzero complete coefficient determinant',determinant!=0 and determinant==terminal['expected_det_mod'])
# Recompute the same finite determinant with plain Python integer arithmetic.
plain=[[int(matrix[i,j]) for j in range(matrix.ncols())] for i in range(matrix.nrows())]
other=1
for k in range(len(plain)):
 pivot=next((i for i in range(k,len(plain)) if plain[i][k]),None)
 assert pivot is not None
 if pivot!=k:plain[k],plain[pivot]=plain[pivot],plain[k];other=-other
 value=plain[k][k];other=other*value%prime;inverse=pow(value,-1,prime);pivotrow=plain[k]
 for i in range(k+1,len(plain)):
  row=plain[i]
  if not row[k]:continue
  multiplier=row[k]*inverse%prime;row[k]=0
  for j in range(k+1,len(plain)):row[j]=(row[j]-multiplier*pivotrow[j])%prime
yes('independent Python determinant replay',other%prime==determinant)
# Three-point types divisible by3 are cube families even when Q3 repeats.
for multiplicities in [[9],[6,3],[3,3,3]]:
 yes('covered nonsimple type '+str(multiplicities),sum(multiplicities)==9 and all(m%3==0 for m in multiplicities))
receipt={'result':'PASS','claim_id':'GAOF7CUBE1','checks':len(checks),'python':platform.python_version(),'sympy':S.__version__,'python_flint':flint.__version__,'ijson':ijson.__version__,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((HERE/'cubic-six-input.json').read_bytes()).hexdigest(),'certificate_sha256':certificate_hash(),'source_nodes':source_nodes,'terminal_matrix_shape':[len(monomials),len(columns)],'terminal_prime':prime,'terminal_det_mod':determinant,'terminal_det_engines':['FLINT','plain Python integers'],'terminal_argument':'Nonzero modular determinant of the entire square INTEGER coefficient matrix implies nonzero integer determinant, hence full coefficient-space span over Q; the constant monomial is present.','local_rational_identities':local_identities,'polynomial_terms':polynomial_terms,'relation_terms':relation_terms,'excluded_nonsimple_types':[[9],[6,3],[3,3,3]],'remaining_nonsimple_types':[[7,2],[5,4],[5,2,2],[4,3,2],[3,2,2,2]],'three_point_fibres':'OPEN; the cubic-cube family is excluded, remaining linear-cube/cubic-square cases unresolved','seconds':time.monotonic()-START,'assurance':'Independent source reconstruction, literal FLINT replay of every predecessor identity after exact integer denominator clearing, and two exact finite determinant replays covering every coefficient row; no independent end-to-end or formal proof review, external review or novelty priority.'}
receipt['checker_mode']='integer-cleared predecessor identities plus full coefficient-space determinant'
(HERE/'receipt-integer.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
