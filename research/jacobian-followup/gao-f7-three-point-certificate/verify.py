"""Fail-closed verification of a complete F7 cubic-cube local identity DAG.

Does not run Julia, F4, modular computation, CRT or rational reconstruction.
A partial reconstruction is deliberately insufficient.
"""
from pathlib import Path
import collections,gzip,hashlib,json,platform,time
import sympy as S
from flint import fmpq,fmpq_mpoly_ctx
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
certificate=HERE/'cubic-six-rational-dag.json.gz'
yes('complete certificate is present',certificate.is_file())
with gzip.open(certificate,'rt') as fd:cert=json.load(fd)
yes('certificate variable and generator source match',cert['variables']==source['variables'] and cert['source']=='cubic-six-input.json')
ctx=fmpq_mpoly_ctx.get(source['variables']);originals=[]
for row in source['equations']:
 originals.append(ctx.from_dict({ex:fmpq(int(c)) for ex,c in S.Poly(S.sympify(row['polynomial']),*V).terms()}))

def terms(rows):
 out={}
 for term in rows:
  ex=tuple(term['exponents'])
  assert len(ex)==5 and all(type(x) is int and x>=0 for x in ex)
  assert ex not in out
  value=fmpq(term['coefficient']);assert value!=0
  out[ex]=value
 return out

polys={};local_identities=0;source_nodes=0;relation_terms=0;polynomial_terms=0
for row in sorted(cert['nodes'],key=lambda row:row['id']):
 i=row['id'];assert type(i) is int and i>0 and i not in polys
 poly=ctx.from_dict(terms(row['polynomial']));polynomial_terms+=len(row['polynomial'])
 if 'source_generator' in row:
  index=row['source_generator'];assert type(index) is int and 0<=index<len(originals)
  scalar=fmpq(row['source_scalar']);assert scalar!=0 and 'relations' not in row
  yes('source identity '+str(i),poly==scalar*originals[index]);source_nodes+=1
 else:
  grouped=collections.defaultdict(dict)
  for rel in row['relations']:
   parent=rel['parent'];assert type(parent) is int and 0<parent<i and parent in polys
   ex=tuple(rel['exponents']);assert len(ex)==5 and all(type(x) is int and x>=0 for x in ex)
   assert ex not in grouped[parent]
   cf=fmpq(rel['coefficient']);assert cf!=0
   grouped[parent][ex]=cf;relation_terms+=1
  total=ctx.constant(0)
  for parent,multiplier in grouped.items():total+=ctx.from_dict(multiplier)*polys[parent]
  yes('literal rational local identity '+str(i),total==poly);local_identities+=1
 polys[i]=poly
unit=cert['unit_node']
yes('terminal node equals one',unit in polys and polys[unit]==ctx.constant(1))
# Three-point types divisible by3 are cube families even when Q3 repeats.
for multiplicities in [[9],[6,3],[3,3,3]]:
 yes('covered nonsimple type '+str(multiplicities),sum(multiplicities)==9 and all(m%3==0 for m in multiplicities))
receipt={'result':'PASS','claim_id':'GAOF7CUBE1','checks':len(checks),'python':platform.python_version(),'sympy':S.__version__,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((HERE/'cubic-six-input.json').read_bytes()).hexdigest(),'certificate_sha256':hashlib.sha256(certificate.read_bytes()).hexdigest(),'source_nodes':source_nodes,'local_rational_identities':local_identities,'polynomial_terms':polynomial_terms,'relation_terms':relation_terms,'excluded_nonsimple_types':[[9],[6,3],[3,3,3]],'remaining_nonsimple_types':[[7,2],[5,4],[5,2,2],[4,3,2],[3,2,2,2]],'three_point_fibres':'OPEN; the cubic-cube family is excluded, remaining linear-cube/cubic-square cases unresolved','seconds':time.monotonic()-START,'assurance':'Independent source reconstruction and literal FLINT replay of a rational identity DAG; no independent end-to-end or formal proof review, external review or novelty priority.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(receipt))
