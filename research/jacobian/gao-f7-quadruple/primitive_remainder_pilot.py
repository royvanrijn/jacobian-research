"""Bounded fraction-free reduction of residual discriminant modulo G.
Every removed univariate coefficient content is retained, including constants.
No cancelled factor is silently discarded from an exceptional projector.
"""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
HERE=Path(__file__).resolve().parent;START=time.monotonic()
def emit(x):print(json.dumps(x,default=int),flush=True)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def worker():
 import sympy as S
 from flint import fmpz_mpoly_ctx,fmpz_poly
 a,s=S.symbols('a s');ctx=fmpz_mpoly_ctx.get(['a','s']);av,sv=ctx.gens()
 data=json.loads((HERE/'residual-discriminant-candidate.json').read_text())
 assert all(ex[0]==0 for ex,c in data['terms'])
 R=ctx.from_dict({(ex[1],ex[2]):int(c) for ex,c in data['terms']})
 src=json.loads((HERE/'generic-triple-input.json').read_text());G=ctx.from_dict({ex:int(c) for ex,c in S.Poly(S.sympify(src['equations'][0]['polynomial']),a,s).terms()})
 assert G.degrees()[0]==8
 def coeff(P,n):return ctx.from_dict({(0,int(j)):int(c) for (i,j),c in P.to_dict().items() if i==n})
 ell=coeff(G,8);removed=[]
 def primitive(P):
  td=P.to_dict();rows={}
  for (i,j),c in td.items():rows.setdefault(int(i),{})[int(j)]=int(c)
  content=fmpz_poly([])
  for row in rows.values():
   u=fmpz_poly([row.get(j,0) for j in range(max(row)+1)]);content=content.gcd(u)
   if content==1:break
  cc=list(map(int,content.coeffs()));C=ctx.from_dict({(0,j):c for j,c in enumerate(cc) if c})
  Q,rem=divmod(P,C);assert not rem
  return Q,cc
 R,cc=primitive(R);removed.append({'stage':'initial','content':list(map(str,cc))})
 emit({'stage':'initial_primitive','degrees':list(R.degrees()),'terms':len(R.to_dict()),'content_degree':len(cc)-1,'seconds':time.monotonic()-START})
 steps=0
 while R and R.degrees()[0]>=8:
  n=int(R.degrees()[0]);R=ell*R-coeff(R,n)*av**(n-8)*G
  R,cc=primitive(R);steps+=1
  removed.append({'stage':steps,'cancelled_a_degree':n,'content':list(map(str,cc))})
  emit({'stage':'remainder','step':steps,'degrees':list(R.degrees()),'terms':len(R.to_dict()),'content_degree':len(cc)-1,'seconds':time.monotonic()-START})
  if steps%10==0:
   (HERE/'primitive-remainder-checkpoint.json').write_text(json.dumps({'steps':steps,'degrees':list(map(int,R.degrees())),'removed_contents':removed,'terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(R.to_dict().items())]},separators=(',',':'))+'\n')
 out={'status':'UNVERIFIED_EXACT_CANDIDATE','input_sha256':digest(HERE/'residual-discriminant-candidate.json'),'curve_sha256':digest(HERE/'generic-triple-input.json'),'variables':['a','s'],'steps':steps,'degrees':list(map(int,R.degrees())),'removed_contents':removed,'terms':[[list(map(int,ex)),str(c)] for ex,c in sorted(R.to_dict().items())],'identity':'ell*R_i - leading_a(R_i)*a^(deg_a(R_i)-8)*G = content_i(s)*R_(i+1); initial discriminant = initial_content*R_0','boundary':'All contents, leading coefficient zeros and source-coordinate poles must be retained in any global candidate set.'}
 (HERE/'primitive-remainder-candidate.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
 emit({'stage':'saved','steps':steps,'degrees':list(R.degrees()),'terms':len(R.to_dict()),'seconds':time.monotonic()-START})
if '--worker' in sys.argv:worker()
else:
 def limits():
  resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(100*1024**2,100*1024**2))
 receipt={'name':'primitive-remainder-pilot','cap_seconds':120,'address_cap_GiB':4,'per_file_cap_MiB':100,'script_sha256':digest(Path(__file__)),'input_sha256':digest(HERE/'residual-discriminant-candidate.json'),'theorem_exclusion':False}
 (HERE/'primitive-remainder-pilot.json').write_text(json.dumps({**receipt,'outcome':'RUNNING'},indent=2)+'\n')
 with (HERE/'primitive-remainder-pilot.log').open('wb') as fd:
  try:
   proc=subprocess.run([sys.executable,str(Path(__file__)),'--worker'],stdout=fd,stderr=subprocess.STDOUT,timeout=120,preexec_fn=limits);receipt.update(outcome='completed',returncode=proc.returncode)
  except subprocess.TimeoutExpired:receipt.update(outcome='timeout',returncode=None)
 receipt['seconds']=time.monotonic()-START
 for name in ['primitive-remainder-candidate.json','primitive-remainder-checkpoint.json']:
  artifact=HERE/name
  if artifact.exists():receipt[name]={'bytes':artifact.stat().st_size,'sha256':digest(artifact)}
 (HERE/'primitive-remainder-pilot.json').write_text(json.dumps(receipt,indent=2)+'\n');emit(receipt)
 print((HERE/'primitive-remainder-pilot.log').read_text()[-4000:])
