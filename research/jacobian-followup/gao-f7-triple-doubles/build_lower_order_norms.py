"""Bounded lower-order minor norms for complete pole-chart stratification.
Order1 tests residual gcd>=2; order0 records the residual discriminant.
"""
from pathlib import Path
import hashlib,json,resource,subprocess,sys,time
from build_pole_norms import tropical,interpolate
HERE=Path(__file__).resolve().parent;START=time.monotonic()
def emit(x):print(json.dumps(x,default=int),flush=True)
def model(row,order):
 from flint import fmpz_poly
 coeff=[[[0]*9 for _ in range(3)] for _ in range(10)]
 for (i,j,k),c in row['terms']:coeff[i][k][j]=int(c)
 F=[[fmpz_poly(v) for v in c] for c in coeff];D=[[p*(j+1) for p in F[j+1]] for j in range(9)];zero=[fmpz_poly([])]*3;n=17-2*order;cols=17-order;M=[[zero for _ in range(cols)] for i in range(n)]
 for i in range(8-order):
  for j,v in enumerate(reversed(F)):M[i][i+j]=v
 for i in range(9-order):
  for j,v in enumerate(reversed(D)):M[8-order+i][i+j]=v
 return M

def matrix(M,cols,value,cubic):
 rows=[[[int(v(value)) for v in row[j]] for j in cols] for row in M];n=len(rows)
 if not cubic:return [[v[0] for v in row] for row in rows]
 out=[[0]*(3*n) for _ in range(3*n)]
 for i,row in enumerate(rows):
  for j,(c0,c1,c2) in enumerate(row):
   block=[[3*c0,2*c2,2*c1],[3*c1,3*c0,2*c2],[3*c2,3*c1,3*c0]]
   for u in range(3):
    for v in range(3):out[3*i+u][3*j+v]=block[u][v]
 return out

def worker():
 from flint import fmpz_mat,fmpz_poly,nmod_poly
 src=json.loads((HERE/'pole-source-candidate.json').read_text());records=[]
 for order in [1,0]:
  for row in src['pole_charts']:
   M=model(row,order);n=len(M);cubic=row['name']=='s-cubic';polys=[]
   for last in range(n-1,len(M[0])):
    cols=list(range(n-1))+[last];trop=tropical(M,cols);bound=trop*(3 if cubic else 1);emit({'stage':'matrix','order':order,'chart':row['name'],'last_column':last,'size':n*(3 if cubic else 1),'degree_bound':bound,'seconds':time.monotonic()-START})
    values=[int(fmpz_mat(matrix(M,cols,a,cubic)).det()) for a in range(bound+1)];co=interpolate(values);F=fmpz_poly(co);content=int(F.content());R=F//content if content else F;power=0
    while R:
     Q,rem=divmod(R,fmpz_poly([1,2]))
     if rem:break
     R=Q;power+=1
    polys.append({'last_column':last,'tropical_degree_bound':trop,'degree_bound':bound,'coefficients':list(map(str,co)),'integer_content':str(content),'denominator_factor_power':power,'stripped_coefficients':list(map(str,R.coeffs()))})
   gcd=fmpz_poly([])
   for x in polys:gcd=gcd.gcd(fmpz_poly(list(map(int,x['stripped_coefficients']))))
   record={'name':row['name'],'order':order,'minors':polys,'stripped_gcd_coefficients':list(map(str,gcd.coeffs()))}
   if order==0:
    R=fmpz_poly(list(map(int,polys[0]['stripped_coefficients'])));p=1000000007;rm=nmod_poly([int(c)%p for c in R.coeffs()],p);record['stripped_squarefree_gcd_degree_mod']=int(rm.gcd(rm.derivative()).degree())
   records.append(record);emit({'stage':'record','name':row['name'],'order':order,'degrees':[len(x['coefficients'])-1 for x in polys],'stripped_degrees':[len(x['stripped_coefficients'])-1 for x in polys],'gcd_degree':gcd.degree(),'squarefree_gcd':record.get('stripped_squarefree_gcd_degree_mod'),'seconds':time.monotonic()-START})
   (HERE/'lower-order-norms-checkpoint.json').write_text(json.dumps({'charts':records},separators=(',',':'))+'\n')
 out={'status':'UNVERIFIED_EXACT_CANDIDATE','source_sha256':hashlib.sha256((HERE/'pole-source-candidate.json').read_bytes()).hexdigest(),'charts':records}
 (HERE/'lower-order-norms-candidate.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
if __name__=='__main__':
 if '--worker' in sys.argv:worker()
 else:
  def limits():
   resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(50*1024**2,50*1024**2))
  rec={'cap_seconds':60,'address_cap_GiB':4,'per_file_cap_MiB':50,'theorem_exclusion':False,'motivation':'Order2 norms completed quickly with coprime residuals. Test the stronger order1 rank bound, and compute order0 discriminant to count the pole-curve exceptions, in a new finite bounded interpolation run.'}
  with (HERE/'lower-order-norms.log').open('wb') as fd:
   try:
    proc=subprocess.run([sys.executable,str(Path(__file__)),'--worker'],stdout=fd,stderr=subprocess.STDOUT,timeout=60,preexec_fn=limits);rec.update(outcome='completed',returncode=proc.returncode)
   except subprocess.TimeoutExpired:rec.update(outcome='timeout',returncode=None)
  rec['seconds']=time.monotonic()-START
  for fn in ['lower-order-norms-candidate.json','lower-order-norms-checkpoint.json']:
   p=HERE/fn
   if p.exists():rec[fn]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
  (HERE/'lower-order-norms-run.json').write_text(json.dumps(rec,indent=2)+'\n');emit(rec);print((HERE/'lower-order-norms.log').read_text()[-4000:])
