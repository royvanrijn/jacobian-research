"""Bounded exact integer determinant interpolation for residual gcd>=3.
The rectangular order2 Sylvester map has13 rows and15 columns.
Every13x13 minor vanishes when gcd(Q9,Q9_prime) has degree>=3.
"""
from pathlib import Path
from fractions import Fraction
import json,time,hashlib,resource,subprocess,sys
HERE=Path(__file__).resolve().parent;START=time.monotonic()
def emit(x):print(json.dumps(x,default=int),flush=True)
def model(row):
 from flint import fmpz_poly
 coeff=[[[0]*9 for _ in range(3)] for _ in range(10)]
 for (i,j,k),c in row['terms']:coeff[i][k][j]=int(c)
 F=[[fmpz_poly(v) for v in c] for c in coeff];derivative=[[p*(i+1) for p in F[i+1]] for i in range(9)];zero=[fmpz_poly([])]*3
 mat=[[zero for j in range(15)] for i in range(13)]
 for i in range(6):
  for j,c in enumerate(reversed(F)):mat[i][i+j]=c
 for i in range(7):
  for j,c in enumerate(reversed(derivative)):mat[6+i][i+j]=c
 return mat

def tropical(mat,cols):
 n=len(mat);dp={0:0}
 for row in mat:
  nxt={}
  for mask,total in dp.items():
   for j,col in enumerate(cols):
    if mask>>j&1:continue
    deg=max(int(v.degree()) for v in row[col])
    if deg<0:continue
    key=mask|1<<j;value=total+deg
    if value>nxt.get(key,-1):nxt[key]=value
  dp=nxt
 return dp[(1<<n)-1]

def evaluated(mat,cols,value,cubic):
 simple=[[[int(p(value)) for p in row[c]] for c in cols] for row in mat]
 if not cubic:return [[v[0] for v in row] for row in simple]
 out=[[0]*39 for _ in range(39)]
 for i,row in enumerate(simple):
  for j,(c0,c1,c2) in enumerate(row):
   block=[[3*c0,2*c2,2*c1],[3*c1,3*c0,2*c2],[3*c2,3*c1,3*c0]]
   for u in range(3):
    for v in range(3):out[3*i+u][3*j+v]=block[u][v]
 return out

def interpolate(values):
 from flint import fmpq_poly
 out=fmpq_poly([]);basis=fmpq_poly([1]);diff=list(values)
 for k in range(len(values)):
  out+=basis*diff[0];diff=[b-a for a,b in zip(diff,diff[1:])];basis=basis*fmpq_poly([-k,1])/(k+1)
 coeff=[Fraction(str(v)) for v in out.coeffs()];assert all(v.denominator==1 for v in coeff);return [v.numerator for v in coeff]

def worker():
 from flint import fmpz_poly,fmpz_mat
 source=json.loads((HERE/'pole-source-candidate.json').read_text());records=[]
 for row in source['pole_charts']:
  mat=model(row);cubic=row['name']=='s-cubic';polys=[]
  for last in [12,13,14]:
   cols=list(range(12))+[last];trop=tropical(mat,cols);bound=trop*(3 if cubic else 1);emit({'stage':'matrix','chart':row['name'],'last_column':last,'matrix_size':39 if cubic else 13,'degree_bound':bound,'seconds':time.monotonic()-START})
   values=[int(fmpz_mat(evaluated(mat,cols,j,cubic)).det()) for j in range(bound+1)];co=interpolate(values);poly=fmpz_poly(co);content=int(poly.content());primitive=poly//content if content else poly;lin=fmpz_poly([1,2]);power=0
   while primitive:
    q,r=divmod(primitive,lin)
    if r:break
    primitive=q;power+=1
   polys.append({'last_column':last,'tropical_degree_bound':trop,'degree_bound':bound,'coefficients':list(map(str,co)),'integer_content':str(content),'denominator_factor_power':power,'stripped_coefficients':list(map(str,primitive.coeffs()))})
  gcd=fmpz_poly([])
  for v in polys:gcd=gcd.gcd(fmpz_poly(list(map(int,v['stripped_coefficients']))))
  records.append({'name':row['name'],'minors':polys,'stripped_gcd_coefficients':list(map(str,gcd.coeffs()))});emit({'stage':'chart_saved','name':row['name'],'gcd_degree':gcd.degree(),'degrees':[len(v['coefficients'])-1 for v in polys],'seconds':time.monotonic()-START})
  (HERE/'pole-norms-checkpoint.json').write_text(json.dumps({'status':'UNVERIFIED_EXACT_CANDIDATE','charts':records},separators=(',',':'))+'\n')
 out={'status':'UNVERIFIED_EXACT_CANDIDATE','source_sha256':hashlib.sha256((HERE/'pole-source-candidate.json').read_bytes()).hexdigest(),'charts':records,'boundary':'Exact determinant interpolation candidate only. Independent full-grid determinant, source, degree and saturation replay required.'}
 (HERE/'pole-norms-candidate.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
if __name__=='__main__':
 if '--worker' in sys.argv:worker()
 else:
  def limits():
   resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(50*1024**2,50*1024**2))
  receipt={'cap_seconds':120,'address_cap_GiB':4,'per_file_cap_MiB':50,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'theorem_exclusion':False,'motivation':'New residual-order2 Sylvester minors on a complete linear pole chart, replacing the retained six-variable factorization searches. At most three degree312 determinant norms per chart, with exact interpolation.'}
  with (HERE/'pole-norms.log').open('wb') as fd:
   try:
    proc=subprocess.run([sys.executable,str(Path(__file__)),'--worker'],stdout=fd,stderr=subprocess.STDOUT,timeout=120,preexec_fn=limits);receipt.update(outcome='completed',returncode=proc.returncode)
   except subprocess.TimeoutExpired:receipt.update(outcome='timeout',returncode=None)
  receipt['seconds']=time.monotonic()-START
  for fn in ['pole-source-candidate.json','pole-norms-candidate.json','pole-norms-checkpoint.json']:
   p=HERE/fn
   if p.exists():receipt[fn]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
  (HERE/'pole-norms-run.json').write_text(json.dumps(receipt,indent=2)+'\n');emit(receipt);print((HERE/'pole-norms.log').read_text()[-3000:])
