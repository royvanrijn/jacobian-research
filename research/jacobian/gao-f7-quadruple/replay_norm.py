"""Exact CRT-height and complete two-engine modular interpolation replay.
No multivariate resultant, Groebner calculation or prime search is performed.
"""
from pathlib import Path
import hashlib,json,os,platform,resource,subprocess,time,math
import sympy as S
HERE=Path(__file__).resolve().parent;START=time.monotonic()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def yes(c):
 if not c:raise AssertionError('norm replay failed')
src=json.loads((HERE/'generic-triple-input.json').read_text());r=json.loads((HERE/'primitive-remainder-candidate.json').read_text());n=json.loads((HERE/'crt-norm-candidate.json').read_text());a,s=S.symbols('a s');gt={ex:int(c) for ex,c in S.Poly(S.sympify(src['equations'][0]['polynomial']),a,s).terms()};rt={tuple(ex):int(c) for ex,c in r['terms']};nc=list(map(int,n['coefficients']))
yes(n['input_sha256']==sha(HERE/'primitive-remainder-candidate.json') and n['curve_sha256']==sha(HERE/'generic-triple-input.json'))
yes(max(i for i,j in gt)==8 and max(i for i,j in rt)==7)
degree_bound=7*max(j for i,j in gt)+8*max(j for i,j in rt);yes(degree_bound==11669 and n['degree_bound']==degree_bound and len(nc)-1==n['degree']<=degree_bound)
B=sum(abs(c) for c in gt.values())**7*sum(abs(c) for c in rt.values())**8;M=math.prod(n['primes']);yes(M>2*B and all(abs(c)<=B for c in nc));yes(len(set(n['primes']))==len(n['primes'])==345)
rows=[{} for _ in range(18)]
for (i,j),c in gt.items():rows[i][j]=c
for (i,j),c in rt.items():rows[9+i][j]=c
rows[17]={j:c for j,c in enumerate(nc)}
with (HERE/'norm-grid-input.txt').open('w') as fd:
 fd.write(str(degree_bound+1)+' '+str(len(n['primes']))+'\n');fd.write(' '.join(map(str,n['primes']))+'\n')
 for row in rows:
  fd.write(str(max(row)+1)+'\n');fd.write(' '.join(str(row.get(j,0)) for j in range(max(row)+1))+'\n')
# Use a matched installed FLINT header/library and its existing BLAS dependency.
# FLINT_ROOT and FLINT_BLAS_LIB can select portable replacements.
if os.environ.get('FLINT_ROOT'):root=Path(os.environ['FLINT_ROOT'])
else:root=next(z for z in Path.home().joinpath('.julia/artifacts').glob('*') if (z/'include/flint/flint.h').exists() and '#define __FLINT_VERSION_MINOR 3' in (z/'include/flint/flint.h').read_text())
blas=Path(os.environ['FLINT_BLAS_LIB']) if os.environ.get('FLINT_BLAS_LIB') else next(Path.home().joinpath('.julia/artifacts').glob('*/lib/libopenblas.0.3.33.so')).parent
command=['cc','-O3','-Wall','-Wextra','-I'+str(root/'include'),str(HERE/'verify_norm_modular.c'),'-L'+str(root/'lib'),'-Wl,-rpath,'+str(root/'lib'),'-Wl,-rpath-link,'+str(blas),'-lflint','-o',str(HERE/'verify_norm_modular')];subprocess.run(command,check=True,timeout=30)
env=os.environ.copy();env['LD_LIBRARY_PATH']=str(blas)+':'+str(root/'lib')+':'+env.get('LD_LIBRARY_PATH','')
def limits():
 resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3));resource.setrlimit(resource.RLIMIT_FSIZE,(10*1024**2,10*1024**2))
with (HERE/'norm-grid-replay.log').open('w') as fd:
 proc=subprocess.run([str(HERE/'verify_norm_modular'),str(HERE/'norm-grid-input.txt')],env=env,stdout=fd,stderr=subprocess.STDOUT,timeout=180,preexec_fn=limits)
yes(proc.returncode==0);lines=(HERE/'norm-grid-replay.log').read_text().splitlines();result=json.loads(lines[-1]);yes(result['result']=='PASS' and result['complete'] and result['primes']==345 and result['comparisons']==2*345*(degree_bound+1))
receipt={'result':'PASS','scope':'Exact integer identity N(s)=Res_a(G,R_final), by complete determinant interpolation modulo345 distinct trial-proved primes and M>2*coefficient_bound. No source/discriminant or chart-exclusion claim by this component alone.','degree_bound':degree_bound,'actual_degree':len(nc)-1,'row_l1_coefficient_bound_bits':B.bit_length(),'crt_modulus_bits':M.bit_length(),'primes':len(n['primes']),'points_per_prime':degree_bound+1,'determinant_comparisons':result['comparisons'],'determinant_engines':['FLINT3.3.1 nmod_mat_det','separate elementary C Gaussian elimination'],'multi_evaluation_engine':'FLINT3.3.1 univariate product-tree evaluation','checker_sha256':sha(Path(__file__)),'c_checker_sha256':sha(HERE/'verify_norm_modular.c'),'input_text_sha256':sha(HERE/'norm-grid-input.txt'),'curve_input_sha256':sha(HERE/'generic-triple-input.json'),'remainder_sha256':sha(HERE/'primitive-remainder-candidate.json'),'norm_sha256':sha(HERE/'crt-norm-candidate.json'),'compiler_command':command,'flint_root':str(root),'blas_lib':str(blas),'python':platform.python_version(),'seconds':time.monotonic()-START}
(HERE/'norm-grid-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
