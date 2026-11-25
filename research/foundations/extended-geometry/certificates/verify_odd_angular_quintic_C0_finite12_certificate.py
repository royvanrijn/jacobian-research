"""Exact normalized seven-term C=0 finite12 obstruction and sharpness.
Reconstructs native moments and verifies complete QQ chart certificates,
small exact resultants, Gauss coprimeness and the full smooth p103 point.
No Singular, SymPy, network, raw-moment input or CRT is required.
"""
from pathlib import Path
from fractions import Fraction
from collections import defaultdict
from math import factorial,gcd,lcm,isqrt
import gzip,json,hashlib,time,argparse
from flint import fmpq,fmpq_poly,fmpq_mpoly,fmpq_mpoly_ctx
EXPECTED='8acaba9157f3c38eccc9d67a038955e4d7908eb2b2af6f68fb9977a1db2812b8'
C=fmpq_mpoly_ctx.get(('a1','a3','b1','q0','lam'),'lex')
a1,a3,b,q,l=C.gens()
def comps(n,weights):
 if len(weights)==1:
  if n%weights[0]==0:yield(n//weights[0],)
  return
 for j in range(n//weights[0]+1):
  for tail in comps(n-j*weights[0],weights[1:]):yield(j,)+tail
def native(m):
 if m%2:return C.constant(0)
 out=defaultdict(Fraction)
 for row in comps(m//2,(1,2,3,1,1)):
  n1,n3,n5,nb,nl=row;r=n1+3*n3+5*n5+nb+nl
  for nq in range(r+1):
   nx=r-nq;s=nb+2*nl+nq+2*nx
   value=Fraction(factorial(m),factorial(nq)*factorial(nx))
   for count in row:value/=factorial(count)
   value*=Fraction(2**r*factorial(r)*factorial(2*s)*factorial(r+s),factorial(s)*factorial(2*r+2*s+1))
   out[(n1,n3,nb,nq,nl)]+=value
 return C.from_dict({e:fmpq(v.numerator,v.denominator)for e,v in out.items()if v})
def prim(n,d):
 if not n:return n,C.constant(1)
 common=n.gcd(d);n//=common;d//=common
 den=lcm(*(int(c.q)for c in n.to_dict().values()));con=gcd(*(abs(int(c.p)*(den//int(c.q)))for c in n.to_dict().values()))
 scale=fmpq(den,con)
 if n.leading_coefficient()<0:scale=-scale
 return n*scale,d*scale
def substitute(f,idx,n,d):
 degree=max((int(e[idx])for e in f.to_dict()),default=0);out=C.constant(0)
 for e,c in f.to_dict().items():
  ex=list(map(int,e));k=ex[idx];ex[idx]=0
  out+=C.from_dict({tuple(ex):c})*n**k*d**(degree-k)
 return out,d**degree
def coefficient(f,idx,power):
 out={}
 for e,v in f.to_dict().items():
  if int(e[idx])==power:
   ex=list(map(int,e));ex[idx]=0;out[tuple(ex)]=v
 return f.context().from_dict(out)
def univariate(f,idx):
 out=[0]*(max((int(e[idx])for e in f.to_dict()),default=-1)+1)
 for e,v in f.to_dict().items():
  assert all(int(e[j])==0 for j in range(len(e))if j!=idx)
  out[int(e[idx])]=v
 return fmpq_poly(out)
def inverse(f,P,idx):
 g,u,v=univariate(f,idx).xgcd(univariate(P,idx));assert g==1
 dd={}
 for k,c in enumerate(u):
  if c:
   ex=[0]*f.context().nvars();ex[idx]=k;dd[tuple(ex)]=c
 ans=f.context().from_dict(dd);assert(ans*f)%P==1;return ans
def verify_complementary(data):
 M={m:native(m)for m in range(1,13)}
 assert all(not M[m]for m in range(1,13,2))
 P1=7*q+3;P3=1615*q**3+2261*q**2+1197*q+231
 assert coefficient(M[2],0,1)==2*P1/105
 assert coefficient(M[4],1,1)==8*P3/4849845
 sol={s:tuple(fmpq_mpoly(data['solutions'][s][k],C)for k in ['numerator','denominator'])for s in ['a1','a3']}
 assert all(d for n,d in sol.values())
 N={}
 for m in range(2,13,2):
  n,d=substitute(M[m],0,*sol['a1']);nn,dd=substitute(n,1,*sol['a3']);nn,dd=prim(nn,d*dd)
  if m<=4:assert not nn
  else:N[m]=nn
  rem=dd
  for factor in [P1,P3]:
   while rem and not rem%factor:rem//=factor
  assert len(rem)==1 and all(all(v==0 for v in e)for e in rem.to_dict())
 U=fmpq_mpoly_ctx.get(('z','q0'),'lex');z,qq=U.gens();F={}
 for m in [6,8,10,12]:
  dd={};par=(m//2)%3
  for e,c in N[m].subs({'lam':0}).to_dict().items():
   assert e[0]==e[1]==e[4]==0 and int(e[2])>=par and (int(e[2])-par)%3==0
   dd[((int(e[2])-par)//3,int(e[3]))]=c
  f=U.from_dict(dd);den=lcm(*(int(v.q)for v in f.to_dict().values()));con=gcd(*(abs(int(v.p)*(den//int(v.q)))for v in f.to_dict().values()))
  f*=fmpq(den,con)
  if f.leading_coefficient()<0:f=-f
  F[m]=f;assert f==fmpq_mpoly(data['lambda0_F'][str(m)],U)
 R=[]
 for m in [8,10]:
  aa,bb=coefficient(F[6],0,1),coefficient(F[6],0,0);cc,dd=coefficient(F[m],0,1),coefficient(F[m],0,0)
  assert F[6].degrees()[0]<=1 and F[m].degrees()[0]<=1
  rr=aa*dd-cc*bb;assert -cc*F[6]+aa*F[m]==rr
  R.append(univariate(rr,1))
 target=fmpq_poly([3,7])**3*fmpq_poly([231,1197,2261,1615]);target/=target.leading_coefficient()
 checks=0
 for name,first,second in [('lambda0_nonzero_b_Bezout',*R),('lambda0_zero_b_Bezout',univariate(F[6].subs({'z':0}),1),univariate(F[12].subs({'z':0}),1))]:
  row=data[name];gg,uu,vv=(fmpq_poly([fmpq(x)for x in row[k]])for k in ['g','u','v'])
  assert gg==target and uu*first+vv*second==target
  assert (uu+1)*first+vv*second!=target;checks+=1
 # P1: derive the actual a3 solution directly from the native M4.
 A=fmpq_mpoly_ctx.get(('a1','a3','lam'),'lex');aa,cc,ll=A.gens();bb=-10*ll/11;qq=A.constant(-3)/7
 f4=M[4].compose(aa,cc,bb,qq,ll);pivot=coefficient(f4,1,1)
 assert len(pivot)==1 and all(all(v==0 for v in e)for e in pivot.to_dict())
 s3=-f4.subs({'a3':0})/pivot.leading_coefficient()
 small=fmpq_mpoly_ctx.get(('a1','lam'),'lex');res={}
 for m in range(2,13,2):
  f=M[m].compose(aa,s3,bb,qq,ll);assert all(e[1]==0 for e in f.to_dict())
  res[m]=small.from_dict({(int(e[0]),int(e[2])):v for e,v in f.to_dict().items()})
 assert not res[2] and not res[4]
 total=small.constant(0)
 for raw,m in zip(data['P1_cofactors'],[6,8,10,12]):total+=fmpq_mpoly(raw,small)*res[m]
 assert total==1 and total+res[12]!=1;checks+=1
 # P3: no polynomial division by the a3 pivot; M4 compatibility stays.
 A=fmpq_mpoly_ctx.get(('a3','b1','lam','q0'),'lex');cc,bb,ll,qq=A.gens();P=1615*qq**3+2261*qq**2+1197*qq+231
 args=(A.constant(0),cc,bb,qq,ll);n,d=sol['a1'];s1=n.compose(*args)*inverse(d.compose(*args),P,3)%P
 assert not M[2].compose(s1,cc,bb,qq,ll)%P
 total=A.constant(0)
 for raw,m in zip(data['P3_cofactors'],[4,6,8,10,12]):
  f=M[m].compose(s1,cc,bb,qq,ll)%P;total=(total+fmpq_mpoly(raw,A)*f)%P
 assert total==1 and(total+f)%P!=1;checks+=1
 return{'native_moments_reconstructed_through':12,'all_odd_native_moments_zero':True,'generic_lambda0_both_b_cases_exact_Bezout_common_factor':'monic(P1^3*P3)','P5_not_inverted':True,'P1_literal_native_unit':True,'P3_literal_native_unit_modulo_full_P3':True,'arithmetic_corruption_controls':checks}

def primitive(f):
 """Remove only a nonzero rational scalar; keep all variable factors."""
 if not f:return f
 den=lcm(*(int(c.q)for c in f.to_dict().values()))
 con=gcd(*(abs(int(c.p)*(den//int(c.q)))for c in f.to_dict().values()))
 out=f*fmpq(den,con)
 return -out if out.leading_coefficient()<0 else out
def prime(n):
 return n>=2 and all(n%d for d in range(2,isqrt(n)+1))
def conv(a,b,p):
 out=[0]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):out[i+j]=(out[i+j]+int(x)*int(y))%p
 return out
def q_content(f,index):
 ans=fmpq_poly([])
 groups={}
 for e,c in f.to_dict().items():
  key=tuple(int(v)for j,v in enumerate(e)if j!=index)
  groups.setdefault(key,{})[int(e[index])]=c
 for row in groups.values():
  a=[0]*(max(row)+1)
  for i,c in row.items():a[i]=c
  ans=ans.gcd(fmpq_poly(a))
 return ans
def modular_value(f,point,mod):
 total=0
 for e,c in f.to_dict().items():
  assert gcd(int(c.q),mod)==1
  term=int(c.p)*pow(int(c.q),-1,mod)%mod
  for x,k in zip(point,e):term=term*pow(x,int(k),mod)%mod
  total=(total+term)%mod
 return total
def modular_derivative(f,j,point,p):
 total=0
 for e,c in f.to_dict().items():
  if not e[j]:continue
  term=int(c.p)*pow(int(c.q),-1,p)*int(e[j])%p
  for i,(x,k)in enumerate(zip(point,e)):term=term*pow(x,int(k)-(1 if i==j else 0),p)%p
  total=(total+term)%p
 return total
def det_mod(mat,p):
 a=[row[:]for row in mat];value=1
 for i in range(len(a)):
  k=next(j for j in range(i,len(a))if a[j][i]%p)
  if k!=i:a[k],a[i]=a[i],a[k];value=-value
  d=a[i][i];value=value*d%p
  for j in range(i+1,len(a)):
   s=a[j][i]*pow(d,-1,p)%p
   a[j]=[(x-s*y)%p for x,y in zip(a[j],a[i])]
 return value%p
def verify_all(data,with_controls=True):
 result=verify_complementary(data)
 M={m:native(m)for m in range(1,13)}
 sol={s:tuple(fmpq_mpoly(data['solutions'][s][k],C)for k in ['numerator','denominator'])for s in ['a1','a3']}
 assert all(d for n,d in sol.values())
 H=fmpq_mpoly_ctx.get(('h','q0','z'),'lex');h,qq,z=H.gens()
 P1=7*qq+3;P3=1615*qq**3+2261*qq**2+1197*qq+231
 P5=930465*qq**5+2225025*qq**4+2314026*qq**3+1285570*qq**2+376805*qq+46189
 G={}
 for m in [6,8,10,12]:
  n,d=substitute(M[m],0,*sol['a1']);nn,dd=substitute(n,1,*sol['a3']);nn,dd=prim(nn,d*dd)
  dic={};par=(m//2)%3
  for e,c in nn.to_dict().items():
   assert e[0]==e[1]==0 and int(e[2])+int(e[4])>=par and (int(e[2])+int(e[4])-par)%3==0
   dic[(int(e[2]),int(e[3]),(int(e[2])+int(e[4])-par)//3)]=c
  G[m]=H.from_dict(dic)
 assert [tuple(map(int,G[m].degrees()))for m in [6,8,10,12]]==[(3,11,1),(4,16,1),(5,18,1),(6,25,2)]
 gen=data['generic_nonzero_lambda'];R={}
 for m,exps in zip([8,10,12],gen['first_stage_content_exponents']):
  full=G[6].resultant(G[m],2);content=P1**exps[0]*P3**exps[1]
  assert full and full%content==0
  assert q_content(full,1)==univariate(content,1)/univariate(content,1).leading_coefficient()
  expected_R=fmpq_mpoly(gen['primitive_R'][str(m)],H)
  assert primitive(full//content)==primitive(expected_R)
  R[m]=expected_R
 F=[]
 for other,exps,text,degree in zip([10,12],gen['outer_pivot_exponents'],gen['residual_F'],gen['residual_degrees']):
  full=R[8].resultant(R[other],0);content=P1**exps[0]*P3**exps[1]*P5**exps[2]
  assert full and full%content==0
  f=fmpq_mpoly(text,H);assert primitive(full//content)==primitive(f)
  ff=univariate(f,1);assert ff.degree()==degree
  # The retained vectors are primitive integer polynomials.
  assert all(int(c.q)==1 for c in ff)
  assert gcd(*(abs(int(c))for c in ff))==1
  F.append(ff)
 gauss=gen['Gauss'];p=int(gauss['modulus']);assert prime(p)
 assert all(int(f.leading_coefficient())%p for f in F)
 first=[int(c)%p for c in F[0]];second=[int(c)%p for c in F[1]]
 u=gauss['bezout_a'];v=gauss['bezout_b'];left=conv(u,first,p);right=conv(v,second,p);size=max(len(left),len(right))
 total=[((left[i]if i<len(left)else 0)+(right[i]if i<len(right)else 0))%p for i in range(size)]
 assert total[0]==1 and not any(total[1:])
 # Full P5 quotient unit: no h-pivot or leading coefficient is inverted.
 P5row=data['P5_nonzero_lambda'];A=fmpq_mpoly(P5row['A'],H)
 assert G[6]%P5==z*A
 unit=P5row['unit'];IA=fmpq_mpoly(unit['A_cofactor'],H);I8=fmpq_mpoly(unit['G8_cofactor'],H);I10=fmpq_mpoly(unit['G10_cofactor'],H)
 k=H.constant(int(unit['constant']));assert k
 assert (IA*A+I8*(G[8]%P5)+I10*(G[10]%P5))%P5==k
 # P1 and P3 are units at every P5 root, with exact polynomial inverses.
 inverse(P1,P5,1);inverse(P3,P5,1)
 result.update(generic_nonzero_lambda_finite12=True,exact_z_and_h_resultants_recomputed=True,
  generic_residual_degrees=[120,148],Gauss_prime=p,full_leading_degrees_preserved=True,
  P5_nonzero_lambda_finite10_literal_unit=True,complete_chart_partition=True)
 # Smooth point gives hypotheses for the separately written Newton proof.
 sharp=data['sharpness'];p=int(sharp['prime']);point=list(map(int,sharp['point']));assert prime(p)and len(point)==5
 assert all(gcd(int(c.q),p)==1 for f in M.values()for c in f.to_dict().values())
 assert all(modular_value(M[m],point,p)==0 for m in range(1,12))
 assert modular_value(M[12],point,p)==int(sharp['M12'])!=0
 orders=[2,4,6,8,10];J=[[modular_derivative(M[m],j,point,p)for j in range(5)]for m in orders]
 assert J==sharp['Jacobian']and det_mod(J,p)==int(sharp['Jacobian_determinant'])!=0
 # A literal first Newton correction checks the matrix/equation orientation.
 vals=[modular_value(M[m],point,p*p)for m in orders];assert all(v%p==0 for v in vals)
 aug=[row[:]+[(-v//p)%p]for row,v in zip(J,vals)]
 for i in range(5):
  k=next(j for j in range(i,5)if aug[j][i]);aug[k],aug[i]=aug[i],aug[k]
  fac=pow(aug[i][i],-1,p);aug[i]=[v*fac%p for v in aug[i]]
  for j in range(5):
   if i!=j:
    fac=aug[j][i];aug[j]=[(x-fac*y)%p for x,y in zip(aug[j],aug[i])]
 delta=[row[-1]for row in aug];lift=[x+p*d for x,d in zip(point,delta)]
 assert all(modular_value(M[m],lift,p*p)==0 for m in orders)
 assert modular_value(M[12],lift,p)==int(sharp['M12'])
 result.update(sharpness_prime=p,first11_native_moments_zero_mod103=True,M12_residue=int(sharp['M12']),
  full_native_Jacobian_determinant=int(sharp['Jacobian_determinant']),first_Newton_step_verified=True)
 if with_controls:
  import copy
  controls=[]
  def control(name,mutate):
   bad=copy.deepcopy(data);mutate(bad)
   try:verify_all(bad,False)
   except(AssertionError,ValueError,ZeroDivisionError,StopIteration):controls.append(name)
   else:raise AssertionError('Corruption accepted: '+name)
  control('a1_solution',lambda d:d['solutions']['a1'].__setitem__('numerator','('+d['solutions']['a1']['numerator']+')+1'))
  control('a3_solution',lambda d:d['solutions']['a3'].__setitem__('numerator','('+d['solutions']['a3']['numerator']+')+1'))
  control('lambda0_F6',lambda d:d['lambda0_F'].__setitem__('6','('+d['lambda0_F']['6']+')+1'))
  for key in ['lambda0_nonzero_b_Bezout','lambda0_zero_b_Bezout']:
   control(key,lambda d,k=key:d[k]['u'].__setitem__(0,str(fmpq(d[k]['u'][0])+1)))
  for key in ['P1_cofactors','P3_cofactors']:
   control(key,lambda d,k=key:d[k].__setitem__(0,'('+d[k][0]+')+1'))
  control('generic_R8',lambda d:d['generic_nonzero_lambda']['primitive_R'].__setitem__('8','('+d['generic_nonzero_lambda']['primitive_R']['8']+')+1'))
  control('generic_F120',lambda d:d['generic_nonzero_lambda']['residual_F'].__setitem__(0,'('+d['generic_nonzero_lambda']['residual_F'][0]+')+1'))
  control('generic_Gauss_Bezout',lambda d:d['generic_nonzero_lambda']['Gauss']['bezout_a'].__setitem__(0,(int(d['generic_nonzero_lambda']['Gauss']['bezout_a'][0])+1)%32003))
  control('Gauss_composite_modulus',lambda d:d['generic_nonzero_lambda']['Gauss'].__setitem__('modulus',32005))
  control('P5_A',lambda d:d['P5_nonzero_lambda'].__setitem__('A','('+d['P5_nonzero_lambda']['A']+')+1'))
  control('P5_unit_A_cofactor',lambda d:d['P5_nonzero_lambda']['unit'].__setitem__('A_cofactor','('+d['P5_nonzero_lambda']['unit']['A_cofactor']+')+1'))
  control('P5_unit_G8_cofactor',lambda d:d['P5_nonzero_lambda']['unit'].__setitem__('G8_cofactor','('+d['P5_nonzero_lambda']['unit']['G8_cofactor']+')+1'))
  control('sharpness_point',lambda d:d['sharpness']['point'].__setitem__(0,(int(d['sharpness']['point'][0])+1)%103))
  control('sharpness_Jacobian',lambda d:d['sharpness'].__setitem__('Jacobian_determinant',95))
  result['corruption_controls']=controls
 return result

def main():
 if not __debug__:raise SystemExit("Optimized Python is unsupported: assertions are required")
 parser=argparse.ArgumentParser();parser.add_argument('--certificate',type=Path,default=Path(__file__).with_name('odd-angular-quintic-C0-finite12-certificate.json.gz'));parser.add_argument('--output',type=Path);args=parser.parse_args()
 if args.output and args.output.exists():raise SystemExit('Refusing to overwrite receipt')
 start=time.monotonic();blob=args.certificate.read_bytes()
 if hashlib.sha256(blob).hexdigest()!=EXPECTED:raise SystemExit('Certificate SHA256 mismatch')
 result=verify_all(json.loads(gzip.decompress(blob)));result.update(status='PASS',elapsed_seconds=time.monotonic()-start,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),certificate_sha256=EXPECTED,scope='Complete seven-term C=0 normalized finite12 obstruction plus exact smooth-prefix11 sharpness hypotheses; written Newton and local finite-Q residue arguments justify characteristic-zero witnesses. No eventual or global quintic theorem.')
 out=json.dumps(result,indent=2)+'\n'
 if args.output:args.output.write_text(out)
 print(out)
if __name__=='__main__':main()
