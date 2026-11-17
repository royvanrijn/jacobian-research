#!/usr/bin/env python3
from pathlib import Path
from fractions import Fraction as Q
import hashlib,itertools,json,math,platform
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];COUNTS={}
def check(key,n=1):COUNTS[key]=COUNTS.get(key,0)+n
def modq(a):
 a=Q(a);assert a.denominator%7
 return a.numerator%7*pow(a.denominator,-1,7)%7
S=s.Symbol('S')
rows=[
[-1-S,3+S,-S,1+S,4-S],
[2*S,3+S,3,5+5*S,3+6*S],
[0,6*S,3+6*S,5+3*S,6+6*S],
[0,0,3*S,2+5*S,5+6*S],
[0,0,0,S,2+6*S]]
def poly7(f):return s.Poly(f,S,modulus=7)
for b in range(5):
 f=[(j+1)*pow(6-j,b,7)%7 for j in range(7)]
 D=[sum((-1)**j*math.comb(a,j)*f[j] for j in range(a+1))%7 for a in range(7)]
 assert D[6]==0 and D[5]==(-int(b==4))%7 and sum(D[:6])%7==0
 check('universal_difference_degree_and_hyperplane_controls',3)
 for a in range(1,6):
  assert poly7(rows[a-1][b]-a*(D[a]-S*D[a-1])).is_zero
  check('universal_pencil_difference_identities')
M=s.Matrix(rows);g=sum(S**i for i in range(6))
assert poly7(M.det()+g).is_zero;check('universal_main_minor_polynomial_identity')
R0=s.Matrix([[6,1,6,1,6]])
R1=s.Matrix([[0,0,0,0,S]])
E0=M.copy();E0[4,:]=R0
E1=M.copy();E1[4,:]=R1
assert poly7(E0.det()-4).is_zero
assert poly7(E1.det()-3*(g-1)).is_zero
assert poly7(E1.det()+3*M.det()-4).is_zero
check('universal_explicit_affine_minor_coefficient_identities',3)
kernels={2:[1,0,3,0,1],3:[0,1,2,5,1],4:[2,2,6,4,1],5:[0,6,2,2,1],6:[2,5,6,3,1]}
for a,q in kernels.items():
 T=M.subs(S,a)
 assert all(int(x)%7==0 for x in T*s.Matrix(q))
 assert 6*sum(x*pow(6,j,7) for j,x in enumerate(q))%7==a
 assert any(int(T.extract(I,J).det())%7 for I in itertools.combinations(range(5),4) for J in itertools.combinations(range(5),4))
 check('complete_finite_field_kernel_states',3)
bad=kernels[2][:];bad[0]=(bad[0]+1)%7
assert any(int(x)%7 for x in M.subs(S,2)*s.Matrix(bad))
check('corrupted_kernel_rejected')
def vp(a):
 if not a:return None
 a=Q(a);ans=0
 while a.numerator%7==0:a/=7;ans+=1
 while a.denominator%7==0:a*=7;ans-=1
 return ans
def act(f,c):
 out=[0]*(len(f)+1)
 for i,x in enumerate(f):out[i]+=(c-i)*x;out[i+1]+=x
 return out
def rawcols(n,c=None):
 B=[1];c=n+5 if c is None else c
 for ell in range(n):
  out=act(B,c)
  for i,x in enumerate(B):out[i]+=6*ell*x
  B=out
 cols=[]
 for a in range(5):
  if a:B=act(B,c)
  cols.append(B)
 return cols
def fpdigits(j):
 d=[]
 while j:d.append(j%7);j//=7
 return d or [0]
profiles=[]
v=s.Symbol('v')
for r in (1,2,3):
 L=7**r;n=L+1;c=n+5;den=6**n*math.factorial(n);cols=rawcols(n)
 trace=1+sum((-1)**i*v**((7**i-1)//6) for i in range(1,r+1))
 for b,B in enumerate(cols):
  f=[(j+1)*pow(6-j,b,7)%7 for j in range(7)]
  D=[sum((-1)**j*math.comb(a,j)*f[j] for j in range(a+1))%7 for a in range(7)]
  for j,x in enumerate(B):
   digits=fpdigits(j);a=sum(digits)
   if a>6:continue
   dj=Q(math.factorial(j)*x,den);assert not dj or vp(dj)>=0
   high=sum(digits[1:]);j0=digits[0]
   if not a:expected=0
   elif high==0:expected=a*(D[a]-D[a-1])%7
   elif high==1:expected=-D[a-1]%7
   else:expected=0
   assert modq(dj)==expected,(r,b,j,digits)
   check('literal_all_small_digit_sum_Lucas_difference_controls')
  for a in range(1,6):
   actual=sum(modq(Q(math.factorial(a)*7**ell*B[a+6*ell],den))*v**ell
              for ell in range(max(0,(len(B)-1-a)//6+1)))
   expected=rows[a-1][b].subs(S,trace) if hasattr(rows[a-1][b],'subs') else rows[a-1][b]
   assert s.Poly(actual-expected,v,modulus=7).is_zero
   check('literal_complete_universal_boundary_row_identities')
  divided=sum(modq(Q(7**ell*B[6*ell],7*den))*v**ell for ell in range((len(B)-1)//6+1))
  expected=6*pow(6,b,7)+(v*trace if b==4 else 0)
  assert s.Poly(divided-expected,v,modulus=7).is_zero
  check('literal_complete_saturated_row_zero_identities')
 profiles.append({'r':r,'terminal':n,'all_column_coefficients_checked':sum(len(B) for B in cols)})
# r=0 terminal2, an actual degree2 Gamma-moment polynomial with strict variance.
mu1=math.factorial(7)//math.factorial(1);mu2=math.factorial(13)//math.factorial(1)
assert mu2>mu1**2
check('r_zero_quadratic_variance_control')

gamma_profiles=[]
for r in (1,2):
 n=7**r+1;delta=7**(r+2);den=6**n*math.factorial(n)
 base=rawcols(n);pert=rawcols(n,n+5+delta)
 for b in range(5):
  for a in range(6):
   for ell in range(max(0,(len(base[b])-1-a)//6+1)):
    j=a+6*ell
    diff=Q(math.factorial(a)*7**ell*(pert[b][j]-base[b][j]),den)
    assert not diff or vp(diff)>=2
    if a==0:assert not diff or vp(diff/7)>=1
    check('literal_Gamma_disk_complete_boundary_coefficient_controls')
 for k in range(1,n+1):
  z=Q(1)
  for j in range(k):z*=Q(delta,6)-j
  z/=math.factorial(k)
  assert vp(z)==r+2-(vp(k) or 0)>=2
  check('literal_Gamma_Vandermonde_valuation_controls')
  phys=7*math.comb(n,k)*math.prod(range(n,n+6*k))
  changed=7*math.comb(n,k)*math.prod(range(n+delta,n+delta+6*k))
  assert phys%7**k==changed%7**k==0
  hp,hg=phys//7**k,changed//7**k
  assert (hp-hg)%7==0
  if k in (1,n):assert hg%7
  check('literal_Gamma_scaled_moment_congruences')
 gamma_profiles.append({'r':r,'terminal':n,'beta':1+delta})
# Ramified positive algebraic beta=1+343sqrt(7), one7-adic place.
def padd(a,b):return (a[0]+b[0],a[1]+b[1])
def pmul(a,b):return (a[0]*b[0]+7*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def pact(f,c):
 out=[(0,0)]*(len(f)+1)
 for j,x in enumerate(f):
  out[j]=padd(out[j],pmul((c[0]-j,c[1]),x));out[j+1]=padd(out[j+1],x)
 return out
n=8;c=(13,343);B=[(1,0)];ram=[]
for ell in range(n):
 out=pact(B,c)
 for j,x in enumerate(B):out[j]=padd(out[j],(6*ell*x[0],6*ell*x[1]))
 B=out
for a in range(5):
 if a:B=pact(B,c)
 ram.append(B)
base=rawcols(n);den=6**n*math.factorial(n)
for b in range(5):
 for a in range(6):
  for ell in range(max(0,(len(base[b])-1-a)//6+1)):
   j=a+6*ell;weight=Q(math.factorial(a)*7**ell,den)
   aa=weight*(ram[b][j][0]-base[b][j]);bb=weight*ram[b][j][1]
   va=vp(aa);vb=vp(bb)
   assert va is None or va>=2
   assert vb is None or vb+Q(1,2)>=2
   check('ramified_Gamma_disk_boundary_coefficient_controls')
# The general Gamma bridge must use6/gamma, rather than6/n.
gamma=Q(8+343);n=8
for k in range(n+1):
 rising=lambda a,m:math.prod(a+j for j in range(m))
 left=math.comb(n,k)*rising(gamma+1,6*k)
 right=math.comb(n,k)*rising(gamma,6*k)*(1+Q(6*k,gamma))
 assert left==right
 check('literal_correct_Gamma_upper_bridge_coefficients')
assert math.comb(n,1)*rising(gamma+1,6)!=math.comb(n,1)*rising(gamma,6)*(1+Q(6,n))
check('incorrect_ordinary_bridge_in_Gamma_rejected')

imports=['PHYSICAL_PRIME_TOTAL_DEGREE_REDUCTION.tex','gap_three_hurwitz_structure.tex','ORDINARY_DEGREE_SEVEN_GAP_SIX_FRONTIER.tex']
out={'passed':True,'scope':'All-r proof identities and explicit affine scalar minor; literal profiles only verify normalization.',
 'uniform_terminal_family':'N=7^r+1 for every integer r>=0','all_index_gap_six_SFC_claimed':False,
 'controls':COUNTS,'literal_profiles':profiles,'Gamma_disk_profiles':gamma_profiles,
 'Gamma_disk_scope':'positive real algebraic beta at one embedding with v7(beta-1)>=r+2,r>=1',
 'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (HERE/'GAP_SIX_THIN_DIGIT_SFC.tex',Path(__file__).resolve())},
 'import_sha256':{name:hashlib.sha256((ROOT/'research/factorial/image-rigidity'/name).read_bytes()).hexdigest() for name in imports},
 'python':platform.python_version(),'sympy':s.__version__,
 'independent_replay':False,'formal_verification':False,'external_review':False,
 'novelty_claim':False,'priority_claim':False,
 'command':'.venv/bin/python research/factorial/image-rigidity/gap_six_thin_digit_sfc.py'}
out['source_sha256'].update(out.pop('import_sha256'))
(HERE/'gap_six_thin_digit_sfc_verification.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({'passed':True,'controls':COUNTS,'total_controls':sum(COUNTS.values())},indent=2))
