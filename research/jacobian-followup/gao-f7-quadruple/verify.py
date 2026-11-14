"""Fail-closed replay for F7 quadruple incidence: generic degree8, finite degenerations and all pole points.

No Groebner basis, symbolic resultant, rational reconstruction or factorizer
is run. Modular unit tests are confined to specified finite FREE algebras,
whose defining leading coefficients are checked invertible at the prime.
"""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,platform,time
import sympy as S
import flint
from flint import fmpz_mpoly_ctx,fmpz_poly,fmpz_mat,nmod_poly,nmod_mat
if not __debug__:raise RuntimeError('Run without -O')
HERE=Path(__file__).resolve().parent;START=time.monotonic();checks=[]
def yes(name,condition):
 assert condition,name
 checks.append(name)
def announce(stage,**extra):print(json.dumps({'stage':stage,'seconds':time.monotonic()-START,**extra}),flush=True)
cert=json.loads((HERE/'exceptional-certificate.json').read_text());exc=json.loads((HERE/'exceptional-source.json').read_text());source=json.loads((HERE/'quadruple-input.json').read_text())
s,a,b,z=S.symbols('s a b z');V=[s,a,b]
ctx=fmpz_mpoly_ctx.get(['x','s','a','b','xi']);x,ss,aa,bb,xi=ctx.gens()
def native(expression):return ctx.from_dict({(0,)+ex+(0,):int(c) for ex,c in S.Poly(expression,*V,domain=S.ZZ).terms()})
G2=a*a+a+s*b+s**3;G3=2*s*G2+a+s**3+s*s*b;G4=3*s*s*G2+b+s**4
X1=S.expand(-S.det(S.Matrix([G2,G3,G4]).jacobian(V)))
G2n,G3n,G4n,X1n=map(native,[G2,G3,G4,X1]);T=ss+x

def eliminant(Y):
 Y1,Y2,Y3,Y4=Y;bt=Y4-3*Y2*T*T+2*Y1*T**3-T**4
 at=Y3-2*Y2*T+Y1*T*T-T**3-T*T*bt
 return at*at+at+T*bt+T**3+Y1*T-Y2,at,bt
full=[xi,G2n+ss*xi,G3n+ss*ss*xi,G4n+ss**3*xi]
Rfull,afull,bfull=eliminant(full)
def x_coefficient(poly,j):return ctx.from_dict({(0,)+ex[1:]:int(c) for ex,c in poly.to_dict().items() if ex[0]==j})
yes('root a reconstruction',x_coefficient(afull,0)==aa)
yes('root b reconstruction',x_coefficient(bfull,0)==bb)
yes('root equation reconstruction',x_coefficient(Rfull,0)==0)
yes('all multiple roots covered without localization',x_coefficient(Rfull,1)==xi-X1n)
X=[X1n,G2n+ss*X1n,G3n+ss*ss*X1n,G4n+ss**3*X1n];R,acrit,_=eliminant(X)
yes('critical source constant and linear coefficients',x_coefficient(R,0)==x_coefficient(R,1)==0)
yes('source variable order',source['variables']==['s','a','b'])
yes('two necessary derivative equations',[r['derivative'] for r in source['equations']]==[2,3])
f=[];fexpr=[]
for j,row in enumerate(source['equations'],2):
 expression=S.sympify(row['polynomial']);F=native(expression);C=x_coefficient(R,j);terms=F.to_dict();ex=next(iter(terms));ratio=Fraction(int(C.to_dict().get(ex,0)),int(terms[ex]))
 yes('fresh necessary derivative '+str(j),ratio!=0 and ratio.denominator*C==ratio.numerator*F)
 f.append(F);fexpr.append(expression)
Q=ctx.constant(0)
yes('all eight residual quotient coefficients',[r['coefficient'] for r in source['quotient_coefficients']]==list(range(9)))
for row in source['quotient_coefficients']:Q+=native(S.sympify(row['polynomial']))*T**row['coefficient']
yes('residual quotient monic',native(S.sympify(source['quotient_coefficients'][8]['polynomial']))==1)
remainder=R-x**4*Q
yes('literal source quotient transport',remainder==sum((x_coefficient(R,j)*x**j for j in [2,3]),ctx.constant(0)))


def trim(q):
 q=list(q)
 while q and q[-1]==0:q.pop()
 return q

def peval(q,value,mod=None):
 out=0
 for c in reversed(q):out=out*value+c;out=out%mod if mod else out
 return out

def integer_bareiss(matrix):
 work=[r[:] for r in matrix];N=len(work);sign=1;previous=1
 for k in range(N-1):
  pivot=next((i for i in range(k,N) if work[i][k]),None)
  if pivot is None:return 0
  if pivot!=k:work[k],work[pivot]=work[pivot],work[k];sign=-sign
  value=work[k][k]
  for i in range(k+1,N):
   left=work[i][k]
   for j in range(k+1,N):
    numerator=work[i][j]*value-left*work[k][j];assert numerator%previous==0
    work[i][j]=numerator//previous
   work[i][k]=0
  previous=value
 return sign*work[-1][-1]

def sylvester(f,g):
 m,n=len(f)-1,len(g)-1;N=m+n;M=[[0]*N for _ in range(N)]
 for i in range(n):
  for j,c in enumerate(reversed(f)):M[i][i+j]=c
 for i in range(m):
  for j,c in enumerate(reversed(g)):M[n+i][i+j]=c
 return M

def add(f,g):
 out=[0]*max(len(f),len(g))
 for i,c in enumerate(f):out[i]+=c
 for i,c in enumerate(g):out[i]+=c
 return trim([c%prime for c in out])

def mul(f,g):
 if not f or not g:return []
 out=[0]*(len(f)+len(g)-1)
 for i,c in enumerate(f):
  for j,d in enumerate(g):out[i+j]+=c*d
 return trim([c%prime for c in out])

def rem(f,h):
 work=trim([c%prime for c in f]);assert h[-1]==1
 while len(work)>=len(h):
  shift=len(work)-len(h);value=work[-1]
  if value:
   for j,c in enumerate(h):work[shift+j]=(work[shift+j]-value*c)%prime
  work=trim(work)
 return work

def npoly(f):return nmod_poly([int(c)%prime for c in f],prime)

def field(h):
 degree=h.degree();coeff=list(map(int,h.coeffs()));yes('finite presentation leading coefficient '+str(degree),degree>0 and coeff[-1]%prime!=0)
 inverse=pow(coeff[-1]%prime,-1,prime);monic=[c*inverse%prime for c in coeff];assert len(monic)==degree+1
 return npoly(coeff),monic

def inverse(f,h,monic):
 g,u,v=f.xgcd(h);yes('modular inverse exists',g.degree()==0 and int(g[0])!=0)
 value=u*pow(int(g[0]),-1,prime)%h;plain=list(map(int,value.coeffs()))
 yes('plain polynomial inverse replay',rem(mul(list(map(int,f.coeffs())),plain),monic)==[1])
 return value

def det_mod_plain(M):
 work=[row[:] for row in M];N=len(work);det=1
 for k in range(N):
  pivot=next((i for i in range(k,N) if work[i][k]),None)
  if pivot is None:return 0
  if pivot!=k:work[k],work[pivot]=work[pivot],work[k];det=-det
  value=work[k][k];det=det*value%prime;inverse=pow(value,-1,prime);pivotrow=work[k]
  for i in range(k+1,N):
   row=work[i]
   if not row[k]:continue
   multiplier=row[k]*inverse%prime;row[k]=0
   for j in range(k+1,N):row[j]=(row[j]-multiplier*pivotrow[j])%prime
 return det%prime

# The complementary slope chart is an exact birational coordinate change.
# It is recorded as a verified working system, not a global exclusion.
d,v=S.symbols('d v');Delta=ss*(3*ss**3-2)
Tslope=2*(3*ss**3-1)*(aa*aa+aa)+6*ss**6+4*ss**5-2*ss**3-3*ss*ss
yes('exact local slope coordinate',x_coefficient(acrit,1)==Tslope+2*Delta*bb)
def b_coefficient(poly,j):return ctx.from_dict({ex[:3]+(0,)+ex[4:]:int(c) for ex,c in poly.to_dict().items() if ex[3]==j})
for index,coordinate in enumerate(X):
 yes('critical targets linear in b '+str(index),all(ex[3]<=1 for ex in coordinate.to_dict()))
 numerator=2*Delta*b_coefficient(coordinate,0)+b_coefficient(coordinate,1)*(xi-Tslope)
 yes('slope targets at most quadratic in a '+str(index),all(ex[2]<=2 for ex in numerator.to_dict()))
raw=json.loads((HERE/'residual-triple-input.json').read_text());slope=json.loads((HERE/'slope-triple-input.json').read_text())
yes('five raw and transformed triple constraints',len(raw['equations'])==len(slope['equations'])==5)
raw_native=[]
for index,(row,expected) in enumerate(zip(raw['equations'],f+[Q,Q.derivative('x'),Q.derivative('x').derivative('x')])):
 original=ctx.from_dict({(ex[0],ex[3],ex[2],ex[1],0):int(c) for ex,c in S.Poly(S.sympify(row['polynomial']),d,b,a,s,domain=S.ZZ).terms()});terms=original.to_dict();ex=next(iter(terms));ratio=Fraction(int(expected.to_dict().get(ex,0)),int(terms[ex]))
 yes('raw residual-triple source transport '+str(index),ratio!=0 and ratio.denominator*expected==ratio.numerator*original);raw_native.append(original)
for index,(original,row) in enumerate(zip(raw_native,slope['equations'])):
 yes('raw source quadratic in b '+str(index),all(ex[3]<=2 for ex in original.to_dict()))
 homogeneous=sum((b_coefficient(original,j)*(xi-Tslope)**j*(2*Delta)**(2-j) for j in range(3)),ctx.constant(0))
 transformed=ctx.from_dict({(ex[0],ex[3],ex[2],0,ex[1]):int(c) for ex,c in S.Poly(S.sympify(row['polynomial']),d,v,a,s,domain=S.ZZ).terms()})
 quotient,remainder=divmod(homogeneous,transformed)
 yes('literal slope source divisibility '+str(index),remainder==0 and bool(quotient))
 yes('only s-dependent slope factors '+str(index),all(ex[0]==ex[2]==ex[3]==ex[4]==0 for ex in quotient.to_dict()))
 coefficients=[0]*(max(ex[1] for ex in quotient.to_dict())+1)
 for ex,c in quotient.to_dict().items():coefficients[ex[1]]=int(c)
 factor=fmpz_poly(coefficients)
 for allowed in [fmpz_poly([0,1]),fmpz_poly([-2,0,0,3])]:
  while True:
   q,r=divmod(factor,allowed)
   if r:break
   factor=q
 yes('no undeclared slope denominator factor '+str(index),factor.degree()==0 and factor[0]!=0)
# The whole denominator-zero carrier is retained, including all four complex
# marked-root values, rather than specializing only a convenient one.
carrier=3*s**4-2*s
yes('specified exceptional carrier',S.expand(S.sympify(exc['carrier'])-carrier)==0)
yes('source exceptional variable order',exc['variables']==['s','a','b'])
yes('two specialized source polynomials',len(exc['equations'])==2)
F=[]
for j,(original,text) in enumerate(zip(fexpr,exc['equations'])):
 remainder=S.rem(S.Poly(original,s),S.Poly(carrier,s)).as_expr();expression=S.sympify(text);ratio=S.cancel(remainder/expression)
 yes('exceptional source remainder '+str(j),ratio.is_Rational and ratio!=0);q=S.Poly(expression,b);yes('exceptional equation linear in b '+str(j),q.degree()==1);F.append(q)
raw=F[0].nth(1)*F[1].nth(0)-F[1].nth(1)*F[0].nth(0)
residual=S.rem(S.Poly(S.expand(raw),s),S.Poly(carrier,s)).as_expr();H=S.sympify(exc['linear_resultant']);ratio=S.cancel(residual/H)
yes('literal linear-resultant transport',ratio.is_Rational and ratio!=0)
yes('a-degree of the exceptional source',S.Poly(H,a).degree()==6)
Ht=S.Poly(S.expand(H.subs(a,z-s)),s,z,domain=S.ZZ);Dc=S.Poly(carrier,s,z,domain=S.ZZ)
def coefficients_in_s(poly):
 rows=[[0]*(poly.degree(z)+1) for _ in range(poly.degree(s)+1)]
 for (i,j),c in poly.terms():rows[i][j]=int(c)
 return rows
fc,gc=map(coefficients_in_s,[Dc,Ht]);m,n=len(fc)-1,len(gc)-1;N=m+n
row=cert['cofactor_row'];yes('specified Sylvester matrix and cofactor row',N==cert['sylvester_size']==13 and row==n==9)
U=fmpz_poly(cert['resultant_coefficients']);P=fmpz_poly(cert['projection_coefficients']);P23=fmpz_poly(cert['nonorigin_projection']);num=fmpz_poly(cert['s_numerator']);den=fmpz_poly(cert['s_denominator'])
bound=n*Dc.degree(z)+m*Ht.degree(z);co_bound=n*Dc.degree(z)+(m-1)*Ht.degree(z)
yes('exact determinant degree bounds',bound==24 and co_bound==18)
yes('stored polynomial degrees',U.degree()<=bound and num.degree()<=co_bound and den.degree()<=co_bound)
for value in range(-12,13):
 matrix=sylvester([peval(q,value) for q in fc],[peval(q,value) for q in gc]);expected=int(U(value))
 yes('FLINT integer projection determinant '+str(value),int(fmpz_mat(matrix).det())==expected)
 yes('plain integer projection determinant '+str(value),integer_bareiss(matrix)==expected)
for value in range(-9,10):
 matrix=sylvester([peval(q,value) for q in fc],[peval(q,value) for q in gc])
 for column,polynomial in [(N-2,num),(N-1,den)]:
  minor=[[c for j,c in enumerate(line) if j!=column] for i,line in enumerate(matrix) if i!=row];sign=(-1)**(row+column);expected=int(polynomial(value))
  yes('FLINT integer cofactor '+str([value,column]),sign*int(fmpz_mat(minor).det())==expected)
  yes('plain integer cofactor '+str([value,column]),sign*integer_bareiss(minor)==expected)
k=next(i for i,c in enumerate(P.coeffs()) if c);ratio=Fraction(int(U[k]),int(P[k]));yes('literal primitive projection normalization',ratio!=0 and ratio.denominator*U==ratio.numerator*P)
yes('actual projection degrees and origin split',P.degree()==24 and P==fmpz_poly([0,1])*P23 and P23.degree()==23 and P23[0]!=0)
# Every retained root is an actual point, proved over Z rather than inferred
# from a modular necessary equation or a bounded point search.
yes('both forward source identities',[r['name'] for r in cert['forward_identities']]==['carrier','linear_resultant'])
for record,polynomial in zip(cert['forward_identities'],[Dc,Ht]):
 degree=polynomial.degree(s);Npow=[fmpz_poly([1])];Dpow=[fmpz_poly([1])]
 for _ in range(degree):Npow.append(Npow[-1]*num);Dpow.append(Dpow[-1]*den)
 homogeneous=fmpz_poly([])
 for j in range(degree+1):
  coefficient=S.Poly(S.Poly(polynomial.as_expr(),s).nth(j),z);value=fmpz_poly([int(coefficient.nth(k)) for k in range(coefficient.degree()+1)]) if not coefficient.is_zero else fmpz_poly([])
  homogeneous+=value*Npow[j]*Dpow[degree-j]
 yes('full homogeneous forward identity '+record['name'],list(map(int,homogeneous.coeffs()))==record['homogeneous_coefficients'])
 quotient=fmpz_poly(record['quotient_coefficients']);yes('integer forward source divisibility '+record['name'],homogeneous==quotient*P)
 left=list(map(int,quotient.coeffs()));right=list(map(int,P.coeffs()));product=[0]*(len(left)+len(right)-1)
 for i,c in enumerate(left):
  for j,d in enumerate(right):product[i+j]+=c*d
 yes('plain integer forward coefficient product '+record['name'],trim(product)==list(map(int,homogeneous.coeffs())))
announce('full source, projection and positive coverage passed')
# All modular unit tests are in fixed finite FREE algebras over Z_(prime).
prime=cert['prime'];yes('fixed prime',prime==1000000007)
yes('prime by complete trial division',all(prime%d for d in range(2,math.isqrt(prime)+1)))
norm_records={row['name']:row for row in cert['norms']};yes('five specified norm witnesses',set(norm_records)=={'s_denominator','projection_squarefree','b_denominator','residual_squarefree','marked_root_nonorigin'})
verified_norms=[]
def unit_norm(name,h,monic,value):
 d=h.degree();current=value%h;columns=[];X=npoly([0,1])
 for _ in range(d):columns.append(list(map(int,current.coeffs())));current=current*X%h
 matrix=[[columns[j][i] if i<len(columns[j]) else 0 for j in range(d)] for i in range(d)]
 plain=[];current=rem(list(map(int,value.coeffs())),monic)
 for _ in range(d):plain.append(current+[0]*(d-len(current)));current=rem([0]+current,monic)
 yes('plain multiplication matrix '+name,matrix==[list(r) for r in zip(*plain)])
 det=int(nmod_mat(matrix,prime).det());other=det_mod_plain(matrix);expected=norm_records[name]
 yes('nonzero FLINT finite free norm '+name,det!=0 and det==expected['det_mod'] and d==expected['dimension'])
 yes('plain finite free norm '+name,other==det);verified_norms.append(expected)
def evaluate(expression,h,monic,values,variables):
 terms={ex:int(c)%prime for ex,c in S.Poly(expression,*variables,domain=S.ZZ).terms()}
 def recurse(data,level,native_mode):
  if not data:return npoly([]) if native_mode else []
  answer=npoly([]) if native_mode else [];value=values[level] if native_mode else list(map(int,values[level].coeffs()))
  for j in range(max(ex[level] for ex in data),-1,-1):
   sub={ex[:level]:c for ex,c in data.items() if ex[level]==j}
   if level==0:
    scalar=sum(sub.values())%prime;q=npoly([scalar]) if native_mode else trim([scalar])
   else:q=recurse(sub,level-1,native_mode)
   answer=(answer*value+q)%h if native_mode else rem(add(mul(answer,value),q),monic)
  return answer
 native_value=recurse(terms,len(values)-1,True);plain_value=recurse(terms,len(values)-1,False)
 yes('plain polynomial source evaluation',list(map(int,native_value.coeffs()))==plain_value)
 return native_value
hp,mp=field(P);unit_norm('s_denominator',hp,mp,npoly(den.coeffs()));unit_norm('projection_squarefree',hp,mp,hp.derivative());Z=npoly([0,1]);S0=npoly(num.coeffs())*inverse(npoly(den.coeffs()),hp,mp)%hp;A0=(Z-S0)%hp
Blead=evaluate(F[0].nth(1),hp,mp,[S0,A0],[s,a]);Bconst=evaluate(F[0].nth(0),hp,mp,[S0,A0],[s,a]);unit_norm('b_denominator',hp,mp,Blead);B0=-Bconst*inverse(Blead,hp,mp)%hp
quotient=[evaluate(S.sympify(r['polynomial']),hp,mp,[S0,A0,B0],[s,a,b]) for r in source['quotient_coefficients']]
yes('finite residual coefficients match witness',[list(map(int,q.coeffs())) for q in quotient]==cert['quotient_coefficients_mod'])
yes('monic finite residual octic',quotient[-1]==1)
qc=[list(map(int,q.coeffs())) for q in quotient];derivative=[[((j+1)*c)%prime for c in qc[j+1]] for j in range(8)]
finite_bound=15*max(len(q)-1 for q in qc);delta=cert['quotient_resultant_mod'];yes('finite resultant bound',finite_bound==345 and delta==trim(delta) and 0<len(delta)<=finite_bound+1)
yes('enough distinct finite interpolation points',finite_bound+1<prime)
for value in range(finite_bound+1):
 matrix=sylvester([peval(q,value,prime) for q in qc],[peval(q,value,prime) for q in derivative]);yes('nominal derivative Sylvester size '+str(value),len(matrix)==15);expected=peval(delta,value,prime)
 yes('FLINT finite residual determinant '+str(value),int(nmod_mat(matrix,prime).det())==expected)
 yes('plain finite residual determinant '+str(value),det_mod_plain(matrix)==expected)
unit_norm('residual_squarefree',hp,mp,npoly(delta));marked=npoly([])
for q in reversed(quotient):marked=(marked*S0+q)%hp
h23,m23=field(P23);unit_norm('marked_root_nonorigin',h23,m23,marked)
# The projection root z=0 is the known exact fivefold origin.
yes('origin s recovery',num[0]==0 and den[0]!=0)
yes('origin b recovery',F[0].nth(0).subs({s:0,a:0})==0 and F[0].nth(1).subs({s:0,a:0})!=0)
t=S.symbols('t');H7=t**7-2*t**4+2*t-1;Qzero=sum(S.sympify(r['polynomial']).subs({s:0,a:0,b:0})*t**r['coefficient'] for r in source['quotient_coefficients'])
yes('exact origin octic',S.expand(Qzero-t*H7)==0)
yes('origin residual simple and exactly fivefold',H7.subs(t,0)==-1 and S.discriminant(H7,t)==713)

# Global structure of the quadruple-root incidence. A good fibre proves
# nonzero GENERIC norms over Q(s), not units at every complex s-value.
generic=json.loads((HERE/'generic-certificate.json').read_text());old_curve=json.loads((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_text())
fb=[S.Poly(q,b) for q in fexpr];Alead,Bmiddle,Cconstant=[fb[0].nth(j) for j in [2,1,0]]
yes('quadratic source leading coefficient',S.expand(Alead-4*s*s*(3*s**3-2)**2)==0)
L=S.Poly(S.sympify(old_curve['linear_equations'][0]),b)
yes('generic source linear cancellation',native(4*s*(3*s**3-2))*native(L.as_expr())==native(Alead)*f[1]-native(fb[1].nth(2))*f[0])
P1,Q1=L.nth(1),L.nth(0);E1=S.expand(Alead*Q1*Q1-Bmiddle*P1*Q1+Cconstant*P1*P1)
yes('generic quadratic norm identity',native(E1)==native(P1*P1)*f[0]+native(L.as_expr())*native(Alead*Q1-Alead*P1*b-Bmiddle*P1))
curve_row=old_curve['equations'][3];curve=S.Poly(S.sympify(curve_row['polynomial']),a)
factor=S.Integer(1)
for term in curve_row['removed_exception_factors']:
 yes('only declared generic curve normalization',term in ['s','3*s**3 - 2']);factor*=S.sympify(term)
ratio=S.cancel(E1/(curve.as_expr()*factor));yes('exact generic plane curve transport',ratio.is_Rational and ratio!=0)
yes('generic plane curve degree',curve.degree()==8)
yes('exact generic leading coefficient',S.expand(curve.LC()-2560*s**4*(3*s**3-2)**2*(3*s**3-1)**2*(3*s**3+1))==0)
coefficients=[]
for j in range(9):
 q=S.Poly(curve.nth(j),s);coefficients.append(fmpz_poly([int(q.nth(k)) for k in range(q.degree()+1)]) if not q.is_zero else fmpz_poly([]))
first=coefficients[0];second=sum(coefficients[1:],fmpz_poly([]));yes('vertical coefficient witness source',list(map(int,first.coeffs()))==generic['coefficient_first'] and list(map(int,second.coeffs()))==generic['coefficient_second'])
vertical=sylvester(list(map(int,first.coeffs())),list(map(int,second.coeffs())));yes('full actual-degree coefficient Sylvester matrix',[len(vertical),len(vertical)]==generic['coefficient_sylvester_shape']==[54,54])
vertical_mod=[[c%prime for c in row] for row in vertical];det=int(nmod_mat(vertical_mod,prime).det());yes('nonzero integer coefficient resultant via FLINT',det==generic['coefficient_det_mod']!=0)
yes('plain coefficient resultant replay',det_mod_plain(vertical_mod)==det)
# Thus at no complex s are all a-coefficients zero. On D0!=0, at most
# eight a-values and two b-values each occur. D0=0 was fully checked above.
point=generic['sample_s'];yes('specified good rational root parameter',point==1 and point*(3*point**3-2)!=0)
P8=S.Poly(curve.as_expr().subs(s,point),a,domain=S.ZZ);sample_coefficients=[int(P8.nth(j)) for j in range(P8.degree()+1)]
yes('actual specialized degree8 polynomial',P8.degree()==8 and sample_coefficients==generic['sample_curve_coefficients'])
sample_poly=fmpz_poly(sample_coefficients);ha,ma=field(sample_poly)
for row in generic['norms']:
 yes('fresh sample norm name',row['name'] not in norm_records);norm_records[row['name']]=row
expected_sample_names={'sample_curve_squarefree','sample_b_denominator','sample_residual_squarefree','sample_marked_root'}
yes('four specified good-fibre norms',{r['name'] for r in generic['norms']}==expected_sample_names)
unit_norm('sample_curve_squarefree',ha,ma,ha.derivative())
Sa=npoly([point]);Aa=npoly([0,1]);Pa=evaluate(P1,ha,ma,[Sa,Aa],[s,a]);Qa=evaluate(Q1,ha,ma,[Sa,Aa],[s,a]);unit_norm('sample_b_denominator',ha,ma,Pa);Ba=-Qa*inverse(Pa,ha,ma)%ha
sample_quotient=[evaluate(S.sympify(row['polynomial']),ha,ma,[Sa,Aa,Ba],[s,a,b]) for row in source['quotient_coefficients']]
yes('good-fibre residual coefficients match',[list(map(int,q.coeffs())) for q in sample_quotient]==generic['quotient_coefficients_mod'])
yes('good-fibre residual monic',sample_quotient[-1]==1)
sample_q=[list(map(int,q.coeffs())) for q in sample_quotient];sample_derivative=[[((j+1)*c)%prime for c in sample_q[j+1]] for j in range(8)]
sample_bound=15*max(len(q)-1 for q in sample_q);sample_delta=generic['quotient_resultant_mod'];yes('good-fibre resultant degree bound',sample_bound==105 and sample_delta==trim(sample_delta) and 0<len(sample_delta)<=sample_bound+1)
for value in range(sample_bound+1):
 matrix=sylvester([peval(q,value,prime) for q in sample_q],[peval(q,value,prime) for q in sample_derivative]);assert len(matrix)==15;expected=peval(sample_delta,value,prime)
 yes('FLINT good-fibre determinant '+str(value),int(nmod_mat(matrix,prime).det())==expected)
 yes('plain good-fibre determinant '+str(value),det_mod_plain(matrix)==expected)
unit_norm('sample_residual_squarefree',ha,ma,npoly(sample_delta));sample_marked=npoly([])
for q in reversed(sample_quotient):sample_marked=(sample_marked*Sa+q)%ha
unit_norm('sample_marked_root',ha,ma,sample_marked)
receipt={'result':'PASS','claim_id':'GAOF7FOUR1','checks':len(checks),'scope':'For Gao v2 F7 on c!=0, quadruple-root incidence has finite fibres over every marked root parameter s, with at most16 points on D0=s*(3s^3-2)!=0 and exactly6 at each D0=0 value. It has generic degree8 and all generic normalized targets have one quadruple root and eight other simple roots, hence8 source points. Its nonordinary normalized-target locus is finite; in particular the still-undecided4+3+2 three-point stratum is finite. The full D0=0 locus is exactly24 normalized targets: one known fivefold origin with7 source points and23 quadruple targets with8 source points. No global4+3+2 exclusion or finite exceptional count is claimed. The3+2+2+2 case, global three-point existence, full strata,F5 and JC2 remain open.','generic_degree':8,'all_marked_parameter_fibres_finite':True,'nonordinary_normalized_locus':'FINITE; explicit complete list not computed','uniform_fibre_bound_on_D0_nonzero':16,'pole_chart_target_count':24,'pole_chart_fibre_counts':{'7':1,'8':23},'good_sample_s':point,'good_sample_target_count':8,'coefficient_coprimality_shape':[54,54],'coefficient_coprimality_det_mod':det,'integer_resultant_evaluation_points':25,'integer_cofactor_evaluation_points':19,'integer_determinant_engines':['FLINT','plain Python Bareiss'],'finite_free_norms':verified_norms,'pole_resultant_degree_bound':finite_bound,'pole_resultant_evaluation_points':finite_bound+1,'sample_resultant_degree_bound':sample_bound,'sample_resultant_evaluation_points':sample_bound+1,'finite_determinant_engines':['FLINT','plain Python Gaussian'],'finite_source_engines':['FLINT univariate polynomials','plain Python coefficient lists'],'positive_coverage':'Exact homogeneous identities on the pole chart; exact generic quadratic-linear source identities with unit denominators and squarefree specialized degree8 presentation.','slope_chart':'Exact derivative coordinate and all five transformed generic equations replayed, with only s and3s^3-2 factors removed. No generic triple-root exclusion.','checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'certificate_sha256':hashlib.sha256((HERE/'exceptional-certificate.json').read_bytes()).hexdigest(),'generic_certificate_sha256':hashlib.sha256((HERE/'generic-certificate.json').read_bytes()).hexdigest(),'source_sha256':hashlib.sha256((HERE/'quadruple-input.json').read_bytes()).hexdigest(),'exceptional_source_sha256':hashlib.sha256((HERE/'exceptional-source.json').read_bytes()).hexdigest(),'raw_triple_sha256':hashlib.sha256((HERE/'residual-triple-input.json').read_bytes()).hexdigest(),'slope_triple_sha256':hashlib.sha256((HERE/'slope-triple-input.json').read_bytes()).hexdigest(),'generic_curve_source_sha256':hashlib.sha256((HERE.parent/'gao-f7-linear-cube/high-root-bivariate-input.json').read_bytes()).hexdigest(),'python':platform.python_version(),'sympy':S.__version__,'python_flint':flint.__version__,'seconds':time.monotonic()-START,'assurance':'Exact arithmetic replay plus written source, adjugate, finite-free unit and generic-specialization arguments. A good specialization proves generic nonvanishing, not a global unit at every s. No independent end-to-end review, formal verification, external review or novelty priority.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
