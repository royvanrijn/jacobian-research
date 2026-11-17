"""Fail-closed sixfold-root exclusion by source transport, determinants and identities.

Runs no Groebner basis, multivariate resultant algorithm or reconstruction.
Every resultant is verified as a polynomial determinant by exact interpolation.
"""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,platform,time
import sympy as S
import flint
from flint import fmpz,fmpz_mat,fmpz_mpoly_ctx,nmod_mat
if not __debug__:raise RuntimeError('Run without -O')
HERE=Path(__file__).resolve().parent;START=time.monotonic();checks=[]
def yes(name,condition):
 assert condition,name
 checks.append(name)
def equal(name,expr):yes(name,S.expand(expr)==0)
def primitive_ratio(name,f,g):
 ratio=S.cancel(f/g)
 yes(name,ratio.is_Rational and ratio!=0)
 return ratio
s,a,b,t,xi=S.symbols('s a b t xi')
Y1,Y2,Y3,Y4=S.symbols('Y1 Y2 Y3 Y4');YY=[Y1,Y2,Y3,Y4]
bt=Y4-3*Y2*t*t+2*Y1*t**3-t**4
at=Y3-2*Y2*t+Y1*t*t-t**3-t*t*bt
R=S.expand(at*at+at+t*bt+t**3+Y1*t-Y2)
G2=a*a+a+s*b+s**3;G3=2*s*G2+a+s**3+s*s*b;G4=3*s*s*G2+b+s**4
X1=S.expand(-S.det(S.Matrix([G2,G3,G4]).jacobian([s,a,b])))
X=[X1,G2+s*X1,G3+s*s*X1,G4+s**3*X1]
full=[xi,G2+s*xi,G3+s*s*xi,G4+s**3*xi]
full_sub=dict(zip(YY,full));critical_sub=dict(zip(YY,X))
equal('root-incidence a reconstruction',at.subs(t,s).subs(full_sub)-a)
equal('root-incidence b reconstruction',bt.subs(t,s).subs(full_sub)-b)
equal('root-incidence equation',R.subs(t,s).subs(full_sub))
equal('multiple-root coverage without variable division',S.diff(R,t).subs(t,s).subs(full_sub)-(xi-X1))
for j in [0,1]:equal('critical derivative '+str(j),S.diff(R,t,j).subs(t,s).subs(critical_sub))
old=json.loads((HERE.parent/'gao-f7-fibres/sixfold-input.json').read_text())
yes('retained sixfold source variable order',old['variables']==['s','a','b'])
yes('exactly four sixfold derivative constraints',len(old['polynomials'])==4)
f=[]
for j,text in enumerate(old['polynomials'],2):
 expression=S.expand(S.diff(R,t,j).subs(t,s).subs(critical_sub)/S.factorial(j))
 source=S.sympify(text);primitive_ratio('fresh necessary derivative '+str(j),expression,source)
 f.append(S.Poly(source,b))
yes('all necessary equations quadratic in b',all(q.degree()==2 for q in f))
A,B,C=[f[0].nth(i) for i in [2,1,0]]
equal('common leading quadratic coefficient',A-4*s*s*(3*s**3-2)**2)
contents=[4*s*(3*s**3-2),S.Integer(4),4*s*s]
source=json.loads((HERE/'high-root-bivariate-input.json').read_text())
yes('bivariate variable order',source['variables']==['a','s'])
yes('three linear equations retained',len(source['linear_equations'])==3)
linear=[]
for j,text in enumerate(source['linear_equations'],1):
 L=S.Poly(S.sympify(text),b);yes('linear b equation '+str(j),L.degree()==1)
 equal('exact quadratic cancellation '+str(j),contents[j-1]*L.as_expr()-(A*f[j].as_expr()-f[j].nth(2)*f[0].as_expr()))
 linear.append(L)
rows=source['equations'];yes('six bivariate necessary constraints',len(rows)==6)
for index,row in enumerate(rows):
 if row['type']=='linear_cross':
  i,j=row['indices'];yes('cross indices '+str(index),0<=i<j<3)
  expression=S.expand(linear[i].nth(1)*linear[j].nth(0)-linear[j].nth(1)*linear[i].nth(0))
  equal('cross identity '+str(index),expression-(linear[i].nth(1)*linear[j].as_expr()-linear[j].nth(1)*linear[i].as_expr()))
 else:
  yes('quadratic resultant type '+str(index),row['type']=='quadratic_resultant')
  i,=row['indices'];yes('quadratic indices '+str(index),0<=i<3)
  P,Q=linear[i].nth(1),linear[i].nth(0)
  expression=S.expand(A*Q*Q-B*P*Q+C*P*P)
  equal('necessary quadratic identity '+str(index),expression-P*P*f[0].as_expr()-linear[i].as_expr()*(A*Q-A*P*b-B*P))
 factor=S.Integer(1)
 for term in row['removed_exception_factors']:
  yes('only retained exceptional factors '+str(index),term in ['s','3*s**3 - 2'])
  factor*=S.sympify(term)
 primitive_ratio('bivariate coefficient transport '+str(index),expression,S.sympify(row['polynomial'])*factor)
# Both exceptional loci are checked below; none is silently discarded.
certificate=json.loads((HERE/'high-root-certificate.json').read_text())
yes('three specified resultant pairs',[r['pair'] for r in certificate['resultants']]==[[0,1],[0,3],[0,2]])
def trim(coefficients):
 coefficients=list(coefficients)
 while coefficients and coefficients[-1]==0:coefficients.pop()
 return coefficients
def evaluate(coefficients,x):
 result=0
 for c in reversed(coefficients):result=result*x+c
 return result
def multiply(f,g):
 result=[0]*(len(f)+len(g)-1)
 for i,a0 in enumerate(f):
  for j,b0 in enumerate(g):result[i+j]+=a0*b0
 return trim(result)
def power(f,n):
 result=[1]
 for _ in range(n):result=multiply(result,f)
 return result
def bareiss(matrix):
 work=[row[:] for row in matrix];N=len(work);sign=1;previous=1
 for k in range(N-1):
  pivot=next((i for i in range(k,N) if work[i][k]),None)
  if pivot is None:return 0
  if pivot!=k:work[k],work[pivot]=work[pivot],work[k];sign=-sign
  value=work[k][k]
  for i in range(k+1,N):
   left=work[i][k]
   for j in range(k+1,N):
    numerator=work[i][j]*value-left*work[k][j]
    assert numerator%previous==0
    work[i][j]=numerator//previous
   work[i][k]=0
  previous=value
 return sign*work[-1][-1]
def sylvester(f,g):
 m,n=len(f)-1,len(g)-1;N=m+n;matrix=[[0]*N for _ in range(N)]
 for i in range(n):
  for j,c in enumerate(reversed(f)):matrix[i][i+j]=c
 for i in range(m):
  for j,c in enumerate(reversed(g)):matrix[n+i][i+j]=c
 return matrix
def coefficients_in_a(poly):
 q=S.Poly(S.sympify(poly),a,s,domain=S.ZZ);m=q.degree(a);N=q.degree(s)
 result=[[0]*(N+1) for _ in range(m+1)]
 for (i,j),c in q.terms():result[i][j]=int(c)
 return result,m,N
reduced=[];point_counts=[];resultant_bounds=[]
for record in certificate['resultants']:
 i,j=record['pair'];fc,m,dsf=coefficients_in_a(rows[i]['polynomial']);gc,n,dsg=coefficients_in_a(rows[j]['polynomial'])
 bound=n*dsf+m*dsg
 coefficients=record['resultant_coefficients']
 yes('integer resultant coefficients '+str([i,j]),all(type(c) is int for c in coefficients))
 yes('resultant degree bounded '+str([i,j]),coefficients==trim(coefficients) and 0<len(coefficients)<=bound+1)
 start=-(bound//2)
 for x in range(start,start+bound+1):
  fv=[evaluate(q,x) for q in fc];gv=[evaluate(q,x) for q in gc]
  # Preserve nominal a-degrees even when a leading coefficient vanishes here.
  matrix=sylvester(fv,gv);assert len(matrix)==m+n
  expected=evaluate(coefficients,x)
  yes('FLINT exact resultant evaluation '+str([i,j,x]),int(fmpz_mat(matrix).det())==expected)
  yes('plain integer Bareiss evaluation '+str([i,j,x]),bareiss(matrix)==expected)
 q=record['primitive_reduced_coefficients'];yes('primitive resultant integer coefficients '+str([i,j]),all(type(c) is int for c in q) and q==trim(q) and len(q)>0)
 e=record['removed_factors'];yes('only two exceptional resultant factors '+str([i,j]),set(e)=={'s','3s^3-2'} and all(type(k) is int and 0<=k<=1000 for k in e.values()))
 product=multiply(power([0,1],e['s']),multiply(power([-2,0,0,3],e['3s^3-2']),q))
 yes('same resultant degree after factor transport '+str([i,j]),len(product)==len(coefficients))
 k=next(k for k,c in enumerate(product) if c)
 scalar=Fraction(coefficients[k],product[k]);yes('nonzero constant normalization '+str([i,j]),scalar!=0)
 yes('literal univariate factor transport '+str([i,j]),all(Fraction(c)==scalar*d for c,d in zip(coefficients,product)))
 reduced.append(q);point_counts.append(bound+1);resultant_bounds.append(bound)
 print(json.dumps({'stage':'resultant replay','pair':[i,j],'degree_bound':bound,'exact_integer_evaluations':bound+1,'seconds':time.monotonic()-START}),flush=True)
# A nonzero modular determinant of this ACTUAL-DEGREE integer Sylvester
# matrix proves a nonzero integer resultant, hence coprimality over C.
yes('fixed combination',certificate['combination']==[1,1])
g=[0]*max(len(reduced[1]),len(reduced[2]))
for q in reduced[1:]:
 for i,c in enumerate(q):g[i]+=c
g=trim(g);matrix=sylvester(reduced[0],g);N=len(matrix)
yes('full actual Sylvester dimensions',certificate['coprimality_sylvester_shape']==[N,N]==[623,623])
prime=certificate['prime'];yes('fixed prime',prime==1000000007)
yes('prime by complete trial division',all(prime%d for d in range(2,math.isqrt(prime)+1)))
modular=[[c%prime for c in row] for row in matrix]
det=int(nmod_mat(modular,prime).det())
yes('nonzero complete Sylvester determinant',det==certificate['det_mod']!=0)
other=1
for k in range(N):
 pivot=next((i for i in range(k,N) if modular[i][k]),None);assert pivot is not None
 if pivot!=k:modular[k],modular[pivot]=modular[pivot],modular[k];other=-other
 value=modular[k][k];other=other*value%prime;inverse=pow(value,-1,prime);pivotrow=modular[k]
 for i in range(k+1,N):
  row=modular[i]
  if not row[k]:continue
  multiplier=row[k]*inverse%prime;row[k]=0
  for j in range(k+1,N):row[j]=(row[j]-multiplier*pivotrow[j])%prime
yes('plain Python modular Sylvester determinant',other%prime==det)
# Literal exceptional-case unit identities, in two exact integer engines.
yes('two exceptional identities',len(certificate['exceptional_identities'])==2)
exception_terms=[]
for case_index,identity in enumerate(certificate['exceptional_identities']):
 name=['high-root-szero-integral-input.json','high-root-cubic-exception-integral-input.json'][case_index]
 yes('specified exceptional source '+str(case_index),identity['source']==name)
 case=json.loads((HERE/name).read_text());V=[a,b] if case_index==0 else [a,b,s]
 yes('exceptional variables '+str(case_index),case['variables']==list(map(str,V)))
 expected=[q.as_expr().subs(s,0) for q in f] if case_index==0 else [S.rem(S.Poly(q.as_expr(),s),S.Poly(3*s**3-2,s)).as_expr() for q in f]+[3*s**3-2]
 yes('exceptional equation count '+str(case_index),len(case['equations'])==len(expected)==len(identity['multipliers']))
 originals=[]
 for j,(row,e) in enumerate(zip(case['equations'],expected)):
  expression=S.sympify(row['polynomial']);primitive_ratio('exceptional source transport '+str([case_index,j]),e,expression)
  originals.append({ex:int(c) for ex,c in S.Poly(expression,*V,domain=S.ZZ).terms()})
 denominator=1;multipliers=[]
 for terms in identity['multipliers']:
  coeffs={}
  for term in terms:
   ex=tuple(term['exponents']);assert len(ex)==len(V) and all(type(k) is int and k>=0 for k in ex) and ex not in coeffs
   value=Fraction(term['coefficient']);assert value!=0
   coeffs[ex]=value;denominator=math.lcm(denominator,value.denominator)
  multipliers.append(coeffs)
 ctx=fmpz_mpoly_ctx.get(list(map(str,V)));native=ctx.constant(0);plain={}
 for original,multiplier in zip(originals,multipliers):
  integer={ex:c.numerator*(denominator//c.denominator) for ex,c in multiplier.items()}
  native+=ctx.from_dict(integer)*ctx.from_dict(original)
  for ex,c in integer.items():
   for ey,d in original.items():
    ez=tuple(x+y for x,y in zip(ex,ey));plain[ez]=plain.get(ez,0)+c*d
 yes('FLINT exceptional integer-cleared unit identity '+str(case_index),native==denominator)
 yes('plain integer exceptional unit identity '+str(case_index),{ex:c for ex,c in plain.items() if c}=={(0,)*len(V):denominator})
 exception_terms.append(sum(map(len,multipliers)))
# Sharpness is already attained by the exact rational target Y=0.
# The exploratory finite-field seed was the reduction of this known witness;
# no Hensel or priority claim is needed.
witness=json.loads((HERE/'fivefold-witness.json').read_text())
yes('exact sharpness target',witness['target_Y']==[0,0,0,0] and witness['root']==0)
H=S.sympify(witness['quotient'])
equal('exact rational fivefold witness',R.subs(dict(zip(YY,[0,0,0,0])))-t**5*H)
yes('residual factor degree seven',S.Poly(H,t).degree()==7)
yes('residual factor nonzero at marked root',H.subs(t,0)==-1)
yes('exact residual discriminant',S.discriminant(H,t)==witness['quotient_discriminant']==713)
yes('residual squarefreeness by gcd',S.Poly(H,t).gcd(S.Poly(S.diff(H,t),t)).degree()==0)
yes('sharp witness has seven simple roots',witness['fibre_cardinality']==7 and witness['root_multiplicity']==5)
receipt={'result':'PASS','claim_id':'GAOF7MULT1','checks':len(checks),'scope':'Every Gao v2 F7 eliminant root has multiplicity at most5, on all root and target values. This is sharp at Y=0, where R=t^5*(t^7-2t^4+2t-1), with seven other simple roots. No leading-coefficient exceptional locus discarded. Combined with GAOF7CUBE1, four three-point nonsimple lists remain.','source_transport':'Fresh root-incidence reconstruction and critical-derivative coverage, followed by exact derivative2..5 reconstruction.','resultant_degrees':[len(r['resultant_coefficients'])-1 for r in certificate['resultants']],'resultant_degree_bounds':resultant_bounds,'exact_integer_evaluation_counts':point_counts,'resultant_engines':['FLINT integer determinant','plain Python Bareiss'],'coprimality_sylvester_shape':[N,N],'coprimality_det_mod':det,'prime':prime,'coprimality_engines':['FLINT modular determinant','plain Python integer Gaussian elimination'],'exceptional_multiplier_terms':exception_terms,'exceptional_engines':['FLINT integer convolution','plain Python integer convolution'],'additional_excluded_three_point_types':[[7,2]],'remaining_three_point_types':[[5,4],[5,2,2],[4,3,2],[3,2,2,2]],'python':platform.python_version(),'sympy':S.__version__,'python_flint':flint.__version__,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'certificate_sha256':hashlib.sha256((HERE/'high-root-certificate.json').read_bytes()).hexdigest(),'input_sha256':hashlib.sha256((HERE/'high-root-bivariate-input.json').read_bytes()).hexdigest(),'seconds':time.monotonic()-START,'assurance':'Exact certificate replay and written resultant/coverage argument. No formal verification, independent end-to-end review, external review or novelty priority. F7 three-point existence and JC2 remain open.'}
(HERE/'high-root-receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
