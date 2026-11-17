"""Fail-closed exact replay of the F7 triple-root linear pole chart.
No symbolic discriminant, resultant, gcd search, or Groebner search is run.
Fresh source identities and all norm polynomials are checked literally.
"""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,platform,time,sympy as S
import flint
from flint import fmpz_mpoly_ctx,fmpz_poly,fmpz_mat,nmod_mat
HERE=Path(__file__).resolve().parent;START=time.monotonic();checks=0
if not __debug__:raise RuntimeError('Run without -O')
def yes(c,name='check'):
 global checks
 assert c,name;checks+=1
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def emit(stage,**kw):print(json.dumps({'stage':stage,'seconds':time.monotonic()-START,**kw}),flush=True)
def bareiss(M):
 W=[row[:] for row in M];prev=1;sign=1
 for k in range(len(W)-1):
  pivot=next((i for i in range(k,len(W)) if W[i][k]),None)
  if pivot is None:return 0
  if pivot!=k:W[pivot],W[k]=W[k],W[pivot];sign=-sign
  value=W[k][k]
  for i in range(k+1,len(W)):
   left=W[i][k]
   for j in range(k+1,len(W)):
    q,r=divmod(W[i][j]*value-left*W[k][j],prev);assert r==0;W[i][j]=q
   W[i][k]=0
  prev=value
 return sign*W[-1][-1]
def sylv(f,g):
 m,n=len(f)-1,len(g)-1;M=[[0]*(m+n) for _ in range(m+n)]
 for i in range(n):
  for j,c in enumerate(reversed(f)):M[i][i+j]=c
 for i in range(m):
  for j,c in enumerate(reversed(g)):M[n+i][i+j]=c
 return M

def det_plain_mod(M,p):
 W=[row[:] for row in M];det=1
 for k in range(len(W)):
  pivot=next((i for i in range(k,len(W)) if W[i][k]),None)
  if pivot is None:return 0
  if pivot!=k:W[pivot],W[k]=W[k],W[pivot];det=-det
  v=W[k][k];det=det*v%p;vi=pow(v,-1,p)
  for i in range(k+1,len(W)):
   z=W[i][k]*vi%p
   if not z:continue
   W[i][k]=0
   for j in range(k+1,len(W)):W[i][j]=(W[i][j]-z*W[k][j])%p
 return det%p

# Rebuild the multiple-root incidence, then the degree9 residual.
source=json.loads((HERE.parent/'gao-f7-quadruple/quadruple-input.json').read_text());data=json.loads((HERE/'pole-source-candidate.json').read_text());norms=json.loads((HERE/'lower-order-norms-candidate.json').read_text());base=json.loads((HERE.parent/'gao-f7-quadruple/receipt.json').read_text())
yes(base['result']=='PASS' and base['claim_id']=='GAOF7FOUR1' and base['checker_sha256']==sha(HERE.parent/'gao-f7-quadruple/verify.py') and base['source_sha256']==sha(HERE.parent/'gao-f7-quadruple/quadruple-input.json'))
yes(data['source_sha256']==base['source_sha256'] and norms['source_sha256']==sha(HERE/'pole-source-candidate.json'))
s,a,b=S.symbols('s a b');ctx=fmpz_mpoly_ctx.get(['d','a','s','b']);d,an,sn,bn=ctx.gens()
def native(e):return ctx.from_dict({(0,ex[1],ex[0],ex[2]):int(c) for ex,c in S.Poly(e,s,a,b).terms()})
def coefficient(F,k,j):return ctx.from_dict({ex[:k]+(0,)+ex[k+1:]:int(c) for ex,c in F.to_dict().items() if ex[k]==j})
G2=a*a+a+s*b+s**3;G3=2*s*G2+a+s**3+s*s*b;G4=3*s*s*G2+b+s**4;X1=S.expand(-S.det(S.Matrix([G2,G3,G4]).jacobian([s,a,b])));g2,g3,g4,x1=map(native,[G2,G3,G4,X1]);X=[x1,g2+sn*x1,g3+sn**2*x1,g4+sn**3*x1];t=sn+d;Y1,Y2,Y3,Y4=X;bt=Y4-3*Y2*t*t+2*Y1*t**3-t**4;at=Y3-2*Y2*t+Y1*t*t-t**3-t*t*bt;R=at*at+at+t*bt+t**3+Y1*t-Y2
yes(coefficient(at,0,0)==an and coefficient(bt,0,0)==bn);yes(coefficient(R,0,0)==coefficient(R,0,1)==0)
f2=native(S.sympify(source['equations'][0]['polynomial']));c2=coefficient(R,0,2);ratio=Fraction(data['f2_ratio_in_Taylor_coefficient']);yes(ratio!=0 and ratio.denominator*c2==ratio.numerator*f2)
Q=ctx.from_dict({tuple(ex):int(c) for ex,c in data['Q9_terms']});yes(R==d**3*Q+d*d*c2);yes(coefficient(Q,0,9)==1 and Q.degrees()[3]==2)
A,B,C=[coefficient(f2,3,j) for j in [2,1,0]];D=sn*(3*sn**3-2);yes(A==4*D*D)
# The coefficient B is nonzero at every actual point on D=0.
Bb=S.Poly(S.sympify(source['equations'][0]['polynomial']),b).nth(1);yes(S.rem(S.Poly(Bb-2*(2*a+1)*(6*s**3-1),s,a),S.Poly(s*(3*s**3-2),s,a)).is_zero)
f2expr=S.sympify(source['equations'][0]['polynomial']);W=27-181*s-60*s*s-144*s**3;remainder=S.rem(S.Poly(108*f2expr.subs(a,-S.Rational(1,2)),s,b),S.Poly(s*(3*s**3-2),s,b)).as_expr();yes(S.expand(remainder-W)==0)
M=sylv([0,-2,0,0,3],[27,-181,-60,-144]);integer_det=int(fmpz_mat(M).det());yes(integer_det!=0 and bareiss(M)==integer_det)
# B=0 on this carrier means a=-1/2; the preceding nonzero determinant
# proves no source point exists there. Therefore b=-C/B is complete.
U=B*B*coefficient(Q,3,0)-B*C*coefficient(Q,3,1)+C*C*coefficient(Q,3,2);chart_models={}
for row in data['pole_charts']:
 h=sn if row['name']=='s-zero' else 3*sn**3-2
 F=ctx.from_dict({(ex[0],ex[1],ex[2],0):int(c) for ex,c in row['terms']});den=int(row['denominator_clearing']);content=int(row['integer_content_removed']);yes(den!=0 and content!=0)
 quotient,rem=divmod(den*U-content*F,h);yes(not rem,'exact pole residual source transport')
 lead=(2*an+1)**2*(1 if row['name']=='s-zero' else 6561);yes(coefficient(F,0,9)==lead);yes(F.degrees()[0]==9 and F.degrees()[1]<=8 and F.degrees()[2]<(1 if row['name']=='s-zero' else 3))
 chart_models[row['name']]=F
emit('source_and_full_pole_coverage',forbidden_coordinate_resultant=integer_det)

# Independent coefficient construction of the rectangular Sylvester maps.
def rectangular(F,order):
 rows=[];top=16-order;terms=F.to_dict();derivative=F.derivative('d')
 for polynomial,count in [(F,8-order),(derivative,9-order)]:
  for shift in range(count-1,-1,-1):
   row=[{} for _ in range(top+1)]
   for (i,j,k,z),c in polynomial.to_dict().items():
    col=top-int(i)-shift;yes(0<=col<=top);row[col][(int(j),int(k))]=int(c)
   rows.append(row)
 yes(len(rows)==17-2*order)
 return rows

def entry_value(entry,value):
 v=[0,0,0]
 for (j,k),c in entry.items():v[k]+=c*value**j
 return v

def matrix_at(M,columns,value,cubic):
 rows=[[entry_value(row[j],value) for j in columns] for row in M];n=len(rows)
 if not cubic:return [[u[0] for u in row] for row in rows]
 out=[[0]*(3*n) for _ in range(3*n)]
 for i,row in enumerate(rows):
  for j,(x,y,z) in enumerate(row):
   # Column k is multiplication of x+y*s+z*s^2 by s^k,
   # reduced by3s^3-2. Multiply the entire regular representation by3.
   block=[[3*x,2*z,2*y],[3*y,3*x,2*z],[3*z,3*y,3*x]]
   for u in range(3):
    for v in range(3):out[3*i+u][3*j+v]=block[u][v]
 return out

verified=[];stripped={};point_total=0
for rec in norms['charts']:
 name,order=rec['name'],rec['order'];yes(order in [0,1]);F=chart_models[name];M=rectangular(F,order);n=len(M);cubic=name=='s-cubic';size=n*(3 if cubic else 1);bound=size*8
 for row in rec['minors']:
  columns=list(range(n-1))+[row['last_column']];yes(len(set(columns))==n and row['last_column'] in range(n-1,len(M[0])))
  polynomial=fmpz_poly(list(map(int,row['coefficients'])));yes(polynomial.degree()<=bound)
  # Full degree_bound+1 INTEGER evaluations prove the polynomial identity.
  offset=bound//2
  for index,av in enumerate(range(-offset,bound-offset+1)):
   mat=matrix_at(M,columns,av,cubic);yes(len(mat)==size);value=int(polynomial(av));yes(int(fmpz_mat(mat).det())==value);yes(bareiss(mat)==value);point_total+=1
   if size>=45 and (index+1)%128==0:emit('norm_grid',chart=name,order=order,last_column=row['last_column'],points=index+1)
  H=fmpz_poly(list(map(int,row['stripped_coefficients'])));power=row['denominator_factor_power'];yes(polynomial==int(row['integer_content'])*fmpz_poly([1,2])**power*H);yes(H and sum(int(c)*(-1)**j*2**(int(H.degree())-j) for j,c in enumerate(H.coeffs()))!=0)
  stripped.setdefault((name,order),[]).append(H);verified.append({'chart':name,'order':order,'last_column':row['last_column'],'matrix_size':size,'degree_bound':bound,'degree':int(polynomial.degree()),'points':bound+1,'denominator_power':power,'stripped_degree':int(H.degree())})
 emit('norm_record_pass',chart=name,order=order)
# Exact nonzero integer resultants via full actual-degree Sylvester matrices.
prime=1000000007;yes(all(prime%d for d in range(2,math.isqrt(prime)+1)));det_records=[]
for name in chart_models:
 pair=stripped[(name,1)];yes(len(pair)==2);M=sylv(list(map(int,pair[0].coeffs())),list(map(int,pair[1].coeffs())));red=[[c%prime for c in row] for row in M];det=int(nmod_mat(red,prime).det());yes(det!=0 and det_plain_mod(red,prime)==det);det_records.append({'name':name,'test':'order1_coprimality','shape':len(M),'det_mod':det})
 H=stripped[(name,0)][0];yes(H.degree()==(72 if name=='s-zero' else 216));M=sylv(list(map(int,H.coeffs())),list(map(int,H.derivative().coeffs())));red=[[c%prime for c in row] for row in M];det=int(nmod_mat(red,prime).det());yes(det!=0 and det_plain_mod(red,prime)==det);det_records.append({'name':name,'test':'projection_squarefree','shape':len(M),'det_mod':det})
yes(stripped[('s-zero',0)][0][0]==0);yes(base['pole_chart_target_count']==24 and base['pole_chart_fibre_counts']=={'7':1,'8':23})
receipt={'result':'PASS','claim_id':'GAOF7TRP1','scope':'For normalized F7 over C with a marked triple-or-higher eliminant root s satisfying D0=s*(3s^3-2)=0, the complete locus is four rational curves parametrized by a!= -1/2. The residual Q9 and its derivative have gcd degree at most1 at EVERY point. Exactly288 normalized targets have seven-point fibres (72 at each carrier root), exactly23 have eight-point fibres, all others have nine-point fibres. The seven-point set contains the known fivefold origin and287 triple-plus-double targets. The eight-point set is the23 ordinary quadruple targets from GAOF7FOUR1. All311 nonordinary targets are explicitly projected by the degree72/216 discriminant polynomials and the old pole parameterization. Every point has a unique triple-or-higher root, so no target is counted twice. In particular no3+2+2+2 case lies on D0=0. The generic D0!=0 case and JC2 remain open.','checks':checks,'integer_norm_evaluation_points':point_total,'integer_determinant_engines':['FLINT','separate Python Bareiss'],'modular_determinant_engines':['FLINT','separate Python Gaussian'],'forbidden_coordinate_sylvester_determinant':integer_det,'norms':verified,'coprimality_and_squarefree':det_records,'normalized_counts':{'7':288,'8':23,'nonordinary_total':311,'9':'infinite ordinary open on each of four rational curves'},'remaining_three_point_gate':'3+2+2+2 with every marked triple root satisfying D0!=0','checker_sha256':sha(Path(__file__)),'hashes':{f:sha(HERE/f) for f in ['pole-source-candidate.json','lower-order-norms-candidate.json']},'source_sha256':sha(HERE.parent/'gao-f7-quadruple/quadruple-input.json'),'base_checker_sha256':sha(HERE.parent/'gao-f7-quadruple/verify.py'),'python':platform.python_version(),'sympy':S.__version__,'python_flint':flint.__version__,'seconds':time.monotonic()-START,'assurance':'Exact coefficient identities, complete integer norm interpolation, full actual-degree integer Sylvester nonvanishing and written source/rank/norm/counting arguments. No independent end-to-end, formal, external or novelty review.'}
(HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt),flush=True)
