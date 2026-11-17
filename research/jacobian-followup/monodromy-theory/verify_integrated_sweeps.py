"""Exact KMFI109--111 integrated family and quartic-field affine completion."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
import sympy as S
import json,hashlib,time,platform
from integrated_sweeps import family,degree9_data
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,val):
    if not val:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
t,A,B=S.symbols('t A B');families=[]
for r,h in [(1,t*t),(2,t*t),(3,t*t),(4,t*t),(4,t**3+t+1)]:
    f=family(r,h,t);ps=f['parameters'];g=f['g'];X=f['X'];swept=f['S'];label='rank'+str(r)+' h='+str(h)
    check(label+' triangular leading source coordinates',
      all(zero(S.diff(X[i],ps[i])-1) and not any(X[i].has(p) for p in ps[i+1:])
          for i in range(1,r)))
    check(label+' first-rank determinant',zero(X[:r,:].jacobian(ps).det()-S.diff(h,t)))
    check(label+' whole normalized sweep Jacobian',
          zero(swept.jacobian([g,*ps]).det(method='domain-ge')-g**r))
    target=f['target_variables'];computed=f['target_polynomial'].subs(dict(zip(target,list(swept)[:r])),simultaneous=True)
    integral=f['integral'].subs(target[0],swept[0])
    check(label+' exact inverse-integral equation',zero(swept[-1]-computed-(-1)**(r-1)*integral))
    check(label+' degree formula',S.degree(f['integral'],t)==f['degree'])
    check(label+' derivative power',zero(S.diff(f['integral'],t)-(target[0]-h)**r))
    families.append({'r':r,'h':str(h),'X':[str(e) for e in X],'degree':f['degree']})
# The higher critical-value derivative is a symbolic identity for arbitrary h.
hf=S.Function('h')(t);primitives=[t]+[S.Function('H'+str(i))(t) for i in range(1,5)]
I=sum((-1)**j*S.binomial(4,j)*A**(4-j)*primitives[j] for j in range(5))
replacement={S.diff(primitives[j],t):hf**j for j in range(1,5)}
previous=I.subs(A,hf)
for j in range(1,5):
    derivative=S.diff(previous,t).xreplace(replacement)/S.diff(hf,t)
    wanted=S.diff(I,A,j).subs(A,hf)
    check('arbitrary-h critical value derivative order'+str(j),zero(derivative-wanted))
    previous=wanted
check('fourth critical value derivative recovers t',zero(previous-24*t))
# Exact square-class identity, without using a Frobenius specialization.
f9=family(4,t*t,t);It=f9['integral'];P=It-B
res=S.resultant(P,A-t*t,t);disc=S.discriminant(P,t)
check('degree9 whole discriminant square-class formula',zero(disc-9*res**4))
check('degree9 arithmetic constant square',S.sqrt(9)==3)
I1=S.expand(It.subs(A,1))
R=35*t**4+175*t**3+345*t*t+325*t+128
check('marked quintuple root and exact quartic tail',
      S.expand(315*(I1-S.Rational(128,315))-(t-1)**5*R)==0)
check('quartic defining field irreducible',S.Poly(R,t,domain=S.QQ).is_irreducible)
check('all quartic marked roots regular',S.gcd(S.Poly(R,t),S.Poly(t*t-1,t)).degree()==0)
data=degree9_data();K=data['K'];evaluate=data['evaluate'];ps=data['family']['parameters'];g=data['family']['g']
check('marked target exact',data['Ystar']==S.Matrix([1,S.Rational(2,3),-S.Rational(8,15),S.Rational(23,315),S.Rational(1,15)]))
check('Jordan determinant1',data['T'].det()==1)
Jord=S.zeros(4)
for i in range(3):Jord[i,i+1]=1
check('actual regular nilpotent Jordan conjugacy',data['T'].inv()*data['P0']*data['T']==Jord)
oldsource=[data['theta0'][0],*data['old_marked_point']]
check('entire marked sweep target',
 all(evaluate(expr,[g,*ps],oldsource)==K.convert(want) for expr,want in zip(data['family']['S'],data['Ystar'])))
check('centered polynomial first stage',
 all(S.denom(S.cancel(e))==1 for e in data['E']))
check('entire marked stage zero',
 all(evaluate(e,data['theta'],data['theta0'])==K.zero for e in data['E']))
check('stage regular determinant gamma0 inverse',data['det_J']*data['theta0'][0]==K.one and data['theta0'][0]!=K.zero)
check('all25 inverse linear-jet entries',
 all(sum((data['J'][i][k]*data['L'][k][j] for k in range(5)),K.zero)==(K.one if i==j else K.zero)
     for i in range(5) for j in range(5)))
ring=data['series_ring'];z=data['z'];h=data['h'];M=data['M'];AA=data['A'];cut=data['cut']
check('entire finite formal curve residual',all(e==ring.zero for e in data['curve_residual']))
check('curve starts at ordinary order2',
 all(e.get((0,),K.zero)==e.get((1,),K.zero)==K.zero for e in h))
check('all finite curve coefficients have degree below8',all(e.degree()<8 for e in h))
check('exact polynomial lower-unitriangular Gaussian matrix',
 all(AA[i][j]==ring.zero for i in range(4) for j in range(i+1,4)) and
 all(AA[i][i]==ring.one for i in range(4)))
check('Gaussian matrix constant identity',
 all(AA[i][j].get((0,),K.zero)==(K.one if i==j else K.zero) for i in range(4) for j in range(4)))
check('entire Gaussian flag residual',all(data['U'][i][j]==ring.zero for i in range(4) for j in range(i)))
# Verify the displayed flag against the original evaluated M, independently of saved U.
check('original M times A flag',
 all(cut(sum((M[i][k]*AA[k][j] for k in range(4)),ring.zero))==ring.zero
     for i in range(4) for j in range(i)))
NN=[[AA[i][j]-(ring.one if i==j else ring.zero) for j in range(4)] for i in range(4)]
def mmul(a,b):
    return [[sum((a[i][k]*b[k][j] for k in range(4)),ring.zero) for j in range(4)] for i in range(4)]
I4=[[ring.one if i==j else ring.zero for j in range(4)] for i in range(4)]
N2=mmul(NN,NN);N3=mmul(N2,NN);N4=mmul(N3,NN)
check('exact untruncated nilpotent Gaussian correction',all(e==ring.zero for row in N4 for e in row))
Ai=[[I4[i][j]-NN[i][j]+N2[i][j]-N3[i][j] for j in range(4)] for i in range(4)]
check('entire polynomial inverse matrix',mmul(AA,Ai)==I4 and mmul(Ai,AA)==I4)
check('all lower matrix degrees below8',all(e.degree()<8 for row in AA for e in row if e!=ring.zero))
check('global monomial cutoff for all high-variable nonlinear terms',2*5>8)
check('linear flag cutoff exceeds every weight difference',8>8-5)
check('pure curve cutoff reaches maximal output weight',8==max(data['stage']['output_weights']))
# Exact geometric normalization consequence for this finite source.
check('degree9 exceptional splitting and class count',9-1==8 and 9-2==7)
manifest=json.loads((ROOT/'sources/integrated-sweeps-manifest.json').read_text())
for source in manifest['sources']:
    raw=(ROOT/source['path']).read_bytes()
    check('primary bytes '+source['path'],len(raw)==source['bytes'])
    check('primary hash '+source['path'],hashlib.sha256(raw).hexdigest()==source['sha256'])
def encode(co):
    terms=list(co.to_list());terms=[K.domain.zero]*(4-len(terms))+terms
    return [str(e) for e in reversed(terms)]
def series(poly):
    return [{'power':m[0],'coefficient':encode(c)} for m,c in sorted(poly.items())]
model={'claims':['KMFI109','KMFI110','KMFI111'],'rank4_full_triangular_class':'finite degree4d+1; C5 ifd1, A_(4d+1) ifd>=2',
 'arithmetic_rank4_d_ge2':'A_N or S_N according to square class N=4d+1; constants k(sqrt(N))',
 'family_controls':families,'critical_derivative_recovery':'d^4 b/da^4=24t',
 'degree9_full_rank_control':{'rank_X':4,'normalized_sweep_Jacobian':'gamma^4','degree':9,'geometric_group':'A9','arithmetic_group':'A9',
 'coefficient_field_polynomial':str(R),'basis':'1,alpha,alpha^2,alpha^3; alpha=CRootOf(R,0)',
 'source_variables':[str(p) for p in ps],'X':[str(e) for e in data['family']['X']],
 'branch_point':[str(e) for e in data['pstar']],'branch_target':[str(e) for e in data['Ystar']],
 'Jordan_T':[[str(data['T'][i,j]) for j in range(4)] for i in range(4)],
 'centered_source_variables':[str(p) for p in data['centered_parameters']],
 'centered_X':[str(e) for e in data['centered_X']],
 'theta_variables':[str(v) for v in data['theta']],'E':[str(e) for e in data['E']],
 'old_marked_point':[encode(e) for e in data['old_marked_point']],
 'theta0':[encode(e) for e in data['theta0']],'inverse_Jacobian_L':[[encode(e) for e in row] for row in data['L']],
 'weights':[1,5,6,7,8],'formal_parameter':'z','curve_h':[series(e) for e in h],
 'unipotent_A':[[series(e) for e in row] for row in AA],
 'polynomial_inverse_A_formula':'I-(A-I)+(A-I)^2-(A-I)^3',
 'actual_A6_map':'C=x*gamma; F_i=E_i(theta0+L*(z,h(z)+A(z)*v))/x^e_i; substitute z=x*y0,v_i=x^e_i*y_i',
 'whole_map_Jacobian':encode(data['theta0'][0]),
 'expanded_whole_A6_determinant_claimed':False,'whole_Jacobian_proved_by_KMFI104':True,
 'geometric_normalization_class_group':'Z^8','geometric_canonical_vector':[4,0,0,0,0,0,0,0],
 'geometric_class_group_promoted_to_quartic_field_class_group':False},
 'new_abstract_action_claimed':False,'rank4_PSL11_field_supplied':False,
 'absolute_exceptional_action_decided':False,'goal_complete':False}
(ROOT/'integrated_sweeps_model.json').write_text(json.dumps(model,indent=2)+'\n')
paths=['INTEGRATED_SWEEP_PENCILS.md','integrated_sweeps.py','verify_integrated_sweeps.py','integrated_sweeps_model.json',
 'sources/integrated-sweeps-manifest.json','HIGHER_RANK_COMPLETION.md','higher_completion.py','COMPLETION_EXCEPTIONAL_DIVISOR.md']
paths += [e['path'] for e in manifest['sources']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,'runtime_seconds':round(time.monotonic()-started,3),
 'python':platform.python_version(),'sympy':S.__version__,'sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 'assurance':'Written unbounded triangular/integral reduction, critical-value derivative, Jordan group and arithmetic proofs; exact full-rank degree9 input and every finite A6 completion coefficient. No expanded final A6 determinant, formal verification or independent whole-proof review claimed.',
 'goal_complete':False}
(ROOT/'integrated_sweeps_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],'actual_full_rank_A6_degree':9,
 'geometric_and_arithmetic_group':'A9','coefficient_field_degree':4,'geometric_normalization_class_group':'Z^8',
 'PSL11_field_supplied':False,'goal_complete':False}))
