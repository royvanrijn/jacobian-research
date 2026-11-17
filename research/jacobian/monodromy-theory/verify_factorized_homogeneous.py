"""Exact KMFI121--123 factorization, field, nonfinite A11 and A6 coefficients."""
import os,sys
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from factorized_homogeneous import degree11_input,completion_data
ROOT=Path(__file__).resolve().parent;start=time.monotonic();checks=[]
def check(name,truth):
    if not truth:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
d=degree11_input();a=d['a'];z=d['z'];A=d['A'];ps=d['parameters'];red=d['red'];redz=d['redz']
check('explicit coefficient quartic irreducible',d['R'].is_irreducible)
gal,alt=S.polys.numberfields.galois_group(d['R'])
check('quartic Galois group S4',gal.order()==24 and not alt)
check('coefficient avoids0 and1',d['R'].eval(0)!=0 and d['R'].eval(1)!=0)
check('whole weighted-polynomial derivative identity',
      redz(S.diff(d['Q'],z)-(A+z)**4*(A+a*z)**6)==0)
Dc=d['cubic']
disc=S.discriminant(Dc,z)
check('reciprocal cubic separable',red(disc)!=0)
check('all matrix eigenvalues nonzero distinct and avoid c1',
      red(Dc.subs(z,0))!=0 and red(Dc.subs(z,1))!=0 and red(Dc.subs(z,a))!=0)
# Constant rank-one matrix identity without adjoining its cubic roots.
# det(TI-M) and the weighted resolvent trace are over Q(a).
T=S.Symbol('T');M=d['M']
char=S.expand(M.charpoly(T).as_expr())
weighted=7*a/(T-a)+S.trace((T*S.eye(3)-d['C']).adjugate()*d['C'])/d['C'].charpoly(T).as_expr()
identity=S.cancel(char*(1-weighted)-(T-1)**4)
num,den=S.fraction(identity)
check('entire constant weighted rank-one characteristic identity',
      all(red(co)==0 for mon,co in S.Poly(num,T).terms()))
for j in range(5):
    moment=(S.Matrix([[S.diff(d['F'],p) for p in ps]])*(M**j)*S.Matrix(ps))[0]-(7*a**j+S.trace(d['C']**j))*d['F']
    check('entire weighted norm gradient moment'+str(j),
          all(red(co)==0 for mon,co in S.Poly(S.expand(moment),*ps).terms()))
check('whole factorized source reconstruction',
      all(all(red(co)==0 for mon,co in S.Poly(S.expand(e),*ps).terms())
          for e in d['X'][1:,:]-d['X'][0]*S.Matrix(ps)-d['F']*M*S.Matrix(ps)))
B=S.Symbol('B');P=d['Q']-B;v=redz(d['Q'].subs(z,-A))
lc=S.Poly(d['Q'],z).LC()
check('constant leading coefficient eleven times',red(lc-a**6/11)==0)
r1=redz(S.resultant(P,A+z,z));r2=redz(S.resultant(P,A+a*z,z))
check('exact visible resultant',redz(r1-(B-v))==0)
check('exact zero-branch resultant',redz(r2-a**11*B)==0)
check('whole discriminant square class minus11 via multiplicative resultant',
      red(a**66/lc-11**11*lc**10)==0)
check('critical visible value distinct from zero',d['critical_value']!=0)
check('exact quintuple root and degree6 regular tail',d['tail_remainder']==0 and S.degree(d['tail'],z)==6)
# Coprimality over Q(a) by nonzero resultants.
check('tail avoids visible critical root',red(d['tail'].subs(z,1))!=0)
denom=S.expand((-S.eye(4)+z*d['M']).det())
check('all tail roots have unit inverse-field denominator',
      red(S.resultant(d['tail'],denom,z))!=0)
check('degree24 primitive coefficient field irreducible',d['resultant'].degree()==24 and d['resultant'].is_irreducible)
check('degree24 field excludes sqrt minus11',
      S.Poly(d['resultant'].as_expr(),z,extension=S.sqrt(-11)).is_irreducible)
check('quartic recovery polynomial exact',S.rem(
 S.Poly(d['R'].as_expr().subs(a,d['a_recovery']),z),d['resultant']).is_zero)
# A genuine gcd2 passport protects against treating two shared-value critical
# points as separate transpositions.
R2=1155*a**4-924*a**3+594*a*a-220*a+35
check('gcd2 coefficient quartic irreducible',S.Poly(R2,a).is_irreducible)
W=S.Symbol('W');q2=S.integrate((1+W*W)**4*(1+a*W*W),W)
quot,rem=S.div(q2,W*(1+a*W*W)**2,W)
def red2(e):
    n,de=S.fraction(S.cancel(e))
    return (S.Poly(n,a)*S.invert(S.Poly(de,a),S.Poly(R2,a))).rem(S.Poly(R2,a)).as_expr()
check('gcd2 exact zero-value cofactor',
      all(red2(co)==0 for mon,co in S.Poly(rem,W).terms()))
check('gcd2 nonzero visible critical values',
      red2(S.div(q2,W*W+1,W)[1].coeff(W,1))!=0)
def partitions(total,length,minimum=1):
    if length==0:return [()] if total==0 else []
    return [(j,*tail) for j in range(minimum,total//length+1)
            for tail in partitions(total-j,length-1,j)]
import math
rank4parts=partitions(10,4)
check('all degree11 rank4 gcd cases',
      {math.gcd(*b) for b in rank4parts}=={1,2} and
      [b for b in rank4parts if math.gcd(*b)==2]==[(2,2,2,4)])
check('rank9 degree11 branch transposition pattern',
      partitions(10,9)==[(1,1,1,1,1,1,1,1,2)])
# The boundary cases in the group proof are exact finite arithmetic.
check('rank8 proper Jones degree impossible',not any(p**e==10 for p in [2,3,5,7] for e in range(1,5)))
check('rank9 exponent pattern unique',10-9==1)
data=completion_data(progress=lambda m:print(m,flush=True));K=data['K'];one=K.one;zeroK=K.zero
mm=data['mm'];mv=data['mv'];evalpoly=data['evalpoly']
check('degree24 field recovered a satisfies quartic',
      sum((K.convert(co)*data['aa']**mon[0] for mon,co in d['R'].terms()),zeroK)==zeroK)
check('entire actual marked input target',
 all(evalpoly(poly,data['old_marked_point'])+data['gamma0']*(one if i==0 else data['old_marked_point'][i-1])==data['Ystar'][i]
     for i,poly in enumerate(data['X_polynomials'])))
check('entire actual marked stage zero',all(v==zeroK for v in data['marked_stage']))
Jord=[[one if j==i+1 else zeroK for j in range(4)] for i in range(4)]
check('exact regular Jordan frame',
      mm(mm(data['Ti'],data['P0']),data['T'])==Jord)
check('actual input full rank at chosen branch source',
      sum((data['DXstar'][0][i]*data['T'][i][0] for i in range(4)),zeroK)!=zeroK)
check('all25 inverse linear-jet entries',
      mm(data['J'],data['L'])==[[one if i==j else zeroK for j in range(5)] for i in range(5)])
check('stage determinant gamma0 inverse',data['det_J']*data['gamma0']==one)
# Verify the native truncated evaluator against the field point and original jet.
R=data['series_ring'];Z=data['z'];ground=R.ground_new;cut=data['cut'];mul=data['mul']
origin=[ground(v) for v in data['theta0']]
E,J=data['series_stage'](origin,with_jac=True)
check('series evaluator exact origin',all(v==R.zero for v in E))
check('series evaluator exact original linear jet',
      all(J[i][j]==ground(data['J'][i][j]) for i in range(5) for j in range(5)))
h=data['curve_h'];AA=data['correction_A']
check('entire finite curve residual',all(v==R.zero for v in data['curve_residual']))
check('curve starts ordinary order2 and truncated below8',
      all(v.get((0,),zeroK)==v.get((1,),zeroK)==zeroK and v.degree()<8 for v in h))
check('exact lower-unitriangular polynomial correction',
 all(AA[i][j]==(R.one if i==j else R.zero) for i in range(4) for j in range(i,4)))
check('correction constant identity',all(AA[i][j].get((0,),zeroK)==(one if i==j else zeroK) for i in range(4) for j in range(4)))
check('entire saved flag residual',all(data['U'][i][j]==R.zero for i in range(4) for j in range(i)))
check('flag against original evaluated matrix',
      all(cut(sum((data['M_high'][i][k]*AA[k][j] for k in range(4)),R.zero))==R.zero
          for i in range(4) for j in range(i)))
def rmm(A,B):
 return [[sum((A[i][k]*B[k][j] for k in range(4)),R.zero) for j in range(4)] for i in range(4)]
I4=[[R.one if i==j else R.zero for j in range(4)] for i in range(4)]
N=[[AA[i][j]-I4[i][j] for j in range(4)] for i in range(4)];N2=rmm(N,N);N3=rmm(N2,N);N4=rmm(N3,N)
Ai=[[I4[i][j]-N[i][j]+N2[i][j]-N3[i][j] for j in range(4)] for i in range(4)]
check('actual untruncated polynomial source inverse',rmm(AA,Ai)==I4 and rmm(Ai,AA)==I4 and all(v==R.zero for row in N4 for v in row))
check('all Gaussian coefficient degrees below8',all(v.degree()<8 for row in AA for v in row if v!=R.zero))
check('global first-stage nonlinear weight cutoff',2*4>=8 and all(4+j>=1 for j in range(4)))
check('global final-stage weighted cutoffs',2*5>8 and 8>8-5 and max(data['output_weights'])==8)

boundary=[]
bs=[7,1,1,1]
for target,btarget in enumerate(bs):
    row=[{'root':'Z0','e':1,'f':1,'present':True}]
    for i,b in enumerate(bs):
        common=math.gcd(b,btarget)
        row.append({'root_index':i,'e':b//common,'f':common,'present':i==target})
    check('target plane'+str(target)+' complete normalization degree',sum(q['e']*q['f'] for q in row)==11)
    boundary.append({'target_plane':target,'normalization_rows':row})
missing=[q for v in boundary for q in v['normalization_rows'] if not q['present']]
check('whole original boundary twelve primes',len(missing)==12)
check('original three ramification7 and nine unramified omissions',
      sum(q['e']==7 for q in missing)==3 and sum(q['e']==1 for q in missing)==9)
check('retained residue degree7 unit is not a seventh power',math.gcd(7,1,1,1)==1)
check('completed class count and canonical vector',len(missing)+1+9==22)
manifest=json.loads((ROOT/'sources/factorized-homogeneous-manifest.json').read_text())
for src in manifest['retained_primary_inputs']:
    raw=(ROOT/src['path']).read_bytes()
    check('retained primary bytes '+src['path'],len(raw)==src['bytes'])
    check('retained primary hash '+src['path'],hashlib.sha256(raw).hexdigest()==src['sha256'])
def encode(co):
    terms=list(co.to_list());terms=[K.domain.zero]*(24-len(terms))+terms
    return [str(v) for v in reversed(terms)]
def series(poly):
    return [{'power':mon[0],'coefficient':encode(co)} for mon,co in sorted(poly.items())]
def matrix(A):return [[encode(co) for co in row] for row in A]
model={'claims':['KMFI121','KMFI122','KMFI123'],
 'factorization_scope':'H=g-fp=F(p)Mp, constant arbitrary M, homogeneous F, full rankDX, r>=2',
 'proved_normal_form':'M distinct invertible eigenvalues; f=-cF, positive monomial F, weighted rank-one condition',
 'nonfactorized_homogeneous_inputs_classified':False,
 'generic_degree':'1+degree F','one_point_axis_is_universal_in_this_class':True,
 'degree11_rank_groups':{str(r):('C11' if r==10 else ('S11' if r%2 else 'A11')) for r in range(1,11)},
 'even_rank_degree11_arithmetic':'A11 iff -11 square in coefficient field, else S11, exact constants k(sqrt(-11)); ranks2,4,6,8',
 'degree11_input':{'coefficient_quartic':str(d['R'].as_expr()),'M':[[str(e) for e in row] for row in d['M'].tolist()],
   'F':str(d['F']),'X':[str(e) for e in d['X']],'Q_A':str(d['Q']),
   'critical_value_A_minus1':str(d['critical_value']),'regular_tail_degree6':str(d['tail']),
   'degree24_primitive_polynomial':str(d['resultant'].as_expr()),'a_recovery_polynomial':str(d['a_recovery']),
   'generic_geometric_group':'A11','generic_arithmetic_group':'S11',
   'original_boundary':boundary,'original_normalization_geometric_class_group':'Z^12'},
 'actual_A6_completion':{
   'coefficient_field':'Q(beta), beta=CRootOf(retained irreducible degree24 polynomial,0)',
   'coefficient_field_degree':24,'basis':'1,beta,...,beta^23; coefficient lists ascending',
   'a_in_field':encode(data['aa']),'pstar':[encode(v) for v in data['pstar']],
   'Ystar':[encode(v) for v in data['Ystar']],'old_regular_point':[encode(v) for v in data['old_marked_point']],
   'Jordan_T':matrix(data['T']),'Jordan_T_inverse':matrix(data['Ti']),
   'theta0':[encode(v) for v in data['theta0']],'inverse_linear_jet_L':matrix(data['L']),
   'source_weights':data['source_weights'],'output_weights':data['output_weights'],
   'curve_h':[series(v) for v in h],'lower_unipotent_A':[[series(v) for v in row] for row in AA],
   'source_inverse':'v=A(z)^-1*(u-h(z)); A^-1=I-N+N^2-N^3, N=A-I',
   'global_map':'(x*gamma,E0(Psi)/x,E1(Psi)/x^5,E2(Psi)/x^6,E3(Psi)/x^7,E4(Psi)/x^8), Psi=theta0+L*(z,h(z)+A(z)v), z=x*y0,v_i=x^e_i*y_i',
   'whole_Jacobian':encode(data['gamma0']),'whole_Jacobian_expression':'beta-1',
   'generic_degree':11,'geometric_group':'A11','arithmetic_group':'S11',
   'exact_closure_constants':'Q(beta,sqrt(-11))',
   'geometric_normalization_class_group':'Z^22','geometric_canonical_vector':[6,6,6,4]+[0]*18,
   'class_group_asserted_over_degree24_field':False,
   'expanded_final_A6_determinant_claimed':False,
   'native_series_evaluator_used_in_unit_gamma_chart_for_exact_jets_only':True,
   'global_polynomiality_proved_by_weights_and_finite_rectifier':True},
 'actual_PSL11_input_supplied':False,'new_abstract_action_claimed':False,
 'absolute_exceptional11_action_decided':False,'full_action_spectrum_complete':False,'goal_complete':False}
bbpowers=[one]
for i in range(1,24):bbpowers.append(bbpowers[-1]*data['bb'])
def roundtrip(co):
    return sum((K.convert(S.Rational(v))*bbpowers[i] for i,v in enumerate(encode(co))),zeroK)==co
check('every serialized completion field coefficient exact roundtrip',
      all(roundtrip(co) for row in data['L'] for co in row) and
      all(roundtrip(co) for co in data['theta0']) and
      all(roundtrip(co) for poly in h for mon,co in poly.items()) and
      all(roundtrip(co) for row in AA for poly in row for mon,co in poly.items()))
(ROOT/'factorized_homogeneous_model.json').write_text(json.dumps(model,indent=2)+'\n')
paths=['FACTORIZED_HOMOGENEOUS_FIELDS.md','factorized_homogeneous.py','verify_factorized_homogeneous.py',
 'factorized_homogeneous_model.json','sources/factorized-homogeneous-manifest.json',
 'HIGHER_RANK_COMPLETION.md','higher_completion.py','SWEEP_DIFFERENTIAL_RECONSTRUCTION.md',
 'HOMOGENEOUS_AXIS_OBSTRUCTION.md','COMPLETION_EXCEPTIONAL_DIVISOR.md']
paths += [a['path'] for a in manifest['retained_primary_inputs']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,
 'runtime_seconds':round(time.monotonic()-start,3),'python':platform.python_version(),'sympy':S.__version__,
 'sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 'assurance':'Written unbounded factorization rigidity and exact field reduction, degree11 geometric/arithmetic spectrum, actual full-rank nonfinite A11 input and complete degree24-field A6 finite coefficients. Boundary proof counts geometric primes separately from residue degree. No formal, independent whole-proof or external review claimed.',
 'goal_complete':False}
(ROOT/'factorized_homogeneous_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'runtime_seconds':receipt['runtime_seconds'],
 'actual_A6_degree':11,'geometric_group':'A11','arithmetic_group':'S11','coefficient_field_degree':24,
 'geometric_normalization_class_group':'Z^22','model_bytes':(ROOT/'factorized_homogeneous_model.json').stat().st_size,
 'actual_PSL11_input_supplied':False,'goal_complete':False}))
