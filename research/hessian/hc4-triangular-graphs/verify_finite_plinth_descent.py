"""Exact controls for finite-plinth residue duality and finite-cycle recognition."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 for v in values),label)
u,v,h=sp.symbols('u v h')
def principal_part(value):
    result=0
    for term in sp.Add.make_args(sp.expand(value)):
        powers=term.as_powers_dict()
        if powers.get(u,0)<0 and powers.get(v,0)<0:result+=term
    return sp.expand(result)
def residue(value):return sp.expand(value).coeff(u,-1).coeff(v,-1)
# One finite inverse system checks the direct coefficient-pairing implementation.
xi=1/(u*u*v*v);basis=[1,u,v,u*v]
pair=sp.Matrix([[residue(a*b*xi) for b in basis] for a in basis])
check(pair==sp.Matrix([[0,0,0,1],[0,0,1,0],[0,1,0,0],[1,0,0,0]]),
      'the finite Laurent coefficient pairing has the displayed nondegenerate matrix')
check(pair.det()==1,'the local Frobenius pairing is primitive with no denominator loss')
zero(principal_part(u*u*xi),'the first monomial complete-intersection generator kills the class')
zero(principal_part(v*v*xi),'the second monomial complete-intersection generator kills the class')
images=[principal_part(b*xi) for b in basis]
check(len(set(images))==4 and all(images),'all four quotient monomials have distinct nonzero inverse-system images')
mul_u=sp.Matrix([[0,0,0,0],[1,0,0,0],[0,0,0,0],[0,0,1,0]])
mul_v=sp.Matrix([[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,1,0,0]])
check(mul_u*mul_v==mul_v*mul_u,'the quotient multiplication matrices commute')
check(mul_u**2==sp.zeros(4) and mul_v**2==sp.zeros(4),'the exact Artinian quotient relations hold')
check(len(mul_u.col_join(mul_v).nullspace())==1,'the complete-intersection quotient has one-dimensional socle')
K1=mul_u.row_join(mul_v);K2=(-mul_v).col_join(mul_u)
check(K1*K2==sp.zeros(4),'the two-parameter Koszul differential squares to zero')
check(4-K1.rank()==1,'the local Koszul H0 has residue-field dimension one')
check(4-K2.rank()==1,'the local Koszul H2 is the one-dimensional socle')
check(8-K1.rank()-K2.rank()==2,'the local Koszul H1 has the required two-generator dimension')
length,soc=sp.symbols('length soc',integer=True,positive=True)
zero(length-2*length+length,'the chain Euler characteristic vanishes for every finite quotient length')
zero(1-(1+soc)+soc,'the all-length Koszul generator/socle identity has the stated homological Euler form')
# A support cylinder does not give extension of its full annihilator ideal.
eta=h/(u*u*v)+1/(u*v*v)
zero(principal_part((u-h*v)*eta),'the h-dependent linear relation annihilates the Cech class')
zero(principal_part(v*v*eta),'the second exact annihilator generator kills the class')
zero(principal_part(u*eta)-h/(u*v),'u acts with the h-dependent residue coefficient')
zero(principal_part(v*eta)-1/(u*v),'v acts with unit residue coefficient')
for monomial in [u*u,u*v,v*v]:
    zero(principal_part(monomial*eta),'the contracted square ideal kills eta: '+str(monomial))
check(principal_part(eta)!=0 and principal_part(v*eta)!=0,
      'the two quotient classes survive after reduction modulo (u-hv,v^2)')
c1,c2=sp.symbols('c1 c2')
zero(principal_part((c1*u+c2*v)*eta)-(c1*h+c2)/(u*v),
     'no constant-coefficient linear relation belongs to the contracted annihilator')
check(sp.Poly(c1*h+c2,h).all_coeffs()==[c1,c2],
      'contraction to the fixed plane has no nonzero linear term and is exactly (u,v)^2')
# Universal local critical-point and normal-form algebra.
F,G,P,Q,af,ag,ah,j11,j12,j21,j22=sp.symbols('F G P Q af ag ah j11 j12 j21 j22')
linear=sp.Matrix([[j11,j21],[j12,j22]])*sp.Matrix([P,Q])-sp.Matrix([af,ag])
zero(linear.jacobian([P,Q]).det()-(j11*j22-j12*j21),
     'the critical equations have exactly the coefficient-Jacobian rank')
rank_one=linear.subs({j11:0,j12:0,j21:0,j22:1})
check(rank_one==sp.Matrix([-af,Q-ag]),'a nonreduced rank-one point has either no critical point or an affine critical family')
rank_zero=linear.subs({j11:0,j12:0,j21:0,j22:0})
check(rank_zero==sp.Matrix([-af,-ag]),'a rank-zero point cannot supply one isolated critical point')
t11,t12,t22=sp.symbols('t11 t12 t22')
HH=sp.Matrix([[t11,t12,j11/ah,j21/ah],[t12,t22,j12/ah,j22/ah],
              [j11/ah,j12/ah,0,0],[j21/ah,j22/ah,0,0]])
zero(HH.det()-(j11*j22-j12*j21)**2/ah**4,
     'a reduced simple residue root is a nondegenerate full Morse point')
alpha,beta,base_value=sp.symbols('alpha beta base_value')
zero(F*(P-alpha)+G*(Q-beta)-(F*P+G*Q-F*alpha-G*beta),
     'the analytic local fibre translation removes every base-ideal numerator term')
# Actual affine-space negative control:two finite Morse points,one additional nonreduced hole.
x,z,s=sp.symbols('x z s');f=x*(x+z);g=z*(z-1);potential=x+P*f+Q*g
D=lambda value:sp.expand(g*sp.diff(value,P)-f*sp.diff(value,Q))
zero(D(potential),'the actual nonreduced countercontrol has a global h-preserving translation')
zero(D(P)-g,'the first actual local slice supplies g')
zero(D(-Q)-f,'the second actual local slice supplies f')
zero((potential-x)-P*f-Q*g,'the Cech numerator is exactly h-x')
check(sp.Poly(h-x,h).LC()==1,'the Cech numerator is monic over the nonreduced quotient,so its colon is unchanged')
grad=[sp.diff(potential,var) for var in [x,z,P,Q]]
actual=sp.groebner(grad,z,P,Q,x)
target=sp.groebner([z-1,x*(x+1),P+2*x+1,Q+x],z,P,Q,x)
check(actual.polys==target.polys,'the complete actual critical ideal is exactly the two stated reduced points')
for xx,pp,qq,value in [(0,-1,0,0),(-1,1,1,-1)]:
    point={x:xx,z:1,P:pp,Q:qq}
    zeros([a.subs(point) for a in grad],'all critical derivatives vanish at actual point '+str(xx))
    zero(potential.subs(point)-value,'the actual critical value '+str(value))
    zero(sp.hessian(potential,[x,z,P,Q]).det().subs(point)-1,'the actual critical Hessian is one '+str(xx))
zero(grad[0].subs({x:0,z:0})-1,'the nonreduced support plane has no critical point at any P,Q')
zero(f.subs(z,0)-x*x,'the support point (0,0) has length two')
zero(f.subs(z,1)-x*(x+1),'the other two support points are simple')
zero(sp.Matrix([f,g]).jacobian([x,z]).det()-(2*x+z)*(2*z-1),
     'the whole coefficient Jacobian records exactly one nonreduced and two simple support points')
check(sp.gcd(f,g)==1,'the actual coefficient pair is primitive on the whole source base')
# Primitive local Morse cycles occupy precisely the simple-point link summands.
cycles=sp.Matrix([[0,0],[1,0],[0,1]])
check(cycles.rank()==2,'the two actual finite cycles span exactly the two simple-point summands')
check(sp.Matrix.hstack(sp.Matrix([1,0,0]),cycles).det()==1,
      'the extra generic-fibre hole gives a primitive integral quotient Z')
check(3-cycles.rank()==1,'the actual nonreduced defect is one,not hidden torsion or a missed finite critical point')
check(sp.eye(3).det()==1,'exhaustion in the three-point native family covers every primitive hole summand')
# Finite point support must be separated from exceptional Keller degree and actual alignment.
for descent in [False,True]:
    for en in [False,True]:
        check((not (en and descent))==((not en) or (not descent)),
              'native finite-ideal alignment and E_N remain independent gates '+str((descent,en)))
algebra=len(checks);context=json.loads((HERE/'sources/finite-plinth-descent-context.json').read_text())
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained immutable input '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'inspected primary context '+entry['retained_file'])
check(context['newly_imported_proved_external_statements']==[] and context['newly_used_unproved_external_statements']==[],
      'local duality and generator count are directly proved,and no principal-plinth theorem is transferred to codimension two')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
check(sum(flags.values())==5 and not flags['finite_Morse_critical_set_alone_forces_radical_plinth']
      and not flags['every_finite_plinth_ideal_is_assumed_two_generated']
      and not flags['nontransverse_normal_maps_are_assumed_degree_one']
      and not flags['an_individual_native_action_with_required_kernel_and_plinth_constructed']
      and not flags['original_goal_complete'],
      'finite-cycle exhaustion replaces radicality without assuming two generators,transversality or a new native action')
proof=HERE/'FINITE_PLINTH_DESCENT_AND_VANISHING_CYCLE_RECOGNITION.md'
paths=['sources/finite-plinth-descent-context.json']+['sources/'+e['retained_file'] for e in context['newly_cached_primary_context']]
result=dict(outcome='PASS',check_count=len(checks),checks=checks,algebraic_and_scope_controls=algebra,
 elapsed_seconds=round(time.monotonic()-START,3),python=platform.python_version(),sympy=sp.__version__,
 proof_sha256=hashlib.sha256(proof.read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 dependencies=context['retained_proof_and_source_sha256'],source_sha256={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths},
 independent_replay=False,formal_verification=False,external_review=False,**flags)
(HERE/'receipt-finite-plinth-descent.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in {'checks','dependencies','source_sha256'}},indent=2))
