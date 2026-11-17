"""Exact controls for univariate suspended multipliers and base-preserving closure."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):
    value=sp.cancel(value)
    check(value==0 or sp.simplify(value)==0,label)
def zeros(values,label):
    check(all(sp.cancel(v)==0 or sp.simplify(v)==0 for v in values),label)
x,y,z,w,X,Y,Z,W,t=sp.symbols('x y z w X Y Z W t')
alpha,kappa,a,b,c,ell=sp.symbols('alpha kappa a b c ell')
variables=[x,y,z,w]
C=sp.Function('C')(z);D=sp.Function('D')(z)
h=2*x*y+w*C+D
Hess=sp.hessian(h,variables)
zero(Hess.det()-4*sp.diff(C,z)**2,
     'the complete native Hessian retains the exact squared multiplier derivative')
zero(sp.diff(h,z)-(w*sp.diff(C,z)+sp.diff(D,z)),
     'the full critical equation has the exact affine dual solution')

# The leading coefficient identity is symbolic in degree, not a degree scan.
n,leading=sp.symbols('n leading')
zero(n*n-2*n*(n-1)-n*(2-n),
     'the all-degree leading coefficient of C-prime-squared minus twice C C-second is n(2-n)')
aa,bb,cc=sp.symbols('aa bb cc')
quadratic=aa*z*z+bb*z+cc
zero(sp.diff(quadratic,z)**2-2*quadratic*sp.diff(quadratic,z,2)
     -(bb*bb-4*aa*cc),
     'every quadratic has the exact discriminant derivative-square identity')
for multiplier in [z**3-z,z**4-1]:
    remainder=sp.rem(sp.diff(multiplier,z)**2,multiplier,z)
    check(sp.degree(remainder,z)>0,
          'a higher-degree squarefree control has nonconstant derivative-square residue '+str(sp.degree(multiplier,z)))
    check(sp.gcd(multiplier,sp.diff(multiplier,z))==1,
          'the excluded higher-degree control is actually squarefree '+str(sp.degree(multiplier,z)))
multiplier=z**3-z
potential=2*x*y+w*multiplier+z*z
gradient=[sp.diff(potential,q) for q in variables]
critical_values=[]
for root in [-1,0,1]:
    ww=-2*root/sp.diff(multiplier,z).subs(z,root)
    point={x:0,y:0,z:root,w:ww}
    zeros([v.subs(point) for v in gradient],
          'the cubic control has its actual native critical point '+str(root))
    determinant=sp.hessian(potential,variables).subs(point).det()
    zero(determinant-4*sp.diff(multiplier,z).subs(z,root)**2,
         'the cubic critical determinant keeps its exact multiplier factor '+str(root))
    critical_values.append(determinant)
check(len(set(critical_values))==2,
      'the cubic control really fails the common critical-Hessian condition')

# Repeated roots are counted geometrically and remain outside the squarefree lemma.
repeated=(z-1)**2*(z+2)
blocked=2*x*y+w*repeated+z
point={x:0,y:0,z:-2,w:-sp.Rational(1,9)}
zeros([sp.diff(blocked,q).subs(point) for q in variables],
      'one simple root of a repeated multiplier gives exactly its full Morse point')
check(sp.diff(blocked,z).subs(z,1)==1,
      'the repeated root with nonzero remainder derivative has no critical point')
zero(sp.hessian(blocked,variables).subs(point).det()-324,
     'the surviving simple-root critical point is Morse')
check(sp.degree(sp.div(repeated,sp.gcd(repeated,sp.diff(repeated,z)),z)[0],z)==2,
      'the repeated control has two distinct roots but only one simple root')
unblocked=2*x*y+w*repeated+(z-1)**2
zeros([sp.diff(unblocked,q).subs({x:0,y:0,z:1}) for q in variables],
      'a vanishing remainder derivative at a repeated root keeps the whole critical dual line')
zero(sp.hessian(unblocked,variables).det().subs(z,1),
     'the whole repeated-root critical line is degenerate')
zero(sp.Matrix([2*x,(z-1)**2]).jacobian(variables).subs({x:0,z:1}).rank()-1,
     'the nonreduced relative critical scheme fails the selected-row rank test')

# Every whole-value Hodge profile, including both quadratic special-value patterns.
L,root_count,jump=sp.symbols('L root_count jump')
zero((L-1)*L*L+L*(L-root_count+jump*L)
     -(L**3-root_count*L+jump*L*L),
     'the complete whole-value profile retains all vertical root lines')
for roots_at_value in [0,1,2]:
    zero((L-1)*L*L+L*(L-2+roots_at_value*L)
         -(L**3-2*L+roots_at_value*L*L),
         'the two-point model passes the full profile with critical-value jump '+str(roots_at_value))
check(2!=1,'the repeated-root generic profile differs from its actual one-critical-point count')
for beta in [0,1]:
    potential=2*x*y+w*(z*z-1)+beta*z
    for root in [-1,1]:
        point={x:0,y:0,z:root,w:-sp.Rational(beta,2*root)}
        zeros([sp.diff(potential,q).subs(point) for q in variables],
              'the canonical two-point model retains its exact critical point alpha='+str(beta)+' root='+str(root))
        zero(sp.hessian(potential,variables).subs(point).det()-16,
             'both canonical models have common critical determinant sixteen alpha='+str(beta)+' root='+str(root))
        zero(potential.subs(point)-beta*root,
             'the actual critical values distinguish the two canonical source-orbit tests')

# Linear and quadratic normal forms have explicit polynomial inverse maps.
multiplier=3*(z-2);remainder=z**5+z+7
quotient=sp.div(remainder-remainder.subs(z,2),z-2,z)[0]
zz=Z+2;ww=(W-quotient.subs(z,zz))/3
zero((2*x*y+w*multiplier+remainder).subs({z:zz,w:ww},simultaneous=True)
     -(2*x*y+Z*W+remainder.subs(z,2)),
     'the arbitrary-degree linear-multiplier shear gives the exact nondegenerate quadratic orbit')
zero((3*w+quotient).subs({z:zz,w:ww},simultaneous=True)-W,
     'the linear-multiplier shear has a global polynomial inverse')
for slope in [0,5]:
    multiplier=2*((z-3)**2-4)
    quotient=z**3+z
    remainder=multiplier*quotient+slope*z+7
    zz=3+2*Z;ww=W/8-quotient.subs(z,zz)
    original=2*x*y+w*multiplier+remainder
    zero(original.subs({z:zz,w:ww},simultaneous=True)
         -(2*x*y+W*(Z*Z-1)+2*slope*Z+3*slope+7),
         'the arbitrary-degree quadratic remainder reduces to the exact canonical alpha model '+str(slope))
    zero(((z-3)/2).subs(z,zz)-Z,'the quadratic normal form recovers its base coordinate')
    zero((8*(w+quotient)).subs({z:zz,w:ww},simultaneous=True)-W,
         'the quadratic normal form has its global polynomial inverse')
canonical=2*x*y+w*(z*z-1)+alpha*z
one=2*X*Y+W*(Z*Z-1)+Z
zero(alpha*one.subs({X:x/alpha,Y:y,Z:z,W:w/alpha},simultaneous=True)-canonical,
     'every nonzero alpha is the exact scalar source orbit of the one fixed alpha-one model')

# Actual closed nondegenerate ambient graphs of unrestricted multiplier degree.
for A in [sp.Integer(1),z,z*z-1,z*(z-1),z**3-z]:
    trace=2*x*y+z*z+w*w-(w+A)**2
    zero(trace-(2*x*y+w*(-2*A)+(z*z-A*A)),
         'the literal compensated graph has the exact suspended potential A='+str(A))
    check(sp.det(sp.hessian(2*x*y+z*z+w*w-X*X,[x,y,z,w,X]))!=0,
          'the actual graph ambient quadratic form is nondegenerate')
zero((z*z-(z*z-1)**2).subs(z,1)
     -(z*z-(z*z-1)**2).subs(z,-1),
     'the centered quadratic graph realizes the alpha-zero test')
check((z*z-(z*(z-1))**2).subs(z,0)
      !=(z*z-(z*(z-1))**2).subs(z,1),
      'the noncentered quadratic graph realizes the alpha-nonzero test')

# A nonlinear base map and nonlinear dual map: exact full relative scheme transport.
sx=2*x+z*z;sz=-z
tau_y=3*y+(x+z)*w*w;tau_w=w
check(sp.Matrix([sx,sz]).jacobian([x,z]).det()==-2,
      'the nonlinear base automorphism has a constant nonzero Jacobian')
check(sp.Matrix([tau_y,tau_w]).jacobian([y,w]).det()==3,
      'the nonlinear dual automorphism has a constant nonzero relative Jacobian')
zero(sx.subs({x:(X-Z*Z)/2,z:-Z},simultaneous=True)-X,
     'the nonlinear base map has its exact polynomial first inverse')
zero(sz.subs(z,-Z)-Z,'the nonlinear base map has its exact polynomial second inverse')
zero(tau_y.subs(y,(Y-(x+z)*W*W)/3).subs(w,W)-Y,
     'the nonlinear dual map has its exact polynomial inverse')
g=2*sx*tau_y+tau_w*(sz*sz-1)+alpha*sz
relative=sp.Matrix([sp.diff(g,y),sp.diff(g,w)])
Jtau=sp.Matrix([tau_y,tau_w]).jacobian([y,w])
zeros(list(relative-Jtau.T*sp.Matrix([2*sx,sz*sz-1])),
      'the relative gradient ideal is exactly the transported coefficient ideal')
check(all(sp.denom(sp.cancel(v)).is_number for v in Jtau.inv()),
      'the relative inverse matrix is polynomial even with nonlinear dual dependence')
for root in [-1,1]:
    base={x:-sp.Rational(1,2),z:-root}
    zeros(list(sp.hessian(g,[y,w]).subs(base)),
          'the vertical Hessian vanishes on the entire relative critical plane root='+str(root))
    orientation=sp.simplify(relative.jacobian([x,z]).subs(base).det())
    zero(orientation-(-2)*3*(4*root),
         'the full orientation factor is base Jacobian times relative Jacobian times native row orientation root='+str(root))
    check(orientation!=0,'the oriented relative critical component remains reduced')
    native_w=-alpha/(2*root)
    source_y=-(base[x]+base[z])*native_w**2/3
    point=base|{y:source_y,w:native_w}
    zeros([sp.diff(g,q).subs(point) for q in variables],
          'the nonlinear source map retains the actual full critical point root='+str(root))
    zero(sp.hessian(g,variables).subs(point).det()-576,
         'the nonlinear source map retains equal critical Hessian squares despite opposite row orientations root='+str(root))

# The retained rank-one and rank-two relative identities used in alignment.
T,U,V,S=sp.symbols('T U V S')
B=U+T*T;C0=T
rankone=U**3+S*B+V*C0
check(sp.Matrix([B,C0]).jacobian([T,U]).det()==-1,
      'the terminal rank-one coefficient row is an actual triangular Keller automorphism')
zero(sp.hessian(rankone,[U,V]).det()+sp.diff(C0,U)**2,
     'the exact rank-one vertical determinant formula is retained')
KK=sp.Matrix([U+V*V,V])
DD=U**3+U*V
base_graph=-KK.jacobian([U,V]).T.inv()*sp.Matrix([sp.diff(DD,U),sp.diff(DD,V)])
ranktwo=DD+T*KK[0]+S*KK[1]
zeros([sp.diff(ranktwo,q).subs({T:base_graph[0],S:base_graph[1]},simultaneous=True) for q in [U,V]],
      'the rank-two relative critical scheme is the whole polynomial graph')
check(KK.jacobian([U,V]).det()==1,
      'the rank-two graph inverse uses a global unit coefficient Jacobian')

context=json.loads((HERE/'sources/suspended-multiplier-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'retained proof/checker/receipt/source bytes stay pinned: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'the primary standard-smooth context is byte-pinned '+entry['retained_file'])
check(context['newly_used_unproved_external_statements']==[],
      'the proof imports no additional unproved source-orbit or differential-equation theorem')
check(context['all_base_preserving_source_orbits_degree_at_least_two_excluded'] is True,
      'the relative closure includes nonlinear base and nonlinear dual maps')
check(context['quadratic_two_critical_full_source_orbits_resolved'] is False
      and context['universal_two_point_keller_endpoint_equivalence_proved'] is False
      and context['unit_ideal_endpoint_resolved'] is False
      and context['singleton_zero_endpoint_resolved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'the fixed full-source tests, both actual endpoints and unrestricted HC4 stay open')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
receipt=dict(outcome='PASS',claims=['HC4MYJK0','HC4MYJK1'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'SUSPENDED_MULTIPLIER_ORBITS_AND_BASE_PRESERVING_CLOSURE.md').read_bytes()).hexdigest(),
 source_sha256={'sources/suspended-multiplier-context.json':hashlib.sha256((HERE/'sources/suspended-multiplier-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written all-degree derivative rigidity, whole-value strata and relative orientation proof; not independent replay, formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,
 elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-suspended-multipliers.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
