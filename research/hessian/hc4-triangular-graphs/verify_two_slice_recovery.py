"""Exact controls for the written all-degree recovery and orientation arguments."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(value)==0 or sp.simplify(value)==0,label)
def zeros(values,label):check(all(sp.cancel(v)==0 or sp.simplify(v)==0 for v in values),label)
x,y,z,w,h,P,Q,U,V=sp.symbols('x y z w h P Q U V');coords=[x,z,y,w]
F=x;G=z*z-1
PP=y+w*w;QQ=w;hh=F*PP+G*QQ
delta=[0,0,G+2*x*w,-x]
def deriv(f):return sp.expand(sum(a*sp.diff(f,b) for a,b in zip(delta,coords)))
zero(deriv(PP)-G,'the first global slice has derivative G')
zero(deriv(QQ)+F,'the second global slice has derivative minus F')
zero(deriv(hh),'the value generator is an actual invariant')
alpha=hh*hh+x;beta=z*hh
p=PP+alpha;q=QQ+beta;a=sp.expand(F*p+G*q)
zero(deriv(p)-G,'an invariant correction preserves the first slice')
zero(deriv(q)+F,'an invariant correction preserves the second slice')
zero(a-hh-F*alpha-G*beta,'the exact affine-socle identity has explicit ideal witnesses')
zero(F*(p-alpha)+G*(q-beta)-hh,'the corrected global slices recover the literal value identity')
for sign in [-1,1]:
    zero(a.subs({x:0,z:sign})-hh.subs({x:0,z:sign}),
         'the socle restriction is independent of the chosen global slices,point='+str(sign))
zero(p/G-(-q/F)-a/(F*G),'the two localization charts have exactly the recorded torsor transition')
source_inverse={x:x,z:z,y:P-Q*Q,w:Q}
zeros([PP.subs(source_inverse,simultaneous=True)-P,QQ.subs(source_inverse,simultaneous=True)-Q],
      'the recovered coordinate frame has actual polynomial reverse identities')
zero((P-Q*Q).subs({P:PP,Q:QQ},simultaneous=True)-y,'the native y inverse identity holds as a polynomial')
check(sp.Matrix([x,z,PP,QQ]).jacobian(coords).det()==1,'the complete recovered coordinate map has constant nonzero Jacobian')
for m,n in [(1,2),(3,1)]:
    native=P-Q*Q
    numerator_F=sp.expand(F**m*native);numerator_G=sp.expand(G**n*native)
    zero(G**n*numerator_F-F**m*numerator_G,'the two chart numerators agree in the common fraction field,m='+str(m))
    quotient,remainder=sp.div(numerator_F,F**m,x,z,P,Q)
    zero(remainder,'the chart denominator cancels by exact polynomial division,m='+str(m))
    zero(quotient-native,'the cancelled chart produces the actual inverse coordinate,m='+str(m))
# A repeated root forces the explicit hypersurface singularity; a simple root
# has a unit h derivative. The global normality argument is written in the note.
for power in [1,2,3]:
    relation=x*P+(z*z-1)*Q-(h-z)**power
    for sign in [-1,1]:
        point={x:0,z:sign,h:sign,P:0,Q:0}
        zero(relation.subs(point),'the point lies on the reconstructed hypersurface,m='+str(power)+',point='+str(sign))
        gradient=[sp.diff(relation,b).subs(point) for b in [x,z,h,P,Q]]
        check(all(v==0 for v in gradient)==(power>1),
              'precisely repeated residue roots make the source singular,m='+str(power)+',point='+str(sign))
lambda_poly=z;row=sp.Matrix([[x,z*z-1,-lambda_poly]])
e=sp.Matrix([0,-1,-z]);v1=sp.Matrix([1,x,x*z]);v2=sp.Matrix([0,-z,1-z*z])
M=sp.Matrix.hstack(e,v1,v2)
check(M.det()==1,'a genuinely nonconstant unimodular row has a global polynomial completion')
zeros(row*M-sp.Matrix([[1,0,0]]),'the completed columns split the row into its section and kernel')
zeros(v1.cross(v2)-row.T,'the cross product retains the exact common orientation factor')
dd=z*z+1;param=dd*e+U*v1+V*v2
zero((row*param)[0]-dd,'the section plus free kernel basis solves the entire affine hyperplane')
recovered=M.inv()*param
zeros(recovered-sp.Matrix([dd,U,V]),'the projective-kernel frame has full polynomial inverse identities')
f=param[2].coeff(U);g=param[2].coeff(V)
Jfg=sp.Matrix([f,g]).jacobian([x,z]).det();JFG=sp.Matrix([F,G]).jacobian([x,z]).det()
for sign in [-1,1]:
    zero((Jfg+M.det()*JFG/lambda_poly).subs({x:0,z:sign}),
         'the recovered coefficient orientation is -mu times J(F,G)/lambda_i,point='+str(sign))
    zero((lambda_poly/JFG).subs({x:0,z:sign})-sp.Rational(1,2),
         'opposite slice slopes can still have one common oriented residue when the original Jacobians also change sign,point='+str(sign))
# The determinant repair is symbolic and has no degree or support bound.
a0,b0,c0,d0,rr,ss,f1,f2=sp.symbols('a0 b0 c0 d0 rr ss f1 f2')
initial=sp.Matrix([[a0,b0],[c0,d0]]);ff=sp.Matrix([f1,f2]);kk=initial*ff
repair=initial+sp.Matrix([rr,ss])*sp.Matrix([[-f2,f1]])
zeros(repair*ff-kk,'adding the polynomial syzygy preserves both ordered generators')
zero(repair.det()-initial.det()-ss*kk[0]+rr*kk[1],
     'the determinant correction is exactly s k1-r k2')
initial=sp.diag(1+x*G,1);kk=initial*sp.Matrix([F,G])
repair=initial+sp.Matrix([x,0])*sp.Matrix([[-G,F]])
check(repair.det()==1,'a globally noninvertible presentation matrix is repaired to a true GL2 matrix')
zeros(repair*sp.Matrix([F,G])-kk,'the actual determinant repair preserves the nonlinear pair')
check(sp.Poly(initial.det(),x,z).total_degree()>0,'the repair is needed globally,not merely at the finite zero set')
bad_pair=sp.Matrix([z*x,G]);badJ=bad_pair.jacobian([x,z]).det()
ratios=[sp.cancel((badJ/JFG).subs({x:0,z:sign})) for sign in [-1,1]]
check(ratios==[-1,1],'equal reduced ideals can carry incompatible unsquared orientations')
check([r*r for r in ratios]==[1,1],'squaring the orientation loses the exact obstruction')
zero(z*(z*x)-x*G-x,'the differently oriented pair nevertheless generates the very same finite-point ideal')
oldpair=sp.Matrix([2*x,G]);oldJ=oldpair.jacobian([x,z]).det()
check([oldJ.subs(z,sign) for sign in [-1,1]]==[-4,4],
      'the old h0/h1 natural-base obstruction is the mixed unsquared orientation')
# General coefficient functions give an exact all-degree Hessian identity.
af=sp.Function('a')(z);Cf=sp.Function('C')(z);Df=sp.Function('d')(x,z)
potential=af*x*P+Cf*Q+Df
zero(sp.hessian(potential,[x,z,P,Q]).det()-(af*sp.diff(Cf,z))**2,
     'every coefficient degree has the same literal mixed-block Hessian determinant')
ell=sp.symbols('ell',nonzero=True)
K1=x+z*z;K2=z
unit_potential=(P*K1+Q*K2-Df)/ell
zero(sp.hessian(unit_potential,[x,z,P,Q]).det()-ell**-4,
     'the elementary two-slice certificate gives the actual unit full Hessian')
# Controls for the general inverse-modulo-C construction,not a degree search.
for N in [2,3]:
    C=sp.prod(z-i for i in range(1,N+1));aa=sp.invert(sp.diff(C,z),C,z)
    zero(sp.rem(aa*sp.diff(C,z)-1,C,z),'the canonical coefficient is the exact inverse derivative modulo C,N='+str(N))
    check(sp.gcd(aa,C,z)==1,'the canonical pair generates the radical graph ideal,N='+str(N))
    T=aa*x*P+C*Q
    zero(sp.hessian(T,[x,z,P,Q]).det()-(aa*sp.diff(C,z))**2,
         'the canonical model uses its literal native Hessian,N='+str(N))
    check(sp.degree(aa*sp.diff(C,z),z)>0,'point orientation one is not a globally constant native coefficient Jacobian,N='+str(N))
    for i in range(1,N+1):
        point={x:0,z:i,P:0,Q:0}
        zeros([sp.diff(T,b).subs(point) for b in [x,z,P,Q]],'each prescribed base zero is a full critical point,N='+str(N)+',point='+str(i))
        zero(sp.hessian(T,[x,z,P,Q]).det().subs(point)-1,'all canonical full critical Hessians equal one,N='+str(N)+',point='+str(i))
        zero((aa*sp.diff(C,z)).subs(z,i)-1,'all canonical unsquared coefficient orientations equal one,N='+str(N)+',point='+str(i))
# A nonaffine finite-point matching is realized by three exact shear inverses.
points=[(0,2),(0,5),(0,9)];targets=[1,2,3]
shear1=sp.interpolate([(v,i) for (_,v),i in zip(points,targets)],z)
shear2=sp.interpolate([(i,i-v) for (_,v),i in zip(points,targets)],x)
for (_,oldv),i in zip(points,targets):
    xx=shear1.subs(z,oldv);zz=oldv+shear2.subs(x,xx);xx=xx-zz
    check((xx,zz)==(0,i),'the three actual interpolating shears match arbitrary distinct points,point='+str(i))
zero((x+shear1).subs(x,x-shear1)-x,'the first polynomial shear has its exact inverse')
zero((z+shear2).subs(z,z-shear2)-z,'the second polynomial shear has its exact inverse')
context=json.loads((HERE/'sources/two-slice-recovery-context.json').read_text());algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,'retained source/checker/receipt bytes: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'inspected primary algebraic import: '+entry['retained_file'])
check(len(context['newly_imported_proved_external_statements'])==5 and all(e['hypotheses_checked'] for e in context['newly_imported_proved_external_statements']),
      'all five proved imports have explicit application hypotheses')
check(context['newly_used_unproved_external_statements']==[],'no action classification or unproved cancellation is imported')
check(context['elementary_two_slice_affine_socle_certificate_equivalent_to_unit_hessian_source']
      and not context['quotient_type_or_annihilator_alone_assumed_sufficient']
      and not context['global_Keller_generators_for_any_N_ge_2_constructed']
      and not context['auxiliary_T_N_assumed_in_original_affine_graph_frontier']
      and not context['unrestricted_HC4_or_JC2_resolved'],
      'coordinate recovery is closed under its exact hypotheses while action existence and global Keller generators remain open')
flags={k:v for k,v in context.items() if isinstance(v,bool)}
sources={'sources/two-slice-recovery-context.json':hashlib.sha256((HERE/'sources/two-slice-recovery-context.json').read_bytes()).hexdigest()}
sources.update({'sources/'+e['retained_file']:e['sha256'] for e in context['newly_cached_primary_context']})
receipt=dict(outcome='PASS',claims=['HC4MYJV0','HC4MYJV1','HC4MYJV2'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'TWO_SLICE_RECOVERY_AND_ORIENTED_COEFFICIENT_MODELS.md').read_bytes()).hexdigest(),
 source_sha256=sources,dependencies=context['retained_proof_and_source_sha256'],
 proved_import_count=5,assurance='Finite exact controls of written all-degree two-slice recovery,normal hypersurface,projective-module and oriented-pair arguments; not independent replay,formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-two-slice-recovery.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['outcome','check_count','algebraic_and_scope_controls','proved_import_count','elapsed_seconds','original_goal_complete']},indent=2))
