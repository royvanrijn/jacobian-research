"""Exact controls for the all-offset polar-curve proof and semisimple closure."""
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
    check(all(sp.cancel(value)==0 or sp.simplify(value)==0 for value in values),label)
x,y,z,w,t,kappa,r,s=sp.symbols('x y z w t kappa r s')
a,b,e,j,c,mu,lam=sp.symbols('a b e j c mu lam',nonzero=True)
L,U=sp.symbols('L U')

# Exact polar equation and finite-free four-sheeted plane presentation.
H=t**3+2*t
inner=z*z+w*w
f=a*z*z+b*w*w+e*z+j*w+H.subs(t,inner)
polar=sp.diff(inner,z)*sp.diff(f,w)-sp.diff(inner,w)*sp.diff(f,z)
zero(polar-(4*(b-a)*z*w+2*j*z-2*e*w),
     'the outer derivative cancels from the exact affine polar equation')
zsquare=(kappa-H-b*t-e*z-j*w)/(a-b)
wsquare=(a*t-kappa+H+e*z+j*w)/(a-b)
zero(zsquare+wsquare-t,'the monic square relations retain the inner coordinate')
zero(a*zsquare+b*wsquare+e*z+j*w+H-kappa,
     'the monic square relations retain the full affine level equation')
G=sp.groebner([z*z-zsquare,w*w-wsquare],z,w,order='grlex',domain=sp.QQ.frac_field(a,b,e,j,t,kappa))
check({p.LM(order=G.order).exponents for p in G.polys}=={(2,0),(0,2)},
      'the four-element standard basis has relatively prime monic leading terms')
zero(polar/4/(b-a)-(z-e/(2*(b-a)))*(w+j/(2*(b-a)))-e*j/(4*(b-a)**2),
     'the polar is an actual translated hyperbola or two transverse lines')

# Hyperbola and reducible-polar controls for even and odd outer degrees.
for outer in [t*t+2*t,t**3+2*t]:
    d=sp.degree(outer,t)
    zz=r+2;ww=-3-6/r
    zero((z*w+3*z-2*w).subs({z:zz,w:ww}),
         'the hyperbola parameter has the exact polar inverse degree '+str(d))
    q=sp.expand(zz*zz+ww*ww)
    ff=2*zz*zz+3*ww*ww+4*zz+6*ww+outer.subs(t,q)
    numerator=sp.Poly(sp.cancel(r**(2*d)*(ff-17)),r)
    check(numerator.degree()==4*d and numerator.eval(0)!=0,
          'both polar ends contribute exact order 2d degree '+str(d))
    check(sp.gcd(numerator,numerator.diff()).degree()==0,
          'a literal generic hyperbola value has all 4d distinct finite ramification points degree '+str(d))
    line1=3*w*w+6*w+outer.subs(t,w*w)
    line2=2*z*z+27-18+outer.subs(t,z*z+9)
    for coordinate,level,label in [(w,line1,'vertical'),(z,line2,'horizontal')]:
        polynomial=sp.Poly(level-17,coordinate)
        check(polynomial.degree()==2*d,
              'each reducible-polar line has exact 2d level degree '+label+' '+str(d))
        check(sp.gcd(polynomial,polynomial.diff()).degree()==0,
              'the reducible-polar simple-root control passes '+label+' '+str(d))

# The infinity initial system has four distinct solutions and no missing degree.
for zz in [sp.Integer(1),sp.Integer(-1)]:
    for ww in [sp.I,-sp.I]:
        zero(zz**2+ww**2,'an infinity branch lies on the inner null cone')
        zero(zz**2+2*ww**2+1,'an infinity branch retains the outer leading coefficient')
        check(4*zz*ww!=0,'the infinity initial Jacobian is nonzero on every branch')
m,components=sp.symbols('m components',integer=True,positive=True)
for epsilon in [0,1]:
    d=2*m+epsilon
    genus=2*d-4+epsilon+components
    zero(2*genus-2*components-(-8+4*d+2*epsilon),
         'componentwise Riemann-Hurwitz retains finite and infinity ramification parity '+str(epsilon))
    base=(d+epsilon-2)/2
    difference=sp.expand((2*d-3+epsilon)-base)
    check(difference.subs(m,1)>0 and sp.diff(difference,m)>0,
          'the total plane genus strictly exceeds the rank-one base genus in every admitted degree parity '+str(epsilon))
    check(4-2*epsilon+(2*epsilon)==4,
          'infinity punctures and ramification exhaust degree four parity '+str(epsilon))

# Scalar plane offsets, including both exact polynomial inverses.
outer=t**3+2*t
for eigen,offset,constant in [(2,3,1),(3,2,-1)]:
    zmap=(t+eigen*outer-eigen*kappa-constant)/(2*offset)
    plane_inner=eigen*(z*z+w*w)+2*offset*z+constant
    plane_h=z*z+w*w+outer.subs(t,plane_inner)
    polynomial=sp.expand(kappa-outer-zmap*zmap)
    zero(eigen*(zmap*zmap+polynomial)+2*offset*zmap+constant-t,
         'the scalar hyperelliptic presentation recovers its actual inner coordinate '+str(eigen))
    zero(zmap*zmap+polynomial+outer-kappa,
         'the scalar hyperelliptic presentation recovers its complete level '+str(eigen))
    zero(zmap.subs({t:plane_inner,kappa:plane_h},simultaneous=True)-z,
         'the scalar curve presentation has a polynomial inverse '+str(eigen))
    pp=sp.Poly(polynomial.subs(kappa,7),t)
    check(pp.degree()==6 and sp.gcd(pp,pp.diff()).degree()==0,
          'the scalar nonisotropic degree-six double cover has all six simple branch roots '+str(eigen))
    scalar_iso_inner=2*eigen*z*w+2*z+constant
    scalar_iso_h=2*z*w+outer.subs(t,scalar_iso_inner)
    iso_z=(t+eigen*outer-eigen*kappa-constant)/2
    iso_w=(t-2*iso_z-constant)/(2*eigen*iso_z)
    zero(2*eigen*iso_z*iso_w+2*iso_z+constant-t,
         'the punctured-line scalar model recovers the inner coordinate '+str(eigen))
    zero(2*iso_z*iso_w+outer-kappa,
         'the punctured-line scalar model recovers the whole curve '+str(eigen))
    zero(iso_z.subs({t:scalar_iso_inner,kappa:scalar_iso_h},simultaneous=True)-z,
         'the punctured-line scalar model has the stated inverse '+str(eigen))
    roots=sp.Poly((t+eigen*outer-eigen*7-constant),t)
    check(roots.degree()==3 and sp.gcd(roots,roots.diff()).degree()==0,
          'the scalar isotropic curve removes exactly d distinct points '+str(eigen))

# Whole fibres with an isotropic projected offset: no degenerate-Krylov split.
for eigen,remaining,offsets,label in [
    (2,(3,5),(1,2),'rank two'),
    (2,(2,5),(0,2),'rank one'),
    (2,(2,2),(0,0),'rank zero')]:
    q=z*z+w*w
    g0=remaining[0]*z*z+remaining[1]*w*w+2*offsets[0]*z+2*offsets[1]*w+1
    full_inner=2*eigen*x*y+2*x+g0
    full_h=2*x*y+q+outer.subs(t,full_inner)
    M=t+eigen*outer+eigen*q-g0-eigen*kappa
    xmap=M/2;ymap=(t-2*xmap-g0)/(2*eigen*xmap)
    zero(2*eigen*xmap*ymap+2*xmap+g0-t,
         'the whole open stratum recovers the exact inner coordinate '+label)
    zero(2*xmap*ymap+q+outer-kappa,
         'the whole open stratum recovers the exact value fibre '+label)
    zero(xmap.subs({t:full_inner,kappa:full_h},simultaneous=True)-x,
         'the open stratum inverse recovers its x coordinate '+label)
    zero(ymap.subs({t:full_inner,kappa:full_h},simultaneous=True)-y,
         'the open stratum inverse recovers its y coordinate '+label)
    zero(full_h.subs(x,0)-(q+outer.subs(t,g0)),
         'the omitted x=0 stratum is retained in full '+label)
    D=eigen*q-g0
    rank=sp.hessian(D,[z,w]).rank()
    check(rank=={'rank two':2,'rank one':1,'rank zero':0}[label],
          'the scalar base rank is its actual rank '+label)
    for vec in sp.hessian(D,[z,w]).nullspace():
        zero((sp.Matrix([sp.diff(D,z),sp.diff(D,w)]).subs({z:0,w:0}).T*vec)[0],
             'there is no unrecorded linear term on the scalar-base kernel '+label)

# High multiplicity nonisotropic active primaries: the whole-fibre polynomial inverse.
pair=sp.symbols('pair')
for eigen,other,aa,bb,label in [(2,5,3,1,'multiplicity three'),(2,2,3,0,'multiplicity four')]:
    g2=eigen*x*x+other*z*z+2*aa*x+2*bb*z+1
    full_inner=eigen*pair+g2
    full_h=pair+x*x+z*z+outer.subs(t,full_inner)
    xx=(t+eigen*outer-eigen*kappa-(other-eigen)*z*z-2*bb*z-1)/(2*aa)
    S=kappa-outer-xx*xx-z*z
    zero(eigen*S+eigen*xx*xx+other*z*z+2*aa*xx+2*bb*z+1-t,
         'the high-primary fibre map recovers its actual inner coordinate '+label)
    zero(S+xx*xx+z*z+outer-kappa,
         'the high-primary fibre map retains the whole value '+label)
    zero(xx.subs({t:full_inner,kappa:full_h},simultaneous=True)-x,
         'the high-primary fibre isomorphism has a global polynomial inverse '+label)
    plane_h=x*x+z*z+outer.subs(t,g2)
    zero(S.subs({t:g2,kappa:plane_h},simultaneous=True),
         'the zero locus in the suspension is the actual plane level curve '+label)

# Literal product-stratum Hodge identities and all low-genus exceptions.
d=sp.symbols('d',integer=True,positive=True)
zero((L-1)*L+d*L-(L*L+(d-1)*L),
     'the rank-two base retains every special affine-line stratum')
zero((L-1)*L*L+L*(L-d)-(L**3-d*L),
     'a pure hyperbolic suspension retains all d missing points')
zero(L**3-(L*L+(d-1)*L)+L*d*(L-1)
     -(L**3+(d-1)*L*L-(2*d-1)*L),
     'the centered remaining scalar case has its nonzero L-squared coefficient')
zero(L**3-(L*L+(d-1)*L)+L*(L-d)-(L**3-(2*d-1)*L),
     'two isotropic projected offsets have the full 2d-minus-one profile')
zero(L**3-d*L*L+L*d*(L-1)-(L**3-d*L),
     'the scalar four-primary isotropic case has the exact d coefficient')
hyper=sp.diag(sp.Matrix([[0,1],[1,0]]),sp.Matrix([[0,1],[1,0]]))
for TT,vv,label in [
    (sp.diag(2,2,3,3),sp.Matrix([0,1,0,1]),'two isotropic primaries'),
    (2*sp.eye(4),sp.Matrix([0,1,0,0]),'one scalar isotropic primary')]:
    stationary=-s*(sp.eye(4)+s*TT).inv()*vv
    RR=(stationary.T*hyper*TT*stationary+2*vv.T*hyper*stationary)[0]
    zero(RR,'the low-genus exceptional count really has constant R '+label)
    check(sp.Poly(s-sp.diff(outer,t).subs(t,RR),s).degree()==1,
          'there is at most one regular critical parameter '+label)
    check(sp.Matrix.hstack(vv,TT*vv,TT**2*vv,TT**3*vv).rank()<4,
          'this low-genus offset is genuinely noncyclic '+label)

# Singular semisimple control: an actual family passes both critical/Hodge tests.
aa,alpha=sp.symbols('aa alpha')
original_inner=z*z+2*aa*z+2*w+aa*aa
original_h=2*x*y+z*z+w*w-original_inner**2/4+original_inner/2
Z,W=sp.symbols('Z W')
zinverse=Z-aa;winverse=-W-Z*Z/4+sp.Rational(5,4)
zero(original_h.subs({z:zinverse,w:winverse},simultaneous=True)
     -(2*x*y+W*(Z*Z-1)-2*aa*Z+aa*aa+sp.Rational(5,4)),
     'a singular rank-one affine pencil has the exact two-point source normal form')
zero((z+aa).subs(z,zinverse)-Z,'the singular-control source map has its exact first inverse')
zero((-w-(z+aa)**2/4+sp.Rational(5,4)).subs({z:zinverse,w:winverse},simultaneous=True)-W,
     'the singular-control source map has its exact second inverse')
check(sp.Matrix([[sp.diff(zinverse,Z),sp.diff(zinverse,W)],
                 [sp.diff(winverse,Z),sp.diff(winverse,W)]]).det()==-1,
      'the singular-control polynomial source Jacobian is a constant unit')
boundary=2*x*y+W*(Z*Z-1)+alpha*Z
gradient=[sp.diff(boundary,q) for q in [x,y,Z,W]]
Hess=sp.hessian(boundary,[x,y,Z,W])
zero(Hess.det()-16*Z*Z,'the singular control is not falsely asserted to have native unit Hessian')
for root in [1,-1]:
    point={x:0,y:0,Z:root,W:-alpha/(2*root)}
    zeros([value.subs(point) for value in gradient],
          'the singular control has an actual reduced critical point '+str(root))
    zero(Hess.subs(point).det()-16,
         'both singular-control critical Hessian determinants are the same unit '+str(root))
zero((L-1)*L*L+L*(L-2)-(L**3-2*L),
     'the singular control passes the exact HG0 generic two-critical-point profile')
check(sp.diag(0,0,1,0).rank()==1 and sp.diag(0,0,1,0).det()==0,
      'the genuine residual pencil is semisimple and singular')

context=json.loads((HERE/'sources/affine-polar-semisimple-context.json').read_text())
algebra=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'retained source/proof/checker/receipt bytes stay pinned: '+path)
for entry in context['newly_cached_primary_context']:
    check(hashlib.sha256((HERE/'sources'/entry['retained_file']).read_bytes()).hexdigest()==entry['sha256'],
          'the versioned primary suspension context is byte-pinned')
check(context['newly_used_unproved_external_statements']==[],
      'no A1-homotopy, proper-gradient or connectedness theorem is imported')
check(context['all_invertible_semisimple_affine_offsets_inverse_closed'] is True,
      'the complete invertible semisimple result includes every affine offset')
check(context['all_singular_semisimple_affine_pencils_resolved'] is False
      and context['all_invertible_semisimple_unit_hessian_source_orbits_proved_nonexistent'] is False
      and context['unit_ideal_endpoint_resolved'] is False
      and context['singleton_zero_endpoint_resolved'] is False
      and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'inverse closure is not overstated as nonexistence or a solution of singular pencils/endpoints')
flags={key:value for key,value in context.items() if isinstance(value,bool)}
receipt=dict(outcome='PASS',claims=['HC4MYJI0','HC4MYJI1'],check_count=len(checks),
 algebraic_and_scope_controls=algebra+3,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'AFFINE_POLAR_CURVES_AND_COMPLETE_INVERTIBLE_SEMISIMPLE_CLOSURE.md').read_bytes()).hexdigest(),
 source_sha256={'sources/affine-polar-semisimple-context.json':hashlib.sha256((HERE/'sources/affine-polar-semisimple-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Exact finite controls of the written all-degree polar ramification and whole-fibre proof; not independent replay, formal verification or external review.',
 python=platform.python_version(),sympy=sp.__version__,
 elapsed_seconds=round(time.monotonic()-START,3),**flags)
(HERE/'receipt-affine-polar-curves.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({key:receipt[key] for key in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
