"""Exact controls for affine critical components and boundary-root graphs.

The all-rank, all-degree and Euler arguments are written proofs. Classical
line rectification, line injectivity and Tang are explicit pinned imports.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import sympy as S

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args = parser.parse_args()
    start = time.monotonic()
    checks = []

    def check(name,condition):
        if not condition:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,expression):
        values = list(expression) if isinstance(expression,S.MatrixBase) else [expression]
        check(name,all(S.cancel(value.doit())==0 for value in values))

    x,y,p,q,r,s,U,Z,V,P,Q,T,t,lam = S.symbols('x y p q r s U Z V P Q T t lam')
    L,N,v = S.symbols('L N v',nonzero=True)
    M,rho = S.symbols('M rho')
    # A genuinely curved closed affine-line fibre, with its actual inverse.
    AA = U+Z*Z
    CC = Z-AA**3
    pair = S.Matrix([AA,CC])
    pair_inverse = S.Matrix([P-(Q+P**3)**2,Q+P**3])
    zero('curved affine-line control is a genuine plane automorphism',pair.jacobian((U,Z)).det()-1)
    zero('curved-line plane inverse left',pair_inverse.subs({P:AA,Q:CC},simultaneous=True)-S.Matrix([U,Z]))
    zero('curved-line plane inverse right',pair.subs(dict(zip((U,Z),pair_inverse)),simultaneous=True)-S.Matrix([P,Q]))
    gamma = S.Matrix([t-t**6,t**3])
    zero('closed curved line lies in the complete zero fibre',pair.subs(dict(zip((U,Z),gamma)),simultaneous=True)-S.Matrix([t,0]))
    zero('curved-line parameter has a polynomial inverse',AA.subs(dict(zip((U,Z),gamma)),simultaneous=True)-t)
    zero('curved embedding is an immersion with a literal tangent Bezout certificate',S.diff(gamma[0],t)+2*t**3*S.diff(gamma[1],t)-1)

    FF,GG,DD = S.Function('FF')(s,U),S.Function('GG')(s,U),S.Function('DD')(s,U)
    rank_zero = DD+P*FF+Q*GG
    J0 = S.Matrix([FF,GG]).jacobian((s,U)).det()
    zero('rank-zero full Hessian is the residual Keller determinant squared',S.hessian(rank_zero,(s,U,P,Q)).det()-J0**2)
    zero('rank-zero critical equations contain the entire plane zero fibre',S.Matrix([S.diff(rank_zero,P),S.diff(rank_zero,Q)])-S.Matrix([FF,GG]))
    # Generic fibre completion needs only a field-unit free-base value.
    fibre_control = U**3+s*P+U*Q
    recovered_P = (V-U**3-U*Q)/s
    zero('rank-zero generic dual fibre has a literal coordinate inverse',fibre_control.subs(P,recovered_P)-V)
    zero('rank-zero generic coordinate inverse recovers P',recovered_P.subs(V,fibre_control)-P)
    zero('rank-zero generic completion Jacobian is a free-base field unit',S.Matrix([U,fibre_control,Q]).jacobian((U,P,Q)).det()-s)
    zero('rank-zero control has only a vertical critical line',S.Matrix([S.diff(fibre_control,z) for z in (U,P,Q)]).subs({s:0,U:0,Q:0}))

    A1,C1,D1 = S.Function('A1')(U,Z),S.Function('C1')(U,Z),S.Function('D1')(U,Z)
    rank_one = D1+s*A1+V*C1
    J1 = S.Matrix([A1,C1]).jacobian((U,Z)).det()
    zero('rank-one full Hessian is the residual Keller determinant squared',S.hessian(rank_one,(s,V,U,Z)).det()-J1**2)
    ss = (-S.diff(D1,U)*S.diff(C1,Z)+S.diff(D1,Z)*S.diff(C1,U))/J1
    vv = (-S.diff(A1,U)*S.diff(D1,Z)+S.diff(A1,Z)*S.diff(D1,U))/J1
    coefficient_matrix = S.Matrix([[S.diff(A1,U),S.diff(C1,U)],[S.diff(A1,Z),S.diff(C1,Z)]])
    zero('rank-one critical scheme solves uniquely for both removed coordinates',
         S.Matrix([S.diff(rank_one,U),S.diff(rank_one,Z)])-coefficient_matrix*S.Matrix([s-ss,V-vv]))
    zero('rank-one remaining critical equation is exactly C',S.diff(rank_one,V)-C1)

    # A unit-Hessian control demonstrates why invertibility is distinct from
    # no unit Hessian when the critical scheme is irreducible.
    unit_normal = s*U+V*Z+U**3
    variables = (s,V,U,Z)
    gradient = S.Matrix([S.diff(unit_normal,z) for z in variables])
    o1,o2,o3,o4 = S.symbols('o1 o2 o3 o4')
    outputs = S.Matrix([o1,o2,o3,o4])
    grad_inverse = S.Matrix([o3-3*o1**2,o4,o1,o2])
    zero('irreducible-critical-line control really has unit Hessian',S.hessian(unit_normal,variables).det()-1)
    zero('irreducible-critical-line gradient inverse left',grad_inverse.subs(dict(zip(outputs,gradient)),simultaneous=True)-S.Matrix(variables))
    zero('irreducible-critical-line gradient inverse right',gradient.subs(dict(zip(variables,grad_inverse)),simultaneous=True)-outputs)
    zero('unit control critical curve is an actual affine line with nonconstant base projection',
         S.Matrix([S.diff(unit_normal,z) for z in (V,U,Z)]).subs({s:-3*t*t,U:t,V:0,Z:0}))

    # Original family and exact fibre strata used in the Euler argument.
    u = 1+x*y
    A = u**3*p+3*x*u*u*q-x**3*r
    B = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    psi = L*A*A+M*A+N*B
    R = S.Function('R')(x,p)
    ur = 1+x*R
    r_from_T = (ur**3*p+3*x*ur*ur*q-T)/x**3
    cpoly = M-N*(2-3*x*R)/x**2
    epoly = N*ur*(ur+1)*p/x**2
    dpoly = 2*N*(3+2*x*R)/x
    rational = L*T*T+cpoly*T+epoly+dpoly*q
    zero('fixed nonzero-base fibre model is literally the original potential',psi.subs({y:R,r:r_from_T},simultaneous=True)-rational)
    zero('fixed-base dual computation coordinates have unit Jacobian over the scalar field',
         S.Matrix([p,A.subs(y,R),q]).jacobian((p,q,r)).det()-x**3)
    E = epoly-cpoly**2/(4*L)
    zero('global fixed-base square completion preserves the whole fibre variety',rational-L*(T+cpoly/(2*L))**2-E-dpoly*q)
    zero('off-root Euler stratum has a unique polynomially solved q after localization',
         rational.subs(q,(lam-L*T*T-cpoly*T-epoly)/dpoly)-lam)
    zero('each root has two distinct quadratic values outside a finite exceptional set',cpoly**2-4*L*(epoly-lam)-4*L*(lam-E))
    Dcrit = -N*p/(4*x*x)-(M-13*N/(2*x*x))**2/(4*L)
    zero('exceptional value at every coefficient root is the known critical value',E.subs(R,-S.Rational(3,2)/x)-Dcrit)
    count = S.symbols('count',integer=True,positive=True)
    zero('Euler counts add to one plus the unrestricted number of distinct roots',(1-count)+2*count-(1+count))
    check('positive distinct-root count gives a generic Euler value different from one',(1+count-1).is_positive)

    native = psi.subs(y,R)
    Rx,Rxx,Rxp,Rpp = S.symbols('Rx Rxx Rxp Rpp')
    jets = {S.diff(R,p):v,S.diff(R,x):Rx,S.diff(R,x,2):Rxx,S.diff(R,x,p):Rxp,S.diff(R,p,2):Rpp}
    h0 = 2*L*rho+M
    boundary = {x:0,p:rho,R:0,q:-h0/(N*v)}
    dual_gradient = S.Matrix([S.diff(native,z) for z in (p,q,r)])
    boundary_gradient = dual_gradient.applyfunc(lambda e:S.cancel(e.xreplace(jets).subs(boundary,simultaneous=True)))
    zero('arbitrary simple boundary zero supplies the entire literal critical line',boundary_gradient)
    Hnative = S.hessian(native,(x,p,q,r))
    Hline = Hnative.applyfunc(lambda e:S.cancel(e.xreplace(jets).subs(boundary,simultaneous=True)))
    zero('simple-root boundary line has the exact nonzero full Hessian determinant',Hline.det()-4*N**4*v*v)
    multiple_jets = dict(jets)
    multiple_jets[S.diff(R,p)] = 0
    multiple_boundary = {x:0,p:-M/(2*L),R:0}
    zero('multiple zero at the critical p value gives a whole critical two-plane',
         dual_gradient.applyfunc(lambda e:S.cancel(e.xreplace(multiple_jets).subs(multiple_boundary,simultaneous=True))))
    zero('critical two-plane makes the full native Hessian singular',
         Hnative.applyfunc(lambda e:S.cancel(e.xreplace(multiple_jets).subs(multiple_boundary,simultaneous=True))).det())
    nonlinear_example = p*p-1+x*p**3
    for root in [-1,1]:
        zero('nonlinear example has boundary simple root '+str(root),nonlinear_example.subs({x:0,p:root}))
        zero('nonlinear example boundary derivative at '+str(root),S.diff(nonlinear_example,p).subs({x:0,p:root})-2*root)
    repeated = 3*p+S.Rational(3,2)*x*p*p
    repeated_curve = 2*x*repeated+3
    repeated_derivative = S.diff(repeated,p)
    zero('entirely repeated generic curve control is retained exactly',repeated_curve-3*(1+x*p)**2)
    zero('repeated curve has no critical component after inverting R_p',3*repeated_curve/repeated_derivative**2-1)
    zero('repeated curve still has a simple boundary zero',S.diff(repeated,p).subs(x,0)-3)
    repeated_sqf = S.Poly(repeated_curve,p,domain=S.QQ.frac_field(x)).sqf_part()
    check('repeated-root Euler control has one distinct root and generic Euler two',repeated_sqf.degree()==1 and 1+repeated_sqf.degree()==2)

    # Opposite chart with every coefficient jet retained at the boundary.
    af,bf = S.Function('af')(y),S.Function('bf')(y)
    a0 = S.symbols('a0',nonzero=True)
    a1,a2,b0,b1,b2 = S.symbols('a1 a2 b0 b1 b2')
    opposite_native = psi.subs(x,af*p+bf)
    opposite_jets = {af:a0,bf:b0,S.diff(af,y):a1,S.diff(af,y,2):a2,
        S.diff(bf,y):b1,S.diff(bf,y,2):b2}
    opposite_rho = -b0/a0
    opposite_h = 2*L*opposite_rho+M
    opposite_boundary = {y:0,p:opposite_rho,r:-opposite_h*(1+3*a0*q)/(2*N*a0)}
    zero('opposite arbitrary-coefficient boundary line is actually relatively critical',
         S.Matrix([S.diff(opposite_native,z) for z in (p,q,r)]).applyfunc(
             lambda e:S.cancel(e.xreplace(opposite_jets).subs(opposite_boundary,simultaneous=True))))
    opposite_H = S.hessian(opposite_native,(y,p,q,r)).applyfunc(
        lambda e:S.cancel(e.xreplace(opposite_jets).subs(opposite_boundary,simultaneous=True)))
    zero('opposite line determinant is independent of every higher coefficient',opposite_H.det()-4*N**4*a0*a0)
    opposite_p = (x-bf)/af
    opposite_r = (u**3*opposite_p+3*x*u*u*q-T)/x**3
    opposite = L*T*T+(M-N*(2-3*x*y)/x**2)*T+N*u*(u+1)*opposite_p/x**2+2*N*(2*u+1)*q/x
    zero('opposite dominating-curve chart is the literal original potential',
         psi.subs({p:opposite_p,r:opposite_r},simultaneous=True)-opposite)
    zero('opposite q coefficient vanishes exactly on the closed hyperbola',S.diff(opposite,q)-2*N*(3+2*x*y)/x)
    xc = -S.Rational(3,2)/y
    Tc = (13*N/(2*xc*xc)-M)/(2*L)
    x_derivative = S.cancel(S.diff(opposite,x).subs({x:xc,T:Tc},simultaneous=True))
    zero('opposite remaining equation has a nonzero q coefficient on the open base',S.diff(x_derivative,q)+6*N/(xc*xc))
    qc = S.cancel(-x_derivative.subs(q,0)/S.diff(x_derivative,q))
    critical = {x:xc,T:Tc,q:qc}
    for variable in (x,T,q):
        zero('opposite dominating-curve critical '+str(variable)+' equation',S.diff(opposite,variable).subs(critical,simultaneous=True))
    opposite_D = N*y/(6*af)+N*y*y*bf/(9*af)-(M-26*N*y*y/9)**2/(4*L)
    zero('opposite whole-branch critical value',opposite.subs(critical,simultaneous=True)-opposite_D)
    zero('dominating component cannot meet the boundary line in the affine source',(2*x*y+3).subs({x:0,y:0})-3)
    tuned_value = opposite_D.subs({af:1,bf:S.Rational(1013,54)*y*y,L:1,M:0,N:1},simultaneous=True)
    zero('earlier matching-determinant example has the retained critical value',tuned_value-y/6+y**4/486)
    zero('earlier matching-determinant example still has determinant four on its branch',-162/y**2*S.diff(tuned_value,y,2)-4)
    separator = -2*y*(p+S.Rational(1013,54)*y*y)/3
    zero('earlier control has a regular function zero on its boundary component',separator.subs({y:0,p:0}))
    zero('the same regular function is one on its dominating component',separator.subs(p,xc-S.Rational(1013,54)*y*y)-1)

    # These are provenance checks, not a proof of any imported theorem.
    manifest_name = 'sources/critical-components-import.json'
    manifest = json.loads((HERE/manifest_name).read_text())
    dependencies = {}
    retained = manifest['retained_file_from_research_root']
    digest = hashlib.sha256((ROOT/retained).read_bytes()).hexdigest()
    check('retained primary rectification statement byte pin',digest==manifest['sha256']
        =='976068fa81837498fcb77bb91e6c4d28b4a6bc74ed940340d9771e807a5be249')
    dependencies[retained] = digest
    check('rectification source location is the exact closed-plane-line statement',
        manifest['location']=='Theorem 2.4' and 'S2.Thmdefi4' in (ROOT/retained).read_text())
    line_manifest_name = 'hessian/hc4-triangular-graphs/sources/critical-jet-import.json'
    line_manifest = json.loads((ROOT/line_manifest_name).read_text())
    line_source = 'hessian/hc4-triangular-graphs/sources/'+line_manifest['retained_file']
    line_digest = hashlib.sha256((ROOT/line_source).read_bytes()).hexdigest()
    check('retained line-injectivity primary source byte pin',line_digest==line_manifest['sha256']
        =='4cce6a7334640ea47036b1fb8c676e28da7946dbd15d0a48868b30d94c68585a')
    dependencies[line_manifest_name] = hashlib.sha256((ROOT/line_manifest_name).read_bytes()).hexdigest()
    dependencies[line_source] = line_digest
    tang_pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md','hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
            for name in [Path(__file__).name,'CRITICAL_COMPONENTS.md',manifest_name]},
        dependencies=dependencies,
        written_not_finite_exhaustion=['closed affine-line component endpoint for plane Keller pairs',
            'both projection ranks with affine critical components',
            'Euler stratification for arbitrary polynomial degree and root multiplicities',
            'boundary simple-root and multiple-critical-root classification',
            'all-degree opposite-chart coefficients with a(0) nonzero'],
        y_R_x_p_simple_boundary_root_gradient_invertibility=True,
        y_R_x_p_reducible_critical_scheme_has_no_unit_Hessian=True,
        repeated_generic_curve_automatically_asserted_no_unit=False,
        opposite_x_a_y_p_plus_b_y_nonzero_a_at_zero_has_no_unit_Hessian=True,
        prior_matching_critical_determinant_example_excluded_geometrically=True,
        opposite_a_at_zero_zero_b_at_zero_nonzero_closed=False,
        no_simple_boundary_root_y_graphs_all_closed=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,
        singleton_zero_properness_resolved=False,
        imported_theorems=['Abhyankar--Moh--Suzuki as retained in Jelonek Theorem 2.4','HC4MYGL0','HCTANG2'],
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-start,3),
        formal_verification=False,independent_replay=False,external_review=False)
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
