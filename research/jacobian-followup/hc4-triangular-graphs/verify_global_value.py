"""Exact controls for whole-value profiles and full source orbits.

Written degree-independent arguments and retained external results are
not certified by this finite list of identities.
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
    started = time.monotonic()
    checks = []

    def check(name,value):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    def zero(name,value):
        values = list(value) if isinstance(value,S.MatrixBase) else [value]
        check(name,all(S.cancel(item.doit())==0 for item in values))

    s,U,P,Q = S.symbols('s U P Q')
    variables = (s,U,P,Q)
    a,b,c,d = S.symbols('a b c d')
    aa,bb,cc = S.symbols('aa bb cc')
    jac = S.Matrix([[a,b],[c,d]])
    top = S.Matrix([[aa,bb],[bb,cc]])
    cotangent_hessian = S.BlockMatrix([[top,jac.T],[jac,S.zeros(2)]]).as_explicit()
    zero('global cotangent Hessian determinant is the plane Jacobian square',
         cotangent_hessian.det()-jac.det()**2)
    Du,Dv = S.symbols('Du Dv')
    solved = -jac.T.inv()*S.Matrix([Du,Dv])
    zero('full critical equations solve both dual variables exactly',
         jac.T*solved+S.Matrix([Du,Dv]))

    F,G,D,value = S.symbols('F G D value')
    raw = D+P*F+Q*G
    zero('first global value stratum has an actual affine-line inverse',
         raw.subs(P,(value-D-Q*G)/F)-value)
    zero('second disjoint stratum has its actual affine-line inverse',
         raw.subs(F,0).subs(Q,(value-D)/G)-value)
    zero('the zero-row exceptional value leaves both dual coordinates free',
         raw.subs({F:0,G:0,D:value},simultaneous=True)-value)
    hu,hv,n,nv,degree,nab = S.symbols('hu hv n nv degree nab')
    t = hu*hv
    global_E = (t*t-n)*t+nv*t*t
    zero('entire global value profile retains every exceptional affine plane',
         global_E-(t**3-n*t+nv*t*t))
    zero('generic whole value profile has coefficient minus the critical count',
         global_E.subs(nv,0)-(t**3-n*t))
    zero('special value jump counts all critical points at that value',
         global_E-global_E.subs(nv,0)-nv*t*t)
    zero('linear-perturbation coefficient difference is exactly the missing length',
         ((t**3-nab*t)-(t**3-degree*t))-(degree-nab)*t)
    zero('the coordinate orbit has the empty-critical profile at every value',
         global_E.subs({n:0,nv:0})-t**3)
    zero('the quadratic orbit has the singleton generic profile',
         global_E.subs({n:1,nv:0})-(t**3-t))
    zero('the quadratic critical value has its entire extra plane',
         global_E.subs({n:1,nv:1})-(t**3+t*t-t))

    # A full four-variable word mixes both quotient and dual coordinates.
    phi = S.Matrix([s+P**2,U+Q**2,P+(U+Q**2)**2,Q+(s+P**2)**2])
    z0,z1,z2,z3 = S.symbols('z0:4')
    backP,backQ = z2-z1**2,z3-z0**2
    back = {s:z0-backP**2,U:z1-backQ**2,P:backP,Q:backQ}
    zero('the fully mixed source word has an actual polynomial inverse',
         phi.subs(back,simultaneous=True)-S.Matrix([z0,z1,z2,z3]))
    zero('the fully mixed source word has unit Jacobian',phi.jacobian(variables).det()-1)
    native = s*P+U*Q
    pulled = phi[0]*phi[2]+phi[1]*phi[3]
    native_gradient_at_phi = S.Matrix([phi[2],phi[3],phi[0],phi[1]])
    zero('the full gradient ideal transforms by an invertible Jacobian matrix',
         S.Matrix([S.diff(pulled,z) for z in variables])-phi.jacobian(variables).T*native_gradient_at_phi)
    extra = S.zeros(4)
    for ii in range(4):
        extra += native_gradient_at_phi[ii]*S.hessian(phi[ii],variables)
    zero('the full Hessian chain rule retains the gradient-weighted correction',
         S.hessian(pulled,variables)-phi.jacobian(variables).T*S.hessian(native,variables)*phi.jacobian(variables)-extra)
    origin = dict.fromkeys(variables,0)
    zero('the critical source point retains its Hessian determinant',
         S.hessian(pulled,variables).subs(origin).det()-S.hessian(native,variables).det())
    off = {s:1,U:0,P:0,Q:0}
    check('off the full critical scheme the Hessian determinant need not be preserved',
          S.hessian(pulled,variables).subs(off).det()!=S.hessian(native,variables).det())

    # Exact native normal jets, before imposing the critical value of p.
    x,y,p,q,r,L,M,N,T0 = S.symbols('x y p q r L M N T0')
    u = 1+x*y
    A = u**3*p+3*x*u**2*q-x**3*r
    B = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    potential = L*A*A+M*A+N*B
    base_origin = {x:0,y:0}
    Dx = S.diff(potential,x).subs(base_origin,simultaneous=True)
    Dy = S.diff(potential,y).subs(base_origin,simultaneous=True)
    zero('native x normal jet at the base origin',Dx-(3*(2*L*p+M)*q+2*N*r))
    zero('native y normal jet at the base origin',Dy-N*q)
    ps = -M/(2*L)
    zero('arbitrary zero-boundary y graph has the exact remaining full derivative',
         (Dx+T0*Dy).subs(p,ps)-N*(q*T0+2*r))
    zero('arbitrary zero-boundary opposite graph has its exact remaining derivative',
         (T0*Dx+Dy).subs(p,ps)-N*(q+2*r*T0))
    Tfun = S.Function('Tfun')(q,r)
    zero('the y-graph critical curve has a nonzero transverse derivative at its fixed point',
         S.diff(N*(q*Tfun+2*r),r).subs({q:0,r:0},simultaneous=True)-2*N)
    zero('the opposite critical curve has a nonzero transverse derivative at its fixed point',
         S.diff(N*(q+2*r*Tfun),q).subs({q:0,r:0},simultaneous=True)-N)
    normal2,normalp = S.symbols('normal2 normalp')
    for qcross,rcross,label in [(N*T0,2*N,'y'),(N,2*N*T0,'opposite')]:
        critical_hessian = S.Matrix([[normal2,normalp,qcross,rcross],
            [normalp,2*L,0,0],[qcross,0,0,0],[rcross,0,0,0]])
        zero('all off-boundary jets retain a singular full critical Hessian '+label,
             critical_hessian.det())

    bfun = S.Function('bfun')(Q)
    split = s*s+U*U+P*P+bfun
    zero('split one-variable critical Hessian has the exact second-derivative factor',
         S.hessian(split,variables).det()-8*S.diff(bfun,Q,2))
    cubic = Q**3-3*Q
    check('a cubic with two Morse points has incompatible critical Hessian values',
          [8*S.diff(cubic,Q,2).subs(Q,rho) for rho in (-1,1)]==[-48,48])
    check('a multiple derivative root gives a degenerate full critical point',
          S.diff(Q**4,Q).subs(Q,0)==0 and S.diff(Q**4,Q,2).subs(Q,0)==0)
    k = S.Symbol('k',integer=True,positive=True)
    check('the derivative comparison uses a strict degree gap in every degree',S.simplify((k-1)-(k-2))==1)

    # An honest coordinate with nonunit Hessian is not promoted to an endpoint solution.
    coord = s*P+(s*U-1)*Q
    Wcoord = P+U*Q
    Vtarget,Wtarget = S.symbols('Vtarget Wtarget')
    pb = U*Vtarget-(s*U-1)*Wtarget
    qb = -Vtarget+s*Wtarget
    zero('a unit-ideal row gives a genuine coordinate and polynomial companion inverse',
         S.Matrix([coord,Wcoord]).subs({P:pb,Q:qb},simultaneous=True)-S.Matrix([Vtarget,Wtarget]))
    zero('the coordinate control really has nonunit full Hessian',
         S.hessian(coord,variables).det()-s*s)

    zrow = S.Matrix([[s,U]])
    generic_matrix = S.Matrix([[a,b],[c,d]])
    grow = zrow*generic_matrix
    mu,nu = S.symbols('mu nu')
    corrected = generic_matrix+S.Matrix([-U,s])*S.Matrix([[nu,-mu]])
    zero('the row correction preserves the entire output row',zrow*corrected-grow)
    zero('the row correction removes any prescribed ideal determinant error',
         corrected.det()-generic_matrix.det()+mu*grow[0]+nu*grow[1])
    singleton = s*P+(U+s*s)*Q
    zero('the singleton control has an actual relative source change to the quadratic',
         native.subs(P,P+s*Q)-singleton)
    zero('the singleton control has unit full Hessian',S.hessian(singleton,variables).det()-1)
    Dterm = s**3+U**2
    zero('every maximal-ideal offset can be absorbed, explicit control',
         Dterm-(s*(s*s-s*U)+(U+s*s)*U))
    zero('the fixed p-graph coefficient base change is nondegenerate',
         S.Matrix([M*x+N*y,2*N*x]).jacobian((x,y)).det()+2*N*N)

    dependencies = {}
    for proof,receipt in [('INTRINSIC_AFFINE_PLANE.md','receipt-intrinsic-plane.json'),
                          ('HODGE_TRACE_CLOSURE.md','receipt-hodge-trace.json'),
                          ('RELATIVE_COORDINATE_ENDPOINT.md','receipt-relative-coordinate.json'),
                          ('P_GRAPH_SINGLETON_ENDPOINT.md','receipt-square-axis.json')]:
        retained = json.loads((HERE/receipt).read_text())
        digest = hashlib.sha256((HERE/proof).read_bytes()).hexdigest()
        check('retained proof pin '+proof,digest==retained['source_sha256'][proof])
        dependencies[str((HERE/proof).relative_to(ROOT))] = digest
    row_proof = ROOT/'jacobian-followup/hc4-relative-row-orbit/PROOF.md'
    # The retained p-graph receipt pins the row proof; the original row
    # receipt instead pins its primary sources by URL.
    row_receipt = json.loads((row_proof.parent/'receipt.json').read_text())
    row_digest = hashlib.sha256(row_proof.read_bytes()).hexdigest()
    p_pins = json.loads((HERE/'receipt-square-axis.json').read_text())['dependencies']
    check('retained row-orbit and finite-normalization proof pin',row_digest==p_pins[str(row_proof.relative_to(ROOT))])
    dependencies[str(row_proof.relative_to(ROOT))] = row_digest
    for filename,url in [('stacks-zmt-02LR.html','https://stacks.math.columbia.edu/tag/02LR'),
                         ('stacks-flatness-00R3.html','https://stacks.math.columbia.edu/tag/00R3')]:
        source = row_proof.parent/'sources'/filename
        digest = hashlib.sha256(source.read_bytes()).hexdigest()
        check('retained primary normalization input '+filename,digest==row_receipt['primary_source_sha256'][url])
        dependencies[str(source.relative_to(ROOT))] = digest
    tang = json.loads((ROOT/'jacobian-followup/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['jacobian-followup/tang-hc4-audit/STRUCTURAL_REVIEW.md',
                 'jacobian-followup/tang-hc4-audit/source-v1.pdf']:
        digest = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,digest==tang[name])
        dependencies[name] = digest
    manifest = json.loads((HERE/'sources/hodge-trace-import.json').read_text())
    source = HERE/'sources'/manifest['retained_file']
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    check('retained primary Hodge source',digest==manifest['sha256'])
    dependencies[str(source.relative_to(ROOT))] = digest
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in [Path(__file__).name,'GLOBAL_VALUE_ORBITS.md']},
        dependencies=dependencies,
        written_not_finite_exhaustion=['whole value fibres and full critical scheme under all source automorphisms',
            'nonautomorphic cotangent whole value profile at every scalar value',
            'linear-perturbation coefficient defect equals omitted finite-normalization length',
            'all zero-boundary base graphs excluded under arbitrary full source automorphisms',
            'all split one-variable potentials of degree at least three excluded in the full source orbit',
            'all origin-avoiding base-only graphs retain exactly the unit-ideal coordinate endpoint',
            'the fixed p graph retains exactly the singleton-zero quadratic endpoint under full source mixing'],
        source_automorphisms_required_to_preserve_free_base=False,
        zero_boundary_graph_exclusion_requires_Tang=False,
        split_high_degree_exclusion_requires_Tang=False,
        full_Hessian_asserted_invariant_off_critical_scheme=False,
        global_profile_asserted_for_every_automorphic_gradient=False,
        Hodge_profile_asserted_to_prove_topological_local_triviality=False,
        whole_linear_perturbation_defect_asserted_zero=False,
        source_orbit_unit_Hessian_assumed_automatic=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
        unrestricted_HC4_or_JC2_resolved=False,
        formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
