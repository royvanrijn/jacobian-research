"""Exact controls for the written relative critical-graph classification."""
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
        check(name,all(S.cancel(S.expand(item.doit()))==0 for item in values))

    s,t,u,v,U,V = S.symbols('s t u v U V')
    source = (s,t,u,v)
    # Rank-one projection: two different vertical-Hessian alternatives.
    A,C = t,U+t*t
    D = t*U+U**3
    f = D+s*A+V*C
    zero('rank-one projection has the claimed vertical Hessian determinant',
        S.hessian(f,(U,V)).det()+S.diff(C,U)**2)
    zero('the rank-one projection control has actual unit full Hessian',
        S.hessian(f,(s,t,U,V)).det()-1)
    zero('nonzero C_U gives the unique V solution on C=0',
        S.Matrix([S.diff(f,U),S.diff(f,V)]).subs({U:-t*t,V:-t-3*t**4},simultaneous=True))
    A,C = U+t*t,t
    f = t*U+s*A+V*C
    zero('C independent of U gives the triangular Keller branch',
        S.Matrix([A,C]).jacobian((t,U)).det()+1)
    zero('its U derivative vanishes on the relative critical scheme',S.diff(C,U))
    zero('a whole V line appears at the special base point',
        S.Matrix([S.diff(f,U),S.diff(f,V)]).subs({t:0,s:0},simultaneous=True))

    # A rank-two graph with an actual polynomial section below three affines.
    H = u+v*v
    G = v+H*H
    K = S.Matrix([H,G])
    JK = K.jacobian((u,v))
    D = -(H*H+G*G)/2
    f = D+s*H+t*G
    P = -JK.T.inv()*S.Matrix([S.diff(D,u),S.diff(D,v)])
    zero('the two nonlinear free coordinates have unit Jacobian',JK.det()-1)
    zero('the critical graph projection is exactly the free coordinate pair',P-K)
    graph = {s:P[0],t:P[1]}
    HE = S.hessian(f,(u,v)).subs(graph,simultaneous=True)
    zero('the relative critical graph satisfies both directional equations',
        S.Matrix([S.diff(f,u),S.diff(f,v)]).subs(graph,simultaneous=True))
    zero('the vertical Hessian and graph differential obey the full identity',HE+JK.T*P.jacobian((u,v)))
    envelope_on_graph = f.subs(graph,simultaneous=True)
    zero('the graph envelope has the claimed differential',
        S.Matrix([S.diff(envelope_on_graph,u),S.diff(envelope_on_graph,v)])-
        P.jacobian((u,v)).T*K)
    section = {v:t-s*s,u:s-(t-s*s)**2}
    zero('the section has a literal polynomial inverse over the retained base',
        K.subs(section,simultaneous=True)-S.Matrix([s,t]))
    zero('the section control has unit actual full Hessian',S.hessian(f,source).det()-1)
    zero('the envelope in retained coordinates is a binary nondegenerate quadratic',
        f.subs(section,simultaneous=True)-(s*s+t*t)/2)
    check('the two remaining free pullbacks have distinct degrees two and four',
        [S.Poly(H,u,v).total_degree(),S.Poly(G,u,v).total_degree()]==[2,4])
    Z,X,Y = S.symbols('Z X Y')
    Q = s*X+t*Y-(X*X+Y*Y)/2+Z*Z
    check('the strict two-affine example has nondegenerate full ambient quadratic',
        S.hessian(Q,(s,t,Z,X,Y)).rank()==5)
    zero('its ambient pullback has the exact source potential',Q.subs({Z:0,X:H,Y:G},simultaneous=True)-f)

    # The rank-one relative critical product and its curvature identity.
    fline = s*v+t*t+u*u
    zero('the relative curve-product control has nonzero constant full Hessian',
        S.hessian(fline,source).det()+4)
    zero('its relative critical locus lies over a straight base curve',
        S.Matrix([S.diff(fline,u),S.diff(fline,v)])-S.Matrix([2*u,s]))
    check('its vertical Hessian has rank one on the whole critical product',
        S.hessian(fline,(u,v)).rank()==1)
    rho,sigma,h = S.symbols('rho sigma h',nonzero=True)
    aa,bb,B1,B2 = [S.Function(name)(rho) for name in ['aa','bb','B1','B2']]
    pair = S.Matrix([h*sigma*S.diff(bb,rho)+B1,-h*sigma*S.diff(aa,rho)+B2])
    curvature = S.diff(aa,rho,2)*S.diff(bb,rho)-S.diff(aa,rho)*S.diff(bb,rho,2)
    remainder = S.det(S.Matrix([[S.diff(B1,rho),S.diff(bb,rho)],
        [S.diff(B2,rho),-S.diff(aa,rho)]]))
    zero('the generic Keller determinant has exactly the curvature coefficient',
        pair.jacobian((rho,sigma)).det()-h*h*sigma*curvature-h*remainder)
    zero('the tangential output row is independent of the free parameter',
        S.diff(S.diff(aa,rho)*pair[0]+S.diff(bb,rho)*pair[1],sigma))
    curved_pair = S.Matrix([2*h*rho*sigma+rho,-h*sigma])
    zero('a nonlinear parabolic base image leaves an uncancellable sigma coefficient',
        S.diff(curved_pair.jacobian((rho,sigma)).det(),sigma)+2*h*h)
    straight_pair = S.Matrix([rho,sigma+rho**3])
    zero('the straight-image endpoint is an actual triangular automorphism',
        straight_pair.jacobian((rho,sigma)).det()-1)
    tau,uu,vv = S.symbols('tau uu vv')
    comps,genus,punctures = S.symbols('comps genus punctures')
    curve_E = comps*uu*vv-genus*(uu+vv)+comps-punctures
    solutions = S.solve(S.Poly(curve_E*uu*vv-(uu*vv)**2,uu,vv).coeffs(),
        (comps,genus,punctures),dict=True)
    check('all smooth curve coefficients uniquely force one genus-zero component and one puncture',
        solutions==[{comps:1,genus:0,punctures:1}])

    # Relative coordinate chain rule and why it only preserves vertical
    # Hessians on the relative critical scheme.
    X0 = u+v*v
    Y0 = v+X0*X0
    psi = S.Matrix([X0,Y0])
    Jpsi = psi.jacobian((u,v))
    native = X*X+s*Y+t*t
    pulled = native.subs({X:X0,Y:Y0},simultaneous=True)
    correction = S.hessian(pulled,(u,v))-Jpsi.T*S.hessian(native,(X,Y))*Jpsi
    zero('the vertical Hessian correction is precisely the relative-gradient chain term',
        correction-(2*X0*S.hessian(X0,(u,v))+s*S.hessian(Y0,(u,v))))
    zero('that correction vanishes on the entire relative critical product',correction.subs({u:-v*v,s:0},simultaneous=True))
    check('the correction is nonzero away from relative criticality',correction!=S.zeros(2,2))
    zero('a rank-zero relative potential has zero vertical Hessian',
        S.hessian(s*X+t*Y,(X,Y)))
    zero('a rank-one B=0 critical fibre has one free coordinate',
        S.Matrix([S.diff(X*X+s*Y,X),S.diff(X*X+s*Y,Y)])-S.Matrix([2*X,s]))
    out = S.symbols('out')
    zero('a constant nonzero Y coefficient gives a literal relative coordinate',
        (X*X+Y+s*X+t).subs(Y,out-X*X-s*X-t)-out)

    # Orientation at coefficient zeros is stronger than a squared determinant.
    p = t*t-1
    zero('the two-root boundary has opposite unsquared derivatives',S.diff(p,t).subs(t,1)+S.diff(p,t).subs(t,-1))
    check('those unsquared derivatives differ',S.diff(p,t).subs(t,1)!=S.diff(p,t).subs(t,-1))
    zero('their squared critical determinants nevertheless agree',
        S.diff(p,t).subs(t,1)**2-S.diff(p,t).subs(t,-1)**2)
    coefficient_pair = S.Matrix([s,t+s*t*t])
    transformed = S.Matrix([coefficient_pair[0]+t*coefficient_pair[1],coefficient_pair[1]])
    zero('an oriented free frame preserves the Jacobian value at the coefficient zero',
        transformed.jacobian((s,t)).det().subs({s:0,t:0})-
        coefficient_pair.jacobian((s,t)).det().subs({s:0,t:0}))
    # Old endpoint constructions are reused, with explicit relative inverses.
    F,G0 = s,s*t-1
    a,b = t,-1
    value = u*F+v*G0
    companion = -b*u+a*v
    zero('the unit-ideal endpoint control has a base Bezout certificate',a*F+b*G0-1)
    zero('its free coordinate completion has a literal global inverse',
        S.Matrix([a*value-G0*companion,b*value+F*companion])-S.Matrix([u,v]))
    zero('the endpoint graph is a genuine relative-coordinate graph',
        S.Matrix([value-s*t-companion**2,companion]).jacobian((u,v)).det()-1)
    zero('the displayed unit-ideal control is explicitly non-Keller',
        S.Matrix([F,G0]).jacobian((s,t)).det()-s)
    frame = S.Matrix([[1,s],[t,1+s*t]])
    zero('the singleton endpoint frame has constant determinant',frame.det()-1)
    free = frame*S.Matrix([u,v])
    row = S.Matrix([[s,t]])*frame
    zero('the singleton graph has the exact cotangent potential',
        s*free[0]+t*free[1]-(row*S.Matrix([u,v]))[0])

    pins = {}
    for receipt_name,name in [
        ('receipt-global-value.json','GLOBAL_VALUE_ORBITS.md'),
        ('receipt-projection-profile.json','MIXED_PROJECTION_PROFILE.md'),
        ('receipt-two-affine-geometry.json','TWO_AFFINE_ENDPOINT_GEOMETRY.md')]:
        retained = json.loads((HERE/receipt_name).read_text())
        pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
    retained = json.loads((HERE/'receipt-intrinsic-plane.json').read_text())
    for name,digest in retained['dependencies'].items():
        if 'tang-hc4-audit' in name or 'ricolfi-hodge' in name:
            pins[name] = digest
    record = json.loads((HERE/'sources/relative-critical-graph-import.json').read_text())
    for item in record['sources']:
        pins[str((HERE/'sources'/item['retained_file']).relative_to(ROOT))] = item['sha256']
    for name,digest in pins.items():
        check('retained mathematical or primary source pin '+name,
            hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
            ['verify_relative_critical_graph.py','RELATIVE_CRITICAL_GRAPH_CLASSIFICATION.md','sources/relative-critical-graph-import.json']},
        dependencies=pins,external_theorem_inputs=['HC2, n=2 assertion of de Bondt Theorem 2.1'],
        written_not_finite_exhaustion=[
            'both relative critical-scheme criteria in every projection rank of Tang plane',
            'curve-product straightness from the Keller syzygy and curvature identity',
            'all constant relative quadratic ranks and all intrinsic-plane positions in the base-only graph family',
            'at most one full critical point and exact disjunction of the two old endpoints',
            'no relative unitization of native graph members with multiple full critical points'],
        cotangent_alignment_assumed=False,
        generic_critical_section_asserted_global=False,
        abstract_product_asserted_base_compatible=False,
        nonlinear_source_Hessian_invariance_used=False,
        arbitrary_four_variable_source_mixing_asserted_covered=False,
        all_two_affine_embeddings_asserted_base_only_graphs=False,
        non_Keller_unit_ideal_control_asserted_Keller=False,
        endpoint_pair_constructed=False,unit_ideal_endpoint_resolved=False,
        singleton_zero_endpoint_resolved=False,unrestricted_HC4_or_JC2_resolved=False,
        new_external_normal_form_classification_imported=False,
        formal_verification=False,independent_replay=False,external_review=False,novelty_priority_asserted=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-relative-critical-graph.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=len(checks),elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
