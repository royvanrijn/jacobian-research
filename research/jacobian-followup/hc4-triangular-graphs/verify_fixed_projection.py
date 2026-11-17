"""Exact controls for fixed-potential mixed projection obstructions.

Written proofs establish the all-degree conclusions; these controls
separate necessary boundary data from a complete graph realization.
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

    s,U,w1,w2,e1,e2,L,M,N = S.symbols('s U w1 w2 e1 e2 L M N')
    Ps,Pu,Qs,Qu,Pv,Qv,c = S.symbols('Ps Pu Qs Qu Pv Qv c')
    jac_formal = S.Matrix([[Ps,Pu],[Qs,Qu]])
    adj_times_K = jac_formal.adjugate()*S.Matrix([Pv,Qv])
    zero('the second pulled-back Euler component has the exact adjugate formula',
         adj_times_K[1]-(Ps*Qv-Qs*Pv))
    zero('its zero-P residue is precisely the derivative factor times the Keller mate',
         adj_times_K[1].subs(Pv,0)-Ps*Qv)

    K = S.Matrix([s,U+s*s])
    jac = K.jacobian((s,U))
    b = jac.inv()*K
    D = s**3+U*U
    gradD = S.Matrix([S.diff(D,s),S.diff(D,U)])
    f = D+w1*K[0]+w2*K[1]
    eta = S.Matrix([S.diff(f,s),S.diff(f,U)])
    zero('the pulled-back Euler vector is polynomial and pushes to the radial target',jac*b-K)
    zero('cotangent reconstruction gives an exact polynomial identity for f',
         D-(b.T*gradD)[0]+(b.T*eta)[0]-f)
    zero('eta has an actual polynomial source inverse',
         jac.T.inv()*(eta-gradD)-S.Matrix([w1,w2]))
    zero('the control has actual unit full Hessian',S.hessian(f,(s,U,w1,w2)).det()-1)
    zero('the pulled-back Euler derivative is the identity on the zero fibre',
         b.jacobian((s,U)).subs({s:0,U:0},simultaneous=True)-S.eye(2))

    # Positive rank-one boundary control; no full Meng--Yang embedding inferred.
    triangular_K = S.Matrix([U,s])
    triangular_D = -L*s*s
    triangular_f = triangular_D+w1*U+w2*s
    triangular_eta1 = S.diff(triangular_f,s)
    zero('a genuine triangular Keller pair has the required pure zero component',
         triangular_K.jacobian((s,U)).det()+1)
    zero('the rank-one common zero scheme leaves exactly a quadratic composite',
         triangular_f.subs({U:0,w2:2*L*s+M},simultaneous=True)-(L*s*s+M*s))
    zero('the arbitrary-offset base equation is solved in the same polynomial fibre',
         (triangular_eta1-M).subs(w2,2*L*s+M))
    zero('the positive boundary control retains unit full Hessian',
         S.hessian(triangular_f,(s,U,w1,w2)).det()-1)
    a,bb,g,t1,t2 = S.symbols('a bb g t1 t2',nonzero=True)
    pp,qq = a*U+bb,g*s+U**4+U
    Ub = (t1-bb)/a
    sb = (t2-Ub**4-Ub)/g
    zero('the forced triangular pair has an explicit source inverse',
         S.Matrix([sb,Ub]).subs({t1:pp,t2:qq},simultaneous=True)-S.Matrix([s,U]))
    zero('the other triangular inverse composition is exact',
         S.Matrix([pp,qq]).subs({s:sb,U:Ub},simultaneous=True)-S.Matrix([t1,t2]))
    zero('its unit Jacobian has the asserted sign',S.Matrix([pp,qq]).jacobian((s,U)).det()+a*g)
    h0,h1,Htop,T = S.symbols('h0 h1 Htop T')
    for degree in (2,3,5):
        composite = Htop*(h0+h1*e2)**degree+(h0+h1*e2)
        check('a nonconstant affine free-variable input has composite degree '+str(degree),
              S.degree(composite,e2)==degree)
        zero('its leading coefficient cannot disappear in a component domain, degree '+str(degree),
             S.Poly(composite,e2).nth(degree)-Htop*h1**degree)
    formal_poly = s**4+U*s*s+U**3+1
    check('a nonzero source derivative has strictly smaller degree than its polynomial',
          S.Poly(S.diff(formal_poly,s),s,U).total_degree()<S.Poly(formal_poly,s,U).total_degree())
    check('a repeated Keller-component candidate violates the unit Jacobian hypothesis',
          S.Matrix([s*s,U]).jacobian((s,U)).det()==2*s)
    # An offset pure fibre changes the exact divisibility statement.
    alpha = S.symbols('alpha')
    zero('a translated pure fibre has the additional derivative term and is not silently covered',
         adj_times_K[1].subs(Pv,alpha)-(Ps*Qv-alpha*Qs))

    # The native common-zero-base boundary and its complete first jet.
    x,y,p,q,r = S.symbols('x y p q r')
    u = 1+x*y
    A = u**3*p+3*x*u*u*q-x**3*r
    B = y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*r
    Psi = L*A*A+M*A+N*B
    origin_base = {x:0,y:0}
    zero('the exact native boundary is quadratic in p alone',
         Psi.subs(origin_base,simultaneous=True)-(L*p*p+M*p))
    native_grad = S.Matrix([S.diff(Psi,z) for z in (x,y,p,q,r)])
    zero('the entire native first jet on the common-zero-base scheme is exact',
         native_grad.subs(origin_base,simultaneous=True)-
         S.Matrix([3*(2*L*p+M)*q+2*N*r,N*q,2*L*p+M,0,0]))
    pstar,lstar = -M/(2*L),-M*M/(4*L)
    zero('the distinguished critical value is the completed-square value',
         L*p*p+M*p-lstar-L*(p-pstar)**2)
    zero('the distinguished ambient point is actually critical',
         native_grad.subs({x:0,y:0,p:pstar,q:0,r:0},simultaneous=True))
    zero('its value is exactly lambda star',
         Psi.subs({x:0,y:0,p:pstar,q:0,r:0},simultaneous=True)-lstar)
    zero('at the critical boundary p the only remaining differential is in the two bases',
         native_grad.subs({x:0,y:0,p:pstar},simultaneous=True)-S.Matrix([2*N*r,N*q,0,0,0]))
    z0,z1,z2,z3 = S.symbols('z0:4')
    gz = S.Matrix([z0,z1,z2,z3])
    phi = S.Matrix([z0+z1*z1,z2+z3**3])
    check('a nonlinear gradient-centered base control vanishes at zero with rank two',
          phi.subs(dict.fromkeys(gz,0))==S.zeros(2,1) and
          phi.jacobian(gz).subs(dict.fromkeys(gz,0)).rank()==2)
    check('a translated gradient base need not be centered and is not silently included',
          S.Matrix([z0+1,z2]).subs(dict.fromkeys(gz,0))!=S.zeros(2,1))

    # The retained fixed p-graph has this exact unit-Hessian normal form.
    TT,WW = S.symbols('TT WW')
    graph_f = (M*x+N*y)*TT+2*N*x*WW
    graph_grad = S.Matrix([S.diff(graph_f,z) for z in (x,y,TT,WW)])
    zero('the retained fixed graph normal form has its exact constant full Hessian',
         S.hessian(graph_f,(x,y,TT,WW)).det()-4*N**4)
    zero('its physical base is a homogeneous rank-two gradient projection',
         S.Matrix([graph_grad[3]/(2*N),(graph_grad[2]-M*graph_grad[3]/(2*N))/N])-S.Matrix([x,y]))
    zero('the retained normal form has exactly the zero critical fibre',
         S.Matrix([graph_grad[3]/(2*N),(graph_grad[2]-M*graph_grad[3]/(2*N))/N,
                   graph_grad[1]/N,(graph_grad[0]-M*graph_grad[1]/N)/(2*N)])-S.Matrix([x,y,TT,WW]))
    zero('the fixed physical p-graph has p zero on its whole base-zero fibre',
         (-2*x*((x*y)**2+1)*q-2*x**3*r).subs(origin_base,simultaneous=True))
    check('the ambient critical value equals its single source critical value exactly when M is zero',
          lstar.subs(M,0)==0 and S.cancel(lstar)!=0)

    # Necessary square equation on the remaining full-eta-rank stratum.
    Cmat = S.Matrix([[1,2],[3,7]])
    Amat = S.Matrix([[2,-1],[1,4]])
    tmat = Cmat.inv()*Amat
    V = D-(b.T*gradD)[0]-(b.T*tmat*K)[0]
    formal_f = D-(b.T*gradD)[0]+(b.T*S.Matrix([e1,e2]))[0]
    zero('solving the full-rank physical zero base yields the exact global square-equation left side',
         formal_f.subs({e1:-(tmat*K)[0],e2:-(tmat*K)[1]},simultaneous=True)-V)
    unmixed_V = D-(b.T*gradD)[0]
    check('the negative square-equation control with zero K mixing has odd total degree',
          S.Poly(unmixed_V,s,U).total_degree()%2==1)
    zero('its gradient vanishes at the actual zero fibre of K',
         S.Matrix([S.diff(V,z) for z in (s,U)]).subs({s:0,U:0},simultaneous=True))
    zero('the square-equation value at a Keller zero is the original D value',
         (V-D).subs({s:0,U:0},simultaneous=True))
    positive_D = -s*s
    positive_V = positive_D-(S.Matrix([s,U]).T*S.Matrix([S.diff(positive_D,s),S.diff(positive_D,U)]))[0]
    zero('the necessary square equation also has genuine positive cotangent data',positive_V-s*s)

    dependency_pins = {}
    for receipt_name,proofs in [
        ('receipt-projection-profile.json',['MIXED_PROJECTION_PROFILE.md']),
        ('receipt-stable-hypersurface.json',['STABLE_HYPERSURFACE_OBSTRUCTION.md']),
        ('receipt-square-axis.json',['P_GRAPH_SINGLETON_ENDPOINT.md']),
        ('receipt-first-slice.json',['FIRST_SLICE_INVERTIBILITY.md']),
        ('receipt-global-value.json',['GLOBAL_VALUE_ORBITS.md'])]:
        retained = json.loads((HERE/receipt_name).read_text())
        for name in proofs:
            dependency_pins[str((HERE/name).relative_to(ROOT))] = retained['source_sha256'][name]
    tang = json.loads((HERE/'receipt-stable-hypersurface.json').read_text())
    for name,expected in tang['dependencies'].items():
        if 'tang-hc4-audit' in name:
            dependency_pins[name] = expected
    for name,expected in dependency_pins.items():
        check('retained written input pin '+name,
              hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==expected)
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                       for name in ['verify_fixed_projection.py','FIXED_POTENTIAL_PROJECTION.md']},
        dependencies=dependency_pins,
        written_not_finite_exhaustion=[
            'a nonlinear univariate boundary closes the rank-one mixed cotangent projection',
            'fixed Meng--Yang gradient-base realizations must be transverse to the cotangent output plane',
            'centered rank-two gradient-dependent physical bases transfer the ambient critical value to at most one source point',
            'remaining transverse projections satisfy a global pulled-back Euler square identity'],
        new_external_classification_imported=False,
        rank_one_boundary_argument_requires_closed_embedding=False,
        shifted_pure_Keller_fibres_asserted_closed=False,
        variable_eta_coefficient_rows_asserted_closed=False,
        rank_one_abstract_embeddings_asserted_impossible=False,
        all_transverse_fixed_potential_projections_asserted_closed=False,
        arbitrary_physical_bases_asserted_gradient_projections=False,
        necessary_square_control_asserted_full_embedding=False,
        critical_value_count_asserted_as_singleton_zero_properness=False,
        nonlinear_Hessian_invariance_used=False,marked_collision_transfer_asserted=False,
        universal_SL2_row_equation_solved=False,unit_ideal_plane_Keller_endpoint_resolved=False,
        singleton_zero_properness_resolved=False,unrestricted_HC4_or_JC2_resolved=False,
        novelty_priority_asserted=False,formal_verification=False,independent_replay=False,external_review=False,
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-started,3))
    output = args.output or HERE/'receipt-fixed-projection.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(outcome=result['outcome'],check_count=result['check_count'],elapsed_seconds=result['elapsed_seconds'])))


if __name__=='__main__':
    main()
