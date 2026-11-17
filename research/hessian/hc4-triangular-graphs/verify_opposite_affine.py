"""Exact controls for the all-degree opposite affine-p closure.

Pole orders and the Tang/Chau/Orevkov exclusion are written arguments,
not inferred from these finite controls. Imports are explicitly pinned.
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

    s,x,p,q,r,T,t,tau = S.symbols('s x p q r T t tau')
    L,N,c,rho,e0 = S.symbols('L N c rho e0',nonzero=True)
    M,b0,b1,A0,B0 = S.symbols('M b0 b1 A0 B0')
    a,b = S.Function('a')(s),S.Function('b')(s)
    u = 1+x*s
    AA = u**3*p+3*x*u**2*q-x**3*r
    BB = s*s*u*(4+3*x*s)*p+(s+3*x*s*s*(4+3*x*s))*q+(2*x-3*x*x*s)*r
    psi = L*AA**2+M*AA+N*BB
    pc = (x-b)/a
    rr = (u**3*pc+3*x*u**2*q-T)/x**3
    cc = M-N*(2-3*x*s)/x**2
    ee = N*u*(u+1)*(x-b)/(a*x**2)
    dd = 2*N*(3+2*x*s)/x
    F = L*T*T+cc*T+ee+dd*q
    zero('rational computation chart reproduces the original polynomial potential',
         psi.subs({p:pc,r:rr},simultaneous=True)-F)
    # The chart is used only on its stated localization, not as an embedding.
    inv_duals = S.Matrix([pc,q,rr])
    zero('relative inverse-chart determinant is exactly reciprocal a*x^3',
         inv_duals.jacobian((x,T,q)).det()-1/(a*x**3))
    zero('at x=0 the native q derivative is N*s',S.diff(psi,q).subs(x,0)-N*s)
    zero('coefficient d has the unique root x=-3/(2s)',dd.subs(x,-S.Rational(3,2)/s))
    zero('remaining critical equation has nonzero coefficient -6N/x^2',S.diff(dd,x)+6*N/x**2)
    xc = -S.Rational(3,2)/s
    Tc = (13*N/(2*xc**2)-M)/(2*L)
    qc = 17*Tc/(12*xc)+(xc-2*b)/(24*a*xc)
    critical = {x:xc,T:Tc,q:qc}
    zero('all three critical equations vanish on the entire explicit curve',
         S.Matrix([S.diff(F,z) for z in (x,T,q)]).subs(critical,simultaneous=True))
    D = N*s*(3+2*s*b)/(18*a)-(M-26*N*s*s/9)**2/(4*L)
    zero('critical value is the claimed rational function',F.subs(critical,simultaneous=True)-D)
    zero('base-gradient on the whole critical curve is D prime',
         S.diff(F,s).subs(critical,simultaneous=True)-S.diff(D,s))
    relative_H = S.hessian(F,(x,T,q))
    zero('relative critical Hessian determinant in the computation chart',
         relative_H.det()+72*L*N*N/x**4)
    zero('native relative Hessian determinant factor',
         (-72*L*N*N/x**4)*(a*x**3)**2+72*L*N*N*a*a*x*x)
    zero('native full critical determinant has the -162 factor',
         (-72*L*N*N*a*a*x*x).subs(x,xc)+162*L*N*N*a*a/s**2)
    # Verify the Schur identity using the tangent of the explicit critical
    # curve; no expansion of the original high-degree native determinant.
    rel_on = relative_H.subs(critical,simultaneous=True)
    cross_on = S.Matrix([S.diff(F,s,z) for z in (x,T,q)]).subs(critical,simultaneous=True)
    tangent = S.Matrix([S.diff(value,s) for value in (xc,Tc,qc)])
    zero('critical-curve derivative solves the Schur linear system',rel_on*tangent+cross_on)
    zero('Schur complement is exactly D double prime',
         S.diff(F,s,2).subs(critical,simultaneous=True)+(cross_on.T*tangent)[0]-S.diff(D,s,2))

    # A genuine s-dependent relative automorphism and its inverse. The
    # native example s*p+q*r makes f_s nonzero on its critical line.
    z = S.Matrix([p,q,r])
    native = s*p+q*r+p**3
    phi = S.Matrix([p+s*q*q,q+s*r*r,r])
    inverse = S.Matrix([p-s*(q-s*r*r)**2,q-s*r*r,r])
    zero('s-dependent relative control inverse left',inverse.subs(dict(zip(z,phi)),simultaneous=True)-z)
    zero('s-dependent relative control inverse right',phi.subs(dict(zip(z,inverse)),simultaneous=True)-z)
    zero('s-dependent relative control has Jacobian one',phi.jacobian(z).det()-1)
    composed = native.subs(dict(zip(z,phi)),simultaneous=True)
    chain_correction = sum(S.diff(native,zi).subs(dict(zip(z,phi)),simultaneous=True)*S.diff(vi,s)
                           for zi,vi in zip(z,phi))
    zero('base-gradient change belongs literally to the pulled critical ideal',
         S.diff(composed,s)-S.diff(native,s).subs(dict(zip(z,phi)),simultaneous=True)-chain_correction)
    # Native critical points: s=-3t^2,p=t,q=r=0, so f_s=t.
    param = {s:-3*t*t,p:t,q:0,r:0}
    pulled = inverse.subs(param,simultaneous=True)
    pulled_param = {s:-3*t*t,**dict(zip(z,pulled))}
    zero('whole critical target-line map survives the relative automorphism',
         S.Matrix([S.diff(composed,s)-t]+[S.diff(composed,zi) for zi in z]).subs(pulled_param,simultaneous=True))

    # The all-degree pole argument is in the proof. These identities check
    # its coefficients and each necessary packet without a degree search.
    E = -(M-26*N*s*s/9)**2/(4*L)
    D3 = S.cancel(D.subs({a:c*s**3,b:b0+b1*s},simultaneous=True))
    zero('c*s^3 packet has unavoidable nonzero order-three D second-derivative term',
         S.limit(s**3*(S.diff(D3,s,2)-N/(c*s**4)),s,0)-2*N*b0/(9*c))
    bI = b0+b1*s+169*c*N*s**3/(9*L)
    q2 = N*(b1+13*c*M/L)/(9*c)
    q1 = N*b0/(9*c)
    q0 = N/(6*c)-M*M/(4*L)
    DI = S.cancel(D.subs({a:c*s,b:bI},simultaneous=True))
    zero('packet I cancels every higher critical-value power',DI-(q2*s*s+q1*s+q0))
    zero('packet I exact critical target-line function',S.diff(DI,s)-(2*q2*s+q1))
    zero('packet I constant critical determinant',
         -162*L*N*N*c*c*S.diff(DI,s,2)+324*L*N*N*c*c*q2)
    zero('packet I unique inverse except its excluded source zero',
         (2*q2*s+q1).subs(s,(tau-q1)/(2*q2))-tau)
    check('packet I omitted target value forces precisely the deleted source point',
          S.cancel(((tau-q1)/(2*q2)).subs(tau,q1))==0)
    BII = -M*M/(4*L)+(N/(6*c)-e0)/rho**2
    numerator = (18*c/N)*(e0+(A0*s+BII)*(s-rho)**2-E*(s-rho)**2)-3
    zero('packet II b numerator is divisible by s',numerator.subs(s,0))
    bII = S.cancel(numerator/(2*s))
    check('packet II b is a polynomial of degree at most five',bII.is_polynomial(s) and S.degree(bII,s)<=5)
    zero('packet II b at zero is the claimed exact coefficient',
         bII.subs(s,0)-(9*c*rho*rho*(A0+2*e0/rho**3)/N-3/rho))
    DII = S.cancel(D.subs({a:c*s*(s-rho)**2,b:bII},simultaneous=True))
    zero('packet II critical value has no omitted polynomial terms',DII-(e0/(s-rho)**2+A0*s+BII))
    zero('packet II exact critical target-line function',S.diff(DII,s)-(A0-2*e0/(s-rho)**3))
    zero('packet II constant critical determinant',
         -162*L*N*N*c*c*(s-rho)**4*S.diff(DII,s,2)+972*L*N*N*c*c*e0)
    fibre_poly = (tau-A0)*(s-rho)**3+2*e0
    zero('packet II fibre equation is exactly the cubic numerator',
         (tau-S.diff(DII,s))*(s-rho)**3-fibre_poly)
    zero('packet II cubic discriminant records only the omitted target value',
         S.discriminant(fibre_poly,s)+108*e0**2*(tau-A0)**2)
    zero('packet II pole point is never a fibre root',fibre_poly.subs(s,rho)-2*e0)
    zero('packet II omitted target value gives no root at all',fibre_poly.subs(tau,A0)-2*e0)
    zero('packet II zero source point is removed at exactly one other target',
         fibre_poly.subs(s,0)+rho**3*(tau-A0-2*e0/rho**3))
    check('packet II special attained target differs from the omitted target',2*e0/rho**3!=0)

    # Abstract cotangent factor: every residual plane preimage supplies
    # exactly one pair of cotangent coordinates for every other output.
    U,V,P,Q = S.symbols('U V P Q')
    FF,GG,HH = (S.Function(name)(U,V) for name in ('FF','GG','HH'))
    cotangent = HH+P*FF+Q*GG
    residual_J = S.Matrix([FF,GG]).jacobian((U,V))
    zero('Tang cotangent Hessian is the square of the plane Keller determinant',
         S.hessian(cotangent,(U,V,P,Q)).det()-residual_J.det()**2)
    target1,target2 = S.symbols('target1 target2')
    unique = residual_J.T.inv()*S.Matrix([target1-S.diff(HH,U),target2-S.diff(HH,V)])
    zero('each residual plane preimage has a unique cotangent lift',
         S.Matrix([S.diff(cotangent,U),S.diff(cotangent,V)]).subs({P:unique[0],Q:unique[1]},simultaneous=True)
         -S.Matrix([target1,target2]))

    manifest_name = 'sources/opposite-affine-import.json'
    manifest = json.loads((HERE/manifest_name).read_text())
    dependencies = {}
    for imported in manifest['imports']:
        name = imported['retained_file_from_research_root']
        digest = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained primary byte pin '+imported['author'],digest==imported['sha256'])
        dependencies[name] = digest
    check('Chau import is specifically Theorem 4',manifest['imports'][0]['location']=='Theorem 4')
    check('degree-three import retains exact Orevkov location and prior registry',
          'Theorem 1.1' in manifest['imports'][1]['location']
          and manifest['imports'][1]['existing_registered_degree_three_result']=='PCC3')
    tang_pins = json.loads((ROOT/'hessian/tang-hc4-audit/structural-proof-receipt.json').read_text())['source_proof_and_replay_sha256']
    for name in ['hessian/tang-hc4-audit/STRUCTURAL_REVIEW.md','hessian/tang-hc4-audit/source-v1.pdf']:
        dependencies[name] = hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        check('retained Tang input '+name,dependencies[name]==tang_pins[name])
    result = dict(outcome='PASS',check_count=len(checks),checks=checks,
        source_sha256={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
                      for name in [Path(__file__).name,'OPPOSITE_AFFINE_ENDPOINT.md',manifest_name]},
        dependencies=dependencies,
        written_not_finite_exhaustion=['whole critical scheme and its preserved target-line fibres',
            'arbitrary-degree finite-pole and infinity-order classification',
            'line projection avoids the plane nonproper locus generically by Chau Theorem4',
            'actual geometric degrees one and three excluded by classical plane theorems',
            'whole opposite affine-p family reduces exactly to the base-only unit-ideal endpoint'],
        every_nonzero_a_opposite_affine_p_graph_forbids_unit_Hessian=True,
        retained_three_dual_automorphism_degrees_unrestricted=True,
        base_only_endpoint_is_exactly_unit_ideal_plane_Keller_existence=True,
        critical_determinant_packets_alone_asserted_sufficient=False,
        unit_ideal_plane_Keller_endpoint_resolved=False,singleton_zero_properness_resolved=False,
        all_no_simple_boundary_root_y_graphs_closed=False,
        unrestricted_HC4_or_JC2_resolved=False,
        imported_theorems=['HCTANG2','Chau arXiv:0711.3894v1 Theorem4',
                           'Orevkov Theorem1.1 / PCC3','classical birational Keller criterion'],
        runtime=dict(python=platform.python_version(),sympy=S.__version__),
        elapsed_seconds=round(time.monotonic()-start,3),
        formal_verification=False,independent_replay=False,external_review=False)
    rendered = json.dumps(result,indent=2)+'\n'
    if args.output:
        args.output.write_text(rendered)
    print(rendered,end='')


if __name__=='__main__':
    main()
