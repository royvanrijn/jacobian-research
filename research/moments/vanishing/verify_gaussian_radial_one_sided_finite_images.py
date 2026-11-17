#!/usr/bin/env python3
"""Bounded original Wick and relative-source controls for radial one-sided phases."""
from math import factorial
import argparse
import hashlib
import json
import platform
from pathlib import Path
import sympy as s
if not __debug__:
    raise RuntimeError("Run without -O: exact assertions must remain enabled")

S,Y,t=s.symbols('S Y t')
N=4

def tr(f,n=N):
    return s.expand(sum(co*t**int(power[0]) for power,co in s.Poly(s.expand(f),t).terms()
                        if int(power[0])<=n))
def multiply(f,g,n=N):
    return tr(f*g,n)
def exponential(f,n=N):
    answer=s.Integer(1)
    power=s.Integer(1)
    for j in range(1,n+1):
        power=multiply(power,f,n)
        answer+=power/s.factorial(j)
    return tr(answer,n)
def gamma(f):
    return s.expand(sum(co*s.factorial(power[0]) for power,co in s.Poly(s.expand(f),S).terms()))
def radial_wick(f):
    return s.expand(sum(co*s.factorial(i)*s.factorial2(j-1)
                        for (i,j),co in s.Poly(s.expand(f),S,Y).terms() if j%2==0))

def controls(g,B,a):
    d=1-2*a*t
    P=g+B*Y+a*Y**2
    unit=s.Integer(1)
    power=s.Integer(1)
    for j in range(1,N+1):
        power=s.expand(power*P)
        unit+=radial_wick(power)*t**j/s.factorial(j)
    invd=sum((2*a*t)**j for j in range(N+1))
    H=sum(s.binomial(2*j,j)*a**j*t**j/2**j for j in range(N+1))
    E=exponential(tr(t*g+multiply(t*t*B*B,invd)/2))
    completed=multiply(H,gamma(E))
    assert s.expand(unit-completed)==0
    # A nontrivial scalar gauge checks the actual relative relations and
    # parameter connection without presuming its normalized value is one.
    p=1+t
    psi=2*t+t*t
    invp=sum((-t)**j for j in range(N+1))
    Hg=multiply(H,invp)
    Eg=exponential(tr(t*g+multiply(t*t*B*B,invd)/2-psi))
    phi=-S+tr(t*g+multiply(t*t*B*B,invd)/2-psi)
    def L(q,n=N):
        return multiply(Hg,gamma(multiply(q,Eg,n)),n)
    for j in range(1,5):
        f=S**j
        assert tr(L(s.diff(f,S)+s.diff(phi,S)*f))==0
    endpoint=L(s.diff(phi,S))
    endpoint_expected=-multiply(Hg,Eg.subs(S,0))
    assert s.expand(endpoint-endpoint_expected)==0
    assert endpoint!=0  # The absolute quotient would incorrectly kill this.
    logarithmic=tr(sum(a*(2*a*t)**j for j in range(N+1))-invp)
    for q in [s.Integer(1),S,S**3,t*S+S**2]:
        lhs=tr(s.diff(L(q),t),N-1)
        rhs=L(tr(s.diff(q,t)+(s.diff(phi,t)+logarithmic)*q,N-1),N-1)
        assert s.expand(lhs-rhs)==0
    return {'g':str(g),'B':str(B),'a':a,'original_unit_jet':str(unit),
            'relative_spatial_relations':4,'relative_parameter_relations':4,
            'nonzero_endpoint_jet':str(endpoint),'absolute_endpoint_erasure_rejected':True}

rows=[controls(S**2,2*S,1),controls(S**3+S,2+S**2,0),
      controls(3+2*S,s.Integer(4),1),controls(s.Integer(3),s.Integer(4),0)]
a,gD,CD=s.symbols('a gD CD')
d=1-2*a*t
top=t*gD+t*t*CD/(2*d)
assert s.cancel(top-t*(gD+t*(CD/2-2*a*gD))/d)==0
# The cancellation face has exactly order0 at infinity.
assert s.cancel(top.subs({a:1,gD:1,CD:4})-t/(1-2*t))==0

Z,W=s.symbols("Z W")
def original_wick(value):
    answer=0
    for (i,j,k),co in s.Poly(s.expand(value),Z,W,Y).terms():
        if i==j and k%2==0:
            answer+=co*s.factorial(i)*s.factorial2(k-1)
    return s.expand(answer)

def one_sided_control(b,g0):
    remainder=Z**5*W**3*Y**4+Z**2
    phase=g0+b*Y+remainder
    order=6
    markers=[s.Integer(1),W**2,Y*W**2,W**4]
    powers=[s.Integer(1)]
    for _ in range(order):
        powers.append(s.expand(powers[-1]*phase))
    factor=exponential(g0*t+b*b*t*t/2,order)
    controls=[]
    for marker in markers:
        D=s.Poly(marker,Z,W,Y).total_degree()
        normalized=0
        for n in range(D//2+1):
            shifted=marker.subs(Y,Y+b*t)*remainder.subs(Y,Y+b*t)**n
            normalized+=t**n*original_wick(shifted)/s.factorial(n)
        normalized=s.expand(normalized)
        actual=sum(original_wick(marker*powers[m])*t**m/s.factorial(m)
                   for m in range(order+1))
        predicted=multiply(factor,normalized,order)
        assert s.expand(actual-predicted)==0
        bound=D//2 if b==0 else D+5*(D//2)
        assert normalized==0 or s.degree(normalized,t)<=bound
        if marker==1:
            assert normalized==1
        if marker==W**2:
            expected=362*t+720*b*b*t**3+120*b**4*t**5
            assert s.expand(normalized-expected)==0
        controls.append({'marker':str(marker),'normalized_exact_period':str(normalized),
                         'sufficient_degree_bound':bound})
    return {'core_B':b,'core_g':g0,'gamma':2,'degree_Y_R':4,
            'R':'Z^5 W^3 Y^4+Z^2','original_orders':[0,order],
            'markers':controls}

one_sided=[one_sided_control(0,0),one_sided_control(2,3)]
result={
    'status':'PASS',
    'scope':'Bounded exact controls only; endpoint comparison and all-degree rigidity are written proofs',
    'radial_orders':[0,N],'radial_cases':rows,
    'relative_spatial_orders':[0,N],'relative_parameter_orders':[0,N-1],
    'generic_top_identity':True,'cancelling_top_order':0,
    'one_sided_original_controls':one_sided,
    'relative_endpoint_retained':True,
}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,help='Write the exact JSON receipt to this path')
    args=parser.parse_args()
    here=Path(__file__).resolve()
    repository=here.parents[3]
    names=[
        'GAUSSIAN_RADIAL_ONE_SIDED_FINITE_IMAGES.md',
        'GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md',
        'GAUSSIAN_AFFINE_ISOTROPIC_CHART.md',
        'GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md',
        'GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md',
        'MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md',
    ]
    canonical=[here]+[here.with_name(name) for name in names]
    imports=[repository/'unpublished/gaussian-boundary'/name for name in [
        'polynomial_stationary_balls.tex','seed_period_comparison.tex','analytic_core.tex',
        'sources/laurent_primary_manifest.json',
    ]]
    paths=canonical+imports
    pins={str(path.relative_to(repository)).replace('\\','/'):hashlib.sha256(path.read_bytes()).hexdigest()
          for path in paths}
    receipt={
        'status':'passed','checked_on':'2026-11-17',
        'environment':{'python':platform.python_version(),'sympy':s.__version__},
        'sha256':pins,'checks':result,
        'independent_replay':False,'formal_verification':False,'external_review':False,
        'unrestricted_result_from_finite_prefix':False,
        'scope':'Bounded exact original Wick and relative-source controls; the all-degree theorem and relative endpoint frame are in the pinned written proof with retained comparison imports.',
    }
    encoded=json.dumps(receipt,indent=2)+'\n'
    if args.output is not None:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(encoded,encoding='utf-8')
    print(encoded,end='')

if __name__=='__main__':
    main()
