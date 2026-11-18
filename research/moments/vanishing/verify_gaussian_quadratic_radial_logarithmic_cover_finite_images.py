"""Exact controls for the fixed-log quadratic radial Gaussian theorem.

These controls do not computationally prove the complete marked-model
comparison or its chain-level estimates. Those are the written proof.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import platform
import sympy as sp

if not __debug__:
    raise RuntimeError('Exact assertions must be enabled')

s,y,t,z,w=sp.symbols('s y t z w')
r=1/t
A=sp.Function('A')(y)
H=sp.Function('H')(y)
G0=sp.Function('G0')(y)
psi=sp.Function('psi')(t)

def zero(expr):
    expr=sp.cancel(sp.expand(expr))
    if expr!=0:expr=sp.simplify(expr)
    assert expr==0,expr

def source_controls():
    P=A*s*s+H*s+G0
    q=1-t*H-2*t*A*s
    def DS(value):return sp.diff(value,s)+(-1+t*(2*A*s+H))*value
    def DY(value):return sp.diff(value,y)+(-y+t*sp.diff(P,y))*value
    def DT(value):return sp.diff(value,t)+(P-sp.diff(psi,t))*value
    reductions=[sp.Integer(1),s]
    for degree in range(2,7):
        reductions.append(sp.expand(((1-t*H)*reductions[-1]
                           -(degree-1)*reductions[-2])/(2*t*A)))
    def reduce(value):
        return sp.expand(sum(coefficient*reductions[power]
          for (power,),coefficient in sp.Poly(sp.expand(value),s).terms()))
    def basis(value):
        value=reduce(value)
        constant,linear=value.coeff(s,0),value.coeff(s,1)
        return (sp.cancel(constant+linear*(1-t*H)/(2*t*A)),
                sp.cancel(-linear/(2*t*A)))
    for degree in range(1,6):
        zero(reduce(DS(s**degree)))
        zero(reduce(q*s**degree)-degree*reductions[degree-1])
    G=t*G0-y*y/2-psi
    x=sp.sqrt(t)*(H-r)/(2*sp.sqrt(A))
    F=G-x*x
    m=sp.sqrt(t)*sp.sqrt(A)
    zero(reduce(DY(q)-sp.diff(G,y)*q))
    zero(reduce(DT(q)-sp.diff(G,t)*q))
    ay,ky=basis(DY(sp.Integer(1)))
    at,kt=basis(P-sp.diff(psi,t))
    zero(ay+sp.diff(A,y)/(2*A)-sp.diff(F,y))
    zero(at+1/(2*t)-sp.diff(F,t))
    zero(m*ky+sp.diff(x,y))
    zero(m*kt+sp.diff(x,t))
    MY=sp.Matrix([[sp.diff(F,y),-sp.diff(x,y)],[0,sp.diff(G,y)]])
    MT=sp.Matrix([[sp.diff(F,t),-sp.diff(x,t)],[0,sp.diff(G,t)]])
    for entry in sp.diff(MY,t)-sp.diff(MT,y)+MY*MT-MT*MY:zero(entry)
    return dict(A_H_G0='arbitrary symbolic functions of Y',psi='arbitrary symbolic time function',
                exact_relative_relations=5,full_Y_and_t_connection=True,
                endpoint_retained=True,curvature_entries=4)

def wick(value):
    total=0
    for (ns,ny),coefficient in sp.Poly(sp.expand(value),s,y).terms():
        if ny%2==0:total+=coefficient*sp.factorial(ns)*sp.factorial2(ny-1)
    return sp.expand(total)

def original_controls():
    fixtures=[(y*y+1,y**3-y,y**5+y*y+2),
              ((y*y-1)**2,(y-1)**3,y**4-y),
              (y*(y*y-1),sp.Integer(0),y**6+y),
              (y*y*(y**3-y-1),y+1,y**5+2)]
    markers=[sp.Integer(1),s,y,y*y,s*y,s*s]
    count=0
    for aa,hh,gg in fixtures:
        pp=aa*s*s+hh*s+gg
        powers=[sp.Integer(1)]
        for _ in range(4):powers.append(sp.expand(powers[-1]*pp))
        def coefficient(a0,a1,order):
            value=wick(a0*powers[order])/sp.factorial(order)
            if order:value+=wick(a1*powers[order-1])/sp.factorial(order-1)
            return value
        for marker in markers:
            for order in range(5):
                zero(coefficient(sp.diff(s*marker,s)-s*marker,
                                 (2*aa*s+hh)*s*marker,order))
                zero(coefficient(sp.diff(marker,y)-y*marker,
                                 sp.diff(pp,y)*marker,order))
                count+=2
        for marker in [sp.Integer(1),y,y*y]:
            for order in range(5):
                zero(coefficient(marker,-(2*aa*s+hh)*marker,order)
                     -wick(marker*gg**order)/sp.factorial(order))
                count+=1
    assert wick(-1)==-1
    return dict(phases=len(fixtures),original_markers=len(markers),orders=[0,4],
                exact_controls=count,nonaffine_boundary=True,
                positive_genus_and_disconnected_cover_fixtures=True,
                rejected_absolute_source='D_S(1) has nonzero endpoint')

def eligible(polynomial):
    polynomial=sp.Poly(polynomial,y)
    if polynomial.is_zero or polynomial.degree()<2:return False
    return all(multiplicity<=2 for factor,multiplicity in sp.factor_list(polynomial.as_expr(),y)[1])

def geometry_controls():
    fixtures=[y*y+1,(y*y-1)**2,y*(y*y-1),y*y*(y**3-y-1),
              (y*y-1)**2*(y**3-y-1)]
    report=[]
    for aa in fixtures:
        assert eligible(aa)
        factors=sp.factor_list(aa,y)[1]
        squarefree=sp.prod(f for f,mul in factors if mul==1)
        doubled=sp.prod(f for f,mul in factors if mul==2)
        B_degree=sp.degree(squarefree,y)
        D_degree=sp.degree(doubled,y)
        assert sp.degree(aa,y)==B_degree+2*D_degree
        components=2 if B_degree==0 else 1
        genus=0 if components==2 else (B_degree-1)//2
        report.append(dict(A=str(aa),degree=int(sp.degree(aa,y)),
                           normalized_components=components,genus_each=int(genus),
                           B_degree=int(B_degree),D_degree=int(D_degree)))
    for bad in [sp.Integer(0),sp.Integer(3),y+1,y**3,(y*y+1)**3]:
        assert not eligible(bad),bad
    j=sp.symbols('j',integer=True,nonnegative=True)
    even_infinity_order=j-1  # degree A=2j+2
    odd_infinity_order=2*j  # degree A=2j+3
    assert (even_infinity_order+1).is_nonnegative
    assert odd_infinity_order.is_nonnegative
    # Exact rational uniformization for the two-simple-root quadratic.
    g,a,m=sp.symbols('g a m',nonzero=True)
    YY=m+a*(z+1/z)/2
    VV=g*a*(z-1/z)/2
    zero(VV*VV-g*g*((YY-m)**2-a*a))
    zero(sp.diff(YY,z)/VV-1/(g*z))
    # Elliptic example: at infinity and at its simple zero the volume is regular.
    YY=w**-2
    VV=w**-3*sp.sqrt(1-w**4)
    zero(sp.diff(YY,w)/VV+2/sp.sqrt(1-w**4))
    YY=w*w
    VV=w*sp.sqrt(w**4-1)
    zero(sp.diff(YY,w)/VV-2/sp.sqrt(w**4-1))
    # A perfect square has two normalized components with only logarithmic roots.
    zero(sp.diff(y,y)/(y*y-1)-sp.Rational(1,2)*(1/(y-1)-1/(y+1)))
    # The volume may have a REGULAR-PHASE logarithmic puncture.
    # It must be avoided by finite-value detours, not treated as a polar tail.
    aa=y*y
    hh=aa*y*y
    gg=aa*y**4/4
    ff=t*gg-y*y/2-t*(hh-1/t)**2/(4*aa)
    zero(ff+1/(4*t*aa))
    assert sp.limit(ff,y,sp.oo)==0
    assert sp.diff(ff,y)!=0
    return dict(normalizations=report,
                all_degree_infinity_orders=['even degree2j+2: j-1>=-1','odd degree2j+3: 2j>=0'],
                finite_root_orders={'simple':'regular on branch coordinate','double':'logarithmic'},
                rejection_controls=['zero/constant/affine excluded from this branch','triple root rejected'],
                exact_genus_zero_uniformizer=True,positive_genus_volume_fixtures=True,
                regular_phase_log_puncture_fixture='A=Y^2,H=Y^4,G0=Y^6/4 gives F=-r/(4Y^2), regular at infinity')

def sector_controls():
    x=sp.symbols('x')
    EE=sp.Function('E')(x)
    jp=sp.exp(-x*x)*(sp.I*sp.sqrt(sp.pi)/2-EE)
    jm=sp.exp(-x*x)*(-sp.I*sp.sqrt(sp.pi)/2-EE)
    for jj in [jp,jm]:zero((sp.diff(jj,x)+2*x*jj+1).subs(sp.diff(EE,x),sp.exp(x*x)))
    zero(jp-jm-sp.I*sp.sqrt(sp.pi)*sp.exp(-x*x))
    zero(jm.subs(x,-x).subs(sp.Function('E')(-x),-EE)+jp)
    degree=sp.symbols('degree',integer=True,nonnegative=True)
    assert (sp.Rational(1,2)+degree).is_positive
    return dict(ODE='J_prime=-2xJ-1',exact_jump='i sqrt(pi) exp(-x²)',
                deck_descent='J_minus(-x)=-J_plus(x)',
                arbitrary_polynomial_unit_normalization_order='1/2+deg p>0')

def finite_original_marker_controls():
    z,w = sp.symbols("z w")
    remainder = z*z*y+z**3*w
    gamma = 2
    markers = [sp.Integer(1),z,w,y,w*w,w**3,y*w*w]
    order = 6

    def full_wick(value):
        total = 0
        for (nz,nw,ny),coefficient in sp.Poly(sp.expand(value),z,w,y).terms():
            if nz==nw and ny%2==0:
                total += coefficient*sp.factorial(nz)*sp.factorial2(ny-1)
        return sp.expand(total)

    report = []
    for cc,bb in [(0,0),(3,2)]:
        phase = remainder+cc+bb*y
        powers = [sp.Integer(1)]
        for _ in range(order):
            powers.append(sp.expand(powers[-1]*phase))
        exponential = [sp.Integer(1)]
        for degree in range(1,order+1):
            value = cc*exponential[degree-1]
            if degree>=2:
                value += bb*bb*exponential[degree-2]
            exponential.append(value/degree)
        fixture = []
        for marker in markers:
            D = sp.Poly(marker,z,w,y).total_degree()
            normalized = sp.expand(sum(t**n*full_wick(
                (marker*remainder**n).subs(y,y+bb*t))/sp.factorial(n)
                for n in range(D//gamma+1)))
            for degree in range(order+1):
                predicted = sum(exponential[degree-j]*normalized.coeff(t,j)
                                for j in range(degree+1))
                zero(full_wick(marker*powers[degree])/sp.factorial(degree)-predicted)
                if bb==cc==0 and degree>D//gamma:
                    zero(full_wick(marker*powers[degree]))
            if marker==1:
                zero(normalized-1)
            bound = D//gamma if bb==0 else D+2*(D//gamma)
            assert normalized==0 or sp.degree(normalized,t)<=bound
            fixture.append({"marker":str(marker),"normalized_period":str(normalized),
                            "valid_degree_bound":bound})
        report.append({"c":cc,"b":bb,"gamma":gamma,"degree_Y_R":1,
                       "R":"Z^2Y+Z^3W","original_moment_orders":[0,order],
                       "markers":fixture})
    return report


def original_support_controls():
    report = {}
    for degree in [5,6]:
        support = {(i,j) for i in range(degree//2+1)
                   for j in range(degree-2*i+1) if i>0 or j<=1}
        expected = ({(0,0),(0,1)}
                    | {(1,j) for j in range(4 if degree==5 else 5)}
                    | {(2,j) for j in range(2 if degree==5 else 3)}
                    | ({(3,0)} if degree==6 else set()))
        assert support==expected
        report[str(degree)] = [f"S^{i}Y^{j}" for i,j in sorted(support)]
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    here=Path(__file__).resolve()
    repository=here.parents[3]
    names=[
      'GAUSSIAN_QUADRATIC_RADIAL_LOGARITHMIC_COVER_FINITE_IMAGES.md',
      'GAUSSIAN_QUADRATIC_RADIAL_AFFINE_COEFFICIENT_FINITE_IMAGES.md',
      'GAUSSIAN_ANGULAR_CORE_FULL_IMAGE_INVARIANCE.md',
      'GAUSSIAN_AMPLITUDE_IDEALS_AND_WEIGHT_INVARIANCE.md',
      'GAUSSIAN_INVARIANT_CORE_BOUNDARY_REDUCTION.md',
      'GAUSSIAN_AFFINE_RADIAL_INVARIANT_CORE_FINITE_IMAGES.md',
      'GAUSSIAN_VARIABLE_RADIAL_QUADRATIC_FIBRE_FINITE_IMAGES.md',
      'GAUSSIAN_AFFINE_ISOTROPIC_CHART.md',
      'GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md',
      'MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md',
    ]
    files=[here]+[here.with_name(name) for name in names]
    primary_base=repository/'results/gaussian-boundary'
    files += [primary_base/name for name in [
      'analytic_core.tex','seed_period_comparison.tex',
      'sources/laurent_primary_manifest.json']]
    manifest=json.loads((primary_base/'sources/laurent_primary_manifest.json').read_text())
    for record in manifest['sources']:
        path=primary_base/record['file']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==record['sha256'],record['file']
        files.append(path)
    report=dict(status='passed',checked_on=datetime.now(timezone.utc).date().isoformat(),
       environment=dict(python=platform.python_version(),sympy=sp.__version__),
       checks=dict(source=source_controls(),original=original_controls(),
                   geometry=geometry_controls(),sectors=sector_controls(),
                   finite_original_markers=finite_original_marker_controls(),
                   original_support=original_support_controls(),
                   primary_manifest_consistency=dict(status='passed',sources=len(manifest['sources']))),
       scope=('Exact arbitrary-function source identities and bounded original-Wick, logarithmic '
              'normalization, Gaussian-sector, regular-puncture and finite-branch controls. '
              'The full fixed-curve chain comparison, open-chain filling and unrestricted '
              'scalar exclusion are the written proof; no finite-prefix extrapolation.'),
       independent_replay=False,formal_verification=False,external_review=False,
       sha256={path.relative_to(repository).as_posix():hashlib.sha256(path.read_bytes()).hexdigest()
               for path in files})
    out=json.dumps(report,indent=2)+'\n'
    if args.output:args.output.write_text(out,encoding='utf-8',newline='\n')
    print(out,end='')

if __name__=='__main__':main()
