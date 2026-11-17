#!/usr/bin/env python3
"""Exact transcription controls; no recurrence fitting or infinite-series proof."""
from pathlib import Path
import hashlib,json,math
import sympy as S

HERE=Path(__file__).resolve().parent
x,y,w,X,t,z,v,s,U=S.symbols('x y w X t z v s U')

def action(lam,f):
    lam=S.sympify(lam);f=S.sympify(f)
    assert not lam.has(S.Float) and not f.has(S.Float)
    assert S.Poly(f,x,y,w).total_degree()<=24
    return S.expand(sum(c*S.diff(f,x,a[0]) for a,c in S.Poly(lam,X).terms()))

def iterate(lam,f,k):
    assert isinstance(k,int) and 0<=k<=12
    for _ in range(k):f=action(lam,f)
    return f

def height(lam,f):
    j=0
    while True:
        f=action(lam,f)
        if f==0:return j
        j+=1
        assert j<=12

def main():
    phi=X**2+X**6+X**24
    profiles=[];F=1
    for m in range(1,7):
        h=height(phi,x**m);assert h==m//2
        assert iterate(phi,x**m,m)==0
        assert iterate(phi,x**(2*m),m)==math.factorial(2*m)
        R=sum(iterate(phi,x**m,k)*z**(m-k)/S.factorial(k) for k in range(h+1))
        F+=R*t**m/S.factorial(m)
        profiles.append(dict(m=m,height=h,deficit=m-h,sharp_marker_degree=m))
    F=S.expand(F)
    formula=S.series(S.exp(t*z*x+z*t*t+z**5*t**6),t,0,7).removeO()
    assert S.expand(F-formula)==0
    assert all(i<=2*j for (i,j,k),c in S.Poly(F,t,z,x).terms())
    # A later coefficient actually participates in the tested source.
    assert action(phi,x**6)==30*x**4+720
    unit=1+X**4+X**22;inverse=1-X**4+X**8-X**12
    assert action(unit,action(inverse,x**12))==x**12
    assert [height(X**2,x**m) for m in range(1,7)]==[v['height'] for v in profiles]
    good=S.exp(t*z*x+z*t*t)
    assert S.simplify(S.diff(good,t)-z*(x+2*t)*good)==0
    # The concise degree-three phase still gives the same actual heights.
    P=x+y**3+w**3;higher=[]
    for m in [2,3,4]:
        h=height(phi,P**m);assert h==m//2
        assert iterate(phi,x**m*P**m,m)==math.factorial(2*m)
        higher.append(dict(m=m,height=h,deficit=m-h))
    # Finite Bell leading terms used by the universal z-valuation proof.
    f=v*v+v**6;bell=S.Integer(1);bells=[]
    for j in range(1,5):
        bell=S.expand(S.diff(bell,v)+S.diff(f,v)*bell/z)
        scaled=S.Poly(S.expand(z**j*bell),z)
        assert scaled.coeff_monomial(1)==S.expand(S.diff(f,v)**j)
        bells.append(dict(order=j,leading_z_power=-j))
    # An exact characteristic equation in the algebraic quadratic control.
    E=S.exp(v*v/z)
    assert S.simplify(z*S.diff(E,v)-2*v*E)==0
    assert S.diff(v*v,v)-2*v==0
    # Multiple tied leading valuations give a NONZERO characteristic polynomial.
    b0=S.Integer(1);b1=S.diff(f,v)/z;b2=S.diff(b1,v)+S.diff(f,v)*b1/z
    expression=S.expand(z*(v+1)*b0+z*z*(v*v-1)*b1+z**3*(v**3+2)*b2)
    characteristic=(v+1)+(v*v-1)*U+(v**3+2)*U**2
    assert characteristic!=0 and S.expand(expression.coeff(z,1)-characteristic.subs(U,S.diff(f,v)))==0
    # The proposed geometric pullback is not multiplicative.
    phase=x*x
    a1=S.simplify(S.exp(-s*phase)*S.diff(S.exp(s*phase),x))
    a2=S.simplify(S.exp(-s*phase)*S.diff(S.exp(s*phase),x,2))
    assert S.expand(a2-a1*a1)==2*s
    guards=0
    for lam,f in [(X,S.Float(1)),(X,x**25)]:
        try:action(lam,f)
        except AssertionError:guards+=1
        else:raise AssertionError('guard failure')
    deps=['verify_scalar_time_gaps.py','verify_linear_suspension.py','verify_formal_symbol_source_rigidity.py','verify_algebraic_formal_symbols.py']
    report=dict(completed=True,checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        supporting_checker_hashes={f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in deps},
        actual_profiles=profiles,reversed_source_time_coefficients_checked=7,
        higher_jet_participating_output=str(action(phi,x**6)),sharp_cutoff='1',height_gap_fraction='1/2',
        unit_height_comparison_pass=True,unit_scalar_equation_pass=True,
        concise_cubic_profiles=higher,bell_leading_controls=bells,
        characteristic_tie_control_pass=True,geometric_pullback_multiplicativity_defect=str(S.expand(a2-a1*a1)),
        guard_refusals=guards,
        infinite_lacunarity_proved_by_prefix=False,nonexistence_proved_by_failed_fitting=False,
        unmodified_reversed_equation_follows_from_eventuality=False,
        formal_unit_invariance_of_reversed_equations=False,
        height_preserving_representative_finiteness_disproved=False,
        general_actual_gap_decided=False,full_exceptional_MZ_decided=False,
        unrestricted_frontier_endpoint_changed=False,independent_replay=False,
        external_review=False,formal_verification=False,novelty_claimed=False)
    assert guards==2
    (HERE/'reversed_time_obstruction_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(completed=True,actual_profiles=6,source_coefficients=7,concise_profiles=3,bell_controls=4,pullback_defect='2*s',guards=guards)))

if __name__=='__main__':main()
