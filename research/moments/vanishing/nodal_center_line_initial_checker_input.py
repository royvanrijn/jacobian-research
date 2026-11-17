#!/usr/bin/env python3
"""Exact auxiliary algebra; does not certify the analytic current argument."""
from pathlib import Path
from math import factorial
import hashlib
import json
import sympy as S


def wick(poly, variables):
    result = 0
    for powers, coefficient in S.Poly(S.expand(poly), *variables).terms():
        if any(power % 2 for power in powers):
            continue
        weight = 1
        for power in powers:
            weight *= factorial(power)//(2**(power//2)*factorial(power//2))
        result += coefficient*weight
    return S.expand(result)


def main():
    p = Path(__file__).resolve().parent
    u,v,z,t,r,U,V,K,s,alpha,ci,eta = S.symbols('u v z t r U V K s alpha ci eta')
    a,d,m = S.symbols('a d m', nonzero=True)
    b,c,k,A,B,C,D,E,lu,lv,lz,q0 = S.symbols('b c k A B C D E lu lv lz q0')
    P = z*u*v+a*u**3+d*v**3+b*u**2+c*v**2+k
    Q = (A*u*u+B*v*v+m*z*z+2*C*u*v+2*D*u*z+2*E*v*z)/2+lu*u+lv*v+lz*z+q0
    phase = t*P-Q
    grad = [S.diff(phase,x) for x in (u,v,z)]
    # The nonescaping nodal chart of the coercivity contradiction.
    scaled = [S.expand(eq.subs({u:r*U,v:r*V,z:K/r,t:1/r}, simultaneous=True))
              for eq in grad]
    assert S.expand(scaled[2]*r).subs(r,0) == -m*K
    assert S.expand(scaled[0]*r).subs(r,0) == K*(V-D)
    assert S.expand(scaled[1]*r).subs(r,0) == K*(U-E)
    # Critical values of P itself, including the nonisolated axis.
    critical_values = [k]
    for point,value in [({u:-2*b/(3*a),v:0,z:0},k+4*b**3/(27*a*a)),
                        ({u:0,v:-2*c/(3*d),z:0},k+4*c**3/(27*d*d))]:
        assert all(S.simplify(S.diff(P,x).subs(point))==0 for x in (u,v,z))
        assert S.simplify(P.subs(point)-value)==0
        critical_values.append(value)
    assert P.subs({u:0,v:0}) == k
    assert all(S.diff(P,x).subs({u:0,v:0})==0 for x in (u,v,z))
    # On a normalized phase cut, this forces every fixed-band limit into P=alpha.
    identity = (ci-alpha*t+Q-eta*s)/t
    assert S.simplify((P-alpha-identity).subs(ci,t*P-Q+eta*s))==0
    # A holomorphic top form annihilates each regular complex level surface.
    level_T = (alpha-P.subs(z,0))/(u*v)
    parametric_jac = S.Matrix([u,v,level_T]).jacobian((u,v,s))
    assert parametric_jac.det()==0
    # Area-formula volume estimate: a unit k-vector has some coordinate minor
    # of size >=1/sqrt(binomial(n,k)); this is a Cauchy-Binet identity.
    x,y,w = S.symbols('x y w', real=True)
    tangent = S.Matrix([[1,0],[0,1],[x,y],[w,0]])
    minors = [tangent[list(rows),:].det() for rows in __import__('itertools').combinations(range(4),2)]
    assert S.expand(sum(z*z for z in minors)-(tangent.T*tangent).det())==0
    # Exact majorant for the uniform unit-band tails, with Re(eta)>=1/2.
    X = S.Symbol('X')
    tail_constants = []
    for power in range(7):
        # Sum_(j>=0) (1+j)^power X^j, evaluated at X=1/2.
        rational = 1/(1-X)
        for _ in range(power):
            rational = S.factor(X*S.diff(rational,X)+rational)
        val = rational.subs(X,S.Rational(1,2))
        assert val.is_positive
        tail_constants.append(str(val))
    # Closed-half-plane rays, including the opposite-action boundary.
    ray_cases = []
    for left,right in [(0,0),(1,0),(1,-1),(1,S.I),(2,3*S.I),(-2,2)]:
        left,right = S.sympify(left),S.sympify(right)
        if not left and not right:
            ray = S.Integer(1)
        elif not left or not right:
            value = left or right
            ray = -S.conjugate(value)/S.Abs(value)
        else:
            z1,z2 = left/S.Abs(left),right/S.Abs(right)
            ray = S.I*S.conjugate(z1) if S.simplify(z1+z2)==0 else -S.conjugate(z1+z2)/S.Abs(z1+z2)
        assert S.simplify(S.Abs(ray)-1)==0
        assert S.simplify(S.re(ray*left))<=0 and S.simplify(S.re(ray*right))<=0
        ray_cases.append({'actions':[str(left),str(right)],'ray':str(S.simplify(ray))})
    triangle = [-1-S.I,2-S.I,-1+2*S.I]
    assert sum(triangle)==0
    assert S.det(S.Matrix([[3,0],[0,3]]))!=0
    # A shifted homogeneous nodal form is not an instance of the old centered
    # homogeneous theorem. The new argument covers it, including a prefix trap.
    tau = S.Symbol('tau')
    translated = (z+tau)*u*v+u**3+v**3
    moments = [wick(translated**n,(u,v,z)) for n in range(1,4)]
    assert moments == [0,tau*tau+31,54*tau]
    # Bounded critical points alone are insufficient: circular pure-null Z^2.
    Z,W,T = S.symbols('Z W T')
    control = t*Z*Z-Z*W-T*T/2
    control_gradient = [S.diff(control,q).subs({Z:1,W:2*t,T:0}) for q in (Z,W,T)]
    assert control_gradient == [0,-1,0]
    prior = json.loads((p/'nodal_gaussian_critical_verification.json').read_text())
    assert prior['completed'] and prior['finite_phase_critical_length']==7
    assert prior['center_line_boundary']['escaping_critical_points']==0
    sources = {}
    for name in ['hardt_lambrechts_turchin_volic_0806_0476v3.pdf',
                 'basu_vorobjov_definable_2006.pdf','simon_gmt_2018.pdf',
                 'pham_asterisque_1985_130_11.pdf']:
        sources[name] = hashlib.sha256((p/'sources'/name).read_bytes()).hexdigest()
    result = {'completed':True,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'prior_checker_sha256':prior['checker_sha256'],'source_sha256':sources,
              'phase_critical_length':7,'center_profile':'h=e=f=0, full affine Gaussian Q, m!=0',
              'critical_values':[str(x) for x in critical_values],
              'closed_class':'0 is not in the interior in R^2 of the convex hull of the critical values; in particular k=0, arbitrary b,c and Gaussian means',
              'ray_examples':ray_cases,'half_plane_method_boundary':{'critical_values':[str(x) for x in triangle],
                   'origin_barycentric_coordinates':['1/3']*3,'pure_null_witness':False},
              'translated_nodal_M1_M2_M3':[str(x) for x in moments],
              'translated_prefix_trap':'tau^2=-31 gives M1=M2=0 and M3=54*tau!=0',
              'bounded_critical_negative_control':'P=Z^2, Q=ZW+T^2/2 has only the origin critical, but grad(tP-Q)(1,2t,0)=(0,-1,0); uniform coercivity fails',
              'geometric_tail_majorants_at_half':tail_constants,
              'analytic_proof_certified_by_checker':False,
              'proof_boundary':'Uniform complex gradient coercivity, Hardt-compatible compactified cut chains, uniform area bounds, integral-current compactness, holomorphic-volume annihilation, the constant-period bridge and all-order exclusion remain internal written mathematics with the imported sources. Symbolic identities do not formally certify that analysis.',
              'remaining':'Finite-critical nodal profiles and center-line critical values whose triangle strictly contains0; arbitrary cuspidal/reducible-leading cubics, numerical prefixes, coefficient schemes and lower-power GVC.'}
    (p/'nodal_center_line_decay_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({key:result[key] for key in ['completed','phase_critical_length','closed_class','analytic_proof_certified_by_checker']}))


if __name__=='__main__':
    main()
