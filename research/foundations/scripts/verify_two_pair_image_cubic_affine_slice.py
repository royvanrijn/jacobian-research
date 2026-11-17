#!/usr/bin/env python3
"""Exact eight-moment elimination for the constant-W cubic sphere chart.

Scope: p=a*w+(v0+v1*z)*t**2+(c+d*z+e*z**2)*t+sum(rj*z**j),
on 2*z*w+t**2=1. The a!=0 branch normalizes a=1,c=0. This
checker proves polynomial identities and a univariate resultant obstruction;
it does not assert that every pure cubic admits this affine chart.
Requires the repository SymPy environment (tested with SymPy 1.14.0).
"""
from functools import lru_cache
from pathlib import Path
import hashlib
import json
import sympy as S


@lru_cache(None)
def sphere(a, b, c):
    if a != b or c % 2:
        return S.Integer(0)
    return S.factorial(a) * (S.factorial2(c-1) if c else 1) / S.factorial2(2*a+c+1)


def sparse_moments(base, cutoff):
    power = {(0,0,0): S.Integer(1)}
    rows = []
    for unused in range(cutoff):
        product = {}
        for left, lv in power.items():
            for right, rv in base.items():
                key = tuple(a+b for a,b in zip(left,right))
                product[key] = product.get(key,0) + lv*rv
        power = {key:S.expand(value) for key,value in product.items()}
        rows.append(S.expand(sum(value*sphere(*key) for key,value in power.items())))
    return rows


def main():
    assert S.__version__ == "1.14.0", "Use the recorded SymPy 1.14.0 environment"
    z,w,t,a,c,k = S.symbols('z w t a c k')
    v0,v1,d,e = S.symbols('v0 v1 d e')
    rr = S.symbols('r0:4')
    q=2*z*w+t*t
    assert S.expand(q.subs({t:t+k*z,w:w-k*t-k*k*z/2}, simultaneous=True)-q)==0
    original=a*w+(v0+v1*z)*t*t+(c+d*z+e*z*z)*t+sum(rr[j]*z**j for j in range(4))
    shifted=S.Poly(S.expand(original.subs({t:t+k*z,w:w-k*t-k*k*z/2}, simultaneous=True)),w,t,z)
    assert all(sum(powers)<=3 for powers,coef in shifted.terms())
    assert S.expand(shifted.as_expr()).coeff(t,1).coeff(z,0).coeff(w,0)==c-a*k
    assert S.expand(shifted.as_expr()).coeff(w,1)==a
    # Monomial integrals are checked directly against the uniform t/phase rule.
    for exponent in range(9):
        for even in range(0,9,2):
            integral=sum((-1)**j*S.binomial(exponent,j)/S.Integer(even+2*j+1) for j in range(exponent+1))/2**exponent
            assert S.simplify(integral-sphere(exponent,exponent,even))==0
    # If a=0, only phase zero contributes; its first three moments
    # force v0=c=r0=0, leaving strict positive z phase weight.
    zero_phase=v0*t*t+c*t+rr[0]
    zero_moments=[]
    for order in range(1,4):
        zero_moments.append(S.expand(sum(coef/(powers[0]+1) if powers[0]%2==0 else 0 for powers,coef in S.Poly(S.expand(zero_phase**order),t).terms())))
    zero_reduced=[S.expand(row.subs(rr[0],-v0/3)) for row in zero_moments]
    assert zero_reduced[0]==0
    assert S.expand(zero_reduced[1]-(4*v0*v0/45+c*c/3))==0
    assert S.expand(zero_reduced[2].subs(c*c,-4*v0*v0/15)+S.Rational(256,4725)*v0**3)==0
    base={(0,1,0):S.Integer(1),(0,0,2):v0,(1,0,2):v1,(1,0,1):d,(2,0,1):e}
    base.update({(j,0,0):rr[j] for j in range(4)})
    moments=sparse_moments(base,8)
    solution={}
    expected=[-v0/3,
      -(2*v0**2+3*v1)/15,
      -4*v0*(31*v0**2+54*v1)/945,
      -(14175*d*d+4624*v0**4+10152*v0*v0*v1+3024*v1*v1)/28350]
    for n in range(4):
        equation=S.expand(moments[n].subs(solution))
        value=S.factor(S.solve(equation,rr[n])[0])
        assert S.expand(value-expected[n])==0
        assert S.Poly(equation,rr[n]).degree()==1
        assert equation.coeff(rr[n])==[S.Integer(1),S.Rational(2,3),S.Rational(2,5),S.Rational(8,35)][n]
        assert not (set(rr[n+1:]) & equation.free_symbols)
        solution[rr[n]]=value
    reduced=[S.factor(row.subs(solution)) for row in moments[4:]]
    N5=155925*d*d*v0+155925*d*e+35504*v0**5+94992*v0**3*v1+52704*v0*v1*v1
    N6=77395500*d*d*v0*v0-116093250*d*d*v1+77395500*d*e*v0-58046625*e*e+30907712*v0**6+63343584*v0**4*v1+3500496*v0*v0*v1*v1-12573792*v1**3
    N7=226800*d*d*v0**3-297675*d*d*v0*v1+226800*d*e*v0*v0+127575*d*e*v1-212625*e*e*v0+53504*v0**7+17328*v0**5*v1-211968*v0**3*v1*v1-108864*v0*v1**3
    N8=4819500*d*d*v0**4-216877500*d*d*v0*v0*v1+4819500*d*e*v0**3-65063250*d*e*v0*v1-75907125*e*e*v0*v0+7986496*v0**8-27544704*v0**6*v1-90351504*v0**4*v1*v1-16632864*v0*v0*v1**3+9144576*v1**4
    numerators=[N5,N6,N7,N8]
    scales=[S.Rational(8,9823275),S.Rational(-8,13408770375),S.Rational(32,273648375),S.Rational(-256,1465387048125)]
    for value,num,scale in zip(reduced,numerators,scales):
        assert S.expand(value-scale*num)==0
    # The v0=0 branch gives v1=0 from M8, then e=0 from M6.
    assert S.expand(N8.subs(v0,0))==9144576*v1**4
    assert S.expand(N6.subs({v0:0,v1:0}))==-58046625*e**2
    assert [S.simplify(solution[x].subs({v0:0,v1:0,e:0})) for x in rr]==[0,0,0,-d*d/2]
    # For v0 != 0, weighted normalization gives four linear equations in D,E,F.
    b,D,E,F=S.symbols('b D E F')
    normalized=[]
    for n,num in enumerate(numerators,5):
        expression=0
        for powers,coef in S.Poly(num,v0,v1,d,e).terms():
            i,j,h,l=powers
            assert i+2*j+2*h+3*l==n
            replacements={(0,0):1,(2,0):D,(1,1):E,(0,2):F}
            expression += coef*b**j*replacements[(h,l)]
        normalized.append(S.expand(expression))
    matrix=S.Matrix([[row.coeff(var) for var in (D,E,F)]+[row.subs({D:0,E:0,F:0})] for row in normalized])
    assert matrix[:,0]-matrix[:,1]-2*b*matrix[:,2]==S.zeros(4,1)
    f=2485161*b**3+6647562*b*b+4908168*b+1097984
    g=10216206*b**4+24356619*b**3-25372278*b*b-41158800*b-11711872
    minor1=S.Poly(matrix[[0,1,2],[0,2,3]].det(),b).primitive()[1].as_expr()
    minor2=S.Poly(matrix[[0,1,3],[0,2,3]].det(),b).primitive()[1].as_expr()
    assert S.expand(minor1-f)==0
    assert S.expand(minor2+g)==0
    resultant=S.resultant(f,g,b)
    assert resultant==-3179537325694882741702679455938048000000000
    u,v,unit=S.gcdex(f,g,b)
    assert unit==1 and S.expand(u*f+v*g)==1
    # Relevant negative control: first four vanishing moments alone do not suffice.
    false_values={v0:0,v1:1,d:0,e:0}
    assert S.simplify(reduced[3].subs(false_values))!=0
    artifact={
      'format':'two-pair-image-cubic-affine-slice-v1',
      'scope':'Constant W coefficient (including zero); no W*t or W*z term. Ten parameters; a!=0 branch normalizes a=1,c=0.',
      'claim':'Sphere moments 1..8 vanish iff either a=0,v0=c=r0=0 (strict positive phase), or a!=0 and p=a*w+(c+d*z)*t-z*(c+d*z)^2/(2*a).',
      'coefficient_domain':'Characteristic zero; algebraic closure for coordinate normalization if needed.',
      'sphere_relation':'2*z*w+t^2=1',
      'proof_source':'TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md',
      'moment_cutoff':8,
      'sympy_version':S.__version__,
      'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'r_eliminations':[str(value) for value in expected],
      'zero_a_branch':'First three moments force v0=c=r0=0.',
      'zero_a_reduced_second_moment':str(zero_reduced[1]),
      'zero_a_reduced_third_moment':str(zero_reduced[2]),
      'reduced_moment_factors':[str(value) for value in reduced],
      'resultant':str(resultant),
      'bezout_coefficients':[str(u),str(v)],
      'full_sixteen_parameter_coverage':False,
      'independent_replay':False,
      'formal_verification':False,
      'external_review':False,
    }
    output=Path(__file__).resolve().parents[1]/'extended-geometry'/'two_pair_image_cubic_affine_slice_verification.json'
    output.write_text(json.dumps(artifact,indent=2)+'\n')
    print('PASS exact constant-W ten-parameter cubic slice: moments 1..8 force one-sided support or the shear family')
    print('PASS all normalization identities, four r-eliminations, degenerate branch, coprime-minor certificate, negative control')
    print('Full sixteen-parameter coverage remains unproved.')


if __name__=='__main__':
    main()
