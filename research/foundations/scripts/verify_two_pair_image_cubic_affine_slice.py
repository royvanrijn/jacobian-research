#!/usr/bin/env python3
"""Exact eight- and ten-moment elimination for constant-W cubic sphere slices.

Scope: p=a*w+(v0+v1*z)*t**2+(c+d*z+e*z**2)*t+sum(rj*z**j),
on 2*z*w+t**2=1. The a!=0 branch normalizes a=1,c=0. This
checker proves polynomial identities and a univariate resultant obstruction;
it does not assert that every pure cubic admits this affine chart.
The same replay covers the additional h*t**3 term as an eleven-parameter
slice: moments 1..10 force h=0, using a retained integer Bezout certificate
and a separate four-moment univariate-cubic certificate when a=0.
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


def verify_eleven_slice():
    """Check every branch of p=a*w+h*t^3+(v0+v1*z)*t^2+(c+d*z+e*z^2)*t+r(z)."""
    from math import factorial, gcd, lcm
    from sympy.polys.rings import ring
    from sympy.polys.domains import QQ
    import gzip

    z,w,t,a,h,v,c,k=S.symbols('z w t a h v c k')
    vv0,vv1,dd,ee=S.symbols('v0 v1 d e')
    rr_symbols=S.symbols('r0:4')
    original=a*w+h*t**3+(vv0+vv1*z)*t*t+(c+dd*z+ee*z*z)*t+sum(rr_symbols[j]*z**j for j in range(4))
    shifted=S.Poly(S.expand(original.subs({t:t+k*z,w:w-k*t-k*k*z/2}, simultaneous=True)),w,t,z)
    assert all(sum(powers)<=3 for powers,coef in shifted.terms())
    assert shifted.coeff_monomial(w)==a
    assert shifted.coeff_monomial(t)==c-a*k
    assert shifted.coeff_monomial(t**3)==h
    assert shifted.coeff_monomial(z*t*t)==vv1+3*h*k
    lam=S.symbols('lam',nonzero=True)
    scaled=S.Poly(S.expand(lam*original.subs({z:lam*z,w:w/lam},simultaneous=True)),w,t,z)
    assert scaled.coeff_monomial(w)==a
    assert scaled.coeff_monomial(t**3)==lam*h

    # a=0: phase averaging leaves a univariate cubic. Its first four
    # moments force h=v0=c=r0=0, with all other terms of positive weight.
    u=h*t**3+v*(t*t-S.Rational(1,3))+c*t
    def interval_average(poly):
        return S.expand(sum(coef/S.Integer(powers[0]+1) for powers,coef in S.Poly(S.expand(poly),t).terms() if powers[0]%2==0))
    M2,M3,M4=[interval_average(u**j) for j in (2,3,4)]
    assert S.expand(M2-(105*c*c+126*c*h+45*h*h+28*v*v)/315)==0
    assert S.expand(S.rem(M3,M2,v)-8*v*(56*c*c+105*c*h+45*h*h)/2205)==0
    R4=28028*c**4+84084*c**3*h+93275*c*c*h*h+46410*c*h**3+8955*h**4
    assert S.expand(S.rem(M4,M2,v)+16*R4/735735)==0
    assert S.expand(S.rem(M3.subs(h,0),M2.subs(h,0),c)+S.Rational(256,4725)*v**3)==0
    zz=S.symbols('zz')
    R=R4.subs({h:1,c:zz})
    P=56*zz*zz+105*zz+45
    B=105*zz*zz+126*zz+45
    U=-1197196*zz**3-2954952*zz**2-2305927*zz-635895
    V=2392*zz+3213
    assert S.expand(U*P+V*R)==157140
    U0=-468468*zz**3-1845844*zz**2-2151045*zz-739960
    V0=1755*zz+3756
    assert S.expand(U0*B+V0*R)==336780

    # a!=0: normalize a=1, then the null rotation makes c=0.
    # h=0 is exactly the already checked ten-parameter slice.
    # For h!=0 the torus plus scalar normalization makes h=1.
    PR,*variables=ring('h,v0,v1,d,e,r0,r1,r2,r3',QQ)
    hh,v0,v1,d,e,*rr=variables
    weights=(1,1,2,2,3,1,2,3,4)
    phase=(0,0,1,1,2,0,1,2,3)
    tpower=(3,2,2,1,1,0,0,0,0)
    @lru_cache(None)
    def df(n):
        return 1 if n<=0 else n*df(n-2)
    def exponent_vectors(n,index=0,prefix=()):
        if index==len(weights):
            if n==0:
                yield prefix
            return
        for count in range(n//weights[index]+1):
            yield from exponent_vectors(n-count*weights[index],index+1,prefix+(count,))
    def moment(n):
        result={}
        for powers in exponent_vectors(n):
            tp=sum(i*j for i,j in zip(powers,tpower))
            if tp%2:
                continue
            y_count=sum(i*j for i,j in zip(powers,phase))
            numerator=factorial(n)*df(tp-1)
            denominator=df(2*y_count+tp+1)
            for count in powers:
                denominator*=factorial(count)
            result[powers]=QQ(numerator,denominator)
        return PR.from_dict(result)
    expected=[-v0/3,
      -3*hh**2/14-2*v0**2/15-v1/5,
      -13*hh**2*v0/21-3*hh*d/7-124*v0**3/945-8*v0*v1/35,
      -29*hh**4/91-5318*hh**2*v0**2/3465-173*hh**2*v1/385-26*hh*v0*d/21-hh*e/3-2312*v0**4/14175-188*v0**2*v1/525-8*v1**2/75-d*d/2]
    solutions=[]
    moments=[]
    for n in range(1,11):
        equation=moment(n)
        for variable,value in zip(rr,solutions):
            equation=equation.compose(variable,value)
        if n<=4:
            coefficient=(QQ(1),QQ(2,3),QQ(2,5),QQ(8,35))[n-1]
            assert equation==coefficient*(rr[n-1]-expected[n-1])
            solutions.append(expected[n-1])
        else:
            moments.append(equation)

    # Check small direct expansions independently of weighted enumeration.
    normalized=original.subs({a:1,c:0})
    for n in range(1,5):
        direct=sum(coef*sphere(i,j,k) for (i,j,k),coef in S.Poly(S.expand(normalized**n),z,w,t).terms())
        assert S.expand(moment(n).as_expr()-direct)==0

    QR,vv0,vv1,ddd,eee=ring('v0,v1,d,e',QQ)
    primitive=[]
    for equation in moments:
        normalized=equation.compose(hh,PR.one)
        denominator=lcm(*(int(coef.denominator) for coef in normalized.values()))
        integers={powers:int(coef*denominator) for powers,coef in normalized.items()}
        content=gcd(*(abs(coef) for coef in integers.values()))
        projected={powers[1:5]:QQ(coef//content) for powers,coef in integers.items()}
        assert all(powers[0]==0 and not any(powers[5:]) for powers in integers)
        primitive.append(QR.from_dict(projected))
    certificate_path=Path(__file__).resolve().parents[1]/'extended-geometry'/'two_pair_image_cubic_eleven_slice_certificate.json.gz'
    with gzip.open(certificate_path,'rt',encoding='utf-8') as source:
        certificate=json.load(source)
    assert certificate['format']=='two-pair-cubic-eleven-slice-bezout-v1'
    assert certificate['variables']==['v0','v1','d','e']
    assert certificate['moment_orders']==list(range(5,11))
    assert certificate['normalization']=='a=1,c=0,h=1'
    def decode(rows):
        result={}
        for powers,coefficient in rows:
            key=tuple(powers)
            assert len(key)==4 and all(isinstance(x,int) and x>=0 for x in key)
            assert key not in result
            result[key]=QQ(int(coefficient))
        return QR.from_dict(result)
    supplied=[decode(row) for row in certificate['inputs']]
    assert supplied==primitive
    multipliers=[decode(row) for row in certificate['multipliers']]
    assert len(multipliers)==len(primitive)==6
    denominator=QR(int(certificate['denominator']))
    assert denominator!=0
    def valid(rows):
        return sum((a*b for a,b in zip(rows,primitive)),QR.zero)==denominator
    assert valid(multipliers)
    # Altering one retained coefficient must fail the exact identity.
    corrupted=list(multipliers)
    corrupted[0]=corrupted[0]+QR.one
    assert not valid(corrupted)
    return {
      'family':'p=a*y+h*t^3+(v0+v1*x)*t^2+(c+d*x+e*x^2)*t+sum(rj*x^j,j=0..3)',
      'parameter_dimension':11,
      'moment_cutoff':10,
      'claim':'Moments 1..10 vanish iff h=0 and the ten-parameter classification holds: either a=0,v0=c=r0=0, or a!=0,v0=v1=e=0 and r(x)=-x*(c+d*x)^2/(2*a).',
      'normalizations':'a!=0: scalar a=1, null rotation c=0 (h unchanged, v1 changes by 3*h*k); h!=0: torus plus scalar h=1.',
      'a_zero_branch':'Univariate uniform cubic cutoff four; h=0 quadratic branch plus h!=0 branches v0=0,v0!=0 verified by two Bezout identities.',
      'univariate_bezout_constants':[157140,336780],
      'r_eliminations':[str(p.as_expr()) for p in solutions],
      'primitive_moment_numerators':[str(p.as_expr()) for p in primitive],
      'certificate':certificate_path.name,
      'certificate_sha256':hashlib.sha256(certificate_path.read_bytes()).hexdigest(),
      'certificate_identity':'Six exact integer polynomial products sum to a nonzero integer denominator.',
      'certificate_negative_control':'Adding 1 to the first multiplier rejects the identity.',
      'full_sixteen_parameter_coverage':False,
    }

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
    eleven=verify_eleven_slice()
    artifact={
      'format':'two-pair-image-cubic-affine-slice-v1',
      'scope':'Constant W coefficient (including zero); no W*t or W*z term. Ten parameters; a!=0 branch normalizes a=1,c=0.',
      'claim':'Sphere moments 1..8 vanish iff either a=0,v0=c=r0=0 (strict positive phase), or a!=0 and p=a*w+(c+d*z)*t-z*(c+d*z)^2/(2*a).',
      'coefficient_domain':'Characteristic zero; algebraic closure for coordinate normalization if needed.',
      'sphere_relation':'2*z*w+t^2=1',
      'proof_source':'TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md',
      'moment_cutoff':8,
      'eleven_parameter_extension':eleven,
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
    print('PASS eleven-parameter extension: moments 1..10 force h=0 and the ten-parameter classification; exact retained Bezout certificate and all branches checked')
    print('This checker covers the eleven-parameter slice; the full-family proof is in the canonical note.')


if __name__=='__main__':
    main()
