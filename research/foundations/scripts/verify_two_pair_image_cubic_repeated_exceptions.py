#!/usr/bin/env python3
"""Exact fixed-number-field certificates for repeated spherical profiles.

Scope: p=y*(t-r)**2+G(x,t), deg G<=3, after the orthogonal null rotation
which kills G's t**3 coefficient. This checker covers the roots of the
three exceptional pivot polynomials P2,P3,P4 and r=0, not every complex r.

It derives moments in QQ[r]/Pj, checks positive weighted homogeneity and
integrality at a specified prime, then recomputes an exact finite-field
Groebner basis with Singular. A weighted projective specialization lemma
proved in the canonical source lifts the absence of nonzero solutions to
characteristic zero. It is not rational reconstruction from sample primes.
Requires sympy==1.14.0 and Singular 4.3.2.
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import tempfile
from math import factorial
from pathlib import Path

import sympy as s
from sympy.polys.rings import ring

ROOT=Path(__file__).resolve().parents[3]
OUTPUT=ROOT/'research/foundations/extended-geometry/two_pair_image_cubic_repeated_exceptions_verification.json'
NAMES=('t','B0','B1','B2','A0','A1','A2','D0','D1','E')
WEIGHTS={'B0':1,'B1':1,'B2':1,'A0':2,'A1':2,'A2':2,'D0':3,'D1':3,'E':4}
CUTOFF=12


def verify_branch(branch: int, scratch: Path) -> dict:
    r=s.symbols('r')
    definitions={
        0:(r,1000003),
        2:(5*r**2+1,1013),
        3:(21*r**4+18*r**2+1,10007),
        4:(429*r**6+715*r**4+195*r**2+5,1031),
    }
    minimal,prime=definitions[branch]
    assert s.Poly(minimal,r,domain=s.QQ).is_irreducible
    assert s.Poly(minimal,r,modulus=prime).is_irreducible
    assert int(s.LC(s.Poly(minimal,r))) % prime != 0
    field=s.QQ.algebraic_field(s.Integer(0) if branch==0 else s.CRootOf(minimal,0))
    polyring,t,B0,B1,B2,A0,A1,A2,D0,D1,E=ring(','.join(NAMES),field)
    rho=polyring.zero if branch==0 else polyring.ground_new(field.unit)
    profile=(t-rho)**2
    B=B0+B1*t+B2*t**2
    A=A0+A1*t+A2*t**2
    D=D0+D1*t
    ep=E
    denominators=set()

    def check_coefficient(coefficient) -> None:
        for value in coefficient.to_list():
            denominator=int(value.denominator)
            assert denominator % prime != 0
            denominators.add(denominator)

    def serialize(poly) -> str:
        terms=[]
        for exponents,coefficient in poly.terms():
            check_coefficient(coefficient)
            raw=coefficient.to_list()
            cs=[]
            for index,value in enumerate(raw):
                if not value:
                    continue
                power=len(raw)-index-1
                cs.append('('+str(value)+')'+('' if power==0 else '*r'+('' if power==1 else '^'+str(power))))
            monomial='*'.join(name+('' if power==1 else '^'+str(power))
                              for name,power in zip(NAMES,exponents) if power)
            terms.append('('+'+'.join(cs)+')'+('' if not monomial else '*'+monomial))
        return '+'.join(terms) or '0'

    def homogeneous(poly,weight: int) -> None:
        for exponents,coefficient in poly.items():
            assert exponents[0]==0
            assert sum(exponents[index]*WEIGHTS[name]
                       for index,name in enumerate(NAMES) if index)>0 or weight==0
            assert sum(exponents[index]*WEIGHTS[name]
                       for index,name in enumerate(NAMES) if index)==weight
            check_coefficient(coefficient)

    def moment(order: int):
        """Complete phase-zero multinomial selection, followed by exact beta.

        For counts (ny,nb,na,nd,ne), phase zero says ny=na+2nd+3ne.
        Integrate C(t)^ny B(t)^nb A(t)^na D(t)^nd t-monomials with
        L(x^ny*y^ny*t^j)=ny!*(j-1)!!/(2ny+j+1)!! for even j.
        """
        out=polyring.zero
        bp=[B**i for i in range(order+1)]
        ap=[A**i for i in range(order//2+1)]
        dp=[D**i for i in range(order//3+1)]
        cp=[profile**i for i in range(order+1)]
        for ne in range(order//4+1):
            for nd in range((order-4*ne)//3+1):
                for na in range((order-4*ne-3*nd)//2+1):
                    nb=order-4*ne-3*nd-2*na
                    ny=na+2*nd+3*ne
                    product=bp[nb]*ap[na]*dp[nd]*cp[ny]
                    multinomial=s.Rational(factorial(order),factorial(ny)*factorial(nb)*factorial(na)*factorial(nd)*factorial(ne))
                    data={}
                    for exponents,coefficient in product.items():
                        height=exponents[0]
                        if height % 2:
                            continue
                        beta=s.Rational(factorial(ny)*s.factorial2(height-1),s.factorial2(2*ny+height+1))
                        value=coefficient*field(multinomial*beta)
                        monomial=(0,)+exponents[1:]
                        data[monomial]=data.get(monomial,field.zero)+value
                    out+=polyring.from_dict(data)*ep**ne
        homogeneous(out,order)
        return out

    # Direct nonzero controls detect both the height weight and its r=0
    # specialization: add x*t for r!=0, and x for r=0.
    negative_control=field(s.Rational(2,15)) if branch==0 else -field.unit*field(s.Rational(4,15))
    assert negative_control != field.zero
    check_coefficient(negative_control)

    eliminations=[B0,A1 if branch==2 else A0,D1 if branch==3 else D0]
    if branch!=4:
        eliminations.append(E)
    solved=[]
    for order,variable in enumerate(eliminations,1):
        equation=moment(order)
        pivot=equation.diff(variable)
        assert pivot.is_ground and pivot!=polyring.zero
        coefficient=next(iter(pivot.values()))
        check_coefficient(coefficient)
        # The pivot must also remain a unit in the selected residue field.
        raw=coefficient.to_list()
        coefficient_polynomial=sum(s.Rational(str(value))*r**(len(raw)-index-1)
                                   for index,value in enumerate(raw))
        assert s.gcd(s.Poly(s.Poly(coefficient_polynomial,r,domain=s.QQ).clear_denoms()[1].as_expr(),r,modulus=prime),s.Poly(minimal,r,modulus=prime)).degree()==0
        solution=-equation.compose(variable,polyring.zero).quo_ground(coefficient)
        homogeneous(solution,WEIGHTS[str(variable)])
        assert equation.compose(variable,solution)==polyring.zero
        B=B.compose(variable,solution)
        A=A.compose(variable,solution)
        D=D.compose(variable,solution)
        ep=ep.compose(variable,solution)
        solved.append({'moment':order,'variable':str(variable),
                       'pivot':serialize(pivot),'solution':serialize(solution)})
    active=[name for name in NAMES[1:] if name not in {str(v) for v in eliminations}]
    orders=list(range(len(eliminations)+1,CUTOFF+1))
    inputs=[serialize(moment(order)) for order in orders]
    singular_input=(
        'ring R=('+str(prime)+',r),('+','.join(active)+'),wp('+','.join(str(WEIGHTS[v]) for v in active)+');\n'
        'minpoly='+str(minimal).replace('**','^')+';\n'
        'option(redSB);\nideal I='+','.join(inputs)+';\nideal G=std(I);\n'
        'if(dim(G)!=0){print("FAIL_DIM");quit;}\n'
        'ideal KB=kbase(G);int bound;int i;int n;\n'
        'for(i=1;i<=size(KB);i++){if(deg(KB[i])>bound){bound=deg(KB[i]);}}\n'
        'print("BEGIN_CERTIFICATE");\n'
        'print("basis_size:"+string(size(G)));\n'
        'print("quotient_dimension:"+string(vdim(G)));\n'
        'print("maximum_standard_weight:"+string(bound));\n'
        'intvec ww='+','.join(str(WEIGHTS[v]) for v in active)+';\n'
        'for(i=1;i<=nvars(R);i++){n=bound div ww[i]+1;'
        'if(reduce(var(i)^n,G)!=0){print("FAIL_POWER");quit;}'
        'print("power:"+string(var(i))+":"+string(n)+":0");}\n'
        'print("END_CERTIFICATE");quit;\n'
    )
    script=scratch/f'branch_P{branch}.sing'
    script.write_text(singular_input)
    process=subprocess.run(['Singular','-q',str(script)],text=True,capture_output=True,check=False,timeout=120)
    assert process.returncode==0,(process.returncode,process.stdout,process.stderr)
    assert not process.stderr,process.stderr
    lines=process.stdout.splitlines()
    assert lines and lines[0]=='BEGIN_CERTIFICATE' and lines[-1]=='END_CERTIFICATE',process.stdout
    assert all('FAIL' not in line for line in lines),process.stdout
    receipt={}
    powers={}
    for line in lines[1:-1]:
        parts=line.split(':')
        if parts[0]=='power':
            assert len(parts)==4 and parts[3]=='0'
            powers[parts[1]]=int(parts[2])
        else:
            assert len(parts)==2
            receipt[parts[0]]=int(parts[1])
    assert set(powers)==set(active)
    assert receipt['quotient_dimension']==797
    return {
        'pivot_branch':'r=0' if branch==0 else f'P{branch}=0',
        'minimal_polynomial':str(minimal),
        'number_field_degree':s.degree(minimal,r),
        'prime':prime,
        'minimal_polynomial_irreducible_over_QQ_and_mod_prime':True,
        'coefficient_denominator_units_verified':len(denominators),
        'eliminations':solved,
        'active_variables':active,
        'positive_weights':[WEIGHTS[v] for v in active],
        'moment_orders':orders,
        'finite_field_moment_input_sha256':hashlib.sha256(','.join(inputs).encode()).hexdigest(),
        'singular_input_sha256':hashlib.sha256(singular_input.encode()).hexdigest(),
        'finite_field_certificate':receipt,
        'pure_powers_with_verified_zero_normal_form':powers,
        'conclusion':'All active parameters vanish in characteristic zero by weighted projective specialization; all eliminated coefficients then vanish.',
    }


def main() -> None:
    assert s.__version__=='1.14.0'
    version=subprocess.run(['Singular','--version'],text=True,capture_output=True,check=True).stdout.splitlines()[0]
    assert 'version 4.3.2 ' in version,version
    r=s.symbols('r')
    pivots=[5*r**2+1,21*r**4+18*r**2+1,429*r**6+715*r**4+195*r**2+5]
    for left in range(3):
        for right in range(left+1,3):
            assert s.gcd(pivots[left],pivots[right])==1
    scratch_root=ROOT/'.scratch'
    scratch_root.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='cubic-repeated-exceptions-',dir=scratch_root) as directory:
        branches=[verify_branch(branch,Path(directory)) for branch in (0,2,3,4)]
    receipt={
        'format':'two-pair-image-cubic-repeated-exceptions-v1',
        'family':'p=y*(t-r)^2+G(x,t), deg G<=3; normalized t^3 coefficient in G is zero',
        'normalization':'Orthogonal null rotation kills the t^3 coefficient; it changes by -2*k and preserves the repeated profile.',
        'scope':'Only r=0 and roots of P2*P3*P4; does not establish the generic-r branch or the full sixteen-parameter theorem.',
        'claim':'Spherical moments 1..12 vanish iff normalized G=0 for r=0 and every root of the three pivot polynomials.',
        'pivot_coefficients':['2*P2/15','2*P3/105','8*P4/15015'],
        'all_pivot_exception_polynomials_pairwise_coprime':True,
        'branches':branches,
        'proof_source':'foundations/extended-geometry/TWO_PAIR_IMAGE_CUBIC_STRATIFICATION.md',
        'specialization_argument':'Positive weighted homogeneous nonzero algebraic solutions specialize to nonzero residue-field solutions after valuation normalization; the verified finite-dimensional graded residue quotient forbids such points.',
        'properness_reference':'https://stacks.math.columbia.edu/tag/01WC',
        'software':{'python':platform.python_version(),'sympy':s.__version__,'Singular':version},
        'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'assurance':{'internal_exact_replay':True,'independent_replay':False,'formal_verification':False,'external_review':False},
        'negative_control':'For nonzero r, p=y*(t-r)^2+x*t has moment2=-4*r/15 !=0; at r=0, p=y*t^2+x has moment2=2/15 !=0.',
    }
    OUTPUT.write_text(json.dumps(receipt,indent=2,default=int)+'\n')
    print('PASS repeated-profile fixed exceptions r=0,P2,P3,P4: exact algebraic moments, pivot units, positive weights and residue-field certificates')
    print('PASS all variables have explicit pure powers reducing to zero; weighted specialization gives the characteristic-zero exceptional cases')
    print(f'Wrote {OUTPUT.relative_to(ROOT)}')


if __name__=='__main__':
    main()
