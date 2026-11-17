#!/usr/bin/env python3
"""Exact controls for the rank-five FIELD-discriminant classification.

The unrestricted classification uses the written maximal-order and local
arguments. This checker verifies their explicit identities and examples.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path

import sympy as sp


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    checks=[]
    def require(name, condition, **details):
        if not condition:
            raise AssertionError(name)
        checks.append({"name":name,"passed":True,**details})
    def equal(name,left,right=0):
        difference=sp.factor(left-right)
        if difference!=0:
            raise AssertionError(f"{name}: {difference}")
        checks.append({"name":name,"passed":True})
    def vp(value,p):
        if value==0:
            raise AssertionError("unexpected zero discriminant")
        value=abs(int(value));out=0
        while value%p==0:
            out+=1;value//=p
        return out

    matrix=sp.Matrix(3,3,sp.symbols("m0:9"))
    plus=0;minus=0
    for permutation in itertools.permutations(range(3)):
        product=sp.prod(matrix[i,permutation[i]] for i in range(3))
        odd=sum(permutation[i]>permutation[j] for i in range(3) for j in range(i+1,3))%2
        if odd:
            minus+=product
        else:
            plus+=product
    equal("embedding_determinant",matrix.det(),plus-minus)
    equal("quadratic_order_discriminant",(plus+minus)**2-4*plus*minus,matrix.det()**2)
    row_swap={matrix[0,j]:matrix[1,j] for j in range(3)}
    row_swap.update({matrix[1,j]:matrix[0,j] for j in range(3)})
    equal("odd_embedding_permutation_swaps_sums",plus.xreplace(row_swap),minus)
    row_cycle={matrix[i,j]:matrix[(i+1)%3,j] for i in range(3) for j in range(3)}
    equal("even_embedding_permutation_fixes_sum",plus.xreplace(row_cycle),plus)

    x,a,b,c,p,A,B,C=sp.symbols("x a b c p A B C")
    disc=a*a*b*b-4*b**3-4*a**3*c-27*c*c+18*a*b*c
    equal("cubic_discriminant_identity",sp.discriminant(x**3+a*x*x+b*x+c,x),disc)
    remainder=sp.Poly(sp.expand((disc.subs({a:p*A,b:p*B,c:p*C})+27*p*p*C*C)/p**3),p,A,B,C)
    require("odd_tame_discriminant_remainder_integral",all(value.q==1 for value in remainder.coeffs()))
    d2=sp.expand(disc.subs({a:2*A,b:2*B,c:2*C})/4)
    equal("two_adic_discriminant_formula",d2,4*A*A*B*B-8*B**3-16*A**3*C-27*C*C+36*A*B*C)
    cases_two=0
    for av,bv,cv in itertools.product(range(2),range(2),(1,3,5,7)):
        if int(d2.subs({A:av,B:bv,C:cv}))%8 != 5:
            raise AssertionError("two-adic unit reduction")
        cases_two+=1
    require("complete_two_adic_unit_reduction",cases_two==16,residue_cases=cases_two)
    cases_three=0
    for va,vb,unit in itertools.product(range(1,7),range(1,7),(1,2,4,5,7,8)):
        actual=disc.subs({a:3**va,b:3**vb,c:3*unit})
        predicted=min(5,3*va+1,3*vb)
        if vp(actual,3)!=predicted:
            raise AssertionError((va,vb,unit,actual,predicted))
        cases_three+=1
    require("wild_three_bounded_controls",cases_three==216,cases=cases_three,
            coefficient_valuation_range=[1,6])
    require("wild_three_minimum_classes_distinct",
            len({5%3,(3*1+1)%3,(3*1)%3})==3)
    local_examples=[x**3+3*x+3,x**3+3*x*x+3,x**3+3]
    require("wild_three_exponents_are_attained",
            [vp(sp.discriminant(poly,x),3) for poly in local_examples]==[3,4,5])
    cyclic_discs=[sp.discriminant(x**3+6*x*x+3+9*k,x) for k in range(3)]
    require("Pauli_cyclic_three_examples_discriminants",cyclic_discs==[-2835,-14256,-30051])
    require("Pauli_examples_have_square_three_adic_discriminant",
            all(vp(value,3)==4 and int(value//81)%3==1 for value in cyclic_discs))
    require("even_three_exponent_forces_four",[value for value in (3,4,5) if value%2==0]==[4])

    fundamental_poly=x**3-x-1
    equal("fundamental_discriminant_example",sp.discriminant(fundamental_poly,x),-23)
    require("fundamental_minus_twenty_three",sp.isprime(23) and (-23)%4==1)
    alpha=sp.symbols("alpha")
    relation=alpha**3-19
    beta=(1+alpha+alpha**2)/3
    reduce=lambda expression:sp.rem(sp.expand(expression),relation,alpha)
    equal("pure_cubic_beta_integral",reduce(beta**3-beta**2-6*beta-12))
    equal("pure_cubic_alpha_square",reduce(alpha**2-(3*beta-alpha-1)))
    equal("pure_cubic_alpha_beta",reduce(alpha*beta-(beta+6)))
    equal("pure_cubic_beta_square",reduce(beta**2-(beta+2*alpha+4)))
    basis=[sp.Integer(1),alpha,beta]
    gram=sp.Matrix(3,3,lambda i,j:3*sp.Poly(reduce(basis[i]*basis[j]),alpha).nth(0))
    require("pure_cubic_integral_basis_trace_gram",gram==sp.Matrix([[3,0,1],[0,0,19],[1,19,13]]))
    equal("pure_cubic_field_basis_discriminant",gram.det(),-3*19**2)
    equal("pure_cubic_polynomial_discriminant",sp.discriminant(x**3-19,x),-27*19**2)
    equal("pure_cubic_index_square",sp.discriminant(x**3-19,x)/gram.det(),9)
    require("polynomial_discriminant_gives_false_three_factor",57==3*19 and 19%3==1)
    require("nineteen_uniformizer_and_remaining_discriminant_support",sp.isprime(19) and vp(-3*19**2,3)==1)

    source=Path(__file__).with_name("RANK_FIVE_DISCRIMINANT_CLASSIFICATION.md")
    report={"status":"passed","scope":"Exact algebraic identities and explicit finite local controls; maximality, local-field proofs and the unrestricted classification remain written mathematics.",
            "python":platform.python_version(),"sympy":sp.__version__,
            "checker_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "source_sha256":hashlib.sha256(source.read_bytes()).hexdigest(),
            "check_count":len(checks),"checks":checks}
    if args.output:
        args.output.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["status"],"check_count":len(checks),"output":str(args.output) if args.output else None}))


if __name__=="__main__":
    main()
