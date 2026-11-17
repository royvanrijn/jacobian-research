#!/usr/bin/env python3
"""Certify dense-open closure of the two degree-42 Kuranishi branches."""

from __future__ import annotations

import shutil
import subprocess

import sympy as sp

from verify_degree42_transported_27_normal_jets import (
    serialize,
    transformed_problem,
)


def main() -> None:
    normals, bases, residuals, _defect = transformed_problem()
    e1, e2, translation, w0, w1, w2 = bases
    a, u, v = sp.symbols("sync42p_a sync42p_U sync42p_V")
    variables = normals + bases + (a, u, v)
    zero_normals = ",".join(f"{variable},0" for variable in normals)
    serialized_residuals = ",".join(serialize(item) for item in residuals)
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required for the certificate"

    factor_a = (
        -e1**2 * e2**2
        + e2**3 * translation
        + 4 * e1**3
        - 6 * e1 * e2 * translation
    )
    factor_b = (
        e1**4 * e2**4
        - 2 * e1**2 * e2**5 * translation
        + e2**6 * translation**2
        - 8 * e1**5 * e2**2
        + 20 * e1**3 * e2**3 * translation
        - 12 * e1 * e2**4 * translation**2
        + 16 * e1**6
        - 48 * e1**4 * e2 * translation
        + 27 * e1**2 * e2**2 * translation**2
        + 9 * e2**3 * translation**3
        + 36 * e1**3 * translation**2
        - 54 * e1 * e2 * translation**3
        + 27 * translation**4
    )
    discriminant_substitution = {
        e1: a**2,
        translation: 2 * a**3 - a**2 * e2,
    }
    assert sp.factor(factor_a.subs(discriminant_substitution)) == (
        a**2 * (2 * a - e2) ** 2 * (a**2 - 2 * a * e2 - e2**2)
    )
    assert sp.factor(factor_b.subs(discriminant_substitution)) == (
        a**4
        * (2 * a - e2) ** 4
        * (
            37 * a**4
            - 22 * a**3 * e2
            - 7 * a**2 * e2**2
            + 4 * a * e2**3
            + e2**4
        )
    )
    assert sp.factor(factor_a) == factor_a
    assert sp.factor(factor_b) == factor_b

    program = f"""
ring q=0,({",".join(map(str, variables))}),(dp(5),dp(9));
ideal I={serialized_residuals};

proc evd(poly p)
{{
  return(subst(p,{zero_normals},{w0},0,{e1},{a}^2,
    {translation},2*{a}^3-{a}^2*{e2}));
}}
proc Dd(poly p)
{{
  return(({e2}-{a})*diff(p,{normals[0]})+diff(p,{normals[1]}));
}}
poly dh3=evd(Dd(Dd(I[5])))/2;
poly dh4=evd(Dd(Dd(I[11])))/2;
poly dh5=evd(Dd(Dd(I[17])))/2;
poly da34=evd(diff(I[5],{normals[3]}));
poly da35=evd(diff(I[5],{normals[4]}));
poly da45=evd(diff(I[11],{normals[4]}));
poly dz5=dh5;
poly dz4=dh4+da45*dz5;
poly dz3=dh3+da34*dz4+da35*dz5;
poly dc18=evd(Dd(Dd(Dd(I[18]))))/6
  +dz3*evd(Dd(diff(I[18],{normals[2]})))
  +dz4*evd(Dd(diff(I[18],{normals[3]})))
  +dz5*evd(Dd(diff(I[18],{normals[4]})));
poly dc19=evd(Dd(Dd(Dd(I[19]))))/6
  +dz3*evd(Dd(diff(I[19],{normals[2]})))
  +dz4*evd(Dd(diff(I[19],{normals[3]})))
  +dz5*evd(Dd(diff(I[19],{normals[4]})));
poly expected19=-5/8*{w2}*{a}^2*({e2}-2*{a})^2;
int discriminant_arc=(
  dc19==expected19 and dc18==(7*{e2}-{a})*expected19
);

proc evw(poly p)
{{
  return(subst(p,{zero_normals},{w0},0,{w1},0));
}}
proc Dw(poly p)
{{
  return({u}*diff(p,{normals[0]})+{v}*diff(p,{normals[1]}));
}}
poly wh3=evw(Dw(Dw(I[5])))/2;
poly wh4=evw(Dw(Dw(I[11])))/2;
poly wh5=evw(Dw(Dw(I[17])))/2;
poly wa34=evw(diff(I[5],{normals[3]}));
poly wa35=evw(diff(I[5],{normals[4]}));
poly wa45=evw(diff(I[11],{normals[4]}));
poly wz5=wh5;
poly wz4=wh4+wa45*wz5;
poly wz3=wh3+wa34*wz4+wa35*wz5;
poly wc18=evw(Dw(Dw(Dw(I[18]))))/6
  +wz3*evw(Dw(diff(I[18],{normals[2]})))
  +wz4*evw(Dw(diff(I[18],{normals[3]})))
  +wz5*evw(Dw(diff(I[18],{normals[4]})));
poly wc19=evw(Dw(Dw(Dw(I[19]))))/6
  +wz3*evw(Dw(diff(I[19],{normals[2]})))
  +wz4*evw(Dw(diff(I[19],{normals[3]})))
  +wz5*evw(Dw(diff(I[19],{normals[4]})));
poly cubic_resultant=resultant(
  subst(wc18,{v},1)/{w2},subst(wc19,{v},1)/{w2},{u}
);
poly expected_resultant=-15625/262144
  *({serialize(factor_a)})*({serialize(factor_b)});
int cubic_pair=(cubic_resultant==expected_resultant);

poly factorA={serialize(factor_a)};
poly factorB={serialize(factor_b)};
poly commonA={e2}*{u}+(2*{e1}-{e2}^2)*{v};
poly alphaB=
  8*{e1}^5-2*{e1}^4*{e2}^2-20*{e1}^3*{e2}*{translation}
  +4*{e1}^2*{e2}^3*{translation}+6*{e1}^2*{translation}^2
  +12*{e1}*{e2}^2*{translation}^2
  -2*{e2}^4*{translation}^2-9*{e2}*{translation}^3;
poly betaB=
  -4*{e1}^5*{e2}+{e1}^4*{e2}^3+4*{e1}^4*{translation}
  +9*{e1}^3*{e2}^2*{translation}
  -2*{e1}^2*{e2}^4*{translation}
  -9*{e1}^2*{e2}*{translation}^2
  -5*{e1}*{e2}^3*{translation}^2
  +{e2}^5*{translation}^2+6*{e2}^2*{translation}^3;
poly commonB=alphaB*{u}+betaB*{v};
poly denominatorA={e2}*(6*{e1}-{e2}^2);
ideal commonAI=std(ideal(factorA,commonA,1-{a}*denominatorA));
ideal commonBI=std(ideal(factorB,commonB,1-{a}*alphaB));
int cubic_common_factors=(
  reduce(wc18/{w2},commonAI)==0 and
  reduce(wc19/{w2},commonAI)==0 and
  reduce(wc18/{w2},commonBI)==0 and
  reduce(wc19/{w2},commonBI)==0
);

print("DEGREE42_KURANISHI_BRANCHES");
print(discriminant_arc);
print(cubic_pair);
print(cubic_common_factors);
"""
    process = subprocess.run(
        [singular, "-q"],
        input=program,
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    compact = process.stdout.split()
    marker = compact.index("DEGREE42_KURANISHI_BRANCHES")
    assert compact[marker + 1 : marker + 4] == ["1", "1", "1"], (
        process.stdout + process.stderr
    )
    print(
        "PASS: the discriminant branch has a nonzero common-tangent "
        "cubic on D(w1*w2*t), and the w1=0 branch has coprime cubic "
        "Kuranishi forms on D(w2*A*B); on A=0 and B=0 their explicit "
        "common linear factors are certified"
    )


if __name__ == "__main__":
    main()
