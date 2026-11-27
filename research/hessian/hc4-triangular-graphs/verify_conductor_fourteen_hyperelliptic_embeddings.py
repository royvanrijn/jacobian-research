"""Replay the coefficient certificates for the two degree-two approximate roots."""

from subprocess import run


def coefficient_certificate(inner_degree, outer_degree, chart):
    """Use the first nonzero nonleading inner coefficient as the scaling chart."""
    inner = [f"p{j}" for j in range(1, inner_degree + 1)]
    outer = [f"c{j}" for j in range(1, outer_degree - 1)]
    if chart < inner_degree:
        inner_terms = [f"T^{inner_degree}"]
        inner_terms.append(f"T^{inner_degree - chart - 1}")
        inner_terms.extend(f"{inner[j]}*T^{inner_degree-j-1}"
                           for j in range(chart + 1, inner_degree))
    else:
        inner_terms = [f"T^{inner_degree}"]
    ring_variables = ",".join(["T", "z", *inner[chart + 1:], *outer])
    outer_terms = [f"q^{outer_degree}"]
    outer_terms.extend(f"{name}*q^{j}" for j, name in enumerate(outer, 1))
    caption = f"Certified ({inner_degree},{outer_degree}) chart {chart + 1}"
    return f"""
ring r=0,({ring_variables}),dp;
proc coefficient(poly f,int j)
{{
    matrix c=coeffs(f,T);
    if (nrows(c)<j+1) {{ return(0); }}
    return(c[j+1,1]);
}}
poly q={'+'.join(inner_terms)};
poly F={'+'.join(outer_terms)};
poly g=T^7;
for (int j=1; j<=7; j++)
{{
    g=g+coefficient(F-T*g^2,15-j)*T^(7-j)/2;
}}
poly R=T*g^2-F;
ideal I=coefficient(R,2),coefficient(R,3),coefficient(R,4),coefficient(R,5),
        coefficient(R,6),coefficient(R,7),z*coefficient(R,1)-1;
matrix witness=lift(I,ideal(1));
matrix product=matrix(I)*witness;
if (product[1,1]==1) {{ "{caption}"; }}
quit;
""", caption


def main():
    for inner, outer in ((3, 5), (5, 3)):
        for chart in range(inner + 1):
            program, caption = coefficient_certificate(inner, outer, chart)
            result = run(
                ["Singular", "-q"], input=program, capture_output=True,
                text=True, check=True, timeout=30,
            )
            if result.stdout.strip() != caption or result.stderr:
                raise RuntimeError(result.stdout + "\n" + result.stderr)
            print(caption)


if __name__ == "__main__":
    main()
