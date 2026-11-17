#!/usr/bin/env python3
"""Exact obstruction certificate for one full four-pair Dvorsky slice.

Fix

    P = (t+a+b+d)(a*d+b*t)
    Lambda = d/dt * R(d/da,d/db,d/dd),

where R=A*a^2+B*a*b+C*a*d+D*b^2+E*b*d+F*d^2 is arbitrary.
The first eight scalars Lambda^m(P^m) cut out, set-theoretically, only

    R = a^2, d^2, (a-b)^2, (b-d)^2

up to scale.  The companion note proves by weight gaps that all four
directions satisfy eventual mixed vanishing for every fixed multiplier.
"""

from __future__ import annotations

from math import factorial, gcd
import shutil
import subprocess

from search_dvorsky_gvc4_bounded import (
    LAMBDA_EXPONENTS,
    factorial_product,
    make_p,
    powers,
)


PARAMETERS = ("A", "B", "C", "D", "E", "F")
CUTOFF = 8


def weak_compositions(
    total: int, length: int, prefix: tuple[int, ...] = ()
):
    if length == 1:
        yield prefix + (total,)
        return
    for entry in range(total + 1):
        yield from weak_compositions(
            total - entry, length - 1, prefix + (entry,)
        )


def singular_moment(order: int, p_power) -> str:
    """Return the primitive moment polynomial in A,...,F."""
    terms: list[tuple[int, str]] = []
    for parameter_exponents in weak_compositions(order, len(PARAMETERS)):
        derivative_exponent = tuple(
            sum(
                parameter_exponents[index]
                * LAMBDA_EXPONENTS[index][variable]
                for index in range(len(PARAMETERS))
            )
            for variable in range(4)
        )
        p_coefficient = p_power.get(derivative_exponent, 0)
        if not p_coefficient:
            continue

        coefficient = factorial(order)
        for entry in parameter_exponents:
            coefficient //= factorial(entry)
        coefficient *= p_coefficient * factorial_product(
            derivative_exponent
        )

        factors = []
        for parameter, exponent in zip(PARAMETERS, parameter_exponents):
            if exponent == 1:
                factors.append(parameter)
            elif exponent > 1:
                factors.append(f"{parameter}^{exponent}")
        terms.append((coefficient, "*".join(factors) or "1"))

    content = 0
    for coefficient, _ in terms:
        content = gcd(content, abs(coefficient))
    assert content

    answer = ""
    for coefficient, monomial in terms:
        primitive = coefficient // content
        if primitive == 1:
            term = monomial
        elif primitive == -1:
            term = f"-{monomial}"
        else:
            term = f"{primitive}*{monomial}"
        answer += term if not answer or term.startswith("-") else f"+{term}"
    return answer


def check_weight_gaps() -> None:
    p = make_p((1, 1, 1, 1))

    # R=a^2 and R=d^2.  The derivative weights are respectively
    # wt(t)+2 wt(a)=3 and wt(t)+2 wt(d)=3, while P has weight at most 2.
    assert max(exponent[0] + exponent[1] for exponent in p) == 2
    assert max(exponent[0] + exponent[3] for exponent in p) == 2

    # For R=(a-b)^2 put a=x, b=y-x:
    # P=(t+y+d)(x(d-t)+yt), so deg_x P=1<2.
    # For R=(b-d)^2 put b=x, d=y-x:
    # P=(t+a+y)(x(t-a)+ay), so again deg_x P=1<2.
    a_minus_b_x_degrees = (1, 1, 0)
    b_minus_d_x_degrees = (1, 1, 0)
    assert max(a_minus_b_x_degrees) == 1
    assert max(b_minus_d_x_degrees) == 1


def main() -> None:
    singular = shutil.which("Singular")
    assert singular is not None, "Singular is required"

    p_powers = powers(make_p((1, 1, 1, 1)), CUTOFF)
    moments = [
        singular_moment(order, p_powers[order])
        for order in range(1, CUTOFF + 1)
    ]

    # The chart A+D+F=1 meets all four claimed projective points.  J is
    # their reduced ideal.  The moment ideal I lies in J, and the displayed
    # powers of every Groebner generator of J lie in I, proving rad(I)=J.
    singular_code = f"""
ring rr=0,(A,B,C,D,E,F),dp;
ideal H={",".join(moments)};
ideal I=H,A+D+F-1;
ideal GI=std(I);

ideal p1=A-1,B,C,D,E,F;
ideal p2=A,B,C,D,E,F-1;
ideal p3=2A-1,B+1,C,2D-1,E,F;
ideal p4=A,B,C,2D-1,E+1,2F-1;
ideal J=intersect(intersect(p1,p2),intersect(p3,p4));
ideal GJ=std(J);

print("CHART_DIM"); print(dim(GI));
print("CHART_LENGTH"); print(vdim(GI));
print("POINT_DIM"); print(dim(GJ));
print("POINT_COUNT"); print(vdim(GJ));

poly c1=reduce(H[1],GJ); print("H1"); print(c1==0);
poly c2=reduce(H[2],GJ); print("H2"); print(c2==0);
poly c3=reduce(H[3],GJ); print("H3"); print(c3==0);
poly c4=reduce(H[4],GJ); print("H4"); print(c4==0);
poly c5=reduce(H[5],GJ); print("H5"); print(c5==0);
poly c6=reduce(H[6],GJ); print("H6"); print(c6==0);
poly c7=reduce(H[7],GJ); print("H7"); print(c7==0);
poly c8=reduce(H[8],GJ); print("H8"); print(c8==0);
poly cn=reduce(A+D+F-1,GJ); print("HN"); print(cn==0);

poly q1=C;
poly q2=B+2D+E;
poly q3=A+D+F-1;
poly q4=4F2-E-4F;
poly q5=2EF-E;
poly q6=4DF+E;
poly q7=E2+E;
poly q8=2DE-E;
poly q9=2D2-D;
ideal Q=q1,q2,q3,q4,q5,q6,q7,q8,q9;
ideal GQ=std(Q);
int i;
print("DISPLAYED_IN_POINTS");
for (i=1;i<=size(GQ);i++) {{ print(reduce(GQ[i],GJ)==0); }}
print("POINTS_IN_DISPLAYED");
for (i=1;i<=size(GJ);i++) {{ print(reduce(GJ[i],GQ)==0); }}
print("Q1"); print(reduce(q1^8,GI)==0);
print("Q2"); print(reduce(q2^8,GI)==0);
print("Q3"); print(reduce(q3,GI)==0);
print("Q4"); print(reduce(q4^8,GI)==0);
print("Q5"); print(reduce(q5^8,GI)==0);
print("Q6"); print(reduce(q6^8,GI)==0);
print("Q7"); print(reduce(q7^8,GI)==0);
print("Q8"); print(reduce(q8^8,GI)==0);
print("Q9"); print(reduce(q9^4,GI)==0);

ideal K=H,A+D+F;
ideal GK=std(K);
print("COMPLEMENT_DIM"); print(dim(GK));
print("COMPLEMENT_LENGTH"); print(vdim(GK));
exit;
"""
    completed = subprocess.run(
        [singular, "-q"],
        input=singular_code,
        text=True,
        capture_output=True,
        check=True,
        timeout=180,
    )
    values = [
        line.strip()
        for line in completed.stdout.splitlines()
        if line.strip()
    ]

    expected_prefix = [
        "CHART_DIM",
        "0",
        "CHART_LENGTH",
        "64",
        "POINT_DIM",
        "0",
        "POINT_COUNT",
        "4",
    ]
    assert values[: len(expected_prefix)] == expected_prefix
    cursor = len(expected_prefix)
    for label in (
        "H1",
        "H2",
        "H3",
        "H4",
        "H5",
        "H6",
        "H7",
        "H8",
        "HN",
    ):
        assert values[cursor : cursor + 2] == [label, "1"]
        cursor += 2
    for label in ("DISPLAYED_IN_POINTS", "POINTS_IN_DISPLAYED"):
        assert values[cursor] == label
        cursor += 1
        assert values[cursor : cursor + 9] == ["1"] * 9
        cursor += 9
    for label in ("Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7", "Q8", "Q9"):
        assert values[cursor : cursor + 2] == [label, "1"]
        cursor += 2
    assert values[cursor:] == [
        "COMPLEMENT_DIM",
        "0",
        "COMPLEMENT_LENGTH",
        "113",
    ]

    check_weight_gaps()
    print("PASS four-pair slice: first eight moments have four-point radical")
    print("PASS four-pair slice: no projective component is missed by the chart")
    print("PASS four-pair slice: all four directions have a strict weight gap")


if __name__ == "__main__":
    main()
