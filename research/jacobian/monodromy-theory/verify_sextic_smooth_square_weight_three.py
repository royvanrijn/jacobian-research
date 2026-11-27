"""Check the smooth-common-factor weight-three reduction over QQ."""

from fractions import Fraction
from subprocess import run


def transverse_ties(axial_weight):
    """List ties after the lower terms at weights (2,3,4) vanish."""
    ordinary = [
        {(2, 2, 2), (1, 4, 1), (0, 6, 0)},
        {(3, 2, 1), (2, 4, 0)},
        {(4, 2, 0)},
    ]
    degrees = [3 + 3 * axial_weight, 4 + 2 * axial_weight, 5 + axial_weight]
    ties = {}
    for component, degree in enumerate(degrees):
        for y in range(1, 7):
            for x in range(7 - y):
                for z in range(7 - x - y):
                    if x + y + z == 6:
                        if (x, y, z) not in ordinary[component]:
                            continue
                    elif 2 * x + 3 * y + 4 * z >= (18, 16, 14)[component]:
                        continue
                    weight = (degree - x - axial_weight * z) / y
                    if 1 < weight <= (1 + axial_weight) / 2:
                        ties.setdefault(weight, set()).add((component, x, y, z))
    return ties


assert transverse_ties(Fraction(5, 2)) == {
    Fraction(3, 2): {(1, 0, 1, 3), (2, 1, 1, 2)},
    Fraction(5, 3): {(2, 0, 3, 1)},
    Fraction(7, 4): {
        (0, 0, 6, 0), (0, 1, 4, 1), (0, 2, 2, 2),
        (1, 2, 4, 0), (1, 3, 2, 1), (2, 4, 2, 0),
    },
}
assert transverse_ties(Fraction(3)) == {
    Fraction(5, 3): {(2, 0, 3, 1)},
    Fraction(2): {
        (0, 0, 3, 2), (0, 0, 6, 0), (0, 1, 1, 3), (0, 1, 4, 1), (0, 2, 2, 2),
        (1, 0, 2, 2), (1, 0, 5, 0), (1, 1, 3, 1), (1, 2, 1, 2), (1, 2, 4, 0), (1, 3, 2, 1),
        (2, 0, 1, 2), (2, 0, 4, 0), (2, 1, 2, 1), (2, 2, 3, 0), (2, 3, 1, 1), (2, 4, 2, 0),
    },
}

SINGULAR = r"""
// Run with Singular -q and require the single success line below.
ring source = 0,(x,y,z,a,b,c,d,e,f,g,h,i,j,k,l,m,n),dp;
poly P = x*z+y^2;
poly A = P^3+a*z^4+b*y^3*z^2+c*x*y*z^3;
poly B = x^2*P^2+d*x*z^3+e*y^2*z^2+f*y^5+g*x*y^3*z+h*x^2*y*z^2;
poly C = x^4*P+i*x^2*z^2+j*y*z^2+k*y^4+l*x*y^2*z+m*x^2*y^3+n*x^3*y*z;
matrix derivatives[3][3] =
    diff(A,x), diff(A,y), diff(A,z),
    diff(B,x), diff(B,y), diff(B,z),
    diff(C,x), diff(C,y), diff(C,z);
poly determinant = det(derivatives);
matrix table = coef(determinant,x*y*z);

ring parameters = 0,(a,b,c,d,e,f,g,h,i,j,k,l,m,n),dp;
matrix rows = imap(source,table);
ideal equations;
for (int column=1; column<=ncols(rows); column++){
    equations[column] = rows[2,column];
}
ideal claims = b-2*h+n,c-2*h+n,2*f-m,2*g-m-n,
    (m-n)^4,j^3,(i+k-l)^3,(2*h-n)^2,e^3,d^3,a^2,(n^2-8*k+4*l)^3;
matrix witnesses = lift(equations,claims);
if (matrix(equations)*witnesses == matrix(claims)){
    print("Verified all twelve ideal membership identities over QQ.");
} else {
    ERROR("The claimed coefficient relations failed");
}
quit;
"""
result = run(["Singular", "-q"], input=SINGULAR, text=True, capture_output=True, check=True)
expected = "Verified all twelve ideal membership identities over QQ."
assert result.stdout.strip() == expected and not result.stderr.strip(), result.stdout + result.stderr
print("Verified both transverse tie tables and all twelve coefficient ideal identities over QQ.")
