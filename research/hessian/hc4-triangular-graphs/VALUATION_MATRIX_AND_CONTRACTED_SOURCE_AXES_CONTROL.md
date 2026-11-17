# A negative valuation matrix can contract the marked source axes

HC4MYVAC0. 17 November 2026. Exact comparison control for the
common-coordinate bridge; NOT an actual Keller example.

Two primitive divisorial valuations at infinity satisfy

    v1(x)=-1, v1(y)=2,
    v2(x)=2, v2(y)=-1.

Both give order zero to the source volume, both admit regular local
volume primitives, and their exact unnormalized Xie dual matrix is

    M=[[-2,1],[1,-2]],

which is strictly negative definite. Their simultaneous nonnegative
polynomial ring nevertheless contains NO nonconstant polynomial
SOURCE coordinate.

The missing actual-source condition is explicit: the induced birational
map contracts both coordinate axes. Every polynomial pair in this
intersection ring has Jacobian divisible by x*y and is therefore not
a Keller pair. The control isolates an invalid valuation-only inference;
it does not disprove common-coordinate existence for an actual Keller
completion.

## 1. The valuations and exact volume charts

For a nonzero polynomial sum c_ij*x^i*y^j define

    v1=min_ij(-i+2j), v2=min_ij(2i-j).

There is no cancellation between distinct monomials with equal order:
their leading coefficients are distinct powers of a transcendental
residue coordinate. These valuations extend to C(x,y).

They are primitive divisorial valuations. Near their generic divisors
use the rational substitutions

    v1: x=s^-1, y=s^2*t,
    v2: x=s^2*t, y=s^-1,

with t transcendental along s=0. Then

    v1: dx wedge dy=-ds wedge dt,
    v2: dx wedge dy= ds wedge dt.                       (1)

In particular both volume orders are zero. Even a regular local
primitive is retained: the source polynomial one-form

    alpha=(x*dy-y*dx)/2, d alpha=dx wedge dy

pulls back respectively to

    (3t/2)*ds+(s/2)*dt,
    -(3t/2)*ds-(s/2)*dt.                                (2)

Thus merely adding regular exact-volume data at the two divisorial
charts does not remove this control.

## 2. The complete polynomial intersection ring

Put

    T=x*y, U=x^2*y, V=x*y^2.

A monomial is nonnegative for both valuations exactly when

    i,j>=0, i<=2j, j<=2i.

This cone has integer Hilbert basis (1,1),(2,1),(1,2). Explicitly, if
i>=j its monomial is U^(i-j)*T^(2j-i); if j>=i it is
V^(j-i)*T^(2i-j). All exponents are nonnegative.
Because equal-order leading terms do not cancel, the FULL intersection is

    A=C[x,y] intersect O_v1 intersect O_v2
      =C[T,U,V]=C[T,U,V]/(U*V-T^3).                     (3)

The notation on the right denotes the abstract presentation. It has
transcendence degree two. Its fraction field is C(x,y), since

    x=U/T, y=T^2/U.

The binomial relation is irreducible and is the full presentation kernel.

Every nonconstant monomial of A contains both x and y. Therefore any
f in A has the form

    f=c+x*y*g(x,y), c in C, g in C[x,y].                 (4)

If f were a nonconstant polynomial source coordinate, its fibre f=c
would be an irreducible reduced affine line. Equation (4) instead
makes that fibre contain BOTH source coordinate axes as components.
This is impossible. Thus A contains no nonconstant source coordinate.

## 3. The exact dual intersection matrix

For completeness the matrix is computed from a smooth projective
blowup model, rather than inferred from a normalization convention.
Start from P2 with its line at infinity. At each of the two distinct
points giving the x-infinity and y-infinity directions, make THREE
successive free blowups, taking zero for the two finite centers.
The last divisors give precisely the substitutions in Section1.

Order the seven boundary components as

    L,E11,E12,E13,E21,E22,E23.

Their intersection matrix I has diagonal

    (-1,-2,-2,-1,-2,-2,-1),

and off-diagonal entries 1 on exactly

    L-E11-E12-E13 and L-E21-E22-E23.

Its determinant is one. The divisor dual to E13 is the column of I^-1
whose intersection with every component is the corresponding Kronecker
delta; similarly for E23. Consequently the dual-intersection submatrix
for the two retained valuations is

    [I^-1]_(E13,E23)=[[-2,1],[1,-2]].                    (5)

Its eigenvalues are -1 and -3, so it is strictly negative definite.
Both primitive valuations have min(v(x),v(y))=-1, so normalized
and unnormalized matrices coincide in this example.

This is exactly the dual-divisor convention in Xie,
[Intersections of valuation rings in k[x,y], Theorem1.2, PDFpage3]
(https://archive.ymsc.tsinghua.edu.cn/pacm_download/423/9433-affinegeo20140928.pdf).
The primary definition and theorem were inspected. That theorem
characterizes transcendence degree two of the polynomial intersection;
it supplies no polynomial-coordinate or whole-source assertion.
Our direct computation of (3),(5) is independent of its proof.

## 4. The actual Keller condition fails visibly

The natural map

    A2_(x,y) -> Spec(A),
    (x,y) -> (T,U,V)=(xy,x^2y,xy^2)

is birational but contracts EACH entire source axis to the same vertex.
It is not the required open immersion of the whole marked source into
a finite normalization/completion.

More strongly, if P,Q in A, equation (4) writes

    P=c+x*y*g, Q=d+x*y*h.

A direct expansion gives

    Jac(P,Q)=x*y*[
       x*(g_x*h-g*h_x)+y*(g*h_y-g_y*h)
                         +x*y*(g_x*h_y-g_y*h_x)].       (6)

Thus every pair in A has a Jacobian vanishing along both axes.
There is no constant nonzero Jacobian pair in this intersection ring.

The control therefore does NOT satisfy actual Keller source geometry.
Its role is precise: zero volume orders, local exact primitives,
transcendence degree two, and strict negative definiteness alone
do not prove a common regular source coordinate. A proof of the
divisorial common-coordinate theorem must use an additional actual
whole-source constraint that excludes such contractions. The
Hamiltonian bridge separately shows that, once a common source
coordinate exists, whole normalized boundary root charts are not
needed to extend its translation field.

No HESS-COMMON-PENCIL-INPUT, bundled HESS-MP, HC4 or JC2 conclusion
is upgraded here.

## 5. Exact checks and replay

Replay with

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_valuation_matrix_source_axes.py

The checker verifies the two substitutions and primitives, the Hilbert
basis decomposition, presentation and fraction-field identities, the
complete blowup intersection model and its dual matrix, the contracted
axes and the unrestricted Jacobian factorization. Finite monomial
controls support the explicit cone proof; they do not replace it.
Internal exact replay is not formal verification or external human
review.
