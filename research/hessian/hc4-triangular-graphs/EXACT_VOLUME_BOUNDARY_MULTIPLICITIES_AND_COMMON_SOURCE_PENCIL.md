# Exact volume, boundary multiplicities and a common source pencil

17 November 2026. Smooth affine complex surface completions of a marked
polynomial plane; all boundary components are retained.

Let X=Spec C[x,y] be a dense open subset of a smooth integral affine
complex surface U. Suppose U has a regular unit two-form Omega restricting
to dx wedge dy, and a regular one-form theta with dtheta=Omega.

Assume the actual whole component opens M_i=X union D_i are smooth
affine, and that the individual-coordinate and whole-chart theorems in
[the individual boundary theorem](INDIVIDUAL_SOURCE_COORDINATES_AND_UNCONDITIONAL_BOUNDARY_CHARTS.md)
apply to every component D_i of U minus X. These hypotheses are supplied
for the ACTUAL full canonical plane Keller completion by the cited
theorem and
[canonical topology](CANONICAL_ETALE_COMPLETION_BOUNDARY_TOPOLOGY_AND_FREE_DECK_ACTION.md).

Then there is ONE polynomial SOURCE coordinate y regular on all of U.
The associated source translation extends to U by the
[unit-volume Hamiltonian proof](HAMILTONIAN_COMMON_COORDINATE_AND_ADDITIVE_COMPLETION.md).

The proof below supplies the missing simultaneous alignment argument.
It does not prove the Keller inverse, remove ramified primes outside U,
or make a target coordinate into a source coordinate. Every classical
surface-resolution input is stated at its actual scope; no Xie theorem
or classification of arbitrary locally nilpotent derivations is used.

## 1. The individual whole-chart data and the strict pole inequality

For each D_i choose source coordinates (x_i,y_i), scaled to have the
same source volume. The individual-coordinate theorem makes y_i
regular on M_i. The whole normalized root-chart theorem gives

    y_i-lambda_i=z^r_i,
    x_i=z^(1-r_i)t+h_i(z),                 r_i>=1,

where h_i is a finite Laurent polynomial after the allowed polynomial
t-shear. Exactness of the SAME volume on M_i kills its z^-r_i
coefficient. Separatedness of the marked completion forces the residual
principal part to be nonzero. Therefore, writing its largest pole p_i,

    v_i(x_i)=-p_i,   v_i(y_i-lambda_i)=r_i,   p_i>=r_i+1.       (1)

Here v_i is the PRIMITIVE divisor valuation at the original D_i.
The root cover is etale at D_i; it has ramification index one at every
chosen lift, so its z-orders are these original orders. A generic t
cannot cancel the leading h_i term, because p_i>r_i.

For later use, it is the STRICT inequality p_i>r_i that matters.
Unit volume alone does not give it. Both the shared primitive and
the actual separated whole source/component chart enter this step.

## 2. Boundary duals and the negative definite simultaneous matrix

Take a smooth projective surface W obtained by point blowups of the
source P2, extracting all v_i, and resolve all rational maps used below.
This is possible by resolution of birational maps between smooth
complex surfaces and of a rational map from a smooth surface to P2.
Every centre lies over the original line at infinity. Extra blowups
do not change the dual pairings of already extracted prime divisors.

Let E_0,...,E_N be the full boundary W minus X. Their classes are a
basis of Pic(W) tensor R: start with the line at infinity and add one
exceptional class per blowup. In the total-transform basis H,F_1,...,F_N,

    H^2=1, F_j^2=-1, H.F_j=F_j.F_k=0 for j!=k.                (2)

Define E_i^vee by E_i^vee.E_j=delta_ij. The matrix of their pairings
is the inverse of the full boundary intersection matrix.

There are two algebraically independent functions a,b in O(U),
since U is affine and its function field is C(x,y). Resolve [a:b:1]
on W and put L=phi^*(line at infinity). It is effective and

    L^2=[C(x,y):C(a,b)]>0.                                   (3)

Its support contains NONE of the retained D_i: both a and b are
regular at their generic points, so [a:b:1] is in the affine target
there. Thus L is a linear combination of all OTHER boundary classes.
In particular

    L.D_i^vee=0 for every i.                                 (4)

A real vector of positive square in the form (2) has negative
definite orthogonal complement. Explicitly, after writing
L=(A,B), A^2>||B||^2, any orthogonal vector (c,d) has
c=(B.d)/A, hence c^2-||d||^2<0 unless d=0. The independent
D_i^vee lie in that complement. Consequently their full matrix is
NEGATIVE DEFINITE.

This proves the required simultaneous sign using the ACTUAL affine
ring O(U). It does not infer a coordinate from a valuation matrix.
For two distinct retained components, if

    (D_i^vee)^2=-n_i,  D_i^vee.D_j^vee=c_ij,

then

    n_i>0,  n_j>0,  n_i*n_j-c_ij^2>0.                         (5)

## 3. The exact local multiplicity bound n_i<=p_i^2

Work first with one valuation and omit i. In its adapted coordinates,
the first source P2 centre is [1:0:0]. After blowing it up the chart
has coordinates

    u=1/x,  w=y-lambda.

The next centre is q=(u,w)=(0,0), a FREE point on the first
exceptional curve, away from the strict line at infinity. All remaining
ancestral centres lie over q. We have v(u)=p and v(w)=r with p>r>0.

Let m_1,...,m_l be the multiplicities of the successive point centres
over q for this divisor valuation. Equivalently, m_j is the coefficient
of the final divisor in the total transform of the exceptional curve
created by the j-th local blowup. These integers are positive; the
last is one. The local relative canonical divisor and the unit volume
give

    ord_v(du wedge dw)=sum_j m_j,
    du wedge dw=-u^2*(dx wedge dy),
    sum_j m_j=2p.                                            (6)

For clarity, the log discrepancy is 1+sum m_j=2p+1.
There is no additional one in the middle equality of (6).

In the total-transform basis (2), the dual of D has coefficient p
on H and on the first exceptional class F_0, and coefficients -m_j
on the subsequent local exceptional classes. Thus

    D^vee=p*(H-F_0)-sum_j m_j*F_j,
    n=-(D^vee)^2=sum_j m_j^2.                                (7)

One can also derive this from proximity: the coefficient of a total
exceptional pullback at D is the valuation of its centre ideal.
Other blowups have coefficient zero. Formula (7) is unchanged by
extra blowups elsewhere.

Put d=gcd(p,r). Successively blowing up the intersection points
forced by the two values (p,r) subtracts the smaller value from the
larger until (d,d) is reached. Including the final blowup at equal
values, its multiplicity block is the subtraction Euclidean sequence.
It satisfies

    sum_block m=p+r-d,       sum_block m^2=p*r.                (8)

These identities follow by induction: remove min(p,r); in the
square sum p*r-min(p,r)^2 is the product of the new pair.
The equal-value base has the one term d.

Until that block ends the centres are forced by unequal values.
At equality, the ratio of the local parameters either has
transcendental residue, in which case the divisor has been reached,
or has a nonzero constant residue and one continues above a point
of the new exceptional divisor. Each later centre multiplicity is
at most d. Indeed its ideal contains a local equation of that
exceptional divisor with value d, and successive point multiplicities
never increase: in blowup coordinates their minimum is at most the
previous minimum. This also covers every later satellite centre.

By (6)-(8), the remaining multiplicities have total

    2p-(p+r-d)=p-r+d.

Therefore the ENTIRE, unbounded extraction sequence obeys

    n<=p*r+d*(p-r+d)<=p^2.                                   (9)

For the second inequality use d<=r and d<=p-r:

    d*(p-r+d)<=d*p<=p*(p-r).

No bound on the extraction length or on the number of characteristic
Puiseux exponents is used. If the Euclidean block were terminal,
the same estimate holds with an empty tail; in the present
unit-volume p>r setting (6) actually requires a nonempty tail.

## 4. A polynomial coordinate is monic in a transverse coordinate

Let (x,y) be source coordinates and g another polynomial source
coordinate. If g is not affine in y, its x-degree m is positive
and its leading x-coefficient is a nonzero CONSTANT.

Here is a proof that does not invoke a normal form for plane
automorphisms. Choose any c. The fibre g=c is a closed affine line,
with a polynomial parameter t furnished by a complementary source
coordinate. The restriction y(t) is nonconstant. Otherwise this
entire irreducible closed curve is the line y=lambda, and equality
of their prime ideals gives g-c=constant*(y-lambda), contrary to
the assumption.

Since y(t) is nonconstant, C[t] is finite over C[y(t)], so x(t)
is integral over that ring. A monic polynomial relation M(x,y)=0
on this fibre belongs to its ideal (g-c) in C[x,y]. In

    M=(g-c)*N,

comparison of the leading x-coefficients gives
1=lead_x(g)*lead_x(N) in C[y]. Hence lead_x(g) is a unit in C[y],
which proves the assertion.

Consequently, for any valuation with v(x)=-p<0 and v(y)>=0,

    v(g)=-m*p,                                               (10)

because every lower x-power, with coefficients polynomial in y,
has strictly larger valuation. Constants in g or shifts of y
do not affect the conclusion.

## 5. Two different source pencils force an excessive cross pairing

Consider two components with adapted coordinates from Section1.
If their coordinate pencils differ, y_1 is not affine in y_2.
Apply Section4 in the source frame (x_2,y_2):

    -v_2(y_1)=m*p_2=:a>0,  m=deg_(x_2)(y_1)>=1.              (11)

We calculate the dual cross pairing in the source P2 frame (x_1,y_1).
The extraction of D_1 starts at A=[1:0:0], with first multiplicity
p_1, and then enters the finite-y_1 free point y_1=lambda_1.
Let b be D_2's coefficient on the total line at infinity.

If D_2 starts at a different point, then b=a: either y_1 has larger
pole than x_1 or their poles are equal and their ratio has a finite
nonzero residue. The two clusters share no point. Formula (2) gives

    D_1^vee.D_2^vee=p_1*b=p_1*a.

If D_2 starts at A, then b=-v_2(x_1)>a and the first multiplicity
is min(v_2(1/x_1),v_2(y_1/x_1))=b-a. Its next centre is the
intersection of the first exceptional curve with the strict line
at infinity, since y_1 has a POLE there. Thus it does not enter
D_1's finite-y_1 free point; the clusters share exactly their
first point. Now

    D_1^vee.D_2^vee=p_1*b-p_1*(b-a)=p_1*a.

Both cases, including any later extra blowups, give the exact formula

    c_12=m*p_1*p_2>=p_1*p_2.                                 (12)

An equivalent way to see it is that no D_1 tail centre is an ancestor
of D_2 when y_1 has a pole, so only H and the first exceptional
class contribute to their dual pairing.

## 6. Alignment, regularity, and the remaining ramified boundary

Combining (5), (9) and (12) under the supposition of different pencils,

    0<n_1*n_2-c_12^2
       <=p_1^2*p_2^2-m^2*p_1^2*p_2^2<=0,

a contradiction. Hence every pair of the chosen source-coordinate
pencils is equal. Two polynomial source coordinates in the same
pencil differ by an affine transformation; alternatively, the
Section4 monicity statement says that non-affine dependence would
already have positive transverse degree and force a pole.

Choose y=y_1 if the boundary is nonempty. For each i,
y_i=a_i*y+b_i with a_i!=0, so y is regular on M_i. The M_i together
with X cover U. Equivalently every boundary valuation of y is
nonnegative, so normality extends y everywhere. Thus y in O(U)
is ONE polynomial SOURCE coordinate. If the boundary is empty,
U=X and any source coordinate works.

The conditional whole-root-chart result now applies in one common
source frame. The Hamiltonian argument gives the nonzero LND
partial_x on O(U), kernel C[y], preserving the entire source X and
each boundary component. Each polynomial replica q(y)partial_x
also preserves U.

For an actual full canonical Keller completion, this discharges the
specific common-source-coordinate existence input and its whole
algebraic boundary/root-chart consequences. It does NOT say that U
contains the omitted RAMIFIED primes of the full finite normalization,
that the LND preserves B or F, that first velocities are B-integral,
that the generic extension is Galois, or that HC4/JC2 is solved.
Any further classification claims, unique-kernel assertions or deck
conclusions require their own complete proof dependencies.

## 7. Verification scope

The checker verifies the subtraction Euclidean identities exactly,
checks the local discrepancy/dual formulas on ancestral free and
satellite blowup trees, exercises coprime and noncoprime multiplicity
blocks, checks coordinate monicity for genuine polynomial
automorphism words, and checks the full two-pencilled determinant
contradiction. The p=r and p<r controls show why the strict
whole-chart pole inequality is indispensable. Lorentz matrix
controls distinguish the actual positive-pullback argument from a
valuation-only coordinate inference.

The proof above is unrestricted. Finite controls support the
transcription; they do not prove completeness for all blowup trees.
Point-blowup extraction, resolution of a surface rational map and
the elementary intersection/proximity identities are adopted
foundations. The individual-coordinate theorem retains its precise
Borisov source import. No formal verification, independent
whole-proof replay, external review or novelty is claimed.

Replay:

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_exact_volume_common_pencil.py
