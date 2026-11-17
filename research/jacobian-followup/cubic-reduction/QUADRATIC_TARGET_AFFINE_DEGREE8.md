# No affine frame of the first nonlinear target layer can improve raw13 (KCUB43)

Subsequent [KCUB44](POSITIVE_SOURCE_QUADRATIC_FIELD_OBSTRUCTION.md) proves
that every nonzero positive-axis quadratic source shear has no rational
Keller repair through source degree eight. It also extends the original
target-field intersection unchanged through degree nine. The transformed
field at degree nine remains unclassified;prior positive-source-gate
proposals below are stage history.

Let H be the saved degree-seven map and set

    H_(L,M)=H+(0,L,M)*H1^2.

Over every characteristic-zero field, every independent invertible affine
source/output frame of a nonzero shear (L,M)!=(0,0) has at least **sixteen
nonlinear monomial occurrences**. The saved degree-eight control
(H1,H2+H1^2,H3) has nineteen, so the minimum over this nonlinear family lies
in **[16,19]**. Sixteen is not claimed attained and the exact minimum in
that interval is not settled here.

This completes the raw13-improvement gate posed by
[KCUB42](RAW13_TARGET_FUNCTION_FIELD_RIGIDITY.md). Its complete rational-target
classification through source degree eight reduces every such target
operation, even with arbitrary target degree or denominators, to this
family followed by an affine output map. Allowing affine source maps as
well cannot reduce the count below thirteen. The entire class has sharp
minimum **thirteen**, attained only by the zero shear in the equality
frames classified by [KCUB41](RAW13_AFFINE_ORBIT_MINIMUM.md): original source
centre and monomial source/output linear parts. No nonzero shear retains13.

There is no new smaller witness, global simplicity optimum or unrestricted
whole-profile exclusion. Positive bounds and the 62 excluded / 12 remaining
strict profiles do not change. Nonlinear source changes, source/target
coordination beyond affine sources and changed constructions remain open.

## Complete centre and parameter reduction

Write X=(a,b,c), t=(alpha,beta,gamma), u=(0,L,M), and P=H1(X+t), P0=H1(t).
An invertible output shear removes the constant-dependent cross term:

    H(X+t)+u*P^2
      = (I+2*P0*u*e1^T)*(H(X+t)+u*(P-P0)^2)-u*P0^2.   (1)

The matrix has determinant one and inverse I-2*P0*u*e1^T, because e1^T*u=0.
Since all output affine frames are allowed, we can study the centered-square
map on the right. This normalization drops no parameter or centre branch.
It is an output equivalence, not an assertion that the individual frame's
term count is preserved by that shear.

Let C_d be the homogeneous output coefficient matrix of H(X+t), let
Q_d=((P-P0)^2)_d, and let r_d be the output-span rank of

    C'_d=C_d+u*Q_d.

Every source linear substitution acts invertibly on homogeneous polynomials,
and output mixing preserves their span. Thus, just as in KCUB41,

    nonlinear occurrences >= sum_(d=2..8) rank(C'_d).    (2)

The exact rank profiles, for all centres and every parameter-zero branch,
are:

| Shear | Centre | r2 | r3 | r4 | r5 | r6 | r7 | r8 | Sum |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| L=M=0 | alpha=beta=0 | 2 | 3 | 3 | 2 | 2 | 1 | 0 | 13 |
| L=M=0 | off that axis | 3 | 3 | 3 | 2 | 2 | 1 | 0 | 14 |
| L=0,M!=0 | alpha=beta=0 | 2 | 3 | 3 | 2 | 2 | 1 | 1 | 14 |
| L=0,M!=0 | off that axis | 3 | 3 | 3 | 2 | 2 | 1 | 1 | 15 |
| L!=0 | alpha=beta=0 | 2 | 3 | 3 | 2 | 2 | 2 | 1 | 15 |
| L!=0 | off that axis | 3 | 3 | 3 | 2 | 2 | 2 | 1 | 16 |

No coefficient grid, generic-point assumption or polynomial ideal solver is
used to prove these profiles. The complete translated matrices have483
entries with all five parameters formal; independent exact reconstruction
checks them and every identity below.

## Rank stability in degrees two through six

KCUB41 supplies the full base ranks of C_d. Two simple observations make
all added-row arguments uniform: H is affine in c and its a exponent is
at most three, even after translation. Therefore a c^2 term or an a
exponent at least four separates a new square row from every old row.
For alpha!=0 the following coefficients of Q_d are nonzero:

| d | Source monomial | Coefficient |
| --- | --- | --- |
|2|c^2|alpha^6|
|3|a*c^2|6*alpha^5|
|4|a^2*c^2|15*alpha^4|
|5|a^3*c^2|20*alpha^3|
|6|a^4*c^2|15*alpha^2|

Thus Q_d is outside the old output span in each case. Independence of
three old rows in degrees2,3,4 cannot be destroyed by adding constant
multiples of an independent fourth polynomial. In degrees5,6 the first
row is zero and remains zero (the first entry of u is zero); independence of the other two rows
likewise survives. This proves all alpha!=0 ranks in those degrees.

At alpha=0, the complete square rows are

    Q2=a^2,
    Q3=-6*beta*a^3,
    Q4=(9*beta^2-2*gamma)*a^4-6*a^3*b,
    Q5=6*beta*gamma*a^5+18*beta*a^4*b-2*a^4*c,
    Q6=gamma^2*a^6+6*gamma*a^5*b+6*beta*a^5*c+9*a^4*b^2.

For degree3 the base determinant on (a^2*b,a*b^2,a*b*c) is -432. Q3
vanishes on those columns, so it is outside their span if beta!=0; if
beta=0 it is zero and the rank is unchanged. For degree4 the base
determinant on (a^3*c,a^2*b*c,a*b^3) is -336. Q4 vanishes there and is
never zero, so is outside the span. In degrees5,6 the fixed coefficients
-2 on a^4*c and9 on a^4*b^2 separate the square row from the old rows.

In degree2 with beta!=0 the first base row is exactly -3*beta*a^2. Adding
L*Q2,M*Q2 to the other rows is just an invertible row operation, retaining
rank three. If beta=0 the first row is zero, and the other two have minor24
on (a*c,b^2); Q2 vanishes on these and has a nonzero a^2 column outside
their span. Their rank remains exactly two. This includes all gamma,L,M
values without divisions by unknown parameters.

## Degrees seven and eight

The complete top forms are

    H3_7=8*a^3*b^3*c,
    Q7=6*a^5*b*c+2*gamma*a^6*c+6*alpha*a^5*c^2,
    Q8=a^6*c^2,
    (H_(L,M))_7=(0,L*Q7,H3_7+M*Q7),
    (H_(L,M))_8=(0,L*Q8,M*Q8).

Q7 is independent of H3_7, since its a^5*b*c coefficient is6 and H3_7
has only its distinct monomial. Hence rank7 is two for L!=0 and one for
L=0. Rank8 is one exactly when the shear is nonzero. The rank table follows.
It already excludes every improvement on13 by a nonzero shear. The next
three short arguments strengthen its uniform lower bound to16.

## Necessary low-count conditions from the top and cubic terms

First, a nonmonomial polynomial square has at least three monomials in
characteristic zero. A two-term Laurent polynomial is squarefree: apply a
logarithmic partial derivative in a coordinate where the two exponent
vectors differ; modulo that binomial the derivative is a nonzero scalar
Laurent unit, so it is coprime to the binomial. A nonunit square cannot be
such a binomial. A polynomial which is not a monomial is not a Laurent
unit, proving the assertion.

The top scalar form after a source matrix is

    (ell_a^3*ell_c)^2.

If it has at most two monomials it must therefore have one. Unique
factorization makes ell_a and ell_c scalar coordinate variables; their
independence makes those coordinates distinct. Permuting and scaling source
coordinates, we may use a,c themselves as coordinate axes and write the
remaining source row as B=b+eta*a+zeta*c.

Second, if L=0,M!=0, the degree-seven scalar form in these coordinates is

    a^3*c*(8*B^3+6*M*a^2*B+2*M*gamma*a^3+6*M*alpha*a^2*c).  (3)

It cannot be one monomial. Its b^3 coefficient forces a purported unique
monomial to be proportional to a^3*b^3*c. The b^2 coefficients in the
bracket are24*eta on a*b^2 and24*zeta on b^2*c, forcing eta=zeta=0.
Its remaining a^2*b coefficient is6*M!=0, a contradiction. This holds
at every centre, not only on the c-axis.

Third, at a c-axis centre the square has no degree-three part, so all
shears have the same cubic output span as H. In the restricted source
frame this is spanned by

    -gamma*a^3-3*a^2*B,
    12*gamma*a^2*B+24*a*B^2,
    6*a*B*c.                                          (4)

There are always three nonzero source columns a^2*b,a*b^2,a*b*c. If exactly
three cubic occurrences are allowed after arbitrary output mixing, there
can be no other nonzero source column: invertible output matrices cannot
annihilate a nonzero column. The a^2*c coefficient in the first row is
-3*zeta, forcing zeta=0. Its a^3 coefficient forces gamma=-3*eta.
The a^3 coefficient in the second row is then -12*eta^2, forcing eta=0
and gamma=0. Thus a three-occurrence cubic layer, together with coordinate
axes a,c, requires a monomial source matrix and the original centre.
All displayed equality coefficients are independently reconstructed from
the full translated original map.

At that centre and a monomial source frame, every nonzero shear has at
least nineteen distinct source columns: thirteen from H and six disjoint
columns from H1^2. The new columns are nonzero multiples of (0,L,M), so
none can be killed by invertible output mixing. This is the KCUB42 fixed-
source argument, here checked directly and used only after rigidity has
forced the source frame to be monomial.

## Exhausting every count below sixteen

Assume a nonzero shear has at most15 nonlinear occurrences in an affine
frame.

For L!=0 off the c-axis, the rank sum is already16. On the c-axis it is15,
so every degree must attain its rank bound. In particular degree8 has one
occurrence, forcing coordinate axes a,c, and degree3 has three, forcing
the original centre and a monomial source matrix by (4). The nineteen-column
bound contradicts the assumed total.

For L=0,M!=0 off the c-axis the rank sum is15. Attainment requires one
occurrence in both degrees8 and7. The first forces axes a,c, and (3) then
forbids the second.

For L=0,M!=0 on the c-axis the rank sum is14, leaving at most one extra
occurrence. Degree8 thus has at most two; the square argument again forces
axes a,c. Equation (3) forces at least two degree-seven occurrences,
consuming the one extra. Therefore degree8 has one and degree3 exactly
three. Equation (4) forces the original monomial source frame, again
contradicting nineteen source columns. This exhausts every field and every
zero/nonzero parameter branch, proving the lower bound16.

The upper19 is the previously verified map (H1,H2+H1^2,H3), with degree8,
identity tangent map, determinant1 and the saved distinct integer collision.
Its exact full determinant and collision are replayed. The interval[16,19]
for the nonlinear subfamily is not collapsed to a sharp minimum here,
since the whole family is already worse than raw13 for the original goal.

## Reproduction, scope and next route

Run `python quadratic_target_affine8.py`, then
`python verify_quadratic_target_affine8.py`. The producer writes complete
matrices and two generic equality models without a solver. The independent
checker instead translates the original map using Fraction composition,
checks all483 shear entries, all base entries and square rows, exact
centre-dependent output equivalence, every rank-separation coefficient and
zero-centre pivot, both coupled equality obstructions, nineteen original-
centre columns and the matching upper control. See
[artifact index](quadratic-target-affine8/README.md) and
[validation](quadratic-target-affine8/validation.json).

The proof imports KCUB41's already verified base ranks and rank-invariance
argument, and KCUB42's complete degree-eight rational-target classification.
Universal span independence, the characteristic-zero Laurent-binomial
argument, UFD and equality exhaustion are written, not formalized or
externally reviewed. No coefficient grid, primary ideal answer, solver
query or new timeout appears in this result. Imported results retain their
own written/external boundaries.

The complete class obtained from H by arbitrary rational target changes,
affine source/output changes and an ordinary-degree<=8 endpoint has sharp
raw minimum13. Equality requires a total affine target and the KCUB41
monomial frames; every total nonlinear target has at least16 occurrences.
This closes the proposed joint affine degree-eight gate without asserting
an unrestricted degree-eight exclusion or any global optimum.

A next route must change the source nonlinearly or the construction. One
specific proposed gate is a quadratic positive-axis source shear
(a,b,c)->(a+q(b,c),b,c), with all three quadratic coefficients, jointly with
uncapped rational target repair. KCUB19 concerns one-axis changes of the
separate twelve-variable cubic and does not settle this three-variable
field problem. Its entire parameter-zero boundary must be retained before
any rank or candidate calculation. Unweighted lower-degree constructions
and the separate cubic dimension/support objectives remain open alternatives.
