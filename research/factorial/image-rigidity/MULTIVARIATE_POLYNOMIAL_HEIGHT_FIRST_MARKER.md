# A real first-marker test for polynomial row heights

Take a finite ordinary two-monomial core and finitely many
additional unequal row types with fixed positive gaps
and orientations. Let their smaller exponents be real
polynomials \(\ell_i(T)\) with rational coefficients,
either nonnegative constants or nonconstant polynomials
with positive leading coefficients. The multiplicities
are arbitrary nonnegative integers. A type with a positive
multiplicity is called active. Only ordinary integer values
of the active rows are used in the conclusions.

**Theorem.** Assume at least one active \(\ell_i\) has odd
degree. With any fixed admissible vector of multiplicities,
only finitely many integer heights have a common root at any
order. Even with independently varying multiplicities, the
number of exceptional integer heights \(0\le T\le H\) is
bounded by \(C(\log(2+H))^\kappa\). At fixed bounds on the
number of types, their gaps and their polynomial degrees,
the exponent and the constant can be taken uniform over the
core and orientations, with a logarithmic dependence on a
bound for the polynomial coefficients specified below.

The same conclusions hold under a more general real marker
condition. For each terminal below the absolute cutoff, and
each allowed nonempty set of active types, require that at
least one active \(\ell_i(T)+1/N\) have a real zero. Its
largest real zero among all active types can be used regardless
of multiplicity of that zero.

## Virtual positive specialization at the first marker

Fix multiplicities and orders \(n,N=n+1\), with \(n\ge2\).
As in the [affine first-marker proof](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md),
strip the factorials at the smaller endpoint in every order
to obtain two polynomials in \(T,z\). At order \(j\), the
moving row of gap \(d_i\) contributes
\[
            \prod_{h=1}^{d_i k_i}(j\ell_i(T)+h)
\]
at index \(k\), where \(k_i=k\) for increasing and
\(k_i=j-k\) for decreasing rows. The fixed core coefficients
are positive. For admissible integer \(T\) these arrays
are the original moments up to nonzero scalar factors.

Let \(T_*\) be the largest real zero of the finitely many
active \(\ell_i(T)+1/N\). A positive-leading odd-degree
polynomial passes through \(-1/N\), ensuring existence in
the stated main case. For \(T>T_*\) all active
\(\ell_i(T)>-1/N\), by their eventual positivity and the
definition of \(T_*\). Thus at \(T_*\) every lower factor
\(n\ell_i(T_*)+h\) is positive, and the only vanishing upper
factors are the first factors \(N\ell_i(T)+1\) in the marked
rows with \(\ell_i(T_*)=-1/N\).

Write \(q_i\) for the positive order of vanishing of such
a first factor at \(T_*\). Its first nonzero Taylor
coefficient there is positive, by positivity on the right.
Let \(u\) and \(v\) be the sums of \(q_i\) times the
multiplicity over increasing and decreasing marked rows,
respectively. The vanishing orders of the upper coefficients
are \(v\) at zero, \(u\) at \(N\), and \(u+v\) at every
interior index. Divide that polynomial by
\((T-T_*)^{\min(u,v)}\). The resulting polynomial at \(T_*\)
is a positive constant or positive multiple of \(z^N\)
when \(u\ne v\), and a positive endpoint binomial when
\(u=v\). The lower specialization has positive endpoints
and still has degree \(n\). This gives coprimality immediately
in the unbalanced cases.

If \(u=v\), some marked row occurs in each orientation.
In each marked row the initial lower Gamma shape is
\(1+n\ell_i(T_*)=1/N\). Reserving one Beta factor from each
orientation as in the
[balanced affine specialization](MULTIVARIATE_AFFINE_BASE_RAY_ALL_DIRECTIONS.md#balanced-cohorts-have-a-positive-geometric-specialization)
turns the lower coefficients into a nonconstant positive
geometric moment sequence, up to a scalar. The other
factors are positive Gamma moments irrespective of \(q_i\).
The lower nonreal arguments lie strictly between consecutive
rays of the upper endpoint binomial by the cited geometric
angle boxes. The opposite parities exclude a shared negative
real argument. Thus the two specialized polynomials are
coprime also in the balanced case.

The fixed-degree resultant of the lower polynomial and
the divided upper polynomial is nonzero at \(T_*\),
even in a degree-drop case, by the product formula with
the lower leading coefficient nonzero. Hence the original
resultant is a nonzero polynomial in \(T\). It has finitely
many admissible zeros. Terminal two follows from real
positivity of an even ordinary moment, and the
[absolute eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
handles all sufficiently large terminals.

For example, at \(n=2,N=3\), take one increasing unit-gap
row with smaller exponent \(((3T+1)^3-1)/3\), and three
decreasing unit-gap rows with smaller exponent \(T\).
At \(T_*=-1/3\) their respective orders of vanishing are
three and one, so \(u=v=3\). The lower specialization is
\(4z^2/9+2z/81+64/729\) and the divided upper one is
\(54z^3+216\).

## Independent multiplicities

The
[freezing and Pfaffian block argument](MULTIVARIATE_AX_FREEZING_POLYLOG_HEIGHTS.md)
uses only two facts about the affine heights. The
factorial-stripped coefficients are products of fixed-degree
positive polynomials to variable powers on the admissible
real domain. No resultant vanishes identically when the
counts are fixed integers at an admissible witness.
Both hold here by the preceding specialization. Split into
the finitely many possible active sets, and in the real
box restrict to \(\ell_i(T)\ge0\) for those types. The base
polynomials of their logarithms are then positive.
Inactive types are omitted.

If \(E\) bounds the polynomial degrees, choose \(A\ge1\)
bounding the sum of absolute values of coefficients of each
\(\ell_i\), plus \(d_i+1\). For \(0\le T\le H\) the maximum
active row endpoint is at most \(A(1+H)^E\).
The [absolute curvature twelve bound](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
therefore gives \(r_i\le12A(1+H)^E\) at any candidate
collision. The proof of the block estimate applies with
this multiplicity box and with auxiliary logarithms of
polynomials of degree at most \(EJ\max_i d_i\), where
\(J\) is the absolute cutoff. Its format bounds depend on
the number of types, degree and gap bounds and \(J\).
Thus the number of exceptions is
\[
                         O\big((\log(2+A+H))^\kappa\big)
\]
with uniform effective constants at fixed number of types,
degree and gap bounds. Coprimality transfers to every
characteristic-zero field through the rational resultant.
The count does not rule out infinitely many sparse exceptional
heights when the multiplicities vary.
