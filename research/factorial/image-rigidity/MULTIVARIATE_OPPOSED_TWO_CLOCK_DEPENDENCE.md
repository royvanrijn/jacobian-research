# Two logarithmic clocks on opposed factorial branches

Fix an ordinary two-monomial core, two fixed positive gaps, and a
fixed signed tail gap. Let the two opposed row heights and their
multiplicities vary independently subject to that signed gap, as in
the [local branch theorem](MULTIVARIATE_OPPOSED_DIFFUSE_PUISEUX_HEIGHT_BRANCHES.md).
At a positive interior Gaussian contact write
\[
 w_*=\exp(-T_* /(2nN)),\qquad v_*=\exp(-\lambda_* /(2nN)).
\]
The clocks \(w_*,v_*\) are positive algebraic numbers and
\(0<w_*<1\). The bias \(\lambda_*/T_*\) can be irrational.

**Theorem.** At an irrational-bias contact, an accumulating sequence
of exact ordinary collisions requires
\(\log w_*,\log v_*\) to be algebraically independent over
\(\overline{\mathbb Q}\). Thus algebraic dependence of these two
specified logarithms excludes all sufficiently high ordinary
collisions in a neighborhood of the contact, including contacts
with multiple reference roots. This is an all-order reduction:
the [absolute cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md)
and [signed-gap bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md)
give finitely many orders and twists for the fixed core and gaps.
Rational-bias interior contacts are already excluded by the
[branch rank theorem](MULTIVARIATE_OPPOSED_AX_BRANCH_RANK.md).

No conjecture about logarithms is assumed. In particular, Schanuel's
conjecture would imply algebraic independence of the two clocks'
logarithms at every irrational-bias contact, because their
exponentials are algebraic and their ratio is irrational.
The present result excludes the alternative of algebraic
dependence and leaves the independent-clock case open.

## Specializing the two clock parameters in a function field

Use \(u=1/r,\ p=r/L,\ q=r/Y\), where the common count \(r\)
parametrizes the two multiplicities. The limiting values
\(p_*,q_*>0\) are invertible rational linear combinations of
\(T_*,\lambda_*\), so
\[
 F=K(p_*,q_*)=K(\log w_*,\log v_*)
\]
for any number field \(K\). Under the algebraic-dependence
hypothesis, \(F\) has transcendence degree one over \(K\).
Indeed at least one logarithm is transcendental by
Hermite–Lindemann. Enlarge \(K\) to contain the finitely many
algebraic constants occurring in the reference coefficient
arrays. Exponentials of rational linear combinations of
\(p_*,q_*\) are algebraic because they are products of
rational powers of the algebraic clocks.

Center \(p,q\) at \(p_*,q_*\), with increments \(P,Q\).
A summand in a normalized row's coefficient logarithm has
form \((c/u+s)\log(1+\alpha u(p_*+P))\) or the corresponding
expression with \(q_*+Q\), with \(c,s,\alpha\) rational.
After removing its constant in \((u,P,Q)\), the coefficient
of \(u^a\) for \(a\ge1\) has total degree at most \(a+1\le2a\)
in \(p_*,q_*\). The coefficient at \(u=0\) then has no
dependence on the two constants. Exponentiating and taking
any inverse with nonzero algebraic constant term preserves
the bound \(2a\). Thus the angular and physical resultants
have their \(u^a\) coefficients in the spans of monomials
\(p_*^i q_*^j\) with \(i+j\le2a\), with coefficients in
\(K[[P,Q]]\).

The same bound survives the two local eliminations. One way
to see this is to lift all the series to power series over
\(K\) in \(P,Q\) and the six placeholders \(U_{ij}\),
\(i,j\ge0,\ i+j\le2\), before specializing
\(U_{ij}=u p_*^i q_*^j\). A monomial of total clock degree
at most \(2a\) splits into \(a\) monomials of degree at most
two, allowing factors of degree zero. This realizes its
coefficient times \(u^a\) by a monomial of total
\(U\)-degree \(a\). The move to centered logarithmic
coordinates \(W,V\) for the clocks is rational linear in
\(P,Q\). Formal Weierstrass preparation in \(W\), the
norm of the physical resultant in the prepared monic
quotient, and a second preparation in \(V\) are therefore
operations over \(K[[U_{ij}]]\), with \(V\) also retained
as a parameter in the first preparation. At vanishing
\(U_{ij}\) the central angular factor is a power of \(W\)
and the norm is a power of the physical resultant at
\(W=0\), which is regular in \(V\). After specialization
these preparations agree with the analytic ones, by
uniqueness. This is the same norm construction as in the
[one-clock auxiliary argument](MULTIVARIATE_OPPOSED_AUXILIARY_POLYNOMIAL_RATIONAL_BIAS.md#one-transcendental-coefficient-parameter).

View \(F\) as a one-variable function field over its
finite number field of constants. The two elements
\(p_*,q_*\) have a fixed finite set of poles on its
projective curve. Hence the coefficient of \(u^a\) in
each prepared polynomial has pole divisor at most \(a\)
times a fixed effective divisor, since its total degree
in the two elements is at most \(2a\).
By the
[linear pole growth lemma for formal roots](MULTIVARIATE_OPPOSED_AUXILIARY_POLYNOMIAL_RATIONAL_BIAS.md#a-pole-bound-for-formal-root-branches),
after finite ramification the local roots for \(V\), and
then for \(W\) after substituting \(V\), have linear pole
growth over a finite extension of \(F\). The substitutions
and the Laurent inversions giving \(L,Y,r\) preserve that
growth bound. The important field degree here is the
transcendence degree of \(F\), irrespective of the number
of initially named generators.

## Integer evaluation and full rank

Fix a branch and its ramified parameter \(s\). For
completeness, choose \(B\) at least the largest pole order
in \(s\) of the three coordinates \(L,Y,r\). Let
\(P_0(L,Y,r)\) have unknown rational coefficients and
degree at most \(D\). In
\(s^{BD}P_0(L(s),Y(s),r(s))\), the coefficient at
order \(j\) has a pole divisor bounded by \(j+D\)
times a fixed divisor on the one-variable coefficient
field. These coefficients therefore vary over a rational
vector space of dimension \(O(j+D+1)\), by counting
principal parts and constants. Requiring vanishing through
order \((B+1)D-1\) imposes \(O(D^2)\) rational conditions
on \(\binom{D+3}{3}\) unknowns. For large enough \(D\),
a nonzero \(P_0\) exists. Its value on the branch is
\(O(s^D)\). On integer samples that value lies in a
fixed rational lattice, so it is zero at all sufficiently
high samples. An accumulating sequence forces an identity
\(P_0(L,Y,r)=0\) on the branch.

On the same branch the common-factor minors vanish
identically if ordinary collisions accumulate. Their
degree is at least two by even-order strict real positivity.
At irrational bias the factorial characters have full
multiplicative rank on the branch: a dependence would
separate into an equality of nonconstant rational functions
of \(L\) and \(Y\) tending to one. Their first nonzero terms
at infinity would make
\(\lim Y/L=(T_*+e\lambda_*)/(T_*-d\lambda_*)\) algebraic,
whereas Gelfond–Schneider makes it transcendental.
The [Ax rank bound](MULTIVARIATE_OPPOSED_AX_BRANCH_RANK.md#a-bound-from-differential-transcendence)
then requires \(\operatorname{trdeg}_{\mathbb C}
\mathbb C(L,Y,r)\ge3\). This contradicts the auxiliary
identity. The finitely many local branches therefore
contain no accumulating integer samples under the stated
clock dependence.

If the two logarithms are algebraically independent, the
natural spaces of bounded-degree polynomials in them have
quadratic rather than linear dimension in the coefficient
order. The unrefined coefficient count for imposing
\(O(D)\) Taylor conditions is then cubic, the same order
as the unknown monomial count. It supplies no dimension
surplus of the kind used above.
