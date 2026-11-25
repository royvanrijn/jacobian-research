# A reduction for all-degree multivariate Strong Factorial

Let \(a,b\) be distinct nonnegative integer exponent columns and
remove rows on which they agree. Throughout, the factorial moments
are the ordinary ones,
\[
 Q_j(t)=\sum_{k=0}^j \binom jk
            \prod_i ((j-k)a_i+kb_i)!\,t^k.
\]
Write \(S=\sum_i(b_i-a_i)^2/\max(a_i,b_i)\).
A row is private when one entry is zero. A profile is fully shared
when there are no private rows.

**Theorem.** For each fixed \(s_0>0\) there is a terminal
\(N_0(s_0)\) such that every fully shared profile with \(S\ge s_0\)
has \(\gcd(Q_{N-1},Q_N)=1\) for \(N\ge N_0(s_0)\).
This is uniform in the dimension, the exponents, and the difference
of the total degrees. The threshold is unevaluated.

For fully shared profiles the midpoint curvature \(W\) satisfies
\(S\le W\le2S\). Suppose an unbounded sequence of terminals had
common roots. The [unbounded-curvature theorem](MULTIVARIATE_UNBOUNDED_CURVATURE_EVENTUAL.md)
allows us to pass to a subsequence with \(W\) in a compact positive
band. Put \(\delta=\max\rho-1\), where
\(\rho(q)=q(1-q)\sum_i(b_i-a_i)^2/(a_i+(b_i-a_i)q)\).
The [strictly subcritical](MULTIVARIATE_UNIFORM_SUBCRITICAL_BANDS.md)
and [strictly supercritical](MULTIVARIATE_UNIFORM_SUPERCRITICAL_SUPPORTS.md)
band theorems force \(\delta\to0\).
The [critical regular comparison](MULTIVARIATE_CRITICAL_REGULAR_COMPLEMENT.md)
then forces both centered fields into the quartic neighborhood.
Bounded \(\sqrt n|\delta|\) is covered by the
[moving quartic theorem](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md).
Otherwise pass to a subsequence with \(\sqrt n|\delta|\to\infty\)
and fixed sign of \(\delta\). The
[supercritical thermal proof](MULTIVARIATE_SHARED_SUPERCRITICAL_THERMAL_APPROACH.md)
excludes the positive sign. The
[subcritical thermal comparison](MULTIVARIATE_SHARED_SUBCRITICAL_THERMAL_REGULAR_FIELDS.md)
places both fields at the same limiting scaled fold, which the
[Airy branch comparison](MULTIVARIATE_SHARED_SUBCRITICAL_THERMAL_FOLD.md)
excludes. This exhausts the sequences and proves the theorem.

**Corollary.** There is one absolute eventual terminal cutoff for
fully shared ordinary profiles. The proof also gives one absolute
eventual cutoff on the union of fully shared profiles and profiles
with private variables on both sides. If degree-unbounded
counterexamples exist, they have a subsequence with bounded \(S\)
and with private variables on exactly one side.

Indeed the [vanishing-curvature theorem](MULTIVARIATE_VANISHING_CURVATURE_EVENTUAL.md)
gives one cutoff throughout a fixed neighborhood of \(W=0\).
Since \(S\le W\le2S\) for shared profiles, this and the theorem
above yield one cutoff on the full shared family. The
[both-private theorem](MULTIVARIATE_BOTH_PRIVATE_UNIFORM_EVENTUAL.md)
handles the second family. Together with the unbounded-curvature
theorem this leaves exactly one private side and bounded \(S\).
There the sum of the private exponents is a bounded positive
integer, since every private row contributes its exponent to \(S\).
The [bounded-private-degree annulus](MULTIVARIATE_BOUNDED_PRIVATE_ROOT_ANNULUS.md)
gives a uniform compact range for endpoint-normalized roots
throughout this family.
For one private linear variable, the
[moving cubic model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md)
gives the exact linked interpolation near the critical endpoint.

## The fixed-terminal obligation

An absolute eventual theorem would leave finitely many terminal
orders starting at six. Each such order still includes arbitrary
dimension and exponents. The
[fixed-terminal compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md)
describes the full closure of the normalized coefficient pairs and
the relative center. A certified separation of that family gives
a finite algebraic certificate for the terminal. The closure includes
a signed center parameter as well as a diffuse variance. A separation
test that varies the two moment arrays independently admits spurious
collisions already at six, as shown by the
[coefficient-link calculation](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md).
A proof at each remaining terminal must retain the links between the
two orders or use the exact arithmetic of the factorial coefficients.
A prescribed finite set of primes and precisions alone cannot give
a universal exclusion. The
[simultaneous prime-power construction](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md#a-simultaneous-obstruction-at-finitely-many-primes)
produces fully shared profiles with both normalized moments congruent
to their binomial polynomials at all the prescribed moduli, within
any open positive curvature interval. In particular, the ordinary
Newton and residual screens at those moduli allow both candidate
common-factor degrees at terminal six.


At each terminal, coprimality over \(\mathbb Q\) gives the assertion
over any characteristic-zero field by the rational Bézout identity.
Both pure monomial coefficients are nonzero at every order, so
the univariate reduction in \(t\) covers all projective coefficient
choices. The two moments at orders one and two already lie in the
low-terminal theorem. An all-degree proof requires both eventual
exclusion and coverage of every terminal it leaves below its cutoff,
or an argument that works directly at every terminal.
