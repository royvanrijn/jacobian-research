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

## The remaining large-order regime

**Corollary.** Every sequence of common roots with terminals tending
to infinity has a subsequence with one private linear variable,
bounded curvature, and shared endpoint curvature \(\sigma\to1\).
Both fields tend to the positive endpoint. The scaled temperatures
and fields in the moving cubic model are bounded. Any limiting
temperature belongs to one of these two finite sets.

- A temperature at which two distinct model zeros have the same
  imaginary coordinate. In this case the subsequence has nonzero
  \(\Delta\), and both \(|\Delta|/n^{1/3}\) and its reciprocal are
  bounded.
- A temperature with a simple zero \(H\) and
  \(2H-T H_T=0\). In this case the subsequence has \(\Delta=0\).

To prove the reduction, first use the preceding support-type
corollary. The
[uniform one-private comparison](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md)
excludes private total degree at least two and all strictly
supercritical limits for private degree one. It localizes critical
limits to the positive endpoint.
The [finite subcritical activation](MULTIVARIATE_UNIFORM_SUBCRITICAL_PRIVATE_ACTIVATION.md)
does the same at strictly subcritical limits, including zero
shared curvature. The
[uniform endpoint comparison](MULTIVARIATE_PRIVATE_ENDPOINT_UNIFORM_THERMAL_TAILS.md)
then excludes all strictly subcritical limits and every escaping
temperature or scaled field at critical curvature.

In the remaining bounded window, the model fields have the same
limiting imaginary coordinate. For nonzero \(\Delta\), their scaled
real drift is \(\ell(n+1)^{2/3}\Delta\log(1+1/n)\), up to an \(O(1/n)\)
change due to rescaling the first field. If both model fields had the
same limit, the
[double-zero and simple-branch comparison](MULTIVARIATE_PRIVATE_CUBIC_DOUBLE_LIMITS.md)
would exclude the pair. If \(\Delta=0\), the two model limits coincide.
The same comparison forces a simple zero with vanishing order
derivative. The [zero-separation theorem](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md)
and its [negative-temperature extension](MULTIVARIATE_HALF_LINE_CUBIC_NEGATIVE_TEMPERATURE.md)
show that both exceptional temperature sets are finite. Their
membership remains unevaluated outside the interval \([-2,2]\)
excluded by the [extended interval certificate](MULTIVARIATE_HALF_LINE_CUBIC_EXTENDED_TEMPERATURES.md).
The [first linked corrections](MULTIVARIATE_PRIVATE_CUBIC_LINKED_SIMPLE_LIMITS.md)
give additional necessary conditions at simple candidates. They
constrain the integer degree drift at a same-height crossing and
the possible first-order row coefficients at a speed zero.

The curvature rescaling by itself does not discretize the limiting
temperature. For example, append a private linear row and choose
shared rows \((A,A+1)\), \((A,A-1)\) with \(A=n^2\), in a total
number \(m=A(1-2^{-2/3}T_0 n^{-1/3})+O(1)\). The number \(m\)
can have either required parity. Choose the orientation counts so
that \(\Delta-1\) is their signed difference. Then any prescribed
finite \(T_0\) and any finite limit of \(\Delta/n^{1/3}\) are
attainable. Indeed \(\sigma=m/A\) and \(c_2\to1/2\), which
gives \(T\to T_0\). This is a construction of limiting profile
data.

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
The [joint curvature measure](MULTIVARIATE_JOINT_FACTORIAL_CURVATURE.md)
uses an exact factorial density for each primitive row curve. It
couples all orders and determines the relative center. Relaxing only
these densities to arbitrary positive measures on the same curves
already gives a collision at orders two and three.
The [uniform curvature closure](MULTIVARIATE_JOINT_CURVATURE_PROFILE_CLOSURE.md)
identifies the actual diffuse part as a diagonal density \(r\,dr\).
Its error is uniform in the two curvature indices, and the relative
centers have one additional signed diffuse parameter.
At terminal six, [one-private-side separation](MULTIVARIATE_TERMINAL_SIX_ONE_PRIVATE_SIDE.md)
and the [two-private linear certificate](MULTIVARIATE_TERMINAL_SIX_TWO_PRIVATE_LINEAR.md),
together with the [other private-degree criteria](MULTIVARIATE_TERMINAL_SIX_PRIVATE_ANGLES.md),
cover every profile with any private variable. Thus the remaining
terminal-six family is fully shared. The same linked angular theorem
covers \(S\ge6\), and the
[effective relative-gap bands](MULTIVARIATE_TERMINAL_SIX_DIFFUSE_BANDS.md)
give explicit angular separation for \(5/4\le S\le6\) under the
respective relative-gap bounds. The
[near-flat certificate](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT.md)
covers the order-five outer mean \(u\ge24/25\). In particular
it covers \(S\le1/11\). To see this, if \(\eta=\max |b-a|/r\)
then \(\eta\le S\) by integrality, and the coefficient bound in
the diffuse-band note gives
\(u\ge \exp(-2S/[5(1-\eta)])\ge1-2S/[5(1-S)]\ge24/25\).
Arbitrary fully shared profiles between these proved ranges remain open.

A prescribed finite set of primes and precisions alone cannot give
a universal exclusion. The
[simultaneous prime-power construction](MULTIVARIATE_PADIC_DIGIT_AUTOMATON.md#a-simultaneous-obstruction-at-finitely-many-primes)
produces fully shared profiles with both normalized moments congruent
to their binomial polynomials at all the prescribed moduli, within
any open positive curvature interval. In particular, the ordinary
Newton and residual screens at those moduli allow both candidate
common-factor degrees at terminal six.
The [marked-row Eisenstein criterion](MULTIVARIATE_MARKED_ROW_VARIABLE_PRIME.md)
gives examples in every degree using a prime selected from the
exponent data. It can also extend any fixed profile by two balanced
shared rows with arbitrarily small added curvature, preserving
private degrees and total-degree difference.

At each terminal, coprimality over \(\mathbb Q\) gives the assertion
over any characteristic-zero field by the rational Bézout identity.
Both pure monomial coefficients are nonzero at every order, so
the univariate reduction in \(t\) covers all projective coefficient
choices. The two moments at orders one and two already lie in the
low-terminal theorem. An all-degree proof requires both eventual
exclusion and coverage of every terminal it leaves below its cutoff,
or an argument that works directly at every terminal.

## The discrete row-merger route

For proportional rows the [joint density comparison](MULTIVARIATE_JOINT_FACTORIAL_CURVATURE.md#inflation-and-identification-along-a-primitive-slope)
orders the effect of integer inflation and identification on all
coefficient curvatures. For nonproportional shared rows the
[exact merger sign law](MULTIVARIATE_ROW_MERGER_CURVATURE_SIGN.md)
gives an explicit flattening threshold and a decrease of the leading
stationary curvature. These comparisons leave the corresponding
root transfer open.

There is an [all-order physical binary transfer](SFC_PHYSICAL_BINARY_RESIDUAL_PRIVATE.md)
when untouched independent rows leave at least two units of private
exponent on each side. The proof applies to every positive correlated
base with finite mixed moments, by Hölder's inequality and separation
of the root arguments.

A proof for arbitrary support size cannot rely on zero-freeness
throughout the physical Gaussian interpolation. An
[ordinary three-monomial example](SFC_PHYSICAL_ROW_MERGE_INTERIOR_COLLISION.md)
has regular endpoints for the window of orders two through four
and a common projective zero at a physical squared correlation
strictly between \(21/25\) and \(17/20\). The discrete endpoint
transfer and the general binary interpolation remain open.
