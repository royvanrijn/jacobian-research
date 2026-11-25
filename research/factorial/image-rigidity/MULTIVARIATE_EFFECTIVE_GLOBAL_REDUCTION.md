# Effective reduction to fully shared near-critical profiles

For two distinct ordinary exponent columns use the moment
polynomials \(Q_j\) of the
[all-degree reduction](MULTIVARIATE_ALL_DEGREE_REDUCTION.md).
For a fully shared profile put
\[
 W=\sum_i\frac{(b_i-a_i)^2}{(a_i+b_i)/2},
 \qquad
 \rho(q)=q(1-q)\sum_i\frac{(b_i-a_i)^2}{a_i+(b_i-a_i)q}.
\]

**Theorem.** There is a computable absolute \(W_0>0\) such
that for every rational \(0<\kappa<1\), one can compute a
terminal \(J(\kappa)\) with the following property. For
\(N\ge J(\kappa)\), any common root of \(Q_{N-1},Q_N\) must
come from a fully shared profile with
\[
               W_0<W<128,\qquad
               \left|\max_{q\in[0,1]}\rho(q)-1\right|<\kappa .
\]
The theorem concerns every ordinary exponent dimension and
degree difference. The constants and cutoffs are unevaluated.

For profiles with private variables on both sides take the
maximum cutoff from the
[effective both-private theorem](MULTIVARIATE_EFFECTIVE_BOTH_PRIVATE_LINEAR.md)
and its mixed-degree consequences. For exactly one private
side use the
[effective regular complement](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_REGULAR_COMPLEMENT.md),
combined with the [effective near-critical endpoint and interior](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_INTERIOR.md)
as specified there. This gives a computable cutoff uniform across
the profiles with any private variable.

For fully shared profiles, \(S\ge W/2\), so \(W\ge128\)
is covered by the
[absolute high-curvature theorem](MULTIVARIATE_EFFECTIVE_ABSOLUTE_CURVATURE_BOUND.md).
The [effective vanishing-curvature theorem](MULTIVARIATE_EFFECTIVE_VANISHING_CURVATURE.md)
covers \(0<W\le W_0\). Choose \(W_0<1\) if necessary.
On \([W_0,128]\) use the
[all-field proof for fixed bands](MULTIVARIATE_EFFECTIVE_SHARED_REAL_BOUNDARIES.md)
with the prescribed \(\kappa\), including real boundaries and
the subcritical simple fold. The maximum of the finite thresholds
proves the assertion.

Within the remaining near-critical band the
[bounded moving quartic theorem](MULTIVARIATE_EFFECTIVE_BOUNDED_SHARED_QUARTIC.md)
already gives computable cutoffs when both scaled quartic fields
and temperatures are confined to prescribed bounded rectangles.
The [effective critical entry](MULTIVARIATE_EFFECTIVE_SHARED_CRITICAL_ENTRY.md)
confines both physical fields near the moving quartic center at a
computable order when the profile is sufficiently close to critical.
The [absolute effective hierarchy](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md)
combines the escaping-field and signed thermal comparisons
with this reduction and the entry estimate.

Every terminal through eight is independently covered by the
[complete terminal-eight theorem](MULTIVARIATE_TERMINAL_EIGHT_ALL_PROFILES.md)
and the earlier fixed-terminal certificates. The computable global cutoff still leaves the intervening fixed
terminals starting at nine to treat. At each such terminal the exact coefficient
and center links must be preserved.
