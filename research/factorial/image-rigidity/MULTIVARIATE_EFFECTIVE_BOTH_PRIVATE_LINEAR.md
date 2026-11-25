# Effective consecutive cutoff for two private linear rows

Suppose each of two distinct ordinary monomials has one private
variable of exponent one. All other unequal rows are shared.
Use \(W=\sum_i d_i^2/m_i\) and \(\Delta=\sum_i d_i\) for the shared
rows, with \(d_i=b_i-a_i\), \(m_i=(a_i+b_i)/2\).
The two private factorials cancel the binomial coefficient, as in
the [both-private comparison](MULTIVARIATE_BOTH_PRIVATE_UNIFORM_EVENTUAL.md).

**Theorem.** A finite calculation gives an absolute terminal threshold
for consecutive coprimality for all such profiles, uniformly in
dimension, shared exponents, and the integer \(\Delta\).
The threshold is unevaluated.

The [absolute curvature bound](MULTIVARIATE_EFFECTIVE_ABSOLUTE_CURVATURE_BOUND.md)
covers \(W\ge128\), since the ordinary shared curvature is at least
\(W/2\). For \(W=0\) the moments reduce to consecutive geometric
sums. In the remaining range use the endpoint-normalized moment
\[
 P_j(z)=\sum_{k=0}^j c(j,k)z^k,\qquad
 c(j,k)=e^{L_j(k/j)},\qquad z_j=e^{\mu_j}t .
\]
Each positive real normalization preserves the argument, and at a
common root \(z_N/z_n=e^{\mu_N-\mu_n}\), with \(N=n+1\).

## Computable real coefficient bounds

For \(W\le128\), each \(|r_i|=|d_i|/(2m_i)\) is at most \(64/65\).
Indeed if \(|r_i|\ge1/2\), then \(2|d_i r_i|\le128\) implies
\(|d_i|\le128\), while the smaller positive exponent gives
\(|r_i|\le|d_i|/(2+|d_i|)\). For \(|r_i|<1/2\) the bound is
immediate. Thus \(x_i(u)/m_i=1+r_i(2u-1)\in[1/65,2]\)
throughout \(0\le u\le1\). Use the [varying-profile phase estimates](MULTIVARIATE_VARYING_PROFILE_PHASE.md)
to write
\[
 L_x=xG+H+R_x,\qquad
 \mu_N-\mu_n=\Delta\log(N/n)+O(W/n^2).
\]
Here \(G,H,R_x\) vanish at \(0,1\), their second derivatives
have bounds \(O(W),O(W),O(W/x)\), and
\(\partial_x R_x''=O(W/x^2)\). These constants are computable:
the factors \(1/(x_i/m_i)\) are bounded by \(65\), and
the real Stirling remainder and its derivatives are bounded
by differentiating the positive-argument Binet integral
or the trigamma integral used in the varying-profile proof.
The summed reciprocal factors satisfy \(\sum_i1/m_i\le W\).
The endpoint conditions bound the functions and their first
derivatives by the same second-derivative bounds.
In particular this also gives computable constants in (2),
(4), (5), and (6) of the both-private comparison.
In (4) use \(|r_i|\le W/2\) and expand the single rational
factor \(1/(1+r_i(2u-1))\). All constants can be used for
\(n\) and \(N\) in the same profile.

## The two geometric-tail ranges

For a rational \(w_0>0\) and \(W\ge w_0\), choose a rational
\(R>1\) so \(R e^{-w_0/8}<1\). Split \(P_j\) into its lower
half and its reversed upper half. Their geometric comparison
parameters are \(\rho_L=e^{G'(0)}\), \(\rho_R=e^{-G'(1)}\).
Each is at most \(e^{-w_0/4}\) because \(G''\ge W/2\).
For each fixed index \(k\) in a half-polynomial, the log
coefficient differs from \(k\log\rho_L\) or \(k\log\rho_R\)
by at most \(C(k+k^2)/j\), by the endpoint-subtracted estimates.
For coefficients beyond a chosen index, the bounds (2)
of the both-private comparison give geometric tails on \(|z|\le R\),
including their first two logarithmic derivatives before taking
the logarithms. Thus for any prescribed tolerance the two
halves approximate \((1-\rho_L z)^{-1}\), \((1-\rho_R z)^{-1}\)
and their derivatives to that tolerance above a computable order.
To compute it, bound a common geometric tail first, then use the
finite coefficient error and require that each tested index is
in both orders' relevant halves.

The comparison functions are bounded above and away from zero
on the disk. For accurate enough halves, the usual two endpoint
Rouché inequalities give \(R^{-1}<|z_j|<R\). At a root the
half-polynomial modulus ratio gives \(|\log|z_j||\le C_R/j\).
Approximate the angular derivatives of each half by those of
the corresponding geometric function, evaluated at the same
radius. The geometric functions are unchanged within the pair.
Their derivatives of the logarithmic slopes are bounded on a
slightly smaller disk, and the two root radii differ by
\(O(1/n)\). Choosing the initial tolerance small and then
increasing the order gives the strict phase derivative bound
(7) of the both-private proof, which excludes common roots.
Every inequality here has a computable threshold.

For \(W\) small and \(\mathcal L=nW\) large, use the cutoff
halves of the
[small-curvature comparison](MULTIVARIATE_BOTH_PRIVATE_UNIFORM_EVENTUAL.md#small-curvature-with-nw-unbounded).
A fixed twice differentiable cutoff on the transition interval
can use \(1-10s^3+15s^4-6s^5\), \(s=2u-1/2\), with the stated
constant values on each side. Use \(e^{bW}\) for the disk radius,
taking \(b>0\) strictly below the constant in (6) and below
\(1/4\). The proof of (8) sums geometric series of polynomials
in \(k\), so its constants are bounded effectively using
\(1-e^{-a W}\ge aW/(1+aW)\). The mixed order and index bounds
use the derivatives of the chosen cutoff and the real
coefficient estimates just described.

Once \(\mathcal L\) exceeds a computable constant, these
bounds certify the nonvanishing tilded halves, the two
annulus inequalities, and, at any root,
\[
             j|\log|z_j||\le C\bigl((jW)^{-1}+W\bigr).
\]
For the last step absorb \(C|\log|z_j||/W\) into the left
side by increasing the lower bound for \(jW\). Choose \(w_0\)
small and \(\mathcal L_0\) large so that the right side is
less than \(1/8\) whenever \(W\le w_0\) and \(nW\ge\mathcal L_0\).
A common root in that regime forces \(\Delta=0\) by the
linked center-drift formula and integrality, for every
sufficiently large computable \(n\). Specifically the two
radial bounds give \(n|\mu_N-\mu_n|<1/4\), whereas
\(\Delta\ne0\) gives \(n/N-O(W/n)>1/4\).
With \(\Delta=0\) the real-order bounds in (8) give a
difference of angular derivatives bounded by
\(C'\bigl(\mathcal L^{-2}+(n\mathcal L)^{-1}\bigr)\).
Increasing \(\mathcal L_0\) and the terminal makes it less
than one, again giving (7). The preliminary small \(W\)
cap for all these estimates is absolute and can be chosen
before \(w_0,\mathcal L_0\).

## Bounded \(nW\)

It remains to treat \(\mathcal L=nW\le\mathcal L_0\); at the
upper order \(NW\le2\mathcal L_0\). Put
\[
  D_j(z)=\sum_{k=1}^j [c(j,k)-c(j,k-1)]z^k,\qquad
  (1-z)P_j(z)=1-z^{j+1}+D_j(z).
\]
On this range the first and second coefficient differences
are bounded by \(C/j\) and \(C/j^2\) by the real coefficient
bounds. For \(|z|\le e^{-B/j}\), sum the first differences
geometrically and choose \(B\) large enough so
\(|z|^{j+1}+|D_j(z)|<1\). Reversal gives the computable
annulus \(|\log|z_j||\le B/j\) for every root. A further
first difference in the expression for \(D_j\) gives
\(|(1-z_j)D_j(z_j)|\le C_1/j\) on this annulus.
Choose \(H\) large enough that \(j|1-z_j|\ge H\) implies
\(|D_j(z_j)|\le C_1/H\) and then
\(j|\log|z_j||\le1/8\).

If \(j|1-z_j|\le H\), take \(h=j\log z_j\) using the
branch near one, so \(|h|\le2H\) for \(j>2H\).
On any fixed bounded disk in \(h\), the trapezoidal estimate
and the uniform small \(W\) expansion give
\[
  j^{-1}P_j(e^{h/j})=F_\Lambda(h)+O(1/j),\qquad
  F_\Lambda(h)=\int_0^1e^{-\Lambda u(1-u)/2+hu}\,du,\quad
  \Lambda=jW .
\]
The existing [model proof](MULTIVARIATE_BOTH_PRIVATE_UNIFORM_EVENTUAL.md#small-curvature-with-bounded-nw)
shows that every \(F_\Lambda\) zero is imaginary for
\(\Lambda\ge0\), and that \(F_\Lambda\) and
\(\mathcal D F_\Lambda=(h\partial_h+\Lambda\partial_\Lambda)F_\Lambda\)
have no common zero. For a concrete compact residual, subdivide
\(\Lambda\in[0,2\mathcal L_0]\), \(|h|\le2H\), \(|\Re h|\ge1/8\)
and interval-evaluate \(|F_\Lambda(h)|\).
The integral and its derivatives have computable bounds on
all these rectangles. Strict positivity and compactness
guarantee termination with a positive lower bound.
Increasing \(j\) so the trapezoidal error is smaller forces
\(j|\log|z_j||\le1/8\) also in this range.
The same center-drift argument now gives \(\Delta=0\)
for any common root.

Use this equality and the real-order estimate (5) of the
both-private proof to bound
\[
       |D_N(z_N)-D_n(z_n)|\le C_2/n .
\]
Indeed each of the \(n\) common coefficient differences
changes by \(O(n^{-2})\), the extra difference is \(O(n^{-1})\),
and moving radially by \(|\log(z_N/z_n)|=O(n^{-3})\)
costs at most \(O(n^{-2})\) using the annular bounds.
Subtracting the two root equations gives
\((z_n-1)z_n^{n+1}=O(1/n)\) since the same radial change
in the terminal power costs \(O(n^{-2})\).
The annulus bounds the factor \(|z_n^{n+1}|\) below.
Consequently \(n|1-z_n|\) is bounded by a computable
constant at a common root, with both scaled fields
\(h=n\log z_n,\ h'=N\log z_N=h+h/n+O(n^{-2})\)
in a fixed computable disk.

The trapezoidal comparison at both orders is now linked.
With \(\Lambda=nW\), the integrand exponent changes from
\(L_n(u)+hu\) to \(L_N(u)+h'u\) by
\[
   n^{-1}\bigl(-\Lambda u(1-u)/2+hu\bigr)+O(n^{-2}).
\]
The explicit trapezoidal endpoint terms change by \(O(n^{-2})\),
and the uniform second derivative controls both quadrature errors.
Thus at a common root both \(F_\Lambda(h)\) and
\(\mathcal D F_\Lambda(h)\) have modulus at most \(C_3/n\).
Interval subdivision on the enlarged compact disk certifies
a positive minimum of
\(|F_\Lambda(h)|+|\mathcal D F_\Lambda(h)|\) by the model
proof. The bound \(C_3/n\) is eventually smaller, completing
this range. The maximum of the computed bounds in the
three ranges gives the absolute terminal threshold.
