# Contact equations for changing ordinary factorial cores

A sequence of ordinary two-monomial profiles with bounded curvature
can retain a countably infinite set of rows in the limit. Relative
to this limiting profile, the finite profiles both add new rows and
omit limiting rows. The two operations contribute with opposite signs.
The following tests keep both contributions in the local contact
equations. The omitted curvature tends to zero, but its second-order
factorial correction can be as large as the correction from the
added rows.

## Limits and the two sets of weak rows

Fix adjacent orders \(n,N=n+1\), \(n\ge2\). Pass to a subsequence
of finite ordinary profiles \(A_\ell\) with bounded curvature and
fixed signed total degree gap \(\Delta_0\) on which every unequal
row's multiplicity stabilizes. At a common root the fixed gap can
be obtained on a subsequence by the integer bound in the
[fixed-terminal compactification](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md).
Let \(C\) be the countable multiset of limiting multiplicities.
Use the multiset differences
\[
                   E_\ell^+=A_\ell\setminus C,\qquad
                   E_\ell^-=C\setminus A_\ell .
\]
For a row with signed gap \(d\), midpoint \(m\), write
\(x=d/m,\ w=d^2/m\). Let \(\eta\) denote the supremum of \(|x|\)
in the two differences. For either difference use the sums
\[
 W=\sum w,\quad \lambda=\sum x,\quad K=\sum wx,\quad
 Z=\sum x^2,\quad B=\sum wx^2,\quad L=\sum x^2/d ,
\]
with subscripts \(+\) and \(-\). Both differences eventually contain
only shared rows. In fact at a bounded midpoint there are finitely
many unequal integer rows, and all the multiplicities agree there
for large \(\ell\). At a large midpoint, \(d^2/m\) is bounded by the
total midpoint curvature, so \(|d|/m\) tends to zero uniformly.
Absolute convergence of the limiting profile's curvature gives
\(W_-\to0\). The other sums converge absolutely. After a further
subsequence, \(W_+\to W_*\ge0\) and \(\lambda_+\to\lambda_*\) with
\(|\lambda_*|\le W_*\).

Use the endpoint-normalized core \(P_j^C\) from the compactification
and define
\[
 G_j(z;W)=\sum_{k=0}^j [z^k]P_j^C\,e^{WF_j(k)} z^k,\qquad
 F_j(k)=-k(j-k)/(2j).
\]
For each ordinary row \(a,b\), let \(E(a,b)\) denote its additive
center correction as in that compactification. The logarithm of
the reference center ratio is
\[
       C_*=\Delta_0\log(N/n)+\sum_{(a,b)\in C}E(a,b)
                                     -\frac{\lambda_*}{2nN}.
\]
The signed gaps of the individual differences need not have finite
sums. The center formula uses \(\Delta_0\) for the whole finite
profile and the absolutely convergent corrections for differing
rows. The reference pair in one common coordinate is
\(G_n(z;W_*),G_N(e^{C_*}z;W_*)\).

At a reference common root assume both corresponding roots are
simple. Allow any finite selection of such pairs, using upper
half-plane roots. Adopt the root responses for each selected pair
from the [matched signed jet expansion](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md#root-coefficients-at-an-algebraic-exceptional-clock).
Their definitions use only local derivatives of the specified
polynomials and also apply to the real coefficients of a countable
core. Take upper minus lower differences within each pair. Write
\(a=(\Delta\alpha_i),\, b=(\Delta\kappa_i),\,
g=(\Delta\gamma_i),\, f=(\Delta\sigma_i),\, h=(\Delta\beta_i^H)\)
for these vectors, and \(e=\mathbf1/(2nN)\). A single pair can use
the corresponding complex scalars. Set
\[
 \begin{gathered}
 X=W_+-W_--W_*,\qquad \delta=\lambda_+-\lambda_--\lambda_*,
 \qquad K_0=K_+-K_-,\qquad Z_0=Z_++Z_-,\\
 r=(Z_+-Z_-)/Z_0,\quad u=(B_+-B_-)/Z_0,\quad
 v=K_0^2/Z_0,\quad p=(L_+-L_-)/Z_0 .
 \end{gathered}                                                    \tag{1}
\]
Assume \(Z_0>0\). At \(W_*>0\) this holds automatically for large
\(\ell\). Define \(c=(n^{-2}-N^{-2})/12\).

## A signed region for arbitrary changing profiles

**Theorem.** Every limiting value of \((r,u,v,p)\) belongs to
\[
 \mathcal U=\operatorname{conv}
 \left(\{(1,u,v,p):(u,v,p)\in\mathcal T_{W_*-|\lambda_*|}\}
                 \ \cup\ \{(-1,0,0,-1),(-1,0,0,1)\}\right),       \tag{2}
\]
where \(\mathcal T_M\) is the
[four-vertex tail polytope](MULTIVARIATE_TRANSVERSE_TAIL_RESPONSE_POLYTOPE.md#the-joint-skew-quartic-and-center-polytope).
For the fixed core \(C\) the two vertices at \(r=-1\) may be
restricted to the interval between the limiting lower and upper
bounds for \(-1/d\) among rows of \(C\) at large midpoints.

Suppose the exact finite moments have common roots continuing
the selected reference pairs. If \(a,b,e\) are real-linearly
independent, and \(P\) is projection along their span onto its
real orthogonal complement, then every limiting value satisfies
\[
  v=0,\qquad r P g+u P f=0,\qquad
                   0\le u\le \frac{(1+r)(W_*-|\lambda_*|)}4.     \tag{3}
\]
If \(a,e\) are independent and \(b=\nu a+\mu e\), let \(P\) be
the projection for \(a,e\) and put
\(Q=h-\nu b'+\nu^2 a'/2\), with primes denoting strength
derivatives. Then every limiting value satisfies
\[
                   rP g+uP f+vP Q=0.                            \tag{4}
\]
In (4) the constraint (2) still applies. In particular, if the
projected vectors in (3) or (4) are independent, any such
sequence requires
\[
                 \frac{Z_-}{Z_+}\to1,\qquad
                         \frac{B_+}{Z_+}\to0 .                  \tag{5}
\]
In both cases it also requires \(K_+^2/Z_+\to0\). If \(P g\ne0\) and the
corresponding fixed-core test strictly excludes the \(r=1\)
slice, both \(Z_-/Z_+\to0\) and \(Z_-/Z_+\to\infty\) are impossible,
even along subsequences. Thus the two second-order contributions
must have comparable sizes.

For a single transverse pair there is also a bias interval.
Using the scalar response notation and any fixed real \(s\),
put
\[
 \begin{split}
  \nu&=\Im b/\Im a,\qquad Q=h-\nu b'+\nu^2 a'/2,\\
  \Pi_s(y)&=\Re y-\frac{\Re a+s/(2nN)}{\Im a}\Im y,\qquad
  R_s=\frac{\delta-sX}{2nN Z_0}+\Pi_s(b)\frac{K_0}{Z_0}.
 \end{split}
\]
Then at exact contacts
\[
          R_s=-r\Pi_s(g)-u\Pi_s(f)-v\Pi_s(Q)+cp+O(\eta).
                                                                    \tag{6}
\]
In particular the asymptotic bounds on \(R_s\) come from the
vertices of (2). The finite-core formula is the slice \(r=1\).

### Proof of the region and contact tests

Write \(M_*=W_*-|\lambda_*|\) and \(\theta=Z_+/Z_0\), so \(r=2\theta-1\).
At a positive limit \(Z_+>0\) eventually. Apply the finite-profile
bounds of the four-vertex theorem to \(E_\ell^+\). Up to errors
tending to zero its normalized moments
\((B_+/Z_+,K_+^2/Z_+,L_+/Z_+)\) belong to \(\mathcal T_{M_*}\).
For \(W_*=0\) the same conclusion follows directly from
\(B_+/Z_+\le W_+\), \(K_+^2/Z_+\le W_+^2\), \(|L_+|\le Z_+\),
with an arbitrary choice for the tuple if \(Z_+=0\).
Since \(B_-\le W_- Z_-\) and \(K_\pm^2\le W_\pm^2 Z_\pm\),
the changes in \(u,v\) from omitting the negative terms are
at most \(W_-\) and \(W_-^2+W_+W_-\), respectively. Thus their
limits are \(\theta\) times the corresponding positive normalized
moments. The remaining coordinate is
\(p=\theta L_+/Z_+-(1-\theta)L_-/Z_-\), with any missing summand
set to zero. The quotients \(L_\pm/Z_\pm\) are weighted averages
of \(1/d\). These observations prove (2) and its refinement
for \(C\).

The logarithms of the coefficient multipliers and the center
corrections in the matched jet expansion add for \(E_\ell^+\)
and subtract for \(E_\ell^-\). Their remainders are bounded by
\(O(\eta Z_0)\), including for the infinite sum by absolute
convergence. For example
\(B_\pm\le W_\pm Z_\pm\), and the remainder after the quadratic
term of the midpoint expansion of \(1/(m+td)\) is at most
a constant times \(\eta B_\pm\) when \(|t|\le1/2\).
The other pointwise remainders are bounded by constants times
\(\eta Z_\pm\) as in the matched expansion. The same bounds
give \(Z_0=O(\eta)\), \(K_0=O(\eta)\), \(K_0^2=O(Z_0)\).
The result of moving the Gaussian strength by \(X\) and using
the implicit root expansion at the indicated simple roots is,
for a joint physical contact,
\[
 \begin{split}
 0={}&aX+e\delta+bK_0
     +g rZ_0+f uZ_0+h K_0^2-cpZ_0\mathbf1\\
     &+\tfrac12 a'X^2+b'XK_0+
              O\bigl(\eta Z_0+|X|Z_0+|X|^3+X^2|K_0|\bigr).
 \end{split}                                                       \tag{7}
\]
The sign of the center term follows from subtracting the upper
argument scale in the common physical equation. In the Taylor
expansion, products of \(K_0\) with \(Z_0,B_+-B_-\) and \(L_+-L_-\)
are \(O(\eta Z_0)\), as is \(K_0^3\). Local simplicity bounds the
derivatives of the root map in a fixed coefficient neighborhood.

In the rank-three case a fixed real left inverse for \(a,b,e\)
first gives \(|X|+|\delta|+|K_0|=O(Z_0)\), absorbing the terms
quadratic in \(X\) and \(K_0\) in a small neighborhood. Projection
and division by \(Z_0\) give \(rPg+uPf=O(\eta)\) and \(v=O(Z_0)\).
The \(v=0\) slice of (2) has \(u\le\theta M_*/2\), proving (3).
In the rank-two case a left inverse for \(a,e\) first gives
\(X=O(|K_0|+Z_0)\), then \(X=-\nu K_0+O(Z_0)\). Substitution in
(7) gives \(rPg+uPf+vPQ=O(\eta)\). This proves (4). Independence
forces \(r,u\to0\) in both cases and \(v\to0\) in the second.
The errors in the region calculation then give (5) and the
additional skew condition. The \(r=-1\) slice has \(u=v=0\),
so \(Pg\ne0\) excludes it. A strict fixed-core test excludes
\(r=1\), proving the comparison of sizes. For a single transverse
pair the imaginary equation in (7) gives
\(X=-\nu K_0+O(Z_0)\) with the stated transverse \(\nu\). Project
with \(\Pi_s\) after adding and subtracting \(sX/(2nN)\) in the
bias term. This gives (6).

## Second-order cancellation at every infinite retained core

**Proposition.** For every countably infinite multiset \(C\) of unequal ordinary
rows with finite curvature, every \(W_*>0\), every
\(|\lambda_*|\le W_*\), and every fixed integer \(\Delta_0\),
there are finite ordinary profiles converging to this reference
with
\[
   X,\ \delta,\ K_0,\ Z_+-Z_-,\ B_+-B_-,\ L_+-L_- = o(Z_0).
                                                                    \tag{8}
\]
In particular the origin of (2) is attainable on original profiles
with any such retained core. Both normalized polynomial arrays
and the center ratio approach the reference with error \(o(Z_0)\).
This is a statement about approximation of the arrays.

Choose a distinct row of \(C\) with midpoint \(m\to\infty\), gap
\(d\), and \(z_i=d^2/m^2\). Retain all rows through a sufficiently
large cutoff, omitting one copy of this row. Choose the cutoff so
the curvature of all more distant rows of \(C\) is \(o(z_i)\).
Hence their contributions to each jet statistic are \(o(z_i)\).
The missing collection has \(Z_-=z_i(1+o(1))\).

Replace the omitted row by a row with the same signed gap \(d\)
and midpoint \(m'>m\) whose type is absent from \(C\). To find it,
write \(T(m)=\sum_{m_i\ge m} d_i^2/m_i\to0\). Every occupied
midpoint in \([m,2m]\) of that gap has cost at least \(d^2/(2m)\).
Among the first \(\lfloor 2mT(m)/d^2\rfloor+1\) integer shifts
of \(m\) there must be an unoccupied type. The shifts lie in
\([m,2m]\) for large \(m\), and
\[
 \frac{m'-m}{m}\le \frac{2T(m)}{d^2}+\frac1m,\qquad
 \frac{|d|(m'-m)}m=o(1).
\]
The two rows' \(Z,B,L\) values differ by \(o(z_i)\), since
\(w\) stays bounded and their relative midpoint difference tends
to zero. Their skew values also differ by \(o(z_i)\), using
\(K_i=d z_i\) and the stronger bound with \(|d|\) above.

Let \(w',\lambda'\) be the replacement row's strength and bias.
Set the strength and bias to be supplied by further rows to
\(W_f=W_*+W_--w'\) and \(\lambda_f=\lambda_*+\lambda_--\lambda'\).
Here \(W_f\ge|\lambda_f|\). Indeed replacing the missing row at a
higher midpoint reduces its curvature by a positive amount and
reduces its signed bias by that amount divided by \(d\); all other
missing rows satisfy the cone inequality separately. Thus the
unit-gap strength amounts \((W_f+\lambda_f)/2\) and
\((W_f-\lambda_f)/2\) are nonnegative. Place these amounts to
error \(o(z_i)\) by rounding counts at arbitrarily large
midpoints of gap \(+1\) and \(-1\), respectively, avoiding \(C\)
and the replacement type. There are arbitrarily large unused
heights for each orientation, since otherwise \(C\)'s curvature
would contain a divergent harmonic sum. Their \(Z,|K|,B,|L|\)
can all be made \(o(z_i)\). Finally correct the total degree gap
to \(\Delta_0\) with the required finite number of unit-gap
rows at still higher unused midpoints, of total curvature
\(o(z_i)\). This gives (8) and \(Z_0=2z_i(1+o(1))\). Substitution
in the signed log multiplier and center formulas proves the
coefficient approximation. At a simple contact the same estimate
holds for the two normalized branch displacements in the common
coordinate.

## Integer replacements invisible through arbitrary finite jet order

**Proposition.** There are infinite ordinary retained profiles,
with arbitrarily small total curvature outside any prescribed
finite core, for which the approach in (8) is faster than every
power of \(\eta\) in the normalized coefficient and center arrays.
More precisely, for any specified \(W_*,\lambda_*,\Delta_0\) in
the preceding proposition, one can use such a single fixed \(C\)
and a sequence of finite integer profiles for which, at every
fixed adjacent pair and each fixed \(a>0\), the two arrays and
center ratio approach the reference with error \(o(\eta^a Z_0)\).
One may choose \(C\) independently of \(W_*,\lambda_*,\Delta_0\).
If the prescribed finite core has a simple transverse physical
contact at a positive interior strength and bias, the enlargement
can be made small enough that the countable core has such a
contact as well. The contact can then serve as the reference
for these approximations.

For an integer order \(r\ge3\) and a sufficiently divisible
integer \(D>0\), consider midpoints \(m_i=D/(i+1)\), \(0\le i\le r\).
Use rows with gap \(+2\). Split the signed counts
\((-1)^i\binom ri\) into positive and negative multisets.
The two multisets have the same cardinality and the same sums of reciprocal
midpoint powers of orders \(1,\ldots,r-1\), since the \(r\)-th
forward difference of \((i+1)^k\) vanishes for \(k<r\).
The midpoints are integers, so the row with midpoint \(m_i\)
is \((m_i-1,m_i+1)\) in both multisets.

At fixed gap and sampled moment order, the normalized log
coefficient and center correction per row have finite expansions
in powers of \(1/m\) through every prescribed order, with a
uniform remainder of the next order for all sufficiently large
\(m\). For the log coefficient, use the curvature integral of
the matched jet expansion and expand the trigamma function at
\(m+td\), \(|t|\le1/2\). For the center use its signed integral
over \(m-d/2\) to \(m+d/2\). The expansions follow by integrating
the Taylor formula for \(u/(e^u-1)\) in the Laplace representation
of trigamma. On \(0\le u\le1\) the bounded Taylor remainder is
\(O_r(u^r)\) and on \(u\ge1\) the same bound follows by bounding
the Taylor polynomial and \(0<u/(e^u-1)\le1\). Expanding the
resulting inverse powers of \(m+td\) to the required order has
uniform remainders because \(d=2\). The expansions start with
\(1/m\), using zero endpoints for the chord and the subtraction
of \(d\log(N/n)\) for the center. Consequently the difference
between the two multisets of ordinary rows changes each of these
quantities by \(O_{j,r}(\eta^r)\), where \(\eta=2(r+1)/D\)
for the union of the two blocks.

Build \(C\) by adding the negative blocks, with \(r_\ell=\ell+4\)
and strictly separated bands of heights. At stage \(\ell\) make
the total curvature of all later blocks \(O(\eta_\ell^{r_\ell})\).
This is possible recursively by bounding each new block's
curvature by \(2^{-\ell}\) times all the earlier prescribed
bounds and sending its common multiple \(D_\ell\) to infinity.
The total added curvature can be arbitrarily small. Take the
first \(\ell-1\) negative blocks and the fixed finite core in the
finite profile, omit the \(\ell\)-th negative block, and include
the \(\ell\)-th positive block. The reciprocal moment equalities at
degrees one and two give exact equality of the two blocks'
strengths, biases and \(Z\). Each block has a \(Z\) at least the
largest individual \(Z\) in their union because their sums
coincide. In particular \(Z_0\ge \eta_\ell^2\) for these
differences. The unretained later blocks contribute
\(O(\eta_\ell^{r_\ell})\) to the arrays and center by the
compactification bounds.

Add the two groups of unit-gap rows and the signed gap correction
from the preceding construction, at heights large enough to
supply the desired \(W_*,\lambda_*\) with error
\(O(\eta_\ell^{r_\ell})\) and have all higher jet errors of this
size. None of these new rows belongs to the gap-two blocks of
\(C\), and they avoid the prescribed finite core. Their relative
gaps are smaller than \(\eta_\ell\), as are those of the later
blocks. The signed gap correction is possible for every
\(\Delta_0\). For any fixed order pair the full log arrays and
center now have error \(O_{n,r_\ell}(\eta_\ell^{r_\ell})\).
When choosing \(D_\ell\), take it large enough that the error
bound is at most \(\eta_\ell^{\ell+3}\) for all pairs with
\(N\le\ell\). Later blocks can meet their curvature budgets
independently. This proves \(o(\eta^a Z_0)\) for each fixed
\(a\) and pair.

Finally, at a simple transverse interior contact of the
prescribed finite core, isolate the two branches on a small
strength interval whose endpoints have opposite angular
differences. Their simplicity and those signs persist under
all sufficiently small added curvature in the retained core,
by the coefficient and center bounds. Choose a zero of the
angular difference in the interval, with derivative kept
nonzero by continuity of the simple-root derivatives. Its
unique physical bias obtained from the modulus equation remains
strictly inside the bias cone by making both the interval and
the added curvature sufficiently small. Here use the same
integer total gap of the reference before and after adding the
blocks, and the absolutely convergent center corrections. The
approximation construction therefore applies at the persisting
contact.

The signed equations force comparable corrections when the
fixed-core tests are strict at both ends. Exact integer
replacement can cancel those corrections through arbitrary
finite order in the coefficient arrays. A uniform exclusion
that allows changing cores must address these replacements
when this cancellation occurs.
