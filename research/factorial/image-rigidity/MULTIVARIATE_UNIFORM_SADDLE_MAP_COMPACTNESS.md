# Uniform compactness of stationary fibers for varying ordinary profiles

Let \(a,b\) be distinct nonnegative integer exponent vectors. Remove equal
rows and write
\[
 d_i=b_i-a_i,\qquad m_i=(a_i+b_i)/2,\qquad
 r_i=d_i/(2m_i),\qquad w_i=d_i^2/m_i=2|d_i||r_i|,\qquad W=\sum_i w_i.
\]
Allow the dimension and all exponents to vary. Let \(\mathbb H\) be the
upper half-plane and \(\mathcal T=\{-\pi<\Im x<0\}\). The centered
stationary map, with principal logarithms on \(\mathbb H\), is
\[
 X_c(q)=S(q)+\log(1-q)-\log q,\qquad
 S(q)=\sum_i d_i\log(1+r_i(2q-1)).
 \tag{1}
\]
This subtracts the real constant \(\sum_i d_i\log m_i\) from the usual
stationary map.

**Theorem.** Fix \(0<W_-\le W_+<\infty\) and a compact \(K\subset\mathcal T\).
There exists a compact \(C\subset\mathbb H\) such that every ordinary
profile with \(W\in[W_-,W_+]\) satisfies \(X_c^{-1}(K)\subset C\).
The same compact containment holds for locally uniform limits of such
profiles. Every such limiting map has exactly
\[
 {\bf1}_{\{P=0\}}+{\bf1}_{\{Q=0\}}
 \tag{2}
\]
preimages, counted with multiplicity, over each \(x\in\mathcal T\), where
\(P=\sum_{a_i=0}b_i\) and \(Q=\sum_{b_i=0}a_i\) are the private degrees
of the retained rows. In particular, the varying shared-support closure
preserves the two-sheet count.

## Escape to infinity

Every contribution to \(\Im S(q)\) is nonnegative, since \(d_i\) and
\(\Im(1+r_i(2q-1))\) have the same sign. Also
\(0<\arg q-\arg(1-q)<\pi\), as \(q/(1-q)\) lies in \(\mathbb H\).
If \(X_c(q)\in K\), compactness gives an \(\epsilon>0\) such that
\[
 \Im S(q)\le \pi-\epsilon .
 \tag{3}
\]
Suppose there is a sequence in these fibers with \(|q|\to\infty\).
The real part of \(S\) stays bounded, since
\(\log|(1-q)/q|\to0\). Put \(z=2q-1\) and pass to a subsequence on which
\(\arg z\) converges.

First suppose the limiting angle is in \((0,\pi)\), so that
\(\eta\le\arg z\le\pi-\eta\) for some \(\eta>0\).
Choose \(L>2\) with
\(2\arcsin(1/L)<\epsilon\) and \(\arcsin(1/L)<\eta/2\).
For \(|r_i z|\ge L\), the circular angular difference between
\(r_i z\) and \(1+r_i z\) is at most \(\arcsin(1/L)\), by the tangent
from zero to the unit disk centered at \(r_i z\). The angles also
differ by this amount on the principal branch because they both have
the same imaginary sign and lie on the same angular half-plane. One
positive and one negative \(r_i\) of this size would contribute at
least \(\pi-2\arcsin(1/L)\) to \(\Im S\), using \(|d_i|\ge1\), which
contradicts (3). Thus all such rows have one sign and their total
\(\sum |d_i|\) is bounded by \(2\pi/\eta\), since each corresponding
logarithm has absolute argument at least \(\eta/2\).

For \(|r_i z|<L\), write \(v=r_i z\). The values \(v\) lie in closed
sectors separated from the real directions, with the origin allowed.
Uniformly on this compact set, \(|\log(1+v)|\le A|v|\) for some
finite \(A\), by continuity of \(\log(1+v)/v\) at zero. Moreover
\[
 |\arg(1+v)|\ge |\sin\arg(1+v)|
 =\frac{|\Im v|}{|1+v|}
 \ge\frac{\sin\eta}{1+L}|v|.
\]
By (3), the sum of \(|d_i||r_i z|\) over these smaller factors is
bounded. Their total \(w_i\)-mass is therefore \(O(1/|z|)\), and
their contribution to \(\Re S\) is bounded. For all sufficiently
large \(|z|\), the remaining mass is at least \(W_-/2\).
The bound on \(\sum |d_i|\) in that group yields at least one
\(|r_i|\ge W_- \eta/(8\pi)>0\). Its real logarithm tends to
\(+\infty\), and all other real logarithms in that same group are
positive by \(|1+r_i z|\ge L-1>1\). Their signed sum diverges with
their common sign, contradicting boundedness of \(\Re S\).

If \(\arg z\to0\), write \(z=x(1+i\kappa)\) with \(x\to+\infty\) and
\(\kappa\to0^+\). For a negative row put \(B=|r_i|x\). If \(B\ge3/2\),
then \((B-1)/(\kappa B)\ge1/(3\kappa)\), so the argument of
\(1-B-i\kappa B\) tends to \(-\pi\) uniformly for all such rows.
Condition (3) and \(|d_i|\ge1\) exclude every \(B\ge3/2\) for
sufficiently small \(\kappa\). For \(B<3/2\) and \(\kappa\le1/2\),
\[
 1-|1-B-i\kappa B|^2=B[2-(1+\kappa^2)B]\ge B/8.
\]
Thus the real contribution of the negative row is at least
\(|d_i|B/16=w_i x/32\), by \(-\log t\ge 1-t\) for \(t>0\)
applied to the squared modulus. For a positive row, concavity of
\(\log(1+rx)\) in \(0\le r\le1\) gives the real contribution at least
\(d_i r_i\log(1+x)=w_i\log(1+x)/2\). Since \(W\ge W_-\), the combined
lower bound tends to \(+\infty\). The case \(\arg z\to\pi\) reduces to
this one by replacing \((z,r_i,d_i)\) with
\((-\overline z,-r_i,-d_i)\), which sends \(S\) to
\(-\overline S\). This rules out every escape to infinity.

## Approach to the finite real boundary

Suppose a sequence of fiber points tends to \(q_0\in\mathbb R\).
Fix a small cutoff \(\rho>0\) so that
\(|r(2q-1)|<1/2\) whenever \(|r|<\rho\) and \(q\) belongs to a fixed
disk around \(q_0\). The kernel
\[
 k_r(q)=\frac{\log(1+r(2q-1))}{2r},\qquad k_0(q)=q-\tfrac12
 \tag{4}
\]
and its \(q\)-derivative are uniformly bounded on the disk for those
rows. The kernel is real on its real diameter. Consequently
\(\sum_{|r_i|<\rho}w_i k_{r_i}(q)\) has bounded real part and
imaginary part \(O(\Im q)\), uniformly over the profiles.

There are finitely many possible row types with \(|r_i|\ge\rho\):
\(|d_i|\le W_+/(2\rho)\), \(m_i\le |d_i|/(2\rho)\), and \(m_i\)
is half-integral. Each multiplicity is also bounded. Pass to a
subsequence of fixed multiplicities for those types. If
\(q_0\notin\{0,1\}\) is a logarithmic singularity of these rows,
all its singular coefficients \(d_i\) have the same sign. Indeed, a
row pole is \(q=(1-1/r_i)/2\), to the left of zero for \(0<r_i<1\)
and to the right of one for \(-1<r_i<0\). The real logarithmic
divergence cannot cancel. At any other \(q_0\notin\{0,1\}\), all
finite row arguments tend to zero or signed \(\pi\), and
\(\Im[\log(1-q)-\log q]\) tends to zero or \(-\pi\).
The limiting height of \(X_c\) is then an integer multiple of \(\pi\),
outside \(K\).

At \(q_0=0\) the only row factors singular there have \(r_i=1\),
and at \(q_0=1\) they have \(r_i=-1\). Exactly,
\[
 X_c(q)=S_{\rm sh}(q)+(P-1)\log q+(1-Q)\log(1-q)+(P-Q)\log 2,
\]
with \(S_{\rm sh}\) the shared-row part of \(S\). If the pertinent
logarithm has nonzero coefficient, the same real divergence
excludes the endpoint. Otherwise the limiting height is zero, since
all nonsingular row factors have positive real boundary value at
that endpoint and the small-factor tail has imaginary part tending
to zero. Every possible finite real limit is excluded. Together
with the bound at infinity, this proves the uniform compactness
claim.

## Retained and diffuse limits

Represent a profile by the positive measure
\(\nu=\sum_i w_i\delta_{r_i}\) on \([-1,1]\). At each cutoff
\(|r|\ge\rho\), the finite-type argument above gives a subsequence
of constant multiplicities. Take nested cutoffs down to zero and a
diagonal subsequence, and pass to a limit of the total masses.
Uniform continuity of a test function at \(r=0\) then gives weak
convergence to a measure of the form
\[
 \nu_*=\sum_{\rm retained\ rows} w_i\delta_{r_i}+\tau\delta_0,\qquad
 \tau\ge0.
\]
The converse approximation for any such summable list and diffuse
mass is by longer finite truncations, with the diffuse mass
approximated by copies of \((A,A+1)\) as \(A\to\infty\). The limiting
total mass is the sum displayed above, and the finite masses lie in
a fixed positive bounded band whenever that sum is positive.
This is the same extraction as for the
[coefficient phases](MULTIVARIATE_VARYING_PROFILE_PHASE.md), now
including \(r=\pm1\).

On compact subsets of \(\mathbb H\), the kernel (4) is continuous
jointly in \(q\) and \(r\in[-1,1]\). Therefore the stationary map
converges locally uniformly to
\[
 X_{\nu_*}(q)=\int k_r(q)\,d\nu_*(r)+\log(1-q)-\log q.
 \tag{5}
\]
For example, uniformity follows by approximating the compact family
\(q\mapsto k_\cdot(q)\) in uniform norm by a finite net in
\(C([-1,1])\), and using the common mass bound. The diffuse term
in (5) is \(\tau(q-\tfrac12)\). The same convergence holds for
derivatives by local analyticity. Private row multiplicities are
eventually fixed because all their atoms have \(|r|=1\).

Fix a target \(x\in K\). The function \(X_{\nu_*}-x\) cannot vanish
identically, since \(X_{\nu_*}\) extends across \(1/2\) and takes
value zero there. A zero outside the common compact \(C\) would,
by local uniform convergence and Rouché's theorem on a small disk
disjoint from \(C\), force nearby finite-profile fibers outside \(C\).
All limit and finite-profile zeros therefore lie in \(C\). Choose
a closed rectangle in \(\mathbb H\) containing \(C\) in its interior.
Rouché's theorem on its boundary proves equality of the number of
zeros, with multiplicity. The finite count (2) is the negative-strip
degree in the [private-support sheet calculation](MULTIVARIATE_PRIVATE_SUPPORT_SADDLE_GEOMETRY.md),
and agrees with the [shared-support degree](MULTIVARIATE_SADDLE_MAP_TWO_SHEETS.md)
when both private degrees vanish. This proves (2) for the closure.
When a limit is specified only by local uniform convergence, extract
a measure subsequence as above to obtain (5). The retained private
degrees can also be read off from the coefficients of the endpoint
logarithms in (5), so they are independent of the subsequence.

## The escaping sheet as curvature tends to zero

Suppose \(W\to0\), and let \(x\) range in a fixed compact \(K\subset\mathcal T\).
The profiles eventually have shared support, since a private row has
mass at least two. The map \(X_c\) therefore has two preimages for each
\(x\), counted with multiplicity. One converges uniformly to
\[
 q_{\rm flat}(x)=\frac{1}{1+e^x}.
\]
The other escapes to infinity. For the escaping preimage \(q\), the
rescaled coordinate \(v=W(2q-1)\) stays in a common compact subset of
\(\mathbb H\). Along any weak limit
\[
 \lambda_*=\lim\sum_i \frac{w_i}{W}\delta_{r_i/W}
\]
and target limit \(x\in K\), it converges to the unique solution of
\[
 F_{\lambda_*}(v)=x+i\pi,\qquad
 F_\lambda(v)=\int\frac{\log(1+\rho v)}{2\rho}\,d\lambda(\rho),
 \tag{6}
\]
with integrand \(v/2\) at \(\rho=0\).

To prove the flat branch statement, note that \(\max |r_i|\le W/2\).
Thus \(S(q)=O(W)\) on every bounded subset of \(\mathbb H\), by the
analytic kernel (4), and the same bound holds on bounded neighborhoods
of any finite real point for sufficiently small \(W\).
The flat map \(q\mapsto\log(1-q)-\log q\) maps \(\mathbb H\)
biholomorphically to \(\mathcal T\). A closed rectangle in
\(\mathbb H\) containing \(q_{\rm flat}(K)\) in its interior
has boundary image a positive distance from \(K\). Rouché gives one
preimage inside for all sufficiently small \(W\), with convergence
uniform in \(x\) by compactness. The other preimage has no finite
limit in \(\mathbb H\), since every such limit is on the flat
branch. A finite real limit other than 0 and 1 would give boundary
height 0 or \(-\pi\), and either endpoint would give divergent real
part. Thus the second preimage escapes to infinity uniformly.

The scaled probability measures \(\lambda=\sum_i(w_i/W)\delta_{r_i/W}\)
are supported on \([-1/2,1/2]\). In the coordinate \(v=W(2q-1)\),
\(S(q)=F_\lambda(v)\) exactly, and
\[
 F_\lambda(v)=x+i\pi-\log(1-1/q)\quad\text{if } X_c(q)=x.
 \tag{7}
\]
For this identity, \(1-1/q\) lies in \(\mathbb H\) and its argument
is \(\arg(1-q)-\arg q+\pi\). The final term in (7) tends to zero
uniformly on the escaping branch. In particular \(F_\lambda(v)\)
ranges in a compact part of \(0<\Im y<\pi\) for small \(W\).
The proof excluding infinity above applies to \(v\) with
\(r_i\) replaced by \(\rho_i=r_i/W\), \(w_i\) by
\(2|d_i||\rho_i|\), and total mass one. It uses just the bound
\(|\rho_i|\le1\), the integer gaps, bounded real part and the upper
imaginary budget.

For a possible finite real limit \(v_0\), choose a cutoff
\(|\rho_i|\ge\rho_0>0\) as in the finite-boundary argument.
There are at most \(1/(2\rho_0)\) such rows, and their integer gaps
are bounded by the same number. Pass to a subsequence with fixed
numbers and gaps, and convergent \(\rho_i\) for these rows.
All small factors give bounded real part and uniformly vanishing
imaginary part near \(v_0\). If any retained factor tends to zero,
all factors doing so have coefficients of the same sign, as
\(1+\rho_i v_0=0\) forces the sign of \(\rho_i\). Their logarithms
give a real divergence. Otherwise the imaginary limit is an
integer multiple of \(\pi\). Both conclusions contradict (7), so
\(v\) stays in a common compact subset of \(\mathbb H\).

The kernels in (6) converge locally uniformly under weak convergence
of the measures, by the same compact-family argument as in (5).
Consequently a solution of (6) exists for every target
\(x\in\mathcal T\) and every attainable limit measure, by applying
the compactness argument to that target. To see uniqueness, for
any probability measure on \([-1/2,1/2]\) and \(v\in\mathbb H\) we have
\[
 \Im\!\big(vF_\lambda'(v)\big)
   =\frac{\Im v}{2}\int |1+\rho v|^{-2}\,d\lambda(\rho)>0.
\]
The derivative of \(u\mapsto F_\lambda(e^u)\) thus has positive
imaginary part on the convex strip \(0<\Im u<\pi\). For two distinct
points in that strip, the difference quotient is the average of
this derivative along their segment and has positive imaginary
part. This shows injectivity, hence uniqueness in (6) and
convergence of the escaping branch on every specified subsequence.

For clarity, the attainable probability measures in (6) are precisely
a summable list of atoms \(2|d_i||\rho_i|\delta_{\rho_i}\) with
\(0<|\rho_i|\le1/2\) and nonzero integer \(d_i\) of the same sign as
\(\rho_i\), plus the remaining mass at zero. Necessity is the
fixed-count, fixed-gap subsequence argument at nested cutoffs.
For sufficiency, temporarily truncate the retained list. Choose a
positive scale \(\varepsilon\to0\), and for each retained atom
choose an admissible half-integer midpoint
\(m_i=|d_i|/(2\varepsilon|\rho_i|)+O(1)\), so that
\((m_i-d_i/2,m_i+d_i/2)\) are nonnegative integers. Then
\(r_i/\varepsilon\to\rho_i\) and
\(w_i/\varepsilon\to2|d_i||\rho_i|\).
Add copies of a positive unit-gap row with midpoint of size
\(\varepsilon^{-2}\) to supply the remaining curvature to
\(o(\varepsilon)\), using integer rounding of the copy count.
Their individual \(r_i/\varepsilon\) tend to zero, and
\(W/\varepsilon\to1\). A diagonal sequence with longer retained
truncations gives the asserted probability measure.

The theorem concerns stationary fibers on compact parts of the
negative strip. A uniform comparison of consecutive polynomial zeros
requires contour and amplitude information in addition to these
stationary locations.
