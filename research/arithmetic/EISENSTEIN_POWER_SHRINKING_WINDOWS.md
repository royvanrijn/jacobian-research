# Power-shrinking windows for Eisenstein binary sums

Let \(R=\mathbf Z[\omega]\), \(\omega^2+\omega+1=0\), and put
\[
A=\{a\ge1:p\mid a\Longrightarrow p\equiv1\pmod3\}.
\]
The number \(1\) belongs to \(A\), and \(3\cdot A\) denotes dilation.
The [fixed-interval angular theorem](EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md)
controls a fixed fraction of a target. Here the relative window shrinks as a
power. The uniform Hecke estimate permits the entire Legendre expansion of a
localized triangular weight, while an annulus count keeps the large-prime loss
inside that same window.

**Theorem EisPower.** Fix \(0<\delta<1/3\), and define
\[
\epsilon=\frac{1-3\delta}{8},\qquad
\eta=\frac{1-3\delta}{4}=2\epsilon,\qquad
r=\left\lceil\frac1\epsilon\right\rceil,\qquad
D_r=(2r)^{2^r},
\]
\[
N_\delta=\left\lceil\max\left\{
2^{1/\delta},\ (10000D_r)^{1/\epsilon}\right\}\right\rceil.
\]
For every integer \(n\ge N_\delta\):

- if \(n\equiv2\pmod6\), then \(n=a+b\) with \(a,b\in A\) and
  \(|a-b|\le n^{1-\delta}\);
- if \(n\equiv4\pmod6\), then \(n=a+3b\) with \(a,b\in A\) and
  \(|a-3b|\le n^{1-\delta}\).

In the mixed case the actual summands are \(a\) and \(3b\). In either
prescribed orientation there are at least
\(n^{-\delta}R_o(n)/20\) signed, ordered coordinate witnesses, where
\[
R_o(n)=12\{\sigma(n)-\sigma(n/4)\}\ge18n.
\]
A quantity at \(n/4\) is zero if \(4\nmid n\). For every fixed
\(\nu>0\), there are consequently
\(\gg_{\delta,\nu}n^{1-\delta-\nu}\) distinct ordered value pairs.
Equivalently, every fixed exponent \(2/3<\theta<1\) is available for the
actual summand difference. The separate
[pointwise semigroup-scale theorem EisSemigroupScale](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md)
now gives matching bounds of order \(n^\theta J(n)/\log n\) for distinct
ordered pairs in both families, where
\(J(n)=\prod_{p\mid n,\ p\ge5,\ p\equiv2(3)}(p-1)/(p-2)\).
Its cofactor argument supplies an existential onset. The
[moving-window extension EisLocationWindows](EISENSTEIN_MOVING_WINDOWS.md)
allows a center \(x\) with \(M=\min(x,n-x)\), under
\(H\le M\) and \(H^3\ge n^{1+\varepsilon}M\). In particular a
designated small summand may lie between fixed multiples of \(n^\theta\)
for any \(\theta>1/2\). This is a position-dependent window result; its
central threshold remains \(\theta>2/3\). It supplies no explicit onset. The
[value-count theorem EisValueCounts](EISENSTEIN_POINTWISE_VALUE_COUNTS.md)
strengthens the distinct-pair consequence to
\(\gg_{\theta,\nu}n^\theta/(\log n)^{\log4+\nu}\) for every fixed
\(\nu>0\). Its moment argument supplies no numerical onset for that stronger
value-count bound. The separate
[bounded-factor theorem EisBoundedFactors](EISENSTEIN_BOUNDED_PRIME_FACTORS.md)
also gives \(\gg_\theta n^\theta/(\log n)^2\) pairs for which every prime
factor of each \(A\)-part exceeds \(n^{\eta_\theta}\), for some fixed
\(\eta_\theta>0\). That additional restriction has no numerical onset here.

**Corollary EisPowerThreeQuarters.** Both conclusions hold with difference
at most \(n^{3/4}\) for every admissible \(n\ge2^{2^{40}}\), with at
least \(n^{-1/4}R_o(n)/20\) coordinate witnesses per orientation.

The original \(99/100\) statement and its witness bound remain valid at
\(n\ge2^{2^{411}}\). The general theorem in fact supplies a smaller
explicit onset for that exponent. None of these cutoffs connects the proof to
a practical finite certificate or identifies the complete small-exception lists.
The pointwise semigroup-scale count is proved in the separate cofactor theorem
throughout this exponent range. Its exact leading constant remains open.

## Quaternionic input and notation

Use the definite discriminant-3 order
\(\mathcal O=R+Rj\), with \(j^2=-1\), \(jz=\bar z j\), and
\[
\operatorname{nrd}(z+wj)=|z|^2+|w|^2.
\]
The manuscript establishes its maximality, class number one, its twelve
units \(R^\times\sqcup R^\times j\), and
\(\#\mathcal Q_n=12\sigma(n)\) for \(3\nmid n\), where
\(\mathcal Q_n=\{\alpha\in\mathcal O:\operatorname{nrd}\alpha=n\}\).
[Norm count](../../results/eisenstein-restricted-prime-sums/sections/norm-count.tex).

Let \(\Omega_n\subset\mathcal Q_n\) have both component norms odd.
For even \(n\) with \(3\nmid n\), it has the displayed cardinality
\(R_o(n)\). For \(\alpha=z+wj\) put
\[
t_\alpha=\frac{|z|^2-|w|^2}{n}.
\]
For an odd inert prime \(p\equiv2\pmod3\), let \(b_p^o(n)\) count
points of \(\Omega_n\) with a component norm divisible by \(p\), and set
\(\mu_p=2/[p(p-1)]\).

The scalar orbit estimate and odd subtraction give the uniform bound
\[
b_p^o(n)\le\mu_pR_o(n)+24\tau(n)\sqrt n.\tag{1}
\]
For \(p\nmid n\), the complementary norm of
\(D_n=T_n-T_{n/4}\) is at most \(2\tau(n)\sqrt n\); the orbit
quadratic-form estimate multiplies it by at most 12. If \(p\mid n\),
anisotropy forces the entire quaternion to be divisible by \(p\):
the bad set is empty unless \(p^2\mid n\), and otherwise its relative
size is
\(R_o(n/p^2)/R_o(n)<p^{-2}\le\mu_p\).
[Scalar Hecke operator and averaging](../../results/eisenstein-restricted-prime-sums/sections/prime-bounds.tex),
[odd filter](../../results/eisenstein-restricted-prime-sums/sections/mixed-local.tex).

The lattice estimate supplies
\[
\frac{b_p^o(n)}{R_o(n)}\le\frac{4M(n)}{p^2},\qquad
M(n)=\max_{u\le n}\prod_{q^a\Vert u,\ q\equiv1(3)}(a+1)
\le\max_{u\le n}\tau(u).\tag{2}
\]
[Large-prime bound](../../results/eisenstein-restricted-prime-sums/sections/elementary-estimates.tex).
Also
\[
\sum_{p\ge5,\ p\equiv2(3)}\mu_p<\frac17.\tag{3}
\]
For example, majorize by all \(6k-1\), sum \(k\le10\), and bound the
remaining decreasing series by \(\frac13\log(59/58)<1/174\).
The resulting rational upper bound is
\(1172270798426263/8207288024655720<1/7\), as in the retained
[finite-polynomial argument](EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md#a-direct-finite-polynomial-proof-of-eventual-balance).

## Uniform archimedean weight and finite level

For even \(d\ge2\), let \(\rho_d\) be the unitary representation
\(\operatorname{Sym}^{2d}(\mathbf C^2)\) of \(\mathrm{SU}(2)\), extended
to the projective quaternion group by normalizing \(\alpha\) to
\(u_\alpha=\alpha/\sqrt{\operatorname{nrd}\alpha}\).
It descends through the central signs. Let \(v_d\) be the unit vector
proportional to \(X^dY^d\). The elementary coefficient identity is
\[
\langle v_d,\rho_d(u_\alpha)v_d\rangle=L_d(t_\alpha),\tag{4}
\]
where \(L_d\) is the Legendre polynomial with \(L_d(1)=1\).
Indeed the coefficient is
\(\sum_r(-1)^r\binom dr^2|z/\sqrt n|^{2d-2r}|w/\sqrt n|^{2r}\).
The torus fixes \(v_d\), and the Weyl element \(j\) acts on it by
\((-1)^d=1\). Thus every order unit fixes \(v_d\).

Fix \(p\nmid n\), \(p\ge5\) inert. Write
\[
G_p=\mathrm{PGL}_2(\mathbf F_p),\quad
\Gamma_p=\operatorname{image}(\mathcal O^\times),\quad
Y_p=\Gamma_p\backslash G_p.
\]
Let \(K_p\) be the nonsplit torus normalizer and
\(O_p=\Gamma_p\backslash K_p\). Only the units \(\pm1\) become scalars
modulo \(p\), and they act trivially on \(\rho_d\), so the coefficient
action of \(\Gamma_p\) is well-defined. Use the Hilbert space
\[
\mathcal H_{d,p}=\{F:G_p\to V_d:
 F(\gamma g)=\rho_d(\gamma)F(g)\},\qquad
\|F\|^2=\sum_{\Gamma_pg}\|F(g)\|^2.
\]
This is the coefficient-valued algebraic automorphic space at archimedean
weight \(\rho_d\), with principal reduction level at \(p\), order-unit
level at 3, and maximal level elsewhere. Class number one identifies it
with these finite equivariant functions.
[Greenberg--Voight, Section 2, Lemma 2.4, Hecke formula (4) and level
formulas (5)--(6)](https://jvoight.github.io/articles/orthohecke-240123.pdf).

For \((n,3p)=1\), its all-ideal operator is
\[
(T_nF)(g)=\frac1{12}\sum_{\alpha\in\mathcal Q_n}
 \rho_d(u_\alpha)^{-1}F(\alpha_pg).\tag{5}
\]
For the ideal normalization, first write the equivalent formula
\(\frac1{12}\sum_\alpha\rho_d(u_\alpha)F(\alpha_p^{-1}g)\).
Replacing \(\alpha\) by \(\alpha u\), for an order unit \(u\), cancels
\(\rho_d(u)\) against the left equivariance of \(F\). Thus the sum is
once per right-ideal generator orbit. Conjugating all generators gives
(5), exactly the same all-ideal convention as the scalar operator.
Replacing \(g\) by \(\gamma g\) and \(\alpha\) by
\(\alpha\gamma\) proves that (5) preserves equivariance. Conjugation
also proves self-adjointness, since
\(u_{\bar\alpha}=u_\alpha^{-1}\) and
\(\bar\alpha_p=\alpha_p^{-1}\) projectively. The operator commutes
with right \(G_p\)-translation.

Its spectral bound on the **whole** positive-weight space is
\[
\|T_n\|_{\mathcal H_{d,p}}\le\tau(n)\sqrt n.\tag{6}
\]
The bound is uniform in \(d\) and \(p\). At maximal finite level it
holds on \(V_d^{\mathcal O^\times}\) as well. The cuspidality and
normalization checks are given in the next subsection.
The finite Hilbert norm and the scalar automorphic \(L^2\) norm differ
only by a positive overall factor: Schur orthogonality integrates the
archimedean matrix coefficient to
\(\|v\|^2\|F(g)\|^2/(2d+1)\). The finite arithmetic stabilizers give
the same counting measure as in the scalar quotient. This common factor
cancels from operator norms; it creates no degree loss in (6).

Now take \(F_p(g)=\mathbf1_{O_p}(\Gamma_pg)v_d\).
It is equivariant, right \(K_p\)-invariant, and has squared norm
\(|O_p|\). So is \(T_nF_p\) right \(K_p\)-invariant. Hence
\[
\langle v_d,(T_nF_p)(1)\rangle
=\frac{\langle F_p,T_nF_p\rangle}{|O_p|}.
\]
Multiplying the left side by 12 gives the Legendre moment over the
\(p\)-bad set, by (4); taking the inverse of \(u_\alpha\) leaves the
same real coefficient. Thus
\[
\left|\sum_{\mathcal Q_n,\ p\text{-bad}}L_d(t_\alpha)\right|
\le12\tau(n)\sqrt n.\tag{7}
\]
At maximal level, the unit vector \(v_d\) gives identically
\(\left|\sum_{\mathcal Q_n}L_d(t_\alpha)\right|
\le12\tau(n)\sqrt n\).

Odd filtering subtracts the scaled norm-\(n/4\) shell; multiplication
by 2 leaves both normalized archimedean and finite projective images
unchanged. Consequently
\[
\left|\sum_{\Omega_n}L_d(t_\alpha)\right|\le24\tau(n)\sqrt n,
\qquad
\left|\sum_{\Omega_n,\ p\text{-bad}}L_d(t_\alpha)\right|
\le24\tau(n)\sqrt n.\tag{8}
\]
For the second assertion when \(p\mid n\), anisotropy identifies the
bad set with \(\Omega_{n/p^2}\), or the empty set. Its moment is bounded
by \(24\tau(n/p^2)\sqrt{n/p^2}\), which satisfies (8). Thus the
second estimate is valid for every odd inert prime, including divisors
of the target.

### Cuspidality at positive archimedean weight

Original Jacquet--Langlands Theorem 14.4 requires infinite-dimensional
local factors at split places. Here is a direct verification that covers
the varying archimedean weight.
For every finite \(\ell\ne3\), the group
\(\Gamma_\ell=\mathcal O[1/\ell]^1\) has dense archimedean image in
\(\mathrm{SU}(2)\). Indeed the count
\(\#\mathcal Q_{\ell^{2k}}=12\sigma(\ell^{2k})\), after removing
points divisible by \(\ell\), gives infinitely many distinct norm-one
S-units \(\alpha/\ell^k\). Single-component points number at most
\(12(2k+1)\), so some S-unit lies outside the normalizer of the
Eisenstein torus.

The closure is an infinite compact Lie subgroup. Every proper
positive-dimensional Lie subalgebra of \(\mathfrak{su}_2\) is
one-dimensional, since its bracket is the cross product. A proper
infinite closed subgroup therefore normalizes a torus. It contains the
order-three unit \(\omega\); the normalizer modulo its torus has order
two, so \(\omega\) belongs to that torus, forcing it to be the
Eisenstein torus. The S-unit just exhibited contradicts this. The closure
is \(\mathrm{SU}(2)\). Any finite-index subgroup remains dense because
\(\mathrm{SU}(2)\) is connected.

A one-dimensional local constituent at \(\ell\) is trivial on
\(B^1(\mathbf Q_\ell)\). Impose the finite-index congruence condition
on \(\Gamma_\ell\) belonging to the level away from \(\ell\); when
\(\ell=p\), its \(p\)-component is absorbed by that trivial local
action. Right invariance and left equivariance then force each
evaluation vector to be fixed by an archimedean-dense subgroup. Since
\(\rho_d\) is irreducible and nontrivial, that vector is zero. Applying
this to all right \(G_p\)-translates forces the entire constituent to
vanish. Thus every split local factor is infinite-dimensional.

[Jacquet--Langlands, Section 14, Theorem 14.4, printed p. 248](https://publications.ias.edu/sites/default/files/automorphic-forms-on-gl2_rpl_9.pdf)
therefore transfers every constituent here to a cusp representation.
The archimedean representation
\(\operatorname{nrd}^{-d}\operatorname{Sym}^{2d}\) corresponds to
holomorphic weight \(2d+2\), with trivial central character;
[Section 5, Proposition 5.20](https://publications.ias.edu/sites/default/files/automorphic-forms-on-gl2_rpl_9.pdf)
gives its local parameter. At every prime \(\ell\mid n\) the level
is maximal and the representation unramified. Deligne's
[classical-weight Frobenius relation and purity, Theorems 4.9 and 5.1](https://www.numdam.org/item/SB_1968-1969__11__139_0.pdf),
combined with [Weil I, Theorem 1.6](https://www.numdam.org/item/PMIHES_1974__43__273_0.pdf),
give classical roots of modulus \(\ell^{d+1/2}\).
The unitary transport in (5) divides the classical coefficient \(a_n\)
by \(n^d\). Thus its roots have modulus \(\sqrt\ell\), their product
is \(\ell\), and the all-ideal recurrence is
\[
T_{\ell^{r+2}}=T_\ell T_{\ell^{r+1}}-\ell T_{\ell^r}.
\]
Equivalently the unnormalized degree-\(\ell+1\) spherical correspondence
has eigenvalue \(\sqrt\ell(\alpha+\beta)\) for unitary Satake parameters
of modulus one. The bound \((r+1)\ell^{r/2}\), multiplicativity and
self-adjointness prove (6) with no level or weight factor.

## The triangular weight and its Legendre coefficients

**Lemma EisTriangle.** For \(0<h\le1/2\), let
\[
b_h(t)=\max(1-|t|/h,0),\qquad
c_d=\frac{2d+1}{2}\int_{-1}^1 b_h(t)L_d(t)\,dt.
\]
Then \(c_0=h/2\), all odd coefficients vanish, and
\[
 |c_d|\le\min\left(3h\sqrt d,\frac8{h d^{3/2}}\right)\quad(d\ge1),
 \qquad
 \sum_{d\ge1}|c_d|\le32h^{-1/2}=:B(h).\tag{9}
\]
Moreover
\[
b_h(t)=\sum_{d\ge0}c_dL_d(t)\tag{10}
\]
converges absolutely uniformly on \([-1,1]\).

For the needed interior bound, the elementary Laplace integral is
\[
L_d(t)=\frac1\pi\int_0^\pi
 \bigl(t+i\sqrt{1-t^2}\cos u\bigr)^d\,du.
\]
Summing in \(d\) near the origin gives the generating function
\((1-2tz+z^2)^{-1/2}\), which proves the identity. The integrand's
modulus is at most \(\exp[-d(1-t^2)\sin^2u/2]\). By symmetry,
\(\sin u\ge2u/\pi\) on \([0,\pi/2]\), and the Gaussian integral,
\[
 |L_d(t)|\le\sqrt{\frac{\pi}{2d(1-t^2)}}\le\frac2{\sqrt d}
 \qquad(|t|\le1/2).
\]
Since \(\int b_h=h\), this gives the first coefficient bound.

For the second use the Legendre differential equation. In distributions,
\[
 c_d=-\frac{2d+1}{2d(d+1)}
       \int L_d(t)\,d\{(1-t^2)b_h'(t)\}.
\]
The signed measure in this formula is supported on \([-h,h]\). Its atoms
at \(-h,0,h\) have masses \((1-h^2)/h,-2/h,(1-h^2)/h\), and its
continuous density is \(2|t|/h\) on \((-h,h)\). Its total variation is
exactly \(4/h\). Combining the interior bound with
\((2d+1)/(2d(d+1))\le1/d\) proves the second estimate in (9).
The integration by parts has no exterior boundary terms because \(b_h\)
vanishes off its compact interior support.

Split the sum at \(D=\lceil1/h\rceil\le2/h\). The first bound contributes
at most \(3hD^{3/2}\le6\sqrt2\,h^{-1/2}\); the integral test on the
second contributes at most \(16/(h\sqrt D)\le16h^{-1/2}\).
Their sum is less than \(25h^{-1/2}\), so the envelope \(B(h)\) in (9)
is sufficient and is at least one. Since \(|L_d|\le1\) on the entire
interval, the series converges absolutely uniformly. Its Legendre
coefficients agree with those of \(b_h\). Orthogonality and density of
polynomials in continuous functions identify the sum with \(b_h\).
Thus the uniform individual-degree estimates (8) may be summed over (10)
without an uncontrolled degree limit.

## The fixed condition at 3 for mixed targets

At the ramified place use
\[
G_3=(\mathcal O/3\mathcal O)^\times/\mathbf F_3^\times,
\qquad\Gamma_3=\operatorname{image}(\mathcal O^\times),
\qquad Y_3=\Gamma_3\backslash G_3.
\]
Their cardinalities are \(36,6,6\), respectively. The finite level at
3 is the reduction kernel, with maximal level elsewhere, exactly as in
the [six-point construction](../../results/eisenstein-restricted-prime-sums/sections/mixed-local.tex).
The same coefficient-valued construction and bound (6) hold for
\(3\nmid n\). The density argument above imposes only a finite-index
condition modulo 3 on \(\Gamma_\ell\).

A component norm is divisible by 9 precisely when its component is zero
modulo \(3R\). On an invertible quaternion this is exactly the identity
coset of \(Y_3\). Let \(\delta_3\) indicate that coset, set
\(f_3=1-\delta_3\), and take \(F_3=f_3v_d\). Its squared norm is 5,
so it is at most 6. Evaluation at the identity coset has norm at most
one. Therefore
\[
|\langle v_d,(T_nF_3)(1)\rangle|
\le\sqrt6\,\tau(n)\sqrt n.
\]
Multiplication by 12 gives the moment on the shell with neither component
norm divisible by 9. Odd subtraction preserves this condition, since
2 is a unit scalar modulo 3. Thus, for the 3-good odd shell
\(\mathcal H_n\),
\[
\left|\sum_{\mathcal H_n}L_d(t_\alpha)\right|
\le24\sqrt6\,\tau(n)\sqrt n.\tag{11}
\]
The estimate is uniform in even \(d\). It involves this one fixed
3-level space, not a simultaneous level \(3p\).

For \(n\equiv4\pmod6\), the independently established
[eta identity EisEta](EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md#proposition-eiseta-the-obstruction-at-3)
gives
\[
\left||\mathcal H_n|-\tfrac23R_o(n)\right|
\le12\tau(n)\sqrt n.\tag{12}
\]
Indeed the removed count is
\(R_o/3+8(c_n-c_{n/4})\), where
\(\eta(3\tau)^2\eta(9\tau)^2=\sum c_nq^n\) has
\(|c_n|\le\tau(n)\sqrt n\) for \(3\nmid n\).

## A large-prime tail inside the central window

**Lemma EisCentralTail.** Let \(b_{p,h}^o(n)\) count points of \(\Omega_n\)
with \(|t_\alpha|\le h\) and a component norm divisible by the odd inert
prime \(p\). For \(0<h\le1/2\) and \(P\ge2\),
\[
 \frac1{R_o(n)}\sum_{p>P}b_{p,h}^o(n)
 \le64M(n)\left\{\frac hP+\frac{1+\log n}{\sqrt n}\right\}.
 \tag{CT}
\]
This estimate includes primes dividing the target and all singular coordinate
representations.

If \(p\) divides the norm of \(z\in R\), anisotropy gives \(z=p\zeta\).
In the central window, \(\zeta\) lies in the annulus
\[
 r_-^2=\frac{n(1-h)}{2p^2}\le|\zeta|^2
 \le\frac{n(1+h)}{2p^2}=r_+^2.
\]
The centered fundamental parallelogram of \(R\) has radius at most one
and area \(\sqrt3/2>1/2\). Attach this cell to every point of the annulus.
The cells have disjoint interiors, and their union lies between radii
\(\max(0,r_--1)\) and \(r_++1\). If \(r_-\ge1\), the enclosing
area is \(\pi(r_+^2-r_-^2+2r_++2r_-)\). If \(r_-<1\), discarding
the inner hole costs at most one extra unit in squared radius. In both cases
the lattice count is at most
\[
2\pi\{r_+^2-r_-^2+2(r_++r_-)+2\}
 \le32\left\{\frac{hn}{p^2}+\frac{\sqrt n}p+1\right\}.
\]
Here \(r_+\le\sqrt n/p\) and \(r_-\le r_+\), and \(\pi<4\).
For each \(\zeta\), the positive complementary norm
\(n-p^2|\zeta|^2\) has at most \(6M(n)\) representations. Allowing
either component to be bad and discarding the oddness restriction gives
\[
 b_{p,h}^o(n)\le384M(n)
 \left\{\frac{hn}{p^2}+\frac{\sqrt n}p+1\right\}.
\]
Only \(p\le\sqrt n\) can occur. Sum over integers in place of primes,
using
\(\sum_{p>P}p^{-2}\le2/P\),
\(\sum_{p\le\sqrt n}p^{-1}\le1+\tfrac12\log n\), and at most
\(\sqrt n\) primes. Dividing by \(R_o\ge18n\) gives (CT).
No residue equidistribution or restriction on \(v_p(n)\) is used here.

## A nonnegative weighted sieve

Put \(h=n^{-\delta}\) and \(P=n^\eta\), with the theorem's \(\eta\).
For pure targets let \(\mathcal S_n=\Omega_n\); for mixed targets let
\(\mathcal S_n=\mathcal H_n\), initially retaining both orientations.
Let \(G_n\) be the subset with no odd inert prime in either component norm.
The nonnegative weight \(b_h\) gives
\[
 \#\{\alpha\in G_n:|t_\alpha|\le h\}
 \ge\sum_{\mathcal S_n}b_h
 -\sum_{p\le P}\sum_{\Omega_n,\ p\text{-bad}}b_h
 -\sum_{p>P}b_{p,h}^o(n).\tag{13}
\]
The small-prime subtraction over the full odd shell only enlarges the loss
in the mixed case. For large primes the weight vanishes off the central
window and is at most one there, so (CT) applies. No estimate for a joint
\(3p\)-level space is required.

By (1), (8) and absolute convergence in (10), each small-prime loss is at
most
\[
 \frac h2\mu_pR_o
 +24\{B(h)+h/2\}\tau(n)\sqrt n
 \le\frac h2\mu_pR_o+48B(h)\tau(n)\sqrt n.
\]
The pure candidate error is at most \(24B(h)\tau(n)\sqrt n\).
In the mixed case (11)--(12) give error at most
\((24\sqrt6 B(h)+6h)\tau(n)\sqrt n\le72B(h)\tau(n)\sqrt n\).
There are at most \(P\) small primes; \(B(h),P\ge1\) and
\(R_o\ge18n\). Thus a common upper bound for the total normalized
spectral error is
\[
 E_{\mathrm{sp}}\le8B(h)P\frac{\tau(n)}{\sqrt n}
 =256h^{-1/2}P\frac{\tau(n)}{\sqrt n}.\tag{14}
\]
Since \(\sum\mu_p<1/7\), the main terms in (13) are at least
\[
 \frac37hR_o\quad\text{in the pure case},\qquad
 \frac{11}{42}hR_o\quad\text{in the mixed case, both orientations}.
\]
The normalized large-prime error is bounded by the right side of (CT).

## Explicit exponents and onset

For every integer \(r\ge1\), the elementary divisor bound is
\[
 \tau(u)\le D_r u^{1/r},\qquad M(u)\le D_r u^{1/r},
 \qquad D_r=(2r)^{2^r}.\tag{15}
\]
For \(p\ge2^r\), \(a+1\le2^a\le p^{a/r}\).
For \(p<2^r\), put \(x=(\log p)/r<\log2<1\). The continuous
maximum of \((a+1)e^{-xa}\), \(a\ge0\), is at most
\(1/x\le r/\log2<2r\). Fewer than \(2^r\) primes occur below this
threshold, proving the first bound. Taking the maximum over integers at
most \(u\) proves the second.

For \(\epsilon=(1-3\delta)/8\), \(r=\lceil1/\epsilon\rceil\),
and \(\eta=2\epsilon\), we have \(1/r\le\epsilon\). Dividing
(14) and (CT) by \(h\), (15) gives
\[
 \frac{E_{\mathrm{sp}}}{h}
 \le256D_r n^{-1/2+3\delta/2+\eta+\epsilon}
 =256D_r n^{-\epsilon},\tag{16}
\]
\[
 \frac{E_{\mathrm{tail}}}{h}
 \le64D_r n^{-\epsilon}
 \left\{1+(1+\log n)n^{-(1-\delta)/4}\right\}
 \le256D_r n^{-\epsilon}.\tag{17}
\]
For the last inequality use \((1-\delta)/4>1/6\) and
\[
 \sup_{n\ge1}(1+\log n)n^{-1/6}=6e^{-5/6}<3.
\]
These exponents explain the range \(\delta<1/3\). There is no assertion
at the endpoint.

At \(n\ge N_\delta\), we have \(h\le1/2\),
\(n^\epsilon\ge10000D_r\), and
\(P=n^{2\epsilon}\ge(10000D_r)^2\ge2\).
Each normalized error in (16)--(17) is therefore less than \(h/20\).
The surviving central weighted mass is at least
\[
 \frac{23}{70}hR_o\quad\text{in the pure case},\qquad
 \frac{17}{105}hR_o\quad\text{in the mixed case, both orientations}.
\]
The even weight is invariant under swapping components, which interchanges
the two mixed orientations. Exactly one norm is divisible by 3, so those
orientations are disjoint and have equal weighted totals. Each receives at
least \(17hR_o/210>hR_o/20\). The pure bound also exceeds \(hR_o/20\).

Since \(b_h\le1\), these are lower bounds for actual coordinate witness
counts. The support of \(b_h\) gives the asserted difference bound.
For pure targets \(n\equiv2\pmod3\) forces both norms to be nonzero
modulo 3. For mixed targets the fixed 3-condition makes the unique divisible
norm have exponent exactly one. Oddness and the exclusion of odd inert
primes now give the required members of \(A\) and \(3\cdot A\).
Each ordered value pair has at most
\(36M(n)^2\ll_\nu n^\nu\) coordinate witnesses, proving the distinct-pair
statement. Apply the elementary divisor estimate with exponent \(\nu/2\)
for this last conclusion.

For the explicit \(3/4\) corollary, take \(\delta=1/4\). Then
\(\epsilon=1/32\), \(\eta=1/16\), \(r=32\), and
\(D_r=2^{6\cdot2^{32}}\). At \(n=2^{2^{40}}\),
\[
 \log_2(n^\epsilon)=8\cdot2^{32}
 >6\cdot2^{32}+\log_2(10000)=\log_2(10000D_r).
\]
Also \(h\le1/2\), so this target exceeds \(N_{1/4}\).
For the original \(\delta=1/100\), \(\epsilon=97/800\), \(r=9\),
and \(D_r=18^{512}\). Using \(18<32\), \(10000<2^{14}\), and
\(800/97<9\), the new onset is at most \(2^{2^{15}}\), which preserves
the original \(2^{2^{411}}\) corollary and its full witness count.

## The Jackson--Parseval route also works

**Lemma EisJacksonParseval.** For the same triangular bump one can construct
even polynomials \(W_h,Q_h\) with degree \(O(h^{-2})\),
\[
 W_h\le b_h,\quad Q_h=W_h+h/100,\quad0\le Q_h\le1,\quad
 \|W_h\|_{L^2(dt/2)}\ll\sqrt h.
\]
Their Legendre coefficients satisfy
\(\sum|c_d|\ll h^{-3/2}\). Using the full-shell tail (2), this alone
proves the weaker range \(\delta<1/7\), or \(\theta>6/7\).

Here is a direct positive-kernel construction, so no approximation constant
is left as an imported hypothesis. For an integer \(m\ge1\), put
\[
 D_m(u)=\frac{\sin(mu/2)}{\sin(u/2)},\qquad
 A_m=\frac1{2\pi}\int_{-\pi}^{\pi}D_m(u)^4\,du
     =\frac{2m^3+m}{3},\qquad K_m=D_m^4/A_m.
\]
The value at zero is interpreted by continuity. The fourth power is a
nonnegative trigonometric polynomial of degree \(2m-2\), and \(K_m\)
has integral one against \(du/(2\pi)\). The exact normalization follows
from \(m^2+2\sum_{j=1}^{m-1}j^2\). Since
\(|D_m(u)|\le\min(m,\pi/|u|)\), splitting at \(\pi/m\) gives
\[
 \frac1{2\pi}\int_{-\pi}^{\pi}|u|K_m(u)\,du
 \le\frac{3\pi}{2m}<\frac5m.
\]
Convolve \(K_m\) with \(b_h(\cos\phi)\), whose circle Lipschitz
constant is \(1/h\). This produces an even, \(\pi\)-periodic
trigonometric polynomial, hence an even algebraic polynomial \(Q_h(t)\)
of degree at most \(J=2m-2\). Positivity gives \(0\le Q_h\le1\).
For \(m=\lceil500h^{-2}\rceil\), its uniform error is at most
\(5/(mh)\le h/100\). Set \(W_h=Q_h-h/100\). Then its constant
Legendre coefficient lies between \(12h/25\) and \(h/2\), and
\[
 \|W_h\|_2\le\sqrt{h/3}+h/50<\tfrac35\sqrt h.
\]
Orthogonality, Parseval and Cauchy--Schwarz give
\[
 \sum_{d\ge2}|c_d|
 \le\sqrt{\sum_{d=0}^J(2d+1)}
       \sqrt{\sum_{d=0}^J\frac{c_d^2}{2d+1}}
 =(J+1)\|W_h\|_2\ll h^{-3/2}.
\]
With \(h=n^{-\delta}\), a sieve cutoff \(n^s\), and a divisor exponent
\(\varepsilon>0\), the normalized errors from this construction and the
old full-shell tail are
\[
 O_\varepsilon\left(
 n^{-1/2+\varepsilon+s+3\delta/2}
 +n^{\varepsilon-s}\right).
\]
Both are \(o(h)\) if
\(\delta+\varepsilon<s<1/2-5\delta/2-\varepsilon\), possible for every
\(\delta<1/7\). The nonnegative majorant controls the small-prime union;
the bound \(W_h\le1\) controls the large-prime union. The mixed fixed-3
estimate is unchanged. This verifies the polynomial upgrade directly. The
stronger theorem above uses the actual triangle coefficients and the central
annulus tail instead of this coarser pair of estimates.

## What individual moment bounds can force

The exponent \(2/3\) is a limit of this particular input, rather than of
Eisenstein representations themselves. The triangle in (9) has the optimal
order of absolute Legendre coefficient mass among central minorants with a
comparable positive average.

**Proposition EisMomentBarrier.** Write \(d\sigma=dt/2\) on \([-1,1]\).
For every \(0<h\le1/4\) there is a symmetric probability measure \(\nu_h\)
with \(\nu_h([-h,h])=0\) such that, with one absolute constant \(C\),
\[
 \left|\int L_d\,d\nu_h\right|\le C h^{3/2}\qquad(d\ge1).
 \tag{18}
\]
Consequently, if \(F\le\mathbf1_{[-h,h]}\),
\(\int F\,d\sigma\ge\kappa h\) for a fixed \(\kappa>0\), and
\(F=\sum_{d\ge0}c_dL_d\) is absolutely summable, then
\[
 \sum_{d\ge1}|c_d|\ge \frac{\kappa}{C}h^{-1/2}.\tag{19}
\]
The assertion includes even polynomials and even absolutely summable series.

To construct the measure, push \(\sigma\) through the map that sends
\(0<|t|<h\) to \(\operatorname{sgn}(t)(|t|+h)\) and leaves the other
points fixed. It moves mass \(h\) by distance \(h\), entirely within
\([-2h,2h]\), and removes all central mass. Endpoints and the point zero
have zero measure. Its symmetry makes all odd moments vanish.

On \(|t|\le1/2\), the interior bound proved above and the identity
\[
 (1-t^2)L_d'(t)=d\{L_{d-1}(t)-tL_d(t)\}
\]
give \(|L_d(t)|\ll d^{-1/2}\) and \(|L_d'(t)|\ll\sqrt d\);
\(d=1\) is checked directly. Since \(\int L_d\,d\sigma=0\),
\[
 \left|\int L_d\,d\nu_h\right|
 \ll h\min\{h\sqrt d,d^{-1/2}\}\le C h^{3/2}.
\]
This proves (18). The minorant is nonpositive on the support of \(\nu_h\),
so absolute convergence permits integration term by term and gives
\[
 \kappa h\le\int F\,d(\sigma-\nu_h)
 \le C h^{3/2}\sum_{d\ge1}|c_d|,
\]
which proves (19).

In particular, an individual-moment envelope of size \(\varepsilon\) is
compatible with an empty central interval whenever
\(h\le c\varepsilon^{2/3}\). At the normalized square-root envelope
\(\varepsilon=n^{-1/2}\), positivity, total mass and these separate
moment bounds alone therefore cannot force a difference window smaller in
power than \(n^{2/3}\). The divisor factors in (8) do not improve that
conclusion. This is an obstruction for the stated abstract information; it
is not an arithmetic counterexample. Joint cancellation across degrees or
additional arithmetic information could give smaller actual windows.

## Scope and review

The exponent family and cutoffs are explicit. No computation through an onset
is proposed. The complete small-exception lists remain open. The separate
cofactor theorem proves the expected pointwise pair-count scale in this range,
including fixed central widths and the unrestricted count. The theorem concerns
ordinary Eisenstein arithmetic; the
simultaneous modulus-12 pointwise lower bound and fixed-prime marginals remain
separate open inputs.

The coefficient-valued Hecke normalization, simultaneous weight/level
uniformity, positive-weight cuspidality and mixed fixed-level extension
received renewed internal source/proof audits. The direct triangle expansion,
annulus geometry, Jackson--Parseval construction and effective tuning also
received separate internal derivation checks. The moment-barrier construction
and its coefficient duality received an independent internal check. These are internal mathematical
reviews. No external specialist review, formal verification or established
priority is claimed. No new finite coverage scan is used.
