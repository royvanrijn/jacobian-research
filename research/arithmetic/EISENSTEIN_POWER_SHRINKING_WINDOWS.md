# Power-shrinking windows for Eisenstein binary sums

Let \(R=\mathbf Z[\omega]\), \(\omega^2+\omega+1=0\), and put
\[
A=\{a\ge1:p\mid a\Longrightarrow p\equiv1\pmod3\}.
\]
The number \(1\) belongs to \(A\), and \(3\cdot A\) denotes dilation.
The [fixed-interval angular theorem](EISENSTEIN_SIFTED_ANGULAR_DISTRIBUTION.md)
controls a fixed fraction of a target. Here the allowed relative difference
tends to zero as an explicit power.

**Theorem EisPower.** Put
\[
N_* = 2^{\,2^{411}}.
\]
For every integer \(n\ge N_*\):

- if \(n\equiv2\pmod6\), then \(n=a+b\) with \(a,b\in A\) and
  \(|a-b|\le n^{99/100}\);
- if \(n\equiv4\pmod6\), then \(n=a+3b\) with \(a,b\in A\) and
  \(|a-3b|\le n^{99/100}\).

In the mixed case the actual summands are \(a\) and \(3b\). In either
prescribed orientation there are at least
\(n^{-1/100}R_o(n)/20\) signed, ordered coordinate witnesses, where
\[
R_o(n)=12\{\sigma(n)-\sigma(n/4)\}\ge18n.
\]
As usual, a quantity at \(n/4\) is zero if \(4\nmid n\). For every fixed
\(\varepsilon>0\), there are consequently
\(\gg_\varepsilon n^{99/100-\varepsilon}\) distinct ordered value pairs.

The onset is intentionally coarse and does not connect this theorem to a
practical finite certificate. The proof needs bounds uniform in the harmonic
degree and sifting prime; it never substitutes a growing degree into a
fixed-weight theta estimate with an unspecified implied constant.

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
[Norm count](../../unpublished/eisenstein-restricted-prime-sums/sections/norm-count.tex).

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
[Scalar Hecke operator and averaging](../../unpublished/eisenstein-restricted-prime-sums/sections/prime-bounds.tex),
[odd filter](../../unpublished/eisenstein-restricted-prime-sums/sections/mixed-local.tex).

The lattice estimate supplies
\[
\frac{b_p^o(n)}{R_o(n)}\le\frac{4M(n)}{p^2},\qquad
M(n)=\max_{u\le n}\prod_{q^a\Vert u,\ q\equiv1(3)}(a+1)
\le\max_{u\le n}\tau(u).\tag{2}
\]
[Large-prime bound](../../unpublished/eisenstein-restricted-prime-sums/sections/elementary-estimates.tex).
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

## An elementary shrinking polynomial

Put \(h=n^{-1/100}\), \(\eta=h/100\), and let \(J\) be the least even
integer at least \(10000h^{-4}\). Define
\[
b_h(t)=\max(1-|t|/h,0),\qquad f_h(x)=b_h(2x-1),
\]
and its Bernstein polynomial
\[
\mathcal B_Jf_h(x)=\sum_{r=0}^Jf_h(r/J)\binom Jr x^r(1-x)^{J-r}.
\]
The Lipschitz constant is \(2/h\). If \(S\sim\operatorname{Bin}(J,x)\),
\[
|\mathcal B_Jf_h(x)-f_h(x)|
\le\frac2h\mathbf E|S/J-x|
\le\frac2h\sqrt{\frac{x(1-x)}J}
\le\frac1{h\sqrt J}\le\eta.
\]
Set
\[
W_h(t)=\mathcal B_Jf_h((t+1)/2)-\eta,\qquad Q_h=W_h+\eta.
\]
Then \(W_h\le b_h\le\mathbf1_{|t|\le h}\),
\(W_h\ge-\eta\), and \(0\le Q_h\le1\). The symmetry of \(f_h\)
makes \(W_h\) even. Write
\(W_h=c_0+\sum_{d=2,4,\ldots,J}c_dL_d\). Its normalized mean and
coefficient sum satisfy
\[
\frac{12h}{25}\le c_0\le\frac h2,\quad
c_0+\eta\le\frac{51h}{100},\quad
\sum_{d\ge2}|c_d|\le B_J:=(J+1)^2.\tag{9}
\]
The mean of \(b_h\) is \(h/2\), and uniform approximation followed by
subtraction costs at most \(2\eta\). For the coefficient bound use
\(|L_d|\le1\), \(|W_h|\le1\), and
\(c_d=(2d+1)\int_{-1}^1W_hL_d/2\).
In particular, with \(A_0=10003^2\),
\[
B_J\le A_0h^{-8}.\tag{10}
\]

## The fixed condition at 3 for mixed targets

At the ramified place use
\[
G_3=(\mathcal O/3\mathcal O)^\times/\mathbf F_3^\times,
\qquad\Gamma_3=\operatorname{image}(\mathcal O^\times),
\qquad Y_3=\Gamma_3\backslash G_3.
\]
Their cardinalities are \(36,6,6\), respectively. The finite level at
3 is the reduction kernel, with maximal level elsewhere, exactly as in
the [six-point construction](../../unpublished/eisenstein-restricted-prime-sums/sections/mixed-local.tex).
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

## A weighted union bound with uniform errors

Choose \(P=n^{1/50}\). For pure targets let \(\mathcal S_n=\Omega_n\);
for mixed targets let \(\mathcal S_n=\mathcal H_n\), initially retaining
both orientations. All surviving positive values have the desired prime
factors. The central good count is at least
\[
\sum_{\mathcal S_n}W_h
-\sum_{p\le P}\sum_{\Omega_n,\ p\text{-bad}}Q_h
-\sum_{p>P}b_p^o(n).\tag{13}
\]
To justify this with a signed minorant, split the bad union into the
small-prime union and the remaining large-prime union. On the first use
\(W_h\le Q_h\) and \(Q_h\ge0\), allowing a union bound; on the second
use \(W_h\le1\). In the mixed case subtraction over the full odd shell
only enlarges the loss. No angular estimate for the intersection of a
\(p\)-condition with the condition at 3 is required.

By (1) and (8), each small-prime loss is at most
\[
(c_0+\eta)\mu_pR_o
+24\{B_J+c_0+\eta\}\tau(n)\sqrt n
\le(c_0+\eta)\mu_pR_o+48B_J\tau(n)\sqrt n.
\]
The pure candidate error is at most \(24B_J\tau(n)\sqrt n\).
In the mixed case (11)--(12) give at most
\((24\sqrt6 B_J+12c_0)\tau(n)\sqrt n
\le72B_J\tau(n)\sqrt n\).
There are at most \(P\) small primes. Since \(B_J\ge1\), \(P\ge1\)
and \(R_o\ge18n\), a common upper bound for the total normalized
spectral error is
\[
8B_JP\frac{\tau(n)}{\sqrt n}.\tag{14}
\]
Here (3) and (9) give the pure main term at least
\(57hR_o/140\), and the mixed main term at least
\(173hR_o/700\).

For \(P\ge2\), (2) bounds the normalized large-prime error by
\[
4M(n)\sum_{p>P}p^{-2}\le\frac{8M(n)}P.\tag{15}
\]
All prime sums here are over the odd inert primes.

## An explicit cutoff

The following deliberately loose divisor estimate avoids unspecified
constants. Put
\[
D=1024^{\,2^{400}}=2^{10\cdot2^{400}}.
\]
For every positive integer \(u\),
\[
\tau(u)\le D u^{1/400},\qquad M(u)\le D u^{1/400}.\tag{16}
\]
For \(p\ge2^{400}\), \(a+1\le2^a\le p^{a/400}\).
For \(p<2^{400}\), put \(x=(\log p)/400<\log2<1\).
The continuous maximum of \((a+1)e^{-xa}\), \(a\ge0\), is at most
\(1/x\le400/\log2<1024\). There are fewer than \(2^{400}\) such
primes, proving the first bound. Taking a maximum over integers at most
\(u\) proves the second.

Combining (10), (14)--(16), the two normalized errors are bounded by
\[
8A_0D n^{-159/400},\qquad 8D n^{-7/400}.\tag{17}
\]
After dividing by \(h=n^{-4/400}\), their base-two logarithms, at
\(n\ge N_*=2^{2^{411}}\), are at most
\[
\begin{split}
3+\log_2 A_0-\frac{3918}{5}\,2^{400},\\
3-\frac{134}{25}\,2^{400}.
\end{split}
\]
Both are less than \(-5\); also \(P\ge2\). Each error is therefore
less than \(h/20\). Formula (13) is consequently at least
\[
\frac{43}{140}hR_o\quad\text{in the pure case},\qquad
\frac{103}{700}hR_o\quad\text{in the mixed case, both orientations}.
\]
The weight \(W_h\) is even, so swapping the two components preserves
it and interchanges the two mixed orientations. Exactly one norm is
divisible by 3, so those orientations are disjoint and receive equal
weighted totals. Each orientation thus has at least
\(103hR_o/1400>hR_o/20\) central witnesses. The pure bound is also
larger than \(hR_o/20\).

Every such witness has \(|t_\alpha|\le h\), proving the stated
\(n^{99/100}\) difference bound. In the pure case, \(n\equiv2\pmod3\)
forces both norms to be nonzero modulo 3. In the mixed case the fixed
3-condition makes the unique divisible norm have exponent exactly one.
Together with oddness and the exclusion of odd inert primes, this gives
the required members of \(A\) and \(3\cdot A\).
Each value pair has at most \(36M(n)^2\ll_\varepsilon n^\varepsilon\)
coordinate witnesses, giving the distinct-pair conclusion. For this last
bound use the elementary divisor argument with exponent \(\varepsilon/2\);
its fixed exponent \(1/400\) above serves only to give the explicit onset.

## Scope and review

The exponent and cutoff are explicit, but no computation through the cutoff
is proposed. This does not identify the complete small-exception lists or
give the conjectured pointwise semigroup-scale count. It treats the ordinary
Eisenstein arithmetic, not the simultaneous modulus-12 restriction.
The vector-valued finite quotient, positive-weight cuspidality, operator
normalization, nonnegative union bound and mixed fixed-level extension
received internal proof/source audits. No external specialist review,
formal verification or established priority is claimed. No new large finite
search is used in this argument.
