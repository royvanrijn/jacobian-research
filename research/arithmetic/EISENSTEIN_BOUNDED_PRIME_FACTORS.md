# Eisenstein binary sums with bounded numbers of prime factors

The angular method also gives decompositions in which each permitted part has
only a bounded number of prime factors. The extra step is to sieve out small
split primes as well as inert primes. A composite-level Hecke estimate and an
elementary lower sieve provide enough weighted mass to do this for every
sufficiently large compatible target.

Put
\[
 A=\{a\ge1:p\mid a\Longrightarrow p\equiv1\pmod3\},
 \qquad\Omega(1)=0,
\]
where \(\Omega\) counts prime factors with multiplicity. The number \(1\)
belongs to \(A\).

**Theorem EisBoundedFactors.** For every fixed \(2/3<\theta<1\), there are
constants \(0<\eta_\theta<1\), \(c_\theta>0\), and \(N_\theta\) such
that the following hold for every integer \(n\ge N_\theta\):

- If \(n\equiv2\pmod6\), at least \(c_\theta n^\theta/(\log n)^2\)
  ordered pairs \((a,b)\in A^2\) satisfy
  \[
   a+b=n,\qquad |a-b|\le n^\theta.
  \]
- If \(n\equiv4\pmod6\), at least \(c_\theta n^\theta/(\log n)^2\)
  ordered pairs \((a,b)\in A^2\) satisfy
  \[
   a+3b=n,\qquad |a-3b|\le n^\theta.
  \]

In each case every prime factor of \(ab\) exceeds \(n^{\eta_\theta}\).
Consequently
\[
 \Omega(a),\Omega(b)\le\left\lfloor\frac1{\eta_\theta}\right\rfloor.
\]
The mixed assertion concerns the actual summands \((a,3b)\); the reversed
orientation has the same count. Its fixed factor \(3\) is excluded from the
bound on the prime factors of the \(A\)-part \(b\).

No numerical onset or small value of the bound on \(\Omega\) is asserted.
The theorem does not require a new finite coverage computation. The separate
[cofactor amplification theorem](EISENSTEIN_POINTWISE_SEMIGROUP_SCALE.md)
uses these rough pair counts to obtain the expected semigroup-scale order
\(H J(n)/\log n\), allowing small squarefree factors in the summands.

## The weighted sequence

Use the quaternion order \(\mathcal O=R+Rj\), \(R=\mathbb Z[\omega]\),
and all notation from [EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md).
Write
\[
 \delta=1-\theta,\quad H=n^\theta,\quad h=H/n,\qquad
 t_{z+wj}=\frac{|z|^2-|w|^2}{n},\quad
 b_h(t)=\max(1-|t|/h,0).
\]
All sufficiently large targets have \(h\le1/2\). Let \(\mathcal S_n\)
be the odd-component shell \(\Omega_n\) for \(n\equiv2\pmod6\).
For \(n\equiv4\pmod6\), let it be the shell \(\mathcal H_n\) with
odd component norms and neither norm divisible by \(9\), initially keeping
both mixed orientations. Define
\[
 X_n=\begin{cases}
       (h/2)R_o(n),&n\equiv2\pmod6,\\
       (h/3)R_o(n),&n\equiv4\pmod6.
      \end{cases}
\]
Thus \(X_n\gg H\), since \(R_o(n)\ge18n\).

For squarefree \(d\) supported on primes at least \(5\), let
\[
 M_n(d)=\sum_{\substack{z+wj\in\mathcal S_n\\
                    p\mid |z|^2|w|^2\text{ for every }p\mid d}}
                    b_h(t_{z+wj}).
\]
We will prove that a multiplicative function \(g_n\) on squarefree integers
satisfies, for every fixed \(\varepsilon>0\),
\[
 M_n(d)=X_ng_n(d)
 +O_\varepsilon\!\left(d^{3/2}n^{1/2+\varepsilon}h^{-1/2}\right).
 \tag{1}
\]
The constant is independent of \(d,n\) and all valuations of \(n\).

## The composite level estimate

This subsection proves the extension of the Hecke input used in (1). Take a
raw shell \(\mathcal Q_m\), with \(3\nmid m\), and a squarefree level
\(D\) coprime to \(m\). The local tests will be bad-prime indicators, the
split-prime tests defined below, and optionally the fixed test at \(3\).
They take values in \([0,1]\), are invariant under projective scalar units,
and are invariant under left multiplication by global order units.

Use
\[
 G=\prod_{p\mid D}\mathrm{PGL}_2(\mathbb F_p),
\]
with the additional factor
\(G_3=(\mathcal O/3\mathcal O)^\times/\mathbb F_3^\times\) when
needed. Let \(\Gamma\) be the diagonal image of the order units. Class
number one identifies the scalar automorphic quotient with
\(\Gamma\backslash G\), with counting measure. This is the same
principalization and stabilizer argument as in
[the scalar spectral proof](../../unpublished/eisenstein-restricted-prime-sums/sections/prime-bounds.tex),
with the intersection of the reduction kernels at the indicated places.

The cuspidal complement has all-ideal Hecke norm at most
\(\tau(m)\sqrt m\). Indeed, the scalar proof that a constituent
one-dimensional at a split place is a global reduced-norm character uses
strong approximation for \(B^1\); replacing its finite level changes only
open stabilizers. All other constituents transfer by Jacquet--Langlands to
weight-two cusp forms, and Deligne gives the same unramified bound.
Every prime dividing \(m\) is outside the level.

The exceptional scalar subspace consists exactly of the characters of the
joint reduced-norm squareclasses. At \(p\mid D\), these are the determinant
classes \(\mathbb F_p^\times/\mathbb F_p^{\times2}\); the optional
factor at \(3\) contributes \(\mathbb F_3^\times/\mathbb F_3^{\times2}\).
Projectivization kills squares. At \(3\), scalar units in the reduction
kernel already have norms covering \(1+3\mathbb Z_3\), so no higher unit
norm character survives. Conversely, every product of the stated quadratic
characters is an actual function on the quotient, since order units have
norm one. Directly its Hecke eigenvalue is \(\sigma(m)\) times its value
at the squareclass of \(m\). Additional finite-group characters at \(3\)
that are not norm characters belong to the cuspidal complement.

Orthogonal projection onto this character space, followed by Hecke action
and evaluation at the identity, is therefore \(\sigma(m)\) times the
uniform conditional average at the joint norm squareclass of \(m\).
Uniform quotient averages lift to uniform averages on \(G\), since all
left \(\Gamma\)-orbits have the same size. The conditional average of a
product test is consequently the product of its local conditional means.

For every positive even Legendre degree \(k\), use the coefficient-valued
space in EisPower. Its entire Hecke operator has norm at most \(\tau(m)\sqrt m\):
the positive-weight S-unit density proof remains valid after any finite-index
congruence condition. For a product test \(F\), the section \(Fv_k\)
has norm at most
\[
 \sqrt{|G|}\le6D^{3/2},
\]
and evaluation has norm at most one. If \(D=1\) without the extra factor
at \(3\), use the maximal-level space \(V_k^{\mathcal O^\times}\)
already treated in EisPower. The same estimate holds there.

Now sum the absolutely convergent triangular expansion from EisTriangle:
its constant coefficient is \(h/2\), and its nonconstant coefficient sum
is \(O(h^{-1/2})\). We obtain
\[
 \sum_{\alpha\in\mathcal Q_m}F(\alpha)b_h(t_\alpha)
 =\frac h2\,12\sigma(m)\prod\text{local conditional means}
 +O\!\left(D^{3/2}\tau(m)\sqrt m\,h^{-1/2}\right).
 \tag{2}
\]
The optional \(3\)-good mean is \(2/3\) when \(m\equiv1\pmod3\):
three of the six quotient points have that norm squareclass, and the test
removes one of them. For \(m\equiv2\pmod3\), both component norms are
already nonzero modulo \(3\). Formula (2) uses the same relative parameter
\(h\) even if \(m\) arose by reducing the original target.

## Unit targets at a sieving prime

Suppose first that \(p\nmid n\). At a split prime \(p\equiv1\pmod3\),
write the Eisenstein subalgebra as \(\mathbb F_p\times\mathbb F_p\).
The four coordinate residues occupy the diagonal and off-diagonal positions
of a matrix, up to signs; a component norm is divisible by \(p\) exactly
when an entry vanishes. Among the \(p(p^2-1)\) matrices of any prescribed
nonzero determinant, each entry vanishes with probability \(1/(p+1)\).
Exactly two opposite-entry pairs can vanish, each with probability
\(1/[p(p+1)]\). Thus
\[
 g_n(p)=\frac4{p+1}-\frac2{p(p+1)}
 \qquad(p\equiv1\pmod3,\ p\nmid n).
 \tag{3}
\]
These probabilities do not depend on the nonzero determinant, so they are
also the conditional probabilities in either determinant squareclass.

At an inert prime, the bad set is the nonsplit torus normalizer, of size
\(2(p+1)\) inside a group of size \(p(p^2-1)\). Its norm squareclasses
are equally represented. Consequently
\[
 g_n(p)=\frac2{p(p-1)}
 \qquad(p\equiv2\pmod3,\ p\nmid n).
 \tag{4}
\]

## Singular primes

At an inert prime \(p\mid n\), a bad representation has both Eisenstein
coordinates divisible by \(p\). Scaling them by \(p\) preserves angle
and every other prime test. If \(e=v_p(n)\), this gives
\[
 g_n(p)=\begin{cases}
  0,&e=1,\\
  \sigma(p^{e-2})/\sigma(p^e),&e\ge2.
 \end{cases}
 \tag{5}
\]

For a split prime choose \(p=\pi\bar\pi\), and write \(n=p^em\)
with \(p\nmid m\). Carry along any angular weight and the other prime
tests appearing in the following recursion. Let \(U_e\) be the raw count
at target \(p^em\), unrestricted at \(p\), and let \(G_e\) count its
representations in which both component norms are \(p\)-units. At the
unit target \(m\), let \(L_0\) be the count with the extra test
\[
 \ell_p(z,w)=\frac12\bigl(
  \mathbf1_{\pi\mid z}+\mathbf1_{\bar\pi\mid z}
 +\mathbf1_{\pi\mid w}+\mathbf1_{\bar\pi\mid w}\bigr).
\]
At most one component norm is divisible by \(p\), so this test lies in
\([0,1]\). It is left-unit invariant and has conditional mean \(2/(p+1)\).

If a norm has \(p\)-valuation \(v\), its \(v+1\) allocations between
\(\pi\) and \(\bar\pi\) are equally frequent. Thus the average of
\(\ell_p\) for a pair with one such divisible norm is \(v/(v+1)\).
Allocations at distinct primes are independent in Eisenstein unique
factorization. Multiplication by a factor over \(p\) preserves every
prime-divisibility test at another rational prime, including another
\(\ell_\ell\) test. The fixed \(3\)-good test and oddness depend only
on the rational norms. These facts justify the same allocation average in
every iteration below.

Classify a \(p^em\)-shell representation by
\(t=\min(v_p(a),v_p(b))\), where here \(a,b\) denote its actual norms.
They are positive on the triangular support. If \(t<e\), both valuations
equal \(t\); dividing both norms by \(p^t\) gives a \(G_{e-t}\) pair,
with multiplier \((t+1)^2\). If \(t=e\), the residual target is a unit.
A pair with both residual norms units has multiplier \((e+1)^2\).
If one residual norm has valuation \(v\ge1\), the multiplier is
\[
 \frac{(e+1)(v+e+1)}{v+1}
 =(e+1)^2-e(e+1)\frac v{v+1}.
\]
Therefore the following angularly weighted identity is exact:
\[
 U_e=\sum_{t=0}^{e-1}(t+1)^2G_{e-t}
       +(e+1)^2U_0-e(e+1)L_0.
 \tag{6}
\]
The relative angle, including the original \(h\), remains unchanged.

Put
\[
 A(T)=\frac{1+T}{(1-T)^3}=\sum_{t\ge0}(t+1)^2T^t,
 \qquad C(T)=\frac{2T}{(1-T)^3}.
\]
Formal inversion of (6) gives
\[
 \sum_{e\ge1}G_eT^e=\frac{U(T)}{A(T)}-U_0+
                         \frac{C(T)}{A(T)}L_0.
\]
The coefficients of \(1/A=(1-T)^3/(1+T)\) have magnitude at most eight,
and those of \(C/A=2T/(1+T)\) at most two. Hence the bad functional
\(U_e-G_e\) is a linear combination of unrestricted counts at levels
\(p^{e-j}m\), \(j\ge1\), and the unit-target \(\ell_p\) count.
Its total coefficient variation is at most \(8e+2\). In particular no
level condition at \(p\) remains on a target still divisible by \(p\).

Insert the predicted scalar values
\(U_e=\sigma(p^e)\), \(U_0=1\), \(L_0=2/(p+1)\). The resulting
series is
\[
 \sum_{e\ge1}G_eT^e=\frac{(p-1)^2}{p+1}\frac T{1-pT}.
\]
Thus the singular split factor is
\[
 g_n(p)=1-\frac{(p-1)^2p^{e-1}}{(p+1)\sigma(p^e)}
 \qquad(p\equiv1\pmod3,\ e=v_p(n)\ge1).
 \tag{7}
\]

Successively eliminate the primes of \((d,n)\). The rewriting commutes
with all tests at other primes. Retained unit-target \(\ell_p\) tests
have mean \(2/(p+1)\), independently of the new determinant squareclass;
the bad-prime means have the same independence. Every extracted split
prime power and every extracted inert square is \(1\pmod3\), preserving
the mixed \(3\)-good mean. The scalar terms in (2) therefore factor as
\(12\sigma(n)(h/2)\prod_{p\mid d}g_n(p)\), with the optional factor
\(2/3\).

The total coefficient variation is bounded by
\[
 \prod_{\substack{p\mid d\\p^e\Vert n}}(8e+2)
 \le\tau(n)^4,
\]
because \(8e+2\le(e+1)^4\) for \(e\ge1\). Each final target divides
\(n\), and its coprime level divides \(d\). Formula (2) gives total
error \(O(d^{3/2}\tau(n)^5\sqrt n\,h^{-1/2})\).

Finally impose oddness by subtracting the raw \(n/4\)-shell when
\(4\mid n\). Scaling by two preserves angle, every test at a prime
at least five, and the condition at three. The target valuations at odd
primes do not change, and \(n/4\equiv n\pmod3\). The main factor becomes
\(X_n\), while the error has the same bound. Since
\(\tau(n)^5\ll_\varepsilon n^\varepsilon\), this proves (1).

## An elementary lower sieve

Equations (3)--(7) give a uniform \(c<1\) with \(0\le g_n(p)\le c\).
Uniformly in \(v_p(n)\), they also give
\[
 g_n(p)=\frac4p+O(p^{-2})\quad(p\equiv1\pmod3),\qquad
 g_n(p)=O(p^{-2})\quad(p\equiv2\pmod3).
\]
Fixed-progression prime harmonic Mertens therefore implies
\[
 V_n(z)=\prod_{5\le p\le z}(1-g_n(p))\asymp(\log z)^{-2},
 \tag{8}
\]
with constants independent of the target. On each nonempty block
\[
 z^{2^{-j-1}}<p\le z^{2^{-j}},\qquad p\ge5,
\]
the sum of the \(g_n(p)\) is at most a fixed \(C\), and its product
\(v_j\) of \(1-g_n(p)\) is at least a fixed \(c_0>0\).
There are \(J=O(\log\log(z+3))\) blocks.

Choose \(\xi>0\) with
\[
 \frac52\xi<\frac{1-3\delta}{2}.
\]
Take a sufficiently large fixed integer \(A_0\), and put
\[
 m_j=2\left\lceil\frac{A_0+j}{2}\right\rceil,
 \qquad S=\sum_{j\ge0}m_j2^{-j}<\infty,
 \qquad\eta_\theta=\xi/S,
 \qquad z=n^{\eta_\theta},\quad D=n^\xi.
\]
For a witness let \(U_j,L_j\) be the ordinary inclusion-exclusion sums
for no bad prime in block \(j\), truncated at degrees \(m_j\) and
\(m_j-1\), respectively. Bonferroni gives
\(U_j\ge\mathbf1_j\ge L_j\) and \(U_j\ge0\), where
\(\mathbf1_j\) indicates no bad prime in that block. Telescoping gives
the pointwise lower bound
\[
 \mathbf1_{\mathrm{all\ good}}\ge
 \prod_jU_j-\sum_j(U_j-L_j)\prod_{i\ne j}U_i.
 \tag{9}
\]
Indeed, telescope \(\prod U_j-\prod\mathbf1_j\), enlarge each other
\(\mathbf1_i\) to \(U_i\), and then use
\(U_j-\mathbf1_j\le U_j-L_j\).

Under the product Bernoulli law with bad probabilities \(g_n(p)\),
\[
 \mathbb E U_j\ge v_j\ge c_0,\qquad
 \mathbb E(U_j-L_j)\le C^{m_j}/m_j!.
\]
Choose \(A_0\) so large that
\(\sum_{j\ge0}C^{m_j}/(c_0m_j!)\le1/2\). Independence then shows
that the expectation of the right side of (9) is at least
\(\tfrac12\prod_jv_j=\tfrac12V_n(z)\).

The polynomial in (9) uses only squarefree divisors at most
\(z^S=D\). Each coefficient has magnitude at most \(J+1\), because
prime sets in different blocks are disjoint. Apply (1) term by term against
the nonnegative coordinate weights. The total remainder is
\[
 O_\varepsilon\!\left((J+1)D^{5/2}
                  n^{1/2+\varepsilon}h^{-1/2}\right)
       =o\!\left(\frac{X_n}{(\log n)^2}\right).
\]
To obtain the last equality choose \(\varepsilon\) smaller than the
strict margin in the choice of \(\xi\). Equations (8)--(9) give weighted
mass \(\gg_\theta X_n/(\log n)^2\) with no prime \(5\le p\le z\)
in either actual component norm.

## Remove the inert tail and count values

[EisCentralTail](EISENSTEIN_POWER_SHRINKING_WINDOWS.md#a-large-prime-tail-inside-the-central-window)
with \(P=z\) bounds the remaining inert-prime loss by
\[
 O\!\left(R_o(n)M(n)
       \left\{\frac hz+\frac{1+\log n}{\sqrt n}\right\}\right).
\]
Relative to \(X_n\), this is
\[
 O\!\left(M(n)\left\{n^{-\eta_\theta}
                +(1+\log n)n^{\delta-1/2}\right\}\right)
       =o((\log n)^{-2}).
\]
Here \(M(n)=n^{o(1)}\), \(\eta_\theta>0\), and \(\delta<1/3\).
The divisor exponent can be decreased independently of the preceding
positive spectral margin.

Every surviving prime factor in an \(A\)-part exceeds \(n^{\eta_\theta}\).
Positivity and the congruence modulo three give the required pure or mixed
types; in the mixed case the fixed test removes exponent at least two at
three. Swapping components preserves every sieve condition and the even
weight, so the two mixed orientations receive equal weighted mass.

Each \(A\)-part is at most \(n\), hence has at most
\(1/\eta_\theta\) prime factors with multiplicity. Its Eisenstein
coordinate multiplicity is at most \(6\,2^{1/\eta_\theta}\), using
\(e+1\le2^e\). An ordered pair therefore has at most
\[
 36\,2^{2/\eta_\theta}
\]
coordinate witnesses; the known mixed factor three has multiplicity one.
Since \(b_h\le1\), the weighted lower mass yields the claimed number
of distinct ordered value pairs.

## Compression for coprime coset tests

The general finite-test estimate (2) pays for the size of the residue space.
The particular divisibility tests admit a smaller bound.

**Proposition EisCosetCompression.** Let \(n\equiv2,4\pmod6\), and let
\(d\) be squarefree with
\((d,6n)=1\), and let \(g(d)=\prod_{p\mid d}g(p)\), where
\[
 g(p)=
 \begin{cases}
 (4p-2)/(p(p+1)),&p\equiv1\pmod3,\\
 2/(p(p-1)),&p\equiv2\pmod3.
 \end{cases}
\]
For the odd-component shell \(\Omega_n\), and the central triangle
\(b_h(t)=\max(1-|t|/h,0)\), \(0<h\le1/2\),
\[
 \sum_{\substack{q\in\Omega_n\\d\mid A_1A_2}}b_h(t_q)
   =\frac h2 R_o(n)g(d)
     +O\!\left(6^{\omega(d)}\tau(n)\sqrt n\,h^{-1/2}\right).
 \tag{C1}
\]
The constant is absolute. For mixed targets one may impose the fixed
condition \(9\nmid A_1,A_2\); the main term is multiplied by \(2/3\),
and the error still has an absolute constant. This proposition concerns
coprime Hecke indices. It makes no new assertion at primes dividing \(n\).

To prove the bound, let \(G\) be a finite group, \(\Gamma\le G\), and
\(\rho\) a unitary representation on \(V\). Give all functions on \(G\)
counting measure. The projection onto the subspace
\(F(\gamma g)=\rho(\gamma)F(g)\) is
\[
 (Pf)(g)=\frac1{|\Gamma|}\sum_{\gamma\in\Gamma}
             \rho(\gamma)^{-1}f(\gamma g).
\]
It commutes with right translations. Suppose \(R\) on this subspace
commutes with right translations and \(\|R\|\le\Lambda\). For any
subgroup \(H\le G\), coset \(gH\), and unit vector \(v\), put
\(F=P(\mathbf1_{gH}v)\). Then \(F\) and \(RF\) are right
\(H\)-invariant, and \(\|F\|\le\sqrt{|H|}\). Averaging evaluation at
the identity over \(H\) gives
\[
 \langle v,(RF)(1)\rangle
   =\frac1{|H|}\langle\mathbf1_Hv,RF\rangle,
 \qquad
 |\langle v,(RF)(1)\rangle|\le\Lambda.
 \tag{C2}
\]
Projection can be inserted on the first factor of the inner product.
No inclusion \(\Gamma\subset H\) is required. Replacing quotient
counting measure by full counting measure rescales both squared norms
by \(|\Gamma|\), so it does not change the Hecke norm bound.

For a split prime identify \(G_p=\operatorname{PGL}_2(\mathbf F_p)\).
The component norms are \(ad\) and \(-bc\), up to units, for an
invertible matrix with entries \(a,b,c,d\). Let \(B,B^-\) be the upper
and lower triangular subgroups, \(T=B\cap B^-\), and
\(w=\left(\begin{smallmatrix}0&1\\-1&0\end{smallmatrix}\right)\).
The exact bad-prime indicator is
\[
 \mathbf1_{p\mid A_1A_2}
 =\mathbf1_B+\mathbf1_{B^-}+\mathbf1_{wB}+\mathbf1_{wB^-}
                  -\mathbf1_T-\mathbf1_{wT}.
 \tag{C3}
\]
Indeed the four single events are the four entries being zero. The only
nonempty pairwise intersections are diagonal and antidiagonal matrices;
all triple intersections are empty. Since
\(|G_p|=p(p-1)(p+1)\), \(|B|=p(p-1)\), and \(|T|=p-1\), (C3) gives
the stated split density. Its coefficient absolute sum is six.

At an inert prime a component norm vanishes only when that component is
zero in \(\mathbf F_{p^2}\). The two disjoint events are an anisotropic
torus \(T_{\rm an}\), of order \(p+1\), and its Weyl coset. Their density
is \(2/[p(p-1)]\), with coefficient absolute sum two. Every subgroup just
used has surjective determinant squareclass: for split primes it contains
the diagonal torus; for inert primes use the surjective finite-field norm.
Consequently its coset has the same density in either determinant
squareclass, and every nontrivial scalar norm character averages to zero
on that subgroup.

For a squarefree product \(d\), expand these identities as products of
cosets of product subgroups in \(G_d=\prod_{p\mid d}G_p\), allowing
\(G_1=\{1\}\). Use \(G_d\times G_3\) throughout, including in the pure
case. At the fixed factor
\(G_3=(\mathcal O/3\mathcal O)^\times/\mathbf F_3^\times\),
only the central signs among the order units become scalars. Hence the archimedean unit action factors through its finite
image, also at \(d=1\). For the unrestricted test at three, extend
each coset \(gH\) to \(gH\times G_3\), a coset of \(H\times G_3\).
The total coefficient absolute sum is at most \(6^{\omega(d)}\).
Apply (C2) termwise to the coefficient-valued Hecke operator at each
positive even Legendre degree. Its norm is at most
\(\tau(n)\sqrt n\), independently of the degree and level. The zonal unit
vector at these degrees is invariant under the order units; projecting a
coset indicator therefore preserves the signed coordinate moment, by
reindexing the norm-\(n\) shell. The shell normalization contributes the
fixed factor twelve. At degree zero subtract the scalar norm-character
part first. The preceding determinant observation kills its nonconstant
terms, leaving \(12\sigma(n)g(d)\). Thus both the nonconstant moments and
the scalar error are at most
\(12\cdot6^{\omega(d)}\tau(n)\sqrt n\).

Odd filtering subtracts the \(n/4\) shell, and costs at most a factor two
in this bound. The triangle has constant Legendre coefficient \(h/2\)
and coefficient absolute sum at most \(32h^{-1/2}\), as proved in
[EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md). These facts prove (C1).
For the mixed test let \(f_3\) be the indicator of the fixed good
condition and use right invariance under \(H\times\{1\}\). The two
projected kernels in (C2) are \(P(\mathbf1_H\delta_1v)\) and
\(P(\mathbf1_{gH} f_3 v)\), where \(\delta_1\) is the indicator of the
identity in \(G_3\). Their squared norms in full counting measure are
at most \(|H|/6\) and \(36|H|\), respectively. In the first projection
the six unit translates have disjoint support at three. The second
bound follows by contraction since \(|G_3|=36\). After division by
\(|H|\), Cauchy--Schwarz costs at most \(\sqrt6\). The already
computed six-point scalar decomposition gives main-term proportion
\(2/3\). This introduces no Hecke index divisible by its principal level.

For coefficients \(|\lambda_d|\ll_\nu d^\nu\), summing the absolute
remainders in (C1) for \(d\le D\) gives, for every fixed \(\nu>0\),
\[
 O_\nu\!\left(\tau(n)\sqrt n\,h^{-1/2}D^{1+\nu}\right).
 \tag{C4}
\]
Here the bound \(6^{\omega(d)}\ll_\nu d^\nu\) permits adjustment of
\(\nu\). The factor \(\tau(n)\) is retained even when \(D=1\).
For \(h\asymp n^{-\delta}\) with \(0<h\le1/2\), \(0\le\delta<1/3\), and a fixed positive
margin \(\epsilon\), this is \(o(nh)\) throughout
\[
 1\le D\le n^{(1-3\delta)/2-\epsilon},
 \tag{C5}
\]
when the upper exponent is positive, by choosing \(\nu\) small enough
and using \(\tau(n)=n^{o(1)}\). Thus these coprime tests reach a larger
summed level than the generic finite-test estimate. The level remains
strictly below the square-root endpoint at macroscopic width.

There is a related exact identity at a split prime square. On
\(\operatorname{PGL}_2(\mathbf Z/p^2\mathbf Z)\), write
\(E_r=\mathbf1_{p^r\mid a}\) and \(D_s=\mathbf1_{p^s\mid d}\),
for \(r,s\in\{1,2\}\). Then
\[
 \mathbf1_{v_p(ad)=1}
 =E_1+D_1-E_2-D_2-2E_1D_1+E_2D_1+E_1D_2.
 \tag{C6}
\]
This follows by separating the alternatives \(v_p(a)=1,v_p(d)=0\)
and \(v_p(d)=1,v_p(a)=0\). Each single event is a coset of a congruence
subgroup. Each intersection \(E_rD_s\) is \(w\) times the subgroup
whose off-diagonal entries are divisible by the corresponding powers;
these are subgroups because \(r,s\ge1\), and invertibility forces their
diagonal entries to be units. The seven coset terms in (C6) have total
absolute coefficient eight. Therefore (C2) also removes the generic
prime-square residue-space size cost for this test. Incorporating this
improvement into the singular-prime and full cofactor sums requires
additional bookkeeping and is not asserted here. Neither (C1) nor (C6)
supplies the missing every-target modulus-12 lower bound.

## All target valuations and joint cofactor levels

The preceding compression survives the singular-prime elimination. Its
coefficients must be summed at the square-root size of each reduced shell.
This avoids a divisor loss from deep target valuations.

**Theorem EisJointCosetDistribution.** Let \(m\) be positive and even,
\(3\nmid m\). Let \(u,v\) be coprime squarefree products of split primes,
\((uv,m)=1\), and let \(\ell\) be squarefree with \((\ell,6uv)=1\).
The primes of \(\ell\) may divide \(m\), to arbitrary powers. Let
\(W(\pm1)=0\) and
\[
 W(t)=\sum_{k\ge0}c_kP_k(t),\qquad T_W=\sum_{k\ge0}|c_k|<\infty.
\]
Put \(t=(A_1-A_2)/m\). On the odd-component shell impose
\(v_p(A_1)=1\) at \(p\mid u\),
\(v_p(A_2)=1\) at \(p\mid v\), and \(\ell\mid A_1A_2\). In the mixed
class one may also impose either prescribed good orientation at three,
their union, or no further test. Write its scalar proportion as
\(s_3=1/3,2/3,1\), respectively, and use \(s_3=1\) in the pure class.
The sum \(C_W\) of \(W(t)/2^{\omega(uv)}\) over these coordinates is
\[
 C_W=c_0s_3R_o(m)\prod_{p\mid uv}\frac{t_p}{2}
                         \prod_{p\mid\ell}g_m(p)
 +O\!\left(T_W\tau(m)\sqrt m\,4^{\omega(uv)}
                    \prod_{p\mid\ell}B_p(v_p(m))\right),
 \tag{J1}
\]
with an absolute constant and \(t_p=2(p-1)^2/[p^2(p+1)]\). The densities
\(g_m(p)\) are those already proved in the unit and singular cases above.
For a split prime set
\[
 B_p(0)=6,\qquad B_p(e)=\frac8{\sqrt p-1}\quad(e\ge1).
\]
For an inert prime set
\[
 B_p(0)=2,\qquad B_p(1)=0,\qquad B_p(e)=p^{-1}\quad(e\ge2).
\]
Thus \(\prod B_p\le6^{\omega(\ell)}\). The bound has only one divisor
factor \(\tau(m)\), uniformly in all target valuations.

For the coprime terms the proof of (C2) works with arbitrary unit vectors
\(v_0,w_0\) and the directional projection of
[EisLocationWindows (3)](EISENSTEIN_MOVING_WINDOWS.md#projection-for-directional-tests).
Keep \(G_3\) even at empty variable level. The projected evaluation
kernel \(P(\mathbf1_H\delta_1w_0)\) has squared norm \(|H|/6\), while
\(P(\mathbf1_{gH}f_3v_0)\) has squared norm at most \(36|H|\).
This costs at most \(\sqrt6\) times the operator norm. No invariance of
\(v_0\) under the units is needed. All positive Legendre degrees,
including odd degrees, therefore have the same bound. At degree zero
subtract the scalar norm characters and use the conditional mean of the
original finite test. Each coprime bad-prime test costs six or two as in
(C3). Each prescribed exact prime-square valuation costs eight by (C6),
or four after its weight division. Either component's test has this
coset decomposition, by the Weyl interchange. All conditional means are
independent of determinant squareclass outside the fixed factor at three.
These facts prove (J1) when \((\ell,m)=1\).

For singular primes apply the valuation identities to the original
\(W\)-weighted sum before its Legendre expansion. Its endpoint vanishing
omits zero coordinates, so every allocation valuation is finite. In the
notation of (6), inversion gives the exact split bad count
\[
 U_e-G_e=-\sum_{j=1}^e a_jU_{e-j}-b_eL_0,
 \quad a_j=[T^j]\frac{(1-T)^3}{1+T},\quad
 b_e=[T^e]\frac{2T}{1+T}.
 \tag{J2}
\]
Here \(U_{e-j}\) has target reduced by \(p^j\) and no test at \(p\).
The last term has target reduced by \(p^e\), now a unit at \(p\), and
test
\[
 l_p=\tfrac12(\mathbf1_{\pi\mid z}+\mathbf1_{\bar\pi\mid z}
               +\mathbf1_{\pi\mid y}+\mathbf1_{\bar\pi\mid y}).
\]
Its four single-entry cosets have total coefficient absolute sum two and
mean \(2/(p+1)\). Since \(|a_1|=4\), \(|a_j|\le8\) and \(|b_e|=2\),
the errors, divided by \(\tau(m)\sqrt m\), cost at most
\[
 \sum_{j=1}^e|a_j|p^{-j/2}+2|b_e|p^{-e/2}
 \le8p^{-1/2}+\frac{8p^{-1}}{1-p^{-1/2}}
 =\frac8{\sqrt p-1}<6.
 \tag{J3}
\]
The last inequality uses \(p\ge7\). Every reduced index has divisor count
at most \(\tau(m)\). This is a bound for the coefficients weighted by
actual reduced shell sizes, rather than for their unweighted variation.
At an inert singular prime the event is empty when \(v_p(m)=1\);
otherwise multiplication of both coordinates by \(p\) reduces the target
by \(p^2\), costing \(p^{-1}\) on the error scale.

Distinct singular primes can be eliminated successively. Common scaling
preserves the ratio of the two norms, hence the same \(W\). All other
norm-divisibility predicates, exact prime-square prescriptions and the
specified three-valuations are preserved. Tests \(l_q\) at previously
eliminated split primes are also preserved: multiplying by an Eisenstein
factor over a different rational prime cannot change divisibility by the
primes over \(q\). The allocation averaging at distinct split primes
factors by unique factorization. Split prime powers and inert prime squares
are \(1\bmod3\); the mixed scalar orientation is unchanged. After all
singular elimination, expand \(W\) and apply the coprime coset estimate
to every resulting shell. The square-root index reductions multiply,
so (J3) multiplies over the singular primes. Substitution of the scalar
factors \(\sigma(p^e)\), \(2/(p+1)\) in (J2) gives the previously proved
\(g_m(p)\). Finally subtract the \(m/4\) shell when integral. Its odd-prime
valuations and densities agree and its error scale is smaller. This proves
(J1) for all valuations. No ramified-level Hecke estimate is assumed.

**Corollary EisJointCosetLevel.** Fix a target \(n\) and the cofactor
triples \((d,u,v)\) used in EisSemigroupScale, with \(d\mid n\),
\((uv,n)=1\), pairwise coprime and squarefree. Sum the absolute errors
in (J1), at \(m=n/d\), over \(duv\le B\) and eligible \(\ell\le L\),
with scalar sieve coefficients at most \(K\) in absolute value. For a
common \(W\) the total is
\[
 \ll K T_W\tau(n)\sqrt n\,
             BL\log^7(2B)\log^5(2L),\qquad B,L\ge1.
 \tag{J4}
\]
Indeed \(6^{\omega(\ell)}\le\tau_6(\ell)\) has summatory bound
\(O(L\log^5(2L))\). For fixed \(t=uv\), its assignments to the two
cofactors and their errors total at most
\(2^{\omega(t)}4^{\omega(t)}\le\tau_8(t)\), whose sum to \(B/d\)
is \(O((B/d)\log^7(2B))\). The remaining factor \(\sqrt{n/d}\)
leaves \(\sum_{d\mid n}d^{-3/2}\le\zeta(3/2)\). This proves (J4),
including either budget equal to one.

With the single budget \(duv\ell\le R\), \(R\ge1\), one instead has
\[
 \ll K T_W\tau(n)\sqrt n\,R\log^{13}(2R).
 \tag{J5}
\]
For fixed \(d\), the squarefree product \(uv\ell\) assigns each prime
to \(u,v,\ell\), with total error coefficient at most \(4+4+6=14\).
Use \(\tau_{14}\) and the same convergent \(d^{-3/2}\) sum.

For \(K=n^{o(1)}\), including the block-sieve envelope
\(K=O(\log\log(3n))\), take the moving bump with
\(T_W\ll\sqrt{M/H}\) and
\(H^3\ge n^{1+\varepsilon}M\). Taking \(B=n^\kappa,L=n^\xi\),
the bound (J4) is \(o(H/\log^A n)\) for every fixed \(A\) whenever
\(\kappa+\xi<\varepsilon/2\). For central triangles
\(h\asymp n^{-\delta}\), \(0<h\le1/2\), \(0\le\delta<1/3\),
the condition is \(\kappa+\xi<(1-3\delta)/2\). At fixed macroscopic
width this gives an all-valuation level approaching \(n^{1/2}\) from
below. It improves the available distribution estimate, without changing
the central exponent threshold, proving a prime-summand theorem, or
supplying the modulus-12 pointwise lower bound.

## Dependencies and review

The proof dependencies are
\[
 \begin{gathered}
 \text{quaternionic norm count, strong approximation, JL and Deligne}
 \Longrightarrow\text{composite-level estimate (2)},\\
 \text{(2) and split-prime allocation (6)}
 \Longrightarrow\text{uniform divisibility estimate (1)},\\
 \text{(1) and prime harmonic Mertens}
 \Longrightarrow\text{small-prime lower sieve (9)},\\
 \text{(9), EisCentralTail and bounded coordinate fibers}
 \Longrightarrow\text{EisBoundedFactors}.
 \end{gathered}
\]
The quaternionic inputs and their exact source locators are recorded in
[EisPower](EISENSTEIN_POWER_SHRINKING_WINDOWS.md) and the
[scalar spectral proof](../../unpublished/eisenstein-restricted-prime-sums/sections/prime-bounds.tex).
The prime harmonic input is
[Loughran--Matthiesen, Lemma 2.4](https://arxiv.org/html/1904.12845v2#S2).
The composite-level projection, singular-prime recurrence, lower-sieve signs,
error budget, inert tail and value fibers received independent internal proof
reviews. A separate internal reviewer reconstructed the coset projection,
scalar exceptions, six- and seven-coset identities, prime-square lifting and
coprime summed level. Independent exact finite enumerations checked (C3)
at primes 5, 7 and 11, including both determinant squareclasses, and (C6)
on all 300,000 invertible matrices modulo 25. These checks support the
elementary identities; the displayed argument proves their general forms.
The all-valuation extension and both joint-level sums received a fresh
internal adversarial reconstruction, including directional predicates and
multiple singular primes. Its exact finite replay is documented in the
[retained checker](eisenstein-restricted-prime-sums/verify_joint_coset_distribution.py).
No external specialist review, formal verification or established
priority is claimed. The result does not assert prime summands, at most two
prime factors, the endpoint \(\theta=2/3\), or a modulus-12 binary theorem.
