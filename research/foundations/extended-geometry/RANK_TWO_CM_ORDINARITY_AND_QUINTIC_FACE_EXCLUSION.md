# Ordinary reduction with two balanced pairs and exclusion of two quintic faces

The two supported quadratic-face families (I) and (II) in
[effective prime-face elimination](EFFECTIVE_PRIME_FACE_PURE_PARAMETER_ELIMINATION.md#6-the-concrete-quintic-outputs)
have empty pure parameter sets. More strongly, every phase in either family has
infinitely many nonzero moments, over every characteristic-zero field. The proof
uses the Frobenius traces at both places of a nonordinary rank-two conjugate
pair. This supplies two powers of the residue characteristic in the trace norm.

## Ordinary primes for rank-two multiplication

**Proposition.** Let \(J\) be an abelian variety over a number field with
geometric multiplication by a CM field \(E\) of degree \(2d=\dim J\).
Assume some polarization induces complex conjugation on \(E\) by its Rosati
involution. Suppose at most two conjugate pairs of embeddings of \(E\) have
Hodge multiplicities \((1,1)\), with all other pairs of type \((2,0)\) up to
order. After a finite extension, the good ordinary primes have natural
density one.

Pass to a number field \(L\) over which the endomorphisms and polarization
are defined, which contains a normal closure of \(E\), and over which the
Zariski closure \(G\) of the action on
\(V=H^1_{\mathrm{et}}(J_{\overline L},\mathbf Q_\ell)\) is connected. These
are finite extensions. Ignore bad primes, the prime \(\ell\), and primes
where the endomorphism order differs locally from the maximal order.
Absolute degree-one primes of \(L\) have density one. At such a prime
\(\mathfrak p\) of residue characteristic \(p\), \(p\) splits completely in
\(E\).

Write \(\tau_{\mathfrak p}\in\mathcal O_E\) for the trace of Frobenius
relative to \(E\) on \(H^1\). Here is one justification of compatibility
between its ℓ-adic and crystalline descriptions. For each endomorphism \(s\)
in a rational basis of \(E\) taken from the integral endomorphism order, the
usual trace of \(s\pi_{\mathfrak p}\) on \(H^1\) is a rational integer, the
same in both cohomologies. Indeed \(s\pi_{\mathfrak p}\) is an endomorphism
of the reduction. Its characteristic polynomial is the degree polynomial of
this endomorphism on an abelian variety, in both cohomologies. The trace
pairing of \(E/\mathbf Q\) then identifies the relative traces with the same
element of \(E\), by the identities
\[
 \operatorname{Tr}(s\pi_{\mathfrak p}\mid H^1)
  = \operatorname{Tr}_{E/\mathbf Q}(s\tau_{\mathfrak p}).
\]
Every conjugate of \(\tau_{\mathfrak p}\) is the sum of two Frobenius
eigenvalues. They are algebraic integers of complex absolute value
\(\sqrt p\), which proves integrality and gives
\[
 \left|N_{E/\mathbf Q}(\tau_{\mathfrak p})\right|\le 2^{2d}p^d.       \tag{1}
\]

For the divisibility, decompose the \(p\)-divisible group, up to isogeny,
under the split \(E\otimes\mathbf Q_p\)-action. The crystalline components
have rank two. The slopes are in \([0,1]\), with sum given by the component's
Hodge multiplicity of slope one. Polarization pairs conjugate components
with reciprocal slopes \(s\mapsto 1-s\). These statements also follow
from the height, dimension and duality formulas for the split factors of the
\(p\)-divisible group. At a pair of type \((2,0)\), one component therefore
has both slopes one. The trace at that place of \(E\) is divisible by \(p\).

At a balanced pair, if one component fails to be ordinary, its two slopes
have sum one and both belong to \((0,1)\). Its conjugate likewise has both
slopes strictly positive. The trace belongs to \(\mathbf Z_p\) in each
component, so positivity of both valuations makes each trace divisible by
\(p\). If there are \(b\) balanced pairs and the reduction is nonordinary,
one such pair must fail, and the places just counted are all distinct. Thus
\[
 p^{d-b+2}\mid N_{E/\mathbf Q}(\tau_{\mathfrak p}).                  \tag{2}
\]
For \(b\le2\), (1)–(2) put the normalized value
\(N_{E/\mathbf Q}(\tau_{\mathfrak p})/p^d\) in a fixed finite set of
integers at every nonordinary degree-one prime. If \(b=0\) all components
are already ordinary.

Over \(\overline{\mathbf Q}_\ell\) write \(V=\bigoplus_\sigma V_\sigma\)
for the \(2d\) rank-two components. Let \(\chi\) be the similitude character,
with \(\chi(\operatorname{Frob}_{\mathfrak p})=p\) in the action with
weight-one Frobenius. The function
\[
 B(g)=\chi(g)^{-d}\prod_\sigma\operatorname{tr}(g\mid V_\sigma)          \tag{3}
\]
is regular and conjugacy invariant. Its value at Frobenius is the normalized
trace norm. Faltings's semisimplicity theorem makes \(G\) reductive.
Its derived group has image either trivial or \(\mathrm{SL}_2\) on each
component, since these are the only connected semisimple subgroups of
\(\mathrm{SL}_2\). If the derived group acts nontrivially on at least one
component, take a maximal torus of the derived group and a cocharacter
outside the kernels of all its nontrivial component weights. Along this
cocharacter (3) is a product of \(u^{n_\sigma}+u^{-n_\sigma}\), with at least
one \(n_\sigma\ne0\), because \(\chi\) is trivial on the derived group.
Its highest power is positive and has positive coefficient, proving
nonconstancy.

If the derived group acts trivially on every component, \(G\) is a torus.
Over \(\overline{\mathbf Q}_\ell\) the commuting algebra has one full
matrix block for each weight space, with the sum of the sizes equal to
\(2\dim J\). By Faltings's endomorphism theorem this commuting algebra is
the scalar extension of \(\operatorname{End}^0_L J\) (or its opposite on
cohomology). Write the rational algebra as a product of matrix algebras
over division algebras. Taking the product along the diagonals of maximal
subfields of the division algebras gives a commutative étale subalgebra
of rational dimension \(2\dim J\). Each rational matrix-algebra block
contributes its split matrix size times the degree of its center.
After adjoining the normal closures of these field factors, degree-one
good primes split in this algebra. The action on \(H^1\) is faithful, so
all split factors have rank one by dimension. Their \(p\)-divisible groups
have height one and slopes zero or one. Thus the torus case is ordinary on
a density-one set.

In the nonconstant case the inverse image of the finite set in
(1)–(2) is consequently a proper conjugacy-invariant closed set.
The algebraic Chebotarev argument gives density zero for primes with
Frobenius in this set. Consequently the nonordinary primes have density zero.
This uses the connected-component form of the trace argument in
[Sawin, Theorem 1 proof](https://arxiv.org/html/1506.04784v2#S0.Thmstatement1)
and [Cantoral Farfán–Li–Mantovan–Pries–Tang, §4.1](https://arxiv.org/html/2508.11174v1#S4.SS1).
The trace compatibility and norm are also the relative-trace construction
in §3.2 of the latter paper. The new divisibility (2) uses the rank-two
hypothesis at a split prime and permits two balanced pairs. No use of the
simple-signature hypothesis or of connectedness over the original ground
field is needed.

## The two quadratic-face families

Use the notation and supports of [prime-face parameter rigidity](PRIME_FACE_PARAMETER_RIGIDITY.md#3-two-residual-faces-with-a-common-quadratic-residue-polynomial).
In both families \(BC\ne0\), so normalize \(B=C=1\) over a finite algebraic
extension of the field of coefficients. Write \(\lambda\) for the resulting
coefficient \(A\). In particular \(\lambda^2=A^2/(BC)\) in original
coordinates. Put
\[
 R_\lambda(z)=1-\tfrac{\lambda}{2}z-\tfrac18 z^2.
\]
The top interval terms at order \(2k\) have respective heights \(7k,9k\).
Their angular constant coefficient is \([z^k]R_\lambda(z)^{2k}\).
For \(p=\ell k+1\), with \(\ell=7,9\) respectively, this is the unique
height whose interval denominator is divisible by \(p\). Thus, at a good
prime of absolute degree one,
\[
 p\Lambda(f^{2k})\equiv [z^k]R_\lambda(z)^{2k}\pmod{\mathfrak p}.       \tag{4}
\]
This holds with every lower coefficient present. It uses the supporting
weights specified in the two family definitions.

If \(\lambda=0\), choose \(p\equiv1\pmod{2\ell}\). Then \(k\) is even
and (4) is \(\binom{2k}{k/2}(-1/8)^{k/2}\). If \(\lambda^2=-2\),
the quadratic is \((1-\lambda z/4)^2\), and (4) is
\(\binom{4k}{k}(-\lambda/4)^k\). Both give units at all good primes
under their stated congruences, since \(4k<p\). These account for the
missing middle coefficient and zero discriminant within \(BC\ne0\).

For all other values take the smooth projective normalization of
\[
 w^\ell=z^{\ell-1}R_\lambda(z)^2 .                                  \tag{5}
\]
The roots of the quadratic are distinct and nonzero. The branch exponents
in the raw homogeneous degree \(2\ell\) form are \((6,2,2,4)\) and
\((8,2,2,6)\). The character-one eigensheaf of the normalized covering
algebra is \(\mathcal O(-2)\), since no local branch exponent reaches
\(\ell\). The Hasse–Witt scalar on its Cech generator \(w/z\) is
\([z^{p-1}](z^{\ell-1}R_\lambda(z)^2)^k=[z^k]R_\lambda(z)^{2k}\).
For \(\ell=7\) this line is in the whole Jacobian. For \(\ell=9\) it
is in the primitive-order-nine factor. Those two abelian sixfolds have
\(E=\mathbf Q(\zeta_7)\) or \(\mathbf Q(\zeta_9)\) and Hodge signature
\((1,1,2,0,1,1)\) in the primitive character order. The principal
polarization makes each deck automorphism Rosati-unitary, since it is
induced by a curve automorphism. This gives complex conjugation on \(E\)
also on the polarized order-nine isogeny factor. Exactly two conjugate
pairs are balanced, so the proposition gives ordinary primes of the
sixfold after a finite extension. At these primes (4) is nonzero. The
elliptic quotient in the degree-nine case does not enter the selected line.

For an arbitrary characteristic-zero field, specialize a single finite-type
coefficient algebra containing the required inverses to a number field.
Keep \(BC\ne0\) and, on the distinct-root chart, the discriminant nonzero.
The number-field argument supplies infinitely many nonzero specialized
moments at unbounded orders. Each such moment is then nonzero before
specialization, since its formula is a rational polynomial in the original
coefficients. This proves non-eventual-purity throughout (I) and (II),
including their degenerate face parameters. For either family the exact
squarefree exceptional polynomial for \(\lambda\) and for \(T=A^2/(BC)\)
is \(1\). This conclusion identifies the geometric exceptional set. It does
not compute a finite order \(b\) with \(1\in(M_1,\dots,M_b)\).

For \(v=2\) this closes the last mixed-winding face in the genuine stripped
quintics of [Theorem 1](STRIPPED_QUINTIC_NEWTON_FACE_REDUCTION.md#1-proven-conclusions-and-the-remaining-face).
For \(v=4,\ \deg P\le4\) it closes the last mixed-winding face in that
note's [Theorem 2](STRIPPED_QUINTIC_NEWTON_FACE_REDUCTION.md#9-a-threefold-input-excludes-the-lower-degree-even-perturbation-sector).
The \(v=4,\ \deg P=5\) sector and arbitrary full-sphere quintics require
additional arguments.

## Exact finite moment computations

The [moment generator](compute_quintic_prime_face_moments.py) expands each
power by its multinomial coefficients, retaining only angularly balanced
terms of even height. Its rational integral on each such monomial is
\[
 \Lambda(x^a y^a t^{2h})
  =\frac{2^a a!(2h)!(a+h)!}{h!(2a+2h+1)!}.
\]
The replay also generates the first four moments by independent literal
Laurent multiplication across (I), (II) and the two even charts. Run

~~~sh
python3 research/foundations/extended-geometry/compute_quintic_prime_face_moments.py
~~~

The same replay certifies that the parameter elimination ideal is zero
through moment six for (I) and through moment four for (II). In (I) retain
the six coefficients of \(1,y,\ldots,y^5\) as pivots alongside \(\lambda\)
and the two normalized face coefficients. In (II) retain the four
coefficients of \(1,y,y^2,y^2t\), set the coefficient of \(t\) to one, and
retain \(\lambda\) and the two normalized face coefficients. All other
coefficients in the respective slices are set to zero. At each order the
new pivot appears as a single linear monomial with a nonzero rational
coefficient, and none of the subsequent pivots occurs. Successive
substitution therefore gives a polynomial section over
\(\mathbb Q[\lambda]\) of each moment prefix. The intersection of the
prefix ideal with \(\mathbb Q[\lambda]\), or with the invariant parameter
ring under \(T=\lambda^2\), is zero. The all-order exceptional sets are
empty by the preceding arithmetic proof.
