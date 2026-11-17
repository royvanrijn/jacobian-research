# Eisenstein two-summand problem: spectral core and normalization audit

Working research note, 17 November 2026. This is a derivation, not a claim of independent referee acceptance.

## 1. The exact arithmetic setting

Let

\[
R=\mathbb Z[\omega],\qquad \omega^2+\omega+1=0,
\qquad Q(a,b)=a^2-ab+b^2.
\]

Let \(B=\mathbb Q(\omega)\oplus\mathbb Q(\omega)j\), with
\(j^2=-1\) and \(jz=\bar z j\). Set

\[
\mathcal O=R\oplus Rj.
\]

Quaternion conjugation sends \(z+w j\) to \(\bar z-wj\), and

\[
\operatorname{nrd}(z+w j)=z\bar z+w\bar w.
\]

The order is closed under multiplication and conjugation. Its basis
\(1,\omega,j,\omega j\) has trace-pairing blocks

\[
\begin{pmatrix}2&-1\\-1&-1\end{pmatrix},
\qquad
\begin{pmatrix}-2&1\\1&-2\end{pmatrix},
\]

whose product determinant is \(-9\). Its reduced discriminant is therefore 3.
The algebra is definite and its finite discriminant is 3; hence this order is
maximal. This last conclusion can also be obtained from the standard fact that
the algebra discriminant divides every order's reduced discriminant, together
with parity of the number of ramified places.

Every complex number lies within squared distance \(1/3\) of an Eisenstein
integer: triangulate the plane by unit equilateral triangles and use the
circumradius. Approximating both complex coordinates independently gives

\[
\forall x\in B\otimes\mathbb R\quad
\exists q\in\mathcal O:\quad\operatorname{nrd}(x-q)\le2/3<1.
\]

Thus \(\mathcal O\) is both left and right norm Euclidean. Taking a nonzero
element of smallest norm in an integral one-sided ideal and applying division
proves that every such ideal is principal.

The units are exactly

\[
\mathcal O^\times=R^\times\cup R^\times j,
\qquad |\mathcal O^\times|=12,
\]

because the two nonnegative integral component norms must sum to 1. Each unit
has reduced norm 1.

## 2. Exact representation count

For \(3\nmid N\), put

\[
\mathcal Q_N=\{\alpha\in\mathcal O:\operatorname{nrd}(\alpha)=N\},
\qquad r(N)=|\mathcal Q_N|.
\]

Then

\[
r(N)=12\sigma(N).
\]

An elementary ideal-counting proof avoids importing a theta normalization.
At every prime \(\ell\mid N\), the maximal order is \(M_2(\mathbb Z_\ell)\).
Its integral right ideals of reduced norm \(\ell^e\) correspond to sublattices
of index \(\ell^e\) in \(\mathbb Z_\ell^2\), and their number is
\(1+\ell+\cdots+\ell^e\). Localizing and recombining ideals gives exactly
\(\sigma(N)\) integral right ideals of norm \(N\). Each is principal, and
has exactly 12 integral generators of norm \(N\), differing by right units.

For \(N\equiv2\pmod3\), the two component norms in every representation are
both \(1\pmod3\), because \(Q(a,b)\equiv(a+b)^2\pmod3\) is 0 or 1.
Consequently both component norms are positive, and neither is divisible by 3.

## 3. Bad primes dividing the target

For an inert rational prime \(p\equiv2\pmod3\), including \(p=2\),

\[
p\mid Q(a,b)\quad\Longleftrightarrow\quad p\mid a,b.
\]

For odd \(p\), this follows because \(R/pR\cong\mathbb F_{p^2}\) is a field;
for 2 it also follows directly by evaluating the four residue pairs.

Let \(b_p(N)\) count representations in \(\mathcal Q_N\) in which either
component norm is divisible by \(p\). If \(p\mid N\), then one component
norm is zero modulo \(p\) precisely when both are; thus all four coordinates
are divisible by \(p\). Hence

\[
b_p(N)=
\begin{cases}
0,&p^2\nmid N,\\
12\sigma(N/p^2),&p^2\mid N.
\end{cases}
\]

Writing \(N=p^e u\), \((p,u)=1\), gives

\[
\frac{b_p(N)}{r(N)}
=\frac{p^{e-1}-1}{p^{e+1}-1}<\frac1{p^2}\qquad(e\ge2).
\]

In particular, for every target \(N\equiv2\pmod6\),
\(b_2(N)/r(N)<1/4\), and it is exactly zero when \(v_2(N)=1\).

## 4. Finite quotient and Hecke operator

Fix an odd inert prime \(p\ge5\), with \(p\nmid N\) and \(3\nmid N\).
Reduction identifies \(\mathcal O/p\mathcal O\) with \(M_2(\mathbb F_p)\),
preserving reduced norm as determinant. Set

\[
G_p=\operatorname{PGL}_2(\mathbb F_p),\quad
\Gamma_p=\rho_p(\mathcal O^\times),\quad
Y_p=\Gamma_p\backslash G_p.
\]

Only the units \(\pm1\) reduce to scalars for \(p\ge5\), so

\[
|\Gamma_p|=6,\qquad |Y_p|=p(p^2-1)/6.
\]

Use counting measure on \(Y_p\), and define

\[
(T_N f)(\Gamma_p g)
=\frac1{12}\sum_{\alpha\in\mathcal Q_N}
 f(\Gamma_p\rho_p(\alpha)g).
\tag{H}
\]

Right multiplication by units permutes \(\mathcal Q_N\), so this is
well-defined on left cosets. It commutes with the right action of \(G_p\).
Quaternion conjugation corresponds projectively to inversion and preserves
\(\mathcal Q_N\); therefore \(T_N\) is self-adjoint for counting measure.
One has

\[
T_N1=\sigma(N)1.
\]

The determinant character

\[
\chi_p(\Gamma_p g)=\left(\frac{\det g}{p}\right)
\]

is well-defined: scalars have square determinant, and all integral units have
norm 1. Directly from (H),

\[
T_N\chi_p=\left(\frac Np\right)\sigma(N)\chi_p.
\]

### Spectral lemma

On \(\{1,\chi_p\}^\perp\),

\[
\|T_N\|\le\tau(N)\sqrt N.\tag{S}
\]

Here the non-elementary input is the standard Jacquet--Langlands correspondence
and the Ramanujan bound for holomorphic weight-two cusp forms. The required
application is as follows; it does not depend on the new Euler preprint.

Let \(H=PB^\times\). Away from \(p\), use the projective images of the local
maximal-order unit groups as level. At \(p\), use the kernel of reduction
\(\operatorname{PGL}_2(\mathbb Z_p)\to G_p\). Denote this finite-adelic level
by \(U(p)\). Class number one identifies

\[
H(\mathbb Q)\backslash H(\mathbb A_f)/U(p)\cong Y_p.
\]

To check the stabilizer, suppose a rational projective quaternion is locally
represented by an order unit at every finite place. Its reduced norm has even
valuation everywhere. Rescale by a rational scalar to make all these
valuations zero. The resulting quaternion is an integral unit globally;
its action at \(p\) is exactly \(\Gamma_p\). This also explains why no extra
ramified-prime coset appears.

The unramified Hecke correspondence of index \(N\) counts integral right
ideals of norm \(N\). Principalizing those ideals gives generators
\(\alpha\), modulo right units; at the \(p\)-component the action is by
\(\rho_p(\alpha^{-1})\). Conjugating the generators changes right-unit
classes to left-unit classes, and projectively
\(\alpha^{-1}=\bar\alpha\). Thus the correspondence is exactly (H), including
its factor \(1/12\), not a spherical operator counting only primitive
generators.

The automorphic functions with trivial infinite component decompose into
characters and infinite-dimensional automorphic representations. The latter
transfer to weight-two cuspidal representations on \(\operatorname{GL}_2\).
They are unramified at every \(\ell\mid N\), since \((N,3p)=1\).
Automorphic characters factor through reduced norm. Restricted to the
finite set \(Y_p\), they depend only on determinant modulo squares, so their
span is contained in \(\langle1,\chi_p\rangle\); both displayed characters
occur already by the explicit formulas above.

For an unramified cusp representation, write its normalized Satake parameters
at \(\ell\) as \(u_\ell,v_\ell\), with modulus 1. In normalization (H), the
eigenvalue of \(T_{\ell^e}\) is

\[
\ell^{e/2}\sum_{r=0}^e u_\ell^r v_\ell^{e-r}.
\]

It has modulus at most \((e+1)\ell^{e/2}\). Multiplicativity in relatively
prime indices gives (S). In particular, \(T_{\ell^2}=T_\ell^2-\ell I\)
on this projective central-character-trivial space; its degree is
\(\ell^2+\ell+1\), matching the all-ideal count rather than only the primitive
sphere count \(\ell^2+\ell\).

## 5. Torus-normalizer averaging: exact constants

Since \(p\) is inert, put

\[
C_p=(R/pR)^\times/\mathbb F_p^\times,\qquad
K_p=C_p\cup C_p\rho_p(j).
\]

These have orders \(p+1\) and \(2(p+1)\), respectively. In contrast to the
Hurwitz Gaussian calculation, **every integral unit lies in this normalizer**:
\(\Gamma_p\subset K_p\). Define

\[
O_p=\{\Gamma_p k:k\in K_p\},\qquad f_p=1_{O_p}.
\]

Consequently

\[
|O_p|=(p+1)/3,\qquad
\mu_p=\frac{|O_p|}{|Y_p|}=\frac2{p(p-1)}.
\]

An invertible quaternion \(z+w j\) reduces projectively into \(K_p\) exactly
when \(z=0\) or \(w=0\) in \(R/pR\). This is exactly the bad-prime condition.
Because \(\Gamma_p\subset K_p\), it is constant on the full 12-element
left-unit class. Therefore

\[
b_p(N)=12(T_Nf_p)(\Gamma_p).\tag{B}
\]

The indicator \(f_p\) is invariant under right translation by \(K_p\).
Since \(T_N\) commutes with all right translations, \(T_Nf_p\) is constant
on \(O_p\), and

\[
(T_Nf_p)(\Gamma_p)=\frac{\langle f_p,T_Nf_p\rangle}{|O_p|}.\tag{A}
\]

The norm map \(\mathbb F_{p^2}^\times\to\mathbb F_p^\times\) is surjective.
Thus \(C_p\), then \(K_p\), have equally many square and nonsquare
determinants. Passing to \(\Gamma_p\)-cosets preserves this equality because
all elements of \(\Gamma_p\) have square determinant. Therefore

\[
f_p-\mu_p1\perp1,\chi_p,
\qquad
\|f_p-\mu_p1\|^2=|O_p|(1-\mu_p).
\]

Using (S) in (A) and then (B) proves the uniform estimate

\[
\boxed{
\frac{b_p(N)}{r(N)}
\le \frac2{p(p-1)}+
\frac{\tau(N)\sqrt N}{\sigma(N)}
\left(1-\frac2{p(p-1)}\right)
}\tag{U}
\]

for \(p\ge5\), \(p\equiv2\pmod3\), \(p\nmid N\), \(3\nmid N\).
For primes dividing \(N\), the exact formula in Section 3 is stronger than
the main term in (U). Thus a universally valid, slightly weaker form is

\[
\frac{b_p(N)}{r(N)}
\le \frac2{p(p-1)}+
\frac{\tau(N)\sqrt N}{\sigma(N)}
\qquad(p\ge5,\ p\equiv2\pmod3).
\]

For a target \(N\equiv2\pmod6\), any representation outside the union of
these bad sets, including the separate \(p=2\) bad set, gives the desired
two positive summands composed only of primes \(1\pmod3\).

## 6. Independent finite-group sanity check

The companion `spectral_check.py` explicitly constructs the finite quotient
for \(p=5,11\) and the degree-3 Hecke operator \(T_2\), using the three
left-unit-class representatives

\[
1+j,\quad 1+\omega j,\quad1+\omega^2j.
\]

Exact integer checks establish: quotient sizes 20 and 220; orbit sizes 2 and 4;
symmetric row-sum-3 matrices; determinant-character eigenvalue \(-3\);
balanced orbit determinant characters; and averaging identity (A).
Floating-point eigenvalue checks give largest nontrivial magnitudes 2 and
\((3+\sqrt5)/2\), both less than \(2\sqrt2\). These computations check
normalizations only; they do not prove (S) for all primes.

## 7. Sources for the standard ingredients

- J. S. Ellenberg, P. Michel, and A. Venkatesh, *Linnik's ergodic method and the
  distribution of integer points on spheres*, especially Section 8.4,
  https://www.math.ias.edu/~akshay/research/gg.pdf. The section explicitly
  describes the definite-quaternion automorphic decomposition, transfer to
  weight-two cusp forms, the Ramanujan bound, and the treatment of norm
  characters. The algebra/level in this note is different and was specified
  above rather than copied from that example.
- John Voight, *Quaternion Algebras*, Chapter 41, *Brandt matrices*,
  https://link.springer.com/chapter/10.1007/978-3-030-56694-4_41,
  for ideal-counting Hecke operators and their relationship to weight-two
  modular forms. The specific principal-level quotient here requires the
  adelic explanation in Section 4.
- Robin Chapman, *Representations of integers by the form
  x²+xy+y²+z²+zt+t²*, International Journal of Number Theory 4 (2008), 709–714,
  DOI 10.1142/S1793042108001638. The representation-count normalization can
  be cross-checked there; our count was proved by ideals in Section 2.

## Audit status

No algebraic gap has been found in this derivation. The non-elementary theorem
used is precisely the established weight-two quaternionic Hecke spectral
bound, with its local level and normalization stated above. An independent
specialist should review Section 4 before any publication claims. The
remaining work for a full theorem is an explicit tail bound, an explicit
analytic threshold, the finite computation up to that threshold, and a
literature search establishing what is already known about this modulus-3
problem.
