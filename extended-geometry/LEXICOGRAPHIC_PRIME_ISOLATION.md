# Lexicographic prime isolation in several radial coordinates

## 1. Question and conclusion

For \(q\) circular Gaussian pairs, write
\[
 U_i=Z_iW_i,\qquad
 \mathcal L_q(U_1^{n_1}\cdots U_q^{n_q})
 =\prod_{i=1}^q n_i!.
\]
After angular constant-term extraction, the one-radial-variable proof in
[`TWO_REAL_GMC_LOWER_FACE_THEOREM.md`](TWO_REAL_GMC_LOWER_FACE_THEOREM.md)
uses a componentwise lowest power \(U^n\), Frobenius at moment \(rp\), and
the divisibility of \((jp)!/(np)!\) for \(j>n\).

The proposed multi-radial extension is to replace componentwise order by
either

1. a rank-\(q\) lexicographic valuation, possibly isolated by successive
   primes; or
2. a scalar encoding
   \[
   \nu_M(\mathbf n)=n_1+Mn_2+\cdots+M^{q-1}n_q.
   \]

The algebraic half of this proposal works: a lexicographic or sufficiently
separated weight valuation selects an initial monomial from any finite
support.  The arithmetic half fails for the ordinary Gaussian functional.
At the prime-dilated moment, the \(p\)-adic order of the factorial weight
sees total radial degree, not the lexicographic vector.  Distinct rational
residue characteristics cannot be iterated, and the Gaussian functional is
invariant under permutations of the radial coordinates.

Thus lexicographic prime isolation does **not** directly extend the
one-radial lower-face theorem.  It can still yield genuine results after an
extra hypothesis aligns the filtration with factorial order **and controls
the non-Frobenius cross terms**.  The most promising unrestricted
replacement is a same-prime, higher-\(p\)-adic initial-functional analysis,
not a sequence of different primes.

## 2. Multi-radial setup

Let \(T^\alpha\), \(\alpha\in\mathbb Z^q\), denote angular Laurent
monomials and write
\[
 P=\sum_{\alpha\in S}T^\alpha B_\alpha(\mathbf U).
\]
Then
\[
 \mathbb E(P^m)
 =\mathcal L_q\bigl(\operatorname{CT}_{\mathbf T}P^m\bigr).
\]
The Duistermaat--van der Kallen theorem still handles the angular part: if
zero lies in the convex hull of an angular face, some power of its face
polynomial has a nonzero angular constant term.  The difficulty begins
after that extraction.  In several radial variables the resulting
polynomial
\[
 F_r(\mathbf U)=\operatorname{CT}_{\mathbf T}(P^r)
\]
may have several incomparable Pareto-minimal exponent vectors and no
componentwise least vector.

Higher-rank valuations are the correct language for selecting one of these
vectors.  Weight-matrix valuations and lexicographically ordered value
groups are standard; see Kaveh--Manon and Amini--Iriarte in the references.
But a useful initial term must also be visible to \(\mathcal L_q\).

## 3. What a prime dilation actually measures

Let
\[
 w_{\mathbf n}(p)
 =\mathcal L_q(\mathbf U^{p\mathbf n})
 =\prod_i(pn_i)!.
\]
If \(p>\max_i n_i\), Legendre's formula gives
\[
 v_p(w_{\mathbf n}(p))
 =\sum_i v_p((pn_i)!)
 =\sum_i n_i
 =|\mathbf n|_1.                                      \tag{3.1}
\]
More precisely,
\[
 \frac{\prod_i(pn_i)!}{p^{|\mathbf n|_1}}
 \equiv
 (-1)^{|\mathbf n|_1}\prod_i n_i!
 \pmod p.                                             \tag{3.2}
\]
Indeed, after removing the \(n_i\) multiples of \(p\), each of the \(n_i\)
blocks contributes \((p-1)!\equiv-1\pmod p\), while the removed multiples
contribute \(n_i!\).

Consequences:

- the first \(p\)-adic order is total degree, not \(\nu_M\);
- vectors with the same total degree tie, including every coordinate
  permutation;
- the first normalized residue on a total-degree face is the ordinary
  multivariate factorial functional on that face, not evaluation at a
  lexicographic initial monomial.

For example, with \(M>1\),
\[
 \nu_M(1,1)=1+M<2M=\nu_M(0,2),
\]
but for every sufficiently large prime \(p\),
\[
 v_p((p!)^2)=2=v_p((2p)!).
\]
The scalar order has separated the vectors while the Gaussian factorial
order has not.

The same issue persists at a prime power.  If \(p>\max_i n_i\), then
\[
 v_p\!\left(\prod_i(p^k n_i)!\right)
 =|\mathbf n|_1(1+p+\cdots+p^{k-1}),                 \tag{3.3}
\]
which again forgets the coordinate distribution.

There is also a finite-support caveat to literal mixed-radix encoding.
The formula \(\nu_M\) reproduces lexicographic order only while the relevant
digits stay in a box of width less than \(M\).  This is harmless for one
fixed \(F_r\), but the supports of \(F_{rp}\) grow with \(p\).  A fixed
\(M\) is not a digit encoding on all prime-dilated supports.  Treating
\(\nu_M\) simply as a linear weight avoids that combinatorial issue, but
does not fix the factorial incompatibility.

### 3.1 The hybrid radial pair at the real \(d=3\) boundary

Three real Gaussian variables give one circular pair and one real variable:
\[
 U=ZW,\qquad V=T^2.
\]
The radial functional is
\[
 \mathcal L_{1,\mathrm{real}}(U^nV^k)
 =n!(2k-1)!!.
\]
For an odd prime \(p>\max(n,2k)\),
\[
\begin{aligned}
 v_p\bigl((pn)!\bigr)&=n,\\
 v_p\bigl((2pk-1)!!\bigr)
 &=v_p((2pk)!)-v_p((pk)!)=2k-k=k.
\end{aligned}
\]
Hence
\[
 v_p\!\left(\mathcal L_{1,\mathrm{real}}
        (U^{pn}V^{pk})\right)=n+k.                  \tag{3.4}
\]
The prime again sees total hybrid radial degree, not a lexicographic order.
Thus the mismatch is already present at the real \(d=3\) boundary; it is
not an artifact of moving to two full circular pairs in four variables.

## 4. Why two residue characteristics do not compose

Taking the moment at \(rp_1p_2\) gives separate congruences
\[
 F_{rp_1p_2}\equiv F_{rp_2}^{p_1}\pmod{p_1},
 \qquad
 F_{rp_1p_2}\equiv F_{rp_1}^{p_2}\pmod{p_2}.
\]
Neither congruence is an iteration of the other: the first contains the
unknown polynomial \(F_{rp_2}\), and the second contains \(F_{rp_1}\).
After reduction modulo \(p_1\), the integer \(p_2\) is a unit.  Algebraically,
\[
 (p_1,p_2)=(1)\subset\mathbb Z
\]
for distinct rational primes, so there is no chain of prime ideals giving
successive residue characteristics \(p_1\) and \(p_2\).

One can put several ordinary valuations on the same rational integer, but
they remain coordinate-symmetric:
\[
 v_\ell\!\left(\prod_i(Nn_i)!\right)
 =\sum_i\sum_{a\ge1}\left\lfloor\frac{Nn_i}{\ell^a}\right\rfloor.
                                                               \tag{4.1}
\]
Such a collection may distinguish some unordered partitions of the
coordinates, but it can never label the coordinates, and it does not
combine with Frobenius into a sequential lexicographic isolation.

## 5. The existing counterexample realizes the obstruction

Long's four-real-variable counterexample is
\[
 P_4=(1+Z_2)\bigl(W_1(1-Z_1)+W_2\bigr).
\]
Its angular-weight-zero part already contains
\[
 -U_1+U_2.                                          \tag{5.1}
\]
The exponent vectors \((1,0)\) and \((0,1)\) have no componentwise least
element.  Any non-symmetric lexicographic encoding selects one of them, but
\[
 \mathcal L_2(-U_1^N+U_2^N)=-N!+N!=0               \tag{5.2}
\]
for every \(N\), hence at \(N=p\), \(p_1p_2\), or any other dilation.

This is not merely a limitation of a prospective proof.  The full
polynomial satisfies
\[
 \mathbb E(P_4^m)=0\quad(m\ge1),
\qquad
 \mathbb E(Z_2P_4^m)=m!\ne0.
\]
Therefore no universal multi-radial isolation theorem can exclude all such
faces.  Any positive theorem beyond two real variables must state a
structural hypothesis that rules out this determinant/Lagrange-cancellation
mechanism.

The three-real-variable counterexample realizes the hybrid version.  With
\(V=T^2\),
\[
 P_3=(1+Z)\left(W-\frac{2+Z}{2}V\right)
\]
has weight-zero radial face
\[
 U-V.                                               \tag{5.3}
\]
Its vectors \((1,0)\) and \((0,1)\) are incomparable and have the same
prime-dilated order by (3.4).  The remaining angular terms supply the exact
all-order cancellation.  Therefore a successful \(d=3\) subclass theorem
must exclude or classify this hybrid Lagrange-cancellation face.

## 6. Positive statements that remain available

### 6.1 Componentwise lower orthant

If, after angular extraction, one has
\[
 F_r(\mathbf U)
 =c\mathbf U^{\mathbf n_0}
 +\sum_{\mathbf n>\mathbf n_0}c_{\mathbf n}\mathbf U^{\mathbf n},
\]
where every other exponent satisfies
\(\mathbf n\ge\mathbf n_0\) componentwise and at least one inequality is
strict, then the original argument extends verbatim:
\[
 \frac{\prod_i(pn_i)!}{\prod_i(pn_{0,i})!}
\]
is an integer divisible by \(p\).  This is the known strong hypothesis.

### 6.2 Why lowest total degree is not yet enough

Let \(H\) be the lowest total-radial-degree homogeneous part of \(F_r\), of
degree \(s\).  On the literal Frobenius terms in \(F_r^p\), equations
(3.1)--(3.2) give the candidate first residue
\[
 (-1)^s\mathcal L_q(H)^p\pmod p.                    \tag{6.1}
\]
But in characteristic zero one only has
\[
 F_{rp}=F_r(\mathbf U)^p+pG_p(\mathbf U).            \tag{6.2}
\]
An exponent vector \(\mathbf d\) of \(G_p\) may have
\(|\mathbf d|_1\ge sp\) while
\[
 \sum_i v_p(d_i!)<s.
\]
Indeed,
\[
 \sum_i\left\lfloor\frac{d_i}{p}\right\rfloor
 \ge s-q+1
\]
is the best bound supplied by total degree alone.  The extra coefficient
\(p\) in (6.2) does not uniformly repair the loss when \(q>1\).  Such cross
terms can occur at or below the candidate Frobenius order.

For a concrete example, take \(q=3\), \(p=7\), and
\[
 H=U_1^3+U_2^3+U_3^3.
\]
The pure Frobenius monomials \(U_i^{21}\) have factorial \(7\)-adic order
three.  The cross monomial obtained with multiplicities \((2,2,3)\) is
\[
 \frac{7!}{2!2!3!}U_1^6U_2^6U_3^9
 =210U_1^6U_2^6U_3^9.
\]
Its coefficient has \(7\)-adic order one and its factorial weight
\((6!)^2 9!\) has order one, so its total score is two, **below** the
Frobenius score three.

Therefore the attractive condition
\(\mathcal L_q(H)\ne0\) is **not by itself a theorem**.  It becomes a valid
certificate only after checking that every term of \(pG_p\) has strictly
larger factorial \(p\)-adic score than the Frobenius face.

For a term \(a_{\mathbf d}\mathbf U^{\mathbf d}\), define
\[
 \sigma_p(a_{\mathbf d},\mathbf d)
 =v_p(a_{\mathbf d})+\sum_i v_p(d_i!).              \tag{6.3}
\]
A genuine broader prime certificate is:

1. all terms of \(F_{rp}\) have score at least \(s\);
2. the score-\(s\) sum has nonzero normalized residue modulo \(p\).

This condition can hold without a componentwise least exponent, but it is
a statement about the full prime-dilated polynomial, not just its Newton
face.

Finding Newton-support hypotheses that imply this score condition is the
cleanest immediate higher-dimensional theorem target.

### 6.3 Anisotropically stretched exponent lattices

Suppose radial exponents lie in
\[
 a_1\mathbb N\times\cdots\times a_q\mathbb N,
\]
and write \(U_i^{a_i n_i}=V_i^{n_i}\).  The induced functional is
\[
 \mathcal L_{\mathbf a}(\mathbf V^{\mathbf n})
 =\prod_i(a_i n_i)!,
\]
so for large \(p\),
\[
 v_p\bigl(\mathcal L_{\mathbf a}
          (\mathbf V^{p\mathbf n})\bigr)
 =\sum_i a_i n_i.                                  \tag{6.4}
\]
Choosing \(a_i=1,M,M^2,\ldots\) makes the actual Gaussian factorial order
reflect the proposed mixed-radix degree.  This gives a real theorem engine
for power-stretched supports, or for polynomials obtained by the monomial
substitution
\[
 Z_i\mapsto Z_i^{a_i},\qquad W_i\mapsto W_i^{a_i}.
\]
It is not a theorem for arbitrary polynomials: the stretching changes the
Gaussian moment problem, and pure-moment vanishing is not preserved by this
substitution.  Moreover, (6.4) describes the pure Frobenius monomials only.
A complete stretched-lattice theorem must still show that the cross terms
in (6.2) have larger score.  One especially natural experiment is to take
the \(a_i\) to be widely separated powers of the **same** good prime; then
Legendre's formula itself supplies a mixed-radix scale, without attempting
to compose different residue characteristics.

### 6.4 Parameter-refined factorial functionals

The deformation
\[
 \mathcal L_{\mathbf t}(\mathbf U^{\mathbf n})
 =\mathbf t^{\mathbf n}\prod_i n_i!
\]
does separate coordinate degrees.  A lexicographic specialization such as
\(t_i=p^{M^{i-1}}\) would implement the desired encoding.  But the Gaussian
hypothesis supplies vanishing only at
\(\mathbf t=(1,\ldots,1)\), not an identity in \(\mathbf t\).  This route is
valid only for families whose moment vanishing is known uniformly in the
variance parameters.

## 7. Best next attack: higher \(p\)-adic initial functionals

After ordering every term by the score (6.3), one can retain the same prime
and examine successive normalized \(p\)-adic coefficients.  Schematically,
\[
 p^{-s}\mathcal L_q\bigl(H(\mathbf U^p)\bigr)
 \pmod{p^2},\pmod{p^3},\ldots                       \tag{7.1}
\]
produces secondary linear functionals involving the \(p\)-adic expansion
of factorials.  These refinements can distinguish some factorial-null faces
without pretending that the functional is lexicographic.  Mellit--Vlasenko's
Dwork congruences for constant terms are relevant technology for organizing
prime-power Frobenius corrections, although their theorem does not by
itself supply the required factorial-functional statement.

This program has a precise stopping obstruction: every functional derived
only from \(\mathcal L_q\) and rational primes remains invariant under
permuting radial coordinates.  Hence it can never separate (5.1).  A useful
theorem should instead classify which factorial-null initial faces survive
all \(p\)-adic jets and show that the survivors have a determinant/Lagrange
structure resembling the known counterexamples.

## 8. Recommended theorem targets

The proposal should be split into three targets of increasing depth.

1. **Factorial-initial certificate.** Compute the score (6.3) on
   \(F_{rp}\), including the non-Frobenius cross terms, and formulate
   checkable hypotheses under which the minimal-score residue is nonzero.
   Search for Newton-support conditions weaker than a componentwise lower
   orthant that imply those hypotheses.
2. **Stretched-lattice theorem.** Formulate a lower-face theorem for
   \(\mathcal L_{\mathbf a}\), preferably with \(a_i\) chosen as widely
   separated powers of one good prime.  Prove rather than assume the
   required cross-term score gap, and state clearly that this concerns an
   anisotropic support class.
3. **Factorial-null face classification.** Compute the first two or three
   normalized \(p\)-adic jet functionals and classify their common kernel
   on bounded Newton faces, modulo coordinate permutations.  Test whether
   the persistent kernel is generated by differences such as
   \(U_i-U_j\) and by the determinant identities arising from multivariate
   Lagrange--Good inversion.

The third target pushes most directly against the \(d=3\) boundary while
respecting the known counterexamples.  It seeks a structural classification
of the exceptional higher-dimensional subclasses, not an impossible
universal positive theorem.

## 9. References

- J. J. Duistermaat and W. van der Kallen,
  [*Constant terms in powers of a Laurent polynomial*](https://doi.org/10.1016/S0019-3577(98)80020-7),
  Indag. Math. 9 (1998), 221--231.
- H. Derksen, A. van den Essen, and W. Zhao,
  [*The Gaussian Moments Conjecture and the Jacobian Conjecture*](https://arxiv.org/abs/1506.05192),
  Israel J. Math. 219 (2017), 917--928.
- C. D. Long,
  [*Small Counterexamples to the Gaussian Moments Conjecture*](https://arxiv.org/abs/2607.18186)
  (2026), especially the four-variable example and its
  Lagrange--Good interpretation.
- K. Kaveh and C. Manon,
  [*Khovanskii bases, higher rank valuations and tropical geometry*](https://arxiv.org/abs/1610.00298),
  for weight-matrix and higher-rank valuations.
- O. Amini and H. Iriarte,
  [*Geometry of higher rank valuations*](https://arxiv.org/abs/2208.06237),
  for quasi-monomial valuations with lexicographically ordered value groups.
- A. Mellit and M. Vlasenko,
  [*Dwork's congruences for the constant terms of powers of a Laurent polynomial*](https://arxiv.org/abs/1306.5811),
  for prime-power constant-term congruences and ghost-term methods.
