# Fresh radical growth for norms of several rational directions

A family of high-cardinality tuples cannot make many norm forms share
enough fresh prime support to reduce their total radical cost. The
forms take several independent rational directions as inputs. Each
fixed norm value recovers the complete tuple up to a subpower count.
Resultants determine a short list of possible shared radicals after
any one row is omitted. A projection inequality then fixes one shared
radical on a large family.

## Statement

Use the height-first convention and balanced endpoint boxes of
[the secant transfer](ABC_SECANT_ENDPOINT_TRANSFER.md#the-height-first-setting).
Thus \(n=X^L\), the height limit comes first at fixed \(L\), and
\(L\) tends to zero along the boxes. Fix \(m\ge 2\) and \(s\ge 1\).
For slot \(i\), let \(\mathcal A_i\) be a family in such a box, with
at most \(n^{1/3+o(1)}\) rows and a fixed ordered choice \(U_i,V_i\)
of two different letters. The coordinates in slot \(i\) are the
primitive homogeneous pair \((U_i,V_i)\). The third letter is
determined by their ratio and the positions of the selected letters.

For \(1\le j\le s\), let \(K_j\) have degree \(m+1\) and choose a
\(\mathbb Q\)-basis \(\theta_{j,0},\ldots,\theta_{j,m}\) of algebraic
integers. Use the integer forms
\[
 \begin{split}
 Z_j(\mathbf U,\mathbf V)
    &=\theta_{j,0}\prod_{i=1}^m V_i
      +\sum_{i=1}^m \theta_{j,i}U_i\prod_{k\ne i}V_k,\\
 F_j(\mathbf U,\mathbf V)&=N_{K_j/\mathbb Q}(Z_j(\mathbf U,\mathbf V)).
 \end{split}
\]
All basis elements have conjugates bounded by \(X^C\) for a fixed
\(C\). The fields and the bases can vary with \(X\) and with the
boxes, and are fixed across the tuple family being tested.
Assume that each pair \(F_j,F_k\), \(j\ne k\), has nonzero homogeneous
resultant in \((U_i,V_i)\) over the function field of the other slots,
for **every** slot \(i\). These resultants are integer polynomials in
the other homogeneous coordinates, of uniformly bounded degree.

Let \(\mathcal P\) be any prescribed prime pool with prime product
\(X^{O(1)}\). Write \(R_{\rm new}\) for the radical of
\(\prod_{j=1}^{s}|F_j|\) on primes outside \(\mathcal P\) and all
three primitive letters of all input rows. For each fixed
\(0<\eta<sm/3\) there is \(\delta>0\) such that
\[
  \#\bigl\{\mathbf A\in\textstyle\prod_i\mathcal A_i:
             R_{\rm new}(\mathbf A)\le n^{sm/3-\eta}\bigr\}
         \ll_{\varepsilon,L} n^{m/3-\delta}X^\varepsilon .
 \tag{TN1}
\]
This holds for sufficiently late boxes and every \(\varepsilon>0\),
uniformly in the fields, bases, and pool under the fixed height
bounds. Dependence on \(m,s,C,\eta\) and the fixed pool-height power
is allowed. In particular, a family of exponent \(m/3-o(1)\) pays
\(m/3\) per form in joint fresh radical in the limiting convention.
The count uses no hypothesis on the nonatomic masses of the input
degree measures.

## Fixed norm values and input prime supports

For any nonzero fixed integer \(M\), there are \(X^{o_X(1)}\) tuples
with \(F_j=M\), uniformly in \(K_j\) and its basis. Every conjugate of
\(Z_j\) has size at most \(X^{O(1)}\), hence at least \(X^{-O(1)}\)
by the nonzero integral norm. The principal ideals \((Z_j)\) have
norm \(|M|\), so there are at most \(\tau(|M|)^{m+1}\) of them,
uniformly in the field. For a given ideal, the generators are a
coset of units. Their logarithmic coordinate vectors occupy a box
of side \(O(\log X)\). Nonzero differences of unit log vectors have
uniform positive minimum spacing in the sup norm in bounded field
degree. Indeed a quotient unit with log coordinates in
\([-1,1]\) has all conjugates bounded by \(e\), so its monic
minimal polynomial has bounded integer coefficients. There are
only finitely many such units across these fields. The same
argument bounds the kernel of the logarithmic map. Packing the
logarithmic box yields \(O((1+\log X)^{m+1})\) generators per ideal.
For each generator the basis coefficients give
\(\prod V_i\) and \(U_i\prod_{k\ne i} V_k\). Their ratios recover all
\(U_i/V_i\) and thus all original rows.

A related elementary enumeration applies to values with bounded
radical outside a fixed pool. If \(|M|\le X^D\) and that radical is
at most \(n^\beta\), enumerate it first at cost at most \(n^\beta\).
Its support can be combined with the fixed pool. For any prime set
with product bounded by a fixed power of \(X\), the number of its
supported integers at most \(X^D\) is \(X^{o_X(1)}\), uniformly.
For any fixed \(\sigma>0\) a Rankin bound gives \(X^{D\sigma}\)
times \(\prod_p(1-p^{-\sigma})^{-1}\) over the pool. The latter is
\(X^{o_X(1)}\) since the number of primes is \(o(\log X)\).
Taking \(\sigma\) arbitrarily small proves the enumeration. Together
with the fixed norm count this bounds the tuple count with one
prescribed outside-pool radical cap by \(n^\beta X^{o_X(1)}\).

The radical common to any \(F_j\) and the original letters of the
tuple has exponent \(o(1)\) on a family of near-full exponent,
after removing power-smaller exceptional sets at each positive
tolerance. First control primes shared by two input rows using
the uniform anisotropic overlap bound (DRR4), with the full prime
support of the other rows fixed. For a prime exclusive to row \(i\),
the three possible original letters give the three boundary ratios
\((0:1),(1:0),(\chi_i:1)\), with \(\chi_i=1\) or \(-1\)
according to the signed relation. All discarded scalar coordinates
are units at the prime. By homogeneity, if the prime divides
\(F_j\), it divides the integer obtained by substituting the
corresponding boundary pair in slot \(i\). This nonzero integer
has polynomial height and depends only on the complementary
rows. Nonvanishing follows already for \(Z_j\), by independence of
the basis coefficients. For the boundary \(V_i=0\), the
\(\theta_{j,i}\) term is nonzero. For the other boundaries the
\(\theta_{j,0}\) coefficient is nonzero. Apply (DRR4) with each
of these specialized values as the prescribed prime pool. There
are a fixed number of tests, and their counts are uniform in the
complementary rows.
[Anisotropic overlap inequality (DRR4)](ABC_DIVISOR_RECORD_RECIPROCITY.md#anisotropic-intersections-for-an-arbitrary-record)

## Fixing a common radical by the omitted rows

The following projection estimate applies to any finite tuple family
\(T\) with a record \(g\) on each tuple. Write \(T_g\) for its fibres
and \(\pi_i\) for deletion of slot \(i\). Suppose each value of
\(\pi_i\) permits at most \(D_i\) records. If each slot has at most
\(N_i\) values, then
\[
        \max_g |T_g|\ \ge\
        \frac{|T|^{m}}
             {\prod_i D_i\,(\prod_i N_i)^{m-1}}.             \tag{TN2}
\]
The discrete Loomis–Whitney inequality gives
\(|T_g|^{m-1}\le\prod_i|\pi_i(T_g)|\).
For example, for a uniform random element of \(T_g\), sum the
entropy chain rules of its coordinate deletions. Omitting a
conditioning coordinate can only increase conditional entropy.
Every chain-rule summand for the full tuple then contributes at
least \(m-1\) times to the sum. Bound each projection entropy by
the logarithm of its cardinality to get the inequality.
Taking \(m\)-th roots and summing by Hölder gives
\[
 \sum_g |T_g|^{(m-1)/m}
 \le \prod_i\left(\sum_g|\pi_i(T_g)|\right)^{1/m}
 \le \left(\prod_i D_i\,(\prod_i N_i)^{m-1}\right)^{1/m}.
\]
The left side is at least
\(|T|/(\max_g|T_g|)^{1/m}\), which proves (TN2).

Apply this with \(g\) the product of all primes dividing at least
two of the evaluated forms. For each omitted slot \(i\), these
primes divide the product of the specialized pairwise resultants.
This includes common roots at infinity since the pair is primitive
and the resultants are homogeneous. Remove complementary tuples
where one of the resultants vanishes. A resultant is a nonzero
bounded-degree polynomial in the complementary rational ratios
after dehomogenizing. For finite coordinate sets of distinct
ratios, a nonzero polynomial of bounded degree has at most a
constant times the sum of their single-coordinate-deleted product
sizes as zeros. This follows by induction on the number of
coordinates, separating the zeros of a nonzero leading
coefficient from the bounded number of roots in a free coordinate.
The exceptional full tuples thus number at most
\(n^{(m-1)/3+o(1)}\).

For each remaining complement, \(g\) belongs to a divisor list of
length \(X^{o_X(1)}\), since the specialized resultants have
polynomial height. For \(s=1\) take \(g=1\). By (TN2), a family
of size \(n^{m/3-\nu}\) has a fibre with a single \(g\) of size
at least \(n^{m/3-m\nu-o(1)}\), including inner subpower allowances.
The record \(g\) itself has polynomial height.

To prove (TN1), take such a fibre among tuples satisfying the
asserted radical cap. Outside \(g\), the forms have pairwise
disjoint prime support. Restrict to tuples where each form's
radical on original row primes has exponent at most a fixed
\(\gamma<\eta/(4s)\), using the overlap bounds above before taking
the fibre. Some form must have fresh radical outside \(g\) at
most \(n^{m/3-\eta/s}\), because the product of their disjoint
fresh radicals there is at most \(n^{sm/3-\eta}\). For this
form its radical outside the fixed pool \(\mathcal P\cup
\{p:p\mid g\}\) is therefore at most \(n^{m/3-\eta/s+\gamma}\).
The fixed-norm enumeration bounds the fibre, after summing over
the \(s\) choices, by \(n^{m/3-\eta/s+\gamma}X^{o_X(1)}\).
This contradicts its lower bound for sufficiently small positive
\(\nu\). Choose \(\nu\) below the fixed power gaps in the
original-support and resultant-zero exclusions, and choose the
inner subpower allowances below the remaining fixed margins at
each \(L\). This yields a positive \(\delta\) in (TN1) under the
stated order of limits.

## A concrete family with arbitrarily many forms

Fix \(K\) and an integral \(\mathbb Q\)-basis
\(\theta_0,\ldots,\theta_m\) as above. For form \(j\) replace
\(\theta_0\) by \(\theta_0+t_j\) with \(t_j\) an integer, keeping
the other basis elements. At most one rational choice of \(t_j\)
can break basis independence. In slot \(i\) over a normal closure,
dehomogenizing the complementary ratios \(x_k=U_k/V_k\) gives the
conjugate linear forms
\[
     \sigma(\theta_i) U_i+
       \left(\sigma(\theta_0)+t_j+
              \sum_{k\ne i}\sigma(\theta_k)x_k\right)V_i .
\]
For a fixed earlier \(t_j\), two embeddings, and a slot \(i\),
the determinant of this form with the corresponding conjugate
form at a new \(t\) is affine in \(t\), with nonzero coefficient.
At most one \(t\) can make it the zero polynomial. Avoid the
finitely many such choices at each stage. This gives the full
pairwise resultant hypothesis for every fixed \(s\), in every
slot, with no further hypothesis on the rational inputs.

## Scope and antecedents

The fixed-value argument extends the repository's
[coefficient-uniform binary norm count](ABC_RATIONAL_PROJECTIVE_CARRIERS.md#linear-and-nonlinear-factor-counts-on-a-primitive-carrier).
Classical binary fixed-value inputs include Stewart's Thue bounds
and the uniform diagonal quadratic count quoted in
Browning–Verzobio. Here the basis components recover several
independent rational directions, and the shared-prime record is
fixed using resultants from all omissions.
[Stewart, Theorem 1](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf)
· [Browning–Verzobio, Lemma 2.3](https://arxiv.org/html/2608.24512v1#S2.SS1)

The counting scope differs from two recent results at fullness two.
Pohl's Theorem 1.2 is an asymptotic for primitive binary quadratic
inputs with fixed coprime squarefree coefficients
\(a,b\equiv1\pmod4\) and squarefulness outside the fixed primes
dividing \(ab\). Van Doorn's Theorem 1 constructs three-term
squarefull progressions with common difference \(2\sqrt N+1\)
by a Pell equation. Neither result supplies a small joint radical
cap for the multivariable norm evaluations in (TN1).
[Pohl, Theorem 1.2](https://arxiv.org/html/2610.20316v1#S1)
· [Van Doorn, Theorem 1](https://arxiv.org/html/2605.06697v1#S3.SS1)

For use against the unrestricted \(1/3\) endpoint, the remaining
obligation is an upper bound for the joint fresh radical on
enough tuples arising from the original additive relation, or a
replacement for that upper bound that exploits the norm counts.
