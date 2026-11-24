# Cyclic positive polynomial blocks with a symmetric common monomial

## Exact statement

Let `K` be a characteristic-zero field. Choose an integer `r>=1`. For
each `ell=1,...,r`, choose:

* a nonzero polynomial `P_ell` in finitely many variables, with
  nonnegative rational coefficients and zero constant term;
* a cycle length `q_ell>=1`, a primitive `q_ell`-th root of unity
  `zeta_ell` in `K`, and a positive rational scalar `rho_ell`;
* a vector of nonnegative integers `a_ell=(a_(ell,i))`, with one
  entry for each variable of `P_ell`.

Make `q_ell` copies `P_(ell,j)`, `0<=j<q_ell`, on variable blocks
`X_(ell,j,i)` that are pairwise disjoint across every `ell,j`. The
variables within one polynomial template need not be separated among
its monomials. Put

\[
 R=\prod_{\ell,j,i}X_{\ell,j,i}^{a_{\ell,i}},\qquad
 G=\sum_{\ell=1}^{r}\rho_\ell
           \sum_{j=0}^{q_\ell-1}\zeta_\ell^j P_{\ell,j},
 \qquad F=\lambda R G,\quad \lambda\in K^*.
 \tag{1}
\]

The common monomial overlaps the polynomial blocks. Its exponent on a
given template coordinate is required to be identical in all copies
of that template. Let

\[
 S=\langle q_1,\ldots,q_r\rangle\subseteq\mathbb N_0,
 \qquad q_*=\min_\ell q_\ell,\qquad L=\sum_\ell q_\ell.
\]

**Theorem 1.** For the ordinary factorial functional
`Lambda(X^alpha)=prod_i alpha_i!`, and for every `n>=0`,

\[
 \boxed{\Lambda(F^n)=0\quad\Longleftrightarrow\quad n\notin S.}
 \tag{2}
\]

If `n in S`, then `Lambda(F^n)/lambda^n` is a positive rational number.
Every window of `q_*` consecutive moment orders contains a nonzero
moment. The actual support of `F` has

\[
                s=\sum_\ell q_\ell\,\#\operatorname{supp}(P_\ell)
 \tag{3}
\]

terms, so this proves its support-size Strong Factorial assertion at
every starting index. The theorem is uniform in dimension, polynomial
degree, number of terms, common exponents, and moment order.

This is a constrained coefficient family, not unrestricted SFC. The
positive rational template coefficients, cyclic phases, disjoint
copies, and symmetry of the common monomial are part of the theorem.
Putting every `a_(ell,i)=0` and taking monomial templates recovers the
[retained cyclic monomial theorem](SFC_CYCLIC_MONOMIAL_BLOCKS.md).

## 1. Positive Gamma cumulants for nonnegative polynomial templates

For positive rational numbers `beta_1,...,beta_d`, define the product
Gamma functional by

\[
 \mathcal L_{\boldsymbol\beta}(X^\alpha)
       =\prod_{i=1}^{d}(\beta_i)_{\alpha_i}.
\]

For any nonzero polynomial `P` with nonnegative rational coefficients
and zero constant term, write

\[
 \phi_{\boldsymbol\beta,P}(t)
      =\sum_{k\ge0}\frac{\mathcal L_{\boldsymbol\beta}(P^k)}{k!}t^k,
 \qquad
 \log\phi_{\boldsymbol\beta,P}(t)
      =\sum_{k\ge1}\frac{\kappa_k}{k!}t^k.
 \tag{4}
\]

**Lemma 2.** Every `kappa_k` is a strictly positive rational number.
More precisely, choose any monomial `c X^d` of `P` with `c>0`, and
choose `i_*` with `d_(i_*)>0`. For every `k>=1`,

\[
 \kappa_k\ge (k-1)!c^k
       \beta_{i_*}^{\,k(d_{i_*}-1)+1}
       \prod_{i\ne i_*}\beta_i^{\,kd_i}>0.
 \tag{5}
\]

**Proof.** The elementary identity

\[
                    (\beta)_b=\sum_{\sigma\in S_b}
                                      \beta^{\#\mathrm{cycles}(\sigma)}
 \tag{6}
\]

follows by inserting the last label: it either forms a new cycle,
contributing `beta`, or is inserted after one of the `b-1` existing
labels. The empty permutation has weight one.

To expand `L_beta(P^k)`, regard its `k` factors as labelled vertices.
At each vertex choose a monomial of `P`, retaining its positive
coefficient, and attach the required number of labelled legs of each
coordinate type. For each coordinate `i`, independently choose a
permutation of all legs of that type and give it weight
`beta_i^(number of cycles)`. Equation (6) proves that the sum of all
these weights is exactly `L_beta(P^k)`.

Connect vertices whenever a permutation cycle visits their legs.
Every tuple of permutations has a unique partition of the vertices
into connected components. Its weight factors over these components.
Conversely independent connected tuples on the blocks of a vertex
partition reconstruct exactly one full tuple. Thus, if `C_k` is the
sum of weights of connected tuples on `k` labelled vertices, then

\[
 \mathcal L_{\boldsymbol\beta}(P^k)
                =\sum_{\pi\in\operatorname{Part}([k])}
                                      \prod_{B\in\pi}C_{|B|}.
 \tag{7}
\]

The formal exponential formula, or its coefficient recurrence with
the component of vertex one singled out, gives
`log phi=sum C_k t^k/k!`. Hence `kappa_k=C_k` and these numbers are
nonnegative.

For strict positivity, choose `cX^d` at every vertex. At coordinate
`i_*`, select one fixed leg from each vertex and put those `k` legs
into one cycle. There are `(k-1)!` such cycles. Make all remaining
legs fixed points, at every coordinate. The resulting tuple is
connected and has exactly the weight on the right side of (5),
including `k=1`. These tuples give the stated lower bound. QED.

All sums defining a fixed coefficient are finite. This proof makes
no analytic convergence assumption and applies to arbitrary positive
rational Gamma shapes. Kaluza log-convexity is not being assumed:
already for `P=X` and `beta=3`, the coefficients `(3)_k/k!` have
decreasing consecutive ratios, whereas `log phi=-3 log(1-t)` still
has strictly positive coefficients.

The permutation-component argument is the Gamma-weighted version of
the connected permutation expansion retained in
[UNRESTRICTED_SFC_MECHANISMS.tex](UNRESTRICTED_SFC_MECHANISMS.tex), in its subsection on the remaining
strict-nonvanishing problem. Here it is given in full and only used
before any elimination or resultant is taken.

## 2. Fix the moment order before applying the filter

Fix one integer `n>=0` and put

\[
 \beta_{\ell,i}(n)=na_{\ell,i}+1,\qquad
 C_n=\Lambda(R^n)=\prod_{\ell,j,i}(na_{\ell,i})!.
 \tag{8}
\]

For every polynomial `H`, coefficientwise factorial cancellation gives

\[
          \Lambda(R^nH)=C_n\mathcal L_{\boldsymbol\beta(n)}(H).
 \tag{9}
\]

The Gamma functional in (9) remains a product over all copies. Its
shape vector is identical in all copies of the same cycle. Define
`phi_(ell,n)` by (4), using `P_ell` and these shapes. The multinomial
theorem and disjointness give the coefficient identity

\[
 \frac{\Lambda(F^n)}{\lambda^n n!C_n}
       =[t^n]\prod_{\ell=1}^{r}\prod_{j=0}^{q_\ell-1}
                      \phi_{\ell,n}(\rho_\ell\zeta_\ell^jt).
 \tag{10}
\]

This fixes `n` first. The kernels `phi_(ell,n)` depend on that chosen
moment order; (10) is not asserted to be one generating-function
identity after summing over all `n`.

Write `log phi_(ell,n)=sum_(k>=1) b_(ell,n,k)t^k`. Every coefficient
is positive rational by Lemma 2. Root-of-unity filtering transforms
the product in (10) into

\[
 \exp\!\left(
    \sum_{\ell=1}^{r}q_\ell\sum_{k\ge1}
          b_{\ell,n,q_\ell k}\rho_\ell^{q_\ell k}t^{q_\ell k}
       \right).
 \tag{11}
\]

Every displayed coefficient in the exponent is positive. Its support
contains each generator `q_ell`, and is contained in their semigroup.
Consequently the support of its exponential is exactly `S`: every
element of `S` has at least one positive contribution, and no term
outside `S` can occur. Extracting the chosen coefficient `t^n`
proves (2). The scalars `n! C_n` are positive integers.

Multiples of `q_*` lie in `S`, so every window of that length contains
a nonzero moment. Because the templates have zero constant terms,
monomials from different copies are distinct. Multiplying by the
common monomial translates all exponent vectors and preserves their
distinctness. This proves (3) and the SFC implication.

The argument is an identity over the prime field with roots of unity
adjoined; all surviving normalized coefficients are rational. Thus it
holds in every characteristic-zero field containing the specified
roots. If the polynomial descends from such an extension to a smaller
field, the same nonvanishing holds there, since field extensions are
injective. No ordering or embedding into the reals is required for
the eventual coefficient field.

## 3. Genuine symmetric variable merging and monomial tilting

**Within-copy identifications.** Partition the coordinate variables of
one template into groups and identify the variables in each group.
Perform the identical identification in every copy of that template,
without identifying variables from different copies. The resulting
template still has nonnegative rational coefficients, is nonzero, and
has zero constant term. Coincident monomials add with positive
coefficients, so no cancellation can remove the template. The common
monomial exponents in each group add, still identically across copies.
Theorem 1 therefore gives the **same exact moment zero set** after
every such identification. The actual support may shrink; it remains
at least `L`, so the proved `q_*` window remains within its support size.

**Adding a symmetric common monomial.** Theorem 1 itself shows that
the entire zero pattern is unchanged by multiplying `G` by the
displayed overlapping common monomial `R`.

This can be realized by actual integer row mergers. First multiply
`G` by a monomial in new independent variables, assigning one such
variable with exponent `a_(ell,i)` to every copy coordinate
`(ell,j,i)`. Those independent multipliers do not change any moment
zero. Then identify each new variable with its assigned original
coordinate, in the same pattern throughout every cycle. The result
is `RG`, and (2) proves preservation at the endpoint. These row sums
can be nonuniform on the polynomial support; the proof controls the
symmetric batch, rather than asserting preservation of every
intermediate asymmetric identification.

These are unbounded-support merger classes with explicit coefficient
restrictions. They are not a proof that the generic regularity property
`R(A,n)` survives arbitrary row merging for all coefficients.

**Negative control.** Symmetry of the overlapping multiplier cannot
be dropped from the exact zero-pattern statement. The two-cycle
`G=X-Y` has every odd factorial moment zero. Multiplication by
`R=X` gives

\[
                     \Lambda(X(X-Y))=2!-1!1!=1.
\]

The predicted odd zero has disappeared. Both examples satisfy their
SFC conclusions; this is a counterexample only to the stronger
unsymmetrical zero-pattern-preservation claim.

## 4. Powers and disjoint products, with support collisions controlled

For every `m>=1`,

\[
                  \Lambda((F^m)^n)=0\iff mn\notin S.
 \tag{12}
\]

Thus the moments of `F^m` are nonzero at every multiple of
`q_*/gcd(q_*,m)`. To compare this period with its actual support,
choose in each template a generic positive rational weight vector
that has a unique maximizing monomial. Its maximum is positive
because the template has no constant term. Normalize the vector so
that this maximum is one, and use the same vector in each copy.
Combine these weight vectors on the disjoint variable blocks.

The resulting maximal face of the support of `G` consists of exactly
one selected monomial from each copy, hence `L` monomials. Their
exponent vectors are linearly independent, because their variable
blocks are disjoint and each selected vector is nonzero. The highest
weight part of `G^m` is the `m`-th power of this face polynomial.
Its `binom(L+m-1,m)` exponent vectors are distinct and their
multinomial coefficients times products of nonzero phase coefficients
cannot vanish in characteristic zero. Lower-weight terms cannot
contribute to that face. Multiplication by `R^m` translates it.
Therefore

\[
 \#\operatorname{supp}(F^m)\ge\binom{L+m-1}{m}
                  \ge L\ge q_*\ge q_*/\gcd(q_*,m).
 \tag{13}
\]

This proves SFC for every positive power, despite possible collisions
among other monomials of the expanded polynomial.

Now let `F_1,...,F_d`, `d>=1`, be phases of Theorem 1 on disjoint
variable sets, and take positive powers `m_i`. For
`H=prod_i F_i^(m_i)`, multiplicativity gives the exact condition

\[
 \Lambda(H^n)\ne0\iff m_i n\in S_i\text{ for every }i.
 \tag{14}
\]

In particular moments are nonzero at all multiples of

\[
       T=\operatorname{lcm}_i
                     \left(q_{*,i}/\gcd(q_{*,i},m_i)\right).
\]

Disjointness makes the actual support sizes multiply. By (13) their
product is at least the product of the displayed periods, and hence
at least `T`. This proves SFC for these products. An additional
nonzero monomial in independent variables multiplies every moment
by a nonzero scalar and does not alter support size or the zero set.

No assertion of arbitrary product preservation is used. Cross-copy
variable identifications are outside the proof, because they destroy
the independent Gamma factorization in (10).

## 5. Exact controls and proof boundary

The companion [verify_sfc_cyclic_polynomial_tilts.py](verify_sfc_cyclic_polynomial_tilts.py) uses exact rational and
cyclotomic arithmetic. It compares connected permutation sums with
formal cumulants in bounded examples, checks target-order Gamma
tilting against direct ordinary-factorial expansion, checks the exact
semigroup zero pattern in the stated finite domains, and verifies
the exposed-face support bound in concrete powers. It includes the
asymmetric-multiplier countercontrol and a Gamma example where
Kaluza log-convexity fails but the cumulants remain positive.

The universal statements are the written proofs above. Finite checks
do not establish the arbitrary-order or arbitrary-support quantifiers.
No general positive-kernel averaging, Gram-positivity preservation,
unrestricted row-merger theorem, or unrestricted SFC is asserted.


Replay from the repository root:

```sh
python research/factorial/image-rigidity/verify_sfc_cyclic_polynomial_tilts.py
```

The finite replay receipt is
[sfc_cyclic_polynomial_tilts_verification.json](sfc_cyclic_polynomial_tilts_verification.json).
The checker rejects a changed proof-note hash and execution with disabled
assertions. Both rejection paths are themselves exercised in the replay.
The proof above is self-contained; the links to the retained connected
permutation expansion and earlier cyclic theorem record the mathematical
provenance and overlap.

Internal mathematical review: the root agent read and reconstructed the
proof. The factorial-review agent separately checked the Gamma permutation
component factorization and strict lower bound, the target-order Gamma
tilting, exact support count, symmetric merger realization, and exposed-face
noncancellation in powers and products. Its separate exact replay also
passed. These are internal reviews and reproducibility checks, not external
peer review or formal verification. No literature-priority claim is made.
