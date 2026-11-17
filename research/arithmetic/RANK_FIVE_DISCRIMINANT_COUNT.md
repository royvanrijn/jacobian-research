# A positive discriminant density of rank-five Hasse fibre algebras

Rank-five Hasse-failing finite etale algebras over `Q` are classified by
non-Galois cubic fields whose quadratic resolvent supplies the missing
local points. Ordered by the cubic field discriminant, these fields have
a positive density among all cubic fields of either signature. This is
a count of field and fibre-algebra isomorphism classes; its height and
family differ from the target-height counts for the one fixed map `Phi`.

## 1. The local specification

Let `L` be a non-Galois cubic field, `K` its quadratic resolvent, and
write the **maximal-order field discriminants** as

\[
 D_L=d_K f^2,\qquad f\in\mathbb Z_{>0}.              \tag{1.1}
\]

The rank-five classification in
[FMA2](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md) and its
[field-discriminant form](RANK_FIVE_DISCRIMINANT_CLASSIFICATION.md),
Theorem RDF1, say that `K times L` is Hasse-failing exactly when:

- every prime `p!=3` dividing `f` satisfies `p=1 mod 3`;
- if `3` divides `f`, then `K` splits over `Q_3`.

The support of `f` is exactly the set of totally ramified cubic primes.
Consequently the allowed maximal cubic local algebras are:

| Prime | Allowed local algebras |
| --- | --- |
| `p=1 mod 3` | All etale cubic `Q_p`-algebras |
| `p=2 mod 3` | All except the totally ramified cubic fields |
| `p=3` | All nontotally ramified algebras, together with the cyclic totally ramified cubic fields |

At the real place either cubic signature is allowed. A cubic algebra
over `R` always has a real point.

The use of field discriminants is essential. Replacing `D_L` by the
discriminant of a chosen polynomial or a nonmaximal order can insert
extra square factors and change the apparent support of `f`.

## 2. Local masses, including the wild prime

For a set of maximal cubic local algebras `Sigma_p`, use the mass

\[
 m_p(\Sigma_p)=\sum_{A\in\Sigma_p}
       \frac{1}{|\operatorname{Aut}_{\mathbb Q_p}(A)|
                    \operatorname{Disc}_p(A)}.      \tag{2.1}
\]

Here `Disc_p(A)` is the positive power of `p` given by its maximal
order discriminant. The total and nontotally ramified masses are

\[
 m_p(\mathrm{all})=1+\frac1p+\frac1{p^2},\qquad
 m_p(\mathrm{nontotal})=1+\frac1p.                  \tag{2.2}
\]

These hold also at `2,3`. The exact primary inputs are Bhargava,
Shankar and Tsimerman, *On the Davenport--Heilbronn theorems and second
order terms*, Inventiones Mathematicae **193** (2013), 439--499,
Lemma 19, p. 453, and Remark 4, equation (58), p. 480:
[author-hosted paper](https://www.math.toronto.edu/~ila/BST.pdf),
[retained copy](sources/bst-2013-davenport-heilbronn.pdf).

At `3`, exactly three totally ramified cyclic cubic extensions occur,
and each has discriminant `3^4` and automorphism group of order three.
This is an exact local classification input from S. Pauli,
*Constructing class fields over local fields*, Journal de Theorie des
Nombres de Bordeaux **18** (2006), 627--652, Corollary 4.7,
pp. 638--639, and the `Q_3` table on p. 649, with the discriminant
subscript convention on p. 647:
[published primary paper](https://numdam.org/item/10.5802/jtnb.563.pdf),
[retained copy](sources/pauli-2006-local-class-fields.pdf).
The three defining Eisenstein polynomials are

\[
 X^3+6X^2+3+9k,\qquad k=0,1,2.                     \tag{2.3}
\]

Their polynomial discriminants are `-2835,-14256,-30051`, each of
three-adic valuation four and square class one over `Q_3`. The
uniformizer integral-basis statement in Milne's *Algebraic Number
Theory*, Proposition 7.55 and Remark 7.56, pp. 122--123, identifies
these as the local field discriminant exponents
([primary text](https://www.jmilne.org/math/CourseNotes/ANT300.pdf),
[retained copy](sources/milne-ant300.pdf)); the complete count and
distinctness come from Pauli's classification, not merely these finite
polynomial checks.

It follows that the allowed mass at three is

\[
 m_3(\mathrm{allowed})=\frac43+
                 3\left(\frac1{3\cdot3^4}\right)
             =\frac{109}{81},\qquad
 \frac{m_3(\mathrm{allowed})}{m_3(\mathrm{all})}
             =\frac{109}{117}.                    \tag{2.4}
\]

At every prime `p=2 mod 3`, the corresponding ratio is

\[
 \frac{1+p^{-1}}{1+p^{-1}+p^{-2}}
     =1-\frac1{p^2+p+1}
     =\frac{1-p^{-2}}{1-p^{-3}}.                   \tag{2.5}
\]

## 3. The counting theorem

**Theorem RDC1 (positive discriminant density).** Let `N_HP^+(X)`
and `N_HP^-(X)` count non-Galois cubic fields up to `Q`-isomorphism,
with `0<D_L<X` and `-X<D_L<0`, respectively, for which the rank-five
algebra `K times L` is Hasse-failing. Put

\[
 \boxed{\rho=\frac{109}{117}
       \prod_{p\equiv2\ (3)}
             \left(1-\frac1{p^2+p+1}\right)>0.}     \tag{3.1}
\]

Then

\[
 \boxed{N_{\rm HP}^+(X)\sim\frac{\rho}{12\zeta(3)}X,
 \qquad N_{\rm HP}^-(X)\sim\frac{\rho}{4\zeta(3)}X.} \tag{3.2}
\]

Thus the limiting proportion is `rho` in each signature. In particular
`0.75979<rho<0.75988`, so it is about seventy-six percent.

**Proof.** Impose the maximal local orders from section 1 and one
chosen archimedean signature. This is an acceptable infinite local
specification in Bhargava--Shankar--Tsimerman's Theorem 8, equation
(8), p. 444: at every sufficiently large prime it contains every
maximal order which is not totally ramified. This verifies the
infinite-specification hypothesis; no interchange of an unsupported
infinite product and a counting limit is used.

That theorem gives the leading coefficient

\[
 \left(\frac12\sum_{A_\infty}
             \frac1{|\operatorname{Aut}(A_\infty)|}\right)
       \prod_p(1-p^{-1})m_p(\Sigma_p).              \tag{3.3}
\]

The real factor is `1/12` for `R^3` and `1/4` for `R times C`.
For unrestricted maximal local orders,
`(1-p^{-1})(1+p^{-1}+p^{-2})=1-p^{-3}`. Dividing the local
factors by these unrestricted factors, (2.4)--(2.5) give exactly
`rho`, while their unrestricted product is `1/zeta(3)`.

The counting theorem also includes cyclic cubic fields. Removing them
does not change the leading coefficient: Lemma 22 of the same primary
paper, p. 456, bounds even all cubic orders with a nontrivial cubic
automorphism by `O_epsilon(X^{3/4+epsilon})`, which is `o(X)` for a
fixed `epsilon<1/4`. Negative-discriminant cubic fields are already
non-Galois. Thus (3.2) counts exactly the fields relevant to the
rank-five classification.

The product converges to a positive number because
`sum_p 1/(p^2+p+1)` converges and every factor is positive. For the
displayed numerical enclosure, let `rho_Y` be the partial product
through `p<=Y`, including `109/117`. The elementary product bound and
the telescoping sum `sum_{n>Y}1/(n(n+1))=1/(Y+1)` give

\[
 \rho_Y\frac{Y}{Y+1}\le\rho\le\rho_Y
       \qquad(Y\in\mathbb Z_{\ge2}).               \tag{3.4}
\]

The exact rational partial product at `Y=10000` lies between
`0.75987590` and `0.75987591`; (3.4) proves the stated wider
enclosure. \(\square\)

The local specification also has the primewise shapes in the paper's
definition of strong acceptability, p. 443. This note uses only its
leading Theorem 8; it does not compute the separate local index
integrals needed for a secondary coefficient or claim a quantitative
error term beyond `o(X)`.

**Corollary RDC2 (fibre-algebra classification and realization).** The
same counts (3.2) count rank-five Hasse-failing etale algebra
isomorphism classes `K times L`, with height `|Disc(L)|` and the
specified cubic signature. Every such algebra is realizable as a full
fibre of a determinant-one Keller map on `A^3_Q` of geometric degree
five, allowing
the map to vary with the algebra.

**Proof.** The abstract rank-five classification FMA2 identifies every
such Hasse algebra with exactly this form, and the unique field
component of degree three recovers `L` from its product algebra.
Conversely the quadratic resolvent is intrinsic to `L`. This is a
bijection on isomorphism classes. The full-fibre realization is
Corollary `cor:finite-etale-realization` of
[Common Arithmetic Fibers](../../unpublished/common-arithmetic-fibers/sections/01-keller-fibers.tex),
applied in rank five. \(\square\)

## 4. Fundamental discriminants as a consistency check

**Corollary RDC3 (the `f=1` subfamily).** Every non-Galois cubic
field with fundamental field discriminant gives a rank-five Hasse
algebra by the same construction. The counts of this subfamily are

\[
 \frac{X}{12\zeta(2)}+o(X),\qquad
 \frac{X}{4\zeta(2)}+o(X)                           \tag{4.1}
\]

in the two signatures. Its density among all cubic fields is
`zeta(3)/zeta(2)`, strictly below `rho`.

**Proof.** Fundamental `D_L` is equivalent to `f=1` in (1.1).
The local criterion is then vacuous; equivalently no prime is totally
ramified. Apply the same primary leading theorem with the nontotally
ramified maximal specification at every prime. Its local Euler factor
is `(1-p^{-1})(1+p^{-1})=1-p^{-2}`, giving (4.1).
This specification is contained in the Hasse specification, with a
strict mass increase already at three, proving the strict density
comparison. \(\square\)

Here “fundamental” does not mean squarefree: even fundamental quadratic
discriminants are included. The theorem does not decide which integers
occur as cubic field discriminants. It tests and counts actual cubic
fields. Nor is the height the full algebra discriminant:
`Disc(K times L)=d_K D_L=d_K^2 f^2`, a different quantity.

## 5. Verification and scope

Run `.venv/bin/python research/arithmetic/verify_rank_five_discriminant_count.py`.
The [checker](verify_rank_five_discriminant_count.py) verifies all mass
identities, the three explicit local polynomials and their three-adic
discriminants, the Euler-factor transformations, the convergent-product
enclosure, and the archimedean and fundamental-subfamily constants.
The [receipt](rank_five_discriminant_count_checks.json) pins the retained
primary sources. Completeness of the local classification and the
infinite-specification counting theorem are explicit imports, not
conclusions from finite computations. No fixed-map target-height
asymptotic, formal verification or external review is asserted.
