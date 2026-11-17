# Field-discriminant classification of rank-five Hasse failures

The local tests for a rank-five finite étale algebra can be replaced by one criterion on the **maximal-order field discriminant** of its cubic component. This classifies which actual rank-five algebras fail the Hasse principle. It does not classify which integers occur as cubic-field discriminants, and it does not assert that the one fixed map realizes every such algebra.

Write `D_L=Disc(O_L)` for the signed field discriminant of a cubic number field. If `L` is non-Galois, let `K` be its quadratic resolvent and `d_K=Disc(O_K)` its fundamental quadratic discriminant.

## 1. Theorem RDF1: complete rank-five criterion

Let `A` be a finite étale `Q`-algebra of rank five. Its spectrum has points over every completion and no rational point if and only if

\[
 A\simeq K\times L,
\]

where `L` is a non-Galois cubic field, `K` is its quadratic resolvent, and, on writing

\[
 D_L=d_K f^2,\qquad f\in\mathbb Z_{>0},                 \tag{1}
\]

the following conditions hold:

\[
 \boxed{p\mid f,\ p\ne3\ \Longrightarrow\ p\equiv1\pmod3},
 \qquad
 \boxed{3\mid f\ \Longrightarrow\ K\otimes\mathbb Q_3
                                      \simeq\mathbb Q_3^2}.      \tag{2}
\]

The splitting condition at three means `3` does not divide `d_K` and `d_K=1 mod3`. Equivalently, among actual cubic fields, the complete admissible form is

\[
 \boxed{f=m\quad\text{or}\quad f=9m},                  \tag{3}
\]

where `m` is squarefree and every prime dividing `m` is `1 mod3`; in the `9m` case the resolvent must split at three. The value `m=1` is allowed.

In particular Hasse failure for the algebra `K x L` depends only on `D_L`, once `L` is known to be a non-Galois cubic field and `K` is its resolvent. Several nonisomorphic cubic fields may have the same discriminant; no uniqueness of the algebra is asserted.

**Proof.** [FMA2 of the ambient classification](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md) forces the global product `K x L` and its matching quadratic resolvent. Conversely this product has no rational point. Lemma RDF2 proves (1), and Lemma RDF3 identifies the primes dividing `f` exactly with primes where the cubic local algebra is a totally ramified field.

At a prime where the cubic algebra is not a field, it has a `Q_p` factor and therefore a point. If it is an unramified cubic field, it is cyclic and the quadratic resolvent splits, supplying the point. These cases exhaust the primes not dividing `f`.

At a totally ramified cubic prime the cubic has no local root. For `p!=3`, Lemma RDF3 identifies the quadratic resolvent with `Q_p(sqrt(-3))`; this splits exactly when `p=1 mod3`, including the exclusion of `p=2`. At `p=3`, local solubility is exactly splitting of `K`, with no point available from the totally ramified cubic. This proves both directions of (2). The real completion always has a cubic root.

For the sharper form (3), Lemma RDF3 gives `v_p(f)=1` at every prime `p!=3` dividing `f`. At three, a totally ramified cubic has field-discriminant exponent `3,4`, or `5`. If its resolvent splits, its discriminant is a square locally, so the exponent is even, necessarily four; since `v_3(d_K)=0`, this gives `v_3(f)=2`. If it is not totally ramified, `v_3(f)=0`. These are exactly (3), and the converse is already (2). ∎

Theorem RDF1 is a criterion for **given fields**, not a characterization of all integers `d_K m^2` or `81d_K m^2` that are realizable as cubic-field discriminants. Applying it to an arbitrary integer without a supplied cubic field would change the quantifiers.

## 2. Lemma RDF2: the resolvent quadratic order

For every non-Galois cubic field `L`, the quantity `D_L/d_K` is a positive integer square.

**Proof.** Choose an integral basis `beta_1,beta_2,beta_3` of `O_L` and its three embeddings into a normal closure. In the determinant of the embedding matrix, let `a` be the sum of the three positive-permutation products and `b` the sum of the three negative-permutation products. Then

\[
 \det(\sigma_i(\beta_j))=a-b,\qquad (a-b)^2=D_L.
\]

Every summand is an algebraic integer. The Galois group is `S_3`; even permutations fix `a,b` and odd permutations exchange them. Thus `a,b` lie in the quadratic resolvent, are conjugate there, and `a+b` and `ab` are rational algebraic integers, hence integers. Since `a!=b`, the element `a` generates `K` and has monic quadratic polynomial

\[
 T^2-(a+b)T+ab\in\mathbb Z[T].
\]

Its order `Z[a]` has discriminant `(a-b)^2=D_L`. The change-of-integral-basis determinant gives

\[
 D_L=d_K[\mathcal O_K:\mathbb Z[a]]^2,
\]

which proves (1) with positive integer `f`. In a quadratic field the index is also the conductor of this order. This last terminology is not needed for the criterion. ∎

Using a basis of a proper cubic order in this argument would multiply the resulting `f` by its index. The maximal-order hypothesis is therefore indispensable.

## 3. Lemma RDF3: local support and exponents

A prime divides `f` if and only if `L tensor Q_p` is a totally ramified cubic field. For `p!=3` in this case,

\[
 v_p(D_L)=2,\qquad v_p(d_K)=0,\qquad v_p(f)=1,
 \qquad K\otimes\mathbb Q_p=\mathbb Q_p(\sqrt{-3}),     \tag{4}
\]

where the last expression denotes the split quadratic algebra if `-3` is a square. For a totally ramified cubic at three,

\[
 v_3(D_L)\in\{3,4,5\},\qquad v_3(d_K)\in\{0,1\},
 \qquad v_3(f)>0.                                      \tag{5}
\]

**Proof.** The possible cubic étale local algebras are a cubic field, `Q_p` times a quadratic field, and the split product. A cubic field is either unramified or totally ramified, because its degree is prime. The integral algebra is the product of the maximal integral local algebras; its trace form is block diagonal. Hence in the product `Q_p x E_2` its discriminant exponent is exactly that of `E_2`. The quadratic resolvent there is `E_2`, since the sign action exchanges precisely its two embeddings. This proves equality of the exponents in `D_L` and `d_K` in the product case, including wild quadratic ramification at two. In the split case both exponents are zero. An unramified cubic field has zero discriminant exponent and cyclic group of order three, so its quadratic resolvent is split and again both exponents are zero.

Now suppose the local cubic field is totally ramified. A uniformizer `theta` has an Eisenstein minimal polynomial

\[
 g(T)=T^3+aT^2+bT+c,
 \qquad v_p(a),v_p(b)\ge1,\quad v_p(c)=1.
\]

Its powers `1,theta,theta^2` are an integral basis of the maximal local ring. These exact local inputs are Milne, *Algebraic Number Theory*, version 3.00, Proposition 7.55 and Remark 7.56, printed pp.122–123 (PDF pp.130–131): [author-hosted primary text](https://www.jmilne.org/math/CourseNotes/ANT300.pdf), [retained copy](sources/milne-ant300.pdf). The stated uniformizer and maximal-ring arguments were inspected; they are classical imported local-field results.

The cubic discriminant is

\[
 \Delta=a^2b^2-4b^3-4a^3c-27c^2+18abc.                 \tag{6}
\]

For an odd prime `p!=3`, the term `-27c^2` uniquely has valuation two; every other term has larger valuation. Thus `v_p(Delta)=2`, and its squareclass is `-3`, since its ratio to `-27c^2` is one modulo `p` and hence a local square. At two write `a=2A,b=2B,c=2C`, with `C` odd. Exact division gives

\[
 \Delta/4=4A^2B^2-8B^3-16A^3C-27C^2+36ABC
                  \equiv5\pmod8.                    \tag{7}
\]

Indeed the extra term modulo eight is `4AB(AB+C)`, divisible by eight because `C` is odd. Consequently `v_2(Delta)=2`, with squareclass `5=-3` modulo eight. This gives the nonsplit unramified quadratic resolvent at two.

For all `p!=3` this proves (4): the field discriminant has exponent two, while the quadratic resolvent is unramified or split, so its discriminant exponent is zero and `v_p(f)=1`. The field `Q_p` contains the roots of `T^2+T+1` exactly when its residue field has a nontrivial cube root of unity, since `p!=3`. Thus the split case is precisely `p=1 mod3`.

At three every term of (6) has valuation at least three. The upper bound and exact possibilities require retaining maximality, not replacing a field discriminant by an arbitrary polynomial discriminant. In the normalized local valuation `v(theta)=1`, one has `v(3)=3`, and

\[
 g'(\theta)=3\theta^2+2a\theta+b.
\]

With `v_3(0)=infinity`, the three term valuations are `5`, `3v_3(a)+1`, and `3v_3(b)`. Their finite values lie in different classes modulo three, so the minimum is unique. The identity `Disc(g)=-Norm(g'(theta))` and total ramification of degree three give

\[
 v_3(D_L)=\min\{5,\ 3v_3(a)+1,\ 3v_3(b)\}
                  \in\{3,4,5\}.                     \tag{8}
\]

A quadratic field discriminant at the odd prime three has exponent at most one. Equation (1) then makes `v_3(f)` a positive integer. The exponents found above show that no prime outside the totally ramified set divides `f`, and every prime in that set does, proving all assertions. ∎

## 4. Corollary RDF4: discriminant invariance and the fixed map

Two non-Galois cubic fields with the same signed maximal-order field discriminant give the same answer to the Hasse question for their products with their quadratic resolvents. Their fields and associated algebras need not be isomorphic.

For the one fixed quintic map `Phi`, a rational target is Hasse-failing exactly when its inverse polynomial is squarefree of degree five, factors over `Q` as an irreducible quadratic times an irreducible cubic with matching quadratic resolvent, and that cubic **field** satisfies RDF1. The existing `pi=0` and shorter-fibre exclusions are retained. Thus no separate prime-by-prime cubic root tests are needed once the maximal-order field discriminant is known; computing that maximal-order invariant is an additional algebraic step.

For an integral Hasse-failing target of this fixed map, [ILO3a](FIXED_MAP_INTEGRAL_LOCAL_OBSTRUCTIONS.md) shows that the quadratic resolvent ramifies at three. It is therefore not split over `Q_3`, and RDF1 forces

\[
 3\nmid f,\qquad v_3(D_L)=v_3(d_K)=1.                 \tag{9}
\]

This is necessary, not a construction of an integral Hasse target. The integral-target existence problem remains open.

Conversely every rank-five algebra classified by RDF1 is realized as a full fibre of some Jacobian-one map on `A^3_Q`, by Common Arithmetic Fibers. This varying-map conclusion does not assert that the single fixed `Phi` realizes every such algebra. No map-equivalence or prescribed-generic-monodromy classification is implied.

**Proof.** The fundamental discriminant `d_K` and positive square factor `f` are determined by `D_L`; apply RDF1. The fixed-map statement follows from FMA1–FMA3 and RDF1. Equation (9) combines the integral resolvent theorem with the split-at-three alternative in (2). The varying-map assertion is the rank-five specialization of [FMA7](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md). ∎

## 5. Example RDF5: fundamental discriminants and the maximal-order warning

For `L=Q(alpha)` with `alpha^3-alpha-1=0`, the polynomial discriminant is `-23`. It is squarefree, so `Z[alpha]` is maximal. Its resolvent has fundamental discriminant `-23`, and `f=1`. RDF1 shows that `Q(sqrt(-23)) x L` is Hasse-failing. Thus neither the quadratic resolvent `Q(sqrt(-3))` nor total ramification at any prime is necessary.

The distinction between polynomial and field discriminants is already decisive in the classical pure cubic example. Let `alpha^3=19` and put

\[
 \beta=(1+\alpha+\alpha^2)/3.
\]

The relations

\[
 \alpha^2=3\beta-\alpha-1,\quad
 \alpha\beta=\beta+6,\quad
 \beta^2=\beta+2\alpha+4,
\]

and the monic equation `beta^3-beta^2-6beta-12=0` show integrality and closure of the lattice with basis `1,alpha,beta`. Its trace Gram matrix and discriminant are

\[
 \begin{pmatrix}3&0&1\\0&0&19\\1&19&13\end{pmatrix},
 \qquad D=-3\cdot19^2.                                \tag{10}
\]

This order is maximal at nineteen because `alpha` is an Eisenstein uniformizer there and its powers generate the maximal local ring. At three the discriminant exponent one precludes a nontrivial index, whose square would divide the discriminant. All other primes have unit discriminant. Hence it is the full ring of integers and

\[
 D_L=-3\cdot19^2,\qquad d_K=-3,\qquad f=19.
\]

The criterion passes, since `19=1 mod3`. Using instead the polynomial discriminant `Disc(T^3-19)=-27*19^2` would produce the false value `f_order=3*19` and would incorrectly reject the fibre at three. This is precisely the failure avoided by the maximal-order hypothesis.

## 6. Exact controls and scope

Run

```sh
.venv/bin/python research/arithmetic/verify_rank_five_discriminant_classification.py --output research/arithmetic/rank_five_discriminant_classification_checks.json
```

The checker verifies the determinant/quadratic-order algebraic identity, the local cubic discriminant reductions, the complete finite two-adic unit residue step in (7), bounded controls of the unrestricted wild-three valuation formula, and the integral-basis relations and trace discriminant in Example RDF5. These finite controls do not classify all cubic number fields or replace the written local-field and maximal-order proofs. No bounded search is promoted to an unrestricted theorem.

The primary local input is the exact Milne version identified above. The global forced `2+3` factorization and matching-resolvent input retain FMA2's classical Galois/local correspondence scope. Formal verification, independent whole-proof replay, external review, integer-discriminant realizability and a classification of maps are not claimed.
