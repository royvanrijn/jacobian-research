# The all-`k` transfer-block theorem

> **Status correction.**  The asserted Boolean/ribbon identification is
> false: its norm map misses a second-order class for `k=2` and four classes
> for `k=3`.  Consequently the uniform rank and flatness theorem is not
> proved by this note.  The affine-difference argument remains valid under
> its stated characteristic-zero hypothesis, and the bounded Groebner
> computations remain evidence.  See
> [C22_CONDUCTOR_RIBBON_AUDIT.md](C22_CONDUCTOR_RIBBON_AUDIT.md).
> The replacement direct presentation and current bounded basis results are
> in [DIRECT_TRANSFER_BASIS.md](DIRECT_TRANSFER_BASIS.md).

Let `K` be a characteristic-zero field and let

\[
 S(Z)=Z^k+s_1Z^{k-1}+\cdots+s_k
\]

be the universal monic polynomial.  Write `Z_k` for the formal completion,
along `(U,V)=(S^3,S^2)`, of the scheme of monic polynomials of degrees
`3k,2k` satisfying

\[
                         U^2=V^3.                              \tag{1}
\]

Write `Z_k^aff` for the same completion with (1) weakened to
`U^2-V^3 in K[Z]_(<=1)`.

Equivalently, before completion the strong scheme is the fiber product

\[
 \operatorname{Poly}_{3k}^{\rm mon}
 \mathop{\times}_{\operatorname{Poly}_{6k}^{\rm mon}}
 \operatorname{Poly}_{2k}^{\rm mon},
\]

for the squaring and cubing maps, and `Z_k` is its completion along
`S mapsto (S^3,S^2)`.  Thus `Z_k` is the local square/cube factorization
fiber; see [UNIVERSAL_FACTORIZATION_GEOMETRY.md](UNIVERSAL_FACTORIZATION_GEOMETRY.md)
for the ambient multiplication-map interpretation.

## Former theorem statement (not established)

For every `k>=1`:

1. `Z_k^aff=Z_k` scheme-theoretically.
2. `Z_k -> A^k_S` is finite flat of rank `2^k`.
3. At `S=Z^k`, if `m` is the maximal ideal of the transverse fiber `A_k`,
   then

   \[
   \dim_K m^d/m^{d+1}={k\choose d},\qquad 0\le d\le k.          \tag{2}
   \]

   Thus

   \[
                    \operatorname{Hilb}_{A_k}(t)=(1+t)^k.       \tag{3}
   \]

The fibers need not be Gorenstein.  Formula (3) is a Hilbert-function
statement, not an assertion that the collided algebra is a tensor product of
dual numbers.

The formerly claimed independent proof is given in
[C22_DEFORMATION_AUDIT.md](C22_DEFORMATION_AUDIT.md). It constructs the same
block as the divided-power symmetric product of the conductor ribbon of the
cusp normalization, proves collision flatness by confluent divided
differences, and compares it with the factorization fiber by a
split-surjection/Nakayama argument.  The counteraudit shows that its claimed
surjection fails already for `k=2`; this paragraph records the old strategy,
not a valid proof.

## 1. Affine difference vanishes

Put `D=U^2-V^3` and

\[
                       K_0=2VU'-3UV'.                           \tag{4}
\]

There are polynomial identities

\[
 VD'-3V'D=UK_0,
 \qquad
 UD'-2U'D=V^2K_0.                                \tag{5}
\]

Suppose `D=aZ+b`.  The left side of the first identity in (5) has degree at
most `2k`, whereas `U` is monic of degree `3k`.  Multiplication by a monic
polynomial preserves degree over every coefficient ring, including rings
with nilpotents.  Hence `K_0=0` and

\[
                         Va-3V'(aZ+b)=0.                         \tag{6}
\]

The coefficient of `Z^(2k)` in (6) is `(1-6k)a`, so `a=0`.  The coefficient
of `Z^(2k-1)` is then `-6kb`, so `b=0`.  This proves

\[
                         Z_k^{\rm aff}=Z_k                       \tag{7}
\]

over every `Q`-algebra, without a bounded Groebner calculation.

## 2. Ordered roots and the Boolean thickening

Let

\[
 R=K[r_1,\ldots,r_k],\qquad A=R^{S_k}=K[s_1,\ldots,s_k],
\]

where `S=prod_i(Z-r_i)`.  Introduce commuting variables `epsilon_i` with
`epsilon_i^2=0`, and put

\[
 C=R[\epsilon_1,\ldots,\epsilon_k]/(\epsilon_1^2,\ldots,
                                     \epsilon_k^2).              \tag{8}
\]

The symmetric group permutes the pairs `(r_i,epsilon_i)`.  Define the
Boolean thickening

\[
                         B_k=\operatorname{Spec}(C^{S_k}).       \tag{9}
\]

On ordered roots set

\[
\begin{aligned}
 V&=\prod_{i=1}^k\bigl((Z-r_i)^2+\epsilon_i\bigr),\\
 U&=\prod_{i=1}^k\bigl((Z-r_i)^3+	frac32\epsilon_i(Z-r_i)\bigr).
\end{aligned}                                                    \tag{10}
\]

For `q_i=Z-r_i`,

\[
 (q_i^3+	frac32\epsilon_iq_i)^2
 =(q_i^2+\epsilon_i)^3
\]

because `epsilon_i^2=0`.  Therefore (10) is `S_k`-invariant and defines a
morphism

\[
                         B_k\longrightarrow Z_k.                 \tag{11}
\]

The following elementary factorization lemma identifies it.

### Cusp-factorization lemma

The morphism (11) induces an isomorphism of formal schemes along the common
reduced `S`-space.

**Proof.**  Write uniquely

\[
                         V=S^2+R_0,qquad \deg R_0<k,             \tag{12}
\]

after using the highest `k` coefficients of `V` as the coefficients of `S`.
In the Laurent-series ring at `Z=infinity`, the monic square root is unique:

\[
 V^{3/2}=S^3\left(1+{R_0\over S^2}\right)^{3/2}.                \tag{13}
\]

Equation (1) is equivalent to the vanishing of the negative Laurent part of
(13), with `U` equal to its polynomial part.  Filter its coefficient ring by
the order of `R_0`.  In filtration degree `d`, polarization of the `d`-th
binomial term gives the orbit sums

\[
 \sum_{i_1,\ldots,i_d\ {m distinct}}
 r_{i_1}^{a_1}\cdots r_{i_d}^{a_d}
 \epsilon_{i_1}\cdots\epsilon_{i_d}.             \tag{14}
\]

These are exactly the degree-`d` invariants of (8).  Indeed, every invariant
monomial is such an orbit sum after the powers at indices not carrying an
`epsilon` are removed with symmetric functions in the `r_i`; conversely,
(14) is the product of the polarized power sums

\[
                         \theta_a=\sum_i r_i^a\epsilon_i,         \tag{15}
\]

followed by Mobius inversion on set partitions.  The diagonal terms vanish
because `epsilon_i^2=0`.

Thus (10) identifies the associated graded coefficient algebra of (13) with
`C^(S_k)`, degree by degree.  Both rings are complete and separated for their
transfer ideals.  The filtered morphism (11), being an isomorphism on
associated gradeds, is an isomorphism.  This proves the lemma.  Notice that
the proof uses the entire negative Laurent part; retaining only its first
quadratic equations would lose the higher collision relations. `square`

Consequently

\[
                         Z_k\simeq\widehat B_k                   \tag{16}
\]

along the reduced `S`-space.  Formula (10) is the promised structural model:
the `2/3` transfer block is the symmetric descent of one square-zero cusp jet
at each ordered root.

## 3. Finite flatness and rank

The reflection-group theorem makes `R` a free `A=R^(S_k)`-module of rank
`k!`.  Hence `C` is finite free over `A`.  Since `char K=0`, the Reynolds
operator makes `C^(S_k)` an `A`-direct summand of `C`; it is therefore finite
projective over `A`.

Over the discriminant complement, the ordered-root cover is an `S_k`-torsor.
After that faithfully flat base change, (9) is simply

\[
 K[\epsilon_1,\ldots,\epsilon_k]/(\epsilon_1^2,\ldots,
                                   \epsilon_k^2),                 \tag{17}
\]

which has basis `epsilon_I=prod_(i in I)epsilon_i`, indexed by the `2^k`
subsets of `{1,...,k}`.  Thus `C^(S_k)` has constant rank `2^k`.  Combining
this with (16) proves finite flatness and the rank formula.

## 4. The collided Hilbert series

At `S=Z^k`, base change in (9) gives

\[
 A_k\simeq
 \left(
 {K[r_1,\ldots,r_k]\over(e_1(r),\ldots,e_k(r))}
 \otimes
 {K[\epsilon_1,\ldots,\epsilon_k]\over(\epsilon_1^2,\ldots,
                                               \epsilon_k^2)}
 \right)^{S_k}.                                    \tag{18}
\]

Invariants commute with this base change because the elementary symmetric
functions form an invariant regular sequence and the Reynolds functor is
exact.  The first factor in (18) is the symmetric-group coinvariant algebra,
which is the regular representation of `S_k`.

Grade (18) by epsilon-degree.  The degree-`d` squarefree epsilon space is the
permutation representation on `d`-subsets and has dimension `binom(k,d)`.
For every finite-dimensional representation `M`,

\[
                    \dim(\operatorname{Reg}\otimes M)^{S_k}
                    =\dim M.                                    \tag{19}
\]

Therefore the degree-`d` part of (18) has dimension `binom(k,d)`.

Finally, the polarized power sums (15) generate the invariant algebra over
the symmetric functions: products of the `theta_a` give the distinct-index
orbit sums (14), and Mobius inversion separates every orbit type.  After the
specialization, all positive epsilon-degree elements lie in the maximal ideal
and the algebra is generated in epsilon-degree one.  Epsilon-degree is thus
the maximal-ideal filtration.  Equations (2)--(3) follow.

## 5. Interpretation

Over separated roots, (17) is the ordinary subset algebra: each root carries
one independent square-zero transfer choice.  Collision does not preserve
its multiplication because taking `S_k`-invariants and then specializing
introduces the coinvariant algebra in (18).  It does preserve the free module
rank and the epsilon filtration.  This explains simultaneously:

* the rank `2^k`;
* the Hilbert series `(1+t)^k`;
* the increasingly non-Gorenstein collided fibers; and
* why the Boolean allocation count survives collision without the algebra
  remaining a tensor product of dual numbers.

## Executable regressions

Run

```bash
python scripts/verify_all_k_transfer_block.py
python scripts/verify_c22_deformation_audit.py
```

The script checks the Wronskian identities, the universal one-root cusp jet,
and the previously uncomputed coincident blocks `k=5,6`.  It obtains lengths
`32,64` and Hilbert functions `(1,5,10,10,5,1)` and
`(1,6,15,20,15,6,1)`.  These computations audit the structural theorem; they
are not used to extend a bounded Groebner pattern to arbitrary `k`.
The second script is dependency-free and separately checks the conductor
norms, exact divided differences, triangular norm generation, normalized
compound determinants, and binomial filtration ranks.
