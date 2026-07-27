# Primitive monodromy makes a Keller map compositionally atomic

Let `k` be a characteristic-zero field.  Write `\mathcal K_d(k)` for the
monoid, under composition, of polynomial maps

\[
 F:\mathbb A_k^d\longrightarrow\mathbb A_k^d
\]

with nonzero constant Jacobian determinant.  Its units are the polynomial
automorphisms.  A nonunit `F` is called **atomic** if every factorization

\[
 F=G\circ H,\qquad G,H\in\mathcal K_d(k),
\]

has a unit factor.  This is the noncommutative-monoid notion of an atom; it
should not be confused with a two-sided prime element.

For a dominant generically finite map `F`, its geometric monodromy is the
Galois group of the Galois closure of

\[
 \bar k(F_1,\ldots,F_d)\subseteq\bar k(x_1,\ldots,x_d),
                                                               \tag{1}
\]

acting on the geometric generic fiber.  Recall that a transitive
permutation action is **primitive** if it has no nontrivial block system.

> **Primitive-monodromy atomicity theorem.**  
> Let `F:\mathbb A_k^d\to\mathbb A_k^d` be a Keller map of geometric degree
> `N>1`.  If its geometric monodromy action is primitive, then:
>
> 1. `F` is atomic, including after every characteristic-zero extension of
>    the constant field;
> 2. for every `r\ge0`, the stabilization
>    \[
>    F^{[r]}=F\times\operatorname{id}_{\mathbb A^r}
>    \]
>    is atomic;
> 3. every polynomial left--right equivalent map
>    `A\circ F^{[r]}\circ B`, with `A,B` automorphisms, is atomic.
>
> In particular, geometric monodromy `S_N` implies absolute and stable
> atomicity.

## 1. Primitive actions and intermediate fields

Work first over an algebraically closed field.  Put

\[
 K=k(F_1,\ldots,F_d),\qquad L=k(x_1,\ldots,x_d),
\]

and let `M/K` be the Galois closure.  If

\[
 G_F=\operatorname{Gal}(M/K),\qquad
 P=\operatorname{Gal}(M/L),
\]

then the generic sheets are the cosets `G_F/P`, and `P` is a point
stabilizer.  Galois correspondence gives

\[
 \left\{K\subseteq E\subseteq L\right\}
 \longleftrightarrow
 \left\{P\subseteq J\subseteq G_F\right\}.               \tag{2}
\]

For a transitive permutation group, the action is primitive if and only if
a point stabilizer is maximal.  Indeed, a subgroup `J` between `P` and
`G_F` gives the block consisting of the `J`-orbit of the distinguished
point; conversely, the setwise stabilizer of a block containing that point
is such an intermediate subgroup.

Consequently, primitivity is equivalent to the field-theoretic statement

\[
 \boxed{\text{there is no }K\subsetneq E\subsetneq L.}   \tag{3}
\]

For the natural action of `S_N`, the stabilizer is `S_{N-1}` and is maximal.
A short proof is useful here.  If
`S_{N-1}\subsetneq J\subseteq S_N`, an element of `J` moves the
distinguished point.  The subgroup `S_{N-1}` moves its image through the
other `N-1` points, so `J` is transitive.  Its stabilizer of the
distinguished point is exactly `S_{N-1}`.  Orbit--stabilizer therefore gives

\[
 |J|=N(N-1)!=N!,
\]

and hence `J=S_N`.

## 2. Factors of a Keller map are Keller

Suppose a Keller map factors as polynomial self-maps

\[
 F=G\circ H.                                             \tag{4}
\]

Both `G` and `H` are dominant.  Dominance of `G` follows because the image
of `F` is dense.  If the image of `H` had dimension less than `d`, its image
under `G` would also have dimension less than `d`, contradicting dominance
of `F`.  Since source and target have the same dimension, both factors are
generically finite.

The chain rule gives

\[
 \det DF=(\det DG)\circ H\cdot\det DH\in k^\times.       \tag{5}
\]

The two factors on the right lie in the polynomial ring
`k[x_1,\ldots,x_d]`, whose only units are the nonzero constants.  Hence
`\det DH\in k^\times` and `(\det DG)\circ H\in k^\times`.  Since `H` is
dominant, pullback

\[
 H^*:k[y_1,\ldots,y_d]\hookrightarrow k[x_1,\ldots,x_d]
\]

is injective.  It follows that `\det DG` itself is constant.  Thus every
polynomial factor of a Keller map is again a Keller map; this did not need
to be assumed in the definition of the factorization.

## 3. A degree-one Keller map is an automorphism

We record the standard birational Keller lemma.

> **Lemma.** If a Keller self-map of `\mathbb A_k^d` has geometric degree
> one, it is a polynomial automorphism.

It is enough to prove this after passing to an algebraic closure, because
being an isomorphism descends under a faithfully flat field extension.  The
map is etale and birational.  An etale morphism is quasi-finite, so the
birational form of Zariski's Main Theorem, with normal target
`\mathbb A^d`, identifies it with an open immersion

\[
 \mathbb A^d\simeq U\hookrightarrow\mathbb A^d.          \tag{6}
\]

This open immersion cannot be proper.  If the complement contained a
divisor, factoriality of affine space would give an irreducible polynomial
vanishing exactly on that divisor.  Its restriction would be a nonconstant
unit on `U\simeq\mathbb A^d`, whereas affine space has only constant units.
Thus the complement has codimension at least two.  Normal Hartogs extension
then gives

\[
 \Gamma(U,\mathcal O_U)=\Gamma(\mathbb A^d,\mathcal O)=
 k[x_1,\ldots,x_d].
\]

Since `U` is affine, the inclusion in (6) is the morphism induced by this
isomorphism of coordinate rings and is therefore an isomorphism.  This
proves the lemma.

## 4. Proof of atomicity

Now suppose the monodromy of `F` is primitive and take a factorization (4).
Pullback of rational functions produces an intermediate field

\[
 K
 \subseteq k(H_1,\ldots,H_d)
 \subseteq L.                                           \tag{7}
\]

By (3), one inclusion in (7) is an equality.  If the right inclusion is an
equality, `H` has geometric degree one.  If the left inclusion is an
equality, `G` has geometric degree one.  Section 2 says the degree-one
factor is Keller, and Section 3 says it is an automorphism.  Thus every
factorization of `F` has a unit factor, proving atomicity.

This also proves the claim over the original, possibly nonclosed, field:
base change a proposed factorization to an algebraic closure, apply the
geometric argument, and descend invertibility of the unit factor.

The assertion is absolute under extension of the constant field.  After
embedding the old and new algebraic closures in a common algebraically
closed overfield, constant extension does not change the Galois closure,
its permutation group, or its point stabilizer.  Hence primitivity and the
preceding proof persist.  Equivalently, the finite generic cover and its
block systems are invariant under extension of algebraically closed
constants.

## 5. Stabilization and left--right equivalence

Adjoin independent variables `T_1,\ldots,T_r`.  The function-field extension
of the stabilization is

\[
 K(T_1,\ldots,T_r)\subseteq L(T_1,\ldots,T_r).           \tag{8}
\]

If `M/K` is the old Galois closure, then

\[
 M(T_1,\ldots,T_r)/K(T_1,\ldots,T_r)
\]

is the new Galois closure, with the same Galois group and the same point
stabilizer.  Thus the stabilized monodromy action is still primitive.
Applying Section 4 in dimension `d+r` proves that every `F^{[r]}` is
atomic.

Finally let `A,B` be polynomial automorphisms.  A factorization

\[
 A\circ F^{[r]}\circ B=G\circ H
\]

would give

\[
 F^{[r]}=(A^{-1}\circ G)\circ(H\circ B^{-1}).
\]

Atomicity of `F^{[r]}` makes one parenthesized factor an automorphism, and
therefore makes `G` or `H` an automorphism.  Atomicity is consequently
preserved by polynomial left--right equivalence and by stable polynomial
left--right equivalence.

## 6. Universal weighted consequence

Let `H(W)` be any admissible weighted seed of degree `N\ge3`.  The
[weighted-seed theorem](WEIGHTED_SEED_THEOREM.md) makes its associated map
`F_H:\mathbb A^3\to\mathbb A^3` Keller of geometric degree `N`, with generic
inverse equation

\[
 H(W)-sW+t=0.
\]

The [universal symmetric-monodromy theorem](UNIVERSAL_SYMMETRIC_MONODROMY.md)
gives geometric monodromy `S_N` for every such `H`, without a genericity
hypothesis.  The theorem above therefore yields

\[
\boxed{\text{Every admissible weighted Keller map of degree }N\ge3
\text{ is absolutely and stably atomic.}}               \tag{9}
\]

There is a particularly sparse explicit seed in every degree:

\[
 H_N(W)=\frac{W^2-W^N}{N-2}.                            \tag{10}
\]

Indeed,

\[
 H_N(0)=H_N'(0)=H_N(1)=0,\qquad H_N'(1)=-1,
\]

\[
 H_N''(1)=-(N+1),\qquad
 a_0=-\frac{N}{N-1}.
\]

Thus (10) is admissible for every `N\ge3`.  Put

\[
 u=1+xy,\qquad
 \gamma=1-\frac{N}{N-1}xy+x^2z.
\]

The associated determinant-one map is

\[
\boxed{
 F_N=
 \left(
 \frac{(N-2)u+u^2-(N-1)u^N\gamma^{N-2}}{(N-2)x^2},
 \frac{(N-2)+2u-Nu^{N-1}\gamma^{N-2}}{(N-2)x},
 x\gamma
 \right).
}                                                       \tag{11}
\]

Weighted polynomiality says that the first two displayed numerators are
divisible by `x^2` and `x`, respectively.  Hence (11) is a polynomial map.
For every `N\ge3`,

\[
 \det DF_N=1,\qquad
 \deg_{\mathrm{geom}}F_N=N,\qquad
 \operatorname{Mon}_{\mathrm{geom}}(F_N)=S_N,
\]

and `F_N` and every stabilization of it are atomic.  Taking `N=12` recovers
the explicit composite-degree map in
[the degree-twelve certificate](INDECOMPOSABLE_COMPOSITE_DEGREE.md).
