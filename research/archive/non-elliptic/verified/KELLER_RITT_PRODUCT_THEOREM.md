# A Keller--Ritt theorem for separated products

Work over a characteristic-zero field \(k\).  For \(n\geq3\), let
\(F_n:\mathbb A^3\to\mathbb A^3\) be the sparse determinant-one weighted
map from

\[
 H_n(W)=\frac{W^2-W^n}{n-2}.
\]

The [primitive-monodromy atomicity theorem](PRIMITIVE_MONODROMY_ATOMICITY.md)
gives

\[
 \operatorname{gdeg}F_n=n,\qquad
 \operatorname{Mon}_{\rm geom}(F_n)=S_n,
\]

and makes \(F_n\), every identity stabilization of it, and every stable
left--right equivalent map atomic.

## 1. The separated product

For \(a,b\geq3\), write the coordinates of \(\mathbb A^6\) as
\((X,Y)\in\mathbb A^3\times\mathbb A^3\), and put

\[
 A=F_a\times\operatorname{id}_{\mathbb A^3},\qquad
 B=\operatorname{id}_{\mathbb A^3}\times F_b.
\]

These two stabilized atoms commute.  Hence

\[
\boxed{
 K_{a,b}=F_a\times F_b
 =B\circ A=A\circ B.                                  \tag{1.1}
}
\]

The map \(K_{a,b}\) is Keller with determinant one and geometric degree
\(ab\).  Its two displayed decompositions have ordered degree words
\((a,b)\) and \((b,a)\).

The disjoint-coordinate generic covers have linearly disjoint Galois
closures.  Thus

\[
 \operatorname{Mon}_{\rm geom}(K_{a,b})=S_a\times S_b \tag{1.2}
\]

in the product action on
\(\{1,\ldots,a\}\times\{1,\ldots,b\}\).  A point stabilizer is

\[
 P=S_{a-1}\times S_{b-1}.                              \tag{1.3}
\]

## 2. The exact diamond

> **Separated-product Keller--Ritt theorem.**
> The interval between \(P\) and \(S_a\times S_b\) has exactly four
> elements:
>
> \[
> \begin{array}{ccccc}
> &&S_a\times S_b&&\\
> &\nearrow&&\nwarrow&\\
> S_a\times S_{b-1}&&&&S_{a-1}\times S_b\\
> &\nwarrow&&\nearrow&\\
> &&S_{a-1}\times S_{b-1}.&&
> \end{array}                                          \tag{2.1}
> \]
>
> Consequently every maximal polynomial decomposition of \(K_{a,b}\):
>
> 1. has length two;
> 2. has degree multiset \(\{a,b\}\); and
> 3. has factor-monodromy group multiset \(\{S_a,S_b\}\).
>
> The two decompositions in (1.1) are an explicit Keller Ritt move.  In
> particular, the ordered degree word is not invariant.

### Proof of the group interval

Put \(G_1=S_a\), \(G_2=S_b\), \(P_1=S_{a-1}\), and
\(P_2=S_{b-1}\).  Let

\[
 P_1\times P_2\leq J\leq G_1\times G_2.                \tag{2.2}
\]

Each projection of \(J\) contains \(P_i\).  Since \(P_i\) is maximal in
\(G_i\), the projection is either \(P_i\) or \(G_i\).

If the projections are \(P_1,P_2\), then \(J=P_1\times P_2\).  If they
are \(G_1,P_2\), multiplication by the contained subgroup
\(1\times P_2\) shows that \(G_1\times1\subseteq J\), so
\(J=G_1\times P_2\).  The other mixed case is symmetric.

Suppose both projections are surjective.  The kernel of the second
projection is a normal subgroup of \(G_1\) containing \(P_1\).  The normal
closure of a point stabilizer \(S_{a-1}\) in \(S_a\) is all of \(S_a\);
therefore \(G_1\times1\subseteq J\).  Symmetrically
\(1\times G_2\subseteq J\), and \(J=G_1\times G_2\).  This proves (2.1).

### Passage to polynomial decompositions

Every polynomial decomposition gives an intermediate function field and
hence an element of (2.1).  Strictly comparable polynomial-sandwich rings
have strictly comparable fraction fields: if two adjacent rings had the
same fraction field, the intervening Keller factor would have geometric
degree one and hence would be a polynomial automorphism, forcing the two
rings to be equal.

The two middle fields in (2.1) are both realized by the displayed
coordinatewise polynomial sandwiches in (1.1).  Therefore an endpoint-only
chain is not maximal, while the diamond leaves no room for a chain of more
than two nonunit factors.  Every maximal polynomial decomposition has
length two.

The two middle fields give the index pairs \((a,b)\) and \((b,a)\), proving
the degree-multiset assertion.  Their two field extensions have monodromy
\(S_a\) and \(S_b\), independently of which compatible affine polynomial
model realizes the middle field.  This proves the final assertion.
\(\square\)

## 3. The explicit degree-twelve Ritt diamond

The smallest unequal pair in the explicit atomic family is

\[
\boxed{
 K_{3,4}=F_3\times F_4
 =
 (\operatorname{id}\times F_4)\circ(F_3\times\operatorname{id})
 =
 (F_3\times\operatorname{id})\circ
 (\operatorname{id}\times F_4).
}                                                       \tag{3.1}
\]

Thus one determinant-one map of \(\mathbb A^6\), of geometric degree
twelve, has the two maximal words

\[
 (3,4),\qquad(4,3).                                    \tag{3.2}
\]

All maximal decompositions of this map have:

| datum | invariant value |
|---|---|
| length | \(2\) |
| degree multiset | \(\{3,4\}\) |
| factor monodromy groups | \(\{S_3,S_4\}\) |
| factor-monodromy composition factors | \(C_2^4,C_3^2\) |

This is a positive Keller--Ritt theorem, not merely a bounded search.
It also supplies the first literal Keller Ritt square in the explicit
weighted family.  It does not contradict the full-wreath order rigidity of
\(F_4\circ F_3\): the product map has direct-product monodromy, whereas the
three-dimensional composite has \(S_3\wr S_4\).

## 4. Boolean towers and strict Coxeter coherence

The same construction has a useful all-arity form.  For
\(\mathbf n=(n_1,\ldots,n_r)\), with every \(n_i\geq3\), put

\[
 K_{\mathbf n}=\prod_{i=1}^r F_{n_i}:
 \mathbb A^{3r}\longrightarrow\mathbb A^{3r}.          \tag{4.1}
\]

Its geometric monodromy and point stabilizer are

\[
 G=\prod_{i=1}^r S_{n_i},\qquad
 P=\prod_{i=1}^r S_{n_i-1}.                            \tag{4.2}
\]

> **Boolean-tower proposition.**  The subgroup interval \([P,G]\) is the
> Boolean lattice \(2^{\{1,\ldots,r\}}\).  Its maximal field chains are the
> \(r!\) labeled orders in which the \(r\) coordinate factors are introduced
> (some degree words coincide when degrees repeat).
> Every such chain has length \(r\), degree multiset
> \(\{n_1,\ldots,n_r\}\), and factor-monodromy group multiset
> \(\{S_{n_1},\ldots,S_{n_r}\}\).

Indeed, if \(P\leq J\leq G\), the \(i\)-th projection of \(J\) is either
\(S_{n_i-1}\) or \(S_{n_i}\).  In the second case,
\(J\cap S_{n_i}\) is normal in the full \(i\)-th projection and contains
\(S_{n_i-1}\).  The normal closure of that point stabilizer is \(S_{n_i}\),
so \(J\) contains the whole \(i\)-th factor.  Thus \(J\) is exactly the
product obtained by choosing, independently at each coordinate, either
\(S_{n_i-1}\) or \(S_{n_i}\).

Every Boolean node has a coordinatewise polynomial sandwich.  The canonical
split composition series are therefore indexed by permutations of the
factors.  Adjacent swaps are literal commuting Keller Ritt moves.  They
satisfy the Coxeter relations strictly:

\[
 s_i s_j=s_j s_i\quad(|i-j|>1),\qquad
 s_i s_{i+1}s_i=s_{i+1}s_is_{i+1}.                    \tag{4.3}
\]

For \(F_3\times F_4\times F_5\), the six canonical words form a filled braid
hexagon.  The two half-braids

\[
\begin{aligned}
345&\to435\to453\to543,\\
345&\to354\to534\to543
\end{aligned}                                         \tag{4.4}
\]

are equal as polynomial-map factorizations of the same separated product,
not merely after reduction or normalization.  With four factors, the
canonical Coxeter \(2\)-skeleton has \(24\) vertices, \(36\) edges, six
commuting squares, and eight braid hexagons.  The stable factor-class
multiset on every canonical split series is

\[
 \{[F_{n_1}]_{\rm stLR},\ldots,[F_{n_r}]_{\rm stLR}\}. \tag{4.5}
\]

The next lemma removes the apparent limitation: once one polynomial
Keller model of an embedded intermediate field is known, its compatible
affine reconstruction open is unique.

## 5. Affine-open rigidity

Let

\[
 X\xrightarrow{H_0}U_0\xrightarrow{G_0}Y,\qquad
 X\simeq U_0\simeq Y\simeq\mathbb A^d                 \tag{5.1}
\]

be a factorization of a Keller map, and put

\[
 k(Y)\subset E=k(U_0)\subset k(X).
\]

Let \(Z\) be the normalization of \(Y\) in the embedded field \(E\).
Zariski's Main Theorem identifies \(U_0\) with an open subset of \(Z\).

> **Affine-open rigidity lemma.**
> If a second polynomial Keller factorization
> \[
> X\xrightarrow{H}U\xrightarrow{G}Y,\qquad U\simeq\mathbb A^d,
> \tag{5.2}
> \]
> realizes the same embedded intermediate field \(E\subset k(X)\), then
> its open immersion \(U\hookrightarrow Z\) has the same image as
> \(U_0\hookrightarrow Z\).  Consequently there is a polynomial
> automorphism \(\psi:U\to U_0\) such that
> \[
> \boxed{G=G_0\circ\psi,\qquad H=\psi^{-1}\circ H_0.}  \tag{5.3}
> \]

### Proof

We may extend the constants to an algebraic closure.  A Keller map is
etale, hence open.  Its image cannot omit a divisor: if the irreducible
hypersurface \(V(q)\subset\mathbb A^d\) were disjoint from the image, then
\(q\circ H_0\) would be a nowhere-vanishing polynomial on affine space and
hence a constant unit.  Dominance of \(H_0\) makes pullback injective, so
this would force \(q\) to be constant.  Therefore

\[
 \operatorname{codim}_{U_0}
 \bigl(U_0\setminus H_0(X)\bigr)\geq2.                \tag{5.4}
\]

The inclusion \(E\subset k(X)\) gives a canonical lift

\[
 \ell_E:X\longrightarrow Z.                           \tag{5.5}
\]

Indeed every element of the normalization ring is integral over \(k[Y]\);
after pullback it is integral over the normal ring \(k[X]\), and hence lies
in \(k[X]\).  Both factorizations in (5.1)--(5.2) induce this same lift.
Thus

\[
 \ell_E(X)\subset U_0\cap U.                           \tag{5.6}
\]

By (5.4), every codimension-one point of \(U_0\) belongs to \(U\).

Write the codimension-one boundary of \(U_0\) in \(Z\) as
\[
 D_1\cup\cdots\cup D_s.
\]
There are no additional higher-codimension boundary components.  This is
the standard affine Hartogs argument: after localizing away from the
divisorial components, a hypothetical codimension-at-least-two component
would give a proper affine open in a normal affine scheme with the same
ring of global functions.

The localization sequence for divisor class groups gives

\[
 k[U_0]^\times/k[Z]^\times
 \longrightarrow
 \bigoplus_{i=1}^s\mathbb ZD_i
 \longrightarrow
 \operatorname{Cl}(Z)
 \longrightarrow
 \operatorname{Cl}(U_0)
 \longrightarrow0.                                   \tag{5.7}
\]

Here \(k[U_0]^\times=k^\times\), every unit of \(Z\) restricts to a
constant, and \(\operatorname{Cl}(U_0)=0\).  Hence

\[
 \boxed{\operatorname{Cl}(Z)
 \simeq\bigoplus_{i=1}^s\mathbb Z[D_i].}              \tag{5.8}
\]

Let \(S\) be the set of boundary primes \(D_i\) whose generic points lie in
the competing affine open \(U\).  Since \(U\) already contains every
codimension-one point of \(U_0\), another application of the localization
sequence yields

\[
 \operatorname{Cl}(U)
 \simeq
 \operatorname{Cl}(Z)/
 \langle[D_i]:D_i\notin S\rangle
 \simeq\bigoplus_{D_i\in S}\mathbb Z[D_i].            \tag{5.9}
\]

But \(U\simeq\mathbb A^d\) is factorial, so \(S\) is empty.  Thus
\(U\subset U_0\), and (5.4)--(5.6) show that \(U_0\setminus U\) has
codimension at least two.  Normal Hartogs extension gives
\(k[U]=k[U_0]\).  Since both opens are affine, the inclusion is an
isomorphism, so \(U=U_0\) inside \(Z\).

The two chosen affine-space identifications now differ by a polynomial
automorphism \(\psi\).  Equality of the canonical lifts (5.5), followed by
restriction of \(Z\to Y\), gives (5.3).  The argument descends from the
algebraic closure because equality of the two open subschemes and the
resulting coordinate-ring automorphism are faithfully flat statements.
\(\square\)

The external structural inputs are the algebraic
[Zariski Main Theorem](https://stacks.math.columbia.edu/tag/03GS),
[purity of the complement of a dense affine open](https://stacks.math.columbia.edu/tag/0BCQ),
and the
[Weil-divisor class presentation](https://stacks.math.columbia.edu/tag/0BE0).

This lemma is stronger than a Torelli comparison of finite marked fibres:
it uses the full finite normalization, the common source lift, factoriality,
and the distinguished affine opens.  It proves that the affine-open gate is
automatic only when the embedded field already has one polynomial Keller
model.  It does not turn an arbitrary imprimitive intermediate field into a
polynomial sandwich.

### The two \(K_{3,4}\) normalizations

Write \(\overline X_n\) for the normalization of the target of \(F_n\) in
the source function field of \(F_n\).  The two middle fields of \(K_{3,4}\)
are

\[
\begin{aligned}
 E_{3|4}&=k(F_3(X),Y),\\
 E_{4|3}&=k(X,F_4(Y)).
\end{aligned}                                         \tag{5.10}
\]

Normalization commutes with the smooth affine-space base change, so their
finite normalizations over the target of \(K_{3,4}\) are intrinsically

\[
\boxed{
\begin{aligned}
 Z_{3|4}&\simeq\mathbb A^3\times\overline X_4,\\
 Z_{4|3}&\simeq\overline X_3\times\mathbb A^3.
\end{aligned}}                                        \tag{5.11}
\]

Their canonical reconstruction opens are the displayed products of the
affine targets with the distinguished opens in \(\overline X_4\) and
\(\overline X_3\).  Their ordered boundary primes, ramification labels,
relative Fitting ideals, conductors, and affine-versus-boundary selectors
are exactly the factor packages pulled back along the affine projection.
This follows either directly from (5.11) or from smooth-base-change
functoriality in the
[stable-normalization theorem](STABLE_NORMALIZATION_FUNCTORIALITY.md).

The affine-open rigidity lemma proves that these are the only compatible
affine-space opens.  The full decorated packages therefore restrict to the
canonical stabilized factors, and the
[decorated Torelli proposition](POLYNOMIAL_GAUGE_DECORATED_TORELLI.md#1-the-full-decorated-normalization-is-complete)
turns their package identifications into ordinary polynomial left--right
equivalences.  Thus the proposed normalization, boundary-selector, and
Torelli attack closes without a degree-specific elimination.

## 6. Stable Keller--Ritt theorem

Return to the separated product \(K_{\mathbf n}\).  Every Boolean
intermediate field has its coordinatewise polynomial Keller model.  The
affine-open rigidity lemma says that every other polynomial model of that
same embedded field differs only by a polynomial automorphism of the
intermediate affine space.

Suppose a polynomial decomposition chain skipped a rank of the Boolean
field lattice.  Insert the canonical coordinate model at that rank and
transport the adjacent coordinate factors through the automorphisms (5.3).
This strictly refines the polynomial decomposition, contradicting
maximality.  Hence every maximal polynomial decomposition is field-maximal.

> **Stable Keller--Ritt theorem for separated products.**
> For every \(r\geq2\) and every \(n_i\geq3\), every maximal polynomial
> decomposition of
> \[
> K_{\mathbf n}=F_{n_1}\times\cdots\times F_{n_r}
> \]
> has length \(r\), and its stable factor-class multiset is
> \[
> \boxed{
> \{[F_{n_1}]_{\rm stLR},\ldots,[F_{n_r}]_{\rm stLR}\}.} \tag{6.1}
> \]
> After absorbing intermediate polynomial automorphisms into adjacent
> factors, every maximal decomposition is one of the \(r!\) canonical
> split orders.  Any two are connected by adjacent commuting Keller Ritt
> moves.  On the canonical representatives the commuting squares and braid
> relations hold literally.

For \(K_{3,4}\), this closes the former affine-model gap.  Every maximal
decomposition has the complete invariant table

| datum | invariant value |
|---|---|
| length | \(2\) |
| degree multiset | \(\{3,4\}\) |
| factor monodromy groups | \(\{S_3,S_4\}\) |
| stable factor classes | \(\{[F_3]_{\rm stLR},[F_4]_{\rm stLR}\}\) |
| factor-monodromy composition factors | \(C_2^4,C_3^2\) |

## 7. Exact regression

Run

```bash
.venv/bin/python scripts/verify_keller_ritt_product.py
```

The checker constructs \(F_3,F_4\), verifies their determinants and the two
commuting compositions, enumerates all \(2^{12}\) candidate blocks
containing one sheet in the \(S_3\times S_4\) action, and confirms the exact
four-element diamond and its induced local and quotient actions.  It then
constructs the strict three-factor braid hexagon and a four-factor commuting
square, checks the \(24\)-vertex Coxeter \(2\)-skeleton counts, and enumerates
the canonical Boolean maximal chains through rank six.  The affine-open
rigidity lemma is the abstract divisor-class argument in Section 5; the
finite checker does not substitute for that proof.
