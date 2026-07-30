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

There is one necessary limitation.  For \(r\geq3\), the Boolean
intermediate-field lattice does not prove that every maximal
**polynomial-sandwich** chain is field-maximal.  An alternative affine model
over one Boolean field might fail to nest with the coordinate models at
adjacent ranks.  Therefore (4.2)--(4.5) prove:

- a full Jordan--Hölder theorem for the intermediate-field towers;
- a strict coherent Ritt \(2\)-complex for the canonical split polynomial
  towers; and
- that the maximum polynomial decomposition length is \(r\);

but not that every maximal polynomial decomposition has length \(r\).
The two-factor theorem is stronger because its diamond has no rank available
to skip.

## 5. Stable factor classes: the remaining gate

For the two coordinate decompositions in (3.1), the stable factor-class
multiset is exactly

\[
 \{[F_3]_{\rm stLR},[F_4]_{\rm stLR}\}.                \tag{5.1}
\]

The group diamond does not prove that (5.1) holds for every maximal
polynomial decomposition.  A middle function field can, in principle,
admit more than one compatible polynomial affine-space model inside
\(k[X,Y]\).  Such models have the same degree and monodromy data but need
not be stably polynomially left--right equivalent.

Thus the four questions have the following exact answer for \(K_{3,4}\):

| proposed invariant | status |
|---|---|
| length | proved invariant |
| degree multiset | proved invariant |
| monodromy composition factors | proved invariant |
| stable factor classes | open beyond the two displayed coordinate models |

If each of the two middle fields in (2.1) has a unique compatible
affine-space reconstruction open up to stable polynomial left--right
equivalence, then the stable factor-class multiset is invariant as well.
This is precisely the affine-model uniqueness gate, not a further
permutation-group calculation.

## 6. Exact regression

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
the canonical Boolean maximal chains through rank six.
