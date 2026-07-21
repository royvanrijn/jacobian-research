# Monodromy, deck transformations, and rigid reconstruction opens

This note isolates the general mechanism behind target-fixed cancellation construction rigidity.
The mechanism is broader than cancellation construction, but its conclusion must be stated
carefully: trivial deck group makes a target-fixed identification unique; it
does not prove that two unrelated finite extensions are isomorphic.

## 1. Automorphisms of a finite normalization

Let `Y` be a normal integral variety with function field `K`, let `L/K` be a
finite separable field extension, and let

\[
 \pi:\overline X=\operatorname{Norm}_Y(L)\longrightarrow Y              \tag{1}
\]

be finite.  Choose a Galois closure `M/K` and put

\[
 G=\operatorname{Gal}(M/K),\qquad H=\operatorname{Gal}(M/L).
\]

The generic fiber is the transitive arithmetic-monodromy `G`-set `G/H`, after
the usual choice of an embedding of `L` into `M`.  After extending the
constant field to an algebraic closure, the corresponding normal subgroup is
the geometric monodromy group.

### Theorem 1.1 (deck-centralizer theorem)

There are natural identifications

\[
 \operatorname{Aut}_Y(\overline X)
 \simeq \operatorname{Aut}_K(L)
 \simeq N_G(H)/H.                                           \tag{2}
\]

Up to passage to the opposite group according to left/right coset
conventions, this is also the centralizer of the monodromy action:

\[
 \operatorname{Aut}_G(G/H)\simeq C_{\operatorname{Sym}(G/H)}(G).         \tag{3}
\]

Consequently the following are equivalent:

1. the finite normalization has no nonidentity target-fixed automorphism;
2. the generic cover has no nonidentity deck transformation;
3. the monodromy centralizer is trivial; and
4. the sheet stabilizer `H` is self-normalizing in `G`.

**Proof.**  A `Y`-automorphism of the normal integral scheme `bar(X)` induces
a `K`-automorphism of `L`.  Conversely, a `K`-automorphism of `L` preserves
the elements integral over every affine open of `Y`, and therefore extends
uniquely to the normalization.  This proves the first identification.

An element of `G` carries the fixed field `M^H=L` to itself exactly when it
normalizes `H`; two such elements induce the same automorphism of `L` exactly
when they differ by an element of `H`.  This proves the second identification
in (2).  Finally, every endomorphism of the transitive `G`-set `G/H` is
determined by the image of `H`, and that image must be fixed by `H`.  These
fixed cosets are precisely those represented by `N_G(H)`, giving (3).  QED

### Corollary 1.2 (natural symmetric and alternating monodromy)

If the arithmetic monodromy action is the natural degree-`n` action of `S_n`
with `n>=3`, then

\[
 \operatorname{Aut}_Y(\overline X)=1.                       \tag{4}
\]

The same holds for the natural action of `A_n` when `n>=4`.

**Proof.**  In either case the sheet stabilizer is the subgroup fixing one
letter.  It fixes exactly that letter globally, so an element normalizing it
must preserve that letter and hence already belongs to the stabilizer.  Thus
`N_G(H)=H`, and Theorem 1.1 applies.  The exclusions matter: the natural
`S_2` and regular `A_3` actions have nontrivial centralizer.  QED

The same conclusion follows from **geometric** natural `S_n` monodromy for
`n>=3`: the arithmetic permutation group contains the geometric `S_n` and is
therefore the same `S_n`.  If the geometric monodromy is natural `A_n` for
`n>=4`, the arithmetic group is `A_n` or `S_n`; both have trivial centralizer
in their natural action.  Thus either the arithmetic or geometric version of
the usual full-monodromy hypothesis suffices.

Full `S_n` is only a convenient sufficient condition.  The correct general
hypothesis is self-normality of the sheet stabilizer; many nonsymmetric
primitive actions also satisfy it, while primitivity alone does not.

## 2. Rigidity of distinguished affine opens

Let `U subset bar(X)` be a dense open over which (1) restricts to a
quasi-finite source model `F:U -> Y`.

### Corollary 2.1 (target-fixed rigidity of one source model)

If `N_G(H)=H`, then

\[
 \operatorname{Aut}_Y(U)=1.                                 \tag{5}
\]

**Proof.**  A target-fixed automorphism of `U` induces a `K`-automorphism of
its function field `L`.  Theorem 1.1 makes this field automorphism the
identity.  Two morphisms from the integral scheme `U` to the separated scheme
`U` which agree at the generic point agree everywhere.  QED

There is a stronger comparison formulation.  Let `U_1,U_2 subset bar(X)` be
two dense affine opens, each regarded as a source model over the same `Y`.

### Corollary 2.2 (open-marking criterion)

If `N_G(H)=H`, then

\[
 U_1\simeq_Y U_2\quad\Longleftrightarrow\quad U_1=U_2
 \text{ as open subschemes of }\overline X.                 \tag{6}
\]

When an isomorphism exists, it is the identity restriction.

**Proof.**  A `Y`-isomorphism induces a `K`-automorphism of `L`, hence the
identity.  As a rational map of `bar(X)` it is therefore the identity.  Its
restriction gives `U_1 subset U_2`; applying the inverse gives the reverse
inclusion.  QED

More generally, for two possibly different extensions `L_1/K` and `L_2/K`,
the set of target-fixed isomorphisms between their normalizations is either
empty or a torsor under their deck groups.  Trivial deck groups make such an
isomorphism unique if it exists, but do not make the set nonempty.  This is
the precise limitation of monodromy-only rigidity.

## 3. Boundary interpretation

Under the open-marking criterion, a source model is not extra generic-cover
data: it is a distinguished open inside a rigid finite normalization.  Since
the reduced boundary is the complement `bar(X)-U`, equation (6) is equivalent
to equality of the boundary marking.  Thus a target-fixed equivalence cannot
permute filled and missing boundary primes when the monodromy centralizer is
trivial.

This conclusion has no unresolved boundary-extension problem.  A field
automorphism extends uniquely across the boundary because `bar(X)` is the
integral normalization; conversely, an isomorphism of source opens is already
determined by its generic-field action.  What remains family-specific is
identifying the two opens inside a common normalization.

## 4. Application to cancellation parameter branches

For fixed `(m,r,C)`, every cancellation branch has the same primitive
extension

\[
 \Psi_{P,Q,R}(T)
 =C\int_0^T\{1-t(Q-Pt)^m\}^r\,dt-R.                         \tag{7}
\]

Its natural monodromy is `S_N` for odd `r` and `A_N` for even `r`, where
`N=r(m+1)+1`.  Here `N>=3` in the symmetric case and `N>=5` in the alternating
case, so Corollary 1.2 makes the common normalization target-fixed rigid.

Different parameter roots therefore do not define different generic covers.
They define different affine reconstruction opens—or equivalently, different
choices of the projective `P=0` branch filled by the source divisor `A=0`.
There is only one possible target-fixed birational identification: the one
fixing the common primitive root `s`, `P`, and `Q`.  In source coordinates it
is

\[
 z'=z+\frac{y^{m+1}\{h(A)-h'(A)\}}{A^{r+1}}.                \tag{8}
\]

The parameter-specific calculation is still needed, but only to test the
open-marking criterion: (8) preserves affine three-space exactly when

\[
 h(A)\equiv h'(A)\pmod {A^{r+1}}.                           \tag{9}
\]

Thus monodromy proves uniqueness and boundary rigidity; the valuation or
coordinate calculation proves inequality of the marked opens for distinct
normalized parameter roots.

## 5. Stabilization

After adjoining identity coordinates and fixing the enlarged target, the
Galois closure becomes `M(t_1,...,t_a)/K(t_1,...,t_a)` with the same group and
the same sheet stabilizer.  Theorem 1.1 and Corollaries 2.1--2.2 therefore
remain valid.  This does not address stable equivalences using a nonidentity
target automorphism; those belong to unrestricted stable left--right
equivalence.
