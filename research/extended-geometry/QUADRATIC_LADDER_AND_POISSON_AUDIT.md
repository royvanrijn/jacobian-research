# Quadratic ladder and rank-two Poisson audit

This note audits the two supplied cards separately.  Their evidentiary status
is different: the quadratic marked-root ladder is exactly identifiable from
the displayed formulas, while the rank-two Poisson card omits most of the
polynomials needed to check its claim.

## The ladder is the `m=1` cancellation column

To avoid conflicting indices, write `ell` for the rung on the card.  Set

\[
 D(s)=1-Bs+As^2,
 \qquad
 P_\ell(\lambda)=-J_\ell\int_0^\lambda D(s)^{\ell-1}\,ds-C.
\]

Then

\[
 P_\ell'(\lambda)=-J_\ell D(\lambda)^{\ell-1},
 \qquad \deg P_\ell=2\ell-1.
\]

This is precisely the generalized cancellation inverse polynomial with

\[
 (m_{\rm cancel},r)=(1,\ell-1),
 \qquad N=r(m_{\rm cancel}+1)+1=2\ell-1.
\]

Thus the card supplies a useful geometric interpretation, not a disjoint new
family: a target is a degree-`2 ell - 1` polynomial with one simple root
marked, and its two stationary points have multiplicity `ell - 1`.  For
`ell=2`, substituting `A=3a`, `B=b`, and `C=-3c`, then dividing by three,
gives exactly

\[
 c-2s+bs^2-2as^3,
\]

the repository's foundational marked-root cubic.

This identification has two useful consequences.

1. The generic degree and inequivalence of different rungs are already
   covered by the cancellation degree theorem.
2. The full source-boundary intersection diagram is certified for every
   rung, not merely in a bounded experiment: the `m=1` branch polynomial is
   irreducible, and the critical-contact coefficient is not an associate, so
   their resultant is nonzero.

The exact regression is
[`verify_counterexample_ladder.py`](../scripts/verify_counterexample_ladder.py).

## What is new in the Poisson card

The first card claims a standard rank-two Poisson endomorphism on
`Q[x,q,p,z]`, with four output polynomials `(R,T,D,S)`, and displays only

\[
 R=x(2-3xq).
\]

This posed the problem of improving the repository's cotangent-lift
consequence from three canonical pairs to two.  The independent construction
below now solves that mathematical problem, while the card's provenance and
actual missing outputs remain unaudited.  It also needs a terminology
explanation.  In the usual convention `A_n` has `2n` associated-graded
polynomial generators.  The four displayed Poisson coordinates have the
symbol size of `A_2`.  The abstract does not say that their direct
quantization automatically gives `A_4`: it advertises a separate construction
using the same four polynomials as commuting base-coordinate multiplication
operators, together with four Hamiltonian duals.  Read as an
inverse-Jacobian/cotangent construction, that
produces eight Weyl generators and therefore lands naturally in `A_4`.  The
full appendix is still needed to check that this is the intended construction
and that its operator identities and non-surjectivity proof are correct.

The screenshot alone is not a certificate.  To admit the claim, the audit
needs the exact formulas for `T,D,S` and must independently verify

\[
 \{D,R\}=1,\quad \{S,T\}=1,
 \quad \{R,S\}=\{R,T\}=\{D,S\}=\{D,T\}=0,
\]

plus a displayed collision or another proof of nonsurjectivity.  The proposed
Weyl lift must then be checked at the operator level; principal-symbol
identities alone do not prove all Weyl commutators.

## An exact foundational fingerprint

One part of the provenance question can be answered without guessing the
missing outputs.  Put

\[
 \begin{aligned}
 X&=x,\\
 Z&=3x^2p+(2-6xq)z,\\
 Y&=q-\frac{xZ}{3},\\
 E&=\frac{1+3xq}{2}p-3q^2z.
 \end{aligned}                                           \tag{1}
\]

This is a polynomial automorphism of `Q[x,q,p,z]`, with determinant `-1`.
Indeed, writing `q=Y+XZ/3`, its inverse is

\[
 p=-6EXq+2E+3Zq^2,
 \qquad
 z=-3EX^2+\frac32ZXq+\frac12Z.                         \tag{2}
\]

For the foundational map `F=(F_1,F_2,F_3)`, substitution of (1) gives

\[
 F_3(X,Y,Z)=2X-3X^2Y-X^3Z=x(2-3xq)=R.                 \tag{3}
\]

Thus the displayed `R` is an exact algebraic fingerprint of the foundational
counterexample, not merely a coincident degree or fiber count.  It is strong
evidence that the announced construction is derived from that map.  It is not
yet a proof that the unavailable `T,D,S` are the manuscript's particular
lift, nor does it identify the author or source.

There is also a useful negative result.  The obvious invariant choices

\[
 S_0=\frac12F_1(X,Y,Z),\qquad T_0=F_2(X,Y,Z)
\]

satisfy exactly

\[
 \{S_0,T_0\}=1,
 \qquad \{R,S_0\}=\{R,T_0\}=0,
 \qquad \{E,R\}=1.                                    \tag{4}
\]

Nevertheless there is **no** polynomial correction
`D=E+f(X,Y,Z)` for which both `{D,S_0}` and `{D,T_0}` vanish.  Since (1) is a
polynomial coordinate system and the kernel of `{-,R}` is `Q[X,Y,Z]`, this
covers every possible polynomial partner `D` for this particular pair
`S_0,T_0`.

Here is the exact obstruction.  With `Q=Y+XZ/3`, the first completion
equation can be integrated explicitly; the two remaining equations reduce to

\[
 \left(-3X^2\partial_X+(6XQ-2)\partial_Q\right)h
 =\frac{Q^3}{2}
  \left(54Q^3X^3+189Q^2X^2+222QX+89\right).           \tag{5}
\]

After localizing at `X`, set `v=1/X` and
`rho=2X-3X^2Q`.  The derivation on the left sends `v` to `3` and `rho` to
zero.  Its forced antiderivative, modulo a polynomial in `rho`, becomes

\[
 \frac{54Q^6X^6+234Q^5X^5+375Q^4X^4+270Q^3X^3
       +90Q^2X^2+24QX+4}{60X^4}.                      \tag{6}
\]

A polynomial in `rho` has no negative power of `X`, so it cannot cancel the
term `1/(15X^4)`.  Equation (5) has no polynomial solution.  This proves that
the manuscript's missing formulas, if correct, contain a genuinely different
choice or construction; they cannot responsibly be filled in by the naive
foundational substitution.  The exact certificate is
[`verify_rank_two_poisson_preaudit.py`](../scripts/verify_rank_two_poisson_preaudit.py).

## The completion criterion

The negative result above suggests the right intrinsic problem.  Work on the
invariant quotient in coordinates `(X,Q,Z)`, where

\[
 Q=Y+XZ/3,\qquad R=2X-3X^2Q.
\]

For any Keller triple `A=(S,T,R)` with determinant `-1`, let `w_A` be the
unique polynomial derivation satisfying

\[
 w_A(S)=w_A(T)=0,\qquad w_A(R)=1.
\]

In the adapted four-dimensional coordinates, put `w_E={E,-}`.  Every
polynomial satisfying `{D,R}=1` has the form `D=E+f(X,Q,Z)`, because
`{-,R}=partial_E`.  The two remaining mixed brackets vanish exactly when

\[
 w_A-w_E=\{f,-\}.                                      \tag{7}
\]

Thus polynomial Darboux completion is a relative Hamiltonian, or polynomial
flux, problem.  This criterion is both necessary and sufficient; it avoids
guessing `D` coefficient by coefficient.

Apply it after the volume-preserving, `R`-fixing shear

\[
 Z\longmapsto Z+cQ^2.
\]

After integrating the first equation in (7), use
`v=1/X` and `rho=R` as in (6).  The entire negative-`X` principal part of the
forced remaining primitive is

\[
 (c+9)\left(
 \frac1{135X^4}+\frac{2Q}{45X^3}
 +\frac{Q^2}{6X^2}+\frac{Q^3}{2X}
 \right).                                               \tag{8}
\]

Consequently `c=-9` is the unique pole-free member of this one-parameter
family.  Unlike the naive case, its forced primitive is polynomial.  This is
how the missing completion was found.

### Relation to a degree-by-degree coefficient search

This closes the proposed reverse-engineering problem without requiring a
blind Gröbner search.  Its requested stages appear here as follows.

1. Equations (1)--(2) give polynomial coordinates with
   `ker{-,R}=Q[X,Q,Z]`: the derivation `{-,R}` kills `X,Q,Z`, sends `E` to
   one, and is therefore exactly `partial_E` in this polynomial coordinate
   ring.  The induced quotient bracket is displayed explicitly in the proof
   below.
2. The foundational triple after the `R`-fixing shear `Z -> Z+cQ^2`
   automatically supplies candidates `(S_c,T_c)` with `{S_c,T_c}=1` and both
   commuting with `R`.  This replaces a general coefficient ansatz by a
   geometrically constrained family that already transports the complete
   foundational fiber.
3. Once `(S_c,T_c)` is fixed, writing `D=E+f` turns the last two mixed-bracket
   equations into the linear Hamiltonian equation (7).  Integration leaves
   the principal part (8), so polynomiality forces `c=-9`; substituting it
   gives the polynomial `K` in (11).
4. Equations (12)--(14) then verify all six brackets and carry the exact
   three-point collision.

Thus existence is settled constructively.  What has **not** been proved is a
degree-bounded classification of every possible canonical pair `(S,T)`, or
that this solution equals the unavailable manuscript's outputs.  A general
Gröbner/SAT coefficient search would now address uniqueness and classification,
not the existence of a rank-two completion.

## An independent rank-two Poisson completion

The mathematical existence gap can now be closed without claiming to have
recovered the unavailable manuscript.  Starting from (1), set

\[
 Q=q,\qquad W=Z-9Q^2,\qquad Y=Q-\frac{XW}{3},
 \qquad U=1+XY.                                         \tag{9}
\]

Define

\[
 \begin{aligned}
 R&=2X-3X^2Y-X^3W=x(2-3xq),\\
 S&=\frac12\left(U^3W+Y^2U(4+3XY)\right),\\
 T&=Y+3XU^2W+3XY^2(4+3XY),                              \tag{10}
 \end{aligned}
\]

and put `D=E+K(X,Y,W)`, where

\[
\begin{aligned}
60K={}&-10X^2W^3-90XYW^2-20W^2
       +18X^3Y^5W+90X^2Y^4W\\
     &+180XY^3W-90Y^2W
       +54X^2Y^6+234XY^5+375Y^4.                       \tag{11}
\end{aligned}
\]

All of (9)--(11) are polynomial expressions in the original variables.  In
particular, these are compact exact formulas for all four outputs, not an
existence ansatz.

### Theorem

With the bracket convention `{p,x}={z,q}=1`, the assignment

\[
 (x,q,p,z)\longmapsto(R,T,D,S)                          \tag{12}
\]

satisfies

\[
 \{D,R\}=1,\quad \{S,T\}=1,
 \quad
 \{R,S\}=\{R,T\}=\{D,S\}=\{D,T\}=0.                  \tag{13}
\]

It has Jacobian determinant one, is not injective, and has generic degree
three.  The target `(R,T,D,S)=(0,0,0,-1/8)` has exactly the three points

\[
 \left(0,0,\frac1{24},-\frac18\right),\quad
 \left(1,\frac23,\frac{247}{96},-\frac{89}{64}\right),\quad
 \left(-1,-\frac23,\frac{247}{96},-\frac{89}{64}\right)                 \tag{14}
\]

in its fiber.  Hence (12) is a counterexample to the Poisson conjecture for
two canonical pairs and a noninjective exact polynomial symplectic etale map
of `A^4`.

### Proof

The coordinate system `(X,Q,Z,E)` is polynomial by (1)--(2).  The base change

\[
 (X,Q,Z)\longmapsto(X,Y,W)
\]

in (9) is also a polynomial automorphism, with inverse

\[
 Q=Y+XW/3,\qquad Z=W+9Q^2.                              \tag{15}
\]

Finally `E -> E+K` is triangular.  Therefore

\[
 H:(x,q,p,z)\longmapsto(X,Y,W,D)
\]

is a polynomial automorphism.  If `F=(F_1,F_2,F_3)` is the foundational map,
then (10) gives the exact factorization

\[
 (S,T,R,D)=(F_1/2,F_2,F_3,\operatorname{id})\circ H.    \tag{16}
\]

Thus (12), after a target-coordinate permutation, is polynomially
right--left equivalent to `(F_1/2,F_2,F_3) x id`.  It follows at once that its
generic degree is three and that all its fiber schemes are the corresponding
foundational fiber schemes.  Applying (15), then (2), to the complete
foundational fiber over `(-1/4,0,0)` gives exactly (14), with `D=0`.

It remains to prove that the particular source automorphism `H` makes (12)
canonical.  In `(X,Q,Z)` the induced bracket is

\[
 \{X,Q\}=0,\qquad \{X,Z\}=-3X^2,
 \qquad \{Q,Z\}=6XQ-2,
\]

and

\[
 \{E,X\}=\frac{1+3XQ}{2},\qquad
 \{E,Q\}=-3Q^2,\qquad
 \{E,Z\}=\frac92QZ.
\]

The determinant of `(S,T,R)` with respect to `(X,Q,Z)` is `-1`, which gives
`{S,T}=1` and the two brackets with `R`.  Substitution of (11) into the flux
identity (7) gives the remaining three relations in (13).  These are
polynomial identities; hence `det d(R,T,D,S)=1`.  Preservation of the
nondegenerate Poisson tensor is equivalent to preservation of the symplectic
form.  The difference of the pulled-back and original polynomial canonical
one-forms is closed and therefore exact on affine four-space by the polynomial
Poincare lemma.  This proves the exact symplectic statement.

Every algebraic step, all six brackets after full substitution into
`Q[x,q,p,z]`, the determinant, the inverse coordinate changes, and (14) are
checked by
[`verify_rank_two_poisson_completion.py`](../scripts/verify_rank_two_poisson_completion.py).
A separate standard-library sparse-polynomial implementation,
[`audit_rank_two_poisson_completion_independent.py`](../scripts/audit_rank_two_poisson_completion_independent.py),
rebuilds (9)--(11) over `Fraction` and independently checks all six brackets,
the determinant, expanded term counts, and (14).

## The completion as a four-dimensional cotangent-graph restriction

There is a direct relation with the six-dimensional cotangent lift, stronger
than the left--right equivalence (16).  It is a restriction to two
four-dimensional symplectic graphs, not an invariant hypersurface or a
Marsden--Weinstein quotient.

Write

\[
 F_0=(S_0,T_0,R_0)=(F_1/2,F_2,F_3),\qquad
 J_0=D F_0,
\]

on the base `A^3_b`, and use fiber coordinates `alpha` on its cotangent
bundle.  On the target cotangent bundle use base coordinates `(s,t,r)` and
fiber coordinates `beta`.  Define closed embeddings

\[
\begin{aligned}
 i_{\rm src}(b,\lambda)&=\left(b,J_0(b)^T(0,S_0(b),\lambda)^T\right),\\
 i_{\rm tgt}(s,t,r,\lambda)&=((s,t,r),(0,s,\lambda)).
\end{aligned}                                                    \tag{17}
\]

They are closed because `det J_0=-1`: on the source image one recovers
`(0,S_0,lambda)^T=J_0^{-T}alpha`, while the target image is cut out by
`beta_s=0`, `beta_t=s`.  If `widehat F_0` is the cotangent lift, then

\[
 \boxed{
 \widehat F_0\circ i_{\rm src}
 =i_{\rm tgt}\circ(F_0\times\operatorname{id}_{\mathbb A^1}).}
                                                                  \tag{18}
\]

Indeed, the momentum component of the left side is
`J_0^{-T}J_0^T(0,S_0,lambda)^T=(0,S_0,lambda)^T`.

This restriction is symplectic.  For the cotangent one-form `theta`,

\[
 i_{\rm src}^*\theta=S_0\,dT_0+\lambda\,dR_0,
 \qquad
 i_{\rm tgt}^*\theta=s\,dt+\lambda\,dr.              \tag{19}
\]

Their exterior derivatives are nondegenerate: the square of the first is a
nonzero constant multiple of
`dS_0 wedge dT_0 wedge d(lambda) wedge dR_0`, because `det J_0=-1`; the target
statement is immediate.  Equation (18) preserves the one-forms in (19)
exactly.

Thus every three-variable Keller map has a canonical four-dimensional
cotangent-graph restriction, but generally its source symplectic form need
not admit a polynomial Darboux trivialization.  For the foundational map,
the relative-flux calculation (7)--(11) supplies exactly such a
trivialization: the unique pole-canceling shear `Z -> Z-9Q^2` and the
Hamiltonian `K` turn the source graph form into the standard form on
`A^4_(x,q,p,z)`.  After the target permutation
`(s,t,r,lambda) -> (r,t,lambda,s)`,
equation (18) is precisely the map `(R,T,D,S)` in (12).

This locates the successful route among the proposed reductions:

\[
 \boxed{\text{cotangent graph restriction}
 \; + \; \text{polynomial relative-Hamiltonian trivialization}.}
                                                                  \tag{20}
\]

The same exact checker verifies (18)--(19) symbolically.

## The `DC(4)` consequence and the classical implication direction

Treat the four polynomials in (12) as a four-variable Keller map `G`.  The
inverse-Jacobian construction sends the four multiplication generators of
`A_4` to `G_i` and its four differential generators to

\[
 \delta_i=\sum_r(DG^{-1})_{ri}\partial_r.
\]

The eight images satisfy the Weyl relations.  Exact preservation of
differential order shows that surjectivity would force a polynomial inverse
to `G`; (14) excludes one.  Therefore this is an injective non-surjective
endomorphism of `A_4`.  This explains why the Weyl rank is four, not two: it
is the cotangent/inverse-Jacobian quantization of the four-variable base map,
not an automatic ordering of the four Poisson coordinates inside `A_2`.
The general proof is [Theorem 2 of the symplectic/Weyl note](SYMPLECTIC_WEYL_LIFT.md).

The standard equivalence theorem does not itself produce this rank-two
completion from JC(3).  Theorem 7 of Adjamagbo--van den Essen,
[*On the equivalence of the Jacobian, Dixmier and Poisson Conjectures in any
characteristic*](https://arxiv.org/abs/math/0608009), gives the fixed-dimension
chain

\[
 CJC(2n)\Longrightarrow CPC(n)\Longrightarrow CDC(n)\Longrightarrow CJC(n).
\]

For `n=2`, its first arrow starts from `JC(4)`, not `JC(3)`.  Formulas
(9)--(11) use the special marked-root structure of the foundational map and
are not a formal corollary of that general theorem.

## Provenance scope (`OP-QP-PROV`, parked)

The **mathematical** existence, bracket, collision, symplectic, and `DC(4)`
gaps are now closed by the independent construction above.  The
**provenance** gap is not: the external paper, authors, version, and its actual
`T,D,S` remain unavailable, so no equality between (9)--(11) and the announced
manuscript is asserted.  This provenance question is tracked only as parked
item `OP-QP-PROV` in [STATUS.md](../STATUS.md).

No public source matching the exact text was located on 22 July 2026 after
exact-phrase and formula searches, arXiv metadata searches, author/repository
searches, and inspection of the public Omniscience Project papers index.  Its
separately published `A_3` paper is not this manuscript.  The supplied item
therefore remains an **external announced manuscript under provenance audit**,
not an attributed Long result and not an external review of this repository.
If its source appears, its outputs should be compared term-by-term with
(9)--(11), and its Weyl appendix should be audited independently.

The broader consequence comparison is in
[External consequences and provenance](EXTERNAL_CONSEQUENCES_AND_PROVENANCE.md).
