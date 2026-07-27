# Imprimitive Keller monodromy and polynomial factorization

This note gives explicit decomposable Keller maps, identifies their
monodromy block systems, and isolates the extra algebra needed to turn an
arbitrary intermediate function field into polynomial factors.

Work over a characteristic-zero field.  For `N\ge3`, let `F_N` be the
determinant-one, degree-`N` weighted map constructed from

\[
 H_N(W)=\frac{W^2-W^N}{N-2}.
\]

Explicitly, put

\[
 u=1+xy,\qquad
 \gamma_N=1-\frac{N}{N-1}xy+x^2z
\]

and define

\[
\boxed{
 F_N(x,y,z)=
 \left(
 \frac{(N-2)u+u^2-(N-1)u^N\gamma_N^{N-2}}{(N-2)x^2},
 \frac{(N-2)+2u-Nu^{N-1}\gamma_N^{N-2}}{(N-2)x},
 x\gamma_N
 \right).
}                                                       \tag{1}
\]

The weighted polynomiality theorem proves that the first two apparent
quotients are polynomials.  The
[primitive-monodromy atomicity theorem](PRIMITIVE_MONODROMY_ATOMICITY.md)
gives

\[
 \det DF_N=1,\qquad
 \deg_{\rm geom}F_N=N,\qquad
 \operatorname{Mon}_{\rm geom}(F_N)=S_N,
\]

and says that every `F_N` is absolutely and stably atomic.

## 1. Explicit decomposable maps

For integers `a,b\ge3`, define

\[
 \boxed{C_{a,b}=F_b\circ F_a.}                         \tag{2}
\]

This is an explicit polynomial map: formula (1), applied first with `N=a`
and then with `N=b`, is a compact straight-line generator for all three
coordinates.  The chain rule and multiplicativity of generic degree give

\[
 \det DC_{a,b}=1,\qquad
 \deg_{\rm geom}C_{a,b}=ab.                            \tag{3}
\]

Both displayed factors are noninvertible because their geometric degrees
are `a,b>1`.  Thus `C_{a,b}` is genuinely non-atomic.

The first three composite benchmarks have the decomposable representatives

\[
 C_{3,4}=F_4\circ F_3\quad(\deg=12),
\]

\[
 C_{3,5}=F_5\circ F_3\quad(\deg=15),
\qquad
 C_{4,5}=F_5\circ F_4\quad(\deg=20).                   \tag{4}
\]

In each of these degrees, `F_{ab}` from (1) is atomic while `C_{a,b}` is
decomposable.  Since stable polynomial left--right equivalence preserves
atomicity,

\[
 \boxed{F_{ab}\text{ and }C_{a,b}\text{ are not stably polynomially
 left--right equivalent}.}                            \tag{5}
\]

This supplies a same-degree stable separation which uses only composition
and monodromy, without a boundary-stratum calculation.

## 2. The block system of a composite

Let

\[
 K=k(C_{a,b}),\qquad
 E=k(F_a),\qquad
 L=k(x,y,z).
\]

The composition gives a strict tower

\[
 K\subsetneq E\subsetneq L,\qquad
 [E:K]=b,\quad [L:E]=a.                                \tag{6}
\]

On a geometric generic fiber of `C_{a,b}`, restriction to the intermediate
field partitions the `ab` sheets into `b` blocks of size `a`.  Therefore its
geometric monodromy `G_{a,b}` is transitive and imprimitive, with

\[
 \boxed{G_{a,b}\hookrightarrow S_a\wr S_b}             \tag{7}
\]

in the product action on `b` blocks of size `a`.

The action induced on the set of blocks is the monodromy of the outer map
`F_b`, hence is `S_b`.  The action induced by a block stabilizer on one
chosen block is the monodromy of the inner map `F_a`, hence is `S_a`.
Thus `G_{a,b}` has full top group and full local group:

\[
 G_{a,b}\twoheadrightarrow S_b,\qquad
 (G_{a,b})_{\mathcal B}\twoheadrightarrow S_a.          \tag{8}
\]

Equations (7)--(8) do **not** by themselves prove

\[
 G_{a,b}=S_a\wr S_b.                                    \tag{9}
\]

The missing assertion is independence of the `b` conjugate inner Galois
closures.  Equivalently, one must prove that the kernel of the block action
is the full base group `S_a^b`, rather than a diagonal or parity-coupled
subdirect product.  A sufficient geometric certificate would be inertia
supported in one block: one within-block transposition fixing all other
blocks, together with `S_b` conjugacy, generates the full base group.  No
such divisor certificate follows formally from composition alone.
Accordingly, (7)--(8) are the general theorem and (9) is a separate
case-by-case problem.  The first case is resolved next.

For `C_{3,4}`, the single-block test has a small exact formulation.  On the
chart where the third intermediate coordinate `c` is nonzero, the normalized
repeated-root divisor of the inner cubic is

\[
 \nu_3(r,c)=
 \left(
 \frac{r^2(1-2r)}{c^2},
 \frac{r(2-3r)}c,
 c
 \right).                                               \tag{9a}
\]

Consider

\[
 F_4\circ\nu_3:\mathbb A^1_r\times\mathbb G_{m,c}
 \longrightarrow\mathbb A^3
\]

and write `(A,B,C)` for its target coordinates.  The outer marked-root
equation, cleared by the scalar two, is

\[
 E_4(V)=V^2-V^4-2BCV+2AC^2.                            \tag{9b}
\]

Pull the non-`c=0` cubic discriminant through the outer reconstruction,
clear only powers of the reconstruction units, reduce modulo (9b), and
saturate by `C`.  The reduced condition is cubic in `V`, with 95 terms.
Its exact resultant with (9b) factors as

\[
 \operatorname{Res}_V(E_4,R_3)=C^8Q(A,B,C),             \tag{9c}
\]

where `Q` is irreducible over `\mathbb Q`, has degrees

\[
 (\deg_AQ,\deg_BQ,\deg_CQ)=(21,28,22),
\]

has 1001 terms, and occurs with exponent one.  The factor `C^8` is removed
by the declared saturation.  Since the incidence over the irreducible
factor `Q` has resultant multiplicity one, the quartic and cubic have
exactly one common simple root at the generic point of `Q`.  Therefore
`F_4\circ\nu_3` is generically birational onto its image.

It remains necessary to distinguish the image of the other inner
discriminant component `c=0`.  At `(r,c)=(2,1)`, the intermediate point is
`(-12,-8,1)` and its outer target is

\[
 \left(-\frac{799529969}{3},\,43960408,\,-204\right).
\]

The specialized quartic (9b) is squarefree, has marked root `V=1649`, and
is coprime to the numerator of the reconstructed third intermediate
coordinate.  Hence none of its four intermediate sheets has `c=0`, proving
that the two ramification images are distinct.

Generic `Q`-inertia is consequently a transposition supported in exactly
one of the four three-sheet blocks.  Conjugating by the full local `S_3`
produces the whole `S_3` in that block; conjugating by the top `S_4`
produces `S_3^4`.  Together with the surjection to `S_4`, this proves

\[
\boxed{
 G_{3,4}=S_3\wr S_4,\qquad |G_{3,4}|=6^4\cdot24=31104.
}                                                       \tag{9d}
\]

The exact reduction, saturated factorization, and boundary-separation gcds
are checked by
[`verify_degree_twelve_wreath_elimination.py`](../scripts/verify_degree_twelve_wreath_elimination.py).

## 3. Why imprimitivity is not yet polynomial factorization

Let

\[
 F:\mathbb A^d\longrightarrow\mathbb A^d
\]

be any Keller map, and abbreviate

\[
 A=k[F_1,\ldots,F_d]\subseteq
 B=k[x_1,\ldots,x_d].
\]

Because `F` is dominant, both `A` and `B` are polynomial rings of
transcendence degree `d`.  Let

\[
 \operatorname{Frac}(A)\subsetneq E\subsetneq
 \operatorname{Frac}(B)                                \tag{10}
\]

be an intermediate field, for example one supplied by an imprimitive block
system.

> **Polynomial-sandwich criterion.**  
> The field `E` is realized by a polynomial factorization
> \[
> F=G\circ H,\qquad
> \operatorname{Frac}k[H_1,\ldots,H_d]=E,               \tag{11}
> \]
> if and only if there is a subalgebra
> \[
> \boxed{A\subseteq R\subseteq B,\qquad
> R\simeq k[t_1,\ldots,t_d],\qquad
> \operatorname{Frac}(R)=E.}                            \tag{12}
> \]
> If the inclusions in (10) are strict, both factors in (11) are
> noninvertible Keller maps.

### Proof

If (11) exists, take `R=k[H_1,\ldots,H_d]`.  Dominance of `H` makes its
components algebraically independent, so `R` is a polynomial ring.
The equations `F_i=G_i(H_1,\ldots,H_d)` give `A\subseteq R\subseteq B`,
and its fraction field is `E`.

Conversely, choose an isomorphism

\[
 k[t_1,\ldots,t_d]\xrightarrow{\sim}R
\]

and write `h_i\in B` for the images of `t_i`.  Since `A\subseteq R`, each
`F_i` is a polynomial `g_i(h_1,\ldots,h_d)`.  The maps

\[
 H=(h_1,\ldots,h_d),\qquad G=(g_1,\ldots,g_d)
\]

then satisfy `F=G\circ H` and realize `E`.  The chain-rule argument in the
atomicity theorem forces both `G` and `H` to be Keller.  Their geometric
degrees are the two field degrees in (10), so strictness makes both
noninvertible.  This proves the criterion.

The condition `E\simeq k(t_1,\ldots,t_d)` alone is not enough.  It only says
that the intermediate field is rational.  Polynomial factorization requires
a **single polynomial affine-space model** `R` which is simultaneously
inside the source coordinate ring and contains the target coordinate ring.
This two-sided regularity is the missing algebraization condition.

## 4. Normalization-open interpretation

Let `Z_E` be the normalization of the target `\mathbb A^d` in `E`.  It is
finite over the target.  If the polynomial sandwich (12) exists, the
intermediate affine space

\[
 U=\operatorname{Spec}R\simeq\mathbb A^d
\]

maps birationally and quasi-finitely to `Z_E`; normality and Zariski's Main
Theorem identify it with an open subset

\[
 U\hookrightarrow Z_E.                                 \tag{13}
\]

The rational lift of the source to `Z_E` is regular on all of
`\mathbb A^d` and lands in this affine-space reconstruction open.  Conversely,
such an open (13), together with regularity of the source lift into it,
recovers the ring sandwich (12).

Thus an imprimitive Keller cover factors polynomially precisely when its
intermediate normalization admits a compatible affine-space reconstruction
open.  There are three logically separate gates:

1. **group gate:** an imprimitive block system, giving `E`;
2. **rationality gate:** `E` is a rational function field;
3. **regularity gate:** the polynomial sandwich (12), equivalently the
   compatible reconstruction open (13).

For one-variable polynomial covers, Lüroth's theorem and total ramification
at infinity make the last two gates automatic; this is why imprimitivity and
Ritt decomposition are equivalent there.  In dimensions at least two they
are genuine additional conditions.

## 5. Research consequences

The comparison between `F_{ab}` and `C_{a,b}` gives a controlled laboratory:

| map | degree | monodromy | intermediate field | composition |
|---|---:|---|---|---|
| `F_{ab}` | `ab` | `S_{ab}` | none | atomic, stably atomic |
| `C_{a,b}` | `ab` | imprimitive subgroup of `S_a wr S_b` | `k(F_a)` | `F_b o F_a` |

The next exact questions are now sharply separated.

1. Determine whether `G_{a,b}` is the full wreath product for the remaining
   explicit pairs in (4), starting with `(a,b)=(3,5)`.
2. Search for primitive-to-imprimitive degenerations whose intermediate
   normalization fails the polynomial-sandwich criterion.  Such examples
   would prove concretely that imprimitive monodromy is not sufficient for
   polynomial decomposition in the Keller setting.
3. Classify all polynomial-ring sandwiches between `k[C_{a,b}]` and
   `k[x,y,z]`.  This is the multivariable analogue of Ritt factorization and
   asks whether the displayed block system is the unique compositional
   structure.
4. Compare the two orders `C_{a,b}` and `C_{b,a}`.  Their block sizes differ,
   so their monodromy systems can distinguish them unless a common refinement
   or Ritt-type move exists.

The generator and exact polynomial regressions for `C_{3,4}` are in
[`scripts/verify_composite_degree_twelve.py`](../scripts/verify_composite_degree_twelve.py).
