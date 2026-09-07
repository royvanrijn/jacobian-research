# Rational solubility and residual Selmer theorems

This is the canonical soluble-cover companion to
[`SPECIALIZATION_QUOTIENT_AND_RANK_JUMP_THEOREMS.md`](SPECIALIZATION_QUOTIENT_AND_RANK_JUMP_THEOREMS.md).
That note supplies the rank-jump sandwich and quotient-height identities;
the present note identifies the rationally soluble subspace, constructs
actual covers for point classes, and specifies what their complexity can
measure. The statements are standard descent consequences and explicit
deductions, not a claim of a new rank-jump theorem for the exceptional fibres.

The target panel is ICARM curves
`351,356,376,377,385,398,400,401,542,543,548`. Its displayed points provide
soluble classes. They do not, by themselves, provide the complementary Selmer
classes or insoluble controls. No point search, class-group campaign, or
complete descent is required by the proofs below.

## 1. The exact object, including saturation

Let `E/K` be an elliptic curve over a number field. Write `A=E(K)`, let
`M subset A` be the displayed specialized generic subgroup of rank `r`, and
put

\[
 S=\operatorname{Sel}_2(E/K),\qquad
 G=\delta(M)\subset S,\qquad R_G=S/G.
\]

Here `delta(M)` means the image in `A/2A`, followed by the Kummer injection.
Writing `M/2M` as an injected subspace requires a separate injectivity check.
Define

\[
 W_G=\delta(A)/G\subset R_G,\qquad
 \epsilon_M=\dim_{\mathbf F_2}(A/M)[2].
\]

### Theorem R1: the residual exact sequence

There are natural identifications and an exact sequence

\[
 W_G\simeq A/(M+2A),\qquad
 0\longrightarrow W_G\longrightarrow R_G
 \longrightarrow\Sha(E/K)[2]\longrightarrow0.
 \tag{R1.1}
\]

Consequently, if `Delta=rank(A)-r`, then

\[
 \boxed{\dim W_G=\Delta+\epsilon_M,\qquad
 \dim R_G=\Delta+\epsilon_M+\dim\Sha(E/K)[2].}
 \tag{R1.2}
\]

A residual class belongs to `W_G` if and only if an associated 2-covering
has a `K`-rational point. This does not depend on which representative of
the residual class is chosen.

**Proof.** Quotient the Kummer--Selmer exact sequence
`0 -> A/2A -> S -> Sha[2] -> 0` by `G`. The first quotient is
`A/(M+2A)=(A/M)/2(A/M)`. A finitely generated abelian group `B` has
`dim(B/2B)=rank(B)+dim B[2]`, proving the formulas. The image in `Sha`
is the class of the underlying genus-one torsor. This class vanishes
precisely when the torsor has a rational point. Adding an element of `G`
does not change that torsor class. QED.

Thus the formula without `epsilon_M` holds when `A/M` has no 2-torsion.
Full saturation at odd primes is unnecessary. Unremoved rational 2-torsion
and 2-primary saturation defects both belong in this correction. Equivalently,

\[
 \epsilon_M=\dim A[2]+r-\dim G.
 \tag{R1.3}
\]

For the irreducible 2-division cubics in the existing R17 pressure theorem,
`A[2]=0`; a certified `dim G=r` therefore gives `epsilon_M=0`, exactly the
hypothesis used in that note. For a specialized family, identifying `Delta`
with the jump from the generic curve still requires preservation of generic
rank under specialization, as in Theorem S1.

An exact free quotient `L/M=Z^q` inside a displayed subgroup `L` does not
alone prove `dim(delta(L)/G)=q`: `L` may have even index in a larger subgroup
of `A`. Conversely, a point outside the rational span of `M` can have zero
residual Kummer class when it is twice another point. Rank independence and
independence modulo `M+2A` require their respective certificates. In
particular a certificate using reduction modulo three cannot silently be
relabelled as a mod-two certificate.

## 2. Cassels--Tate is an obstruction, not a solubility test

Identify `1/2 Z/Z` with `F_2`. The Cassels--Tate pairing pulls back from
`Sha[2]` to an alternating pairing on `S`. Its radical is the image of
the natural map `Sel_4(E/K) -> Sel_2(E/K)`. These standard facts are stated
in Sections 2--3 of
[Fisher--Schaefer--Stoll, *The yoga of the Cassels--Tate pairing*](https://www.dpmms.cam.ac.uk/~taf1000/papers/casselspairing.pdf).
No finiteness assumption on the full Tate--Shafarevich group is needed here.

### Theorem R2: the residual radical retains divisible Sha classes

Since `G` pairs to zero with all of `S`, the pairing descends to `R_G`.
Write it as `lambda_G`. Then

\[
 W_G\subset\operatorname{rad}(\lambda_G),\qquad
 \boxed{\operatorname{rad}(\lambda_G)/W_G
        \simeq 2\Sha(E/K)[4].}
 \tag{R2.1}
\]

In particular the equality `rad(lambda_G)=W_G` is equivalent to
`2 Sha[4]=0`, not merely to finiteness of `Sha`.

**Proof.** The Kummer image pairs to zero because its image in `Sha` is
zero. The natural diagram for multiplication by two has rows

```text
0 -> A/4A -> Sel_4 -> Sha[4] -> 0
     |          |          | [2]
0 -> A/2A -> Sel_2 -> Sha[2] -> 0.
```

The left map is surjective. Hence the image of the middle map is exactly
the inverse image of `2 Sha[4]`. Apply the standard radical theorem and
quotient by `G`, and then by the full Kummer image. QED.

Two elementary examples make the distinction sharp. On `F_2^2` with
`lambda(e1,e2)=1`, the line `<e1>` is isotropic, but its nonzero element
is outside the radical. It cannot be a global point class for this pairing.
On the abstract group `(Z/4Z)^2`, the perfect alternating form
`((a,b),(c,d)) -> (ad-bc)/4 mod Z` restricts to zero on the whole
two-torsion subgroup. Thus even a finite perfect Cassels--Tate group can
have two nonzero-dimensional directions of `Sha[2]` invisible to the
restricted pairing. This is a group-theoretic example, not a claimed
Tate--Shafarevich computation for a panel curve.

### Corollary R3: useful pairing certificates

If `dim R_G=s` and the **complete** pairing matrix has rank `2a`, then

\[
 \dim W_G\leq s-2a,\qquad
 \operatorname{rank}E(K)\leq r+s-2a-\epsilon_M.
 \tag{R3.1}
\]

A certified nonzero pairing `lambda_G(x,y)` proves that both classes have
insoluble underlying torsors. A zero pairing proves neither one soluble.
The largest isotropic subspace has dimension `s-a`, whereas the radical has
dimension `s-2a`; replacing the latter by an isotropic subspace weakens the
condition in precisely the wrong way.

For a certified subspace `V subset R_G`, a restricted matrix of rank `2a`
only gives

\[
 \dim(W_G\cap V)\leq\dim V-2a.
\]

It supplies no full rank upper bound without a certified ambient upper
space. Its radical may exceed `V intersect rad(lambda_G)` because classes
outside `V` have not been paired. If `V` consists solely of known point
classes, every row pairs to zero with **all** Selmer classes by R2. This is
an exact global-witness certificate, but it provides no information about
the pairing between unknown complementary classes.

The componentwise Hilbert tensors already studied in the repository are
not automatically Cassels--Tate matrices. The actual Cassels construction
requires auxiliary functions on the covering, suitable local points, and
the prescribed product of local symbols; see
[Shukla--Stoll, *The Cassels--Tate pairing on 2-Selmer groups of elliptic curves*](https://arxiv.org/abs/2302.01640).
A locally insoluble ambient squareclass is outside `Sel_2` and must not
serve as an insoluble **Selmer** control.

## 3. Explicit actual 2-coverings of known point classes

Assume characteristic zero and fix a short model

\[
 E:y^2=x^3+Ax+B,\qquad P=(a,b)\in E(K).
\]

Here `A,B` denote curve coefficients, rather than the group `A` of Section 1.

### Theorem R4: the pointed quartic with its covering map

Let `C_P` be the smooth projective model of

\[
 w^2=t^4-6at^2-8bt-3a^2-4A.
 \tag{R4.1}
\]

The formulas

\[
 x(\phi_P(t,w))=\frac{t^2-a+w}{2},\qquad
 y(\phi_P(t,w))=\frac{t(t^2-3a+w)}2-b
 \tag{R4.2}
\]

extend to a `K`-isomorphism of curves `phi_P:C_P -> E`. Equip this curve
with the **degree-four** map

\[
 \boxed{\pi_P=[2]\circ\phi_P-P:C_P\longrightarrow E.}
 \tag{R4.3}
\]

Then `(C_P,pi_P)` is the 2-covering representing `delta(P)`. It has rational
points at both infinities; under `phi_P` they are `O` and `P`, respectively.
The quartic involution corresponds to `R -> P-R`.

**Proof.** Intersect the line `y=t(x-a)-b` with `E` and remove its known
intersection `-P`. The remaining quadratic in `x` is

\[
 x^2+(a-t^2)x+a^2+A+at^2+2bt=0.
\]

Its discriminant is R4.1, and solving it gives R4.2. This also gives the
inverse `t=(y+b)/(x-a)`, `w=2x+a-t^2` on a dense open subset. The quartic
discriminant is `256 Delta_E`, after using `b^2=a^3+Aa+B`; hence it is
separable. The birational map therefore extends between smooth projective
curves. Its degree-two function `t` has pole divisor `O+P`, identifying
the infinity points and involution.

Choose `H` over an algebraic closure with `2H=P`. The map
`psi=phi_P-H` satisfies `[2] psi=pi_P`. Its descent cocycle is translation
by `H-sigma(H)`, which equals the Kummer cocycle modulo two. This proves
the cover class, and multiplication by two has degree four. QED.

The quartic **alone**, used with its birational map `phi_P`, remains the
point-search chart in Theorem S3. Its interpretation as a labelled Selmer
cover requires R4.3. Solubility does not make a covering the identity
class in `H^1(K,E[2])`: it makes its underlying `H^1(K,E)` torsor trivial.

### Theorem R5: cubic-descent intersection of two quadrics

Suppose `b != 0`, let `F=K[theta]/(theta^3+A theta+B)`, and put

\[
 \alpha=a-\theta,\qquad U=u_0+u_1\theta+u_2\theta^2,\qquad
 \alpha U^2=q_0+q_1\theta+q_2\theta^2.
\]

The coefficients `q_i` are explicit quadratic forms over `K`, obtained by
polynomial reduction. In `P^3` with coordinates `(u0:u1:u2:s)`, set

\[
 D_P:\quad q_2=0,\qquad q_1+s^2=0.
 \tag{R5.1}
\]

This is a smooth genus-one 2-covering with affine map

\[
 \boxed{x=\frac{q_0}{s^2},\qquad
 y=\frac{b\operatorname{Norm}_{F/K}(U)}{s^3}.}
 \tag{R5.2}
\]

It represents `delta(P)=[a-theta]` and has the exact rational witness
`(1:0:0:1)`, whose image is `P`.

**Proof.** On R5.1 the identity is `alpha U^2=q0-theta s^2`. Taking norms,
using `Norm(alpha)=b^2`, gives

\[
 b^2\operatorname{Norm}(U)^2=q_0^3+Aq_0s^4+Bs^6,
\]

so R5.2 lies on `E`. Over a splitting field, the equations read
`(a-e_i)U_i^2=q0-e_i s^2`. These are the standard square-root equations
for halving on `E`, twisted by `a-e_i`: explicitly, for each root `e` of
`f(X)=X^3+AX+B`, the duplication identity is

\[
 x([2]R)-e=
 \left(\frac{(x(R)-e)^2-f'(e)}{2y(R)}\right)^2.
\]

The equation for `y` fixes the
product of the three signs, leaving four choices over a generic point.
After adjoining square roots of `a-e_i` with product `b`, they identify
the map with multiplication by two; consequently its source is smooth
of genus one and its twist class is `[alpha]`. Substitution of the stated
witness proves solubility and its image. QED.

This norm identity also supplies a cheap exact replay: reduce the product,
check both quadrics at the witness, recompute the norm, verify the image,
and verify the curve and Kummer-model transports. It needs no class group.
The same construction works without removing the quadratic term of a
separable monic cubic `f`: take `F=K[theta]/f`, and its norm identity gives
`y^2=f(x)` with the corresponding quadratic term. This is the form used
for the panel's completed-square models.
For `P=O`, use `(E,[2])`; rational two-torsion requires the usual modified
Kummer coordinates, rather than substituting a zero norm into R5.

## 4. What translate and half-lattice machinery actually discovers

### Theorem R6: translations and chart fibres preserve residual labels

Let `q:A -> A/(M+2A)=W_G` be the residual map. For `m in M` and `R in A`,

\[
 q(R+m)=q(R),\qquad q(m-R)=q(R).
 \tag{R6.1}
\]

For the pointed chart with centre `Q in M`, both birationally recovered
points `R,Q-R` therefore represent the same residual class. However the
actual 2-cover `(C_Q,pi_Q)` in R4 represents `delta(Q)`, whose residual
class is **zero**. Its birational coordinate map `phi_Q` can still recover
points in every residual class.

**Proof.** The first two identities follow by quotienting by `M+2A`, where
sign is immaterial. R4 identifies the cover class. Since `phi_Q` is an
isomorphism on the smooth projective models, it is a bijection on rational
points; restricting to an affine chart only removes its explicitly known
boundary points. QED.

More generally `P -> P+2R0` preserves the absolute 2-cover class:
on the copies of `E`, translation `R -> R+R0` satisfies

\[
 2(R+R_0)-(P+2R_0)=2R-P.
\]

Adding an arbitrary element of `M` to `P` preserves its residual cover
class. These identities justify changing representatives for smaller
models, provided every transport and subgroup is recorded exactly.

A further consequence is operationally decisive. For a soluble cover of
class `xi`, all its rational images under the **covering map** lie in one
coset of `2A`, and hence have the same residual Kummer label. Indeed,
choosing one rational point identifies its map with `R -> 2R+P`.
Consequently, searching only the covers of already known classes cannot
enlarge their residual span through their covering maps. New residual
directions require either witnesses on additional Selmer classes or point
discovery through a birational chart followed by a new Kummer calculation.
Further points of the same parity class can still add free rank; Section 1
explains why this is a separate question.

The midpoint decomposition of S3 remains the quantitative height statement:

\[
 \frac14\widehat h(2R-Q)
 =\widehat h_{/M}(\bar R)
  +\widehat h(\operatorname{pr}_M(R)-Q/2).
\]

It explains accessibility of point representatives. It supplies no
distribution of Selmer classes among the old pointed charts. In particular
the `2^r` generic midpoint labels do not enumerate `2^r` distinct residual
torsors.

## 5. A rigorous low-complexity soluble-cover filtration

Freeze a rational model for `E`, the subgroup `M`, and a presentation policy.
An admissible cover presentation must include its equations, map to that
fixed `E`, exact class identification, and every coordinate transport.
Choose explicit integer-valued encoding costs `c_model` and `c_witness`;
for example coefficient and map numerator/denominator bit lengths and
primitive weighted-projective coordinate bit lengths. Include transports
in the model cost. Finite bounds must allow only finitely many encodings.
Every soluble torsor is abstractly isomorphic to `E`, so an invariant of
the unmarked genus-one curve alone loses the Selmer label entirely. The
covering map `pi:C -> E` is an essential part of the complexity object.
One may optimize over all equivalent admissible presentations, but a
recorded minimum over a bounded dictionary is only an upper bound for that
unrestricted optimum.

For `H,T >= 0`, let `F_G(H,T)` be the set of residual classes having such a
presentation of model cost at most `H` and a rational witness of cost at
most `T`. Define its **span**, rather than assuming the set is a subspace,

\[
 W_G(H,T)=\operatorname{span}_{\mathbf F_2}F_G(H,T),\qquad
 d_G(H,T)=\dim W_G(H,T).
 \tag{R7.1}
\]

### Theorem R7: witness filtration and its exact rank implication

The spaces `W_G(H,T)` increase in both bounds, lie in `W_G`, and exhaust
`W_G` as the bounds increase, provided the presentation policy includes
R4 or R5 for every rational point. In particular

\[
 d_G(H,T)\leq\dim W_G=\Delta+\epsilon_M.
 \tag{R7.2}
\]

A collection of `q` exactly independent residual labels with exact witnesses
of the declared costs proves `Delta >= q-epsilon_M`. Under the usual
`epsilon_M=0` hypothesis it proves a jump of at least `q`.

**Proof.** A witness makes the torsor soluble, so R1 puts its label in
`W_G`. Monotonicity and the inequality follow from taking spans. Every
class of `W_G` is represented by some rational point `P`; R4 supplies a
finite presentation and rational witness for it. Since `W_G` is finite,
a basis fits within some finite pair of bounds. QED.

A bounded search produces a certified **lower** space inside `W_G(H,T)`.
Unless it exhausts the declared finite dictionary and witness box, it does
not compute that filtered space exactly. Even an exhaustive bounded miss
does not prove an individual cover insoluble, since witnesses beyond the
box remain possible. The span definition also does not assert that every
linear combination has a witness within the same bounds.

For a single budget `B`, write `d_G(B)=d_G(B,B)` and define successive
complexity thresholds

\[
 \lambda_i(E,M)=\min\{B:d_G(B)\geq i\}.
 \tag{R7.3}
\]

These make “many independent soluble classes have inexpensive witnesses”
a precise assertion. They depend on the declared presentation policy,
unlike the underlying space `W_G`. Set `lambda_i=+infinity` if the set in
R7.3 is empty, including when `i>dim W_G`.

### Proposition R8: a cheap height-to-cover bound

For an integral short model over `Q`, write
`P=(p/d^2,q/d^3)` with `d>0` and integral primitive coordinates. The change
`t=T/d`, `w=W/d^2` in R4 gives the integral quartic

\[
 W^2=T^4-6pT^2-8qT-3p^2-4Ad^4.
 \tag{R8.1}
\]

Put `X=max(1,abs(p),d^2)` and
`C_E=max(1,6,8 sqrt(1+abs(A)+abs(B)),3+4 abs(A))`.
The maximum absolute coefficient of R8.1 is at most `C_E X^2`.

**Proof.** The point equation gives
`q^2=p^3+A p d^4+B d^6`, hence
`abs(q)<=sqrt(1+abs(A)+abs(B))*X^(3/2)`. Bound the four displayed
coefficients separately. QED.

Thus low-height representatives give controlled coefficient sizes for a
soluble quartic. This direction is constructive but does not produce a
representative from an unproved soluble class. Every R4 quartic already
has tiny infinity witnesses; ignoring its coefficients and covering map
would hide the original point's complexity and make the comparison
tautological. Compare model **and** witness costs after the same exact
normalization, and record when a model was constructed from the answer.

## 6. Exact eleven-fibre soluble panel

The theorem-directed replay
[`../elliptic-curves/cas/certify_exceptional_soluble_selmer_panel.sage`](../elliptic-curves/cas/certify_exceptional_soluble_selmer_panel.sage)
constructs **110 independent basis covers**, distributed as follows.
Here `q` is a certified lower bound for `dim W_G`, not a complete residual
Selmer or Mordell--Weil calculation.

| curve | generic rank `r` | certified `q <= dim W_G` | certified rank lower bound `r+q` | primitive quadric coefficient bits, min--max |
|---:|---:|---:|---:|---:|
| 351 | 17 | 8 | 25 | 319--358 |
| 356 | 17 | 12 | 29 | 414--472 |
| 376 | 17 | 5 | 22 | 259--285 |
| 377 | 17 | 6 | 23 | 282--317 |
| 385 | 17 | 12 | 29 | 464--520 |
| 398 | 16 | 14 | 30 | 414--417 |
| 400 | 16 | 12 | 28 | 370--372 |
| 401 | 16 | 11 | 27 | 321--322 |
| 542 | 16 | 10 | 26 | 320--324 |
| 543 | 17 | 12 | 29 | 392--407 |
| 548 | 16 | 8 | 24 | 276--276 |

The exact artifact is
[`../artifacts/generated-results/elliptic-curves/exceptional_soluble_selmer_panel_v1.json`](../artifacts/generated-results/elliptic-curves/exceptional_soluble_selmer_panel_v1.json).
It retains each model, generic subgroup, candidate points, auxiliary-prime
roots and Legendre signatures, selected independent residual basis, both
quadrics, norm formula, witness `(1:0:0:1)`, and map to the point. Every
cubic is irreducible and the generic signatures have full rank `r`, so
R1.3 gives `epsilon_M=0` on every row. The auxiliary-prime cap is 2000,
all actually used primes are at most 409, and each fibre has a checkpoint.

The curve-542 mod-two certificate uses ten already found blind-search
points; together with MW16 they have binary signature rank 26. It does not
promote the different public basis's previously obtained signature rank
25 to 26. The repeated
MW16 parent labels are coordinate presentations of the previously
deduplicated fibrations, not additional observations.

R2 certifies that the known-soluble block of size `q by q` is zero on
each curve and that each of its rows pairs to zero against every actual
Selmer class. The artifact records that proof method explicitly. The full
residual Selmer dimensions, complementary Cassels--Tate entries, full
radical dimensions, and `Sha[2]` dimensions remain `UNKNOWN`; there are
currently **no certified insoluble Selmer controls** in this panel.

The coefficient-bit column is descriptive, depends on the source models
and chosen basis, and has not been minimized over equivalent covers. It
does not include a claim of small `lambda_i` or a rank predictor. Its
visible range also illustrates why a tiny stored rational witness alone
cannot certify low total complexity.

Replay uses existing points only, with zero point searches, complete
descents, or class-group computations:

```bash
sage -python elliptic-curves/cas/certify_exceptional_soluble_selmer_panel.sage --check
```

## 7. The theorem-directed comparison and its remaining gap

The existing
[`R17_KUMMER_CLASSGROUP_PRESSURE_COMPARISON_2026-09-04.md`](R17_KUMMER_CLASSGROUP_PRESSURE_COMPARISON_2026-09-04.md)
places much of the already known residual Kummer information in the
everywhere-even cubic class-group image. It says where those point classes
live. R1--R8 identify the additional data needed to explain their production.

For each of the eleven fibres the finite comparison should retain:

1. the exact generic embedding, displayed subgroup, actual mod-two Kummer
   rank, and independently proved free quotient rank;
2. explicit covers and maps for a certified basis of the known soluble
   residual subspace, with model and witness costs;
3. a separate completeness field for the full residual Selmer space;
4. exact nonzero Cassels--Tate obstructions where available, with unknown
   complementary pairing entries retained as `UNKNOWN`;
5. separate statuses for globally soluble, locally soluble but globally
   unresolved, certified globally insoluble Selmer, and locally excluded
   ambient classes.

The proposed discovery can now be stated precisely: within a fixed
normalization and comparable family, the exceptional fibres have many small
successive complexity thresholds `lambda_i(E,M)` in their **global**
residual spaces. Every such space is contained in the Cassels--Tate radical;
isotropy alone is a weaker necessary condition. Small thresholds on covers
constructed from known points are retrospective certificates. A claim that
they predict jumps requires a construction or comparison whose inputs do
not already contain those exceptional points, with genuinely unresolved
locally soluble classes included.

The point-class panel can therefore establish exact soluble-cover data and
zero Cassels--Tate rows immediately. It cannot yet compare soluble versus
insoluble Selmer classes if no latter classes have been constructed. Neither
a large class group, a large isotropic subspace, a zero pairing matrix, nor
a low-complexity old midpoint chart fills that gap.

<!-- status-consumer: EC-RATIONAL-SOLUBILITY-RESIDUAL-SELMER 431d915185bf3de9 -->

<!-- status-consumer: EC-EXCEPTIONAL-SOLUBLE-SELMER-PANEL 539bd8ec36b36c44 -->
