# The RR double plane and an absolutely simple genus-two Jacobian

## Result and relevance to the requested descent comparison

**Verified application and new deduction.** The historical generic RR net
gives an explicit double-plane model of the determinant1092 surface,

\[
W^2=c\,B(X,Y,Z),\qquad \deg B=6.
\]

The polynomial has25 nonzero terms and exactly one geometric singularity:
an ordinary node `N=[0:1:0]`. The two previously conflated curve families
have different plane descriptions:

| Plane line | Curve on the resolved double plane |
|---|---|
| `X=tau Z`, through `N` | The original elliptic fibre `E_tau`, with the contracted bisection removed |
| `Y=uZ+vX`, not through `N` | The genus-two RR member `C_(u,v)`, when smooth |

**Verified application.** The fixed, source-selected RR member `B`, meaning
`u=v=0`, has absolutely simple Jacobian. Its good reduction at17 has

\[
\#C_B(\mathbf F_{17})=30,\quad
\#C_B(\mathbf F_{289})=280,\quad
\boxed{P_{17}(x)=x^4+12x^3+67x^2+204x+289.}
\]

**New deduction.** Consequently neither this fixed member nor the geometric
generic RR member admits a nonconstant map to *any* elliptic curve, of any
degree. This excludes a generic elliptic-quotient bridge from these
genus-two Jacobians to the marked fibres, not merely the previously rejected
degree-two branch-cover construction.

**Unresolved.** Special RR loci with elliptic factors have not been
fully classified. The degree-two locus now has an explicit necessary
equation, and the entire rational first-witness pencil is excluded from it
below. Other loci could still support an appropriate construction,
but it must be specified and its elliptic quotient identified. Absolute
simplicity is not a rational-solubility obstruction and does not disprove
an empirical incidence correlation. A subsequent
[uniform inherited-point identity](DET1092_UNIVERSAL_RR_DESCENT_PREFLIGHT_2026-09-08.md#decisive-correction-all-smooth-rr-members-are-already-rationally-pointed)
does, however, prove that raw RR-curve solubility is constant on the smooth
rational locus. The proposed *residual* global class characterization remains open.

## Explicit double-plane equation and maps

**Verified application.** Start with the generic seventeen sections and the
historical norm10 centre `C=P_w`. As in the
[universal RR equation](DET1092_UNIVERSAL_RR_DESCENT_PREFLIGHT_2026-09-08.md),
put `r=u+vT` and write

\[
q(T;u,v)=Q(T,u+vT),\qquad
Q(T,r)=\sum_{i,j}q_{ij}T^i r^j.
\]

Exact reconstruction from `A,B` verifies `i+j<=6`, `j<=4`, and no common
polynomial factor in `T` among the coefficients in `r`. Therefore

\[
B(X,Y,Z)=\sum_{i,j}q_{ij}X^iY^jZ^{6-i-j}
\]

is the explicit homogeneous sextic. Its scalar `c` is retained, not removed
as an irrelevant geometric factor. Every coefficient is in the
[construction certificate](../../artifacts/generated-results/elliptic-curves/det1092_rr_plane_jacobian_gate_v2.json).
The independent replay checks the entire polynomial identity

\[
cB(T,u+vT,1)=c_{\rm previous}q_{\rm previous}(T;u,v),
\]

and reconstructs the parent equation identity with denominators cleared.
It does not infer the equation from a finite set of line restrictions.

On the affine surface, let `A(T,x,y)` and `B_RR(T,x,y)` denote the two saved
RR line functions, after cancelling the common trace component. The plane
map on its displayed open set is

\[
r=-B_{\rm RR}(T,x,y)/A(T,x,y),\qquad
[X:Y:Z]=[T:r:1].
\]

For the reverse map, set `f_i=B_i+rA_i` and
`m=-f_1/f_2+a_1/2`. In the retained short elliptic model use

\[
\omega=h^3W/f_2^2,\qquad
x_s=(m^2-c_x+\omega)/2,\qquad
y_s=m(x_s-c_x)-c_y,
\]

followed by the inverse short-model transport. These maps are defined off
the displayed denominators and birational identifications are understood
on the smooth models. Generic invertibility of the slope-to-`r` Mobius map
is already certified. No claim of everywhere-defined affine formulas is
made.

**Established literature, with exact application here.** The divisor
`D=C_min+F` has square2 and `D.C_min=0`. The double-plane interpretation of a
rational bisection plus an elliptic fibre, with the bisection contracted to
a node and the original fibration recovered from the lines through that
node, is classical; see
[Rank jumps and multisections of elliptic fibrations on K3 surfaces, Proposition3.4](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/rank-jumps-and-multisections-of-elliptic-fibrations-on-k3-surfaces/602807FB00055D106E3CEA5418DE08F7).
The coefficient identities above verify this application without asserting
that the double-plane construction itself is new.

## Node and the two kinds of line sections

**Verified application.** At `N=[0:1:0]`, put `Y=1`. Every monomial has
`X,Z`-degree at least2. The degree-two tangent cone has nonzero discriminant,
so `N` is an ordinary node. Its three coefficients are retained in the
certificate. If `q_4(T)` is the coefficient of `r^4` in `Q`, then the replay
also proves

\[
cq_4(T)=\kappa^2q_{\min}(T),\qquad\kappa\in\mathbf Q^\times,
\]

where `q_min` is the previously certified minimal rational bisection's
quadratic branch polynomial. This identifies the exceptional conic over the
node with that existing construction; it is not another independent point
source.

**New deduction using the completed halving theorem.** Over `Qbar(T)`, the
quartic in `r` is Mobius-equivalent to the connected halving quartic. With
the checked absence of vertical polynomial content, the plane sextic is
geometrically integral. Its normalization is the already certified
genus-nine halving curve. An integral plane sextic has arithmetic genus10.
Thus the node's delta-invariant1 exhausts the total singularity defect;
there are no additional geometric singularities. No plane singular-locus
enumeration or Groebner-basis calculation was needed.

Substitution `X=tau Z` gives a factor `Z^2` and a homogeneous quartic in
`Y,Z`. Removing the square factor gives the pointed elliptic model. In
contrast, the line `Y=uZ+vX` misses `N`; its transverse branch divisor has
degree six and gives genus two. This explains the parameter mismatch in the
earlier descent proposal geometrically, not just notationally.

## Exact Jacobian certificate for the fixed member

**Verified application.** No genus-two member was fitted to302. The
constructor uses exactly `u=v=0`, previously fixed by generic RR row
algebra. A protocol recorded before computation permits at most the same
eighteen fixed primes through197 and stops at the first certificate. The
first prime,17, suffices. The reduced curve is

\[
y^2=3T^6+7T^5+T^4+7T^3+6T^2+6\pmod {17}.
\]

Its sextic is squarefree. The constructor counts points by iterating
`F17` and `F289`. The independent replay represents `F289` as integer pairs
`a+b*sqrt(3)`, evaluates the sextic by Horner arithmetic, and uses the
quadratic character of the norm `a^2-3b^2`. Both include the points at
infinity. They agree on30 and280. The point counts determine the displayed
degree-four Frobenius polynomial.

**Established literature applied exactly.** The polynomial is irreducible
over `Q` and ordinary since `67` is prime to17. In the notation
`x^4+ax^3+bx^2+aqx+q^2`, the four exceptional equalities of
[Howe--Zhu, Theorem6](https://arxiv.org/abs/math/0002205) all fail:

\[
a=12\ne0,\qquad a^2=144\notin\{q+b,2b,3b-3q\}
=\{84,134,150\}.
\]

The reduction is therefore absolutely simple. Independently of this
coefficient test, the replay raises a rational companion matrix to powers
`2,3,4,6` and verifies that all four extension Frobenius polynomials are
irreducible. Those are precisely the possible first splitting degrees in
the cited ordinary-surface theorem. Their coefficients are retained in the
[replay certificate](../../artifacts/generated-results/elliptic-curves/det1092_rr_plane_jacobian_gate_replay_v2.json).

**Established specialization principle and new deduction.** Endomorphisms
specialize injectively under good reduction; see
[Conrad, Semistable reduction, Proposition6.4](https://math.stanford.edu/~conrad/DarmonCM/2011Notes/SemistableReduction.pdf).
A nontrivial idempotent expressing geometric splitting in characteristic
zero would specialize to one in the absolutely simple reduction, which is
impossible. Thus `Jac(C_B)` is absolutely simple over `Qbar`. The same
specialization argument in the smooth RR family shows that its geometric
generic Jacobian is absolutely simple. One can apply the argument along
successive DVR specializations `v=0`, then `u=0`; smoothness at the final
member ensures the required good models.

A nonconstant map from a smooth genus-two curve to an elliptic curve
induces a nonzero homomorphism from its Jacobian onto that elliptic curve.
An absolutely simple abelian surface has no such quotient. Hence no such
map exists for the fixed `B` member or for the geometric generic RR member,
even after algebraic field extension. This also rules out a nonzero
Jacobian homomorphism induced by an algebraic correspondence to an elliptic
curve for those members. It does not classify special split-Jacobian loci.

## What this says about descent, and what remains unknown

### A universal degree-two quotient gate and the entire first-witness pencil

**Established covariant construction and new application.** There is now
an exact equation-side gate for special RR members with a degree-two
elliptic quotient. For the binary sextic `F(X,Z)` of any RR member, define
unnormalized transvectants by

\[
(f,g)_k=\sum_{s=0}^k(-1)^s\binom{k}{s}
\frac{\partial^k f}{\partial X^{k-s}\partial Z^s}
\frac{\partial^k g}{\partial X^s\partial Z^{k-s}}.
\]

Put

\[
i=(F,F)_4,\quad Y_1=(F,i)_4,\quad
Y_2=(i,Y_1)_2,\quad Y_3=(i,Y_2)_2,
\qquad
\mathcal I_{15}(F)=\det\begin{pmatrix}
[Z^2]Y_1 &[XZ]Y_1 &[X^2]Y_1\\
[Z^2]Y_2 &[XZ]Y_2 &[X^2]Y_2\\
[Z^2]Y_3 &[XZ]Y_3 &[X^2]Y_3
\end{pmatrix}.
\]

These are the quadratic covariants of coefficient degrees3,5,7 in
[Cardona--Quer, Table1](https://arxiv.org/abs/math/0207015).
The different normalization of transvectants only changes nonzero constant
factors. Our determinant has coefficient degree15 and `GL2` determinant
weight45. In particular its vanishing is independent of a rational
coordinate change or a nonzero scalar twist of `F`.

**New deduction, with a direct proof of the necessary condition.** A
degree-two map from a smooth genus-two curve to an elliptic curve supplies
a nonhyperelliptic involution. Over the algebraic closure it can be written
as `X -> -X` on the hyperelliptic base, with an even sextic equation. For
an even sextic all three `Y_i` are even quadratics, so their `XZ` columns
vanish. Covariance therefore proves

\[
\boxed{\text{degree-two elliptic quotient}\quad\Longrightarrow\quad
\mathcal I_{15}(F)=0.}
\]

Only this necessary implication is used. A zero is not claimed to construct
a quotient over `Q`, identify it with `E_tau`, or certify an independent
point.

Applied to the universal binary RR restriction
`F_(u,v)(X,Z)=B(X,uZ+vX,Z)`, this determinant is an explicit source-only
polynomial recipe for the possible bielliptic locus. Projectively it is a
section of `O(45)` on the dual plane: changing a two-vector basis of a line
multiplies the invariant by its determinant to the45th power. Equivalently,
the `SL2`-invariant polynomial in that frame is a degree45 polynomial in
its Plucker coordinates. Thus the affine expression has total degree at
most45. The nonzero degree45 restriction below also proves that the
universal affine degree is exactly45. This is a necessary-condition locus,
not a complete rational parametrization of its points.

**Retrospective verified application.** Use exactly the immutable pencil
`u=u0` determined by the first historical witness, not a newly fitted
pencil or a prospective selection rule. The exact primitive integer
polynomial

\[
A(v)=\operatorname{primitive}\bigl(\mathcal I_{15}(F_{u_0,v})\bigr)
\]

has degree45. Its complete coefficients and all three intermediate
quadratic covariants are retained in the
[certificate](../../artifacts/generated-results/elliptic-curves/det1092_first_pencil_elliptic_involution_v1.json).
At the thirteenth prime of the unchanged eighteen-prime list,127, its
degree remains45 and

\[
\gcd\bigl(A(v)\bmod127,\;v^{127}-v\bigr)=1.
\]

All earlier trials are retained, including the degree drops at67,89,101.
No prime was added after seeing an outcome. Since `A` has integral
coefficients and leading coefficient a127-adic unit, any rational root
would be127-integral and reduce to a root. Therefore `A` has no rational
root.

**New deduction.** Every rational affine member of this entire pencil
has no degree-two elliptic quotient, even over the algebraic closure.
These members are all smooth of genus two by the previous
[genus-drop obstruction](DET1092_FIRST_WITNESS_PENCIL_GENUS_GATE_2026-09-08.md).
The projective member at `v=infinity` is the known reducible
`C_min+F_0`, not a genus-two construction. Thus the excluded class includes
every genus-two RR member in this fixed centre that passes through the
literal first witness and could use a degree-two elliptic quotient.

**Verified distinction.** These genus-two curves still have rational
points: the constant coefficient `q(0,v)` is independent of `v`, and its
product with the retained squareclass scale is a nonzero rational square.
This is the common first-witness point, so every rational member is
everywhere locally soluble as well. The invariant obstruction concerns
an *elliptic quotient*, not solubility of the curve. Selecting this pencil
using its known point would make rational nonemptiness tautological and is
not an admissible blind selector for the controls.

The [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_first_pencil_elliptic_involution_replay_v1.json)
uses Sage10.9's documented `ubs` implementation of the normalized Mestre
covariants, whose source hash is retained. It reproduces all covariant
coefficients after the explicit normalization factors, verifies the entire
primitive polynomial, and checks all residue root counts by Frobenius gcds.
It also compares the pencil with the universal equation, checks the saved
chart-to-`u0` map, and verifies symbolic even-sextic vanishing and four
`GL2` transformation regressions. No constructor function is imported.

**Scope remaining open.** Elliptic quotients of degree at least three,
other rational points of the universal necessary-condition locus, and
other RR pencils remain unclassified. This calculation does not identify
a Selmer class or establish an arithmetic rank-incidence predictor.

**Superseded solubility status.** The finite-field computation proved local
solubility at17. The subsequent uniform-point certificate now proves
rational and everywhere-local solubility, and nonemptiness of the curve's
fake2-Selmer set, for this member and every smooth rational RR member.
The full Jacobian2-Selmer group and Cassels--Tate data remain unknown.
These conclusions come from the inherited point, not absolute simplicity.

**New deduction.** A generic genus-two RR Jacobian cannot be connected to
the marked elliptic fibre by an elliptic quotient, so a genus-two descent
on arbitrary RR members cannot be justified as that quotient's descent.
An empirical arithmetic comparison does not require an elliptic quotient.
A constructive rank-gain claim does require a specified map to the marked
fibre and a non-inherited point criterion. The degree-two
split-locus route cannot pass through the first witness in its fixed RR
centre, by the whole-pencil obstruction above. No claim is made that
all predictive arithmetic mechanisms are impossible.

**Computation record.** The initialv1 attempt stopped before any finite-field
count because Sage's exact polynomial division returned coefficients in a
fraction field. Its source and frozen protocol remain retained. Version2
adds the exact coercion back to the polynomial ring; it does not change the
member, primes or stopping rule. The v2 constructor and independent replay
each finished in under one second under a25-second process cap. No Selmer,
class-group, rational point search, control sweep or V3 mutation was run.
No background process is active for this calculation.

The invariant supplement computed one exact polynomial on one immutable
diagnostic pencil and13 reductions within an18-prime cap. Its constructor
and independent replay each finished in under one second under a25-second
process cap. The full803KB certificate retains the exact intermediate
covariants; the independent replay certificate is compact. No rational
point, Selmer, class-group or new control search was started.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_rr_plane_jacobian_gate.sage
sage -python research/elliptic-curves/cas/verify_det1092_first_pencil_elliptic_involution.sage
```
