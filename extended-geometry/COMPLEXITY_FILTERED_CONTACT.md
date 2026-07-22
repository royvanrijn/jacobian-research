# Complexity-filtered orbit contact

The unfiltered order of contact with the full polynomial automorphism orbit is
infinite for every fixed-Jacobian arc through a Keller map on `A^3`.  This
note retains the missing information by filtering the automorphism ind-group
by coordinate degree.  It also isolates a computable lower bound and evaluates
that bound for the degree-five family.  The intrinsic continuation is an
associated-graded obstruction for the left--right deformation complex, not an
exact minimization of the numerical slope in one target gauge.

Work over a characteristic-zero field `k`.  Let

\[
 F_t\in k[[t]][x,y,z]^3,
 \qquad F_0=F,
 \qquad \det D_xF_t=\det DF=c\in k^\times.           \tag{1.1}
\]

The coefficients of `t` are polynomial in `x,y,z`; no uniform source-degree
bound is assumed.

Here `Aut_1(A^3)={A:det DA=1}`, while
`SAut(A^3)=<Ga-subgroups>`.  The formal source trivializer is initially only
an `Aut_1` object; finite-jet interpolation genuinely constructs its reduced
representatives in `SAut`.

## 1. The degree filtration

For a polynomial automorphism `A` define its two-sided degree complexity by

\[
 \|A\|_{\deg}=
 \max\{\deg_x A,\deg_x A^{-1}\}.                    \tag{1.2}
\]

Including the inverse is essential: it makes the filtration symmetric and
prevents a low-degree presentation from hiding a high-degree inverse.  It is
submultiplicative in the coarse sense

\[
 \|A\circ B\|_{\deg}
 \le \|A\|_{\deg}\|B\|_{\deg}.                     \tag{1.3}
\]

Let `C_D` be the set of reduced parameter curves

\[
 A(t)\in\operatorname{SAut}_{k[t]}(\mathbb A^3),
 \qquad A(0)=\mathrm{id},
 \qquad \|A(t)\|_{\deg}\le D.                      \tag{1.4}
\]

The coefficient ring may instead be a localization of `k[t]` at a function
nonzero at zero.  The degree in (1.2) is only the degree in source variables,
not the degree or pole order in `t`.

For algebraization questions a stronger filtered profile records

\[
 \kappa_m^{\rm par}=(D_x,D_t,\nu_t),                 \tag{1.5}
\]

where `D_x` is the two-sided source-variable degree, `D_t` bounds parameter
numerator degree after choosing a reduced fraction, and `nu_t` bounds pole
order at the base point.  The source/inverse profile below is its useful first
coordinate, not a complete measure of parameter complexity.

## 2. Filtered contact and approximation complexity

Define the source approximation complexity

\[
 \boxed{
 \kappa_m^{\rm src}(F_t/F)=
 \min\left\{D:\exists A(t)\in\mathcal C_D,
 F_t\equiv F\circ A(t)\pmod {t^{m+1}}\right\}.}     \tag{2.1}
\]

Equivalently, for a fixed degree budget define

\[
 \boxed{
 c_D^{\rm src}(F_t/F)=
 \sup\{m:\kappa_m^{\rm src}(F_t/F)\le D\}.}        \tag{2.2}
\]

These are inverse encodings of the same profile.  The finite-jet
interpolation theorem in
[formal orbit triviality](FORMAL_ORBIT_TRIVIALITY.md) proves that every
`kappa_m^src` is finite.  Unlike raw contact order, the sequence can still
grow.

For left--right equivalence put

\[
 \kappa_m^{\rm LR}(F_t/F)=
 \min_{A,B}\max\{\|A\|_{\deg},\|B\|_{\deg}\},       \tag{2.3}
\]

where `A(0)=B(0)=id` and, with a fixed composition convention,

\[
 B(t)\circ F_t\circ A(t)
 \equiv F\pmod {t^{m+1}}.                           \tag{2.4}
\]

Stable versions allow identity variables before taking the minimum and record
the stabilization dimension as a second filtration parameter.

For a genuinely stable profile one should not minimize over the stabilization
dimension separately at each contact order: that would allow the number of
identity variables to grow with `m`, whereas stable equivalence uses one fixed
finite stabilization.  The appropriate object is the resource spectrum
indexed by `(m,r,D_x,D_t,nu_t)` developed in
[the ind-group constructibility note](IND_GROUP_CONSTRUCTIBILITY.md).

The numerical values depend on coordinates.  Under fixed polynomial changes
of source and target coordinates, however, (1.3) distorts them by bounded
multiplicative constants.  Thus boundedness, unboundedness, and coarse growth
type are the more intrinsic outputs.

## 3. The canonical-jet lower bound

Formal source triviality gives a unique series

\[
 \widehat\alpha(t)=
 \mathrm{id}+\sum_{r\ge1}t^rV_r,
 \qquad F_t=F\circ\widehat\alpha(t),                \tag{3.1}
\]

and a unique inverse series

\[
 \widehat\beta(t)=
 \mathrm{id}+\sum_{r\ge1}t^rW_r.                  \tag{3.2}
\]

Define

\[
 b_m(F_t/F)=
 \max_{1\le r\le m}
 \{\deg_xV_r,\deg_xW_r\}.                          \tag{3.3}
\]

### Proposition 3.1

For every `m`,

\[
 \boxed{b_m(F_t/F)\le\kappa_m^{\rm src}(F_t/F).}   \tag{3.4}
\]

Indeed, every curve matching through order `m` has, by uniqueness, the same
forward jet (3.1).  Its inverse has the jet (3.2).  Neither the curve nor its
inverse can have source degree smaller than the maximum degree of one of
their coefficients.

This turns the recursive formal trivializer into a directly computable
obstruction.  The first coefficients are

\[
 V_1=(DF)^{-1}H_1,                                  \tag{3.5}
\]

\[
 V_2=(DF)^{-1}
 \left(H_2-\frac12D^2F[V_1,V_1]\right),             \tag{3.6}
\]

while

\[
 W_1=-V_1,
 \qquad W_2=(DV_1)V_1-V_2.                          \tag{3.7}
\]

## 4. Algebraization obstruction

If a single source automorphism curve of degree at most `D` satisfies
`F_t=F\circ A(t)` exactly, then

\[
 \kappa_m^{\rm src}\le D
 \quad\hbox{and}\quad b_m\le D
 \qquad\text{for every }m.                          \tag{4.1}
\]

Consequently,

\[
 \boxed{
 \sup_m b_m=\infty
 \Longrightarrow
 \text{the canonical formal source trivializer does not algebraize with
 bounded source degree}.}                           \tag{4.2}
\]

Similarly, unbounded `kappa_m^LR` obstructs exact polynomial left--right
equivalence.  The converse is false without additional hypotheses: a
bounded-source-degree formal series may still fail to be regular or algebraic
in the parameter, and its inverse may fail to descend to the chosen reduced
base.  Thus the degree profile is a sufficient obstruction, not a complete
algebraization criterion.

There is also a distinction between (3.3) and the left--right profile.  The
canonical lower bound is attached to a chosen target gauge.  Minimizing over
all target jets may reduce it.  A stable-moduli application must therefore
either prove the target gauge intrinsic or work directly with (2.3).

Equally importantly, unbounded growth excludes one bounded-degree algebraic
automorphism **family** realizing the arc; it does not exclude left--right
equivalence of each nearby fiber separately.  Pointwise equivalences may have
degrees tending to infinity and need not assemble regularly.  A stable-moduli
theorem from this filtration therefore also needs a boundedness or
constructibility result ensuring that pointwise equivalences in an algebraic
neighbourhood can be chosen in one bounded degree stratum.

The finite-type part of that bridge is now proved in
[ind-group constructibility](IND_GROUP_CONSTRUCTIBILITY.md): bounded
two-sided automorphism degree gives a finite-type incidence scheme, pointwise
stable equivalence over an uncountable irreducible base gives one fixed degree
and stabilization dimension generically, and a generically finite base change
produces a regular equivalence on a dense open.  What constructibility does
**not** give is regularity at the chosen base point; a dominating incidence
branch may have a pole there.  This is the remaining no-escape problem and is
the reason pole order belongs in the full profile.

## 5. Degree-five test

Consider the family `F_lambda` from
[degree-five stable moduli](DEGREE_FIVE_STABLE_MODULI.md), for which

\[
 \det DF_\lambda=\frac{\lambda-1}{30}.
\]

Center at `lambda=2` and write `lambda=2+t`.  Rescale the first target
coordinate by

\[
 \frac{c_2}{c_{2+t}}=\frac1{1+t}.                   \tag{5.1}
\]

This is a based degree-one target automorphism over the local reduced base
`k[t,1/(1+t)]`, and the resulting arc `G_t` has determinant `1/30`.

The exact recursion (3.5)--(3.7) gives

\[
 \begin{array}{c|ccc}
 r&0&1&2\\ \hline
 \deg_x V_r&1&35&69\\
 \deg_x W_r&1&35&69.
 \end{array}                                        \tag{5.2}
\]

In particular,

\[
 \kappa_1^{\rm src}(G_t/F_2)\ge35,
 \qquad
 \kappa_2^{\rm src}(G_t/F_2)\ge69.                 \tag{5.3}
\]

The leading homogeneous vector terms are

\[
 (V_1)_{35}
 =360x^{20}y^8z^6(-x,0,3z),                         \tag{5.4}
\]

\[
 (V_2)_{69}
 =194400x^{40}y^{16}z^{12}(x,0,z).                  \tag{5.5}
\]

These two terms are the beginning of an all-order leading recurrence.  Put

\[
 M=(xy)^8(x^2z)^6=x^{20}y^8z^6.                    \tag{5.6}
\]

Let `W'=W+Delta` be the normalization root which reconstructs the canonical
source trivializer.  More explicitly, put

\[
 S_t=B_{2+t}C=c_{2+t}\gamma+H_{2+t}'(W),
 \qquad r_t=\frac{c_2}{c_{2+t}}=\frac1{1+t}.
\]

The target-normalized first coordinate gives the root equation

\[
 H_2(W')-S_tW'
 +r_t^2\bigl(S_tW-H_{2+t}(W)\bigr)=0.               \tag{5.7}
\]

At `t=0`, its derivative with respect to `W'` is exactly

\[
 H_2'(W)-S_0=-c_2\gamma,                            \tag{5.8}
\]

even though its highest ordinary-degree part has a double root.  If
`Delta_r` is the coefficient of `t^r`, then

\[
 \deg_x\Delta_r=34r-12,
\]

and its leading homogeneous term is

\[
 (\Delta_r)_{34r-12}
 =a_r(xy)^{8r-3}(x^2z)^{6r-2}.                      \tag{5.9}
\]

The first coefficient is `a_1=-12`.  For `r>=2`, the unique top-degree terms
in the root equation are the exact linear derivative `-c_2 gamma Delta_r`
and the quadratic Taylor term of the leading quintic.  Since the leading
quintic coefficient is `1/20`, comparison gives

\[
 a_r=15\sum_{i=1}^{r-1}a_i a_{r-i}.                 \tag{5.10}
\]

All cubic and higher Taylor terms are at least seventeen ordinary degrees
lower.  The parameter-forcing terms are also lower after `r=1`.  Thus no
omitted term can cancel (5.8).  The solution is

\[
 \boxed{
 a_r=15^{r-1}(-12)^r\operatorname{Cat}_{r-1}\ne0.} \tag{5.11}
\]

The reconstructed factor `gamma'/gamma` has leading series

\[
 1-30\sum_{r\ge1}a_r(tM)^r
 =\sqrt{1+720tM}.                                   \tag{5.12}
\]

Consequently the associated leading source transformation is

\[
 \frac{x'}x=(1+720tM)^{-1/2},
 \qquad
 \frac{z'}z=(1+720tM)^{3/2}                         \tag{5.13}
\]

in the associated ordinary-degree leading algebra.  These formulas follow
from

\[
 \gamma'=\frac{S_t-H_2'(W')}{c_2},
 \qquad x'=\frac C{\gamma'},
 \qquad u'=\frac{W'}{\gamma'},                      \tag{5.14}
\]

followed by `y'=(u'-1)/x'` and reconstruction of `z'` from `gamma'`.

Every binomial coefficient in the first series is nonzero.  The inverse
transformation has the corresponding factor
`(1-720tM')^(-1/2)`, so its coefficients are nonzero in the same degrees.
Therefore

\[
 \boxed{
 \deg_xV_r=\deg_xW_r=34r+1\quad(r\ge1),
 \qquad b_m=34m+1.}                                 \tag{5.15}
\]

This proves that the determinant-normalized degree-five arc has no exact
source trivialization of bounded polynomial degree.  In particular, its
canonical formal source trivializer does not algebraize to a polynomial
automorphism curve in this target gauge.

The same conclusion holds after every fixed identity stabilization.  Indeed,
the unique formal source trivializer of
`G_t x id_r` relative to `F_2 x id_r` is
`widehat(alpha)_t x id_r`, so its forward and inverse coefficient degrees
remain `34m+1`.  Thus stabilization introduces no new source-only gap; the
open issue is cancellation by target automorphisms.  See Proposition 5.1 of
[the constructibility note](IND_GROUP_CONSTRUCTIBILITY.md).

The chosen target gauge is not minimal.  With

\[
 \tau=1+\frac{t}{28},\qquad
 T_t(A,B,C)=(\tau^{40}A,\tau^{-23}B,\tau^{-17}C),
\]

the target automorphism has determinant one, target degree one, bounded
parameter height, and is regular at zero.  The canonical source and inverse
degrees for `T_t` composed with `G_t` are `(25,49)` through orders one and two,
rather than `(35,69)`.  Their leading monomial is
`N=x^14 y^6 z^4`, with leading vectors

\[
 (\widetilde V_1)_{25}=\frac{435}{7}N(x,0,-3z),
 \qquad
 (\widetilde V_2)_{49}=\frac{567675}{98}N^2(x,0,z).
\]

Thus target-gauge minimality cannot be proved for the determinant-normalized
gauge as stated.  The viable alternatives are intrinsic minimization over
target strata, or a stronger canonical target slice.  Redoing the root
recurrence over the invariant ring `k[xy,x^2z]` proves the exact all-order law

\[
 \deg_x\widetilde V_m=\deg_x\widetilde W_m=24m+1.
\]

The recurrence seed is `29/14`, with the same Catalan convolution as above;
see the [torus filtered-module note](TORUS_FILTERED_LR_MODULE.md).  This is
still a gauge profile, not an intrinsic target-minimality theorem.

The checker
[`verify_degree_five_contact_profile.py`](../scripts/verify_degree_five_contact_profile.py)
constructs the target-normalized arc, computes both canonical source
coefficients and inverse coefficients through order two, checks (5.2), and
verifies the first two determinant-one identities.  It also checks the
leading terms (5.4)--(5.5), the Catalan recurrence and closed form through
eight orders, and the nonvanishing binomial factors.  The all-order conclusion
uses the degree separation in the proof above, not bounded symbolic expansion.
The separate checker
[`verify_degree_five_torus_module.py`](../scripts/verify_degree_five_torus_module.py)
verifies the two-variable torus root equation, its first two exact
coefficients, the all-order `24m+1` recurrence, and survival of the candidate
class in the invariant-ring-saturated target quotient.

## 6. The associated-graded LR obstruction

The correct linear object is the polynomial LR deformation complex at `F`.
Put

\[
 \mathcal E_F=\Gamma(T_{\mathbb A^3})
       \oplus\Gamma(T_{\mathbb A^3}),\qquad
 \mathcal M_F=\Gamma(F^*T_{\mathbb A^3}),
\]

and, with the source field written first, define

\[
 L_F:\mathcal E_F\longrightarrow\mathcal M_F,
 \qquad L_F(V,W)=DF\cdot V+W\circ F.              \tag{6.1}
\]

Signs change with the convention in (2.4), but the cokernel does not.  Equip
these modules with a filtration for which `L_F` is filtered.  A Newton-polytope
filtration does this without losing the torus weights.  For ordinary degree
one may equivalently use the bifiltration recording source-field and
target-field degree separately; the elementary bounds

\[
 \deg(DF\cdot V)\le \deg F-1+\deg V,
 \qquad \deg(W\circ F)\le(\deg W)(\deg F)          \tag{6.2}
\]

compare it with the usual coordinate-degree resource.

Let `P` denote the selected filtration and introduce the Rees modules

\[
 \mathcal R_P(\mathcal E_F)=\bigoplus_d P_d\mathcal E_F\,h^d,
 \qquad
 \mathcal R_P(\mathcal M_F)=\bigoplus_d P_d\mathcal M_F\,h^d, \tag{6.3}
\]

together with `Rees_P(L_F)`.  The fiber `h=1` is the ordinary deformation
operator and the fiber `h=0` is `gr_P L_F`.  Since `det DF` is a unit, the
source summand makes `L_F` surjective without a filtration.  The Rees
cokernel therefore measures failure of filtered strictness: the obstruction
is special-fiber (equivalently, Rees-torsion) information, exactly what the
unfiltered deformation complex forgets.

### Proposition 6.1 -- reduction to target jets

Let `B_t` be any based formal target automorphism.  There is a unique based
formal source automorphism `delta_F(B_t)` such that

\[
 B_t\circ F=F\circ\delta_F(B_t).
\]

Consequently

\[
 \Lambda_F(B_t):=\delta_F(B_t)^{-1},\qquad
 B_t\circ F\circ\Lambda_F(B_t)=F.                 \tag{6.4}
\]

The map `delta_F` preserves composition, while `Lambda_F` reverses it.  If

\[
 F_t=F\circ\widehat\alpha_t,
\]

then every based formal LR trivialization
`B_t circ F_t circ A_t=F` is uniquely of the form

\[
 \boxed{A_t=\widehat\alpha_t^{-1}\circ\Lambda_F(B_t).}          \tag{6.5}
\]

Indeed, apply the formal source-triviality theorem to the deformation
`B_t circ F` over every Artin quotient `k[[t]]/(t^(m+1))`.  This gives the
unique compatible `delta_F(B_t)`.  Substitution in the LR equation gives

\[
 F\circ\delta_F(B_t)\circ\widehat\alpha_t\circ A_t=F.
\]

Formal source uniqueness applied to the constant deformation `F` forces the
source factor after `F` to be the identity, proving (6.5).  The composition
statements follow from uniqueness as well.

Thus arbitrary lower source jets are not independent variables: after a
target jet is chosen, they are forced by (6.5).  The complexity problem is the
filtered distance from `widehat alpha_t^(-1)` to the lifted target subgroup
`Lambda_F(widehat Aut_Y)`.  Notice that `Lambda_F(B_t)` is only coefficientwise
polynomial a priori; its degrees can grow without bound even when `B_t` has
small coordinate degree.  This is precisely the phenomenon OP-LR measures.

The first two coefficients make the nonlinear forcing explicit.  Write

\[
 B_t=\operatorname{id}+tW_1+t^2W_2+O(t^3),\qquad
 \delta_F(B_t)=\operatorname{id}+tD_1+t^2D_2+O(t^3).
\]

Then

\[
 D_1=(DF)^{-1}(W_1\circ F),
\qquad
 D_2=(DF)^{-1}\left(W_2\circ F
       -\frac12D^2F[D_1,D_1]\right).              \tag{6.6}
\]

If `widehat alpha_t=id+tV_1+t^2V_2+O(t^3)`, formula (6.5) gives

\[
 A_1=-(V_1+D_1),
\]

\[
 A_2=(DD_1)D_1-D_2+(DV_1)(V_1+D_1)-V_2.          \tag{6.7}
\]

Put `ell_F(W)=(DF)^(-1)(W circ F)`.  After quotienting the order-two equation
by the new linear target choice `ell_F(W_2)`, its universal class for a fixed
first target jet is represented by

\[
 \begin{split}
 \Theta_2(W_1)={}&(DD_1)D_1
 +\frac12(DF)^{-1}D^2F[D_1,D_1]\\
 &+(DV_1)(V_1+D_1)-V_2,
 \qquad D_1=\ell_F(W_1).                           \tag{6.8}
 \end{split}
\]

This is the first nontrivial arbitrary-lower-gauge audit: order one is the
linear quotient, while order two asks whether the quadratic map (6.8) can
cancel the torus residue.  Higher orders have the same form: the new `W_m`
enters linearly through `ell_F`, and every other term is a known polynomial in
the lower target jets.

There is a useful logarithmic refinement in characteristic zero.  The based
formal automorphism groups are `t`-adically prounipotent, so write

\[
 \widehat\alpha_t=\operatorname{Exp}(X(t)),\qquad
 B_t=\operatorname{Exp}(Y(t)).                     \tag{6.9}
\]

Because `delta_F` is a group homomorphism with differential `ell_F`, one has
the exact identity

\[
 \delta_F(\operatorname{Exp}Y)=\operatorname{Exp}(\ell_FY),
\qquad
 \log A=\operatorname{BCH}_{\rm map}(-X,-\ell_FY). \tag{6.10}
\]

Here the map-composition bracket is
`{P,Q}_map=(DP)Q-(DQ)P`, the negative of the usual derivation bracket.  If
`X=tX_1+t^2X_2+...` and `Y=tY_1+t^2Y_2+...`, the order-two logarithmic class,
modulo the new linear choice `ell_F(Y_2)`, is

\[
 \boxed{
 \Xi_2(Y_1)=
 \left[-X_2+\frac12\{X_1,\ell_F(Y_1)\}_{\rm map}\right].}       \tag{6.11}
\]

Thus pure target self-interactions disappear in logarithmic coordinates; they
already lie in the lifted target subgroup.  Formula (6.11) is affine linear in
the only lower target variable.  Returning to coordinate coefficients adds
the explicit term `(DZ_1)Z_1/2`, where
`Z_1=-(X_1+ell_F(Y_1))`, and recovers the quadratic representative (6.8).
The next finite audit can therefore be split into a linear bracket-module
calculation and this one known quadratic reconstruction.  Logarithmic degree
is not being substituted for automorphism coordinate degree: the
BCH-to-coordinate comparison must remain part of the Rees filtration.

For a fixed coordinate budget `D`, this gives a finite order-two incidence
problem.  Parameterize logarithmic jets `Y_1,Y_2` subject to the actual
coordinate conditions

\[
 \deg Y_1\le D,\qquad
 \deg\left(\mathord\pm Y_2+\frac12(DY_1)Y_1\right)\le D,
\]

where the two signs impose the bounds for the target map and its inverse.
Impose the analogous degree-`D` condition on the reconstructed first source
coefficient, and project

\[
 \mathord\pm\left(
 -X_2+\frac12\{X_1,\ell_F(Y_1)\}_{\rm map}
 -\ell_F(Y_2)\right)+\frac12(DZ_1)Z_1             \tag{6.11a}
\]

to source degrees above `D`; the two signs impose the bounds for the source
map and its inverse.  Vanishing of both projections is exactly the order-two
two-sided bounded-degree condition.  In the degree-five torus decomposition,
`X_1` has weight zero.  The weight-zero residue of the bracket term therefore
depends only on the weight-zero part of `Y_1`; nonzero weights enter the
quadratic reconstruction only in opposite-weight pairs.  A fixed degree
budget leaves finitely many such pairs.  Hence the next audit is a finite
weight-by-weight elimination over `R=k[v,S]`, not an unbounded
three-variable target search.

### Proposition 6.2 -- the order-two weight-zero leading sector

Let `R_24` be the five-dimensional space

\[
 R_{24}=\langle v^{12},v^9S^2,v^6S^4,v^3S^6,S^8\rangle,
 \qquad \deg v=2,\quad\deg S=3.                   \tag{6.11b}
\]

In the space of ordinary-degree-25, torus-weight-zero diagonal source fields,
the kernel of the saturated target quotient is exactly

\[
 \boxed{R_{24}\,(x,-y,-2z).}                      \tag{6.11c}
\]

Normalize the nonzero first-order torus class to

\[
 N(x,0,-3z),\qquad N=v^6S^4.
\]

Every representative of that class in this leading sector is therefore

\[
 Z_P=N(x,0,-3z)+P(x,-y,-2z),\qquad P\in R_{24}.   \tag{6.11d}
\]

For every `P`, the quadratic reconstruction

\[
 Q_P=\frac12(DZ_P)Z_P                              \tag{6.11e}
\]

has nonzero class in `\mathcal Q_R`.  In fact the third summand `R/(gamma)`
already detects it.

To verify the last assertion, put

\[
 E(P)=\left(v\partial_v-S\partial_S\right)P.
\]

Writing `Q_P=(xA,yB,zC)`, direct differentiation gives

\[
 \begin{aligned}
 2A&=(N+P)^2+2N^2+NE(P),\\
 2B&=P^2-NE(P),\\
 2C&=3N^2+12NP+4P^2-2NE(P).                       \tag{6.11f}
 \end{aligned}
\]

Its normal residue is

\[
 \rho_P=
 \left[-\frac87(A+B)v+(2A+C)S\right]_{\gamma=0}. \tag{6.11g}
\]

Writing `P` in the basis (6.11b), the coefficients of `rho_P` generate the
unit ideal in `k[p_0,\ldots,p_4]`.  Hence `rho_P` is never the zero polynomial,
even after extending the ground field.  This proves the proposition.

The calculation closes the entire homogeneous weight-zero ambiguity at
orders one and two.  In coordinate coefficients, nonzero opposite weights
appear able to pair in `(DZ_1)Z_1/2` and land in weight zero.  Proposition 6.3
explains why the correct next step is to prove Rees descent for the coset,
rather than expand those pairs individually.

The exact checker is
[`verify_order_two_weight_zero_lr_obstruction.py`](../scripts/verify_order_two_weight_zero_lr_obstruction.py).

### Proposition 6.3 -- filtered-coset BCH descent

Let `widehat G_X` and `widehat G_Y` be the `t`-adically prounipotent groups of
based formal source and target automorphisms, and let

\[
 \delta_F:\widehat G_Y\longrightarrow\widehat G_X
\]

be the target-lift homomorphism of Proposition 6.1.  Put

\[
 \widehat H_F=\delta_F(\widehat G_Y),\qquad
 \mathfrak h_F=\ell_F(\mathfrak g_Y)
       \subset\mathfrak g_X.                       \tag{6.11h}
\]

All formal LR source representatives form the single right coset

\[
 \boxed{
 \widehat{\mathcal C}_{\alpha}
   =\widehat\alpha^{-1}\widehat H_F.}              \tag{6.11i}
\]

Thus the source part of intrinsic LR complexity is the least filtered
complexity of a representative of this coset.  The full resource spectrum
uses the parametrized coset map and retains the degree and parameter resources
of `B` as well.

Give polynomial vector fields the excess-degree filtration

\[
 \operatorname{fil}(V)=\deg_xV-1,                 \tag{6.11j}
\]

or its Newton-polytope analogue.  For the map-composition pre-Lie product

\[
 P\triangleright Q=(DP)Q
\]

one has

\[
 \operatorname{fil}(P\triangleright Q)
 \le \operatorname{fil}(P)+\operatorname{fil}(Q). \tag{6.11k}
\]

Consequently `Exp`, `log`, and every BCH Lie polynomial are filtered on every
finite `t`-jet.  Equip the target Lie algebra with a compatible filtration for
which `ell_F` is filtered.

Write

\[
 \widehat\alpha=\operatorname{Exp}X,qquad
 \delta_F(B)=\operatorname{Exp}D,qquad D=\ell_F(Y).
\]

For

\[
 A=\widehat\alpha^{-1}\delta_F(B)^{-1}
   =\operatorname{Exp}(-X)\operatorname{Exp}(-D),
\]

put `Z=log A`.  Then, modulo the lifted target Lie algebra,

\[
 [Z_1]=[-X_1]\quad\text{in }\mathfrak g_X/\mathfrak h_F,      \tag{6.11l}
\]

and

\[
 \boxed{
 [Z_2]=
 \left[-X_2+\frac12\{X_1,D_1\}_{\rm map}\right].}             \tag{6.11m}
\]

In particular, no quadratic expression involving only target directions
survives on the formal homogeneous-space quotient.  Such expressions appear
when `Exp(Z)` is rewritten in coordinate coefficients, but they are the
coordinate expansion of a curve already lying in `widehat H_F`.

More generally, after passage to the formal quotient every surviving
nonlinear BCH word contains at least one `X` and at least one `D`; words made
only from target directions lie in `mathfrak h_F` because `ell_F` is a Lie
homomorphism.

#### Associated-graded form

The same conclusions hold in the associated-graded LR quotient provided the
coset action is **Rees-strict** at the filtered face in question: namely, the
special fiber of the Rees closure of

\[
 B\longmapsto\widehat\alpha^{-1}\delta_F(B)^{-1}  \tag{6.11n}
\]

must equal the orbit of the corresponding special-fiber lifted target group,
scheme-theoretically.  Equivalently, taking the initial filtered face must
commute with passing to the target-lift coset; there may be no additional
Rees torsion caused by cancellation between source and lifted-target terms.

Under this hypothesis, (6.11l)--(6.11m) specialize directly to the
associated-graded normal module.  If a torus preserves the filtration,
`X_1` has weight zero, and

\[
 D_1=\sum_pD_{1,p},
\]

then

\[
 \operatorname{pr}_0[Z_2]
 =\left[-X_2+\frac12\{X_1,D_{1,0}\}_{\rm map}\right].          \tag{6.11o}
\]

Indeed, `{X_1,D_{1,p}}_map` has weight `p`.  Opposite nonzero weights cannot
occur at logarithmic order two: their apparent coordinate-quadratic term is
target-only and disappears by Rees-strict descent.  The first possible mixed
weight-zero word involving an opposite pair occurs at order three, in a term
of the form

\[
 \{\{X_1,D_{1,p}\}_{\rm map},D_{1,-p}\}_{\rm map}.             \tag{6.11p}
\]

### Proof

Equation (6.11i) is Proposition 6.1 and the fact that a subgroup is closed
under inversion.  Inequality (6.11k) follows from
`deg((DP)Q)<=deg P+deg Q-1`; hence the `t`-adically finite BCH formulas pass to
the Rees algebra.  Expanding

\[
 \operatorname{BCH}_{\rm map}(-X,-D)
 =-X-D+\frac12\{X,D\}_{\rm map}+\cdots
\]

gives (6.11l)--(6.11m), since every `D_m` belongs to `mathfrak h_F`.
Target-only Lie words remain in `mathfrak h_F` because it is a Lie subalgebra.
Rees-strictness is precisely the condition allowing this homogeneous-space
calculation to commute with the special-fiber functor.  Finally, torus weights
add under the bracket, which proves (6.11o)--(6.11p).  QED

This proposition identifies the conceptual replacement for the proposed
opposite-weight expansion.  Proposition 6.2 verifies the required behavior in
the full weight-zero degree-25 face.  To close order two for arbitrary weights,
one should prove Rees-strictness of the degree-five target-lift coset at that
face.  If strictness holds, opposite weights disappear by (6.11o) without any
pair-by-pair coefficient calculation.  If it fails, the resulting Rees-torsion
class, rather than a larger raw expansion, is the new obstruction that must be
computed.

### Proposition 6.4 -- the target-lift second fundamental form

The quadratic part of Rees strictness is controlled by one canonical bilinear
tensor.  For target fields `Y_1,Y_2`, put

\[
 D_i=\ell_F(Y_i)=(DF)^{-1}(Y_i\circ F).
\]

Then the symmetric pre-Lie defect of the target lift is

\[
 \begin{split}
 \operatorname{II}_F(Y_1,Y_2)
 :={}&\frac12\Bigl((DD_1)D_2+(DD_2)D_1\\
 &\qquad-\ell_F\bigl((DY_1)Y_2+(DY_2)Y_1\bigr)\Bigr)\\
 ={}&-(DF)^{-1}D^2F[D_1,D_2].                     \tag{6.11q}
 \end{split}
\]

Thus `II_F` is exactly the second fundamental form of the lifted target
subgroup when source and target automorphism groups are embedded by ordinary
coordinate coefficients.  Once the linear target-lift inclusion is strict at
a filtered face, quadratic Rees strictness there is equivalent to vanishing
of the initial normal class

\[
 [\operatorname{in}\operatorname{II}_F(Y_1,Y_2)]
 \in\operatorname{gr}(\mathfrak g_X/\mathfrak h_F)             \tag{6.11r}
\]

for the allowed leading target fields.

If the setup is torus-equivariant and `Y_i` have weights `p_i`, then
`II_F(Y_1,Y_2)` has weight `p_1+p_2`.  Consequently the entire nonzero-weight
order-two strictness problem is the family of bilinear module maps

\[
 \boxed{
 \operatorname{II}_{F,p,-p}:
 \mathfrak g_{Y,p}\otimes\mathfrak g_{Y,-p}
 \longrightarrow(\mathfrak g_X/\mathfrak h_F)_0.}             \tag{6.11s}
\]

It is enough to prove that the relevant initial normal symbol of every map
(6.11s) is zero.  At a fixed Newton or ordinary-degree face these are maps of
finitely generated invariant-ring modules, so the calculation can be made on
module generators rather than on arbitrary target-jet coefficients.

### Proof

Differentiate the identity

\[
 DF\,D_2=Y_2\circ F
\]

in the source direction `D_1`.  The chain rule gives

\[
 D^2F[D_1,D_2]+DF\bigl((DD_2)D_1\bigr)
 =(DY_2\circ F)(Y_1\circ F).                      \tag{6.11t}
\]

The right side is

\[
 DF\,\ell_F((DY_2)Y_1).
\]

Multiply by `(DF)^(-1)`, interchange `1,2`, and average.  Symmetry of `D^2F`
gives (6.11q).  Formula (6.11r) is the standard second-jet criterion for the
Rees graph to have no new normal quadratic term after its strict linear
specialization.  Torus equivariance makes weights additive, proving
(6.11s).  QED

Propositions 6.3--6.4 replace the proposed opposite-weight expansion by a
specific module theorem: establish the vanishing of the initial normal symbols
of `II_(F,p,-p)`.  Higher orders are governed analogously by higher mixed BCH
fundamental forms; unlike a raw coefficient recursion, every surviving term
contains at least one deformation direction `X`.

Now expand arbitrary based source and target gauges and suppose their
coefficients through order `m-1` solve the LR equations through that order.
At order `m` the equation has the form

\[
 L_F(V_m,W_m)=-R_m(g_{<m}),                        \tag{6.12}
\]

where `R_m` is the known nonlinear forcing polynomial in the arc and in all
lower source and target coefficients `g_<m`.  Let `Z_(m-1)` be the coefficient
ind-scheme of all such lower LR jets, exhausted by its finite resource
truncations.  Then `R_m` is a universal section on `Z_(m-1)`.  On each finite
truncation, homogenize it at the relevant filtration level to obtain
`widetilde R_m`; the compatible Rees classes define

\[
 \omega_m=[\widetilde R_m]\in
 \operatorname{coker}(\operatorname{Rees}_P L_F)
       \otimes_k\mathcal O_{Z_{m-1}}.              \tag{6.13}
\]

Its relevant special-fiber leading piece

\[
 \overline\omega_m(g_{<m})
   =[\operatorname{in}_P R_m(g_{<m})]
   \in\operatorname{coker}(\operatorname{gr}_P L_F)             \tag{6.14}
\]

is the order-`m` obstruction.  This formulation handles a small but important
canonicity issue: forcing terms obtained from two lower gauges need not be
literally equal.  By Proposition 6.1 it is enough to let the lower target jet
vary; the lower source jet is then forced.  The intrinsic assertion is that
the universal class (6.13),
or the required leading face of it, is nonzero at **every** lower-jet solution
in the resource range.  It is not enough to evaluate it on the one torus
section used to discover the recurrence.

This module class is the linearized shadow of the Rees coset in Proposition
6.3.  Proving Rees-strictness makes the two descriptions agree.  Without
strictness, an apparent cancellation by lower gauges is measured by Rees
torsion of the coset action and cannot be decided from the linear cokernel
alone.

### Proposed OP-LR theorem

For the degree-five arc, choose an ordinary-degree-compatible Rees filtration.
There are constants `c>0` and `C` such that, for every `m` and every choice of
lower source and target gauges solving through order `m-1`, either a lower
coefficient already has coordinate degree at least `cm-C`, or the order-`m`
class (6.14) cannot be killed by an order-`m` correction below that degree.
Equivalently, the arc leaves every fixed coordinate-degree stratum of the
polynomial LR ind-group.

The torus computation supplies the candidate restriction of this universal
class.  In its invariant coordinates

\[
 v=xy,\qquad S=x^2z,\qquad
 N^m=v^{6m}S^{4m},                                 \tag{6.15}
\]

the source vector with coefficient `N^m` survives the
invariant-ring-saturated equivariant target quotient.  Exact equality with
the torus-gauge profile `24m+1` is stronger than OP-LR needs.

There is already a simple comparison with ordinary source degree on this
candidate quotient.  Modulo

\[
 \gamma=1-\frac87v+S
\]

one has

\[
 N^m\longmapsto
 v^{6m}\left(\frac87v-1\right)^{4m},              \tag{6.16}
\]

of `v`-degree `10m`.  A polynomial of ordinary source degree at most `D`,
written in `v,S`, restricts modulo `gamma` to a polynomial of `v`-degree at
most `D/2`, because `deg_x v=2` and substitution for `S` cannot increase
`v`-degree past that bound.  Hence every representative of this nonzero
residue has ordinary degree at least `20m` (up to the fixed shifts coming from
the vector-field basis).  The extra nonzero linear factor in the residue from
the torus-module calculation does not weaken the estimate.  Thus the
invariant-ring-to-ordinary-degree translation already gives a positive linear
bound **once** (6.15) is identified with the universal class (6.13).  The real
remaining algebraic step is survival after arbitrary lower-order gauges.

Any positive linear bound is sufficient.  By the bounded-box algebraization
theorem, unbounded intrinsic coordinate degree forces escape from every fixed
full resource box and rules out one exact rational LR equivalence in such a
box.  Neither the coefficient `24` nor an exact target-minimal slope is needed.

### A geometric route to no escape

The pole-at-the-base issue is logically separate from the associated-graded
degree obstruction.  A promising valuative formulation uses the proper marked
finite-cover object supplied by the
[Hurwitz--LL compactification](HURWITZ_LL_COMPACTIFICATION.md) and decorated
normalization.  A meromorphic LR equivalence over a punctured DVR should
induce an isomorphism of the two proper marked objects over the generic point.
Coarse affine-mark descent now supplies the unique specialization of each mark
over every DVR limit.  What remains for no escape is to prove that the generic
LR-induced isomorphism belongs to a finite proper `Isom` closure compatible
with those marks.  The valuative criterion would then extend it over the DVR,
and rigidity would identify its special fiber.  This could exclude pole escape
without coefficientwise pole estimates.  The required functorial `Isom`
statement is still a proposed no-escape lemma; uniqueness of the individual
marked specializations alone does not prove it.

## 7. Current conclusion

The filtration repairs the defect of ordinary contact order:

- raw contact is infinite and carries no information;
- `kappa_m` measures how far into the automorphism ind-group one must move to
  realize that contact;
- `b_m` is an exact computable lower bound;
- unbounded growth would be a genuine algebraization obstruction;
- the determinant-normalized degree-five family satisfies the exact law
  `b_m=34m+1`, giving a source-algebraization obstruction.

What is not yet proved is Rees-strictness of the target-lift coset and
subsequent survival of its class after every lower-order target choice, or the
no-escape upgrade from a generic bounded incidence branch to a regular based
family.  Parameter numerator/denominator degree and pole order remain part of
the full resource spectrum, but a linear intrinsic coordinate-degree
obstruction would already force escape from every fixed box.

More precisely, the constructible generic boundedness bridge is now complete,
as is a Noetherian theorem saying that contact in one fixed full resource box
at every order produces one exact rational equivalence.  The remaining
family-theoretic issue is no escape at the base point, and the remaining
degree-five computation is the universal associated-graded nonvanishing in
every fixed stabilization dimension.  Stable moduli do not depend on this
continuation: decorated normalization together with the affine sheet already
supplies them.  `OP-LR`, tracked in [STATUS.md](../STATUS.md), is therefore an
independent theorem about filtered polynomial LR orbits.
