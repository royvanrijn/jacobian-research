# Complexity-filtered orbit contact

The unfiltered order of contact with the full polynomial automorphism orbit is
infinite for every fixed-Jacobian arc through a Keller map on `A^3`.  This
note retains the missing information by filtering the automorphism ind-group
by coordinate degree.  It also isolates a computable lower bound and evaluates
that bound for the degree-five stable-moduli family.

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
target strata, or a stronger canonical target slice.  The new degrees suggest
`24m+1`, but only the first two orders are currently established in this
torus gauge.

The checker
[`verify_degree_five_contact_profile.py`](../scripts/verify_degree_five_contact_profile.py)
constructs the target-normalized arc, computes both canonical source
coefficients and inverse coefficients through order two, checks (5.2), and
verifies the first two determinant-one identities.  It also checks the
leading terms (5.4)--(5.5), the Catalan recurrence and closed form through
eight orders, and the nonvanishing binomial factors.  The all-order conclusion
uses the degree separation in the proof above, not bounded symbolic expansion.
The torus-gauge calculation is an additional exact symbolic audit and is not
yet incorporated into that checker.

## 6. Current conclusion

The filtration repairs the defect of ordinary contact order:

- raw contact is infinite and carries no information;
- `kappa_m` measures how far into the automorphism ind-group one must move to
  realize that contact;
- `b_m` is an exact computable lower bound;
- unbounded growth would be a genuine algebraization obstruction;
- the determinant-normalized degree-five family satisfies the exact law
  `b_m=34m+1`, giving a source-algebraization obstruction.

What is not yet proved is intrinsic target coercivity or the no-escape upgrade
from a generic bounded incidence branch to a regular based family.  Parameter
numerator/denominator degree and pole order must also be controlled in the
actual lower bound.  These are the remaining steps before this filtration
itself becomes a left--right stable-moduli theorem.

More precisely, the constructible generic boundedness bridge is now complete,
as is a Noetherian theorem saying that contact in one fixed full resource box
at every order produces one exact rational equivalence.  The remaining
family-theoretic issue is no escape at the base point, and the remaining
degree-five computation is target-minimal coercivity, preferably in every
fixed stabilization dimension.  This continuation is tracked only as
`OP-LR` in [STATUS.md](../STATUS.md).
