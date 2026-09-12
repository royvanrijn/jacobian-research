# Torus reduction of the filtered left--right problem

This note carries out the two reductions suggested by the residual source
torus for the degree-five arc.  It proves the all-order law `24m+1` in the
specific target torus gauge and computes an `R=k[v,S]`-saturated target-image
quotient in the relevant torus-weight sector.  The quotient residue gives a
linear comparison with ordinary degree.  It does **not** yet prove the
intrinsic filtered left--right lower bound: arbitrary lower source and target
choices still have to be incorporated into one filtered induction.

## 1. Two invariant variables

Put

\[
 v=xy,\qquad S=x^2z,\qquad
 u=1+v,\qquad \gamma=1-\frac87v+S,\qquad W=u\gamma .
\]

The source torus has weights `(1,-1,-2)`, so its invariant ring is

\[
 R=k[v,S]=k[u,\gamma].                              \tag{1.1}
\]

Let

\[
 \tau=1+\frac t{28},\qquad r=(1+t)^{-1}.
\]

After the fixed-Jacobian normalization is followed by

\[
 (A,B,C)\longmapsto
 (\tau^{40}A,\tau^{-23}B,\tau^{-17}C),              \tag{1.2}
\]

write

\[
 \Sigma_t=c_{2+t}\gamma+H'_{2+t}(W).
\]

The normalization root `W'=W+Delta` for the source trivializer is determined
by the single equation over `R[[t]]`

\[
 H_2(W')-\tau^{-40}\Sigma_tW'
 +r^2\tau^6\bigl(\Sigma_tW-H_{2+t}(W)\bigr)=0.      \tag{1.3}
\]

Thus the three-variable reconstruction is not needed for the leading
recurrence.

## 2. The all-order recurrence

Give the invariant variables their ordinary source degrees

\[
 \deg u=2,\qquad \deg\gamma=3.                      \tag{2.1}
\]

At `t=0`, the derivative of (1.3) with respect to `W'` is exactly

\[
 H_2'(W)-\Sigma_0=-\frac{\gamma}{30}.               \tag{2.2}
\]

The first-order forcing at `Delta=0` is

\[
 [t]E(W,0,t)=
 \frac{\gamma^2u(29\gamma u^2-24u-5)}{420}.         \tag{2.3}
\]

Its top part is `29u^3 gamma^3/420`.  Consequently, if

\[
 \Delta_m=a_m u^{6m-3}\gamma^{4m-2}
          +\text{terms of lower weight},            \tag{2.4}
\]

then

\[
 a_1=\frac{29}{14}.                                 \tag{2.5}
\]

For `m>=2`, the unique terms of weight `24m-9` in (1.3) are

\[
 -\frac{\gamma}{30}\Delta_m
 +\frac{W^3}{2}\sum_{i=1}^{m-1}\Delta_i\Delta_{m-i}.
\]

Indeed, a positive-`t` coefficient in the linear term loses at least seven
weights, every cubic or higher Taylor term loses at least seventeen, and the
parameter forcing has weight at most twenty-five.  It follows that

\[
 \boxed{
 a_m=15\sum_{i=1}^{m-1}a_i a_{m-i}},                \tag{2.6}
\]

and hence

\[
 \boxed{
 a_m=15^{m-1}\left(\frac{29}{14}\right)^m
       \operatorname{Cat}_{m-1}\ne0}.               \tag{2.7}
\]

Put

\[
 N=v^6S^4.
\]

The top part of `u^(6m) gamma^(4m)` is `N^m`.  Reconstruction from the root
therefore gives

\[
 \frac{\gamma'}\gamma
 \sim 1-30\sum_{m\ge1}a_m(tN)^m
 =\sqrt{1-\frac{870}{7}tN}.                         \tag{2.8}
\]

The scalar powers of `tau` do not affect the top source degree.  Thus

\[
 \frac{x'}x\sim
 \left(1-\frac{870}{7}tN\right)^{-1/2},\qquad
 \frac{z'}z\sim
 \left(1-\frac{870}{7}tN\right)^{3/2},              \tag{2.9}
\]

while the top term of the `y` coefficient vanishes.  In particular,

\[
 (V_m)_{24m+1}=N^m(\alpha_mx,0,\beta_mz),           \tag{2.10}
\]

where

\[
 \alpha_m={-1/2\choose m}\left(-\frac{870}{7}\right)^m,
 \qquad
 \beta_m={3/2\choose m}\left(-\frac{870}{7}\right)^m. \tag{2.11}
\]

Both coefficients are nonzero for every `m`.  For `m=1,2`, (2.10) is exactly

\[
 \frac{435}{7}N(x,0,-3z),\qquad
 \frac{567675}{98}N^2(x,0,z).
\]

At the same leading level,

\[
 \frac{v'}v\sim(1-q)^{-1/2},\qquad
 \frac{S'}S\sim(1-q)^{1/2},\qquad
 N'\sim\frac{N}{1-q},\qquad q=\frac{870}{7}tN.
\]

Writing `q'=(870/7)tN'` gives `q=q'/(1+q')`.  The inverse `x` and `z`
factors are therefore `(1+q')^(-1/2)` and `(1+q')^(3/2)`, whose coefficients
are again all nonzero.  Hence the forward and inverse source profiles in this
fixed target gauge satisfy

\[
 \boxed{\deg V_m=\deg W_m=24m+1\quad(m\ge1).}        \tag{2.12}
\]

## 3. Target-image module over the invariant ring

For `F_2`, write

\[
 F_2=(x^{-2}a(u,\gamma),x^{-1}b(u,\gamma),x\gamma),
\]

where

\[
 a=u(6\gamma^3u^4-21\gamma^2u^3+17\gamma u^2-3u+1),
\]

\[
 b=\frac{15\gamma^3u^4-56\gamma^2u^3+51\gamma u^2-12u+2}{60}. \tag{3.1}
\]

The target weights are `(-2,-1,1)`.  Over the target invariant ring, the
semi-invariant coefficient modules in those three weights have generators

\[
 (A,B^2),\qquad(B,AC),\qquad(C),                    \tag{3.2}
\]

respectively.  Pull back and enlarge coefficients to `R`.  This enlargement
can only make the target image larger, so survival in the enlarged quotient
implies survival in the actual target quotient.  The resulting diagonal
quotient is

\[
 \boxed{
 \mathcal Q_R=
 R/(a,b^2)\oplus R/(b,a\gamma)\oplus R/(\gamma).}   \tag{3.3}
\]

Equivalently, identify source fields with map perturbations using the
logarithmic differential matrix

\[
 J=
 \begin{pmatrix}
 -2a&a_u&a_\gamma\\
 -b&b_u&b_\gamma\\
 \gamma&0&1
 \end{pmatrix},
 \qquad \det J=\frac1{30}.                          \tag{3.4}
\]

Since its determinant is a unit, `J^(-1)` transports the target-field image
back to the source module without denominators.  Thus no three-variable
Gröbner basis occurs.

## 4. Survival of the candidate class

Apply `DF_2` to the leading source vector (2.10), or equivalently use the
source quotient transported by `J^(-1)`.  Its normalized third component,
reduced modulo the last target summand `(gamma)`, is

\[
 N^m\left[-(2\alpha_m+\beta_m)
       +\frac87(\alpha_m+\beta_m)v\right]           \tag{4.1}
\]

in

\[
 R/(\gamma)=k[v],\qquad S=\frac87v-1.              \tag{4.2}
\]

Here

\[
 N^m=v^{6m}\left(\frac87v-1\right)^{4m}\ne0.
\]

The linear factor in (4.1) could vanish identically only if
`alpha_m+beta_m=0` and `2alpha_m+beta_m=0`, which would force both coefficients
to be zero.  Equation (2.11) excludes this.  Therefore

\[
 \boxed{[N^m(\alpha_mx,0,\beta_mz)]\ne0
        \text{ in }\mathcal Q_R\quad(m\ge1).}       \tag{4.3}
\]

This is stronger than a bounded target-degree test because (3.3) used the
`R`-saturation of the whole equivariant target image.

If identity variables are assigned weight zero, they only tensor the modules
above with a weight-zero polynomial ring and add block summands.  Restriction
to the zero section preserves the nonzero residue (4.1), so those summands do
not kill this particular class.

## 5. Associated-graded interpretation and remaining step

The exact profile `24m+1` is a property of the chosen torus gauge.  The part
relevant to intrinsic LR coercivity is instead the nonzero class (4.3), viewed
as the restriction of an order-`m` class in the associated-graded cokernel of

\[
 (V,W)\longmapsto DF_2\cdot V+W\circ F_2.          \tag{5.1}
\]

There is a direct ordinary-degree estimate on the detected residue.  Give
`v,S` their ordinary source degrees `2,3`.  If a polynomial of weighted degree
at most `D` is reduced modulo `gamma`, then its image in `k[v]` has `v`-degree
at most `D/2`: a monomial `v^aS^b` maps to
`v^a(8v/7-1)^b`, of `v`-degree `a+b`, while

\[
 2(a+b)\le 2a+3b.                                  \tag{5.2}
\]

The factor `N^m` in (4.1) has residue

\[
 v^{6m}\left(\frac87v-1\right)^{4m},              \tag{5.3}
\]

of `v`-degree `10m`.  Multiplication by the nonzero linear factor in (4.1)
leaves the degree at least `10m`.  Therefore any representative detected by
this quotient has ordinary source degree at least `20m`, up to the fixed
degree shift introduced by the logarithmic vector-field basis.  Thus a
positive linear ordinary-degree bound is already present in the torus-weight
quotient; recovering the exact slope `24` is unnecessary.

The remaining filtered step is canonicity under arbitrary lower gauges.  The
formal target-lift reduction shows that every lower source jet is uniquely
forced after the lower target jet is chosen.  Thus if `Z_(m-1)` denotes the
target-jet scheme through order `m-1`, the nonlinear order-`m` forcing defines
a universal Rees-cokernel class over `Z_(m-1)`.  One must prove that its
relevant special-fiber leading face restricts to (4.3), or to another class
with a positive linear degree bound, at every point of `Z_(m-1)` in the
proposed resource range.  At order two this is the explicit quadratic class
`Theta_2(W_1)` in coefficient coordinates, or the affine-linear bracket class
`Xi_2(Y_1)` plus its known quadratic reconstruction in logarithmic
coordinates.  The present torus calculation verifies only its
diagonal-torus section.  See
[complexity-filtered contact](COMPLEXITY_FILTERED_CONTACT.md#6-the-associated-graded-lr-obstruction)
for the Rees-module formulation.

The first extension beyond the diagonal-torus section is now complete.  In
the full degree-25 homogeneous weight-zero source sector, the saturated target
kernel is exactly

\[
 R_{24}(x,-y,-2z),
 \qquad
 R_{24}=\langle v^{12},v^9S^2,v^6S^4,v^3S^6,S^8\rangle.
\]

Adding any such kernel field to `N(x,0,-3z)` leaves a quadratic order-two
reconstruction with nonzero class; its third `R/(gamma)` residue alone cannot
vanish.  Thus arbitrary weight-zero lower target choices do not kill the
order-two leading obstruction.  The filtered-coset BCH descent lemma shows
that opposite nonzero weights would disappear at logarithmic order two under
Rees-strict descent.  Its quadratic obstruction is the second fundamental form
`II_F(Y_1,Y_2)=-(DF)^(-1)D^2F[ell_F(Y_1),ell_F(Y_2)]`, so the nonzero-weight
problem is a family of invariant-module maps `II_(F,p,-p)`, not an arbitrary
coefficient expansion.  The minimal pair is now known to be nonzero: for
constant target fields `(partial_B,partial_C)` of weights `(1,-1)`, the third
saturated normal summand has leading symbol `-146880u^5/7`.  Hence quadratic
Rees strictness fails and this symbol replaces the former vanishing gap with a
canonical invariant; linear target-lift Rees strictness remains separate.  The
exact weight-zero calculation is Proposition 6.2 of the contact note and is checked by
[`verify_order_two_weight_zero_lr_obstruction.py`](../scripts/verify_order_two_weight_zero_lr_obstruction.py).

The exact checker is
[`verify_degree_five_torus_module.py`](../scripts/verify_degree_five_torus_module.py).
