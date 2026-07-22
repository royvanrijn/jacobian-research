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

The values `35,69` and the doubled exponent pattern suggest linear growth
`deg V_r=34r+1`, but two orders do not prove that formula or unboundedness.
The next useful task is to derive the recursion on the leading homogeneous
valuation rather than expand the full third jet.  That would turn the current
bounded certificate into an algebraization obstruction for the entire arc.

The checker
[`verify_degree_five_contact_profile.py`](../scripts/verify_degree_five_contact_profile.py)
constructs the target-normalized arc, computes both canonical source
coefficients and inverse coefficients through order two, checks (5.2), and
verifies the first two determinant-one identities.

## 6. Current conclusion

The filtration repairs the defect of ordinary contact order:

- raw contact is infinite and carries no information;
- `kappa_m` measures how far into the automorphism ind-group one must move to
  realize that contact;
- `b_m` is an exact computable lower bound;
- unbounded growth would be a genuine algebraization obstruction;
- the degree-five family already shows a jump from degree `35` to degree `69`
  between its first two nontrivial jets.

What is not yet proved is the all-order growth law or target-gauge minimality.
Those are the two precise remaining steps before this profile itself becomes
a stable-moduli theorem.
