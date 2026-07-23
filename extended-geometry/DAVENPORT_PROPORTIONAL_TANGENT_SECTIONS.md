# Proportional tangent sections for the Davenport pair

The fixed-mark tangent curve does not admit a nonconstant map from an
affine line, but allowing both tangent points to move changes the answer.
There are four exact proportional polynomial sections.  One is fixed by
the quadratic conjugation and gives simultaneous marks for the two
Davenport covers over the same quadratic base change

\[
T=-s^2.
\]

This removes the finite tangent-mark cover from the relative Sunada--Keller
construction.  It does not yet produce an absolute affine-space Keller
pair: determinant normalization and the distinguished endpoint introduce
two explicit divisors on the \(s\)-line.

Work over

\[
K=\mathbb Q(a),\qquad a^2+a+2=0.
\]

## 1. The proportional square-discriminant locus

For the first Davenport polynomial \(g_T\), put

\[
\mathcal T(T,q,r)=
\frac{g_T(r)-g_T(q)-g_T'(q)(r-q)}{(r-q)^2}.
\]

Write \(r=q+d\) and impose \(d=kq\).  The equation
\(\mathcal T(T,q,q+kq)=0\) is quadratic in \(T\).  After removing its
forced \(q^4\) factor, let \(D_k(q)\) be its discriminant in \(T\).
An exact calculation gives

\[
\begin{aligned}
\operatorname{disc}_q D_k
=-\frac{16}{49}F_1(k)F_2(k),
\end{aligned}
\]

where

\[
F_1=k^5+7k^4+21k^3+35k^2+35k+21
\]

is irreducible over \(K\), while

\[
\begin{aligned}
F_2={}&(11a-7)
 \left(k+1-\frac{3a}{2}\right)
 \left(k+\frac{48}{23}-\frac{9a}{23}\right)\\
&\mathrel{}\cdot
 \left(k+2+\frac a2\right)(k+4+a)^2.
\end{aligned}
\]

For all four roots of \(F_2\), the square multiplier in \(D_k(q)\) is
itself a square in \(K\).  Exactly one root of the quadratic equation in
\(T\) is polynomial.  The four sections are

\[
\begin{array}{c|c}
k& T/q^2\\ \hline
\frac{3a}{2}-1&-\frac{3(5a+2)}{44}\\[1mm]
\frac{9a-48}{23}&\frac{9(16a-1)}{529}\\[1mm]
-\frac a2-2&\frac a4\\[1mm]
-a-4&-1.
\end{array}
\]

Thus the affine-linear-in-\(T\) no-go from the fixed-mark audit is sharp:
the first polynomial escape occurs with marks linear in \(s\) and
\(T\) quadratic in \(s\).

## 2. A simultaneous conjugate section

The last row has rational coefficient \(T/q^2=-1\), so it is preserved by
the nontrivial automorphism of \(K/\mathbb Q\).  Set

\[
\begin{array}{c|cc}
&\text{first point}&\text{second point}\\ \hline
g& s&-(a+3)s\\
h& s&(a-2)s.
\end{array}
\]

At \(T=-s^2\), direct substitution proves

\[
\begin{aligned}
g_{-s^2}(-(a+3)s)-g_{-s^2}(s)
  -g_{-s^2}'(s)(-(a+4)s)&=0,\\
h_{-s^2}((a-2)s)-h_{-s^2}(s)
  -h_{-s^2}'(s)((a-3)s)&=0.
\end{aligned}
\]

Consequently both covers acquire polynomial tangent marks after the same
degree-two base change.  Since \(\mathrm{GL}_3(\mathbb F_2)\) is simple and
has no quotient of order two, this quadratic extension cannot be a
subextension of the common Galois closure.  The generic point/line
monodromy therefore remains \(\mathrm{GL}_3(\mathbb F_2)\).

## 3. Removing the spurious sixth-order collapse

For either row define

\[
H(W)=f_{-s^2}(q+dW)-f_{-s^2}(q)-d f_{-s^2}'(q)W.
\]

Both primitives have a common factor \(s^6\):

\[
H_g=s^6\overline H_g,\qquad H_h=s^6\overline H_h,
\]

with \(\overline H_g,\overline H_h\in K[s,W]\), of degree one in \(s\)
and degree seven in \(W\).  They retain

\[
\overline H(0)=\overline H'(0)=\overline H(1)=0.
\]

Their endpoint constants are

\[
\begin{aligned}
c_g&=-\overline H_g'(1)
 =28(a-10)(2s+1),\\
c_h&=-\overline H_h'(1)
 =-28(a+11)(2s+1),
\end{aligned}
\]

and their normalized second derivatives are

\[
\begin{aligned}
\kappa_g&=
\frac{2as+a-50s-11}{2(2s+1)},\\
\kappa_h&=
-\frac{2as+a+52s+12}{2(2s+1)}.
\end{aligned}
\]

Therefore the weighted construction extends across \(s=0\).  Its remaining
bad divisors are

\[
2s+1=0
\]

and

\[
\begin{aligned}
L_g&=2as+a-42s-7=0,\\
L_h&=2as+a+44s+8=0,
\end{aligned}
\]

where \(L_g\) and \(L_h\) are the numerators of
\(\kappa_g+2\) and \(-(\kappa_h+2)\).  They are conjugate up to sign, and

\[
L_gL_h=-2(928s^2+326s+29).
\]

On the common open

\[
\Omega=\operatorname{Spec}
K\left[s,\frac1{(2s+1)L_gL_h}\right],
\]

take \(b_g=c_g^{-1}\) and \(b_h=c_h^{-1}\).  The universal weighted
suspension theorem gives two relative maps

\[
F_g,F_h:\mathbb A^3_\Omega\longrightarrow\mathbb A^3_\Omega
\]

with relative Jacobian one.

These localizations are forced inside the universal weighted ansatz.  Its
Jacobian is \(b_\bullet c_\bullet\).  Because \(c_\bullet\) has degree one
in \(s\), no polynomial \(b_\bullet\in K[s]\) can make this product a
nonzero constant; determinant normalization necessarily uses
\((2s+1)^{-1}\).

The endpoint parameter is uniquely forced by polynomiality:

\[
a_\bullet=-\frac{1+\kappa_\bullet}{2+\kappa_\bullet}.
\]

For \(g\) and \(h\), respectively, its numerator and denominator resultants
are both

\[
\operatorname{Res}_s(N_\bullet,L_\bullet)=-56.
\]

Moreover

\[
\operatorname{Res}_s(2s+1,L_g)=28,\qquad
\operatorname{Res}_s(2s+1,L_h)=-28.
\]

Thus neither endpoint pole cancels, and neither coincides with the
determinant pole.  Removing these divisors requires leaving the
parameter-preserving universal weighted suspension, not just simplifying
its displayed formulas.

## 4. A two-chart marking atlas

The determinant and endpoint zeros above do not mean that the Davenport
polynomial has no admissible tangent pair at those parameters.  The second
row of the four-section table can also be placed over the same base
\(T=-s^2\).  Indeed, put

\[
\mu=\frac{5+2a}{3},\qquad
k_2=\frac{9a-48}{23},\qquad
\tau_2=\frac{9(16a-1)}{529}.
\]

Then

\[
\tau_2\mu^2=-1,\qquad k_2\mu=-a-4.
\]

Thus \(q_g=\mu s\), \(d_g=k_2\mu s\), together with the conjugate choices
for \(h\), give a second simultaneous tangent chart over \(T=-s^2\).

Up to nonzero constants, the complete bad polynomials of the first and
second charts are

\[
\begin{aligned}
B_4(s)&=(2s+1)(928s^2+326s+29),\\
B_2(s)&=(1516s^2-648s+729)
        (796300s^2+3240s+88209).
\end{aligned}
\]

Exact Euclidean reduction gives

\[
\gcd_{\mathbb Q[s]}(B_4,B_2)=1.
\]

Consequently

\[
D(B_4)\cup D(B_2)=\mathbb A^1_s.
\]

At every value of \(s\), at least one of the two simultaneous marking
charts is admissible.

The overlap is particularly simple.  If \(W_4\) and \(W_2\) denote the
root coordinates in the two charts, then

\[
\begin{aligned}
W_{2,g}&=W_{4,g}+\frac{3a+5}{21},\\
W_{2,h}&=W_{4,h}+\frac{2-3a}{21}.
\end{aligned}
\]

After these constant translations, the two primitives differ only by an
affine polynomial in \(W\).  This changes precisely the slope and intercept
of the inverse pencil and therefore gives a genuine two-chart atlas of the
same marked-root incidence covers.

This clears the existence-of-marks gap over the whole quadratic base.  What
is not yet proved is that the corresponding three-dimensional weighted
Keller suspensions glue through polynomial source and target automorphisms
to a trivial affine-space bundle.

There is already a sharp obstruction to the most natural lift.  Write

\[
\overline H_4(W)-\overline H_2(W+\delta)
=m(s)W+n(s).
\]

The exact formulas give nonzero \(m_g(s)\) and \(m_h(s)\).  Hence the
incidence slopes on the overlap satisfy

\[
\sigma_4=\sigma_2+m_\bullet(s).
\]

In the weighted target, however, \(\sigma=BC\).  Any transition retaining
the canonical boundary coordinate \(C\) is therefore forced to have

\[
B_4=B_2+\frac{m_\bullet(s)}{C}.
\]

Because \(m_\bullet(s)\ne0\), this is not polynomial across \(C=0\).
Thus the marked-root incidence atlas glues, but its canonical
\(C\)-preserving weighted lift is intrinsically Laurent.  A successful
absolute construction must change or modify the \(C=0\) boundary, rather
than merely glue the two standard weighted charts.

Allowing both weighted factors to change does not rescue a
product-compatible transition.  Such a transition would require two
polynomial coordinate functions \(\widetilde B,\widetilde C\) satisfying

\[
\widetilde B\widetilde C=BC+m_\bullet(s).
\]

But \(BC+m_\bullet(s)\) is irreducible in
\(K(s)[A,B,C]\), whereas the left side is a product of two nonunit
coordinate polynomials.  Hence no polynomial automorphism of the weighted
target can lift the overlap while retaining the presentation of incidence
slope as a product of two coordinates.  The remaining construction must
abandon that product presentation or use a genuine affine modification.

## 5. Recovery of the common Sunada fibers

Let \(u\) be the common Davenport target coordinate.  On \(C=1\), put

\[
B_\bullet=-\frac{d_\bullet f_\bullet'(q_\bullet)}{s^6},
\qquad
c_\bullet A_\bullet=
\frac{f_\bullet(q_\bullet)-u}{s^6}.
\]

Then

\[
\overline H_\bullet(W)-B_\bullet W+c_\bullet A_\bullet
=\frac{f_\bullet(q_\bullet+d_\bullet W)-u}{s^6}.
\]

For \(s\ne0\), the affine change
\(Y=q_g+d_gW\), respectively \(Z=q_h+d_hW\), identifies the two inverse
pencils with the original Davenport point and line covers pulled back by
\(T=-s^2\).  Hence every common good finite-field fiber has identical zeta
data, and the two generic covers retain their common
\(\mathrm{GL}_3(\mathbb F_2)\) Galois closure.

The Hessian-divisor comparison is unchanged by the affine changes in the
root coordinate, so the two relative Keller maps remain stably
inequivalent over the generic point.

## 6. What remains

This replaces the unspecified finite marking curve by an explicit quadratic
base change and improves the Keller base from that curve to a punctured
affine line.  It also proves that \(s=0\) is removable for the Keller seeds.

The two-chart atlas ensures that no point of the \(s\)-line lacks an
admissible tangent mark.  It does **not** by itself fill the divisors
\(2s+1=0\), \(L_g=0\), and \(L_h=0\) inside one weighted Keller chart.
The determinant-one normalization is forced to use
\(c_\bullet^{-1}\), while the uniquely determined weighted endpoint
cancellation uses \((\kappa_\bullet+2)^{-1}\).
Accordingly the total spaces are still affine-space bundles over
\(\Omega\), not affine spaces over \(K\).  An absolute construction must
absorb these remaining divisors and the overlap pole along \(C=0\) through
a nontrivial higher-dimensional modification or a different suspension.
Polynomial recoordination within the standard weighted-product model is
now ruled out.

The minimal affine-modification attempt is analyzed separately in the
[weighted-glue obstruction](DAVENPORT_WEIGHTED_GLUE_OBSTRUCTION.md).
The slope quotient alone gives an affine-space hypersurface, but the
intercept quotient has a comaximal center; adjoining both quotients forces
`C` to be a unit.

## 7. Verification

Run

```bash
.venv/bin/python scripts/verify_davenport_proportional_tangent_sections.py
```

The checker verifies the two-quintic factorization, irreducibility of
\(F_1\), the four polynomial sections, the simultaneous conjugate section,
the \(s^6\) division, all endpoint formulas, the two inverse-pencil
identities, the coprime two-chart cover of the \(s\)-line, and the constant
root-coordinate transition maps.
