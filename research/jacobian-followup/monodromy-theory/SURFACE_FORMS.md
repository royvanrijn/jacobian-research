# Global two-forms obstruct both A12-to-PSL2(11) target surfaces

Work over an algebraically closed field k of characteristic zero.
Let M be the Galois closure of either actual A12 plane field pair
from KMFI13 or KMFI15, and let J=PSL2(11)<A12 be the natural
degree-twelve subgroup. The independent threefold coordinate C
is omitted until the conclusion.

**KMFI19.** Each target surface field M^J has a smooth projective
model with a nonzero regular two-form. In particular neither field
is stably rational. No common independent stabilization or birational
change of coordinates can turn either of these subgroup-reduced
field pairs into a polynomial map between ordinary affine spaces.

This closes the two surface routes left open in
[KMFI18](SUBGROUP_REDUCTION.md). Combined with KMFI17 it closes all
five subgroup reductions audited there. It does not exclude an
unrelated PSL2(11) field pair or either abstract action, and does
not complete the action-spectrum classification.

The proof gives explicit J-invariant rational two-forms, resolves
all their possible poles, and proves nonvanishing by exact good-prime
certificates. [The checker](verify_surface_forms.py),
[model](surface_forms_model.json), and [receipt](surface_forms_receipt.json)
retain the branch polynomials, local extension bounds and certificates.
This is a written proof with exact controls, not formal verification
or independent whole-proof review.

The [orbit-module extension](ORBIT_FORMS.md), KMFI20--22,
now supplies a general subset-orbit criterion. For these two
normal closures it excludes every proper transitive subgroup
reduction, including both Mathieu groups, and gives source
obstructions for root-configuration actions.

## 1. Two elementary differential-form principles

First, a regular top form invariant under a finite group on a
smooth surface extends to a resolution of its coarse quotient.
Here is a direct argument which also covers the finite quotient
charts used below. Let U be a smooth local chart, eta an invariant
regular two-form, and Y a resolution of U/H. Normalize U times_(U/H) Y
and resolve it, obtaining V with morphisms to U and Y.
The pullback of eta to V is regular. If its descended rational form
on Y had a pole of order a>=1 on a divisor, its pullback on any
divisor above it with ramification index e would have order
-e*a+(e-1)<0, since ramification is tame in characteristic zero.
This is impossible. Thus the form extends.
The same argument applies to smooth finite-quotient stack charts.

Second, a smooth projective surface with a nonzero regular two-form
is not stably rational. Regular differential forms on smooth
projective varieties are birational invariants; adjoining independent
variables corresponds to taking a product with projective space,
which preserves a nonzero pulled-back two-form. Projective space
has no such form. Thus a nonzero form also obstructs a rational
affine-space target after any independent stabilization.

## 2. A J-invariant cut function on the twelve inverse roots

Label the roots w_0,...,w_11 by P1(F11), with label 11 for infinity,
and let O be the J-orbit of {0,1,2,5}. It has 165 four-element subsets.
For A in O put

\[
 B_A(w)=\prod_{i\in A,\ j\notin A}(w_i-w_j),\qquad
 U_A(w)=\sum_{i\in A}w_i .
\]

Each B_A has degree 4*8=32. Define

\[
 F_{10}(w)=\sum_{A\in O}\frac{U_A(w)^{10}}{B_A(w)},\qquad
 F_0(w)=\sum_{A\in O}\frac1{B_A(w)}.
 \tag{1}
\]

These are invariant under J; no sign character enters the products,
whose orientation always runs from the four-set to its complement.
They have homogeneous degrees -22 and -32 respectively.

The useful local bound is elementary. If r roots form a collision
cluster and A contains k of them, precisely k(r-k) factors of B_A
come from differences within that cluster. For r=3,4,5,7 the maxima,
with |A|=4, are 2,4,6,12. Differences to other roots are units.
If the local root difference has order h, a denominator contributes
at most h*k(r-k). The numerator U_A^10 never introduces a pole.
The proof checks this bound term by term, so cancellation in the
orbit sum is not needed for regularity.

All natural embeddings of J are conjugate in A12. Indeed their
S12 conjugacy class does not split: the projective multiplier
x -> 2x is an odd permutation normalizing J. It also preserves
the 165-element orbit. Consequently choosing the other component
of the ordered-root cover or relabeling the natural action does
not affect the existence of the form constructed below.

## 3. The compactified root equations

Normalize the cone seed aW^2+bW^6 to W^2+W^6 by constant source
and target rescaling; a,b are nonzero. This is legitimate over k
and preserves the target function field up to k-isomorphism.
Use homogeneous coordinates (s,Q,P) of weights (1,6,11).
The cone root equation is

\[
 D_c(X)=Q^2+2s^4QX^2+\frac25QX^6-\frac13s^8X^4
        -\frac27s^4X^8-\frac1{11}X^{12}-PX.
 \tag{2}
\]

Every term has weight twelve, with X of weight one.
On s=1 it is exactly X(Pi(X,Q)-P).
The leading coefficient is the nonzero constant -1/11.

Normalize the shear coefficients c=d=1 by the constant rescaling
from KMFI13. Its compactification has weights (s,Y,P)=(1,2,9):

\[
\begin{split}
 D_s(X)={}&\frac19X^{12}-\frac67YX^{10}+3Y^2X^8
 -\frac{20}3Y^3X^6+15Y^4X^4+6Y^5X^2-\frac13Y^6\\
 &+s^4\left(2Y^3X^2+6Y^2X^4-2YX^6+\frac25X^8\right)
 +s^8X^4-PX^3 .
\end{split}
 \tag{3}
\]

Again the equation has weight twelve, with constant leading
coefficient 1/9, and its s=1 chart is the verified shear root equation.

Take the normalization of each weighted projective plane in M,
or equivalently the appropriate component of its ordered-root cover.
The roots are integral sections of weight one, since the equations
are monic after multiplication by a constant.
Thus they have no finite poles on the weighted local charts.
The only possible poles of (1) come from collisions of roots.
The two-form analysis below uses smooth charts for weighted planes
and weighted blowups; residual stabilizers are finite, so Section 1
handles their coarse quotient resolutions.

## 4. Complete finite branch geometry

Both branch curves admit a rational parametrization for which,
putting z equal to the square of the parameter,

\[
 Q=z^3+z,\qquad P^2=\mu(z).
\]

For the cone,

\[
 \mu_c(z)=z^3\left(\frac83+\frac{144}{35}z^2
                         +\frac{72}{55}z^4\right)^2.
\]

For the shear, replace Q by Y and take

\[
 \mu_s(z)=\frac{64}{315^2}z^3
 (105+630z^2+2205z^4+4452z^6+5040z^8+2880z^{10}
                      +640z^{12})^2 .
\]

Let R(V,Q) be the monic cubic resultant of z^3+z-Q and V-mu(z).
Then Gamma(P,Q)=R(P^2,Q) is the reduced branch curve.
The critical resultants of (2) and (3), on s=1, are respectively
-Q*Gamma and -Y^9*Gamma. The derivative identities of KMFI13 and
KMFI15 therefore give root discriminants, up to nonzero constants,

\[
 \operatorname{disc}(D_c)=\Gamma_c^2,\qquad
 \operatorname{disc}(D_s)=Y^{12}\Gamma_s^2.
 \tag{4}
\]

The shear factor Y=0 is an apparent collision, not generic
ramification of the normalized root field; it is treated below.

The exact cubic resultants give

\[
 R_c(0,Q)=u Q^3 C_c(Q)^2,\quad
 C_c=6751269Q^4+22773960Q^2+154000,
\]
\[
 \operatorname{disc}_V R_c
 =v(27Q^2+4)^5 E_c(Q)^2,\quad
 E_c=62523502209Q^6-132324872400Q^4+639360000Q^2+4000000.
\]

For the shear they give the same pattern, with

\[
\begin{split}
 C_s(Y)={}&262144000Y^{12}+309657600Y^{10}
 +170311680Y^8+54441408Y^6\\
 &+10734885Y^4+1277640Y^2+82320,\\
 E_s(Y)={}&5832000000Y^8+4682957625Y^6+1761706800Y^4
             +331914240Y^2+30118144 .
\end{split}
\]

All packet polynomials are squarefree, nonzero at zero, and
pairwise coprime with each other and 27Q^2+4.
The two depressed-cubic triple-root equations are coprime,
and R(0,Q) and R_V(0,Q) are coprime.
These exact identities are retained and checked over Q.

They imply the following complete finite geometry:

| family | ordinary nodes | cusps (2,5) | origin cusp (2,3) |
|---|---:|---:|---:|
| cone | 4+2*6=16 | 4 | 1 |
| shear | 12+2*8=28 | 4 | 1 |

To justify the local assertions, away from the zeros of the
parameter derivative all branches are analytic functions of Q.
A simple zero of C gives a crossing of the opposite values
+P and -P at P=0. A simple zero of E gives a simple crossing of
two distinct nonzero squared values and hence two ordinary nodes,
one at each sign of P.
The exponent two in the discriminant shows that their slopes
differ. The coprimality checks exclude triple crossings, a
simultaneous zero and repeated squared value, and intersections
with the critical-parameter packet.

At a nonzero zero of the parameter derivative, the original
rational map has derivative equal to a unit times
(q-alpha*z^2)^2 in centered source coordinates. Subtracting an
analytic function of q from P gives initial polynomial
p=z^5/5-2alpha*q*z^3/3+alpha^2*q^2*z, up to nonzero scalings.
Thus the branch cusp is (2,5). There are four such source parameters.
At zero, the branch parametrization starts in orders (2,3), giving
the remaining cusp.

This classification can also be checked globally. The cone branch
curve has weighted degree 66 in P(1,6,11), arithmetic genus 25,
and finite delta sum 16+4*2+1=25. It is smooth at its weighted
boundary. The shear has degree 54 in P(1,2,9), arithmetic genus 64;
its finite delta sum is 28+4*2+1=37, and its sole boundary cusp
has delta 27, as proved in Section 6. The total is again 64.
Both parametrizations are birational onto their branch curves.
For the cone this follows from the coprime parameter degrees six
and eleven. For the shear, the fiber over (P,Y)=(0,0) consists
only of the zero parameter: at the other solutions of Y=0 the
polynomial in the shear parametrization takes value 28.
The local orders at zero are two and three, so the degree of
this finite parametrization divides their gcd and is one.

The two arithmetic genera equal the counts of weighted monomials
of degrees 48 and 42 respectively; no orbifold correction occurs
because the branch curves avoid the stack stabilizer points.

The shear polynomials C_s,E_s also agree with the algebraic
branch-collision packets in the retained Kistner--Shaska source.
This proof does not assume its conditional admissible-cover
Hypothesis 11.3: only the explicit root equations and the exact
resultants above are used.

## 5. Extension at every finite collision

At a smooth branch point, three roots collide and inertia is a
single three-cycle. A local uniformizer u makes their differences
order one and the pulled-back base two-form order two.
The cut denominator has order at most two, so either proposed
form extends.

At an ordinary node there are two disjoint three-root clusters.
The normalized Galois cover has local coordinates u,v with the
two base branch parameters u^3,v^3. Its volume has orders two
in each parameter. The cut denominator has order at most two
in each as well, so the same argument applies.

At a finite (2,5) cusp, make the weighted blowup with weights
(q,p)=(2,5). On its smooth stack charts the five roots in the
cluster differ in order one; the base volume has order
2+5-1=6. The maximum cut is 2*3=6.
The normalized cover is smooth on those charts away from
a simple three-cycle branch: the leading five-root polynomial
is separable on the exceptional curve away from that branch,
and its branch transform is smooth and transverse.
The weighted axis charts have distinct roots too.
The three-cycle calculation then handles the remaining intersection.
Finite stabilizers introduce no further poles by Section 1.

At the cone origin, (2) has a four-root cluster with source-root
weight one and base weights (Q,P)=(2,3).
The eight spectator roots are distinct and nonzero.
The volume order is four and the maximum cut is four.
The same weighted-blowup argument applies, with a simple
three-cycle branch on the transformed chart.

The shear needs an additional factor for its lost roots.
For Y=0, P nonzero, three roots have X=Y^2 Z. Dividing the
root equation by Y^6 leaves the distinct residual roots

\[
 -P Z^3-\frac13=0.
\]

Thus their pairwise differences have Y-order two although the
normalized field is unramified there. Multiplication by Y^4
cancels the worst cut, two edges each of order two.
This still works where Y=0 meets a distinct three-cycle branch:
the two clusters are disjoint and both cancellations apply.

At the shear origin, four roots have weight three when
(Y,P) have weights (2,3). The volume order is four,
and Y^4 has order eight, totaling twelve.
The worst cut has four edges of order three, also twelve.
After this weighted blowup, the axis Y=0 has a three-root
apparent collision again; the residual factor Y^4 on the chart
cancels it exactly as above. The other axis is separable.
All remaining branch points are simple three-cycle points.
The spectator roots are again distinct and nonzero.

These arguments give an explicit smooth finite-quotient
modification of the normalized cover near every finite point,
with the forms regular there. In particular, they address
exceptional divisors over singular points, not only generic
points of the original branch divisor.

## 6. Extension at the projective boundary

Let E be the Euler vector field for the indicated weights and set

\[
 \widetilde\Omega_c
  =s^4 F_{10}(w)\,\iota_E(ds\wedge dP\wedge dQ),\qquad
 \widetilde\Omega_s
  =s^{12}Y^4F_0(w)\,\iota_E(ds\wedge dP\wedge dY).
 \tag{5}
\]

Both expressions have total weight zero and are horizontal, hence
descend from the homogeneous cone.
Indeed their weight sums are respectively
4-22+(1+11+6)=0 and 12+8-32+(1+9+2)=0.
On s=1 they are F_10 dP wedge dQ and Y^4 F_0 dP wedge dY.

For the cone, the branch curve is smooth and transverse to s=0.
Its homogeneous parametrization is

Q=W^6+s^4 W^2,
P=(8/3)s^8 W^3+(144/35)s^4 W^7+(72/55)W^11.

At s=0 one has W nonzero. All roots outside the unique
three-root branch cluster are distinct. The ordinary
three-cycle estimate applies at the crossing with the boundary,
and the extra factor s^4 is nonnegative.
The weighted axis charts are disjoint from the branch and have
distinct roots. They only require finite quotient descent.

For the shear, on the boundary chart Y=1 the leading Laurent
primitive has

\[
 P_0'(X)=(X^2-1)^6/X^4,\qquad P_0(1)=1024/63.
\]

At each of the equivalent points X=+1 and X=-1, exactly seven
roots collide; five spectators are distinct. Center at X=1,
write z=X-1, and subtract P_s(1) from P.
The exact weighted initial form, with wt(s,z,p)=(3,4,28), is

\[
 p=s^8z+4s^4z^4+\frac{64}{7}z^7,\qquad
 p_z=(s^4+8z^3)^2.
 \tag{6}
\]

Higher terms have strictly larger weight.
The critical branch has s of order three and p of order 28,
with nonzero initial coefficients. Since gcd(3,28)=1, this is
a single (3,28) cusp, of delta (3-1)(28-1)/2=27.

Make the weighted blowup (s,p)=(3,28). Root differences in
the seven-cluster have order four.
The volume order is 3+28-1=30; s^12 adds 36.
The largest cut is 4*3=12, so the possible denominator order
is at most 48, below the available 66.
After normalization, the leading seven-root polynomial in (6)
is separable away from a simple three-cycle branch on the
exceptional curve. The transformed branch is smooth and
transverse; both weighted axis polynomials are separable.
Thus the same smooth-chart and three-cycle arguments complete
extension there. The calculation is unchanged at X=-1.

At the boundary axis Y=0,P nonzero, the three small roots still
have X=Y^2 Z with distinct residual roots. The factor Y^4 cancels
their cut denominator just as in Section 5, and s^12 adds zeros.
The other weighted axis is separable. This covers every boundary
point, so (5) extends globally on a smooth finite-quotient model.

## 7. Exact nonzero certificates

Regularity alone would allow the orbit sums to be identically
zero on this special two-dimensional family. We exclude that
with two exact specializations, with all denominator factors units.

For the cone take s=Q=0 and P=-1/11 in (2). The monic root
polynomial is X(X^11-1). Over F_991 take omega=118 of order eleven,
label the first eleven roots by omega^i and the last by zero.
The orbit sum F_10 is

\[
 F_{10}=473\pmod{991}\ne0 .
\]

For the shear take s=P=0,Y=1. Multiplying (3) by 63 gives

\[
 N(X)=7X^{12}-54X^{10}+189X^8-420X^6
       +945X^4+378X^2-21 .
\]

It splits with twelve distinct roots over F_179[h]/(h^2-2).
In the label order, write the roots as pairs (a,b) for a+bh:

(0,4), (0,175), (21,1), (21,178), (34,18), (34,161),
(43,0), (136,0), (145,18), (145,161), (158,1), (158,178).

The checker substitutes every root, verifies distinctness,
and evaluates every cut denominator exactly. It gives

\[
 F_0=173+114h\ne0 .
\]

These are good-prime certificates of nonzero characteristic-zero
rational functions, not numerical estimates. The monic root
schemes are etale at the specialized fibers; the listed roots
lift in characteristic zero, and a rational function with unit
denominator and nonzero reduction cannot vanish identically.
The finite orbit computation and the odd normalizer ensure the
chosen labeling causes no component ambiguity.
Although s=0, the nonzero leading coefficient of the respective
power of s proves that the global rational two-form is nonzero.

## 8. Descent and the affine-space obstruction

The forms (5) are regular by Sections 4--6 and invariant under J
by Section 2. They descend to nonzero regular two-forms on smooth
projective resolutions of the target surfaces M^J, using Section 1.
Thus both have geometric genus at least one and are not stably
rational. The proof does not infer this from the positive genera
of the pencils computed in KMFI18.

The full threefold subgroup target is M^J(C), with C independent,
so it too remains nonrational after every further independent
stabilization. An ordinary affine-space target has a rational
function field. Therefore neither subgroup-reduced field pair
admits a polynomial affine-space model, regardless of subsequent
Jacobian or source-coordinate adjustments.

This proves KMFI19. Together with KMFI17 it excludes all five
audited subgroup-reduction inputs to PSL2(11).
The independent Klein Shabat completion and genuinely different
higher-dimensional field pairs are not addressed by this
obstruction. Their occurrence or action-level exclusion and a
complete classification remain open.
