# Restricting symmetric monodromy to PSL2(11) does not supply an affine target

Work over an algebraically closed field k of characteristic zero.

**KMFI17.** Reducing the known curve-induced A11 and S11 Keller
extensions to an exceptional transitive PSL2(11) subgroup by
Galois-closure base change gives target curves of genera 5678,
7331 and 9343 for the three inputs below. Their function fields
are not stably rational. Consequently none of these subgroup reductions
can yield a polynomial self-map of ordinary affine space, even after
arbitrary birational coordinate changes and common independent
stabilization. This excludes a construction route, not either abstract
degree-eleven PSL2(11) action.

**KMFI18.** For the two known A12 families, analogous subgroup
reductions in the natural degree-twelve PSL2(11) action have generic
target-pencil genera 527826 and 524161. These are genera over the
algebraic closure of the pencil parameter field. They do not establish
that the total two-dimensional target field is nonrational.
The subsequent [global two-form proof](SURFACE_FORMS.md), KMFI19,
now proves both total target surfaces are not stably rational,
closing these two specific affine-space routes.

The [checker](verify_subgroup_reduction.py), [exact model](subgroup_reduction_model.json)
and [receipt](subgroup_reduction_receipt.json) retain all class intersections,
all powers of each inertia generator, cyclic orbit profiles and
Riemann--Hurwitz calculations. Only the 660-element group is enumerated;
the large symmetric and alternating groups need not be enumerated.
Universal geometric conclusions use the written proofs below.
No formal verification or independent whole-proof review is asserted.

## 1. The subgroup-reduction construction and its extra obligation

Let M/K be the geometric Galois closure of a degree-n point extension
L/K, with group G acting faithfully on n inverse sheets.
Write L=M^H, and let J<G be a transitive subgroup on those same sheets.
The base change

\[
 K_J=M^J,\qquad L_J=LK_J=M^{H\cap J}
 \tag{1}
\]

has degree n and normal-closure group J.
Transitivity gives [J:H intersection J]=n, and faithfulness gives
trivial core of H intersection J in J.
In particular, the tensor product L tensor_K K_J is a field, so
no inverse component was silently discarded.

This is a field-theoretic group realization.
An affine-space Keller realization additionally requires K_J to be
a rational function field, as well as rational source coordinates and
a polynomial map with constant nonzero Jacobian.
Adding independent variables to a positive-genus curve field cannot
satisfy even the first of these obligations.

Here is the rationality obstruction used below.
If k(C)(u_1,...,u_r) were rational over k for a smooth projective
curve C of positive genus, there would be a dominant rational map
from projective space to C. Restricting to a suitable general line
would give a nonconstant rational map P1 to C. It extends to a
morphism since C is proper, contrary to Riemann--Hurwitz.
Thus the obstruction survives all independent stabilization.
It is intrinsic to the field, so a birational change of the
prospective affine coordinates cannot remove it.

If the original pair is induced from a fixed curve cover and then
independently stabilized, M itself is the Galois-closure curve
field with the same independent variables adjoined. Hence

\[
 M^J=k(C_J)(u_1,\ldots,u_r).
 \tag{2}
\]

A positive genus for C_J is therefore a complete obstruction to
this specific subgroup-reduction route.

## 2. Exact fixed-point and genus formulas

For an element g in a finite group G, its action on G/J has

\[
 \operatorname{Fix}_{G/J}(g)
   =\frac{|C_G(g)|\,|g^G\cap J|}{|J|}
 \tag{3}
\]

fixed cosets. For m=ord(g), the number of cycles of g is

\[
 c_{G/J}(g)=\frac1m\sum_{i=0}^{m-1}
                      \operatorname{Fix}_{G/J}(g^i).
 \tag{4}
\]

This cyclic Burnside formula handles composite inertia orders too:
using only the fixed points of g would be insufficient.
If N=[G:J], its cycle index is N-c_(G/J)(g). Thus a curve with
branch generators g_b gives

\[
 2g(C_J)-2=-2N+\sum_b\bigl(N-c_{G/J}(g_b)\bigr).
 \tag{5}
\]

For a symmetric-group cycle partition lambda, the centralizer order is
the product of j^(m_j) m_j!. An even class splits in A_n exactly when
all its cycle lengths, including fixed points, are distinct and odd.
In a split class the alternating centralizer is unchanged; otherwise
it is half the symmetric centralizer.

For the eleven-cycle classes used here, every cyclic order-eleven
subgroup of J contributes five generators to each split alternating
class. Indeed conjugation of an eleven-cycle to its a-th power is
even exactly when a is a quadratic residue modulo 11.
The checker verifies that sign on all ten multipliers.
It follows that a 120-element order-eleven census contributes
60 elements to each split class; a 10-element census contributes five.
No arbitrary choice of split class changes the genera below.

The complete cycle dictionaries for J=PSL2(11) are

| order | elements | exceptional degree 11 | natural degree 12 |
|---:|---:|---|---|
| 1 | 1 | 1^11 | 1^12 |
| 2 | 55 | 2^4 1^3 | 2^6 |
| 3 | 110 | 3^3 1^2 | 3^4 |
| 5 | 264 | 5^2 1 | 5^2 1^2 |
| 6 | 110 | 6 3 2 | 6^2 |
| 11 | 120 | 11 | 11 1 |

The checker constructs SL2(F11)/{+/-I}, its natural action and an
A5 subgroup of order 60, then constructs the degree-eleven coset
action. The two exceptional actions have the same cycle dictionary
by the retained [Gassmann calculation](../../archive/non-elliptic/extended-geometry/PSL2_11_KELLER_ACTION_SPECTRUM.md);
the split-class argument applies to either. Thus the conclusions
are independent of which exceptional action is requested.

The stabilizers J intersection H have orders 60 in degree eleven
and 55 in degree twelve. Their nonidentity class counts needed here
are 24 elements of type 5^2 1 for A5, and ten elements of type
11 1 for the natural point stabilizer. Their complete dictionaries
are also retained in the model.

## 3. Three actual degree-eleven inputs

The [higher-power radicand construction](ALL_ALTERNATING.md) gives
an actual A11 Keller map at (p,d)=(3,5), and an actual S11 map at
(p,d)=(2,10). The
[quadratic beta construction](ALTERNATING_REALIZATIONS.md) at p=6
gives another actual S11 map. Their independent curve presentations
are of the common form

\[
 \tau=wH_{p,d}(w)^d,\qquad
 H_{p,d}(w)=\sum_{j=0}^{p-1}\binom{p-1}{j}
       \frac{(-1)^{p-1-j}w^j}{dj+1}.
 \tag{6}
\]

For the p=6 quadratic beta input, the same field reconstruction
used in ALL_ALTERNATING.md applies without requiring p=2 or 3:
put Z=1+dQ, R=delta^p P/(dp)+integral_0^1(theta^d-Z)^(p-1)dtheta,
ell=R/Z^(p-1), tau=ell^d/Z and w=rho^d/Z.
The inverse is Z=ell^d/tau and rho=ell/H_(p,d)(w).
Together with the independent C, these formulas prove exactly
K=k(tau,ell,C) and L=k(w,ell,C).
Thus this input, too, is a fixed curve followed by independent
stabilization; genus is not being inferred from a merely chosen pencil.

The passports are d^(p-1) 1 | p 1^(11-p) | 11.
For J in its exceptional action, equations (3)--(5) give:

| actual input | fixed-curve passport | [G:J] | target genus | source genus |
|---|---|---:|---:|---:|
| A11, (p,d)=(3,5) | 5^2 1; 3 1^8; 11 | 30240 | 5678 | 62493 |
| S11, (p,d)=(6,2) | 2^5 1; 6 1^5; 11 | 60480 | 7331 | 80641 |
| S11, (p,d)=(2,10) | 10 1; 2 1^9; 11 | 60480 | 9343 | 102813 |

For the first row, the target inertia calculations are especially
transparent:

| original branch generator | fixed cosets | complete coset cycles | index |
|---|---:|---|---:|
| 5^2 1 | 10 | 5^6046 1^10 | 24184 |
| 3 1^8 | 0 | 3^10080 | 20160 |
| 11 | 1 | 11^2749 1 | 27490 |

Thus

\[
 2g(C_J)-2=-60480+24184+20160+27490=11354,
\]

so g(C_J)=5678. The corresponding source quotient by A5 has genus
62493. Its degree-eleven map to C_J has ramification degree
124984-11*11354=90, a separate Riemann--Hurwitz control.

For the second row, a single six-cycle and all of its nontrivial
powers have cycle partitions absent from J. It acts freely on
the 60480 cosets; so does the odd 2^5 1 inertia.
The eleven-cycle has two fixed cosets in the symmetric ambient
group. These give genus 7331.

For the third row, the tenth-order generator has a genuine
short-orbit contribution. Its four fifth-order powers each have
20 fixed cosets, although the generator itself fixes none.
Its complete cycles are 10^6046 2^10 and its index is 54424.
The other indices are 30240 and 54980, giving genus 9343.
This illustrates why a prime-order shortcut would give an incorrect
genus for this route.

**KMFI17 follows.** In every row the new target, and also the new
source, has positive curve genus before independent parameters
are adjoined. Equations (1)--(2) construct the desired smaller
generic group at the field level, but the rationality obstruction
in Section 1 prevents an affine-space polynomial realization of
these field pairs in every independently stabilized birational model.

This does not exclude an unrelated genus-zero PSL2(11) input such
as the retained Klein Shabat cover with passport (3^3 1^2)|(2^4 1^3)|11.
It does not exclude a genuinely higher-dimensional PSL2(11) field pair.
There is also no general claim that a subgroup of a realized
monodromy group is itself realized by an affine-space Keller map.

## 4. The A12 reductions give different information

There are two actual A12 families already verified:

* the level-two cone at d=6, whose generic Q-pencil has six single
  three-cycle branch generators and one (11,1) generator;
* the alternating shear, whose generic Y-pencil has six single
  three-cycle generators and one (9,3) generator.

Use the natural degree-twelve embedding J=PSL2(11)<A12,
with point stabilizer of order 55. The same exact class-intersection
calculation gives

| input pencil | [A12:J] | generic target-curve genus | generic source-curve genus | relative ramification degree |
|---|---:|---:|---:|---:|
| level-two cone | 362880 | 527826 | 6333906 | 10 |
| alternating shear | 362880 | 524161 | 6289921 | 0 |

For the cone, each finite three-cycle acts freely on G/J and
has index 241920. The (11,1) element has one fixed coset and
index 329890. Hence 2g-2=-725760+6*241920+329890=1055650.
For the shear, every nontrivial power of (9,3) has cycle partition
absent from J; it acts freely and has index 322560. This gives
2g-2=1048320. The reduced shear curve cover is unramified,
as the source and target genera independently verify.

**KMFI18 is a generic-fiber genus calculation.**
Here the fixed fields M^J are two-dimensional over k before the
independent C is adjoined. We have computed their generic fiber
genera over the algebraic closure of k(Q) or k(Y).
We have not identified these total fields with the independent
stabilization of a fixed positive-genus curve.

Positive generic-fiber genus does not imply that a surface field
is nonrational over k. For example k(x,y) is a rational surface
field, while the function t=y^2-x^3 makes its generic fiber
y^2=x^3+t a genus-one curve. Therefore the argument proving
KMFI17 cannot be applied to these A12 surfaces.
The subsequent [KMFI19 proof](SURFACE_FORMS.md) constructs global
J-invariant regular two-forms and now proves that both total
surfaces are not stably rational. This supplies the previously
missing field-intrinsic obstruction; the generic-fiber genus
calculation alone still does not imply it.

## 5. Arithmetic-parametric statements retain their stronger scope

Geometric occurrence of A5, or of every natural A_n for n>=5,
does not by itself produce one polynomial map over Q whose
generic arithmetic group is that same group and whose complete
fibers realize every point field of every corresponding G-extension.

The original [Beckmann--Black Keller chart result](../../archive/non-elliptic/KELLER_BECKMANN_BLACK_SPECIALIZATION.md)
asks for this stronger uniform arithmetic property.
Its absolute arithmetic-parametric frontier remains open for the
four pilot actions, even though geometric A5 occurrence is now
proved. The current ledger makes that distinction explicit.

The full action spectrum and both specified PSL2(11) realization
problems remain open. KMFI17 closes the three curve routes in Section 3; KMFI19 now
closes both additional surface routes in Section 4.
