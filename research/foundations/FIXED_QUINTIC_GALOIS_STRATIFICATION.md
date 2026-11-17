# Galois incidences in the target of the fixed quintic Keller map

## Status

This note is an exact atlas and a frontier statement.  It proves
pairwise-nonisomorphic specialization theorems for `S_5`, `F_{20}`, and
`D_5`: a generic target family, a rational De Moivre surface, and a rational
Brumer curve, respectively.  It also identifies and sharpens the two
incidence covers that must still be parametrized for `A_5` and `C_5`.
For `A_5` it gives fixed-map anchors in both real signatures, two exact
Mestre source pencils, and a genus-eight obstruction to lifting the
totally-real pencil by affine changes of generator.  It determines the
group-theoretically possible real signatures and unramified splitting
types.

It does **not** claim infinitely many fixed-map fields for `A_5` or `C_5`.
For those groups the missing step is a rational curve or surface on the
relevant incidence cover with the indicated generic group.  In particular,
the five groups do not form five locally closed algebraic strata of the
rational target.  Square classes, rational resolvent roots, resolvent
factorizations, and rational automorphisms are thin arithmetic conditions.

The fixed map and the isolated certificates are in
[`FIXED_QUINTIC_MODULI_DOMINANCE.md`](FIXED_QUINTIC_MODULI_DOMINANCE.md).
The new discriminant-branch and group-theoretic assertions are checked by

```bash
.venv/bin/python scripts/verify_fixed_quintic_galois_stratification.py
```

## 1. The normalized inverse family

On `\Pi\ne0`, use

\[
 f_{\Pi,B,C}(T)
 =
 T^5-5T^3-2\Pi BT^2+4\Pi^3T-2\Pi^5C.
\]

Let

\[
 U=\{(\Pi,B,C):\Pi\Delta(\Pi,B,C)\ne0\},
\]

where

\[
 \operatorname{Disc}(f_{\Pi,B,C})=16\Pi^8\Delta(\Pi,B,C).
\]

The polynomial `\Delta` is the one displayed in Section 9.1 of the main
note.  Every irreducible specialization in `U(Q)` is a connected full fiber
of the fixed Keller map.  Thus the arithmetic problem may be carried out
entirely with this normalized polynomial.

Write `R_{10}^F(X)` for the pulled-back pair-sum resolvent in Section 9.3,
and write

\[
 \mathcal D_6^F(X)
 =
 \mathcal D_6(-5,-2\Pi B,4\Pi^3,-2\Pi^5C;X)
\]

for Dummit's solvability sextic.  This is an explicit polynomial in
`\mathbb Z[\Pi,B,C,X]`.  For orientation, its leading terms are

\[
\begin{aligned}
\mathcal D_6^F(X)
={}&X^6+32\Pi^3X^5\\
&-40\Pi^2(B^2+5BC\Pi^4-16\Pi^4+15\Pi)X^4+\cdots .
\end{aligned}
\]

The full universal formula and its exact quintic certificates are checked
by
[`scripts/verify_universal_quintic_calculator.py`](scripts/verify_universal_quintic_calculator.py).

## 2. Classification on rational points

Restrict to targets `y\in U(Q)` for which `f_y` is irreducible.  The
transitive-subgroup classification of `S_5`, the pair action, and Dummit's
criterion give the following exact decision tree.

| group | arithmetic condition |
|---|---|
| `C_5` | `R_{10}^F` factors as `5+5`, and `Q[T]/(f_y)` has a nonidentity `Q`-automorphism |
| `D_5` | `R_{10}^F` factors as `5+5`, and `Q[T]/(f_y)` has no nonidentity `Q`-automorphism |
| `F_{20}` | `\Delta(y)` is nonsquare and `\mathcal D_6^F(y;X)` has a rational root |
| `A_5` | `\Delta(y)` is a square and `R_{10}^F(y;X)` is irreducible |
| `S_5` | `\Delta(y)` is nonsquare and `\mathcal D_6^F(y;X)` has no rational root |

Here “factors as `5+5`” means a product of two irreducible rational
quintics.  On the irreducible-quintic locus, reducibility of the pair-sum
resolvent already leaves only `C_5` and `D_5`.

For the last distinction, a nonidentity automorphism of a degree-five field
has order five.  It therefore makes the field normal and cyclic.  Conversely,
a cyclic quintic field supplies such an automorphism.

This table is a classification of rational points, not a decomposition into
Zariski strata.  For example, “`\Delta` is a rational square” is represented
geometrically only after passing to the double cover `W^2=\Delta`.

## 3. The four exceptional incidence covers

The conditions above pull back to explicit finite-type incidence varieties.
These are the spaces on which a rational curve or surface must be built.

### 3.1 Alternating incidence

\[
 \boxed{\mathcal I_{A}: W^2=\Delta(\Pi,B,C).}
\]

On the algebraic open `\Pi W\ne0`, restrict to rational points for which
`f` and `R_{10}^F` are irreducible.  Those points give exactly `A_5`
fibers.  The isolated `A_5` certificate proves that this arithmetic subset
is nonempty, but it does not prove infinitude or Zariski density of its
rational points.

### 3.2 Frobenius incidence

\[
 \boxed{\mathcal I_{F}: \mathcal D_6^F(X)=0.}
\]

Remove the ramification divisors of the incidence projection.  Among the
remaining rational points, those for which `f` is irreducible and `\Delta`
is nonsquare give exactly `F_{20}` fibers.  The known Dummit root `X=-13/2`
supplies one rational point.

### 3.3 Dihedral incidence

Introduce two monic quintics

\[
\begin{aligned}
 U_5(X)&=X^5+u_4X^4+\cdots+u_0,\\
 V_5(X)&=X^5+v_4X^4+\cdots+v_0,
\end{aligned}
\]

and equate coefficients in

\[
 \boxed{\mathcal I_D:\quad R_{10}^F(X)=U_5(X)V_5(X).}
\]

This gives ten explicit equations over the target; quotienting by the
involution exchanging `U_5` and `V_5` forgets the chosen pair orbit.  The
rational points for which `f`, `U_5`, and `V_5` are irreducible and which
do not lift to the cyclic automorphism incidence below give exactly `D_5`
fibers.

### 3.4 Cyclic incidence

Let `A=\mathbb Q[T]/(f_{\Pi,B,C})` and introduce

\[
 \sigma(T)=s_0+s_1T+s_2T^2+s_3T^3+s_4T^4.
\]

Reduce `f(\sigma(T))` and
`\sigma^{\circ5}(T)-T` modulo `f(T)` and set their ten remainder
coefficients equal to zero.  After removing `\sigma(T)=T` and the
noninvertible locus, this defines

\[
 \boxed{\mathcal I_C:\quad
 f(\sigma)=0,\qquad \sigma^{\circ5}(T)=T\pmod f,\qquad \sigma\ne1.}
\]

On the irreducible locus, its rational points are exactly cyclic quintic
fibers.  This automorphism incidence is preferable to an absence-of-
Frobenius test: the latter can never certify cyclicity in a family.

## 4. A rational `F_{20}` surface

The classical De Moivre family admits a particularly small non-affine
Tschirnhaus lift.  Let `t,\pi` be independent rational parameters and put

\[
 u=\frac{1-t^2}{1+t^2},\qquad
 v=\frac{2t}{1+t^2},\qquad u^2+v^2=1,
\]

\[
 \boxed{
 h=\frac{\frac45\pi^3-1+3u^2v^2}{uv}.
 }
\]

On the open where `\pi tuv\ne0`, let

\[
 G_h(X)=X^5-5X^3+5X+h
\]

and let `\theta` denote the class of `X`.  Define

\[
 \boxed{\eta=u\theta+v(\theta^2-2).}
\]

Direct trace calculation gives

\[
 \operatorname{Tr}(\eta)=0,\qquad
 \operatorname{Tr}(\eta^2)=10.
\]

Its characteristic polynomial is

\[
 \chi_\eta(T)=T^5-5T^3+qT^2+4\pi^3T+s,
\]

where

\[
 \boxed{q=5uv(vh-2u)}
\]

and

\[
\boxed{
\begin{aligned}
s={}&u^5h+10u^4v+5u^3v^2h-10u^2v^3\\
   &-5uv^4h+2v^5-v^5h^2.
\end{aligned}}
\]

Consequently the rational map

\[
\boxed{
(t,\pi)\longmapsto
(\Pi,B,C)=
\left(\pi,-\frac q{2\pi},-\frac s{2\pi^5}\right)
}
\]

lands in the target of the one fixed Keller map, and its normalized inverse
polynomial is exactly `\chi_\eta`.

This is genuinely a surface in the target.  At `(t,\pi)=(1/2,1)`,

\[
 \det\frac{\partial(\Pi,B)}{\partial(\pi,t)}
 =-\frac{1884}{625}\ne0.
\]

The source discriminant is

\[
 \boxed{\operatorname{Disc}(G_h)=3125(h^2-4)^2
 =5\bigl(25(h^2-4)\bigr)^2.}
\]

The identity

\[
 (r-s)^5+5rs(r-s)^3+5(rs)^2(r-s)=r^5-s^5
\]

shows that `G_h` is solvable by a quadratic equation and radicals.  At

\[
 (t,\pi)=\left(\frac12,1\right)
\]

one has

\[
 (u,v,h)=\left(\frac35,\frac45,\frac{307}{300}\right)
\]

and

\[
 300X^5-1500X^3+1500X+307
\]

is irreducible modulo `11`.  Its discriminant has nonsquare class `5`.
The transitive solvable group is therefore neither `C_5` nor `D_5`, so it is
`F_{20}`.  Specialization embeds this group into the generic group of the
surface, while the De Moivre identity makes the generic group solvable.
Thus the generic group over `Q(t,\pi)` is exactly `F_{20}`.

The corresponding fixed-map target is

\[
 \boxed{
 \left(1,\frac{286}{625},
 \frac{67834669}{140625000}\right),
 }
\]

with primitive normalized polynomial

\[
\begin{aligned}
70312500T^5-351562500T^3-64350000T^2\\
{}+281250000T-67834669.
\end{aligned}
\]

It is also irreducible modulo `11`, giving a certificate directly in the
fixed target.

### 4.1 Signatures and unramified local types on the surface

The identity

\[
 X^5-5X^3+5X=2\cos(5\vartheta)
 \quad\text{when }X=2\cos\vartheta
\]

shows that `G_h` has five real roots for `|h|<2` and one real root for
`|h|>2`.  Both chambers meet the rational surface: with `t=1/2`, the choices
`\pi=1` and `\pi=2` give respectively

\[
 h=\frac{307}{300}\quad\text{and}\quad h=\frac{1269}{100}.
\]

Thus both group-theoretically possible `F_{20}` signatures occur on
infinitely many specializations.

At the displayed `\pi=1` target, exact squarefree factorizations give all
four possible unramified types:

\[
\begin{array}{c|c}
p&\text{factor degrees}\\ \hline
7&(4,1)\\
11&(5)\\
29&(2,2,1)\\
89&(1,1,1,1,1).
\end{array}
\]

Because the discriminant square class on this surface is `5`, the sign of
unramified Frobenius is `(5/p)`.  In particular, at primes where `5` is a
nonsquare the only possible `F_{20}` type on this surface is `(4,1)`.
At split-sign primes the three even types occur.  A finite collection of
these conditions is imposed by taking the corresponding residue opens in
the rational parameter surface.

### 4.2 Pairwise nonisomorphic `F_{20}` fields

Fix `t=1/2`.  Then

\[
 h=\frac53\pi^3-\frac{193}{300}.
\]

Let `p\equiv2\pmod3` avoid the finite set

\[
 \{2,3,5,7,11,13,31,61\}.
\]

The cube map is bijective modulo `p`, so Hensel's lemma gives
`\pi_p\in\mathbb Q_p` with

\[
 \pi_p^3=\frac{793}{500}+\frac35p,
 \qquad h(\pi_p)=2+p.
\]

Modulo `p`,

\[
 G_{2+p}(X)\equiv
 (X+2)(X^2-X-1)^2.
\]

The quadratic is separable and coprime to `X+2`.  At each of its geometric
roots the perturbation by `p` has a unit quadratic term and constant term
of valuation one.  The corresponding local factor is therefore a ramified
quadratic Eisenstein factor.  Hence the specialized quintic algebra is
ramified at `p`.  This remains true on a sufficiently small parameter
neighborhood of `(1/2,\pi_p)`.

Apply Hilbert irreducibility with local conditions on the rational
`(t,\pi)`-surface.  After finitely many fields have been selected, choose a
new prime of the above congruence class that is unramified in all of them,
and impose its ramified parameter neighborhood.  The next specialization
has group `F_{20}` and ramifies at the new prime.  Induction proves:

> **Theorem.**  The fixed Keller map has infinitely many pairwise
> nonisomorphic connected full `F_{20}` fibers with five real roots, and
> infinitely many with one real root.  Any fixed finite collection of
> nonempty local open conditions on the displayed rational surface can be
> imposed simultaneously.

## 5. A rational `D_5` curve

A one-parameter specialization of Brumer's dihedral family has an especially
small affine lift into the fixed target.  Let `r\ne0` and put

\[
 \lambda=20r^3.
\]

Start with

\[
\begin{aligned}
G_\lambda(X)={}&X^5-5X^4+5(2-\lambda^2)X^3\\
&+5(2\lambda^2-1)X^2+5(1-\lambda^2)X-2.
\end{aligned}
\]

If `\theta` denotes the class of `X`, the trace-zero scaled generator

\[
 \eta=\frac{\theta-1}{\lambda}
\]

has characteristic polynomial

\[
\boxed{
\begin{aligned}
f_r(T)={}&T^5-5T^3
-\frac{400r^6-1}{1600r^9}T^2\\
&+\frac1{16000r^{12}}T+\frac1{800000r^{15}}.
\end{aligned}}
\]

It is the normalized inverse polynomial of the fixed map at

\[
\boxed{
(\Pi,B,C)=
\left(
\frac1{40r^4},
\frac{400r^6-1}{80r^5},
-64r^5
\right).
}
\]

The discriminant is the rational square

\[
\boxed{
\operatorname{Disc}(f_r)=
\left(
\frac{
256000000r^{18}+160000r^{12}-800r^6-7
}{
40960000000r^{30}
}
\right)^2.
}
\]

The pair-sum resolvent factors identically.  More precisely, if

\[
\begin{aligned}
U_r(X)={}&3200000r^{15}X^5
-16000000r^{15}X^3-40000r^9X^3\\
&-800000r^{12}X^2-2000r^6X^2
+40000r^9X+2000r^6-7,
\end{aligned}
\]

and

\[
\begin{aligned}
V_r(X)={}&1600000r^{15}X^5
-16000000r^{15}X^3+20000r^9X^3\\
&+2000r^6X^2+40000000r^{15}X
-120000r^9X-50r^3X\\
&-6000r^6+19,
\end{aligned}
\]

then

\[
 \boxed{R_{10}^{f_r}(X)=
 \frac{U_r(X)V_r(X)}{5120000000000r^{30}}.}
\]

At `r=-1/2`, the target becomes

\[
 \left(\frac25,-\frac{21}{10},2\right),
\]

the compact `D_5` witness from Section 6.2 of the main note.  Its quintic is
irreducible modulo `3`, both pair-resolvent factors are irreducible modulo
`3`, and its factor degrees modulo `11` are `(2,2,1)`.  Hence the generic
group on this rational curve is exactly `D_5`, rather than `C_5`.

### 5.1 Signatures and unramified local types on the curve

Exact Sturm counts give five real roots at `r=1/2` and one real root at
`r=1/4`.  The corresponding real chambers are open, so Hilbert
specialization supplies infinitely many `D_5` fields of both possible
signatures.

At `r=-1/2`, the primitive integral polynomial

\[
3125T^5-15625T^3+5250T^2+800T-128
\]

has the three possible unramified `D_5` types

\[
\begin{array}{c|c}
p&\text{factor degrees}\\ \hline
3&(5)\\
11&(2,2,1)\\
23&(1,1,1,1,1).
\end{array}
\]

Their residue neighborhoods give simultaneous local conditions on the
`r`-line.

### 5.2 Pairwise nonisomorphic `D_5` fields

The unscaled Brumer discriminant is

\[
 \operatorname{Disc}(G_\lambda)
 =
 62500\bigl(4\lambda^6+\lambda^4-2\lambda^2-7\bigr)^2.
\]

Put

\[
 H(r)=4(20r^3)^6+(20r^3)^4-2(20r^3)^2-7.
\]

Schur's theorem gives infinitely many primes dividing values of the
nonconstant integral polynomial `H`.  Outside the finite set dividing its
leading coefficient or discriminant, such a prime supplies a nonzero simple
root `r_0` modulo `p`.  Lift it to `\rho\in\mathbb Z_p` and perturb inside
`\rho+p\mathbb Z_p` so that

\[
 v_p(H(r))=1.
\]

At this branch the quintic has two double geometric roots and one simple
root.  For `p\nmid10`, tame inertia is therefore a reflection of `D_5`,
with cycle type `(2,2,1)`, so the quintic field is ramified at `p`.

As in Sections 4.2 and 6, impose one such new ramified neighborhood at each
Hilbert step, choosing `p` unramified in all previously obtained fields.
This proves:

> **Theorem.**  The fixed Keller map has infinitely many pairwise
> nonisomorphic connected full `D_5` fibers with five real roots, and
> infinitely many with one real root.  Any fixed finite collection of
> nonempty local open conditions on the displayed rational curve can be
> imposed simultaneously.

## 6. The completed `S_5` infinitude statement

The existing fixed-map Hilbert theorem allows any nonempty real open and
finitely many nonempty `p`-adic opens in `U`.  Its stated infinitude concerns
targets.  The following fresh-ramification argument upgrades it to fields.

Consider the rational point on the discriminant

\[
 y_0=\left(1,-\frac32,\frac32\right).
\]

At this point

\[
 f_{y_0}(T)=(T-1)^2(T+1)(T^2+T-3),
\]

and

\[
 \nabla\Delta(y_0)=(0,-468,-468).
\]

Thus for every prime `p\notin\{2,3,13\}`,

\[
 y_p=\left(1,-\frac32,\frac32+p\right)
\]

satisfies

\[
 v_p(\Delta(y_p))=1.
\]

The same equality holds throughout a sufficiently small `p`-adic
neighborhood `\Omega_p\subset U(\mathbb Q_p)` of `y_p`.

Fix once and for all an allowed real chamber and any finite compatible
package of local open conditions.  Inductively, after fields
`K_1,\ldots,K_n` have been selected, choose a prime
`p_{n+1}\notin\{2,3,13\}` outside the fixed local package and unramified in
all previous fields.  Add `\Omega_{p_{n+1}}` to the local conditions and
apply the fixed-map Hilbert theorem.  It produces an `S_5` quintic field
`K_{n+1}` whose defining monic polynomial is `p_{n+1}`-integral and has
discriminant
valuation one at `p_{n+1}`.  The index-discriminant formula then gives

\[
 v_{p_{n+1}}(\operatorname{Disc}K_{n+1})=1.
\]

Hence `p_{n+1}` ramifies in `K_{n+1}` and in none of the earlier fields.
The fields are pairwise nonisomorphic.

Therefore:

> **Theorem.**  One fixed explicit Keller map has infinitely many pairwise
> nonisomorphic connected full `S_5` fibers in every quintic real signature.
> Any fixed finite collection of nonempty local open conditions can be
> imposed simultaneously.  In particular, any finite collection of
> available unramified splitting conditions may be imposed.

This proof does not use the density-one theorem for all quintic fields.
That theorem concerns ordering abstract fields by field discriminant, whereas
target-height ordering in this fixed three-parameter family has different
multiplicities and a different discriminant map.

## 7. Real signatures and unramified local types

Complex conjugation has cycle type

\[
\begin{array}{c|c}
\text{number of real roots} & \text{cycle type}\\ \hline
5&(1,1,1,1,1)\\
3&(2,1,1,1)\\
1&(2,2,1).
\end{array}
\]

Intersecting these types with each group gives the only possible
signatures.  The full cycle-type set gives the possible unramified local
factorizations.

| group | possible real-root counts | possible unramified factor degrees |
|---|---|---|
| `C_5` | `5` | `(1,1,1,1,1)`, `(5)` |
| `D_5` | `5,1` | `(1,1,1,1,1)`, `(5)`, `(2,2,1)` |
| `F_{20}` | `5,1` | `(1,1,1,1,1)`, `(5)`, `(2,2,1)`, `(4,1)` |
| `A_5` | `5,1` | `(1,1,1,1,1)`, `(5)`, `(3,1,1)`, `(2,2,1)` |
| `S_5` | `5,3,1` | all seven partitions of five |

This table is group-theoretic.  In the fixed Keller map, all three `S_5`
signatures and both signatures for each of `F_{20}` and `D_5` are proved.
Both `A_5` signatures occur in the fixed Keller map; exact targets are
given in Section 8.1.  A `C_5` quintic is necessarily totally real, and
the compact cyclic anchor realizes both of its unramified cycle types:
factor degrees `(5)` modulo `2` and `(1,1,1,1,1)` modulo `7`.

For an explicit regular `G`-family, a residue class with squarefree reduction
of one of the listed types gives a `p`-adic open imposing that Frobenius
class.  Ramified extensions require separate inertia and decomposition-group
analysis and are not implied by the table.

## 8. What a successful exceptional parameterization must prove

For either `G=A_5,C_5`, it is enough to exhibit a rational
curve or surface `V` and a rational map

\[
 V\dashrightarrow\mathcal I_G
\]

such that:

1. the induced quintic over `Q(V)` is irreducible and has group exactly `G`;
2. the map meets the squarefree full-fiber locus;
3. every desired real chamber and unramified Frobenius class has a local
   point on `V`;
4. `V` has weak approximation after removal of the bad divisors; and
5. the family meets a simple discriminant branch at infinitely many good
   primes.

Hilbert irreducibility with local conditions then gives infinitely many
specializations with group `G`.  The fresh-ramification induction of
Section 6 gives pairwise nonisomorphic fields.  Conditions 3 and 4 give the
signature and local refinements.

### 8.1 Exact `A_5` anchors in both real chambers

The target

\[
 y_+=(\Pi,B,C)=\left(1,0,-\frac25\right)
\]

has inverse polynomial

\[
 P_+(T)=T^5-5T^3+4T+\frac45.
\]

Its integral multiple `5P_+` has discriminant `145000^2`, has five real
roots, and has the following good-prime factor degrees:

\[
\begin{array}{c|cccc}
p&3&17&23&211\\ \hline
\text{degrees}&(5)&(2,2,1)&(3,1,1)&(1,1,1,1,1).
\end{array}
\]

The target

\[
 y_-=(\Pi,B,C)=\left(2,4,-\frac14\right)
\]

has inverse polynomial

\[
 P_-(T)=T^5-5T^3-16T^2+32T+16.
\]

It has discriminant `54176^2`, has one real root, and has good-prime
factor degrees

\[
\begin{array}{c|cccc}
p&3&7&13&389\\ \hline
\text{degrees}&(5)&(3,1,1)&(2,2,1)&(1,1,1,1,1).
\end{array}
\]

In each row the `(5)` factorization proves irreducibility.  The square
discriminant puts the group in `A_5`, while the `(3,1,1)` factorization
excludes the transitive proper subgroups `C_5` and `D_5`.  Thus both
targets have group exactly `A_5`.  Together the rows prove that both
allowed signatures, and every unramified `A_5` cycle type, actually occur
in the fixed map.  This is an attainment statement at exact anchors, not
yet an infinitude or simultaneous-local-conditions theorem.

### 8.2 Two Mestre source pencils and the affine-lift obstruction

The totally-real anchor lies on an especially small square-discriminant
pencil.  Put

\[
\begin{aligned}
P(X)&=X^5-5X^3+4X+\frac45,\\
Q(X)&=-X^4+X^3+2X^2-X-\frac{21}{25},\\
R(X)&=X^4-X^3-X^2+X+\frac85.
\end{aligned}
\]

Then the exact Mestre identity

\[
 \boxed{PQ'-P'Q=R^2}
\]

gives

\[
\operatorname{Disc}_X(P-tQ)
=\frac4{15625}
\left(
427t^4+2335t^3+7925t^2+12125t+14500
\right)^2.
\]

The specialization at `t=0` is the `A_5` polynomial above.  Specialization
therefore forces the generic group over `Q(t)` to contain `A_5`, while the
square discriminant puts it inside `A_5`; hence the generic group is
exactly `A_5`.  A Sturm calculation shows that the displayed quartic has
no real zero.  The number of real roots is consequently constant on the
real parameter line, so every member has five real roots.

There is also an exact one-real-root pencil through `P_-`.  Define

\[
\begin{aligned}
Q_-(X)={}&-X^4+\frac{9499}{1790}X^3
-\frac{36290687}{3204100}X^2
+\frac{23868087}{1602050}X
-\frac{37781321}{3204100},\\
R_-(X)={}&X^4-\frac{9499}{1790}X^3
+\frac{4841}{895}X^2
+\frac{26683}{1790}X
-\frac{22208}{895}.
\end{aligned}
\]

The verifier checks

\[
 P_-Q_-'-P_-'Q_-=R_-^2
\]

and

\[
 \operatorname{Disc}_X(P_--tQ_-)
 =
 \left(\frac{1693\,A_-(t)}{16447056722460500000}\right)^2,
\]

where

\[
\begin{aligned}
A_-(t)={}&599061012975492710t^4
+18908224953293263933t^3\\
&+193635834698041154800t^2
+591172039750728710000t\\
&-526305815118736000000.
\end{aligned}
\]

Again the generic group of `P_--tQ_-` is `A_5`.  Its discriminant has no
zero on `[-1/2,1/2]`, and the member at zero has one real root, so this
interval is a one-real-root chamber.  These are regular `A_5` source
polynomials; they are not yet curves in the three-dimensional fixed-map
target.

For the first pencil the direct affine lift can be analyzed completely.
Center by `X=Z-t/5`.  The cubic and linear coefficients become

\[
 -\frac{2t^2+5t+25}{5},\qquad
 -\frac{3t^4+15t^3-25t^2-125t-500}{125}.
\]

If the centered root is rescaled affinely into the fixed normalized chart,
the necessary and sufficient square and cube conditions are

\[
\lambda^2=\frac{2t^2+5t+25}{25},\qquad
\Pi^3=
-\frac{5(3t^4+15t^3-25t^2-125t-500)}
{4(2t^2+5t+25)^2}.
\]

The conic has the rational parametrization

\[
t=\frac{5(1-2r)}{r^2-2},\qquad
5\lambda=-\frac{5(r^2-r+2)}{r^2-2}.
\]

After substitution, the remaining cover is

\[
\boxed{
4(r^2-r+2)^4\Pi^3=N_8(r),
}
\]

where

\[
\begin{aligned}
N_8(r)={}&4r^8-10r^7-7r^6+160r^5-429r^4\\
&+290r^3-23r^2-60r+59.
\end{aligned}
\]

The polynomial `N_8` is squarefree and coprime to `r^2-r+2`.
The degree-three map from the smooth projective normalization of this
curve to the `r`-line is totally ramified over the eight zeros of `N_8`
and the two zeros of `r^2-r+2`, and is unramified at infinity.
Riemann--Hurwitz gives

\[
 2g-2=3(-2)+10(3-1)=14,\qquad \boxed{g=8}.
\]

By Faltings' theorem this curve has only finitely many rational points.
Thus this Mestre pencil cannot prove infinitely many fixed-map `A_5`
fibers by affine changes of generator.  The compact lift is the rational
point `r=1/2,\Pi=1`.  A successful use of either source pencil must employ
a genuinely non-affine generator.

The first quadratic non-affine ansatz can also be excluded completely.
Let `z` be the centered root above and put

\[
 \eta=z+h\left(z^2-\frac{\operatorname{Tr}(z^2)}5\right).
\]

Then `\operatorname{Tr}(\eta)=0` and

\[
 \operatorname{Tr}(\eta^2)=\frac{2Q_h(t,h)}{125},
\]

where

\[
\begin{aligned}
Q_h={}&(18t^4+90t^3+325t^2+500t+875)h^2\\
&-(60t^3+225t^2+375t)h+50t^2+125t+625.
\end{aligned}
\]

Rational rescaling to second trace `10` requires `Q_h` to be a square in
`\mathbb Q(t)`.  In fact the conic

\[
 w^2=Q_h(t,h)
\]

has no `\mathbb Q(t)`-point.  Here is a local proof.  Set `u=1/t`.  A
comparison of leading `u`-adic terms shows first that any solution must
have `h=uH` with `H(0)=5/3`: if `v_u(h)<1` the leading square class is
`18`, if `v_u(h)>1` it is `50`, and if `v_u(h)=1` the initial form is
`2(3H(0)-5)^2`.

Write `H=5/3+K`.  After multiplying by `u^2`, the first terms are

\[
\begin{aligned}
u^2Q_h(1/u,uH)
={}&18K^2+(90K^2+75K)u\\
&+\left(325K^2+\frac{2125}{3}K+\frac{8125}{9}\right)u^2+\cdots.
\end{aligned}
\]

If `v_u(K)>1`, the leading square class is
`8125/9=(25/3)^2\cdot13`, impossible.  If `v_u(K)=1`, write
`K=cu+\cdots`.  The leading coefficient would have to satisfy

\[
 d^2=18c^2+75c+\frac{8125}{9}.
\]

Equivalently, with `X=36c+75` and `Y=6d`,

\[
 X^2-2Y^2=-59375=-5^5\cdot19.
\]

But `2` is nonsquare modulo `5`, so
`\mathbb Q_5(\sqrt2)/\mathbb Q_5` is unramified and every norm has even
5-adic valuation.  The right side has valuation five, a contradiction.
Thus even this quadratic generator cannot meet the fixed trace
normalization over `\mathbb Q(t)`.  Higher-degree non-affine generators,
or a direct curve on `\mathcal I_A`, remain open.

### 8.3 A smaller explicit `A_5` descent surface

There is a substantially smaller model than the full oriented-discriminant
threefold.  Put

\[
 a=5k^2-1,\qquad
 g_k(X)=X^5+5aX+4a.
\]

Its discriminant is identically a square:

\[
 \operatorname{Disc}(g_k)
 =4\,000\,000\,k^2a^4=(2000ka^2)^2.
\]

At `k=1`, reduction modulo `3` is irreducible and reduction modulo `7`
has factor degrees `(3,1,1)`.  Thus the generic group over `Q(k)` is
`A_5`.  If `\theta` is the class of `X`, take the sparse generator

\[
 \eta=u\theta+v\theta^3.
\]

Exact power traces give

\[
 \operatorname{Tr}(\eta)=0,\qquad
 \operatorname{Tr}(\eta^2)=-40auv.
\]

Consequently the first fixed-map normalization is simply

\[
 \boxed{auv=-\frac14.}
\]

The characteristic polynomial is

\[
\begin{aligned}
\chi_\eta(T)={}&T^5+20auvT^3
{}+20av(u^2-3av^2)T^2\\
&+5a(25a^2v^4+10au^2v^2-16auv^3+u^4)T\\
&+4a(25a^2uv^4+16a^2v^5+10au^3v^2+u^5).
\end{aligned}
\]

The remaining cube descent is therefore the single equation

\[
 \boxed{
 5a(25a^2v^4+10au^2v^2-16auv^3+u^4)=4\Pi^3.
 }
\]

Together with `a=5k^2-1` and `auv=-1/4`, this defines an explicit surface
`\mathcal S_A`.  On its primitive irreducible open it maps to
`\mathcal I_A` by

\[
\begin{aligned}
B&=-\frac{10av(u^2-3av^2)}{\Pi},\\
C&=-\frac{
 2a(25a^2uv^4+16a^2v^5+10au^3v^2+u^5)
}{\Pi^5}.
\end{aligned}
\]

The normalization equation can be eliminated without losing any nonzero
solutions.  Write

\[
 u=\frac{s}{2a},\qquad v=-\frac1{2s}.
\]

Then `auv=-1/4` identically, and the remaining surface is the single
cyclic-cover equation

\[
\boxed{
\Pi^3=
\frac{
5(25a^6+16a^4s^2+10a^3s^4+s^8)
}{64a^3s^4},
\qquad a=5k^2-1.
}
\]

Thus the alternating gap is now a rational-curve problem on an explicit
two-parameter cubic cover.  This smaller equation is also suitable for
local-obstruction and low-degree-section searches.

There is an exact first degree obstruction.  Suppose `s=s(k)` is a
polynomial and `\Pi\in\mathbb Q(k)`.  For `\deg s=1` or `2`, the rational
function on the right side of the boxed equation has order `2` at
`k=\infty`.  It therefore cannot be a cube in `\mathbb Q(k)`.
More generally, for `\deg s=d\ge2` its order at infinity is `4d-6`;
hence `d` must be divisible by `3`.

The constant case is also impossible.  Put `s=c\ne0` and

\[
 P_c(k)=25a^6+16a^4c^2+10a^3c^4+c^8.
\]

If the right side were a rational cube, then

\[
 W=4ac^2\Pi,\qquad W^3=5c^2P_c(k).
\]

Since the right side is a polynomial, `W` must be a polynomial.  It is
even and has degree four, so write `W=A k^4+B k^2+C`.  Exact coefficient
elimination in

\[
 (Ak^4+Bk^2+C)^3-5c^2P_c(k)=0
\]

puts `c^8` in the coefficient ideal.  Thus `c=0`, a contradiction.

For a general nonconstant polynomial `s`, clear denominators by putting

\[
 W=4as^2\Pi,\qquad
 W^3=5s^2(25a^6+16a^4s^2+10a^3s^4+s^8).
\]

If an irreducible factor `q` of `s` does not divide `a`, reduction modulo
`q` leaves `25a^6` in the final factor.  Hence `2v_q(s)` must be divisible
by three.  The same conclusion holds for `q=a=5k^2-1`.  Indeed, if
`e=v_a(s)>0`, the four terms in the final factor have valuations

\[
 6,\quad 4+2e,\quad 3+4e,\quad 8e.
\]

For `e=1`, the two initial terms have residue
`25+16(s/a)^2`, which is nonzero in the real quadratic field
`\mathbb Q(\sqrt5)`; for `e\ge2`, the first term is minimal.  The total
valuation is therefore divisible by three only when `e` is divisible by
three.  Consequently every irreducible factor of `s` has multiplicity
divisible by three:

\[
 \boxed{s=c\,h(k)^3.}
\]

After dividing `W` by `h^2`, the cube identity becomes

\[
 V^3=5c^2(25a^6+16a^4s^2+10a^3s^4+s^8).
\]

If `\deg s=3`, then `h=k-r` after absorbing its leading coefficient into
`c`.  Evaluation at `k=r` gives

\[
 V(r)^3=125c^2a(r)^6=(5a(r)^2)^3c^2,
\]

so `v_5(c)\equiv0\pmod3`.  If `\deg s=6`, then `h` has degree two.
For any irreducible factor `q` of `h`, the same evaluation, after dividing
by `a^6` when `q=a`, shows that `c^2` is a cube in the number field
`\mathbb Q[k]/(q)`.  Every ramification index above `5` in this extension
of degree at most two is `1` or `2`; hence again
`v_5(c)\equiv0\pmod3`.

In both degrees, however, the leading coefficient at infinity on the
right is `5c^{10}` times an evident rational cube.  It would have to be a
cube, forcing

\[
 1+10v_5(c)\equiv0\pmod3,
\]

a contradiction.  Consequently:

\[
\boxed{\text{a polynomial-section search must begin at }\deg s=9.}
\]

At degree nine the argument leaves only `s=c h^3` with `h` cubic, and a
necessary condition is that every prime of `\mathbb Q[k]/(h)` above `5`
have ramification index divisible by three.  Thus `5` must be totally
ramified in the cubic algebra.  This does not exclude rational,
nonpolynomial choices of `s`.

Thus rationally parametrizing `\mathcal S_A`, with a nonconstant `k`, is
now the precise alternating-group gap.  The equations above are an exact
reduction, not a claim that `\mathcal S_A` is rational.

### 8.4 A trace-normalized non-affine `C_5` section

Let `\theta` satisfy Emma Lehmer's cyclic quintic

\[
\begin{aligned}
L_n(X)={}&X^5+n^2X^4-(2n^3+6n^2+10n+10)X^3\\
&+(n^4+5n^3+11n^2+15n+5)X^2\\
&+(n^3+4n^2+10n+10)X+1.
\end{aligned}
\]

Its cyclicity and the automorphism below are classical; see
[Darmon, *Note on a polynomial of Emma
Lehmer*](https://www.math.mcgill.ca/darmon/pub/Articles/Research/03.Lehmer/paper.pdf).
Put

\[
 Q=n^4+5n^3+15n^2+25n+25.
\]

If `r_0=\theta` and

\[
 r_{i+1}=
 \frac{(n+2)+nr_i-r_i^2}{1+(n+2)r_i},
\]

then `r_0,\ldots,r_4` are cyclically ordered roots of `L_n`.  Define

\[
\begin{array}{c|rrrrr}
i&0&1&2&3&4\\ \hline
A_i&2&1&0&-1&-2\\
B_i&5/2&5/2&5&-5/2&-15/2\\
C_i&0&5&10&-10&-5
\end{array}
\]

and

\[
 d_i=A_in^2+B_in+C_i,\qquad
 \boxed{\eta_n=\frac1Q\sum_{i=0}^4d_ir_i.}
\]

The elementary identities

\[
 \sum_i d_i=0,\qquad \sum_i d_i^2=10Q
\]

combine with the cyclic trace-pairing identity

\[
 \operatorname{Tr}\left(\left(\sum_i c_ir_i\right)^2\right)
 =Q\sum_i c_i^2\qquad\left(\sum_i c_i=0\right)
\]

to give

\[
 \boxed{\operatorname{Tr}(\eta_n)=0,\qquad
 \operatorname{Tr}(\eta_n^2)=10.}
\]

Thus a genuinely non-affine normal-basis generator removes the square
obstruction over the entire Lehmer parameter line.  More generally, the
coefficient space

\[
 \sum_i c_i=0,\qquad Q\sum_i c_i^2=10
\]

is a three-dimensional quadric over `Q(n)`.  The displayed coefficients
give it a rational point, so this quadric is rational over `Q(n)`.  Its
The rationality can be made completely explicit.  For

\[
 z=(z_0,\ldots,z_4),\qquad \sum_i z_i=0,\qquad
 Z_2=\sum_i z_i^2\ne0,
\]

put `D_z=\sum_i d_i z_i` and

\[
 \boxed{
 c_i(z)=\frac{d_i}{Q}-\frac{2D_z}{QZ_2}z_i.
 }
\]

Then

\[
 \sum_i c_i(z)=0,\qquad Q\sum_i c_i(z)^2=10.
\]

This is the second-intersection parametrization by lines through the
known point `(d_i/Q)`.  Conversely, for any other point `c` on the
quadric, taking `z=c-(d_i/Q)` recovers `c`; hence it is birational from
`\mathbb P^3_{\mathbb Q(n)}`.

The full remaining cyclic incidence is therefore the explicit
three-dimensional cubic cover

\[
\boxed{
\Pi^3=
\frac{50-\operatorname{Tr}\left(
  \left(\sum_i c_i(z)r_i\right)^4
\right)}{16}
}
\]

over `\mathbb A^1_n\times\mathbb P^3_z`.  The displayed one-parameter
section is obtained by retaining the original point.  Its fourth trace
gives

\[
 \frac{50-\operatorname{Tr}(\eta_n^4)}{16}
 =\frac{N(n)}{64Q^3},
\]

where

\[
\begin{aligned}
N(n)={}&64n^{12}+960n^{11}+6780n^{10}+31800n^9
 +119875n^8+385875n^7\\
&+1055000n^6+2370000n^5+4206250n^4+5600000n^3\\
&+5237500n^2+3125000n+1000000.
\end{aligned}
\]

Consequently this section enters the fixed map exactly on the explicit
cubic cover

\[
 \boxed{w^3=N(n),\qquad \Pi=\frac{w}{4Q}.}
\]

At `n=0`, `N(0)=100^3`, `\Pi=1`, and the characteristic polynomial is

\[
 T^5-5T^3+4T+\frac75,
\]

recovering the compact cyclic target `(1,0,-7/10)`.  This is a nonconstant
trace-normalized section through the known point, but the displayed
one-parameter cubic cover is not rational.  Indeed, `N` is squarefree of
degree twelve.  The degree-three map from the smooth projective
normalization of `w^3=N(n)` to the `n`-line is totally ramified at the
twelve finite zeros of `N` and unramified at infinity.  Therefore

\[
 2g-2=3(-2)+12(3-1)=18,\qquad \boxed{g=10}.
\]

By Faltings' theorem this particular section has only finitely many
rational points.  It cannot by itself prove cyclic infinitude.  A rational
curve or surface must instead use the remaining three rational parameters
on the trace quadric before imposing the cube condition.

For completeness, the integral polynomial at the compact point is

\[
 5T^5-25T^3+20T+7.
\]

It has discriminant `26875^2`, five real roots, is irreducible modulo `2`,
and splits completely modulo `7`.  Thus both possible unramified `C_5`
local types are attained in the fixed map.

For comparison, merely centering and affinely rescaling the original
defining generator leaves both obstructions visible.  Namely,

\[
 \operatorname{Tr}(\eta^2)=\frac{4Q}{5}.
\]

An affine rescaling into the fixed normalized chart would have to satisfy

\[
 \lambda^2=\frac{2Q}{25},\qquad
 \Pi^3=
 -\frac{5(3n^4+15n^3+20n^2-50)}{16Q}.
\]

These independent square and cube conditions explain why the direct affine
presentation does not give a rational parameterization without solving
further covers.  The explicit `\eta_n` above removes the square cover and
isolates the remaining cyclic calculation.

The present state is therefore:

| group | pullback condition | rational family with generic group | pairwise-nonisomorphic theorem |
|---|---|---|---|
| `S_5` | complete | complete | complete |
| `F_{20}` | complete | complete | complete |
| `D_5` | complete | complete | complete |
| `A_5` | complete | open; explicit cubic surface, with affine, quadratic, and degree-at-most-six polynomial obstructions | open |
| `C_5` | complete | open; displayed section has genus ten | open |

The next `A_5` computation should use a non-affine generator over one of
the Mestre pencils, or search directly on `\mathcal I_A`; the affine route
through the totally-real pencil is now closed.  The oriented-discriminant
equation is smaller than the coefficient-factorization and automorphism
incidences, and the two real-signature anchors provide bases for searches
in different real chambers.  The next `C_5` computation must vary the
remaining rational parameters on the trace quadric: the displayed
one-parameter section is now proved incapable of giving infinitude.

## 9. Arithmetic-statistics experiment

A bounded target census should keep presentation multiplicity separate from
field isomorphism.  For primitive projective targets `[W:P:B_0:C_0]` of
height at most `H`, the reproducible pipeline is:

1. apply one declared target symmetry convention and record the number of
   raw and retained presentations;
2. discard `\Pi=0`, zero discriminant, and reducible quintics;
3. compute the exact Galois group with PARI/GP;
4. record field discriminant and signature;
5. use `polredabs` only as a fast reduced-polynomial key, not as a proof of
   isomorphism;
6. within equal degree, signature, and field-discriminant buckets, use
   `nfisisom` to form exact isomorphism classes;
7. record good-prime factorization types and exclude primes dividing the
   field discriminant from unramified equidistribution counts.

The primary outputs should be both

\[
 \frac{\#\{\text{isomorphism classes}\}}
      {\#\{\text{connected target presentations}\}}
 \quad\text{and}\quad
 \frac{\#\{\text{presentations of group }G\}}
      {\#\{\text{connected target presentations}\}}.
\]

Neither ratio should be described as a density of quintic fields ordered by
discriminant.  Report results simultaneously by target height, coefficient
height, and field discriminant; otherwise duplicate presentations can
dominate the apparent group frequencies.
