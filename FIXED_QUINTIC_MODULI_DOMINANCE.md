# One fixed quintic Keller map dominating affine quintic moduli

## Status

This note proves a geometric fixed-seed strengthening of the
finite-etale-fiber construction and gives an explicit arithmetic zoo inside
that same map.  It does **not** assert that every rank-five finite etale
algebra over `Q` occurs at a rational target of the fixed map.  The
distinction between geometric moduli dominance and arithmetic surjectivity
is recorded explicitly in Section 5.

The construction uses the root-engineered quadratic gauge and its
scheme-theoretic fiber reconstruction from
[`papers/common-arithmetic-fibers/sections/02-quadratic-gauge.tex`](papers/common-arithmetic-fibers/sections/02-quadratic-gauge.tex).

## 1. The fixed seed and fixed Keller map

Take the split squarefree seed

\[
 G(S)=S(S^2-1)(S^2-4)=S^5-5S^3+4S.
\]

Thus

\[
 (g_1,g_2,g_3,g_4,g_5)=(4,0,-5,0,1),
 \qquad g_1g_3g_5\ne0.
\]

Put

\[
 t=1+xy,\qquad
 q=t^2z-\frac45y^2(1+3t),
\]

and define

\[
 \boxed{
 F(x,y,z)=
 \left(
 tq,\;
 y-\frac{15}{4}xq+\frac54t^2x^3q^5,\;
 x(5-3t)+\frac54x^3z-\frac34(xq)^5
 \right).
 }
\]

This is the specialization of the general quadratic-gauge formula.  Hence

\[
 \boxed{\det DF=-2,\qquad \operatorname{gdeg}(F)=5.}
\]

The direct symbolic expansion is included in
[`scripts/verify_fixed_quintic_moduli_dominance.py`](scripts/verify_fixed_quintic_moduli_dominance.py).
The coordinate-degree profile is `(7,32,30)`.

Write the target coordinates as `(\Pi,B,C)`.  The inverse polynomial is

\[
 \boxed{
 E_{\Pi,B,C}(S)
 =\Pi^5S^5-5\Pi S^3-2BS^2+4S-2C.
 }
\]

For `\pi\ne0` and squarefree `E_{\pi,b,c}`, the existing localized-fiber
theorem gives the full scheme fiber

\[
 F^{-1}(\pi,b,c)
 \simeq
 \operatorname{Spec}\mathbb Q[S]/(E_{\pi,b,c}).
\]

At `(1,0,0)` the inverse polynomial is the seed itself, so the fixed map has
the complete split fiber with roots `0,\pm1,\pm2`.

## 2. Tiny coefficient-Jacobian certificate

On `\Pi\ne0`, make the affine generator change

\[
 S=\Pi^{-2}T
\]

and multiply the inverse equation by `\Pi^5`.  Then

\[
 \Pi^5E_{\Pi,B,C}(\Pi^{-2}T)
 =
 T^5-5T^3+c_2T^2+c_3T+c_4,
\]

where

\[
 \boxed{
 c_2=-2\Pi B,\qquad
 c_3=4\Pi^3,\qquad
 c_4=-2\Pi^5C.
 }
\]

Consequently

\[
 \boxed{
 \det
 \frac{\partial(c_2,c_3,c_4)}
      {\partial(\Pi,B,C)}
 =-48\Pi^8.
 }
\]

At the rational point

\[
 (\Pi,B,C)=\left(1,\frac52,0\right)
\]

the determinant is `-48`.  The corresponding inverse polynomial is

\[
 S^5-5S^3-5S^2+4S,
\]

whose discriminant is `-1139056`, so the certificate lies in the full
finite-etale locus.

## 3. Geometric moduli consequence

Let `M_5^{aff}` denote the moduli of unordered squarefree five-point
configurations on the affine line modulo affine changes of coordinate.
Equivalently, this is the moduli of squarefree quintic presentations modulo
nonzero polynomial scaling and changes `T\mapsto aT+b`.

On the open set where the centered cubic coefficient is nonzero, every
geometric orbit has a representative

\[
 T^5-5T^3+c_2T^2+c_3T+c_4.
\]

It is unique up to `T\mapsto-T`, which sends

\[
 (c_2,c_3,c_4)\longmapsto(-c_2,c_3,-c_4).
\]

Thus this coefficient cross-section is a finite dominant double cover of
the corresponding open of `M_5^{aff}`.  The nonzero Jacobian above proves
that the target of the single fixed map dominates this cross-section and
hence dominates `M_5^{aff}`.

More explicitly, over an algebraically closed field and on `c_3\ne0`, choose

\[
 \Pi^3=\frac{c_3}{4},\qquad
 B=-\frac{c_2}{2\Pi},\qquad
 C=-\frac{c_4}{2\Pi^5}.
\]

After also removing the quintic discriminant, every resulting inverse
polynomial is squarefree and therefore gives a full degree-five etale fiber.
This proves:

> **Fixed-quintic geometric dominance.**
> One explicit Keller map of `A^3_Q`, of geometric degree five, contains
> among its full fibers a family dominating the moduli of affine quintic
> root configurations.

## 4. A rational quotient chart

The preceding normalization is rational along this Keller family, but a
generic quintic reaches the fixed cubic coefficient only after extracting a
square root.  A rational quotient chart makes the arithmetic descent issue
visible.

After dividing `E_{\Pi,B,C}` by `\Pi^5`, write the centered monic quintic as

\[
 S^5+aS^3+bS^2+cS+d
\]

with

\[
 a=-\frac5{\Pi^4},\quad
 b=-\frac{2B}{\Pi^5},\quad
 c=\frac4{\Pi^5},\quad
 d=-\frac{2C}{\Pi^5}.
\]

On `ab\ne0`, scaling by `\lambda=b/a` gives the unique rational
normalization

\[
 T^5+uT^3+uT^2+vT+w,
\]

where

\[
 \boxed{
 u=-\frac{125}{4\Pi^2B^2},\qquad
 v=\frac{625}{4\Pi B^4},\qquad
 w=-\frac{3125C}{16B^5}.
 }
\]

Its Jacobian is

\[
 \det
 \frac{\partial(u,v,w)}
      {\partial(\Pi,B,C)}
 =
 \frac{732421875}{128\Pi^4B^{12}}.
\]

At `(1,5/2,0)` it equals `96`.  This independently proves dominance on a
rational quotient chart without using a square-root normalization.

## 5. Arithmetic limitation

Geometric dominance does not imply that every rational moduli point has a
rational target above it.  From the small coefficient certificate,

\[
 \Pi^3=\frac{c_3}{4}.
\]

From the rational quotient coordinates one similarly obtains

\[
 \Pi^3=\frac{25}{4}\frac{v}{u^2},
 \qquad
 B^2=-\frac{125}{4\Pi^2u}.
\]

Thus the generic rational lift requires a cubic and then a quadratic
extraction.  The map to rational affine-quintic moduli is generically
six-to-one.  Its rational image is Zariski dense but is a thin subset in the
sense relevant to Hilbert irreducibility.

There is a second distinction.  `M_5^{aff}` parametrizes a finite etale
algebra together with a primitive generator, modulo affine changes of that
generator.  Two primitive generators of the same abstract algebra need not
be affinely related.  Therefore the calculation does not rule out the
possibility that every abstract quintic algebra occurs by some other
primitive generator, but it does not prove it either.

The proved statement should consequently be described as geometric
dominance of affine quintic-presentation moduli, not as rational
surjectivity onto all quintic algebras.

## 6. Arithmetic range inside the same map

The fixed map already contains irreducible `S_5` fibers of all three
possible quintic signatures:

| target `(\Pi,B,C)` | inverse polynomial | real roots | discriminant |
|---|---|---:|---:|
| `(1,-1,-1)` | `S^5-5S^3+2S^2+4S+2` | 1 | `500624` |
| `(1,0,-1)` | `S^5-5S^3+4S+2` | 3 | `-57056` |
| `(1,0,-1/2)` | `S^5-5S^3+4S+1` | 5 | `38569` |

For each row, exact Sturm counting gives the displayed number of real roots.
The checker proves irreducibility and the presence of cycle types `(5)`,
`(4,1)`, and `(3,2)` at good primes.  The transitive subgroup classification
in degree five then gives Galois group `S_5`.

The general symmetric-monodromy theorem for quadratic-gauge seeds and Hilbert
irreducibility already imply infinitely many connected `S_5` fibers of this
fixed map.  Combining Hilbert irreducibility with real local conditions gives
infinitely many in each of the three signature regions above.

### 6.1 One explicit arithmetic zoo

The same fixed map has five complete rational target fibers with sharply
different arithmetic:

| type | target `(\Pi,B,C)` | inverse polynomial |
|---|---|---|
| `Q^5` | `(1,0,0)` | `S(S-1)(S+1)(S-2)(S+2)` |
| irreducible `S_5` | `(1,0,-1/2)` | `S^5-5S^3+4S+1` |
| irreducible `A_5` | `(1,-4/3,6)` | `S^5-5S^3+(8/3)S^2+4S-12` |
| `K_2 x K_3` | `(1,-3/2,-9/2)` | `(S^2+S+1)(S^3-S^2-5S+9)` |
| Hasse failure | `(4,-335/27,4807/20736)` | `(192S^2-72S+19)(55296S^3+20736S^2+1224S-253)/10368` |

Every first coordinate is nonzero and every displayed polynomial is
squarefree of degree five.  The localized reconstruction theorem therefore
identifies each row with the full scheme fiber, not merely a subset of its
affine points.

For the `S_5` row,

\[
 \operatorname{Disc}(S^5-5S^3+4S+1)=38569
\]

is nonsquare.  The polynomial is irreducible modulo `2`, while modulo `79`
its factor degrees are `(2,1,1,1)`.  Thus its transitive Galois group
contains a 5-cycle and a transposition and is `S_5`.

For the alternating row, the monic polynomial has discriminant `984^2`.
After clearing the denominator, it is irreducible modulo `5` and has factor
degrees `(3,1,1)` modulo `7`.  Its transitive group lies in `A_5` and
contains a 3-cycle, forcing the full group `A_5`.

For the product row, the quadratic discriminant is `-3`, the cubic is
irreducible modulo `5`, and the factor resultant is `181`.  Hence this full
fiber is exactly a product of irreducible fields `K_2 x K_3`.

### 6.2 All five transitive quintic groups

The same map realizes every transitive subgroup of `S_5`.  It is convenient
to display the centered monic polynomial

\[
 \widetilde E_{\Pi,B,C}(T)
 =\Pi^5E_{\Pi,B,C}(\Pi^{-2}T)
 =T^5-5T^3-2\Pi BT^2+4\Pi^3T-2\Pi^5C.
\]

The `A_5` and `S_5` rows are those already listed above.  The three solvable
rows are:

| group | target `(\Pi,B,C)` | `\widetilde E_{\Pi,B,C}(T)` |
|---|---|---|
| `C_5` | `(1,-15/11,331/242)` | `T^5-5T^3+(30/11)T^2+4T-331/121` |
| `D_5` | `(5/2,-27/8,-738/3125)` | `T^5-5T^3+(135/8)T^2+(125/2)T+369/8` |
| `F_{20}` | `(31/5,5229/310,9618099/114516604)` | `T^5-5T^3-(5229/25)T^2+(119164/125)T-9618099/6250` |

For the cyclic row, let `\theta` satisfy

\[
 m(\theta)=\theta^5+\theta^4-4\theta^3-3\theta^2+3\theta+1=0
\]

and put

\[
 \eta=\frac{-3\theta^4+\theta^3-4\theta+15}{11}.
\]

Reduction modulo `2` proves that the displayed target polynomial is
irreducible.  Direct reduction in `Q[\theta]/(m)` gives

\[
 m(\theta^2-2)=0,\qquad
 (\theta\mapsto\theta^2-2)^5=\mathrm{id},
\]

with no smaller positive power equal to the identity.  The characteristic
polynomial of `\eta` is the displayed `C_5` polynomial.  Thus its degree-five
field has five automorphisms and is cyclic.  As useful regressions,

\[
 \operatorname{Disc}(\widetilde E)
 =\left(\frac{109\cdot2663}{11^4}\right)^2
\]

and its pair-sum resolvent factors as

\[
\frac1{11^4}
\left(121X^5-1210X^3+715X^2+1364X+23\right)
\left(121X^5-605X^3-385X^2+209X+43\right).
\]

For the dihedral row, irreducibility modulo `11` gives transitivity, while
modulo `7` the factor degrees are `(2,2,1)`.  Its discriminant is

\[
 \left(\frac{2048625}{256}\right)^2,
\]

and its pair-sum resolvent is

\[
\frac1{512}
\left(16X^5+60X^3+500X^2-585X+2196\right)
\left(32X^5-600X^3-460X^2-180X-5553\right),
\]

with both quintic factors irreducible over `Q`.  The two pair orbits put the
transitive group in `D_5` or `C_5`, and the `(2,2,1)` Frobenius element
excludes `C_5`.

For the Frobenius row, start from

\[
 H(U)=U^5-10U^3+20U+20.
\]

It is irreducible modulo `29` and has discriminant
`231200000=2^8\cdot5^5\cdot17^2`.  We use the explicit sextic criterion
from [Dummit, *Solving solvable quintics*](https://site.uvm.edu/ddummit/files/2021/04/Solving_Solvable_Quintics__Math_Comp_57_no195_1991__pp_387_401.pdf):
an irreducible rational quintic has group contained in `F_{20}` exactly when
this resolvent has a rational root.  Here the resolvent is

\[
\begin{aligned}
 \mathcal C_H(X)
 &=(X^3-700X^2+110000X-15880000)^2
   -2^{10}\operatorname{Disc}(H)X\\
 &=(X-500)
 \left(
 X^5-900X^4+260000X^3-55760000X^2
 +6452000000X-504348800000
 \right).
\end{aligned}
\]

The rational Cayley root puts the transitive group inside `F_{20}`; the
nonsquare discriminant excludes `C_5` and `D_5`, so the group is exactly
`F_{20}`.  In `Q[U]/(H)`, the element

\[
 \eta=\frac{22}{5}+\frac95U-\frac{11}{10}U^2-\frac35U^3
\]

has traces

\[
 \operatorname{Tr}(\eta)=0,\qquad
 \operatorname{Tr}(\eta^2)=10,\qquad
 \operatorname{Tr}(\eta^4)=50-16(31/5)^3,
\]

and has the displayed target polynomial as its characteristic polynomial.
This supplies an exact non-affine transport into the fixed inverse pencil.
The bounded trace search that found this element is reproduced by

```bash
.venv/bin/python scripts/search_fixed_quintic_trace_points.py \
  --u -10 --v 20 --bound 18
```

The search is discovery evidence only; the element and all group
certificates above are checked exactly without trusting the search.

Consequently one fixed Keller map contains complete fibers with all five
transitive quintic Galois groups

\[
 \boxed{C_5,\ D_5,\ F_{20},\ A_5,\ S_5.}
\]

### 6.3 The Hasse row

Put

\[
\begin{aligned}
 Q(S)&=192S^2-72S+19,\\
 R(S)&=55296S^3+20736S^2+1224S-253.
\end{aligned}
\]

Then

\[
 E_{4,-335/27,4807/20736}(S)=\frac{Q(S)R(S)}{10368}.
\]

The quadratic is irreducible, the cubic is irreducible modulo `7`, and

\[
\begin{aligned}
 \operatorname{Disc}(Q)&=-3\cdot56^2,\\
 \operatorname{Disc}(R)&=-3\cdot28366848^2,\\
 \operatorname{Res}(Q,R)&=93138374098944\ne0.
\end{aligned}
\]

Thus `R` has Galois group `S_3`, and the quadratic subfield of its splitting
field is exactly the splitting field `Q(sqrt(-3))` of `Q`.  Neither factor
has a rational root.

The displayed integral models can be bad only at

\[
 \{2,3,7,19\}.
\]

Outside those primes, the decomposition group in the common `S_3` splitting
field is cyclic.  An element of order two fixes a cubic root, while an
element of order three acts trivially on the two quadratic roots.  Hence the
product has a root over every unramified completion.

The remaining completions have the following exact witnesses.

- At `2`, put `U=8S`.  The cubic becomes
  \[
   108U^3+324U^2+153U-253,
  \]
  which has the simple root `U=1` modulo `2`.
- At `3`, put `U=3S`.  For
  \[
   f(U)=2048U^3+2304U^2+408U-253
  \]
  one has
  \[
   v_3(f(-7))=7>4=2v_3(f'(-7)).
  \]
- At `7`,
  \[
   v_7(Q(5))=3>2=2v_7(Q'(5)).
  \]
- At `19`, `Q(0)=19` and `Q'(0)=-72` is a unit.

Ordinary or strong Hensel therefore supplies a local root in every
exceptional case.  The cubic supplies a real root.  Consequently this
complete fiber has a point over every completion of `Q` and no rational
point.

For the standard projective height on the target, its primitive coordinates
are

\[
 [20736:82944:-257280:4807],
\]

so the new Hasse target has height `257280`.  This improves the first
certificate `(2,-6741,458080)`, whose height was `458080`; no global
height-minimality is claimed.

The completed bounded search that found the improved row is

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_targets.py
```

It enumerates `|num(Pi)|<=16`, `den(Pi)<=8`, `|Pi|<=4`,
`|num(a)|<=120`, `den(a)<=8`, and quadratic-factor parameter height
`H(b)<=10000`.  It reports only irreducible quadratic-cubic pairs having a
local root at every possibly ramified prime and target height below `458080`.
The two sign-related targets of height `257280` are the only reported rows in
that box.  This is reproducible search evidence, not a global minimality
theorem.

Varying the shared quadratic field produces a smaller second certificate.
For

\[
 (\Pi,B,C)=\left(-7,\frac{387}{14},\frac{400}{2401}\right)
\]

the normalized inverse polynomial is

\[
\widetilde E(T)
=(T^2-4T+32)(T^3+4T^2-21T+175).
\]

Equivalently,

\[
E_{\Pi,B,C}(S)=
-\frac{
(2401S^2-196S+32)
(117649S^3+9604S^2-1029S+175)}
{16807}.
\]

The quadratic discriminant and cubic discriminant are

\[
 -7\cdot4^2,\qquad -7\cdot(5\cdot79)^2,
\]

respectively.  The quadratic is irreducible, while the cubic is irreducible
modulo \(2\), since \(T^3+T+1\) has no root in \(\mathbb F_2\).  Thus the
cubic has Galois group \(S_3\), and its quadratic resolvent is the splitting
field \(\mathbb Q(\sqrt{-7})\) of the quadratic factor.

Outside \(2,5,7,79\), the same unramified common-resolvent argument applies.
At the exceptional primes:

- at \(2\), \(-7\equiv1\pmod8\), so the quadratic splits over
  \(\mathbb Q_2\);
- at \(5\), the cubic has the simple root \(0\) modulo \(5\);
- at \(7\), the cubic has the simple root \(3\) modulo \(7\); and
- at \(79\), the quadratic has the simple root \(31\) modulo \(79\).

The cubic supplies a real root.  Consequently this is another complete
everywhere locally soluble fiber with no rational point.  Its primitive
projective target is

\[
 [4802:-33614:132741:800],
\]

of height \(132741\), improving the preceding \(257280\) certificate.  No
global height-minimality is claimed.  Its independent exact audit is

```bash
.venv/bin/python scripts/verify_fixed_quintic_hasse_minus_seven.py
```

and the bounded varying-discriminant search that found it is

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_discriminants.py
```

A wider run of the same search found a still cleaner row:

\[
 (\Pi,B,C)=\left(5,-\frac{144}{5},-\frac{188}{3125}\right),
\qquad
\widetilde E(T)
=(T^2-8T+47)(T^3+8T^2+12T+8).
\]

Here the factor discriminants are

\[
 -31\cdot2^2,\qquad -31\cdot8^2.
\]

The quadratic is irreducible, and the cubic is irreducible modulo \(5\).
Their common quadratic resolvent is therefore
\(\mathbb Q(\sqrt{-31})\).  Only \(2\) and \(31\) divide the displayed
factor discriminants.  At \(2\), \(-31\equiv1\pmod8\), so the quadratic
splits over \(\mathbb Q_2\).  At \(31\), the cubic has the simple root
\(15\) modulo \(31\).  The unramified common-resolvent argument covers every
other finite prime, and the cubic covers the real place.

Thus this complete fiber is everywhere locally soluble with no rational
point.  Its primitive projective target is

\[
 [3125:15625:-90000:-188],
\]

of height \(90000\).  This improves both preceding certificates, but no
global height-minimality is claimed.  Its independent exact audit is

```bash
.venv/bin/python scripts/verify_fixed_quintic_hasse_minus_thirty_one.py
```

### 6.4 Infinitely many Hasse failures: present frontier

Section 9.2 reduces the common-resolvent condition for a quadratic-cubic
factorization to the explicit square equation

\[
 y^2=-N(\Pi,a,b)(a^2-4b),
\]

where `N` is the displayed cubic-discriminant polynomial used by
`scripts/search_fixed_quintic_hasse_targets.py`.  With `y` included this is
a three-dimensional arithmetic incidence variety.  The Hasse row above is
one rational point on it, and the required local-root conditions are open
at each fixed completion.  These observations make infinitude plausible,
but they do not prove weak approximation or a positive-rank rational curve
through a point satisfying all the required local conditions.

There is nevertheless an exact rational parametrization of the
common-resolvent condition.  Put

\[
 A=\Pi^2a,\qquad V=\Pi^4b,\qquad W=\Pi^3.
\]

On the quadratic-cubic factorization locus the normalized inverse polynomial
has factors

\[
\begin{aligned}
q(T)&=T^2+AT+V,\\
h(T)&=T^3-AT^2+(A^2-V-5)T
       +\frac{4W-V(A^2-V-5)}{A}.
\end{aligned}
\]

Write \(M=-A^2\operatorname{Disc}(h)\), viewed as a quadratic polynomial in
\(W\).  Direct calculation gives

\[
 \operatorname{Disc}_W(M)
 =-256A^2(2A^2-3V-15)^3.
\]

Impose

\[
 V=\frac{A^2+3R^2}{4},\qquad
 H=\frac{5A^2-9R^2-60}{4},
\]

so that \(\operatorname{Disc}(q)=-3R^2\), and put

\[
 L=-\frac{
 5A^4+270A^2R^2+180A^2-243R^4-1620R^2}{2}.
\]

For every nonzero parameter \(\kappa\), set

\[
\boxed{
 W=\frac{8AH(H/\kappa-\kappa)-L}{864},
 \qquad
 Y=\frac{RAH(\kappa+H/\kappa)}{3}.
}
\]

Then the exact identity

\[
 Y^2=-(A^2-4V)M(A,V,W)
\]

holds.  Thus the discriminant-square part of the Hasse construction is
rational; the remaining hard conditions are that \(W\) be a rational cube,
that both factors be irreducible, and that the cubic have a root at every
completion where the quadratic does not.

The proportional slices \(\kappa=cA\) give two especially simple cube
equations.  The coefficient of \(A^4\) in the numerator of \(W\) is
\(-5(c+1)(4c-5)\), and hence it vanishes exactly at
\(c=-1,5/4\).  At these values,

\[
\begin{aligned}
c=-1:\quad&
 W=\frac{(3R^2+5)(3A^2-3R^2-20)}{48},\\
c=\frac54:\quad&
 W=\frac{135A^2R^2-99R^4-420R^2+1600}{960}.
\end{aligned}
\]

For \(c=-1,R=1\), the equation \(W=\Pi^3\) becomes

\[
 (54A)^2=(18\Pi)^3+22356,
\]

and PARI/GP returns the exact rank interval \([2,2]\), with generators
\((-11,145)\) and \((73/4,1349/8)\) on this displayed elliptic curve.  For
\(c=5/4,R=4\), the model

\[
 (12A)^2=(4\Pi)^3+\frac{30464}{15}
\]

has PARI rank interval \([1,1]\).  These positive-rank slices produce
infinitely many points on the common-resolvent incidence, but do not by
themselves establish the required simultaneous local-root conditions.

The reproducible bounded experiment is

```bash
.venv/bin/python scripts/search_fixed_quintic_hasse_curves.py
```

It verifies all identities above, checks the two rank intervals, and searches
\(\operatorname{den}(R),\operatorname{den}(A)\le4\),
\(1\le\operatorname{num}(R)\le30\),
\(|\operatorname{num}(A)|\le30\), and
\(0<|\Pi|\le30\) on all proportional slices.  Of 32 irreducible
common-resolvent presentations, only four have a cubic root over
\(\mathbb Q_2,\mathbb Q_3,\mathbb Q_5\); they are the sign and conic-parameter
presentations of the already known Hasse target.  Since the search tests only
three completions and a bounded box, this is evidence about where the next
search must go, not an exhaustiveness or infinitude result.  No infinitude
theorem is claimed here.

For the particularly clean \(\mathbb Q(\sqrt{-31})\) row, two natural
infinitude attempts admit exact reductions.  First fix its normalized
parameters \(A=-8,R=2\).  Then

\[
 M=16(27W^2-8254W+617811),
\]

and the known point is \((W,Y)=(125,1984)\) on \(Y^2=31M\).  Lines through
that point give the rational conic parametrization

\[
 W=\frac{125\kappa^2-3968\kappa-2419984}
          {\kappa^2-13392}.
\]

The requirement \(W=\Pi^3\) is equivalent to the genus-two curve

\[
 Z^2=1984(27\Pi^6-8254\Pi^3+617811),
\]

which contains the certified point \((\Pi,Z)=(5,3968)\).  Thus this fixed
vertical slice is not an immediate positive-rank elliptic route.  An exact
enumeration of rational \(\Pi\) of height at most \(600\) finds only
\(\Pi=5\) (with the two signs of \(Z\)); this is bounded evidence, not a
determination of all rational points.

There is a complementary fixed-algebra calculation.  Let \(\theta\) satisfy

\[
 \theta^3+8\theta^2+12\theta+8=0
\]

and vary the two generators affinely as

\[
 \eta_2=u\sqrt{-31}+v,\qquad \eta_3=w\theta+s.
\]

Trace zero forces \(s=(8w-2v)/3\), and the second-moment normalization is
the rational quadric

\[
 \boxed{5v^2+28w^2-93u^2=15.}
\]

On this quadric the fourth-moment condition becomes

\[
\Pi^3=-\frac1{216}\left(
25947u^4-5022u^2v^2+35v^4+672v^2w^2
+1504vw^3+2352w^4-675
\right).
\]

The known row is \((u,v,w,\Pi)=(1,4,1,5)\).  These calculations do not prove
that the full fixed-algebra threefold lacks rational curves; they isolate why
the most direct affine generator variation does not immediately prove
infinitude.  The identities and the bounded genus-two point search are
checked by

```bash
.venv/bin/python scripts/analyze_fixed_quintic_hasse_minus_thirty_one.py
```

A fixed-discriminant integral search through
\(R\le100,\ |A|,|\Pi|\le200\) finds only the certified point and its sign
mate.  A rational grid with denominator at most \(4\) and numerator bounds
\((20,30,30)\) for \((R,A,\Pi)\) finds two additional common-resolvent
presentations, but both fail local solubility at \(17\).  These are bounded
experiments, not finiteness results.

There is also a useful exact obstruction to the most obvious attempted
proof.  Begin with the classical intersective algebra

\[
 A_m=\mathbb Q(\sqrt{-3})\times\mathbb Q(\sqrt[3]{m})
\]

and allow independent affine generators

\[
 \eta_2=u\sqrt{-3}+v,\qquad
 \eta_3=w\sqrt[3]{m}+s.
\]

Trace zero forces `2v+3s=0`.  Since a pure cubic generator has first and
second trace moments zero, the normalization
`\operatorname{Tr}(\eta^2)=10` required by the fixed pencil becomes

\[
 5v^2-9u^2=15.
\]

This conic has no `Q_5`-point.  Indeed, for a primitive projective solution

\[
 5V^2-9U^2=15W^2,
\]

reduction modulo `5` gives `5\mid U`; after division by `5`, reduction
modulo `5` gives `V^2=3W^2`, forcing `5\mid V,W` because `3` is a
nonsquare modulo `5`, a contradiction.  Thus the standard infinite
pure-cubic Hasse family cannot enter this fixed normalized trace chart,
even after independent affine changes on its two factors.

The remaining infinitude problem is therefore genuinely about finding
rational curves or a weak-approximation mechanism on the general
common-resolvent threefold, not about importing the classical pure-cubic
family by an affine change of generators.

## 7. Verification

Run

```bash
.venv/bin/python scripts/verify_fixed_quintic_moduli_dominance.py
.venv/bin/python scripts/verify_fixed_quintic_arithmetic_zoo.py
```

The first checker verifies:

1. the seed factorization and discriminant;
2. the full three-variable Jacobian determinant and coordinate degrees;
3. the inverse-polynomial normalization identity;
4. both coefficient-Jacobian certificates;
5. squarefreeness at the rational dominance point; and
6. the three explicit signature and `S_5` specialization certificates.

The second checker verifies the five arithmetic rows, their Galois and
factorization certificates, the common quadratic resolvent, all four local
Hensel witnesses, and the improved target height.

## 8. Fixed-map Hilbert and local-engineering theorem

Let

\[
 U=\{(\pi,b,c)\in\mathbb A^3_{\mathbb Q}:
       \pi\ne0,\ \operatorname{Disc}_S(E_{\pi,b,c})\ne0\}.
\]

The generic splitting field of `E_{\Pi,B,C}` over
`\mathbb Q(\Pi,B,C)` is regular with Galois group `S_5`, by the general
symmetric-monodromy theorem for quadratic-gauge seeds.

> **Theorem (fixed-map arithmetic specialization).**
> Let `\Omega_\infty\subset U(\mathbb R)` be a nonempty Euclidean open set,
> let `S` be a finite set of rational primes, and for every `p\in S` let
> `\Omega_p\subset U(\mathbb Q_p)` be a nonempty `p`-adic open set.  Then
> there are infinitely many
> \[
>  y\in U(\mathbb Q)\cap\Omega_\infty
>       \cap\bigcap_{p\in S}\Omega_p
> \]
> for which `E_y` is irreducible and has splitting-field group `S_5`.
> Every such polynomial represents a connected full fiber of the one fixed
> Keller map `F`.

This is the weak-approximation form of Hilbert irreducibility applied to the
regular `S_5` splitting cover of `U`.  The complement of the corresponding
Hilbert subset is thin, while rational affine three-space satisfies weak
approximation.  Intersecting the Hilbert subset with the specified finite
product of local opens gives infinitely many rational targets.  The full
scheme-fiber assertion then follows from the quadratic-gauge reconstruction
because `y\in U`.

Real-root count is locally constant on `U(\mathbb R)`, so the three explicit
signature points in Section 6 supply admissible choices of `\Omega_\infty`
for all quintic signatures.

If a target has integral coordinates at `p`, `\pi` is a unit, and its
reduction is squarefree with factor-degree partition `\lambda`, then the
whole residue box is a valid `\Omega_p`.  Every fiber selected from that box
has unramified local degree partition `\lambda`.

For this fixed map, all seven partitions of five already occur modulo `7`:

| partition | `(\bar\Pi,\bar B,\bar C)` in `F_7^3` |
|---|---|
| `(5)` | `(1,0,1)` |
| `(4,1)` | `(1,0,3)` |
| `(3,2)` | `(1,0,2)` |
| `(3,1,1)` | `(1,1,3)` |
| `(2,2,1)` | `(3,2,0)` |
| `(2,1,1,1)` | `(1,2,6)` |
| `(1,1,1,1,1)` | `(1,0,0)` |

The checker verifies squarefreeness and the exact factor-degree pattern in
every row.  More generally, finite-field Chebotarev for the regular `S_5`
cover shows that every conjugacy class occurs at all sufficiently large good
primes.  Therefore arbitrary finite collections of unramified quintic
splitting types at sufficiently large primes can be imposed simultaneously,
together with any quintic signature, on infinitely many connected full
fibers of this one fixed map.

## 9. Arithmetic loci in the fixed target

### 9.1 The discriminant and its normalization

The discriminant factors as

\[
 \operatorname{Disc}_S(E_{\Pi,B,C})=16\Pi^8\Delta(\Pi,B,C),
\]

where

\[
\begin{aligned}
\Delta={}&432B^5C\Pi^2-432B^4\Pi^2
+12600B^3C\Pi^3-2000B^3C\\
&+9000B^2C^2\Pi^7+20625B^2C^2\Pi^4
-11520B^2\Pi^3+2000B^2\\
&+18750BC^3\Pi^8-25600BC\Pi^7
+56000BC\Pi^4-45000BC\Pi\\
&+3125C^4\Pi^{12}-40000C^2\Pi^8
+112500C^2\Pi^5-84375C^2\Pi^2\\
&+16384\Pi^7-51200\Pi^4+40000\Pi.
\end{aligned}
\]

On `\Pi\ne0`, the repeated-root incidence gives a rational normalization
chart for this hypersurface.  If `r\ne0` is the repeated root, then

\[
 \boxed{
 B=\frac{5\Pi^5r^4-15\Pi r^2+4}{4r},
 \qquad
 C=-\frac r4(3\Pi^5r^4-5\Pi r^2-4).
 }
\]

The root `r=0` never lies on the discriminant because
`E'_{\Pi,B,C}(0)=4`.

Since `16\Pi^8` is a square, an irreducible specialization has Galois group
contained in `A_5` exactly when `\Delta(\Pi,B,C)` is a rational square.
Thus the oriented-quintic pullback is the explicit double cover

\[
 W^2=\Delta(\Pi,B,C).
\]

For example,

\[
 (\Pi,B,C)=\left(1,-\frac43,6\right)
\]

gives

\[
 S^5-5S^3+\frac83S^2+4S-12,
\qquad
 \operatorname{Disc}=968256=984^2.
\]

It is irreducible modulo `5` and has factorization type `(3,1,1)` modulo
`7`.  Hence its transitive Galois group lies in `A_5`, contains a 5-cycle
and a 3-cycle, and is exactly `A_5`.

### 9.2 Rational-root and `2+3` incidences

A rational root `r` is equivalent to

\[
 \boxed{
 C=\frac{\Pi^5r^5-5\Pi r^3-2Br^2+4r}{2}.
 }
\]

This is a rational threefold incidence over the target.  Its rational image
is the locus of fibers with a rational point; it is an arithmetic thin image,
not a geometric divisor.

For the `2+3` incidence, prescribe a monic quadratic factor

\[
 q(S)=S^2+aS+b.
\]

On the chart `a\ne0`, division by `q` gives the unique target

\[
 \boxed{
\begin{aligned}
B={}&-\frac{\Pi^5a^4-3\Pi^5a^2b+\Pi^5b^2
                  -5\Pi a^2+5\Pi b+4}{2a},\\
C={}&\frac{b(\Pi^5a^2b-\Pi^5b^2-5\Pi b-4)}{2a}.
\end{aligned}}
\]

The cubic quotient is

\[
\begin{aligned}
\Pi^5S^3-\Pi^5aS^2
(\Pi^5(a^2-b)-5\Pi)S
\frac{-\Pi^5a^2b+\Pi^5b^2+5\Pi b+4}{a}.
\end{aligned}
\]

On the omitted chart `a=0`, divisibility requires

\[
 \Pi^5b^2+5\Pi b+4=0,\qquad C=Bb.
\]

These formulas turn the search for intersective `2+3` fibers into explicit
conditions on `(\Pi,a,b)`: irreducibility of both factors, equality of the
quadratic field with the cubic discriminant field, and the ramified-prime
local-root conditions.

### 9.3 The two-subset resolvent

For a normalized quintic

\[
 f(T)=T^5-5T^3+aT^2+bT+c
\]

define the monic pair-sum resolvent

\[
 R_{10}(X)=\prod_{i<j}(X-r_i-r_j).
\]

An exact resultant calculation gives

\[
\begin{aligned}
R_{10}(X)={}&X^{10}-15X^8+aX^7+(75-3b)X^6
 +(-10a-11c)X^5\\
&+(-a^2+10b-125)X^4
 +(-4ab+25a+20c)X^3\\
&+(5a^2+7ac-4b^2+25b)X^2\\
&+(-a^3+4bc-25c)X-a^2b-5ac-c^2.
\end{aligned}
\]

The checker proves the identity

\[
 \operatorname{Res}_Y(f(Y),f(X-Y))
 =2^5f(X/2)R_{10}(X)^2.
\]

Substituting

\[
 (a,b,c)=(-2\Pi B,4\Pi^3,-2\Pi^5C)
\]

pulls this resolvent back explicitly to the Keller target:

\[
\begin{aligned}
R_{10}^{F}(X)={}&X^{10}-15X^8-2B\Pi X^7
+(75-12\Pi^3)X^6\\
&+(20B\Pi+22C\Pi^5)X^5
+(-125+40\Pi^3-4B^2\Pi^2)X^4\\
&+(32B\Pi^4-50B\Pi-40C\Pi^5)X^3\\
&+(20B^2\Pi^2+28BC\Pi^6-64\Pi^6+100\Pi^3)X^2\\
&+(8B^3\Pi^3-32C\Pi^8+50C\Pi^5)X\\
&-16B^2\Pi^5-20BC\Pi^6-4C^2\Pi^{10}.
\end{aligned}
\]

For an irreducible quintic, the action on the ten unordered root pairs has:

- one orbit for `S_5`, `A_5`, and `F_{20}`; and
- two five-element orbits for `D_5` and `C_5`.

The `F_{20}` assertion is important: its natural degree-five action is
sharply two-transitive.  Consequently `R_{10}^{F}` is irreducible in the
generic `S_5/A_5/F_{20}` cases, whereas a factorization into two irreducible
quintics detects the combined `D_5/C_5` region.  A modular `(2,2,1)` cycle
then distinguishes `D_5` from `C_5`; an explicit order-five automorphism
certifies the cyclic row in Section 6.2.  Distinguishing `F_{20}` from
`S_5` requires a different invariant, supplied there by Cayley's sextic
solvability resolvent.

As with rational factorization of `E`, reducibility of `R_{10}^{F}` over
`\mathbb Q` is an arithmetic thin-image condition, not a Zariski-closed
subvariety of the geometric target.

## 10. Descent while varying the primitive generator

### 10.1 Exact obstruction for one presentation

Let

\[
 f(T)=T^5+aT^3+bT^2+cT+d
\]

be a centered monic squarefree quintic with `ac\ne0`.  It is affinely
equivalent over `\mathbb Q` to an inverse polynomial of the fixed Keller map
if and only if

\[
 \boxed{
 -\frac a5\in\mathbb Q^{\times2},
 \qquad
 \frac{25c}{4a^2}\in\mathbb Q^{\times3}.
 }
\]

Indeed, choose `\lambda,\pi\in\mathbb Q^*` with

\[
 \lambda^2=-\frac a5,
 \qquad
 \pi^3=\frac{25c}{4a^2}.
\]

Then

\[
 \frac{f(\lambda T)}{\lambda^5}
 =
 T^5-5T^3+\frac b{\lambda^3}T^2
             +4\pi^3T+\frac d{\lambda^5},
\]

and the required target is

\[
 \boxed{
 \Pi=\pi,\qquad
 B=-\frac{b}{2\pi\lambda^3},\qquad
 C=-\frac{d}{2\pi^5\lambda^5}.
 }
\]

Conversely, the normalized inverse family has exactly these square and cube
classes, proving necessity.

### 10.2 Intrinsic trace form

Let `A` be a rank-five finite etale `\mathbb Q`-algebra and let
`\eta\in A` have trace zero.  Put

\[
 p_2(\eta)=\operatorname{Tr}(\eta^2),
 \qquad
 p_4(\eta)=\operatorname{Tr}(\eta^4).
\]

The characteristic polynomial of `\eta` has the form

\[
 T^5+aT^3+bT^2+cT+d
\]

with

\[
 \boxed{
 a=-\frac{p_2}{2},
 \qquad
 c=\frac{p_2^2-2p_4}{8}.
 }
\]

Therefore the two presentation obstructions become

\[
 \boxed{
 \frac{p_2}{10}\in\mathbb Q^{\times2},
 \qquad
 \frac{25(p_2^2-2p_4)}{8p_2^2}
       \in\mathbb Q^{\times3}.
 }
\]

Translations and nonzero rational scalings of `\eta` do not change either
class.

### 10.3 The quadratic obstruction always disappears

Let `A_0=\ker(\operatorname{Tr}_{A/\mathbb Q})`, a four-dimensional rational
space.  The trace form

\[
 q(\eta)=\operatorname{Tr}(\eta^2)
\]

is nondegenerate on `A_0`: the full trace pairing is nondegenerate and
`\mathbb Q\cdot1` is its nondegenerate orthogonal complement.

Consider the five-variable quadratic form

\[
 q(\eta)-10s^2.
\]

Over the reals it is indefinite for every possible signature of a rank-five
etale algebra.  It is nondegenerate and has dimension five, so Meyer's
theorem makes it isotropic over `\mathbb Q`.  Its projective quadric is then
rational, and its rational points are Zariski dense.  Consequently one may
avoid `s=0` and the proper closed nonprimitive-element locus.

Thus every rank-five finite etale algebra has a primitive trace-zero element
`\eta` with

\[
 \operatorname{Tr}(\eta^2)=10s^2,\qquad s\ne0.
\]

After replacing `\eta` by `\eta/s`, the quadratic obstruction is normalized
uniformly:

\[
 \boxed{\operatorname{Tr}(\eta^2)=10,\qquad a=-5.}
\]

### 10.4 The single remaining descent variety

After this normalization, the cube condition becomes

\[
 \boxed{
 \operatorname{Tr}(\eta^4)=50-16\Pi^3.
 }
\]

Define

\[
 \mathcal V_A:
 \quad
 \operatorname{Tr}(\eta)=0,\qquad
 \operatorname{Tr}(\eta^2)=10,\qquad
 \operatorname{Tr}(\eta^4)=50-16\Pi^3.
\]

Here `\eta\in A` and `\Pi\in\mathbb A^1`; after eliminating the trace-linear
coordinate this is a threefold.  The fixed map realizes `A` if
`\mathcal V_A(\mathbb Q)` contains a point with `\Pi\ne0` and primitive
`\eta`.  Conversely, every realization on the centered nonzero-cubic chart
produces such a point.  Once a point is found, if

\[
 \chi_\eta(T)=T^5-5T^3+bT^2+4\Pi^3T+d,
\]

the target is simply

\[
 \boxed{
 (\Pi,B,C)=
 \left(\Pi,-\frac{b}{2\Pi},-\frac{d}{2\Pi^5}\right).
 }
\]

This removes the quadratic descent obstruction for every algebra and
isolates the rational-point problem on `\mathcal V_A` as the only remaining
obstruction in this fixed-seed strategy.  Universal rational solubility of
these cubic covers is open here.

### 10.5 Kummer-cover and projective-complete-intersection form

There is a useful geometric simplification of the remaining obstruction.
Put

\[
 Q_A=\left\{\eta\in A_0:
        \operatorname{Tr}(\eta^2)=10\right\},
 \qquad
 D_A=\left\{\eta\in Q_A:
        \operatorname{Tr}(\eta^4)=50\right\}.
\]

The nondegeneracy of the trace form makes `Q_A` a smooth affine quadric
threefold.  Section 10.3 gives it a rational point, so `Q_A` is rational
over `\mathbb Q`.  The descent variety is the cyclic cubic cover

\[
 \boxed{
 \mathcal V_A\longrightarrow Q_A,\qquad
 \Pi^3=\frac{50-\operatorname{Tr}(\eta^4)}{16}.
 }
\]

Consequently, over `Q_A\setminus D_A` its `\Pi\ne0` part is a
`\boldsymbol\mu_3`-torsor.  In particular, the only remaining condition is
not an additional quadratic-form problem: it is a Kummer class on a
rational threefold.

The same observation gives a compact geometric model.  Introduce homogeneous
coordinates `[\eta:s:W]` after eliminating `\operatorname{Tr}(\eta)=0`, and
let

\[
 \overline{\mathcal V}_A\subset\mathbb P^5
\]

be cut out by

\[
 \boxed{
 \operatorname{Tr}(\eta^2)=10s^2,\qquad
 \operatorname{Tr}(\eta^4)=50s^4-16sW^3.
 }
\]

This is a `(2,4)` complete intersection of dimension three.  Its dualizing
sheaf is

\[
 \omega_{\overline{\mathcal V}_A}
 \simeq\mathcal O_{\overline{\mathcal V}_A}(2+4-6)
 \simeq\mathcal O_{\overline{\mathcal V}_A}.
\]

It is not smooth as written: the rational boundary point
`[\eta:s:W]=[0:0:1]` is singular.  On the actual realization locus
`sW\ne0`, however, smoothness is immediate.  The Jacobian minor in the
`(s,W)` columns is

\[
 \det
 \begin{pmatrix}
 -20s&0\\
 -200s^3+16W^3&48sW^2
 \end{pmatrix}
 =-960s^2W^2\ne0.
\]

Over a separable closure, writing
`A\simeq\overline{\mathbb Q}^{\,5}` identifies this model with the symmetric
power-sum complete intersection

\[
 \sum_i\eta_i=0,\qquad
 \sum_i\eta_i^2=10s^2,\qquad
 \sum_i\eta_i^4=50s^4-16sW^3.
\]

Thus `\overline{\mathcal V}_A` is the `S_5`-twist of one fixed split
complete intersection by the torsor of ordered geometric factors of `A`.
This packages all algebra-dependence into the twist and all remaining
arithmetic into rational points on its smooth Kummer open.

### 10.6 A non-affine generator succeeds for `T^5-T-1`

Let

\[
 A=\mathbb Q[\theta]/(\theta^5-\theta-1).
\]

The original polynomial has vanishing cubic coefficient and lies outside
the normalization chart.  Nevertheless, take

\[
 \boxed{\eta=\theta^3-\theta^2+2\theta.}
\]

Exact multiplication-matrix calculation gives

\[
 \operatorname{Tr}(\eta)=0,\qquad
 \operatorname{Tr}(\eta^2)=10,\qquad
 \operatorname{Tr}(\eta^4)=-78=50-16\cdot2^3,
\]

and

\[
 \chi_\eta(T)=T^5-5T^3-13T^2+32T-23.
\]

Thus `\Pi=2` and the same fixed Keller map realizes this field at

\[
 \boxed{
 (\Pi,B,C)=
 \left(2,\frac{13}{4},\frac{23}{64}\right).
 }
\]

This example shows concretely why failure of the square/cube test for one
chosen defining polynomial is not an obstruction for the abstract algebra:
a genuinely non-affine change of primitive generator can remove it.

### 10.7 Exact audit on three nonisomorphic `S_5` fields

The preceding small generator is not isolated.  The following three
power-basis fields have squarefree, pairwise distinct discriminants, hence
are pairwise nonisomorphic.  In every row the displayed element has trace
zero and second trace moment `10`; its fourth moment is
`50-16\Pi^3`, and the listed target reconstructs its characteristic
polynomial.

| defining polynomial for `\theta` | field discriminant | `\eta` | `(\Pi,B,C)` |
|---|---:|---|---|
| `T^5-T-1` | `2869` | `\theta^3-\theta^2+2\theta` | `(2,13/4,23/64)` |
| `T^5-4T^3+T+1` | `-55563` | `(-4+16\theta+13\theta^2-5\theta^3-3\theta^4)/7` | `(8/7,209/784,2273/16384)` |
| `T^5-4T^3-T-1` | `-179467` | `\theta-\theta^3/2` | `(-3/4,5/4,-784/243)` |

Each Galois group is certified as `S_5` without a black-box group
calculation.  Reduction modulo the primes shown below gives one irreducible
factorization `(5)` and one factorization `(3,2)`:

\[
\begin{array}{c|cc}
 T^5-T-1 & p=3:(5) & p=2:(3,2)\\
 T^5-4T^3+T+1 & p=13:(5) & p=2:(3,2)\\
 T^5-4T^3-T-1 & p=11:(5) & p=2:(3,2).
\end{array}
\]

The first pattern proves transitivity and supplies a 5-cycle; the second
supplies an odd element of order six.  Among the transitive subgroups of
`S_5`, this forces `S_5`.

For the last two Hermite presentations, the direct affine criterion already
fails because `-u/5=4/5` is not a rational square.  Their successful
generators are therefore genuinely non-affine evidence for rational points
on the twisted Kummer model, rather than points visible in the original
coefficient chart.

The bounded search that found and independently rechecks such points is

```text
.venv/bin/python scripts/search_fixed_quintic_trace_points.py \
  --u -4 --v 1 --bound 18
.venv/bin/python scripts/search_fixed_quintic_trace_points.py \
  --u -4 --v -1 --bound 18
```

This computation is evidence only; the displayed rows themselves are exact
certificates verified by
`scripts/verify_fixed_quintic_moduli_dominance.py`.

## 11. Comparison with the two-parameter generic `S_5` polynomial

The arithmetic issue is clearer when compared with the classical generic
family

\[
 H_{u,v}(T)=T^5+uT^3+vT+v.
\]

This polynomial is generic for `S_5` over `\mathbb Q`: every `S_5`-extension
of `\mathbb Q` is the splitting field of a rational specialization.  A
convenient published statement is
[Ghitza--Yamauchi, Section 4.1](https://www.numdam.org/item/10.5802/jtnb.1291.pdf);
the more general Hermite reduction of a separable quintic to
`T^5+bT^3+cT+d` is explained by
[Kraft](https://arxiv.org/abs/math/0403323).  This two-parameter phenomenon
is compatible with
`\operatorname{ed}_{\mathbb Q}(S_5)=2`
([Buhler--Reichstein, Theorem 6.5](https://www.maths.ed.ac.uk/cheltsov/ed/pdf/BUHLER97.pdf)).
It concerns arbitrary primitive-generator changes, not merely affine
changes of a fixed generator, so it does not contradict the
three-dimensional affine-configuration moduli calculation above.

There is an exact pullback relation with the fixed Keller family.  Put

\[
 \boxed{
 u=-5\lambda^2,\qquad
 v=4\lambda^4\Pi^3,\qquad
 B=0,\qquad
 C=-\frac{2}{\Pi^2\lambda}.
 }
\]

Then

\[
 \boxed{
 \Pi^5E_{\Pi,0,-2/(\Pi^2\lambda)}(\Pi^{-2}T)
 =\frac{H_{-5\lambda^2,\,4\lambda^4\Pi^3}(\lambda T)}
        {\lambda^5}.
 }
\]

Equivalently, a direct affine use of a Hermite generator reaches this
fixed map precisely after solving

\[
 \lambda^2=-\frac u5,\qquad
 \Pi^3=\frac{25v}{4u^2}.
\]

The parameter map

\[
 (\lambda,\Pi)\longmapsto
 (u,v)=(-5\lambda^2,4\lambda^4\Pi^3)
\]

has Jacobian `-120\lambda^5\Pi^2`, and is generically a degree-six cover.
Thus the standard two-parameter generic polynomial does **not** by itself
prove rational universality of the fixed Keller map: its direct comparison
reproduces exactly the quadratic and cubic descent classes already found in
Section 10.1.

This also sharpens the remaining question.  For `S_5` quintic fields it is
enough to study non-affine primitive generators of specializations of
`H_{u,v}`; any successful universal argument must split the above degree-six
cover by a further Tschirnhaus transformation, not by affine normalization
alone.  The Kummer threefold of Section 10.5 is an intrinsic formulation of
that required extra transformation.

## 12. A rational curve of genuinely non-affine `S_5` realizations

The non-affine successes are not isolated.  Let `\tau` be a rational
parameter and put

\[
\begin{aligned}
 u_\tau&=
 \frac{(4-\tau^2)(1+5\tau^2)}{49\tau^2},\\
 v_\tau&=
 \frac{(4-\tau^2)^2(1+5\tau^2)}{343\tau^4},\\
 \beta_\tau&=
 \frac{343\tau^3}
 {(4-\tau^2)^2(1+5\tau^2)},\\
 \Pi_\tau&=\frac7{4-\tau^2}.
\end{aligned}
\]

For `\tau\ne0,\pm2`, let

\[
 A_\tau=
 \mathbb Q[\theta]/
 \left(\theta^5+u_\tau\theta^3+v_\tau\theta+v_\tau\right)
\]

and take the sparse cubic Tschirnhaus generator

\[
 \boxed{\eta_\tau=\beta_\tau\theta^3.}
\]

Exact Newton-sum calculation gives the identities

\[
 \boxed{
 \operatorname{Tr}(\eta_\tau)=0,\qquad
 \operatorname{Tr}(\eta_\tau^2)=10,\qquad
 \operatorname{Tr}(\eta_\tau^4)=50-16\Pi_\tau^3.
 }
\]

The corresponding target of the one fixed Keller map is

\[
\boxed{
\begin{aligned}
 \Pi_\tau&=\frac7{4-\tau^2},\\
 B_\tau&=
 \frac{21\tau(6-5\tau^2)}
 {2(4-\tau^2)(1+5\tau^2)},\\
 C_\tau&=
 -\frac{7\tau^3(4-\tau^2)}
 {2(1+5\tau^2)^2}.
\end{aligned}
}
\]

Indeed, the characteristic polynomial of `\eta_\tau` is

\[
\begin{aligned}
 T^5-5T^3
&+\frac{147\tau(5\tau^2-6)}
 {(4-\tau^2)^2(1+5\tau^2)}T^2\\
&+\frac{1372}{(4-\tau^2)^3}T
+\frac{117649\tau^3}
 {(4-\tau^2)^4(1+5\tau^2)^2},
\end{aligned}
\]

which is exactly

\[
 T^5-5T^3-2\Pi_\tau B_\tau T^2
 +4\Pi_\tau^3T-2\Pi_\tau^5C_\tau.
\]

This curve is genuinely non-affine over the generic point.  The original
Hermite generator would require

\[
 -\frac{u_\tau}{5}
 =
 \frac{(\tau^2-4)(1+5\tau^2)}{245\tau^2}
\]

to be a square in `\mathbb Q(\tau)`.  It is not: the valuations at
`\tau=2` and `\tau=-2` are odd.

Finally, the generic Galois group on this curve is `S_5`.  At `\tau=1`,

\[
 H_1(T)=T^5+\frac{18}{49}T^3+\frac{54}{343}T+\frac{54}{343}.
\]

It is irreducible modulo `47`, while modulo `5` it has factor pattern
`(3,2)`.  Thus this specialization has group `S_5`, forcing the generic
group over `\mathbb Q(\tau)` to be `S_5`.  After removing the finite
discriminant locus, Hilbert irreducibility therefore supplies infinitely
many rational `S_5` specializations on this single explicitly parametrized
non-affine curve of full fibers.

## 13. A dominant rational surface over the Hermite base

The curve above is the slice `\rho=1` of a two-parameter construction.
For independent parameters `\rho,\tau`, abbreviate

\[
 D=4\rho^3-\tau^2,\qquad H=1+20\rho^3
\]

and define

\[
\boxed{
\begin{aligned}
 u&=\frac{9D(1+5\tau^2)}{\tau^2H^2},&
 v&=\frac{27D^2(1+5\tau^2)}{\tau^4H^3},\\
 \beta&=\frac{\tau^3H^3}{27D^2(1+5\tau^2)},&
 \Pi&=\frac{\rho H}{3D}.
\end{aligned}
}
\]

Let `\theta` satisfy the Hermite equation

\[
 \theta^5+u\theta^3+v\theta+v=0
\]

and again take `\eta=\beta\theta^3`.  The general power-trace identities

\[
\begin{aligned}
 \operatorname{Tr}((\beta\theta^3)^2)
 &=-2\beta^2u(u^2-3v),\\
 \operatorname{Tr}((\beta\theta^3)^4)
 &=2\beta^4
 (u^6-6u^4v+9u^2v^2-6uv^2-2v^3)
\end{aligned}
\]

specialize to

\[
 \boxed{
 \operatorname{Tr}(\eta)=0,\qquad
 \operatorname{Tr}(\eta^2)=10,\qquad
 \operatorname{Tr}(\eta^4)=50-16\Pi^3.
 }
\]

Thus this entire rational surface maps into the fixed Keller family.  One
compact form of the target is

\[
\boxed{
 B=-\frac{3\beta^3v(u^2-v)}{2\Pi},
 \qquad
 C=-\frac{\beta^5v^3}{2\Pi^5}.
}
\]

The parameter map to the Hermite plane is dominant.  Its exact Jacobian is

\[
\boxed{
 \det\frac{\partial(u,v)}{\partial(\rho,\tau)}
 =
 \frac{
 29160\rho^2D^3(1+5\tau^2)^2
 }{
 \tau^7H^6
 }.
}
\]

More precisely, away from the displayed exceptional divisors, the inverse
radicals are

\[
\boxed{
 \tau^2=
 -\frac{u(u^2-3v)}{5v^2},
 \qquad
 \rho^3=
 -\frac{(3u+v)(u^2-3v)}{20u^2v}.
}
\]

Hence this is another generically degree-six cover of the two-parameter
generic `S_5` base, but unlike Section 11 it is produced by a genuinely
non-affine cubic Tschirnhaus generator.  Indeed,

\[
 -\frac u5=
 -\frac{9D(1+5\tau^2)}{5\tau^2H^2}
\]

is not a square in `\mathbb Q(\rho,\tau)` because the divisor `D=0`
occurs to odd order.  The `S_5` specialization `(\rho,\tau)=(1,1)` from
Section 12 shows that the generic group of this surface is `S_5`.

This does not yet give arithmetic surjectivity on all rational `(u,v)`:
lifting a rational Hermite point requires the displayed square and cube
classes to vanish.  It does, however, put a dominant rational
two-parameter family of genuinely non-affine `S_5` realizations inside the
one fixed Keller map, and supplies a second explicit descent chart whose
arithmetic image is different from the direct affine chart.
