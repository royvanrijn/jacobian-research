# Virtual norm rows and their affine divisor support

The specialization theorem of
[Gillibert–Levin, Theorems 1.1 and 2.4](sources/gillibert-levin-torsion-line-bundles.pdf)
constructs ideal classes from torsion line bundles. The statements here
go in the other direction for a fixed list of functions in a quadratic
family. Arithmetic virtuality on an unbounded integral set restricts
where their affine divisors may fail to be divisible. The restriction
uses the two valuations at a split prime individually.

## A support bound for simultaneous virtuality

**Theorem VF1.** Let \(F\in\mathbb Q[t]\) be nonconstant and squarefree,
let \(C^\circ:y^2=F(t)\), and fix nonzero rational functions
\(\alpha_1,\ldots,\alpha_h\in\mathbb Q(C^\circ)\). Fix an odd prime \(l\).
For an integer \(n\) with \(F(n)\) a nonzero nonsquare, call a specialized
row virtual if every prime-ideal valuation of \(\alpha_i(n)\) in
\(K_n=\mathbb Q(\sqrt{F(n)})\) is divisible by \(l\).
Exclude the finitely many fibers at zeros or poles of the rows.

Let \(\mathcal H\) be the set of primitive irreducible polynomials
\(g\in\mathbb Z[t]\), up to sign, for which some closed point above \(g\)
has order in the divisor of some \(\alpha_i\) nonzero modulo \(l\).
There is a finite set of rational primes \(S\), depending on these
fixed data, such that the whole list is virtual away from \(S\) at
\(n\) exactly when
\[
           v_p(g(n))\equiv0\bmod l
              \qquad(g\in\mathcal H,\ p\notin S).       \tag{1}
\]
If \(\mathcal H\) contains a polynomial of degree at least two or two
different linear polynomials, there are only finitely many such \(n\),
with an effective bound. If \(\mathcal H\) consists of one linear
polynomial, there are \(O(T^{1/l})\) such \(n\) with \(|n|\le T\).
In particular virtuality on all but finitely many members of one
integer progression of nonsquare specializations forces
\(\mathcal H=\varnothing\).

Multiplication by a rational-function \(l\)th power clears denominators
in \(t\), leaving the virtual condition and the support modulo \(l\)
unchanged. We may take each \(\alpha_i=A_i(t)+B_i(t)y\) regular.
At all but finitely many rational primes, the distinct irreducible
factors of the norms and \(F\) are separately simple modulo that
prime and have nonzero leading coefficients and pairwise nonzero
resultants after reduction. Include denominator primes and the primes
from each of the following local factorizations in \(S\).

Fix a factor \(g\). If the quadratic cover over \(\mathbb Q[t]/g\)
has one unramified point, with row order \(a\), division by \(g^a\)
leaves a function that is regular and nonzero along that fiber.
Its norm has order zero at \(g\), and for good primes dividing \(g(n)\)
the function evaluates to a local unit at every prime over \(p\).
The valuation of a row there is thus \(a\,e(\mathfrak p/p)v_p(g(n))\).
If the cover is ramified at \(g\), divide locally by \(y^a\).
The remaining factor has order zero and again specializes to a unit,
and \(v_{\mathfrak p}(y)=e(\mathfrak p/p)v_p(g(n))/2\) since \(F\)
is squarefree. Because \(l\) is odd, both formulas give exactly
the condition in (1) when \(l\) does not divide \(a\).
The local unit assertions follow by writing the function and its
inverse with denominators in \(\mathbb Q[t]\) coprime to \(g\),
then excluding the fixed denominator and resultant primes.

If \(g\) splits in the quadratic cover, its two \(y\)-coordinates
are \(s(t),-s(t)\) in \(\mathbb Q[t]/g\) with \(s\ne0\).
For good \(p\mid g(n)\) these give two distinct roots modulo \(p\)
and hence two split completions. Let their generic row orders be
\(a,b\). Divide the row by \(g^{\min(a,b)}\). The remaining row is
regular along both branches and a unit at at least one of them.
At a good \(p\) this unit remains a unit, by reduction of its residue
in \(\mathbb Q[t]/g\) and a resultant. Its norm has \(g\)-exponent
\(|a-b|\), with the rest prime to \(g\). The two specialized row
valuations are therefore \(a\,v_p(g(n)),b\,v_p(g(n))\), in the
order specified by the two reductions. This proves (1).
Primes away from \(S\) and all norm factors give unit rows.

For any \(g\in\mathcal H\) we now have \(g(n)=c z^l\), where \(z\)
is an integer and \(c\) is from the fixed finite list of signed
products of primes of \(S\) with exponents \(0,\ldots,l-1\).
For \(\deg g\ge2\), this is the effective superelliptic theorem of
[Bérczes–Evertse–Győry, Theorem 2.1](sources/berczes-evertse-gyory-hyper-superelliptic-excerpt.pdf).
For two different linear polynomials, eliminate \(n\). The result
is \(c_1 z_1^l+c_0=c_2 z_2^l\) with nonzero fixed rational
coefficients and \(c_0\ne0\). Its degree-\(l\) polynomial has
simple roots, so the same theorem applies. For one linear
polynomial just count the \(z\)'s of size \(O(T^{1/l})\).
Alternatively, any nonconstant \(g\) has a simple root modulo
infinitely many good primes, by Chebotarev in its splitting field.
For a prime coprime to the progression modulus, combine the
progression with an exact valuation-one class modulo \(p^2\).
This directly excludes that \(g\) on a full nonsquare progression
as in the theorem.

## The ordinary class rank on a Pell chart

**Corollary VF2.** Suppose in VF1 that \(F=A(t)^2+\sigma\), where
\(A\in\mathbb Z[t]\) is nonconstant and \(\sigma=\pm1\).
Let \(C\) be the smooth projective model, let
\(J=\operatorname{Jac}(C)\), and put \(m=\deg A\).
Assume that the list is virtual at every prime of \(K_n\) for
infinitely many integer parameters \(n\). Write
\(r_\circ=\dim_{\mathbb F_l}\operatorname{Pic}(C^\circ)[l]\).
For all sufficiently large parameters in this common set, the
ordinary class \(l\)-torsion span of the list has dimension at most
\(r_\circ+1\). The additional direction requires the unique
\(g\in\mathcal H\) to be linear and split into two rational points
\(P_+,P_-\). It also requires
\[
                   [P_+]\in l\,\operatorname{Pic}(C^\circ).
                   \tag{2}
\]
If either condition fails, the bound is \(r_\circ\).
On a full common progression the bound is \(r_\circ\), and for
\(l\nmid m\) this equals \(\dim_{\mathbb F_l}J(\mathbb Q)[l]\).

Let \(U=A+y\). Its norm is \(-\sigma\), so its specialization is
an integral unit. The affine unit group is \(\mathbb Q^\times U^{\mathbb Z}\).
Indeed an affine unit \(a+b y\) has constant norm. If \(b\ne0\),
the leading terms cancel in the norm only with
\(\deg a=\deg b+m\). At one of the two rational infinities the
leading terms of the unit then add, giving a pole of order at
least \(m\). The conjugate infinity difference therefore has exact
order \(m\) in \(J(\mathbb Q)\), since \(U\) realizes pole order \(m\).
Taking the quotient of projective divisors by the infinity points
gives
\(\operatorname{Pic}(C^\circ)=J(\mathbb Q)/\langle\infty_+-\infty_-\rangle\).
When \(l\nmid m\), its \(l\)-torsion is the stated Jacobian torsion.

Consider the classes of the rows in
\[
  \frac{\mathbb Q(C)^\times}
       {\mathbb Q(t)^\times U^{\mathbb Z}\,\mathbb Q(C)^{\times l}}.
                   \tag{3}
\]
Every relation among a finite list in this quotient holds on
specialization outside finitely many fibers. Since the rows are
virtual and \(U(n)\) is an integral unit, the rational factor in
the relation is virtual. A rational number virtual in a quadratic
field is a rational \(l\)th power, as its valuations are multiplied
by indices one or two. Hence the class-group span of the list is
bounded by the dimension of its span in (3).

The affine divisor map takes (3) to affine divisors modulo \(l\)
and modulo divisors of rational functions of \(t\). An element
in its kernel can be adjusted by a rational function of \(t\) so
its affine divisor is \(lD\). The class \([D]\) then belongs to
\(\operatorname{Pic}(C^\circ)[l]\). Conversely, the group of
functions with divisible affine divisor modulo affine units and
\(l\)th powers is isomorphic to that Picard torsion by the map
\(\alpha\mapsto[\operatorname{div}(\alpha)/l]\). Consequently the
kernel on (3) has dimension at most \(r_\circ\).

By VF1, the image of the row span is supported over at most one
linear \(g\). A one-point fiber contributes no image because its
divisor is the divisor of \(g\) with coefficient one or two.
In a split fiber the two exponents are taken modulo the common
diagonal and yield at most one dimension. If that dimension is
nonzero, a principal divisor modulo \(l\) has the form
\(aP_++bP_-\) with \(a-b\ne0\bmod l\). Since \(P_++P_-\) is
principal on the affine curve, its image in
\(\operatorname{Pic}(C^\circ)/l\) forces (2).
For a full progression VF1 removes the entire image.

For the genus-zero Pell chart of
[PI3](PELL_PERIOD_RING_CLASS_INERTIA.md) and
[PM3](PELL_NORM_SECTIONS_MORDELL_RANK.md), use
\(y=B\sqrt{n(B^2n+2)}\) and \(A=1+B^2n\). Its affine ring is
\(\mathbb Q[U^{\pm1}]\) by
\[
       n=\frac{U+U^{-1}-2}{2B^2},\qquad y=\frac{U-U^{-1}}2.
\]
Thus \(r_\circ=0\) and every fixed collection of rational functions
can give at most one ordinary direction on any infinite integral
set. The row \(Bn+c+\sqrt{n(B^2n+2)}=(U+Bc-1)/B\) has one
nontrivial split fiber at \(n=-c^2/(2(Bc-1))\), and the genuine
ideal of PI3 shows the bound is attained. The statement applies
to any degree of fixed row on this conic.

For [CP1](QUARTIC_PELL_TWO_CLASS_DIRECTIONS.md) use its
integral progression parameter and \(l=3\). Its two known
ordinary directions exhaust what any fixed list of rows can give
on a common full progression. A thin infinite set can have a
third direction from such rows only through a rational split
fiber satisfying (2).

## Realizing the extra direction by a ramified change of parameter

**Corollary VF3.** In the Pell chart of VF2 suppose the fiber at
\(t=a\in\mathbb Q\) is split, with points \(P_+,P_-\), and (2)
holds. There is a \(\lambda\in\mathbb Q^\times\) such that on the
smooth curve
\[
                   C':Y^2=F(a+\lambda q^l)                \tag{4}
\]
one new rational \(l\)-torsion divisor class is independent of
the pullback of \(J(\mathbb Q)[l]\). The pullback is injective.
If \(l\nmid m\), the new infinity class of order \(l\) is also
independent of these classes and gives a second new direction.
The added finite-fiber direction survives quotienting by the
infinity difference, independently of the pullback of
\(\operatorname{Pic}(C^\circ)[l]\).

Choose \(f\in\mathbb Q(C)^\times\) with
\(\operatorname{div}_{\rm aff}(f)=P_+-lD\). This is possible by
(2), with \(D\) an affine divisor over \(\mathbb Q\).
Then
\[
                 N(f)=c(t-a)G(t)^l
\]
for some \(c\in\mathbb Q^\times\), \(G\in\mathbb Q(t)^\times\).
Take \(\lambda=c^{-1}\). The norm becomes an \(l\)th power
under (4). Because \(F(a)\ne0\), the new polynomial is
squarefree. The map \(C'\to C\) is the extension by a root of
\(q^l=(t-a)/\lambda\); it has degree \(l\), and ramifies
totally above both points of \(t=a\) and both infinities.
In particular the entire pullback divisor of \(f\) is divisible
by \(l\), defining its rational Jacobian class on \(C'\).

The kernel of restriction on classes of functions modulo
\(l\)th powers is the line generated by \((t-a)/\lambda\).
Indeed, adjoining \(\zeta_l\) reduces this to ordinary Kummer
theory, and passage to that extension is injective on power
classes by its degree prime to \(l\). The extension stays
degree \(l\) there by a simple valuation at \(P_+\).
If a pullback \(l\)-torsion class becomes zero, its defining
function becomes an \(l\)th power up to a rational constant.
The kernel line has nonzero orders at both \(P_+,P_-\);
the defining function has orders zero modulo \(l\) there,
so its original class was already zero. A similar comparison
for \(f\) gives the unequal orders \((1,0)\), which neither
the kernel line nor any old divisible divisor can supply.
This proves its independence. The same comparison applies
on the affine chart even after multiplication by powers of
the specialized infinity unit. That unit is the pullback of
\(U\), as follows from the degree proof of VF2 applied to (4).

On the projective chart the pullback \(U\) has divisor
\(lm(\infty'_+-\infty'_-)\), up to sign. It gives a new class
of order \(l\) when \(l\nmid m\). Any relation to old classes
and \(f\) first forces the coefficient of \(f\) and the kernel
line to vanish by the orders over \(t=a\). Then the base
orders of \(U\) at infinity, which are \(\pm m\), force its
coefficient to vanish.

VF1 gives virtuality of the pulled-back rows outside a fixed
finite set of rational primes. An ordinary specialization
result with an extra direction still requires virtuality at
those remaining primes and independence after specialization
and quotienting by the unit line. Corollary VF3 supplies the
geometric condition and the parameter change for that problem.

## An arithmetic realization and its norm-only deformation

**Corollary VF4.** Take the base Pell curve \(y^2=t^4+1\) with
\(U=t^2+y\) and the split fiber \(t=0\).
Its rational Jacobian three-torsion is zero. After \(t=q^3\)
the genus-five curve \(R^2=q^{12}+1\) has rational Jacobian
three-torsion exactly \(C_3^2\), generated by its infinity
direction and the divisor class of \(V=2(U-1)\) divided by
three. On the integral slice
\[
                q\equiv1\bmod 2,3,7,\qquad q\equiv0\bmod13,
\]
the fields \(K=\mathbb Q(\sqrt{q^{12}+1})\), \(q>0\), have an
integral norm-minus-one unit \(U=q^6+R\) and an independent
ordinary class direction supplied by \(V\). Three is inert,
\(U\) is noncubic at the completion over three, and \(V\) is
cubic there. The slice yields \(\gg X^{1/12}\) distinct fields
of discriminant at most \(X\).
In the even-\(q\) deformation the same norm is a cube, but
\(V\) fails virtuality at both places over two.

On the base the affine divisor of \(U-1\) is \(2P_+\), where
\(P_+=(0,1)\), since \(N(U-1)=-2t^2\). Hence \([P_+]\) is
divisible by three in the affine Picard group, as (2) requires.
The choice \(f=V^2\) yields
\(N(f)=64t^4=(4t)^3t\), so the change of parameter in VF3 can
take \(\lambda=1\). The base has good reduction at five and
four points there, excluding rational three-torsion.
On the new curve the two independent classes come from VF3,
since \(m=2\) is prime to three. The point counts over
\(\mathbb F_{5^i}\), \(1\le i\le5\), give Frobenius trace sums
\((2,-14,-22,58,82)\). Newton's identities and reciprocity give
the Jacobian group order \(2880=2^6\cdot3^2\cdot5\).
Good reduction at five gives the matching upper bound.

At every specialization \(N(V)=(-2q^2)^3\).
An odd split prime cannot divide both conjugate values of \(V\),
as this would give \(U\equiv\bar U\equiv1\) at one place and
contradict \(U\bar U=-1\). At odd nonsplit primes the cube norm
gives the desired valuation. At two, odd \(q\) gives
\(q^{12}+1\equiv2\bmod8\), so the unique ramified prime
has valuation divisible by three by the cube norm.
In the even deformation the field splits at two. The two
branches of \(R\) are congruent to \(1,-1\bmod8\).
In the first one
\(v_2(R-1)=2v_2(q^6)-1>v_2(q^6)\), and in the second
\(v_2(R-1)=1\). Thus the valuations of \(V\) are
\(1+6v_2(q),2\), neither divisible by three.

At seven take \(R\equiv3\) and the cube root of unity \(2\).
The residue exponents of \(U,V\) are \(1,0\).
At thirteen take \(R\equiv-1\) and cube root \(3\).
The exponents are \(0,2\). These tests show that \(U,V\) are
independent modulo cubes in \(K^\times\). The positive and
negative real units modulo cubes span only one line, so \(V\)
gives a nonzero class in \(\operatorname{Cl}(K)[3]\).

At three, the sixth power \(q^6\equiv1\bmod9\) gives the inert
quadratic completion and \(R^2\equiv2\bmod9\). A unit in the
unramified quadratic completion is a cube exactly when its
eighth power is \(1\bmod9\), by the Teichmüller decomposition
and the cube map on principal units. There,
\(U^8\equiv1+3R\) and \(V^8\equiv1\bmod9\).
Finally [PF1](PELL_UNIT_PRIMITIVITY_FROM_TRACE_EQUATIONS.md)
applies to \(2q^6\) with negative norm. Almost all \(q\) give
the fundamental positive unit and therefore distinct fields,
with discriminant \(O(q^{12})\).

## The fixed quartic has no rational split fiber

**Corollary VF5.** The fixed \(w=1/48\) curve of CP1 has precisely
its two points at infinity over \(\mathbb Q\). Consequently the
two ordinary directions are the most that any fixed list of rows
on this curve can supply on any infinite integral set, including
thin sets. In the reflected family at \(w=5/24\), the split fibers
from the two roots of \(H\) exist, but neither fiber meets (2)
for \(l=3\).

Write \(P(t)=a t^4+b t^2+c\) and \(U=P(t)+y\).
The involution \(t\mapsto-t\) has the elliptic quotient
\[
 E:\quad Y^2=X^3+(b^2-4ac)X^2-4a^2X,\qquad
 X=2aU,\qquad Y=2a(2at^2+b)U .
 \tag{5}
\]
Indeed \(U-U^{-1}=2P(t)\); completing the square in \(t^2\)
gives (5). The point over the pole of \(U\) maps to the identity
and the point over its zero to \((0,0)\). A rational affine
point must give \(X\ne0\) and \(t^2=(Y/X-b)/(2a)\).

For \(w=1/48\), a rational Weierstrass change gives the model
\[
     E_0:\quad y^2+xy=x^3-81381060x-282436946928.
\]
The change from this model back to (5) has Weierstrass parameters
\((u,r,s,t)=(8/1153,327168/1329409,4/1153,0)\).
The exact [PARI rank bounds](sources/pari-ellrank-rank-bounds.html)
are \([0,0,2,[]]\) on \(E_0\), from the 2-descent and its
Cassels pairing. Both rank bounds agree. The counts at good
primes seven and eleven are \(6,12\), so the rational torsion
order divides six. The six rational points are
\[
 O,\quad (-5112,2556),\quad (-5256,-7812),\quad(-5256,13068),
 \quad(15624,-1511172),\quad(15624,1495548).
\]
The first two are the images of the infinities. On the last
four points \((Y/X-b)/(2a)\) is one of
\[
                   -\frac{4613}{442368},\qquad
                    \frac{593923}{49152}.
\]
The negative value and the positive value of odd three-adic
valuation are both nonsquares. This proves the rational-point
assertion and rules out the extra split fiber in VF2.
The rank input is the unconditional upper bound from the
descent. As a separate computational check, the \(L\)-series
central value has a strict positive interval from its first
1800 coefficients and a geometric tail. The latter check
uses the squarefree conductor and its positive root number.

For \(w=5/24\) the two roots of \(H\) from CP1 are \(269/32\)
and \(-181/96\). They give \(y=\pm16409/1280\) and
\(y=\pm181/180\), respectively. The induced map on affine
Picard groups takes \([P]\) to the quotient elliptic point
modulo \((0,0)\). This boundary point has order two and is
already divisible by three. Thus (2) would imply divisibility
of the quotient point in \(E(\mathbb Q)\). At the first fiber,
both signs in (5) give points outside \(3E(\mathbb F_7)\)
and at the second both signs give points outside
\(3E(\mathbb F_{11})\). The two good reduced groups have order
12. Reflection supplies the same quotient points for the
opposite values of \(t\). This obstruction concerns the split
fibers arising from the norm-section roots at this parameter.

## Verification

The [exact checker](verify_pell_virtual_affine_fibers.py) checks
the quadratic Laurent coordinate, the prime-path divisor
pairing, the two residue tests and exact finite-field traces
for VF4, and the complete elliptic quotient controls for VF5. The finite-support proof uses the stated effective
height theorem and does not enumerate superelliptic solutions.
