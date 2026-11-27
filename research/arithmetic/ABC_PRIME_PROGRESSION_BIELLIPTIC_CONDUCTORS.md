# Bielliptic projections and conductors for ordinary prime progressions

The ordinary \((1,1,1)\) cubic prime progression has two
fixed elliptic quotients. Their Kummer marks are the
quotient and product of the consecutive norm classes.
Only the quotient is a local cube at \(3\). Both still
give independent elliptic twist constructions on the
same quadratic field, with fixed conductor laws and
opposite finite \(3\)-Selmer parities.

## Identifying the two elliptic directions

**Theorem PQ1.** The genus-two progression curve
\[
 C:\quad y^2=-3t^6+24t^4+24t^2+4
\]
has degree-two projections to the fixed curves
\[
 E_1:\ Y^2=X^3+(6X+14)^2,\qquad
 E_2:\ Y^2=X^3+(6X+28)^2.
\]
They generate the two rational elliptic factors of its
Jacobian, and each quadratic twist by the displayed
sextic has generic rank one.
If
\[
 \alpha_1=(t^3-6t^2-2+y)/2,\qquad
 \alpha_2=(t^3+6t^2+2+y)/2,\qquad
 \gamma=\bar\alpha_1\alpha_2,\quad
 \eta=\alpha_1\alpha_2,
\]
then the marked cubic classes on \(E_1,E_2\) pull back
to \(\gamma^{-1},\eta^{-1}\), respectively.

The explicit coordinates are
\[
 (X_1,Y_1)=(-3t^2-4,-3y),\qquad
 (X_2,Y_2)=(4(t^{-2}-1),4y/t^3).
\]
The abscissa maps have degree two. The first projection
is fixed by \(t\mapsto-t\) and the second is negated.
Thus the two elliptic factors span the Jacobian up to
isogeny, as in BT1. Their traces at \(13\) are \(2,-4\),
so the factors are nonisogenous over \(\mathbb Q\).
This gives the generic rank assertion. Direct
multiplication yields
\[
 \begin{split}
 (Y_1+6X_1+14)\gamma&=(y+5t^2+2)^3,\\
 (Y_2+6X_2+28)\eta&=-8(t^2-1)^3.
 \end{split}
\]
Thus the two arithmetic directions of the sign
decomposition are explicit, including their cube factors.

## Specialization to genuine prime progressions

**Theorem PQ2.** In the ordinary coherent progression
family
[PP1–PP3](ABC_COHERENT_PRIME_PROGRESSION_NORMS.md),
impose \(h\equiv0\bmod126\) and \(p\equiv1\bmod7\),
retaining the two auxiliary cube residue tests of PP1.
On either sign interval of the quadratic radicand, let
\(D\) be the resulting fundamental discriminant.
The \(D\)-twists of \(E_1,E_2\) have positive rank, with
conductors \(84D^2,21D^2\), respectively.
The first marked class is a local cube at \(3\), and the
second is a local noncube there. The real or imaginary
mirror of \(K=\mathbb Q(\sqrt D)\) has the unramified cyclic
cubic extension defined by the first mark, split at \(3\).

The two PP1 residue tests give independent classes
\(\alpha_1,\alpha_2\) also on an interval for the real
sign, since they are tests on local power classes.
Here the displayed rational functions on \(C\) are the
homogeneous prime-triangle functions divided by \(h^3\).
Thus both elliptic directions are noncube by PQ1.
Their rational twist evaluations are the general
[Cohen–Pazuki Kummer evaluations](sources/cohen-pazuki-2009-three-isogeny.pdf)
up to cubes and signs, as in BT2. The \(3\)-division
polynomials on the original curves factor as
\[
 \begin{split}
 3X(X^3+48X^2+336X+784),\\
 3X(X+28)(X^2+20X+112).
 \end{split}
\]
The first cubic is irreducible, and on the second curve
the order-three points outside the marked rational line
are defined only over \(\mathbb Q(\sqrt{-3})\) among
quadratic fields. In particular this exceptional
quadratic field is excluded by the two residue tests.
There is no anti-invariant \(3\)-torsion to account for
the noncubic evaluations, so both displayed twist
points have infinite order.

At \(3\), let \(r=\sqrt{-3}\) in the completion of \(K\).
The square root of the radicand divided by \(p^3\)
can be chosen to be \(r\) times an element of
\(1+3^4\mathbb Z_3\), since \(h\equiv0\bmod9\).
The two rational traces divided by \(p^3\) are
\(1\mp6h/p\mp2(h/p)^3\).
Both normalized norm classes are consequently
\((1+r)/2\) modulo local cubes.
Indeed the quotient by this unit differs from one
by an element of valuation greater than two, with
the valuation normalized on \(\mathbb Q_3\).
The unit \((1+r)/2\) is a primitive sixth root of
unity and is noncubic in \(\mathbb Q_3(r)\), which
does not contain \(\mu_9\).
Thus precisely the \(\gamma\)-line in the marked
two-dimensional subspace consists of local cubes.
This is the mirror direction of
[RF3](ABC_COHERENT_NORM_CUBIC_REFLECTION.md#the-fundamental-discriminant-on-the-fermat-slice),
now on the fixed first elliptic quotient.

At \(2\) the radicand is \(5\bmod8\), and at \(3\)
its squareclass is \(-3\). The field has squareclass
one at \(7\), using \(h\equiv0\bmod7\) and
\((-3/7)=1\). In particular \(D\) is odd and \(v_3(D)=1\).
The base conductors are \(756=2^2\cdot3^3\cdot7\)
and \(189=3^3\cdot7\), and the fixed local classes
\(5,-3,1\) at \(2,3,7\) preserve exactly those
exponents. All other primes dividing \(D\) give
exponent two by good reduction of the base curves
and odd ramified twisting. Dividing the base
conductors by \(3^2\) and multiplying by \(D^2\)
gives the two formulas.

## Exact local Selmer parity on both signs

**Proposition PQ3.** For negative \(D\) in PQ2 the
finite \(3\)-Selmer dimensions on the two twists have
parities even, odd. For positive \(D\) the parities
are odd, even. In either case the even dimension is
at least two and the odd dimension at least one.
The \(D\)-twist of \(E_3:Y^2=X^3+(12X+4)^2\), which is
\(3\)-isogenous to the second twist, has conductor
\(21D^2\) and finite \(3\)-Selmer dimension at least
three for negative \(D\) and at least two for
positive \(D\).

Use the \(3\)-isogenies with \(X=0\) as the marked
subgroup. For the model
\(Y^2=X^3+D(aX+b)^2\), the normalized isogenous
model in Cohen–Pazuki has
\(\widehat D=-3D,\ \widehat a=a,\
\widehat b=(27b-4Da^3)/9\).
The rational twist models here have \(a=6,\ b=b_0D\)
for \(b_0=14,28\), after the usual quadratic twist
coordinate change. The isogeny in those coordinates
pulls back the original model differential with
multiplier one. The local ratio on minimal models
is the target Tamagawa number divided by the source
Tamagawa number, times the inverse absolute value
of the differential multiplier.
At \(2,3,7\) the required local ratios are
\[
 \begin{array}{c|rrr}
 &2&3&7\\\hline
 E_1^{(D)}&3&9&1/3\\
 E_2^{(D)}&1&1&1/3 .
 \end{array}
\]
They are computed on \(D=5,-3,1\), respectively,
which are the local squareclasses throughout.
In detail, the source and target Tamagawa pairs
are \((1,3),(1,3),(3,1)\) for the first curve and
\((1,1),(3,1),(3,1)\) for the second.
The minimal differential multiplier at \(3\) has
valuation one in both cases. At \(2,7\) it is a
unit. At the other primes the ratio is one, since
both good-reduction factors and the prime-to-three
component groups on odd ramified twists are unchanged
under a \(3\)-isogeny.

The real ratio is one for negative \(D\) and \(1/3\)
for positive \(D\), from the real kernel and
surjectivity on real points for an odd isogeny.
Thus the global ratios are \(9,1/3\) in the negative
case and \(3,1/9\) in the positive case.
In both cases the twists have trivial rational
\(3\)-torsion. The general \(3\)-isogeny parity and ratio formulas give
the finite Selmer parity as the valuation of the
global ratio modulo two.
[Bhargava–Klagsbrun–Lemke Oliver–Shnidman, Section 9](sources/bhargava-klagsbrun-lemke-oliver-shnidman-three-isogeny.pdf)
The displayed non-torsion points then give the
dimension lower bounds.

On \(E_2\) the other rational subgroup has
\(X=-28\). Recenter the \(D\)-twist at \(X=-28D\).
Its model becomes \(Y^2=X^3-3D(4X-28D)^2\),
so the isogenous model from this subgroup is isomorphic
to \(E_3^{(D)}\). The minimal differential multiplier
is a unit at \(3\) for this isogeny. The local ratios
from \(E_2^{(D)}\) to \(E_3^{(D)}\) at \(2,3,7\)
are \(1,1/3,1/3\), respectively. The real ratio is
\(1/3\) for negative \(D\) and one for positive \(D\),
since this marked kernel has character \(-3D\).
Thus the dual isogeny from \(E_3^{(D)}\) has global
ratio \(27\) or \(9\), respectively. The rational
kernels on both sides are trivial, as \(D\) and
\(-3D\) are nonsquares. By Cassels' ratio formula the
dual-isogeny Selmer group on \(E_3^{(D)}\) has at
least that size. It injects into its \(3\)-Selmer group
by the multiplication-by-three exact sequence.
This gives the dimension bounds for \(E_3^{(D)}\).
Its conductor and positive rank follow by isogeny.

## Counts with conductor and Selmer conditions

**Corollary PQ4.** For either sign there are
\(\gg X^{1/3}/(\log X)^B\) distinct fields \(K\) in
PQ2 of absolute discriminant at most \(X\), where
\(B\) is fixed. The number of the corresponding
conductor pairs with larger conductor at most \(Z\)
is \(\gg Z^{1/6}/(\log Z)^B\).
The imaginary base fields have two independent
ideal classes of order three as in PP1.

The added fixed congruences are compatible with the
auxiliary prime conditions and the three linear forms
are locally admissible. Choose a fixed compact
interval of \(p/h>1\) away from the roots of the
sextic for the requested sign. Apply the PP2–PP3
prime-form count and twist multiplicity bound.
The branch splitting field and its fixed
two-torsion control are independent of the
quadratic sign, so the same argument works on
both intervals. The conductor pair is injective
in \(D\) on each sign.
The [exact checker](verify_prime_progression_bielliptic_conductors.py)
replays the maps and cube identities, the local
Tamagawa and differential inputs, and both signed
prime-progression fixtures.
