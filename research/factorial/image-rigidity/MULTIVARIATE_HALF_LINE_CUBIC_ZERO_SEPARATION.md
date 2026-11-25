# Zero separation for the half-line cubic model

For real \(T\) and complex \(H\) let
\[
 F(T,H)=\int_0^\infty e^{-y^3/3-Ty^2/2+Hy}\,dy.
\]
This is the limiting function in the
[moving one-private cubic window](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md).
At a simple zero, write \(H_T\) for the derivative of its
local zero branch with respect to \(T\).

**Theorem.** On each compact real temperature interval, only
finitely many \(T\) can have a multiple zero, two distinct
zeros at the same imaginary height, or a simple zero with
\(2H-T H_T=0\). For all sufficiently large positive \(T\),
none of these three events occurs. The bounds here are
unevaluated.

In particular the model comparison establishes the bounded
moving one-private consecutive exclusion when the limiting
temperature is outside a locally finite exceptional set.
Membership at individual bounded temperatures is a separate
question.

## An exact Airy decomposition

Put \(a=T/2\), \(Z=H+a^2\), and \(\omega=e^{2\pi i/3}\).
Shifting the integration variable by \(a\) gives the
exponential factor \(e^{-aH-2a^3/3}\) times
\(\int_a^\infty e^{-u^3/3+Zu}\,du\).
Subtract the integral from \(a\) to infinity in direction
\(\omega\). The difference, after the change of variable
\(v=-\omega u\) and reversal of orientation, is the ordinary
Airy contour. The resulting identity is
\[
 F(T,H)=\omega F(T\omega^2,\omega H)
      +2\pi i\omega^2 e^{-aH-2a^3/3}\operatorname{Ai}(\omega^2 Z).
                                                               \tag{1}
\]
The contour identity is valid for complex \(T,H\) by cubic
decay at the ends. For \(T=0\) it is the usual Scorer
connection formula. The Airy contour and its leading
sectorial expansion used here are the classical formulas
[9.5.4](https://dlmf.nist.gov/9.5.E4) and
[9.7.5](https://dlmf.nist.gov/9.7.E5) with the uniform
sector bounds in [9.7(iv)](https://dlmf.nist.gov/9.7.iv).
For the zero-temperature connection and zero asymptotics
see [Gil, Segura and Temme](https://ir.cwi.nl/pub/10698/10698D.pdf).

For bounded complex \(T\), choose a small fixed angle
\(\eta>0\). Uniformly for
\(\pi/3+\eta\le\arg H\le\pi\), rotate the original integration
ray to direction \(\pi/6-\eta/2\). The cubic has strictly
negative real leading term there, while the linear term is
at most \(-c|H|y\). The residual factor after the linear
exponential is one plus \(O(y^2)\) near zero and bounded
on the ray. Thus \(F(T,H)=-H^{-1}(1+O(|H|^{-2}))\) there.
On \(0\le\arg H\le\pi/3+\eta\), the first term of (1)
has this same expansion: \(\Re(\omega H)\le-c|H|\), and
the original integration ray for \(F(T\omega^2,\omega H)\)
suffices. For the second term of (1), let
\[
 S=\tfrac23 Z^{3/2}-aH-\tfrac23 a^3 .
]
On this sector and for large \(|H|\), the principal
determinations give
\((\omega^2 Z)^{3/2}=-Z^{3/2}\).
The second term is
\(\sqrt\pi Z^{-1/4} e^S(1+O(|Z|^{-3/2}))\).
In particular it dominates for \(0\le\arg H\le\pi/3-\eta\),
because \(\Re S\ge c|H|^{3/2}-O(|H|)\).

Consequently the large upper zeros for bounded real \(T\)
lie in an arbitrarily small fixed sector about \(\arg H=\pi/3\).
There (1) writes a zero as
\[
 \Xi(T,H)\in2\pi i\mathbb Z,\qquad
 \Xi=S+\log(\sqrt\pi H Z^{-1/4})
                  +\log(1+R_A)-\log(1+R_E),
 \tag{2}
\]
with \(R_A=O(|Z|^{-3/2})\), \(R_E=O(|H|^{-2})\), both analytic.
Here the first term of (1) is \(-H^{-1}(1+R_E)\) and
the second is \(\sqrt\pi Z^{-1/4}e^S(1+R_A)\).
By taking the expansion on a slightly wider sector and a
larger temperature disk, Cauchy's estimate gives
\[
 \Xi_H=\sqrt H+O(1),\qquad
 \Xi_T=-H/2+O(\sqrt{|H|}).
 \tag{3}
\]
Thus every large zero in a compact temperature band is
simple. Since \(\Re\Xi_H>0\), a horizontal segment in
this sector can contain at most one zero: the real part
of \(\Xi\) strictly increases, whereas it is zero at any
zero of \(F\). In particular every sufficiently large
imaginary height has at most one zero. Also
\(H_T=-\Xi_T/\Xi_H=O(\sqrt{|H|})\) and
\(2H-T H_T\ne0\) for sufficiently large \(|H|\).
Conjugation gives the lower-half-plane assertions.

## All zeros at large positive temperature

For \(T\ge0\), no zero has \(\Re H\le0\).
For \(\Im H>0\) and \(\Re H\le0\) the real density
\(y\mapsto e^{-y^3/3-Ty^2/2+(\Re H)y}\) is strictly
decreasing. Pairing consecutive positive and negative
half periods of its sine transform gives a strictly
positive imaginary part of \(F\). Conjugation and
positivity on the real axis give the rest of the claim.
The same argument for \(T<0\) excludes
\(\Re H\le -T^2/4\), where the density is also strictly
decreasing.

Take \(T\ge1\), \(\Re H\ge0\), \(\Im H\ge0\), and
\(|H|\ge L\sqrt T\), for a sufficiently large fixed \(L\).
To estimate the endpoint term in (1), rotate its
integration ray from \(\omega\) to direction \(7\pi/9\).
The real parts of both the cubic and the quadratic
exponent terms are nonpositive on that ray. The real
linear term is at most \(-c_0|H|y\) throughout the
upper-right quadrant of \(H\). Using
\(|e^v-1|\le |v|\) for \(\Re v\le0\) now gives
\[
 R_E=O(T/|H|^2+|H|^{-3}).
]
The rotation is valid for all finite \(T,H\) by cubic
decay at infinity. Here \(|Z|\asymp a^2+|H|\) and
\(0\le\arg Z\le\pi/2\), so the same Airy expansion in
(2) gives \(R_A=O(|Z|^{-3/2})\) uniformly. These estimates
persist on fixed relative complex neighborhoods of
positive \(T\) and of \(H\) inside the quadrant away
from its boundary rays. For the endpoint ray choose
the relative angle of \(T\) small enough that the
quadratic still has nonpositive real exponent.
For the Airy argument, \(Z\) remains in a fixed larger
sector about the right half-plane and its modulus
has the same comparison. Cauchy's estimate therefore
bounds logarithmic error derivatives after multiplication
by \(|H|\) or \(T\) by
\(O(T/|H|^2+|H|^{-3}+|Z|^{-3/2})\) on closed subcones.

Put \(V=|H|/(a+\sqrt{|H|})^{1/2}\).
The logarithmic amplitude in (2) has real part
\(\log V+O(1)\). The saddle exponent is also
\[
 S=H^2\int_0^1\frac{1-s}{2\sqrt{a^2+sH}}\,ds .
]
If \(\theta=\arg H\), the integrand after multiplication
by \(H^2\) has argument between \(3\theta/2\) and \(2\theta\).
It follows that \(\Re S\ge c_1 V^2\) for
\(0\le\theta\le\pi/6\), and \(\Re S\le-c_1 V^2\)
for \(5\pi/12\le\theta\le\pi/2\).
For example, the lower bound on the absolute scalar
integral uses
\(|a^2+sH|^{1/2}\le(a^2+|H|)^{1/2}\).
With \(L\) large, \(V\) is uniformly large on
\(|H|\ge L\sqrt T\). Hence (2) forces every zero in
this exterior region into the open cone
\(\pi/6<\arg H<5\pi/12\).

On the closed cone,
\[
 S_H=P=\sqrt Z-a=H/(\sqrt Z+a),\qquad
 \Re P\ge c_2 |H|/(a+\sqrt{|H|}).
]
This follows from the argument of \(\sqrt Z+a\) lying
between zero and \(\theta/2\).
The remaining derivative terms of \(\Xi_H\) are
\(O(1/|H|)\). Thus \(\Re\Xi_H>0\) on the cone,
after increasing \(L\). This proves simplicity and
distinct imaginary heights for exterior zeros by the
horizontal segment test. For the order derivative,
\(S_a=-P^2\) gives
\[
 (2H\partial_H+T\partial_T)S=P^2(2\sqrt Z+a).
]
Its modulus is comparable to \(V^2\).
The same operator on the logarithmic amplitude in (2)
gives \(3/2\) and on the remainders gives a bounded term,
by the relative complex neighborhoods above. Hence
\(2H\Xi_H+T\Xi_T\ne0\) there, proving
\(2H-T H_T\ne0\) at exterior zeros.

For the inner fields, put
\(H=\sqrt T\,\zeta\) and \(s=T^{-3/2}\). Then
\[
 \sqrt T F(T,\sqrt T\,\zeta)=
         \int_0^\infty e^{-u^2/2-su^3/3+\zeta u}\,du,
 \qquad E(\zeta)=\int_0^\infty e^{-u^2/2+\zeta u}\,du .
]
The integral and its field derivative converge uniformly
on compact \(\zeta\)-sets to \(E,E'\), and the \(s\)-derivative
is bounded there for positive \(s\). The half-Gaussian
function satisfies \(E'=\zeta E+1\), so its zeros
are simple. They lie in \(\Re\zeta>0\): for \(\Re p>0\),
the Fourier representation of the Gaussian gives
\[
 E(-p)=\frac1{\sqrt{2\pi}}
             \int_{\mathbb R}\frac{e^{-v^2/2}}{p-iv}\,dv ,
]
which has positive real part. On the imaginary axis
the cosine integral of the Gaussian is positive.
Also \(E(\zeta)+E(-\zeta)=\sqrt{2\pi}e^{\zeta^2/2}\).
For \(\Re\zeta>0\) the horizontal derivative of
\(\log|E(-\zeta)|-\Re(\zeta^2)/2\) is
\(-\Re(1/E(-\zeta))<0\). The equality of these two
moduli at a zero therefore permits at most one
zero at each imaginary height.

Rouché's theorem and the implicit function theorem
now give simple inner zeros, with distinct heights
on any fixed compact \(\zeta\)-set for large \(T\).
Each belongs to a branch \(H=\sqrt T\,\zeta(s)\) with
bounded \(\zeta'(s)\). Since the limiting zeros stay
away from zero on such a set,
\[
       2H-T H_T=	frac32\sqrt T\,\zeta(s)
                         +\tfrac32 T^{-1}\zeta'(s)\ne0 .
]
If an exterior zero shares a height with a zero in
\(|H|\le L\sqrt T\), the exterior cone implies
\(|H|\le2L\sqrt T\) for that exterior zero as well.
Thus the compact inner comparison with enlarged radius
covers all mixed pairs. This proves the three assertions
for all sufficiently large positive \(T\).

## Isolated exceptional temperatures

Differentiation and integration by parts give
\[
 F_T=-F_{HH}/2,\qquad F_{HH}+T F_H-HF=1.
 \tag{4}
\]
The integral is entire in the two complex variables.
At a multiple zero its second field derivative is
one, so the zero has multiplicity exactly two, and
\(F_T=-1/2\). The equations \(F=F_H=0\) have
invertible Jacobian in \((T,H)\) there. Hence their
solutions are isolated and finite on compact boxes.
At a simple zero, the order derivative vanishes
exactly when
\[
                    (4H+T^2) F_H-T=0 .
 \tag{5}
\]
This equation cannot hold identically on a local
simple zero branch. Near nonzero \(T\) the identity
\(2H-T H_T=0\) would force \(H=C T^2\); then
\(F(T,C T^2)\) would vanish identically by analyticity,
contrary to \(F(0,0)>0\). At \(T=0\) it would
require \(H=0\), also impossible.
The joint equations \(F=0\) and (5) can meet a
multiple zero only at \(T=0\). There \(F_T\ne0\),
so \(F=0\) is a smooth curve parametrized by \(H\)
and \(dT/dH=0\) at the multiple zero. The derivative
of the left side of (5) along that curve is then
\(4H\ne0\). Thus the zero-speed events are isolated,
including near a multiple zero. By the uniform
large-field result (3), both exceptional sets have
finitely many temperatures on every compact real interval.

To treat same-height events, first exclude any persistent
pair of distinct simple branches at equal imaginary height
on an open real temperature interval. Such an identity
continues along both simple branches toward increasing
\(T\) as long as they remain simple. No finite escape in
\(H\) is possible, by the derivative bound from (3) on
each compact temperature interval. If the continuation
stopped at a finite \(T_0\), each branch would have a unique
limit by compactness and the discreteness of zeros at
\(T_0\). Indeed its cluster set is the intersection of
the nested compact connected closures of the tails.
At least one limit would be a double zero. If both
were simple, either they would be distinct and the
identity would continue, or they would be the same
simple root, contrary to uniqueness of its analytic
branch nearby.

At a double zero \((T_0,H_0)\), put \(t=T-T_0\) and
\(h=H-H_0\). Since \(F_T=-1/2\), the zero set is an
analytic curve with \(t=h^2+O(h^3)\). Its two local
branches in \(s^2=t\) are given by taking the two signs
of \(s\) in
\[
 H=H_0+s-T_0 s^2/3
                  +(H_0/12+T_0^2/36)s^3+O(s^4).
 \tag{6}
\]
For verification, (4) gives the coefficients
\[
 \begin{split}
 F={}&(h^2-t)/2-T_0 h^3/6+T_0 t h/2\
    &+(H_0+T_0^2)(h^4/24-t h^2/4+t^2/8)+O(s^5)
 \end{split}
\]
when \(h=O(s)\) and \(t=s^2\), and substitution yields (6).
Approaching a proposed stopping temperature from below,
write \(s=\pm i\sqrt{T_0-T}\). The imaginary expansions
then distinguish a double branch from a simple one
by their linear terms in \(\sqrt{T_0-T}\). For two
branches tending to the same double zero the signs
are opposite. For branches tending to distinct double
zeros at the same limiting height, the signs would
have to agree, but the cubic coefficients of their
imaginary parts differ by the difference of the
real parts of the centers divided by \(12\).
In each case persistent equality is impossible.
Thus a persistent pair would continue to arbitrarily
large positive \(T\), contradicting the result above.

Finally, on a compact real temperature interval the
same-height events have bounded field by (3) and the
large-field angular localization. If their temperatures
accumulated, both zeros would limit to zeros at the
accumulation temperature. Use the implicit function
theorem at a simple limit and the local square-root
branches (6) at a double limit. On one side of that
temperature a subsequence is covered by two fixed
branches, each analytic in the real square-root
coordinate. Their imaginary difference is real
analytic in that coordinate and vanishes on a
sequence accumulating at zero. It must vanish
identically, giving the persistent pair just ruled
out. This proves local finiteness and completes
the theorem. \(\square\)
