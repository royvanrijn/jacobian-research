# Exact contact/source filtration and three surviving height-eight strata

The subsequent [PF2RULED16 chart-leading obstruction](POLE_HEIGHT_EIGHT_EXCLUSION.md)
excludes all three contact strata below. They are historical necessary
controls, not current realization targets; the full filtration identities
and jet ranks remain valid and are used by that exclusion.

14 November 2026. PF2RULED15 works on the actual marked Wright inclusion,
conditional on a finite quartic Keller target algebra. It narrows the
contact-only height-eight frontier of PF2RULED14 as follows:

| Pole pair | Remaining common-symbol divisor | Ordinary source degrees |
| --- | --- | --- |
| (6,8) | div(sigma)=4p+2q, p!=q, sigma a sextic | (72,108) |
| (7,8) | div(tau)=p+2q or3q, p!=q, tau a cubic | (72,108) |

The former contact multiplicities3,4,6 for(6,8) and1 for(7,8) are
excluded. All three remaining strata have boundary restriction degrees
(18,27) in the fixed parameter v=D. The cusp multiplier is exactly
linear, with slope20/27, and its nonzero value at the cusp is determined
by the first two leading coefficients of Q0. Neither pole pair, the
all-height pole-four problem nor actual compatible completion is closed.

## The two filtrations on the actual ring

Use the source coordinates X,Y and the fixed Cox slice from
[PF2WRIGHT1](WRIGHT_MARKED_MODEL.md):

 alpha=1+X^2Y, t=1+X^2+X^4Y,
 beta=-1-2X^2Y-X^4Y^2-Y, s=X,
 alpha*t+beta*s^2=1.

By [PF2COORD1](SOURCE_COORDINATES.md), a complete basis of V_m is

 alpha^(m-j)*beta^j*s^i*t^(m+j-i),
 0<=j<=m, 0<=i<=m+j.                                  (1)

There are(m+1)(3m+2)/2 sections. In the actual cubic frame g0=beta*t^2
on B, put r=s/t, so q is r=0. The symbol of(1) is

 (-1)^(m-j)*r^(2(m-j)+i).                              (2)

This frame and q remain fixed throughout. No arbitrary base change is
substituted for the actual marking. Write W_(m,d) for the elements of
V_m whose ordinary source degree is at most d.

The [exact certificate](contact_filtration_certificates.json) gives a
rational invertible triangular change from(1) to a polynomial echelon
basis, ordered first by total X,Y degree and then by X exponent. The
[checker](verify_contact_filtration.py) reconstructs every echelon row
from its recorded combination, checks the triangular inverse exists,
and checks distinct leading monomials and each full source polynomial.
Therefore the rows of degree<=d span exactly W_(m,d): a linear
combination involving a higher row has a unique largest source monomial
and cannot lose that monomial. This certifies the complete spaces,
including every lower-pole term, rather than selected lifts.

Near the entire omitted prime use

 X=u^-1, Y=u^2*(u*v-1).

All section polynomials become polynomials in u,v; the checker verifies
cancellation of every negative u power. Write F=F0(v)+u*F1(v)+... and
write s_k for the coefficient of r^k in its leading B-symbol at pole m.
For elements with smaller actual pole all these s_k are zero.

The complete filtration computation gives:

| Space | Dimension | Unique top source form | Max v degree | Max u degree on that v face |
| --- | --- | --- | --- | --- |
| W_(6,72) | 63 | X^54 Y^18 | 18 | 1 |
| W_(7,72) | 77 | X^54 Y^18 | 18 | 2 |
| W_(8,108) | 112 | X^81 Y^27 | 27 | 1 |

"Unique" means the homogeneous image at that exact degree is one
dimensional. Its coefficient equals the u^0 coefficient on the top
v face. In W_(7,72) every leading B-symbol is divisible by r^8.
All entries follow from the full rational echelon basis and exact
substitution; the receipt records this finite computation's scope.
No all-height claim is inferred from these bounded spaces.

## The highest v faces are constants in u

Let P,Q be an actual surviving pole(6,8) or(7,8) pair. PF2RULED11-14
leave precisely ordinary degrees72,108 in this order. Thus P and Q
belong to the corresponding spaces above, with top v coefficients

 P=A(u)*v^18+..., Q=B(u)*v^27+...,
 deg(A)<=1 or2, deg(B)<=1, A(0)*B(0)!=0.

The actual Jacobian identity is J_(u,v)(P,Q)=kappa*u. Its v^44
coefficient is

 27*A'*B-18*A*B'=0.

Consequently A^3/B^2 is constant. Factoring over C shows
A=a*T^2, B=b*T^3 for a polynomial T and nonzero a,b. Since deg(B)<=1,
T is constant. Hence the whole top v faces are nonzero constants.
In particular

 deg(P0)=18, deg(Q0)=27,
 deg(P1)<=17, deg(Q1)<=26.                            (3)

This deduction uses the full Jacobian and the complete support bounds,
not merely the ordinary leading terms.

## Exact first-jet identities and the common cusp multiplier

For every F in W_(6,72), the complete basis calculation gives

 [v^18]F1=-s5,
 [v^17]F1-13*[v^18]F0=-s6,
 [v^16]F1-12*[v^17]F0=5*s5-s7.                       (4)

For every F in W_(8,108), it gives

 [v^27]F1=-s4,
 [v^26]F1-20*[v^27]F0=-s5,
 [v^25]F1-19*[v^26]F0=7*s4-s6,
 [v^24]F1-18*[v^25]F0+7*[v^27]F0=7*s5-s7.            (5)

These are exact linear identities on entire finite section spaces.
They are verified on every echelon row after substitution, retaining
all symbol coefficients and all lower-pole lifts.

The [actual quartic cusp theorem PF2OMIT9](../jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md)
supplies h=v-v*, the monic derivative gcd, and polynomial lambda with
lambda(v*)!=0 such that

 h*P1=lambda*P0', h*Q1=lambda*Q0'.                    (6)

The Q symbol has contact order at least8 for(6,8) and at least8 for
(7,8). Hence s4,...,s7 are zero in(5). Put

 c=[v^27]Q0!=0, d=[v^26]Q0, e=[v^25]Q0.

By(3) and(6), deg(lambda)<=1. Equation(5) then forces

 lambda=(20/27)*(v-v*)+delta,
 delta=lambda(v*)=-7*d/(729*c)!=0.                    (7)

In particular d!=0, and lambda has degree exactly one. The next equation
of(5) also determines the fixed actual cusp parameter:

 v*=(2*e+27*c)/d-26*d/(27*c).                         (8)

Equations(7)-(8) are necessary restrictions on a proposed pair; no global
Q or numerical cusp point is supplied here. The fixed parameter v=D
cannot be translated freely while preserving the actual model.

## Excluding the four contact alternatives

For(6,8), PF2RULED11-14 allow the common sextic sigma to have contact
multiplicity k in{2,3,4,6}, and the P symbol is a nonzero multiple of
sigma^3. Thus s5=0 and its first nonzero order is3k.

Let p=[v^18]P0!=0. Equations(4),(6),(7) give

 13*p-s6=[v^17]P1=18*(20/27)*p,
 s6=-p/3!=0.                                        (9)

Therefore3k=6 and k=2. This excludes k=3,4,6; it leaves exactly
sigma divisor4p0+2q, p0!=q. The normalized coefficient identity(9)
is retained as an additional constraint on the surviving stratum.
The away point may be the marked point r=infinity; it was not removed
by choosing an affine root coordinate.

For(7,8), its P symbol is a nonzero multiple of tau^7. Since every
W_(7,72) symbol vanishes to order>=8, k=ord_q(tau)=1 would give the
impossible order7. Thus k=2 or3, leaving tau divisor p0+2q or3q.
This step needs no new Newton-corner import.

## A precise limit of first-jet searches

The lower-pole space W_(6,72) maps surjectively to

 C[v]_(<=18) + u*C[v]_(<=18) modulo u^2.

The exact matrix rank is38. Thus a W_(7,72) element can have its entire
first jet adjusted arbitrarily within these degree bounds by an element
of V6 without changing its pole-seven symbol or source bound. This is
a negative control against attempting to exclude the remaining(7,8)
symbols from the unrestricted first-jet identity alone.

For comparison, W_(5,72) has first-jet rank35 of38, and W_(7,108)
has rank52 of56; the missing top equations are exactly the source of
(4)-(5). The rank computation proves only linear first-jet freedom,
not freedom for a pair satisfying the full Jacobian, finiteness,
second jets or the global image conditions.

## Assurance and next gates

The [receipt](contact_filtration_receipt.json) covers the complete exact
section-space calculation and rational basis changes. The global
identification, finite-map hypotheses and cusp/Euler input are the
written theorems explicitly cited above. No formal verification,
independent full reproof, external review or novelty priority is claimed.

The actual ordinary leading monomials have now been proved, but no
identification with a published standard Newton-corner case is inferred
without checking the full normalization and coordinate-change hypotheses.
The conditional2026 degree125 upgrade is not used.

The next gates are the three remaining contact strata, with the fixed
multiplier(7), coefficient restrictions(8)-(9), higher cusp jets, and
actual global Jacobian/finiteness. The all-height pole-four problem,
higher poles, arbitrary actual-normalization compatibility, both full
quartic packets and JC2 remain open. The full goal remains active.

The active frontier is now height9, with first reduced pairs
(5,9),(6,9),(7,9),(8,9). Do not repeat a higher-jet realization search on
the three height-eight strata; PF2RULED16 closes them by a full chart-leading
Jacobian obstruction. Arbitrary compatible completion remains open.
