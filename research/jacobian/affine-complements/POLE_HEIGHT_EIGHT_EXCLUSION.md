# A primitive chart-leading form excludes the entire compatible quartic height-eight layer

14 November 2026. PF2RULED16 excludes every remaining contact(6,8) and
(7,8) presentation on the actual marked Wright model. Combined with the
previous exclusions and polynomial target-pole reduction, every compatible
finite quartic Keller map has pole height at least9 in every target
presentation. The first unexcluded reduced pairs are

 (5,9), (6,9), (7,9), (8,9).

The three contact strata of PF2RULED15 are now excluded as realization
targets. Its exact filtration and first-jet calculations remain valid
historical controls and inputs. No height-nine, all-height pole-four,
arbitrary compatible-completion or full quartic-packet theorem follows.
The full research goal remains active.

## A finite-dimensional algebraic obstruction

Let W_(m,d) consist of functions in V_m with ordinary source degree<=d,
as in the [complete double filtration](CONTACT_SOURCE_FILTRATION.md).
All functions are in the specified inclusion

 A=C[Y,XY,X^2Y,X^3Y+X] inside C[X,Y],
 X=u^-1, Y=u^2*(u*v-1).

Let P belong to W_(7,72), and let Q belong to W_(8,108) with exact
ordinary source degree108. In the fixed cubic frame g0=beta*t^2 on B,
write s_k for the r^k coefficient of Q's leading pole-eight symbol;
all s_k are zero if its actual pole is smaller than eight. Suppose

 s4=s5=0.                                                (1)

Then there is no nonzero kappa for which

 J_(u,v)(P,Q)=kappa*u.                                   (2)

This core statement does not assume finiteness, field degree four or the
cusp/Euler identity. Those hypotheses enter the application by forcing
the precise pole/source/contact conditions above. The actual inclusion
and marked chart are essential.

## The complete chart-leading identity

Reconstruct the PF2RULED15 rational echelon bases from their complete
Cox section basis. Their source leading monomials and invertible
triangular changes certify exactly W_(7,72), dimension77, and
W_(8,108), dimension112, including every lower-pole term.

The full W_(7,72) basis has chart total degree at most20:

 deg_(u,v)(P)<=20.                                      (3)

No assumption on its leading contact symbol is needed for this bound.

For every Q in the whole W_(8,108) space, put c=[v^27]Q0. The exact
homogeneous projection to chart total degrees at least27 satisfies

 Q_[>=27] = c*v^7*(u+v)^20
             -s4*u*v^7*(u+v)^20
             -s5*u*v^7*(u+v)^19.                        (4)

There are no terms of degree greater than28. The s4 term has degree28,
and the other two terms have degree27. Equation(4) is verified on all112
echelon rows, with exact rational coefficients, after reconstructing the
full source polynomials and changing to the actual chart. Linearity
therefore proves it on the entire complex section space. No restricted
lift, omitted remainder or sampled root position is involved.

Moreover, the unique ordinary homogeneous form at source degree108 is
X^81Y^27, and

 c=[X^81Y^27]Q=[v^27]Q0.

Since the source degree is exactly108, c!=0. Under(1), equation(4) proves

 deg_(u,v)(Q)=27,
 Q_27=c*v^7*(u+v)^20.                                  (5)

The [certificate](height_eight_certificates.json) retains(4), its expanded
polynomials, the reconstructed-space dimensions and the nonresonance
calculation below. The [checker](verify_height_eight.py) verifies every
basis row and every high-degree projection, rather than relying only on
support hulls or on the one-dimensional top source image.

## The leading Jacobian cannot vanish

If(2) held, P would be nonconstant. Let d=deg_(u,v)(P), so1<=d<=20,
and let P_d be its nonzero homogeneous part. The potential highest
Jacobian degree is d+27-2=d+25>1. As the prescribed Jacobian has degree
one, its highest component must vanish:

 J(P_d,v^7*(u+v)^20)=0.                                 (6)

Use the invertible linear chart w=u+v, keeping v. Write
P_d(u,v)=F_d(w,v). Equation(6) becomes

 v^6*w^19*(7*w*partial_w(F_d)-20*v*partial_v(F_d))=0.

For a degree-d monomial w^i*v^(d-i), the remaining diagonal operator
acts by

 7*w*partial_w-20*v*partial_v:
 w^i*v^(d-i) -> (27*i-20*d)*w^i*v^(d-i).                (7)

If27*i=20*d, coprimality of27 and20 implies27 divides d. This is
impossible for1<=d<=20. Every coefficient of F_d must therefore be
zero, contradicting its definition. This proves the algebraic obstruction.

The bound on d matters: degree27 itself has the nonzero resonant control
v^7*w^20. The checker retains that commuting control as well as every
nonzero eigenvalue in the complete degree range1..20. No assertion that
the centralizer contains only constants in unrestricted degree is made.

## Application to all compatible quartic height-eight presentations

The full compatible-completion and period reductions identify the actual
source open, volume and marked ring used above. Thus a source Keller
identity is precisely J_(u,v)=kappa*u, kappa!=0, in this chart.

[PF2RULED15](CONTACT_SOURCE_FILTRATION.md), using the retained published
2022 degree alternative and the preceding intrinsic gates, leaves:

* Pole(6,8): source degrees72,108 and a Q symbol of contact order8.
* Pole(7,8): source degrees72,108 and a Q symbol of contact order16 or24.

In particular P belongs to W_(7,72), since V6 is contained in V7, and
Q satisfies the exact source degree108 and(1). The obstruction excludes
all three remaining contact strata, including every away-point position
and every lower-pole lift. It needs neither higher cusp-jet elimination
nor construction of a companion or residue calculation.

All smaller pole heights are already excluded by PF2RULED7. At reduced
height8, PF2RULED8 removes coordinate poles<=3; divisible-pole shears
handle(4,8) and equal-pole cases; PF2RULED12 closes(5,8). The new result
closes(6,8),(7,8). These are every reduced height-eight possibility.
The shears preserve the target algebra, actual map degree and source
Jacobian, and do not increase the maximum pole. Therefore every target
presentation has height>=9, not merely a selected reduced one.

At height9, all-height coordinate exclusions and divisible shears leave
minimum coordinate pole>=4. PF2RULED9 excludes the odd companion(4,9),
so the first unexcluded reduced pairs are exactly(5,9),(6,9),(7,9),(8,9).
The gcd3 pair(6,9) is outside the rational gcd1/gcd2 period criterion;
no contact-symbol classification is asserted for it.

## Assurance and remaining obligations

The [receipt](height_eight_receipt.json) seals the proof, exact certificates
and the PF2RULED15 inputs. The checker reconstructs the complete rational
bases and verifies the universal identity(4), source coefficient equality,
chart degree bounds and every eigenvalue in the exhaustive degree range.
It shares exact polynomial arithmetic utilities with PF2RULED15 and is
not an independent whole-proof implementation. The section-space and
compatible-source identifications remain the cited written geometry.

The core obstruction is intrinsic finite algebra. The quartic application
uses the already pinned2022 degree alternative, with its exception
retained globally: this proves its impossibility in these pole spaces,
not at arbitrary pole height. No published standard Newton-corner
identification, conditional2026 upgrade, formal verification, independent
whole-proof replay, external review or novelty priority is claimed.

Height9 and higher, all-height pole-four contact coordinates, the
unrestricted mixed-weight marked-ring problem, actual compatible embedding
of arbitrary quartic normalization packets and both full quartic packets
remain open. An abstract affine-plane complement still does not supply
that finite map or the required compatible completion.
