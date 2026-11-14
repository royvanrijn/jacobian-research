# Global cusp jets exclude every one-root noncontact leading symbol

14 November 2026. PF2RULED14 proves an all-height obstruction on the
actual marked Wright model of PF2WRIGHT1, conditional on a finite
quartic Keller target algebra. No element of that target algebra can
have B-pole m>=1 and leading symbol tau^m when the cubic tau has
divisor3p with p distinct from the contact point q.

This statement includes elements that are not target coordinates.
Their first boundary jets inherit the same target-algebra identity;
we do not infer a new constant-Jacobian pair from a noncoordinate.

Combined with PF2RULED11, every coordinate pair whose pole gcd is one
or two must have its common leading section vanish at q. In particular,
the entire noncontact(6,8) case is excluded. The three PF2RULED13 profiles
are superseded as realization targets, while their exact local curve
controls are retained. The remaining reduced height-eight pairs(6,8)
and(7,8) both require contact symbols and ordinary source degrees72,108.
Neither contact pair, the higher-pole problem nor actual ruled
compatibility is settled.

## The fixed marked model and its boundary-jet identity

Use the actual inclusion

 A=C[Y,XY,X^2Y,X^3Y+X] inside C[X,Y].

Write A1=XY, B1=X^2Y, D=X^3Y+X and C=1+B1+D^2, as in
[PF2WRIGHT1](WRIGHT_MARKED_MODEL.md). Near the entire omitted prime R,

 Y=u^3*v-u^2, A1=u^2*v-u, B1=u*v-1, D=v.

For F in A put F=F0(v)+u*F1(v)+O(u^2). The
[quartic boundary theorem PF2OMIT9](../jc2-omitted-divisors/WRIGHT_IMAGE_RESTRICTION.md)
gives, for an actual finite pair(P,Q),

 h=gcd(P0',Q0')=v-v*,
 P1=lambda*P0'/h, Q1=lambda*Q0'/h,
 lambda(v*)!=0.                                        (1)

Here h is monic, lambda is polynomial, and v is the fixed actual
parameter D, not a freely renormalized branch coordinate. The identity
uses finiteness, quartic degree, the actual source open and the Euler
boundary calculation. Those hypotheses are essential.

By the chain rule, every F in the target algebra C[P,Q] satisfies

 h*F1=lambda*F0'.                                       (2)

No hypothesis that F is a coordinate is needed for(2).

## The symbol-to-jet comparison on the whole pole filtration

Normalize the V1 generators as

 g0=f0, g1=f1, g2=-f2, g3=-f3,
 f0=-Y-C^2-2A1*D-B1*D^2,
 f1=-A1-D-2B1*D-D^3, f2=C, f3=D.

Let r=s/t be the actual base coordinate on B. The contact point is
r=0, and B meets the negative section at r=infinity. Relative to the
g0 cubic frame, the symbols of g0,g1,g2,g3 are1,r,r^2,r^3. This is
the marked frame from the explicit source-coordinate theorem; no
arbitrary automorphism of B is substituted for the actual marking.

Define the corrected first jet

 J(F)=F1-(3/4)*F0'.

It is a derivation to C[v] relative to the restriction F->F0. Direct
substitution in the actual ring gives:

| generator | F0 | J(F) |
| --- | --- | --- |
| g0 | -v^4+v^2 | v/2 |
| g1 | -v^3+v | (v^2+1)/4 |
| g2 | -v^2 | v/2 |
| g3 | -v | 3/4 |

Since Vm is spanned by products of at most m V1 generators,

 deg(F0)<=4m,  deg(J(F))<=4m-2 for m>=1.                (3)

Indeed a product only of g0 has its corrected-jet degree drop by
three. A product with another generator already has restriction
degree at most4m-1, and every corrected-jet term drops degree by at
least one. Fewer factors give smaller bounds; constants have zero
corrected jet. This proves(3) for every height, not by extrapolation
from finitely many monomials.

If the leading symbol is nonzero at q, its g0^m coefficient is
nonzero, while all other binary orders have smaller restriction
degree. The lower-pole remainder has restriction degree at most4m-4.
Consequently deg(F0)=4m exactly.

Assume such a noncontact F belongs to the target algebra. Equations
(2)-(3) force lambda to be linear with leading coefficient3/4:
F1=(3/4)F0'+J(F) has degree4m-1, and F0'/h has degree4m-2.
Thus

 lambda=(3/4)h+delta, delta=lambda(v*)!=0,
 h*J(F)=delta*F0'.                                     (4)

The same h,lambda apply to the whole target algebra.

## A degree gap passes the identity to a V1 lift

Suppose F has leading symbol tau^m, div(tau)=3p, p!=q.
Scale to normalize tau(q)=1 and lift tau to T in V1. Then

 F=T^m+R, R in V_(m-1).

For any G define L(G)=h*J(G)-delta*G0'. It is also a derivation
relative to restriction. Let C(v)=T0(v), of exact degree four. Then

 L(T^m)=m*C^(m-1)*L(T).

For m>=2, (3) bounds deg(L(R)) by4m-5. If L(T) were nonzero, its
product with m*C^(m-1) would have degree at least4m-4 and could not
cancel L(R). Since L(F)=0, we obtain L(T)=0. For m=1, R is constant
and the same conclusion is immediate.

This is only a necessary first-jet identity for the auxiliary T.
It does not put T in the target algebra or an intermediate plane.
The strict degree gap is what makes this argument valid at unbounded m.

## The fixed V1 identity has no solution

If p is the marked point r=infinity, choose T=g0. Then C=-v^4+v^2
and J(T)=v/2. The identity h*J(T)=delta*C' is impossible: the left
has degree two and the right has degree three with delta!=0.

Otherwise write r(p)=a, where a!=0. The normalized cubic and its lift are

 tau(r)=(1-r/a)^3,
 T=g0-(3/a)*g1+(3/a^2)*g2-(1/a^3)*g3.

The exact restriction C and corrected jet L0=J(T) are

 C=-v^4+(3/a)*v^3+(1-3/a^2)*v^2+(-3/a+1/a^3)*v,
 L0=-(3/(4a))*v^2+(1/2+3/(2a^2))*v-3/(4a)-3/(4a^3).

Set h=v-t, t=v*. Comparing the v^3 and v^2 coefficients of

 (v-t)*L0-delta*C'=0

gives

 delta=3/(16a), t=1/(4a)-2a/3.

After this substitution the two remaining coefficients are exactly

 [v] = (4a^2-3)/(12a),
 [1] = -(2a^2-1)/(4a^2).

Their denominators are nonzero. The two numerators cannot vanish,
because the exact polynomial certificate is

 2*(2a^2-1)-(4a^2-3)=1.                                (5)

This excludes every finite noncontact point a as well as the marked
point treated above. The target-algebra obstruction is proved.

## Consequences for the surviving programme

Under compatible completion the period normalization supplies the
actual marked model, so the theorem applies to every such quartic map.
For a pole gcd d=1 or2, PF2RULED11 leaves at most one common-symbol
root away from q. If that symbol did not vanish at q, all its3d
zeros would be at one point p!=q, and it would be a d-th power of a
cubic tau with divisor3p. The component symbols would be tau^m,tau^n,
contradicting the theorem. Hence the common symbol must vanish at q.

At reduced height8, only(6,8),(7,8) remain after PF2RULED12, and their
noncontact branches are now both excluded. The necessary contact lists
from PF2RULED11 remain:

* (6,8): common sextic has contact multiplicity2,3,4 or6.
* (7,8): common cubic has contact multiplicity1,2 or3.

Both require ordinary source degrees(72,108), in the displayed order,
by the pinned2022 alternative already imported in PF2RULED7. This
exception is retained; no source Newton-corner identification or
conditional2026 removal is assumed.

More generally, for a gcd-one/two presentation of height H>=8, the
highest coordinate has contact vanishing at least H/d>=4. The source
filtration gives maximum ordinary degree at most16(H-1). A surviving
pole-four coordinate, whose symbol is tau^4 by PF2RULED10, must likewise
have tau(q)=0; its ordinary source degree is at most48. No theorem
excluding an arbitrary Keller component of degree<=48 is inferred.

The three local curves of PF2RULED13 remain reproducible controls:
they meet local pole/discriminant counts but cannot extend to a finite
quartic Keller target algebra in this actual marked model. Their former
residue/global-Q search is superseded by the obstruction above.

## Assurance and remaining gates

The [checker](verify_cusp_leading.py) replays the actual generator jets,
the complete generator degree-drop bounds, the universal derivation and
power identities, the degree-gap inequality, and the exact unit
certificate(5). It retains [symbolic certificates](cusp_leading_certificates.json)
and a [receipt](cusp_leading_receipt.json). The global all-height
filtration argument and the PF2OMIT9 finite/Euler input are written
proofs. No formal verification, independent whole-proof implementation,
external review or novelty priority is claimed.

The contact72,108 height-eight cases, higher pole gcds and contact
families, the all-height pole-four problem, actual ruled compatibility,
both full quartic packets and JC2 remain open. The full goal remains active.

The subsequent [PF2RULED15 complete contact filtration](CONTACT_SOURCE_FILTRATION.md)
further excludes contact multiplicities3,4,6 for(6,8) and1 for(7,8).
Only sextic4p+2q or cubic p+2q/3q remain, with boundary restriction
degrees18,27 and common cusp multiplier slope20/27. The larger contact
lists above are this proof's earlier necessary table.

The subsequent [PF2RULED16 full height-eight exclusion](POLE_HEIGHT_EIGHT_EXCLUSION.md)
closes both remaining pairs(6,8),(7,8), including every contact stratum.
The earlier height-eight handoffs in this note are historical. The current
minimum height is9; first reduced pairs are(5,9),(6,9),(7,9),(8,9).
