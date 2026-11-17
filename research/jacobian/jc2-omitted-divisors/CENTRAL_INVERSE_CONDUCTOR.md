# Draft: central inverse conductor accounting

14 November 2026. **Unregistered, unverified draft.** The current replay
fails its independent L0 divided-difference resultant check. The formulas
below require correction and review; no PF2OMIT20 claim has been added
to MATH_STATUS.json. Work stopped when the user cleared the goal.

This proposed continuation of PF2OMIT19 concerns an admissible inverse,
conductor polynomials and target representatives. It supplies no actual
inverse, additional image, retained-residue choice or global Keller source.

## 1. An explicit degree24 conductor on the known first image

Use the monic birational normalization (p(t),q(t)) of C and its reduced
equation h=H_C, monic of V-degree6. Then

    c_C(t)=H_C,V(p(t),q(t))/p'(t)
          =Res_s((p(s)-p(t))/(s-t),(q(s)-q(t))/(s-t)),
    H_C,U(p,q)=-q'(t)*c_C(t).                       (I1)

For the resultant identity, factor p(s)-p(t)=(s-t)*Delta_p in the
defining resultant for H_C(p(t),V), differentiate in V at V=q(t),
and use multiplicativity. The two signs from Delta_p's degree5 cancel.

This polynomial generates the conductor ideal of C[p,q] in C[t], up
to a nonzero scalar. On an immersed branch of a multibranch value,
its order is the sum of its contacts with the other branches. At the
ordinary cusp its order is2, as in C[t^2,t^3]. At a smooth projection
critical point it is a unit: the chain rule cancels p', while the
other gradient/parameter derivative is a unit. The complete catalogue
therefore gives the conductor orders

| Image value | Orders at its normalization parameters |
|---|---|
| node | 1,1 |
| ordinary triple | 2,2,2 |
| tacnode | 2,2 |
| A5 | 3,3 |
| ordinary cusp | 2 |

Their sum is24 at every coefficient point. Equivalently the infinity
gradient has pole29, so dividing by p' of degree5 gives exact degree24.

The [models](central_inverse_ledger.json) give every coefficient of
c_C on L and scaled K. Their leading coefficients are-1/6561 and
2z^6/19683 respectively, units on the full domains d in C and z!=0.
No cusp, collision or auxiliary parameter locus is divided out.

For efficient construction, import the exact key equation of
[PF2R7KEY1](../jc2-r7-elliptic-contact/PROOF.md):

    H_C=G_C^2+ell(U)*G_C+H(U)+L(U)*V.

With g=G_C(p,q), (I1)'s numerator is
(2g+ell(p))*G_C,V(p,q)+L(p). The [replay](verify_central_inverse.py)
checks both gradient identities as full rational polynomial identities
on each complete chart. Separate divided-difference resultants at
exact witnesses check the construction independently.

The full C[p] basis from PF2OMIT16 also produces a unique filtered
target representative A_C(U,V), with

    A_C(p,q)=c_C(t), total degree(A_C)<=7.           (I2)

The thirteen available pole degrees at most24 come from p^i*q^j*g^k,
j0..2,k0..1. Their degree bounds give
A0(U)+A1(U)V+A2*V^2+B0(U)G_C+B1*V*G_C, with degrees4,2,0,2,0
in the five coefficients. The resulting two explicit representatives
are recorded with the conductor polynomials. This is a polynomial
identity in the original target coordinates; the finite key map is
not being treated as a Keller-preserving coordinate change.

## 2. Truncate the inverse before normalizing the central curve

Choose a generic admissible source coordinate T=x+c*y as in PF2OMIT3,
and write its primitive degree6 inverse as R6=sum a_i(U,V)*T^i.
In addition choose T to separate the generic retained points over C.
This excludes only finitely many additional slopes: for any pair of
distinct source points over the generic C point, their x,y coordinates
do not both agree, so at most one constant slope makes T agree.

The surviving central polynomial is

    S_C(T;t)=R6(T;p(t),q(t)), deg_T S_C=m,
    m=3 for A,H=C, and4 otherwise.                  (I3)

All a_i(p,q) for i>m vanish identically. Its coefficient ideal in
C[t] is the unit ideal, by pulling back the original coefficient-unit
identity. Its degree-m projective root cover is finite flat; its
normalization is exactly the retained finite curve E_C of PF2OMIT19.
The whole omitted R component and any H=C component have been removed
by this truncation. Their original ranks must not be reinserted into m.

Let p_R,p_S,p_D be the positive pole orders of T along the original
omitted primes; p_D exists only in A. At every finite point of such
a prime, 1/T is its local equation to that power times a unit. This
follows from the global infinity Cartier divisor and disjointness of
the original omitted boundary. Thus restriction to a curve multiplies
that pole order by the curve's local contact with the boundary.

For a primitive polynomial of degree m over the DVR of t-t0, whose
coefficient ideal is a unit, the valuation of its leading coefficient
is the sum of all normalized root pole orders. This follows from
multiplicativity of the Gauss valuation, or the reciprocal-root Newton
polygon. A normalization-index factor cannot add zeros to this charge.

## 3. The exact conductor-weighted leading coefficient

Apply that pole formula to every residual deleted point of E_C.
For a cyclic block w^e=t^r, its gcd(e,r) branches each have T pole
p*r/gcd(e,r), so their total pole charge is p*r. This remains true
when every normalized branch is unramified. At the own cusp of C,
the residual branch (u,w)=(-3t^2,2t) meets the boundary u+3w^2=0
to order2, contributing2p_R, precisely the cusp conductor order.
At B's cusp, a transverse C path has pole2p_S on its index3 branch;
a tangent path has poles p_S,2p_S on its index1,2 branches. Their
totals are again p_S times the respective image contacts2,3.

The other R branches therefore contribute c_C^p_R. If H=C, its
other D branches contribute another c_C^p_D; C is immersive in this
grouping. The other images contribute their reduced equations pulled
back to C. Consequently, with gamma a nonzero complex scalar,

| Grouping | a_m(p(t),q(t)) / gamma |
|---|---|
| A, H distinct | c_C^p_R * H_B(p,q)^p_S * H_H(p,q)^p_D |
| A, H=C | c_C^(p_R+p_D) * H_B(p,q)^p_S |
| A, H=B | c_C^p_R * H_B(p,q)^(p_S+p_D) |
| B,N,C | c_C^p_R * H_B(p,q)^p_S |

These are exact polynomial identities, not just lower bounds. The
only finite zeros are the listed pole points, so equality of divisors
gives a constant ratio. Replacing c_C by A_C turns each row into the
same exact congruence for a_m(U,V) modulo H_C. In particular this is
a concrete coefficient constraint for an actual global inverse.

Write d=deg H_B(p,q), and u=deg H_H(p,q) only when H is distinct.
The exact t-degree of this coefficient is24p_R+d*p_S, plus

    24p_D if H=C; d*p_D if H=B; u*p_D if H distinct. (I4)

Under the separate F2 source-degree(75,125) hypothesis, PF2OMIT14's
pole budget W<=24 and d<=30j, u<=30l give degree at most684, or648
when H=C. Indeed these bounds are30W-36p_R, or
30W-36(p_R+p_D). They are necessary caps, not geometric sharpness or
new realizations of the pole labels.

## 4. Central discriminants still separate ramification and index

Add the clean hypotheses of PF2OMIT19: C has twelve nodes, B avoids
them, and all C/B crossings are transverse smooth crossings. In the
shared-image A rows and in row C, c_C and R_CB=H_B(p,q) are coprime
and squarefree. The normalized finite-curve discriminant of E_C is

    c_C * R_CB^(e_S-1), up to a nonzero scalar.

At a node parameter the residual fold contributes different1; at a
crossing S contributes e_S-1. Extra unramified omitted points contribute
zero different. The trace-index identity of PF2OMIT4, now over C[t],
therefore gives

    Disc_T(S_C)=constant*c_C*R_CB^(e_S-1)*J_C(t)^2.  (I5)

Here J_C is the finite binary-order normalization-index determinant.
The argument applies to disconnected normal curves and product fields
as well. Extra squared zeros can be coordinate collisions. This is
the discriminant of the **truncated** central equation; restricting
the degree6 discriminant to C merely gives zero and does not yield(I5).

The square classes are c_C*R_CB in A and c_C in row C. Thus the
cyclic-index3 image is still omitted and ramified even when absent
from the central sign square class. The global sign equation sqrt(H_C)
cannot simply be restricted to C to recover this central sign cover.

In A,H=C, at each node parameter all three retained-cover sheets are
deleted: a residual fold and an unramified D point. Their reciprocal
T-root valuations are p_R/2,p_R/2,p_D. If these sorted values are
q1<=q2<=q3, then a0 is a unit and

    ord a_i(p,q)>=ceil(q1+...+qi), i1,2,3,
    ord a3(p,q)=p_R+p_D,
    ord Disc S_C>=p_R+4min(p_R/2,p_D),
    ord J_C>=ceil((p_R+4min(p_R/2,p_D)-1)/2).        (I6)

The exact different there is1. In particular unit poles p_R=p_D1
force raw discriminant order at least3 and index order at least1;
the extra order is not an additional ramified prime. Exact local
quadratic-plus-linear controls attain these integer bounds, but are
not global Keller realizations.

## 5. The split residue alternatives become inverse factors

In clean row C, a split3+1 retained cover requires the polynomial
lift of PF2OMIT19. Hence the degree4 S_C has a monic linear factor
T-A(t), A in C[t]. The remaining cubic has the same leading
coefficient c_C^p_R*R_CB^p_S, up to a constant.

In clean A,H=C, the split2+1 section has poles exactly of order p_D
at every one of the24 node parameters and no other finite pole.
Writing c_C up to a scalar as the monic node-preimage polynomial,
its T-value is A(t)/c_C(t)^p_D with gcd(A,c_C)=1. Gauss's lemma gives

    S_C=(c_C^p_D*T-A(t))*Q2(T;t),
    leading_T Q2=constant*c_C^p_R*R_CB^p_S.          (I7)

The quadratic field is C(t)(sqrt(c_C*R_CB)), as before. Its raw
discriminant has that square class; the full cubic discriminant gains
the squared resultant with the linear factor. This resultant is
coordinate-order data, not an extra omitted or ramified component.
The section's coordinates can have poles; replacing this factor by
a monic polynomial section would erase all24 unramified omissions.

The connected cubic case has the same quadratic sign field only in
its Galois closure, not as an intermediate field of its degree3 curve
extension. Likewise row C's connected quartic has its conductor sign
field in the S4 Galois closure. The central residue extensions are
different objects from the primitive degree6 surface or generic-pencil
extension; no forbidden global field tower has been asserted.

Actual a_i,B,H,pole orders,index factors,section numerators and the
global source still need a compatible normalization/source package.
The [replay](verify_central_inverse.py) checks exact complete conductor
models, independent resultant/cusp controls, all pole-charge formulas
and local index controls. The inverse/curve-normalization and trace
arguments are written proofs, not a formalization or external review.
