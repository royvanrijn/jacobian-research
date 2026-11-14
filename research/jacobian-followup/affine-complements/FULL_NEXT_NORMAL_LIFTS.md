# All sixteen marked candidates lift through the full next-normal equation

14 November 2026. PF2RULED24 solves the remaining next-normal lifting
problem of [PF2RULED23](CUSP_HIGH_JACOBIAN_REDUCTION.md): each of the
sixteen marked noncontact pole(6,9) candidates admits actual global
sections with

    Jac_uv(P,Q)=u+O(u^2).

Their boundary maps have a clean singleton ordinary cusp and are
birational onto their images on a nonempty open locus. The eleven
previous lower consistency conditions therefore have simultaneous
solutions. This test removes no marked candidate. These are finite
polynomial sections, not solutions of the full Keller identity or
finite degree-four surface maps. The coefficient of u^2 and all later
coefficients, contact cases and actual compatible embedding remain open.

## Fix P and solve for the entire second coordinate

Keep the actual Wright ring, omitted-prime chart and complete first
cusp-lift spaces of [PF2RULED22](FULL_CUSP_FIRST_JET_LIFTS.md). Write

    F=F0(v)+u*F1(v)+u^2*F2(v)+u^3*F3(v)+...,
    h=v-s, lambda=3h/4+delta, L(F)=h*F1-lambda*F0'.

The marked constants delta and s are fixed by the common nonic H.
All pairs with symbols H^2,H^3 and L(P)=L(Q)=0 form an affine space
of dimension31+85=116. For any such pair the constant coefficient of
the chart Jacobian vanishes, and its u coefficient is

    J1(P,Q)=2*(P2*Q0'-P0'*Q2)+P1*Q1'-P1'*Q1.          (1)

PF2RULED23 proves that its degree in v is at most48 on the complete
family. Fix one first coordinate P. Equation (1) is linear in Q, so
the full next-normal problem is one affine linear system of49
coefficient equations in the85 parameters of Q. This includes all
three high conditions and all eleven lower membership conditions;
no normal controls are frozen before solving it.

For each of the two degree-eight marked candidate fields K, the exact
certificate gives a P for which

    ker(L|V8) -> K[v]_(<=48), U |-> J1(P,U)

has rank49. Hence every right-hand side is available, including
J1(P,Q)=1. The reference Q has the prescribed symbol, and all its
corrections lie in V8, so the resulting pair retains H^2,H^3.

Surjectivity persists on a nonempty open subset of the affine31-space
of P. Over it the fiber with J1=1 has dimension85-49=36. If the
nonzero constant kappa=J1 is allowed to vary, the fiber has dimension37
and the pair chart has dimension31+37=68. Its projection to the data
through u^2 has dimension32: the complete directions invisible through
u^2 have dimensions3 and33 and act freely. These are dimensions of
this surjective chart, not a classification of the rank-drop strata.
Those strata remain in scope for any eventual complete exclusion.

## Exact characteristic-zero certificates from residue units

The new proof uses an exact reduction to certify nonzero determinants,
not a bounded search extrapolation. The arithmetic lifting argument
is as follows.

For each irreducible polynomial P8 from PF2RULED19 put

    K=Q[y]/(P8), p=1009, y0=990.

Normalize P8 to be monic. Every denominator of P8 and x(y) is prime
to1009, and P8(990)=0 modulo1009. The evaluation y->990 defines a
residue map on the corresponding local subring O of K. Concretely,
start with Z_(1009)[y]/(P8) and localize at the kernel of that
residue map. Irreducibility makes this a subring of the field K.
Any element with nonzero residue is a unit of O and a nonzero
member of K. Matrix inversion using such determinants stays in O.

All constructions below use only additions, products and the inverses
of determinants whose nonzero residues are recorded. They therefore
define exact finite arithmetic circuits over K, regular at the chosen
residue point. No characteristic-zero solution is inferred from the
mere existence of a finite-field solution of nonlinear equations.

The circuit is fully specified as follows.

1. Use the complete signed Cox bases of PF2RULED23. Choose, in the
   documented basis order, the first section for each symbol power
   r^i, i=0,...,9, and combine them with the coefficients of H to form
   T in V3. This keeps the marked frame and its required signs.
2. For m=5,8 let M_m be the full matrix of L on V_m. The certificate
   specifies4m pivot columns I_m. Their square minor is a residue
   unit. Embed its inverse on those columns to obtain a right inverse.
   Correct T^2 and T^3 to references P*,Q*. For a free column j use

       U_j=e_j-embed_(I_m)(M_(I_m)^-1*M_j).

   These vectors form the complete kernel, because the chosen minor
   is invertible. All31 and85 affine directions are retained.
3. Put P=P*+sum U_j using all31 first-coordinate kernel vectors with
   integer weight1. Form the49-by85 matrix B with entries

       B_(ell,j)=[v^ell]J1(P,U_j), ell=0,...,48.

   The certificate specifies49 pivot columns I. The determinant of
   B_I has nonzero residue. Set all other Q kernel controls to zero
   and choose the49 controls by

       B_I^-1 * coefficients(1-J1(P,Q*)).

   This constructs the exact global section Q over K. It solves the
   whole polynomial J1=1, because the characteristic-zero degree
   bound48 is already proved on the full family.

The new nonzero residue determinants and singleton-cusp resultants are:

| Common-nonic multiplicities | y0 | x(y0) | det(B_I) mod1009 | singleton resultant mod1009 |
| --- | --- | --- | --- | --- |
| (2,7) | 990 | 912 | 796 | 38 |
| (4,5) | 990 | 264 | 285 | 240 |

Both full cusp inverse minors and all division guards are retained in
the [certificate](full_next_normal_certificates.json). The checker
reconstructs every residue matrix from the complete global sections,
verifies its determinant and the full J1=1 residue identity, and checks
irreducibility of both P8 over Q. The written local-ring argument
above supplies the characteristic-zero implication.

A nonzero element of K remains nonzero at every complex embedding.
Thus each circuit supplies a solution at all eight conjugate marked
points in its sector, not just at one selected numerical root. The two
sectors together cover all sixteen points. Exact characteristic-zero
coefficients are specified by these inverse-matrix circuits; their
potentially large expanded rational representations are not needed
or claimed to have been independently replayed.

## Why the clean cusp and singleton conditions survive

Put p=P0'/h, q=Q0'/h. The full cusp identity gives P1=lambda*p and
Q1=lambda*q, so (1) becomes

    2h*(P2*q-p*Q2)+lambda^2*(p*q'-p'*q)=1.              (2)

If p and q shared any finite root, its substitution into (2) would
give0=1. Hence gcd(p,q)=1. At v=s equation (2) gives

    delta^2*(p*q'-p'*q)(s)=1.

Since delta!=0, the second/third derivative wedge of the boundary map
is nonzero. Its only derivative-gcd root is h, and it is an ordinary
cusp. These conclusions hold for every J1=1 solution in the chart,
not just for the particular circuit witnesses.

For singletonness, translate the target so P0(s)=Q0(s)=0 and write

    P0=h^2*A, Q0=h^2*C.

The circuit witnesses have a nonzero resultant Res(A,C), certified by
the two residues in the table. Because these polynomials and their
resultant belong to O, the resultant is a nonzero element of K. Thus
no second boundary parameter maps to the cusp image. Nonvanishing of
this resultant is open and supplies a nonempty clean singleton cusp
locus in the characteristic-zero solution chart. The existing
normalization-degree argument of PF2RULED22 then proves birationality
onto the boundary image. It proves no surface finiteness statement.

## The following normal coefficient and all its hidden controls

The next required coefficient is

    J2=[u^2]Jac_uv(P,Q)
      =3*(P3*Q0'-P0'*Q3)
       +2*(P2*Q1'-P1'*Q2)+P1*Q2'-P2'*Q1.             (3)

A Keller pair must have J2=0. The circuit witnesses have a nonzero J2,
with degree43 already after reduction modulo1009. This proves that
these actual characteristic-zero witnesses still fail the full Keller
identity. It excludes neither their marked candidates nor other
solutions of (2).

For any fixed data through u^2, the complete remaining lower-pole
freedoms are the kernels computed in PF2RULED22: dimension3 in V5 and
33 in V8. The new exact rational full-Cox calculation gives their
entire third-coefficient images:

    P3 <- P3+A,
    A in span{(v^2-1)^2, v*(v^2-1), v^3*(v^2-1)},

    Q3 <- Q3+B, deg B<=17.

The first map has rank3 and no further kernel. The second has rank18
and a15-dimensional kernel invisible through u^3. Every listed space
is a full image of a complete kernel, not a monomial ansatz. In
particular the P3 image is only three-dimensional and every element
is divisible by v^2-1; it is not the whole degree-five polynomial space.

These controls change (3) by

    3h*(A*q-p*B),                                      (4)

of degree<=40. Since gcd(p,q)=1 and deg p=22 while deg A<=5, the map
in (4) is injective: a zero image forces p|A, hence A=0 and B=0.
Its rank is therefore3+18=21 at every clean noncontact candidate.
The other15 directions must be retained at later normal orders.

The crude cross-term bound is52, but the intrinsic source Jacobian
improves it to43. Put

    D=Jac_XY(P,Q)=-Jac_uv(P,Q)/u.

The first cusp identity makes the chart Jacobian divisible by u.
Therefore D is regular across the actual omitted prime; it is already
regular on the source open, so it belongs to the normalization ring.
At B, poles6 and9 would give a volume pole of at most16. The symbols
H^2,H^3 make its coefficient9*(H^2)'*H^3-6*H^2*(H^3)' vanish.
Thus dP wedge dQ has pole at most15. The actual source volume has
pole3 at B, so

    D in V12.

After J1=1, D+1 vanishes on E and belongs to H0(12B-E). Its first
normal coefficient has degree<=12*(B.E)-E^2=48-5=43 on E. The
normal frame u is regular and nonvanishing as a frame along the entire
affine E; its missing frame point is precisely q. Hence this coefficient
is a polynomial in the actual v of degree<=43. Equivalently,

    deg J2<=43.                                       (5)

For an exact full-space control, the checker uses all247 sections of
V12. The restriction map has kernel dimension198; the entire first-jet
image of that kernel is the44-dimensional space of polynomials of
degree<=43. This controls the actual chart frame without transporting
an abstract line-bundle bound through an unspecified trivialization.
The divisor argument placing D in V12 is a written proof.

Thus a complete next necessary test for fixed data through u^2 is:

- coefficients in degrees41,42,43 vanish (at most3 equations);
- J2(s)=0, so the remaining degree<=40 polynomial is divisible by h;
- its quotient by h belongs to the21-dimensional image
  3*(A*q-p*B) in the40-dimensional space of polynomials of degree<=39.
  This gives19 local matrix-membership conditions.

This is a bound of23 consistency conditions, not a claim that they
are independent or that the high3 coefficients all survive the earlier
equations. Their dependency, solvability and rank-drop strata remain
open. No dimension count certifies a finite map or an exclusion.

## A finite gate for the entire source Jacobian identity

The same V12 bound shows exactly how many further normal coefficients
are sufficient. Write B=C0+2f and E=C0+3f on the marked F1. Then

    12B-9E=3C0-3f

has no global sections: the effective cone is generated by C0 and f.
Equivalently the ruling decomposition has only line bundles of
negative degree. Consequently a section of V12 vanishing to order9
on E is zero. The complete rational matrix of its first nine jets
has rank247, which the checker verifies directly as a second finite
control on this statement.

It follows that, in this noncontact pole(6,9) family, the equations

    J1=kappa!=0, J2=J3=...=J9=0

are sufficient for the full source Keller identity. Indeed
Jac_uv=kappa*u+O(u^10) makes D+kappa vanish to order9 on E, hence
identically. Conversely a constant source Jacobian satisfies them.
There are eight further normal coefficients to settle after the
present lifting theorem; no infinite-order hypothesis is required.
This finite criterion proves neither the consistency of those eight
coefficient equations nor field degree4, surface finiteness or actual
compatible embedding. Those remain separate parts of the original goal.

## Verification and remaining scope

Run from the repository root:

```sh
sage -python research/jacobian-followup/affine-complements/verify_full_next_normal.py
.venv/bin/python research/jacobian-followup/affine-complements/audit_pole_integration.py
```

The [checker](verify_full_next_normal.py),
[arithmetic circuits and exact residue/rational certificates](full_next_normal_certificates.json)
and [receipt](full_next_normal_receipt.json) retain all pivot choices,
unit determinants, complete matrix constructions, cusp controls and
following normal image spaces. The characteristic-zero local-ring
implication, open-chart dimensions and curve geometry are written
proofs. Shared Cox/Jacobian utilities are used. No formal verification,
independent whole-proof implementation, external review, new imported
theorem or novelty priority is claimed.

The first open reduced pairs remain (5,9),(6,9),(7,9),(8,9). The full
next-normal lifting test now passes all sixteen noncontact marked
candidates. All following Jacobian coefficients, the remaining contact
cases, other height-nine pairs, higher heights, all-height pole-four,
field degree four, global surface finiteness and actual compatible
embedding remain unresolved. The full original goal stays active.
