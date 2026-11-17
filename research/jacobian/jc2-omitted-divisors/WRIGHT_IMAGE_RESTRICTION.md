# Image and exceptional-length accounting on the marked Wright model

14 November 2026. PF2OMIT9 is conditional on the **actual** marked model
of PF2WRIGHT1 and a finite target map. It does not supply that map or a
compatible completion for an arbitrary Keller normalization.

## 1. The input and its restriction to the omitted prime

Write A=C[Y,A1,B1,D]=C[Y,XY,X^2Y,X^3Y+X] inside C[X,Y], with the
relations and source open proved in
[PF2WRIGHT1](../affine-complements/WRIGHT_MARKED_MODEL.md). Its sole omitted
prime R=(Y,A1,B1+1) is A1_v. A neighborhood of the entire R has coordinates

    Y=-u^2+u^3 v, A1=-u+u^2 v, B1=-1+uv, D=v, X=1/u.

In particular dX wedge dY=-u du wedge dv. Suppose P,Q belong to A,
their source Jacobian is the nonzero constant c, and
pi=(P,Q):Spec(A)->A2 is finite, of field degree d. Then

    J_(u,v)(P,Q)=kappa u, kappa=-c !=0.                 (W1)

Every ramified prime is omitted because the actual source is Keller.
Here R is the only omitted prime, and (W1) gives e_R=2. PF2OMIT7 gives
f_R=1, including the residue extension rather than just a chosen local
parametrization. Thus its target image C has birational parametrization

    gamma(v)=(p0(v),q0(v))=(P,Q)(0,0,-1,v).           (W2)

The image polynomial is the irreducible generator of the kernel of
C[U,V]->C[v], U->p0,V->q0. It is not specified until P,Q are specified.

## 2. Exact jet restrictions at every pole height

For any representative of P in C[Y,A1,B1,D], take the following partial
derivatives at (Y,A1,B1,D)=(0,0,-1,v):

    p1=-P_A1+v P_B1,
    p2=-P_Y+v P_A1+P_A1A1/2-v P_A1B1+v^2 P_B1B1/2.

Taylor substitution gives P=p0+u p1+u^2 p2+O(u^3). Define q1,q2 likewise.
These expressions do not depend on the representative modulo the defining
ideal: they are coefficients of the same function in the actual boundary
chart. This is an all-height Taylor identity, not a bounded search.
The coefficients of u^0 and u^1 in (W1) give

    p1 q0'-p0' q1=0,
    2(p2 q0'-p0' q2)+p1 q1'-p1' q1=kappa.             (W3)

Finiteness prevents contraction of R, so p0',q0' are not both zero.
Let h be their monic gcd, and put a=p0'/h, b=q0'/h. Since (a,b)=1,
the first equation gives a polynomial lambda with
p1=lambda a, q1=lambda b. The second becomes

    2h(p2 b-a q2)+lambda^2 W=kappa,
    W=a b'-a' b.                                    (W4)

In particular lambda and W are units modulo h. This is a practical
necessary condition on restrictions of any proposed pair in the marked
ring, including pairs at unbounded boundary pole height. It does not
certify finiteness or produce such a pair.

## 3. Degree and all exceptional fibre lengths

PF2OMIT7's simple-boundary formula applies to this actual image:

    ell_v=2+ord_v h, E_R=deg h.

The all-prime Euler identity of PF2OMIT2 now gives

    d-1=2+deg h, hence d=3+deg h.                     (W5)

Thus for a finite **quartic** map h=v-v*. There is exactly one
nonimmersive parameter v*, and PF2OMIT7 makes it an ordinary (2,3) cusp
with ambient algebra C[w]/(w^3). All other points of R are folds with
algebra C[w]/(w^2). In (W4), W(v*)!=0 and

    det(gamma''(v*),gamma'''(v*))=2W(v*),
    lambda(v*)^2 W(v*)=kappa.                         (W6)

This determines the exceptional target value from a proposed pair:
y*=gamma(v*). It is not a value computed for a presently supplied global
pair. If b(y) is the number of normalization parameters of C over y,
finite-flat length conservation and the etale source give

    #F^-1(y)=4-2b(y)-1_(y=y*).                       (W7)

Nonnegativity forces b(y*)=1 and b(y)<=2 elsewhere. In particular the
cusp has no additional image branch through it. The complete conditional
target-value table is:

| Target value | Retained source points | Omitted local factors |
|---|---:|---|
| Outside C | 4 | None |
| One immersed image branch | 2 | C[w]/w^2 |
| Two immersed image branches | 0 | Two separate copies of C[w]/w^2 |
| y*, the unique cusp value | 1 | C[w]/w^3 |

The two-branch row is present only if the actual image has a collision.
No nodality or transversality is assumed there. Values with three
branches, or with a cusp branch and another branch, are impossible.
There are no further exceptional ambient lengths on R. In particular
every such quartic model automatically has a reduced singleton fibre
over its unique cusp value; an off-branch singleton is impossible.

After blowing up that value, the omitted radial packet is the cusp
packet (e,f)=(3,1), with length3 in every exceptional direction, plus the
retained radial packet (1,1). Its principalization canonical vector is
(2,4,5), as in PF2OMIT7. The original ramified R still has generic index2;
these indices belong to different target models.

## 4. Leading coefficient and compulsory index at the cusp

Choose a generic admissible source coordinate T=X+tY by PF2OMIT3. On R
its pole order is exactly one: X has order -1 and Y has order2, so they
cannot cancel. There is only one omitted prime, with f=1. Consequently
the primitive inverse leading coefficient satisfies

    a_d=constant*H_C,                                (W8)

where H_C is the reduced irreducible image equation. In degree4 its
sole positive Newton segment over generic C has width2 and slope1/2.
The leading coefficient's vanishing order1 records pole charge; it does
not replace the missing generic length2. PF2OMIT4 also gives

    Disc(R_d)=constant*H_C J^2.                       (W9)

For the quartic cusp y*, J(y*)=0. To see this, use the reciprocal chart
z=1/T. The boundary chart gives z=u*unit. In the normalized cusp fibre
C[w]/w^3, a defining equation u of the smooth critical curve has order2:
this follows from the regular ordinary-cusp form proved in PF2OMIT7.
Hence z lies in (w^2) and cannot generate C[w]/w^3. The binary-root
order's length3 infinity factor is generated by z. Its normalization
therefore cannot be an isomorphism over y*. If the normalization-index
determinant J were a unit at y*, the inclusion of the two finite locally
free target algebras would be an isomorphism there, a contradiction.
Thus

    ord_(y*) Disc(R_4)>=2+2=4.                       (W10)

Here H_C has ordinary-cusp multiplicity2. This is a necessary index
contribution at an on-branch singleton, distinct from the off-branch
singleton staircase in PF2OMIT4/6. Specialization R_4 at y* has degree1
by PF2OMIT3, but that degree alone does not identify the normal cusp
factor's generator or its index.

An exact local control displays this mechanism. Put

    U=u+w^2, V=u*w+(2/3)w^3, z=u.

Its Jacobian is u, normal fibre at zero is C[w]/w^3, and

    w^3-3U*w+3V=0,
    G(z)=z^3+3U*z^2+9V^2-4U^3=0.

In the normal basis (1,w,w^2), (1,z,z^2) has inclusion determinant -3V.
The discriminants are -27(9V^2-4U^3) and
-243V^2(9V^2-4U^3), respectively. At zero z=-w^2. Thus the normal
length3 algebra and its nonnormal monogenic presentation must be kept
separate. This is a finite cubic local control, **not** a quartic map
on the Wright ring or a global Keller realization.

The [replay](verify_wright_image.py) checks the explicit chart, invariant
Taylor coefficients, universal jet identities, cusp/index control and
finite fibre rows. Finiteness, global Euler accounting and the imported
simple-ramification theorem remain written mathematical arguments.
Actual P,Q, H_C and y* are still absent; no compatible-completion,
quartic-exclusion or unrestricted singleton-zero-properness theorem follows.
