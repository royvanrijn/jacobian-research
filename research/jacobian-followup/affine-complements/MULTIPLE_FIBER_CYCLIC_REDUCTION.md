# Cyclic reduction of the multiple omitted fiber

14 November 2026. PF2RULED37 gives a structural reduction in every
finite field degree. It does not bound the multiplicity by a calculation.

For a smooth one-boundary normalization with generic index two and
omitted fiber qR, cyclic base change separates q reduced boundary lines.
Each single-boundary piece is an affine-line torsor with linear bundle
O(-(q+2)). The global volume test detects the same one translation
coefficient on every piece, independent of q. The lifted finite map
has wreath-product monodromy and an additional branch curve, rather
than the simple symmetric cover used in the quartic exclusion.

An explicit smooth affine family with every q>=2 has K=R and an
exact source volume form. Thus those geometric and differential
conditions cannot alone remove the multiple-fiber gate. No member
is asserted to have a finite Keller map.

## 1. Hypotheses, and coefficients which must be retained

Let N be a smooth affine surface with U=A2_(x,y) and N minus U=R=A1.
Suppose Sigma:N->P1 restricts to x on U and

    Sigma^*(infinity)=qR,  q>=1,
    div_N(dx wedge dy)=(e-1)R,  e>=2.                    (1)

The index e is kept separate from q. In the actual finite normalization
of a Keller map with one generic-index-two omitted prime, e=2;
PF2OMIT7 supplies the required smoothness and residue degree one.

The divisor localization sequence gives

    Pic(N)=Cl(N)=Z[R].                                  (2)

Surjectivity follows from Cl(U)=0. If a rational function had divisor
supported on R, its restriction would be a unit on the polynomial
plane, hence a nonzero constant. This proves injectivity as well.
The same proof works for any normal plane completion with finitely
many omitted divisors: its class group is free on those divisors.
Smoothness is needed only to identify Pic with Cl.

In particular the canonical class in (1) is (e-1)[R]. Knowing just the
ramification support does not set its coefficient to one; no claim
that the primary-module condition forces e=2 is used.

## 2. Normalize the cyclic base change

Put x=z^q on the source plane and let

    r:N' -> N

be the normalization in C(z,y), z^q=x. It is a connected finite
cyclic cover of degree q: x has order one along the retained closed
line L=Sigma^*(0)={x=0} in U and so is not a proper power.
Over U, N' is exactly A2_(z,y). Write L'={z=0}.

Near every point of R a local parameter s gives 1/x=s^q*a, with a
a unit. After an etale unit-root extension, normalized cyclic base
change is q disjoint copies of this smooth chart. Thus r is etale
along R. It is also etale outside L, and is locally z->z^q along L.
Its inverse image of R is a finite etale degree-q cover of A1.
Finite etale covers of complex A1 are trivial. Consequently

    r^*R=R_1+...+R_q,  each R_i=A1 mapping isomorphically to R,
    Ram(r)=(q-1)L',   div(z)=L'-sum R_i.                 (3)

The deck group cyclically permutes the q lines. This is not an
everywhere etale cover of the simply connected N.

For each i remove all R_j with j!=i and call the resulting open
surface X_i. It is affine: those R_j are effective Cartier divisors
on the affine smooth N', so their complement is an affine open
morphism, locally a principal open on its target. The induced
ruling X_i->P1_z has reduced omitted fiber R_i and its other fibers
are A1. Borisov--Gabber--Vasiu Theorem11.2(2), with multiplicity one,
makes it locally affine-line trivial. Its affine transition
defines a torsor under a line bundle on P1.

The primary manuscript and precise import are retained in
[sources/multiple-fiber-cyclic/manifest.json](sources/multiple-fiber-cyclic/manifest.json).
The local splitting calculation above also explains the normalization;
no algorithm for a generic hypothetical Keller polynomial is assumed.

## 3. The torsor index is q+e, not three

Keep the actual common source coordinates z,y, with no separate
coordinate change on each X_i. Its second trivialization has the form

    u=1/z,   v_i=z^m_i*y+g_i(z),   g_i in C[z,z^-1].      (4)

The coefficient of y can be normalized to the displayed monomial by
a constant rescaling. The affine part is allowed to have high powers:
one cannot normalize it independently on each piece while retaining
the common source y and its deck action.

On this chart

    dz wedge dy=-u^(m_i-2) du wedge dv_i,
    r^*(dx wedge dy)=q*z^(q-1) dz wedge dy
                   =-q*u^(m_i-q-1) du wedge dv_i.       (5)

As r is etale at R_i, (1) says that the last exponent is e-1.
Therefore every i has

    m_i=q+e.                                           (6)

For index two these are O(-(q+2))-torsors. Independently,

    Pic(N')=direct sum_i Z[R_i],
    [K_N']=(q+e-2)*sum_i [R_i].                         (7)

Indeed the ramification formula gives
(e-1)*sum R_i+(q-1)L', and div(z) identifies [L'] with sum[R_i].
The unweighted source form dz wedge dy represents the same class.

All finite intersections of distinct X_i are the contractible
source plane U'. Each X_i is smoothly a line bundle over P1,
hence has the homotopy type of S2. Mayer--Vietoris and van Kampen
then give H1(N',Z)=0 and H2(N',Z)=Z^q, with the corresponding
restriction isomorphism H2_deRham(N')->direct sum H2_deRham(X_i).
Also chi(N')=q+1, consistent with the finite-cover calculation
q*chi(N)-(q-1)*chi(L)=2q-(q-1).

These equalities do not identify the end pair or its intersection
pairing; no lens-space conclusion follows from them alone.

## 4. One shared coefficient detects the volume class

Write a_i for the coefficient of z^e in g_i. On U' choose the
primitive theta_U=-q*z^(q-1)*y dz of r^*(dx wedge dy).
On the second plane chart choose theta_i=q*u^(e-1)*v_i du.
Their difference is

    theta_U-theta_i=-q*u^(e-1)*g_i(1/u) du,
    Res_(u=0)(theta_U-theta_i)=-q*a_i.                   (8)

The two charts are A2 and their overlap is Gm x A1. A Laurent
one-form in u is exact on this overlap exactly when its residue
vanishes. The polynomial Poincare lemma and the two-chart de Rham
sequence show that (8) detects the entire H2 class on X_i,
not just a sufficient local obstruction. For example the Cech
primitive of a residue-zero difference can be split across the
charts because X_i is affine and H1(X_i,O)=0.

Choose deck generator z->zeta*z, where zeta^q=1, and index the R_i
cyclically. In the normalized coordinates (4), its action implies

    g_(i+1)(zeta*z)-zeta^(q+e)*g_i(z)

has only nonpositive z powers; the difference comes from a function
regular on the second base chart. Taking the z^e coefficient gives

    a_(i+1)=zeta^q*a_i=a_i.                             (9)

Thus there is one shared scalar a, not q independent period tests.
Changing y by a common polynomial in z adds terms of degree at least
q+e to g_i; changing the second trivialization adds nonpositive
powers. Neither changes vanishing of this coefficient.

The volume form on N is exact if and only if a=0. Necessity follows
by pullback and restriction. Conversely, (8) and the de Rham
restriction isomorphism make the pulled-back class zero when a=0.
N' is affine, so it has a global algebraic one-form primitive.
Average that primitive over the cyclic group. The invariant form
descends regularly: over the etale locus this is ordinary descent;
near L' an invariant f(z,y) dz has f divisible by z^(q-1) and is
a regular multiple of dx, while invariant dy terms descend directly.
This gives a global primitive on N.

In the actual finite Keller setting dx wedge dy=c^-1*dP wedge dQ
is already exact. Hence

    [z^e]g_i=0 for every i; for index two, [z^2]g_i=0.    (10)

The exponent tested is independent of the omitted-fiber multiplicity.
This is the multiple-fiber counterpart of the retained
[volume-period obstruction](VOLUME_PERIOD_OBSTRUCTION.md).

## 5. The lifted finite map acquires full wreath monodromy

Return to the actual finite normalization of a Keller map
F=(P,Q) of field degree d with precisely one index-two omitted prime.
Its boundary is the entire exceptional curve C, with simple
meridians generating S_d. These deductions in PF2RULED32 use only
the one-simple-prime hypotheses, not its Wright chart.

For clarity, C necessarily has a finite nonimmersion. Chau's
Theorem1.1 forces a collision of its finite birational A1
parametrization, since otherwise C is homeomorphic to C and F
would be bijective. Each collided boundary point costs at least
two sheets, so d>=4. The index-two local dictionary and Euler
identity give

    deg gcd(p',q')=d-3>=1.                              (11)

The line L is closed in N and disjoint from R. Its finite image
D=pi(L) is an immersed polynomial curve. Indeed L->D is finite;
its factorization through the normalization D^nu=A1 is a polynomial
map with no critical point, since F is etale on L. That polynomial
is affine linear. Hence D^nu->D is everywhere immersive, and
D cannot equal the nonimmersed C. Self-collisions of D are allowed.

For q>1, the composite N'->A2 is finite of field degree qd. Its
entire finite branch locus is C union D. At a generic point of C
its meridian is q disjoint transpositions. At a generic point of D
it is one q-cycle and fixes the other q(d-1) sheets. Its source
plane restriction is

    G(z,y)=(P(z^q,y),Q(z^q,y)),
    J(G)=c*q*z^(q-1),                                  (12)

not a Keller map for q>1.

The intermediate cyclic cover gives a block system of d blocks
of q sheets, with monodromy contained in C_q wreath S_d.
The D meridian gives a generator of C_q in one block. Conjugating
by lifts of the transitive S_d action gives the generator on each
block, independently, hence the whole base subgroup (C_q)^d.
The quotient is S_d. Therefore the monodromy is exactly

    C_q wreath S_d,  of order q^d*d!.                   (13)

This identifies the new cover problem without treating the qd
cover as a simple symmetric one. Restricting to one X_i also
removes the other q-1 boundary primes and loses finiteness over
the original target when q>1. Neither operation transfers the
marked quartic exclusion directly.

## 6. Explicit affine controls for every multiplicity

For any integer q>=1 glue U=A2_(x,y) to W=A2_(u,t) by

    y=u,  x=u*t+u^(-q),  u!=0.                         (14)

Here is an explicit affine ring for this smooth surface N_q.
On U put

    f=x*y^q-1,  a=x*f,  b=x^2*y^(q-1)*f,
    S_q=C[y,f,a,b] inside C[x,y].                       (15)

All four functions extend to W:

    y=u, f=u^(q+1)*t,
    a=u*t*(1+u^(q+1)*t),
    b=t*(1+u^(q+1)*t)^2.                               (16)

After localization at f, x=a/f recovers the source chart D_U(f).
After localization at f+1, u=y and t=b/(f+1)^2 recover D_W(f+1).
These are precisely the two charts of the glued surface on those
opens. Since (f+1)-f=1, they cover Spec(S_q). Thus the gluing is
separated and affine, is Spec(S_q), and has sole complement
R={u=0}=A1. It is smooth by its two plane charts.

The projection x extends to P1 with infinity fiber qR; at R its
numerator 1+u^(q+1)*t is a unit. Moreover

    dx wedge dy=-u du wedge dt,  K_N_q=R.               (17)

For q>=2 a global regular primitive is

    theta|U=(q*x dy+y dx)/(q-1),
    theta|W=((q+1)*u*t du+u^2 dt)/(q-1),
    d theta=dx wedge dy.                               (18)

These formulas agree exactly. At q=1 the local primitives instead
differ by du/u, with nonzero residue, so the form is not exact.

The normalized x=z^q cover of this control has q components over R.
Writing zeta_i for the q-th roots of unity, its pieces have

    u_i=1/z, v_i=z^(q+2)*y-zeta_i*z^(q+1).              (19)

Near R_i, z*y approaches zeta_i. The defining identity
(zy)^q-1=y^(q+1)*t shows that (19) is regular there and
v_i=zeta_i^2*t/q on R_i. It gives the whole second plane chart,
including u_i=0; on u_i!=0 it is the exact inverse of (19).
Thus g_i=-zeta_i*z^(q+1). Its z^2 coefficient is nonzero only
when q=1, agreeing with (18) and the uniform criterion.

These are geometric controls with the correct smoothness,
one-boundary class, canonical divisor and arbitrary multiple fiber.
For every q>=2 they also pass the exact-volume test. They are
not asserted to arise as finite normalizations of Keller maps,
nor to exhaust the possible normalizations or end pairs. No assertion
here rules out another source-coordinate ruling of a different
multiplicity on N_q; the controls refute a forced multiplicity-one
inference for the given ruling, not an unproved existential statement
about every alternative coordinate system.

## Consequence for the research gate

Cyclic normalization and the single invariant residue replace an
unstructured multiple-fiber question by q torsor charts with fixed
linear degree q+2, cyclic gluing and an auxiliary immersed branch
curve with full wreath monodromy. The residue obstruction alone
does not exclude q>=2; the controls prove that limitation.

Next use the actual finite algebra, the cyclic gluing or the marked
end pair to constrain the remaining translation data. Full quartic
packets and higher finite Wright field degrees remain open. No
coefficient/Groebner search, generic surface classification, end
recognition, formal verification, independent full reproof, external
review or novelty priority is claimed.
