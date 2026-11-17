# Regular source coordinates, canonical zeros and primitive residues

For every n>=2 and q>=1 there is an explicit smooth affine surface N_(n,q)
containing the WHOLE affine plane X=A2_(x,y) as an open, with complement
one affine line E. The ACTUAL source coordinate y is regular on N_(n,q),
and

    div_N(dx wedge dy)=q*(n-1)*E,
    Cl(N)=Z[E].

For q=1 the source volume has no global regular primitive. For EVERY
q>=2 it has an explicit global regular primitive. Thus smoothness,
affineness, a regular source coordinate, positive nontorsion canonical
class and exact source volume coexist. They do not alone settle the
one-integral-coordinate problem for an actual Keller target.

The same surfaces are the full finite normalizations of explicit
non-Keller polynomial plane maps. That finite map is part of the control,
not a claimed Keller counterexample. For q=1 the residue proof excludes
EVERY proposed regular Keller target pair on this entire surface.

## 1. An exact affine-plane completion

Inside R=C[x,y] set

    z=x*y^q,
    a_j=x^j*(z-1), 1<=j<=n,
    S_(n,q)=C[y,z,a_1,...,a_n], N_(n,q)=Spec S_(n,q).

This is a finite-type domain with fraction field C(x,y), because
x=a_1/(z-1). The two principal opens D(z-1) and D(z) cover N.
On the first,

    S[(z-1)^-1]=C[x,y,(x*y^q-1)^-1].                    (1)

On the second write w=a_n. Then

    S[z^-1]=C[y,z,z^-1,w]/(y^(qn)*w-z^n*(z-1)),
    a_j=y^(q(n-j))*w/z^(n-j).                           (2)

The equalities follow by substitution inside C(x,y). The displayed
hypersurface is a domain: its linear-in-w equation is primitive and
irreducible, since y and z^n*(z-1) are relatively prime in C[y,z,z^-1].
It therefore gives the whole localized ring, not only a birational
parametrization.

Both charts are smooth. In (2), if y!=0 the w derivative is nonzero.
If y=0 its equation forces z=1, where the z derivative is -1. Thus
N is a SMOOTH normal affine surface.

The given inclusion S subset R defines an open immersion of the whole
source plane. Its image is D(z-1) union D(y). On D(y) recover x=z/y^q;
together with (1) these are the two inverse source charts. Their complement
is exactly

    E={y=0,z=1,a_1=...=a_(n-1)=0}, w=a_n arbitrary, E=A1. (3)

No isolated point is omitted from this calculation. The fibre y=0 consists
of the whole source line L_0 and E, both reduced and disjoint. In particular
y is an ACTUAL polynomial source coordinate that extends regularly, not
a nonlinear invariant substituted for a coordinate.

## 2. Canonical order and the first varying coefficient

Put C(z)=z^(n-1)*((n+1)*z-n). Differentiating the relation in (2) gives

    C(z)*dz=qn*y^(qn-1)*w*dy+y^(qn)*dw.                   (4)

At E, C(1)=1. Since x=z/y^q,

    dx wedge dy=dz wedge dy/y^q
               =y^(q(n-1))/C(z)*dw wedge dy.             (5)

The right side is regular near the ENTIRE E. It is a unit times
y^(q(n-1)) in these smooth local coordinates. On X it is already a unit,
so its divisor is exactly q(n-1)E.

The free class-group calculation is elementary. Every divisor becomes
principal on the factorial plane open X, so E generates Cl(N). A relation
kE=div(f) makes f a source unit, hence a nonzero complex constant.
Therefore k=0 and Cl(N)=Z[E]. Formula (5) makes the canonical class
positive and nontorsion.

The boundary expansion reads

    x=y^-q+y^(q(n-1))*w+terms of higher y order.           (6)

The leading pole y^-q has constant coefficient along E. The first
coefficient varying with w has POSITIVE order q(n-1). It cannot be
replaced in the volume calculation by the leading pole order.
Here ord_E(y)=1 while ord_E(dx wedge dy)=q(n-1)>0.
Consequently an affine-plane completion with a regular source coordinate
need not have a nonpositive boundary volume order.

## 3. The q=1 family has no regular volume primitive

Assume a regular one-form theta on N has dtheta=dx wedge dy.
On the WHOLE source plane, theta-x*dy is a closed polynomial one-form,
so

    theta|X=x*dy+dF, F in C[x,y].                         (7)

The polynomial Poincare assertion follows by integrating its x
coefficient and then its remaining y-only coefficient.

Use (y,w) as smooth local coordinates at a generic point of E.
For q=1, expansion (6) makes x*dy have a simple logarithmic pole
of residue ONE. The differential of a rational function F has zero
logarithmic residue: in its Laurent series over C(w), differentiating
any power y^k gives k*y^(k-1), whose y^-1 coefficient is zero.
Thus dF cannot cancel that residue. Equation (7) contradicts regularity
of theta.

In particular there are NO P,Q in S_(n,1) with nonzero constant
Jacobian relative to x,y, in ANY degrees. Such a pair would give the
global regular primitive P*dQ divided by that constant. This excludes
a complete realization architecture, with its plane source, boundary,
normality and finite-type algebra fixed explicitly.

## 4. Every q>=2 family DOES have an exact regular source volume

For q>=2 define the rational one-form

    theta=(y*dx+q*x*dy)/(q-1)
          =dz/((q-1)*y^(q-1)).                           (8)

Its source expression is polynomial, and dtheta=dx wedge dy.
At E, (4) gives

    theta=[qn*y^(q(n-1))*w*dy
           +y^(q(n-1)+1)*dw]/((q-1)*C(z)).               (9)

Both exponents are nonnegative, and C(z) is a unit near E.
Thus theta is regular at EVERY point of E as well as throughout X.
It is a global regular primitive on N.

This supplies exact-volume controls with a regular actual source
coordinate, smooth affine normal completion, one affine-line boundary
and positive nontorsion canonical class, for all n>=2,q>=2.
An argument using only those properties cannot force zero boundary.
There is no regular finite Keller target pair supplied by (8).

## 5. An explicit full finite normalization and its failed Keller hypothesis

Take B=C[P,Q] subset S with

    P=y, Q=a_n=x^n*(x*y^q-1).                             (10)

These two elements are algebraically independent. The element z is
integral over B, satisfying the MONIC equation

    z^(n+1)-z^n-P^(qn)*Q=0.

Every a_j is integral over B[z], since

    a_j^n-Q^j*(z-1)^(n-j)=0.                              (11)

Thus S is finite over B. It is normal with fraction field C(x,y), so
it is the FULL finite normalization of B in that source field.
In particular the actual source coordinate y is integral over B
(and even belongs to B).

The induced whole-plane polynomial map is

    (x,y) -> (y,x^n*(x*y^q-1)),
    Jac=-x^(n-1)*((n+1)*x*y^q-n).                        (12)

This Jacobian vanishes on source curves. Hence the map is NOT Keller.
Its generic degree is n+1: after inverting y, the source variable x
satisfies a polynomial of degree n+1 over C(y,Q), irreducible because
Q=x^n*(x*y^q-1) is a degree-(n+1) rational function of the transcendental
x over C(y). Equivalently the monic z polynomial defines that extension.

At generic E, y=P has valuation one and Q=w is a separating tangent
coordinate, so this omitted prime is UNRAMIFIED over P=0. The finite
ramification meets X. Along E the polynomial source Jacobian in (12)
has order -q(n-1), precisely cancelling the zero in (5):

    dP wedge dQ=Jac_source*dx wedge dy

is a unit there. The target and source volume forms are therefore NOT
constant scalar multiples. Replacing (12) by a constant-Jacobian target
pair would be a substantial additional condition, not a consequence
of finiteness, the integral coordinate or exactness of (8).

### 5.1. Finite volume-preserving source translations need not descend

There is a separate all-degree obstruction to using a regular additive
action on the FULL finite normalization to manufacture a target action.
Let H(T) be ANY complex polynomial of degree d>=3, write H_t and H_tt
for its first and second derivatives, and set

    S=C[u,t], P=H_t(t)-u,
    Q=H(t)-t*H_t(t)+u*t=H(t)-P*t, B=C[P,Q].                  (13)

This is a finite map of smooth affine planes, of generic degree d, with
S the full finite normalization of B in C(u,t). Indeed u=H_t(t)-P, so
S=C[P,t]. The function Q is a degree-d polynomial in the transcendental
t over C(P), so P,Q are algebraically independent. The equation

    H(t)-P*t-Q=0                                             (14)

is monic after dividing by the nonzero leading coefficient of H.
It proves finiteness. As a polynomial in C[P,T,Q], its left side is
irreducible: its quotient eliminates Q and is the domain C[P,T]. It
is primitive over C[P,Q], with nonzero constant leading coefficient
in T. Gauss therefore proves irreducibility over C(P,Q), giving the
claimed degree d.

The ACTUAL source coordinate u is integral over B, because it is a
polynomial in the integral element t with coefficients in B. Define
the regular locally nilpotent source derivation D=partial_t, keeping
u fixed. Direct calculation gives

    Jac_(u,t)(P,Q)=-u,
    Omega=dP wedge dQ=-u*du wedge dt,
    Lie_D(Omega)=0, D(u)=0.                                 (15)

Thus the action preserves the entire ramification divisor u=0, the
pulled-back target volume form and the ordinary nowhere-zero source
volume du wedge dt. It is defined on the whole finite normalization,
not just its etale open or a formal boundary chart.

Nevertheless D does NOT preserve even the target FUNCTION FIELD:

    D(P)=H_tt(t) not in F=C(P,Q).                            (16)

The minimal polynomial of t has degree d, while H_tt is a nonconstant
polynomial of degree d-2. If H_tt(t) were in F, subtracting that field
element would give a nonzero polynomial of degree less than d vanishing
at t, a contradiction. Consequently D(B) is not contained in B. The
additive action cannot descend to the target. The degree restriction
matters: for a quadratic H, H_tt is constant and this argument fails;
in fact D(P) and D(Q) then do belong to B.

These covers are NOT whole-plane Keller maps: their ramification line
is inside the retained affine plane. They admit no dense affine-plane
open on which the target map is Keller either. Such an open would have
to avoid u=0, making the nonconstant regular function u a unit on A2;
its pullback would be constant, contrary to dominance of the open.
The etale open u!=0 is Gm times A1 and has nonconstant units.

Thus finiteness over a polynomial target, a smooth full normalization,
an actual integral source coordinate, an algebraic source translation,
and preservation of both volumes and ramification do not imply target
descent. A proof in the actual Keller setting must use its entire-plane
etaleness beyond these properties. This control does not refute the
one-integral-coordinate Keller criterion; its failed hypothesis is
explicit in(15).

## 6. Scope of the exclusion and remaining trace problem

The q=1 nonexact family cannot admit a regular Keller target pair.
The q>=2 family blocks the broader geometric shortcut after exactness
is added. This does not refute, prove or assume that one integral
source coordinate forces inversion for an ACTUAL Keller map.
In the controls, finite ramification lies on the retained source,
whereas the Keller condition would place every ramified prime outside it.
That indispensable difference remains visible in (12).

The construction, charts, smoothness, canonical class, primitive
obstruction and exact primitive have written proofs for all n,q.
The checker verifies independent source Jacobians, defining relations,
the chart volume calculation, positive/negative primitive controls, and
the all-polynomial chain-rule identities in(13)--(16). Finiteness and
failure of target-field descent follow from the written irreducibility
and minimal-degree proof, not a finite sample of H.
No new classification theorem is imported. These controls do not resolve
automatic trace integrality, HC4 or JC2, and no independent mathematical
review, formal verification or external review is claimed.
