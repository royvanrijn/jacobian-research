# Three-point translation quotient fields and a Morse countercontrol

## 1. Results and scope

This note proves HC4MYKG0/1/2. The all-degree result concerns the actual
three-point suspension from [JQ](TWO_NORMAL_MODELS_AND_RELATIVE_RANK_ONE_ALIGNMENT.md):

    R=C[P,V,u,v], h=S_rho=2PV+B_rho(u,v),
    T=uv+rho, B_rho=u(T^2-1)+T.

For every complex rho, every nonzero action in either standard translation
module has a fixed affine line on EVERY noncritical whole value fibre.
After dividing the coefficient gcd, the action is primitive on the entire
source. Consequently [JT0](PRIMITIVE_QUOTIENTS_PLINTH_AND_GLOBAL_ACTION_CERTIFICATES.md)
excludes every action with the same embedded invariant field, including
arbitrary invariant replicas and h-preserving source conjugates. This
excludes these quotient-field classes from pure source alignment with
[KB's R_rho](GLOBAL_ORIENTED_NORMAL_FORMS_AND_SEIFERT_BARRIER.md), independently
of the existence of a three-point Keller fibre E3. No degree bound occurs.

An actual primitive quadric action also shows that a Morse invariant,
polynomial kernel, generic freeness, vanishing linearization and even
length-one plinth fibres away from one value do not imply global base
coordinate recovery. Its plinth support has a component contained in that
value. It does not meet the fixed-cylinder hypothesis of
[KF](FLAT_PLINTH_FAMILIES_AND_ONE_VALUE_RECOGNITION.md).

The three-point action classification and source alignment remain open.
The fixed-cylinder minimum-length question remains open. The countercontrol
does not disprove that more restrictive conjecture or settle either old
universal inversion endpoint. These are structural proofs, not searches.

## 2. An actual tangent-identity Morse action: HC4MYKG0

Work on the actual polynomial source R0=C[x,y,z,w]. Put

    q=xy+zw, k=yz-qw,
    D(y)=0, D(z)=-qy, D(w)=-y^2, D(x)=yz+qw.

Define commuting derivations D1,D2 by

    D1(x)=z, D1(w)=-y, D1(y)=D1(z)=0,
    D2(x)=w, D2(z)=-y, D2(y)=D2(w)=0.

They are triangular LNDs, commute on all four generators, and both fix y,q.
Thus D=yD1+qD2 is locally nilpotent: its iterates are the binomial expansion
in commuting locally nilpotent operators with invariant coefficients.
Alternatively D^2(x)=-2qy^2, D^3(x)=D^2(z)=D^2(w)=0.
Direct differentiation gives Dq=Dk=0. The Hessian determinant of q is 1;
its only critical point is the origin. Every coefficient of D has order
at least two there, so the action has identity differential at that point.

The action is primitive. A common prime factor of D(z),D(w) must divide y,
since gcd(q,y)=1. But D(x) modulo y is zw^2, which is nonzero.

### The actual kernel, including denominator removal

On y!=0 there is an exact chart with polynomial inverse:

    (R0)_y=C[y,y^-1,q,k,w],
    z=(k+qw)/y, x=(q-zw)/y, D=-y^2 partial_w.

Hence its localized kernel is C[y,y^-1,q,k]. Modulo y, q=zw and k=-zw^2
are algebraically independent: their Jacobian in z,w is -zw^2, nonzero.
In particular C[q,k] injects into R0/(y). If an invariant in R0 had a
minimal expression y^-m a(y,q,k), m>0 and y not dividing a in the polynomial
ring C[y,q,k], regularity in R0 would make a vanish modulo y. The injection
contradicts that minimality. Therefore

    ker D=A=C[y,q,k].

This is an equality of embedded rings, not a fraction-field assertion.

### Radical plinth and its value localization

Write J=D(R0) intersect A. The exact slice identities are

    D(-w)=y^2, D(-qx-z^2)=qk.

Therefore V(J) is contained in V(y,qk), a union of two lines. At y=0 the
quotient map has q=zw, k=-qw. No quotient point with q!=0,k=0 is reached,
and no point with q=0,k!=0 is reached. Every point outside V(J) IS reached:
for any slice value g=D(s) in A with g nonzero there, the local-slice
identity (R0)_g=A_g[s/g] supplies its affine-line fibre. Thus both missing
open lines lie in V(J); taking closures gives the reverse inclusion and

    sqrt(J)=(y,qk)=(y,k) intersect (y,q).

The line V(y,q) is contained in the single value q=0. Consequently this
support is not I_Z C[y,k,q] for a fixed finite Z in Spec C[y,k]. Indeed its
image in Spec C[y,k] contains an entire line, whereas a fixed cylinder has
finite image. The order and placement of the kernel generators matter.

After inverting q the further identity D(-z/q)=y gives

    J_q=(y,k) A_q.

Here containment of (y,k) follows from that identity and qk in J, while
the radical already equals (y,k), proving equality. The action is free
on every q=kappa!=0: D(w)=0 forces y=0, and then D(x)=qw cannot vanish
when q=zw!=0. On each such value the ACTUAL ideal fibre A/(J,q-kappa)
has length one. Nevertheless qk belongs to J and k does not, so A/J has
nonzero q-torsion and is not flat over C[q]. The support hypothesis of KF
is exactly what this example lacks; specialization here is of the actual
ideal, without substituting the plinth of a specialized action.

Finally B=C[y,k] is not a coordinate plane in this source. At y=0,k=lambda
with lambda!=0 the whole fibre has ring

    C[x,z,w]/(zw^2+lambda)=C[x,w,w^-1].

It has a nonconstant unit and cannot be A2. If R0=B[U,W], every such fibre
would be A2. Also d(y,k) has rank one at the Morse point. Thus primitive
LND plus Morse invariant does not force a nonzero linear action or a
transverse base projection. The polynomial q itself already has an
invertible gradient; this example is not an HC4 counterexample or a
solution of the universal quadratic endpoint.

## 3. Every standard S_rho module has fixed lines: HC4MYKG1

The two commuting h-preserving LNDs fixing P are

    D_u(u)=2P, D_u(v)=D_u(P)=0, D_u(V)=-B_u,
    D_v(v)=2P, D_v(u)=D_v(P)=0, D_v(V)=-B_v.

Their coefficients and commutation follow directly from differentiating
B; their common kernel is C[P,h], as retained in
[JS](GLOBAL_ADDITIVE_ACTIONS_AND_THE_K2_TRANSLATION_BARRIER.md).
For arbitrary a,b in C[P,h], not both zero, set

    D=a(P,h) D_u+b(P,h) D_v.

Invariant coefficients and commutation again prove local nilpotence.
Interchanging P,V gives the second whole module. It suffices to handle
the first.

Fix ANY kappa outside {1,-1,rho}; these are all critical values, with
repetitions allowed. Set alpha=a(0,kappa), beta=b(0,kappa). At P=0,h=kappa
the whole source is the curve B=kappa times the affine V-line. By the
actual JQ whole-curve isomorphism it has parameter

    t=T in P1 minus {-1,1,kappa,infinity},
    d=t^2-1, u=(kappa-t)/d, v=(t-rho)d/(kappa-t).

For completeness, B=kappa implies u(t^2-1)=kappa-t. If t=+-1 then
kappa=t, excluded. If t=kappa then u=0 and t=rho, also excluded. Conversely
the displayed formulas satisfy T=t and B=kappa and are regular at every
allowed t. Thus this chart omits no affine curve point or component.

On this exceptional cylinder the only possible D velocity is

    D(V)=-(alpha B_u+beta B_v).

A zero of that direction on the base curve makes its entire affine
V-line fixed. Define

    E(t)=2kappa t-t^2-1,
    A(t)=(t^2-1)(kappa-t)+(t-rho)E(t)
        =-2t^3+(3kappa+rho)t^2-2kappa rho t+rho-kappa.

Substitution gives B_u=A/(kappa-t), B_v=(kappa-t)E/d^2. The fixed-point
numerator with nonvanishing denominator (kappa-t)d^2 is

    F(t)=alpha A(t)d^2+beta(kappa-t)^2 E(t).             (1)

If alpha=beta=0, the whole exceptional cylinder is fixed. Otherwise:

* If alpha beta!=0, F has degree seven and leading coefficient -2alpha.
  At t=kappa it is alpha(kappa-rho)(kappa^2-1)^3, nonzero. At r=+-1
  it is beta(kappa-r)^2 E(r), also nonzero, since
  E(1)=2(kappa-1), E(-1)=-2(kappa+1). Every finite root of F is therefore
  allowed, and there is at least one over C.
* If alpha=0,beta!=0, use a root of the quadratic E. Its values at 1,-1
  above and E(kappa)=kappa^2-1 are nonzero, so all its roots are allowed.
* If alpha!=0,beta=0, use the cubic A. Its values are

      A(kappa)=(kappa-rho)(kappa^2-1),
      A(1)=2(kappa-1)(1-rho),
      A(-1)=2(kappa+1)(1+rho).

  If rho!=+-1, every root is allowed. In the two collision cases remove
  the single forbidden linear factor:

      rho=1: A=(t-1) Q+, Q+=-2t^2+(3kappa-1)t+kappa-1,
      Q+(1)=4(kappa-1), Q+(-1)=-2(kappa+1),
      Q+(kappa)=kappa^2-1;

      rho=-1: A=(t+1) Q-, Q-=-2t^2+(3kappa+1)t-kappa-1,
      Q-(1)=2(kappa-1), Q-(-1)=-4(kappa+1),
      Q-(kappa)=kappa^2-1.

  The residual quadratic still has roots and none is forbidden.

In every case there is a fixed affine line on that entire regular fibre.
This proves the theorem for every coefficient degree and every rho,
including both collisions. The fixed numerator degree seven is unrelated
to the generic degree or exceptional fibre cardinality of a Keller map.

## 4. Primitive normalization and whole quotient fields: HC4MYKG2

First suppose a,b are coprime in C[P,h]. On P!=0 the actual source chart is

    R_P=C[P,P^-1,h,u,v], V=(h-B)/(2P).

Here the derivation is 2P(a partial_u+b partial_v). Coprimality of a,b is
preserved under this localization and polynomial extension. Any prime
common divisor of its coefficients in the whole source must therefore
be P. It remains to exclude that divisor; one cannot assume coefficient
coprimality survives an arbitrary source embedding.

If P divided D(V), the polynomial

    a(0,B) B_u+b(0,B) B_v

would vanish identically in C[u,v]. Restrict to every generic B=kappa
curve in Section 3. Formula (1) would be identically zero as a polynomial
in t. Its degree-seven coefficient forces a(0,kappa)=0; the remaining
nonzero degree-four term forces b(0,kappa)=0. Since this holds on an
infinite value open, both coefficient polynomials vanish identically at
P=0. That makes P divide both a,b, contradicting coprimality. Thus D is
primitive on the ENTIRE source.

For arbitrary a,b, remove their gcd g in C[P,h]. The remaining D0 is again
an action in the same module, is primitive by the preceding proof, and
has the same embedded invariant field as D=gD0. By JT0, any other nonzero
LND E on R with that SAME EMBEDDED field is

    E=fD0, f in ker D0.

For clarity, the retained argument is elementary: equality of invariant
fields makes the derivations rationally proportional in the one-dimensional
space of derivations over that field. A denominator of the rational ratio
would divide all primitive coefficients of D0, so the ratio is polynomial.
The LND degree rule applied to E(f)=fD0(f) forces D0(f)=0. Invariant replicas
can only add fixed points and thus retain all of Section 3's fixed lines.
If E is primitive too, f is a nonzero constant. Interchange P,V and
transport by any actual h-preserving polynomial source automorphism to
obtain the same conclusion for both whole conjugacy-saturated classes.

### Why this excludes pure source alignment, before E3

In the explicit KB target

    R_rho=a3(z) X P+C3(z) Q+D_rho(z),
    C3=(z-1)(z-2)(z-3), a3 C3'=1 modulo C3,
    D_rho(1)=1, D_rho(2)=-1, D_rho(3)=rho,

the coefficient translation E_X=C3 partial_P-a3 X partial_Q is an LND.
It is free on every regular value fibre: a fixed point has C3=0 and X=0
(a3 does not vanish at these roots), hence value in {1,-1,rho}. An actual
source equivalence S_rho~R_rho would transport this to an h-preserving LND
free on EVERY noncritical whole fibre. Its embedded invariant field cannot
belong to either class just excluded. This argument does not assume E3
or an actual unit-Hessian source member.

A remaining action must use a genuinely different embedded quotient field;
changing kernel generators, multiplying action speed by an invariant, or
conjugating within the h-preserving source group cannot evade this theorem.
We have not classified every h-preserving LND, every action fixing P or V,
or actions with extra-invariant coefficients outside these standard modules.
An abstractly isomorphic field or polynomial kernel is not the embedded-field
hypothesis. Neither S_rho alignment nor its impossibility has been proved.

## 5. Literature, reproducibility and the remaining gate

The publisher's [Bisi--Furter--Lamy abstract](https://numdam.org/articles/10.5802/jep.8/)
concerns the tame automorphism group of the affine quadric SL2 and its square
complex. Only that abstract and bibliographic page were inspected for this
note; no full-paper automorphism theorem or classification of all additive
actions is imported. Its scope gives no missing all-action classification
for this A4 problem. New proofs here use the displayed identities and
retained JQ/JS/JT/KB/KF statements; there are zero new external proved imports.

The [checker](verify_three_point_translation_barrier.py),
[receipt](receipt-three-point-translation-barrier.json) and
[source context](sources/three-point-translation-barrier-context.json) bind
the retained inputs and check the kernel charts, tangent identity,
primitivity controls, plinth slices, full rational curve, every puncture
test, collision factors and target action. These algebraic controls support
the all-degree arguments above; they do not enumerate actions or certify an
unproved classification. Independent replay, formal verification and external
review are not claimed.

The useful next gate remains one new native quotient: actual embedded
ker(delta)=C[b1,b2,h], fixed finite cylindrical radical plinth support and
one ACTUAL plinth-ideal value fibre of minimum length three. KG2 requires a
different embedded field from all the standard-module classes above. KG0
prevents a blanket appeal to a primitive Morse action's first linear jet;
any minimum-length proof must use the fixed-cylinder hypothesis itself.
E3 remains independent of pure alignment, and the two universal endpoint
inversions, E2, E3, full embedding recovery and unrestricted HC4/JC2 remain open.
