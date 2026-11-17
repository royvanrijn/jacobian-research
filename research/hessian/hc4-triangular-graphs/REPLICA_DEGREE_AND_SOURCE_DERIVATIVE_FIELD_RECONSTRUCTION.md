# Replica degree and first source-derivative field reconstruction

HC4MYMR0. 16 November 2026. ALL actual complex polynomial plane Keller
maps, all generic degrees and all polynomial source-coordinate systems.
Full HC4/JC2 and the trace-dual/marked-energy implication remain OPEN.
No novelty, formal verification or external review is asserted. The
first-column field argument in Section 6 was independently checked internally
on 17 November 2026; this is not an independent replay of the whole note.

The new obstruction is an exact degree divisibility. For a source LND
delta=p(y)*partial_x and an intermediate field E stable under delta,

    [L:E] divides deg(p)+1.                                     (1)

A bare source translation therefore preserves NO proper intermediate
field. In fact its TWO FIRST derivatives of the target generators already
recover the ENTIRE source function field, without a differential-stability
assumption on their generated field (Section 6). This is field reconstruction, not
regular coordinate reconstruction or inversion. Nonconstant replicas
must not be divided out in an unknown intermediate field.

## 1. The actual smooth algebra in the intermediate field

Fix an ACTUAL complex plane Keller inclusion

    B=C[P,Q] subset R=C[x,y], F=Frac(B), L=C(x,y),
    F subset E subset L, delta(E) subset E,
    delta=p(y)*partial_x, p in C[y] nonzero.

Scale a target coordinate so Jac_(x,y)(P,Q)=1; B and F are unchanged.
The source coordinates (x,y) are arbitrary. The already adopted plane
Rentschler input, explicitly retained in KX0, puts EVERY nonzero source
LND in this form after a polynomial source change. In adapted coordinates
the degree k=deg p is intrinsic: a different generator of ker(delta)
is affine in y and a complementary coordinate is a*x+b(y), a!=0.
These changes preserve k.

Use the complete retained KV0 classification INSIDE the actual source:

    C=R intersect E

is a finitely generated NORMAL smooth affine algebra, etale over B.
Its specified source map X=A2->V=Spec(C) is etale, with open image
whose complement has at most finitely many points. KV0 removes the
divisorial gaps before taking this spectrum; no arbitrary subfield
normalization is substituted by a plane.

Since both R and E are delta-stable, C is stable too. The restriction
is a nonzero LND: if it killed P,Q it would kill F and hence its finite
separable extension L, impossible. Let

    D=ker(delta|C)=C intersect C(y) subset C[y].

We will justify D=C[h(y)] for a nonconstant polynomial h, rather than
assume that the source invariant y descends.

## 2. Polynomial invariant and a rational slice

Here is the local-slice argument already written in LO0. Choose a in C
which is a last nonzero derivative of an element of C. Then a in D,
and some s in E satisfies delta(s)=1 and

    C[a^-1]=D[a^-1][s], E=Frac(D)(s).                          (2)

Successively subtracting the highest powers of s proves the ring
identity. It also proves Frac(D)=E^delta and trdeg(D)=1. In particular
D contains a nonconstant polynomial h0(y). Since C[y] is finite over
C[h0], D is a C[h0]-submodule of a finite module, hence is finitely
generated. D is normal in its fraction field: an element integral over
D is integral over normal C and, being delta-invariant, belongs to D.

Its fraction field is a subfield of C(y), so is rational by the already
retained Luroth input. Units of D are scalar because D subset C[y].
A smooth rational affine curve with scalar units has exactly one point
at infinity and is A1. Therefore

    D=C[h(y)], h a nonconstant source polynomial.                (3)

This is the elementary rational-curve/one-infinity argument used in
the retained surface-kernel proofs, not a new quotient classification.

Now delta(s)=1 gives

    s=x/p(y)+beta(y), beta in C(y),
    E=C(h(y),s), L=C(y,s).                                       (4)

The integration constant is rational because the kernel of partial_x
on C(x,y) is C(y). Since s is transcendental over C(y), (4) gives

    m=[L:E]=[C(y):C(h(y))]=deg h.                                (5)

No division by p INSIDE E has occurred. The coordinate expression
for s is an equality in L; p need not belong to E.

## 3. The target volume forces a polynomial derivative factorization

Put u=h(y). The ACTUAL target functions P,Q lie in E=C(u,s), so the
coefficient of dP wedge dQ relative to ds wedge du belongs to E.
On the other hand, from x=p(y)*(s-beta(y)),

    dx wedge dy=p(y)*ds wedge dy
               =(p(y)/h'(y))*ds wedge du.

As dP wedge dQ=dx wedge dy, the nonzero ratio p/h' belongs to E.
It is killed by delta, hence belongs to E^delta=C(h(y)). Thus

    p(y)=h'(y)*r(h(y)), r in C(u)^*.                             (6)

Crucially r has NO finite pole. If r had a pole of order l>=1 at c,
take any root y0 of h-c, of multiplicity e>=1. Characteristic zero
gives ord_(y0)(h')=e-1, so the right side of (6) would have order

    e-1-e*l <= -1,

contradicting polynomiality of p. Therefore r is a NONZERO POLYNOMIAL.
Taking degrees in (6) now gives the exact identity

    k=deg p=m-1+m*deg r,
    k+1=m*(deg r+1).                                             (7)

This proves (1), including m<=k+1. For k=0 it forces E=L.
No Euler, boundary-count, minimality or source-degree cutoff is used.

The retained ALL-coordinate primitive-element lemma in PF53 Section1
gives L=F(y). Since h in E and deg h=m, the tower degrees show

    E=F(h(y)).                                                   (8)

Indeed [L:F(h(y))]<=m, while it equals m*[E:F(h(y))]. This uses
the general lemma, not PF53's later special ramification hypotheses.

For E=F the stronger retained KX0 already forces inversion for EVERY
nonzero LND, including replicas. The present theorem instead controls
descent to arbitrary intermediate fields; it does not reopen that case.

## 4. The equivariant source image has no isolated point gaps

The invariant h defines an A1-fibration V->A1_u by (2)--(3). Reuse
the local structure package retained in MQ0/PF46: fibres on this smooth
affine surface are disjoint A1s with positive multiplicities. This package
keeps its original adopted/proof-sketch assurance boundary.

Every source line L_c={y=c} maps into the fibre u=h(c). It maps
nonconstantly because the source map is etale. Its image is a WHOLE
A1 fibre component, and the map to that component is an isomorphism:
a nonconstant polynomial A1->A1 is finite, and an unramified such map
has degree one. Source etaleness supplies the unramified condition.

KV0 allows only a finite missed set in V. The generic point of every
fibre component is therefore reached. Some L_c then maps onto the
ENTIRE component, so none of its points is missed. All fibres are
covered, proving that

    X=A2 -> V=Spec(R intersect E) is SURJECTIVE ETALE.           (9)

For a fibre at a, write its components as e_j*C_j and let n_j>=1
be the number of distinct roots c of h(c)=a whose source lines map
onto C_j. Etaleness and pullback of u-a give exactly

    ord_c(h-a)=e_j on each such root,
    sum_j n_j*e_j=m.                                             (10)

Thus ALL component multiplicities and fibre cardinalities are tied to
the one polynomial h. They are not arbitrary passports. In particular
a multiple root of h-a forces a nonreduced component in this actual
source image. We do not infer an entire intermediate etale model or
integrality of h over B from (9).

## 5. Monomial replicas give exactly cyclic deck fixed fields

Suppose p(y)=c*(y-a)^k, c!=0. All roots of h' in (6) are at a.
With m=deg h this forces

    h(y)=A*(y-a)^m+B, A!=0, m divides k+1.                       (11)

Put v=y-a. Multiplying the slice in (4) by r(h)=p/h', which DOES
belong to E, and scaling by a constant gives

    E=C(v^m, t), t=x/v^(m-1)+b(y), b in C(y),
    L=C(v,t).                                                    (12)

Hence L/E is GALOIS CYCLIC of order m, generated by v->zeta*v,
t fixed. Its automorphisms fix F and are actual generic rational decks.
For m>1 it is therefore the fixed field of a cyclic deck subgroup
of order m dividing gcd(d,k+1), where d=[L:F]. The induced x action
has multiplier zeta^(m-1)=zeta^-1 and a rational shift; it need not
preserve the source plane.

This is consistent with, and sharp at, MQ0: an actual cyclic subgroup
of order g has source-image kernel C[y^g], and its invariant replica
y^(g-1)*partial_x descends, with g=deg(p)+1. Thus one must NOT replace
the bare-translation theorem by a false prohibition on every replica's
intermediate-field descent.

If the actual generic deck group is TRIVIAL, no monomial replica in
ANY source-coordinate frame can preserve a proper intermediate field.
General nonmonomial replicas can have nongalois polynomial h and are
governed by (6)--(10), not silently made cyclic.

## 6. Finite source-derivative data recovers the whole function field

For EVERY bare source translation, in EVERY polynomial source frame,

    C(P,Q,partial_x P,partial_x Q)=C(x,y).                    (14a)

Here is a direct proof. Normalize Jac(P,Q)=1 and put
E_1=C(P,Q,P_x,Q_x). The Hamiltonian identity gives

    dy=-Q_x*dP+P_x*dQ in image(Omega_(E_1/C)).                (14b)

Let n=[L:E_1]. In characteristic zero every derivation of E_1 extends
uniquely to L and commutes with field trace. To verify this, extend it
to a normal closure: the derivative of the separable minimal polynomial
determines its value on each algebraic element, so the extension commutes
with every embedding; summing those embeddings proves trace compatibility.
Using a transcendence basis of E_1 gives the same compatibility for
Kahler differentials. Apply it to (14b). The differential trace is n times
the identity on forms from E_1, hence

    d(Tr_(L/E_1)(y))=n*dy.

Thus y-Tr(y)/n has zero differential in C(x,y). A rational function in
x,y with both partial derivatives zero is constant in characteristic
zero. Therefore y belongs to E_1.

Here is an elementary completion that needs no intermediate-surface
classification or coordinate-primitive theorem. Take any E_1-embedding
sigma of L into a finite normal closure. It fixes y and P_x,Q_x.
Extend delta=partial_x to that closure. If P_x!=0, chain rule gives

    delta(P(sigma(x),y))=P_x(sigma(x),y)*delta(sigma(x))
                       =P_x*delta(sigma(x)).

Since sigma(P)=P, the left side is P_x. Hence delta(sigma(x))=1.
If P_x=0 use Q_x, which is nonzero by the unit Jacobian. Thus
beta=sigma(x)-x is killed by delta. Differentiate its monic minimal
polynomial over C(y)(x). The resulting relation has lower degree,
so all its coefficient derivatives vanish. Those coefficients lie in
C(y), proving that beta is algebraic over C(y). The variable x stays
transcendental over C(y)(beta). Choose P or Q of positive x-degree r.
In its identity P(x+beta,y)=P(x,y), the x^(r-1) coefficient is
r*leading_x(P)*beta; it must vanish, and beta=0. Every embedding
therefore fixes x,y. Separability gives L=E_1, proving (14a).

The same calculation applies to a replica delta=p(y)*partial_x with
first velocities a=delta P,b=delta Q. If H'=p and
E_1=C(P,Q,a,b), then

    dH(y)=a*dQ-b*dP, H(y) belongs to E_1,
    L=E_1(y), [L:E_1]<=deg p+1.                             (14c)

The trace proof is identical. After adjoining y to E_1, division by
p(y) recovers P_x,Q_x, so (14a) gives E_1(y)=L. The polynomial
H(T)-H(y) over E_1 gives the bound. This is an upper bound; the stronger
divisibility statement (7) retains its differential-stability hypothesis.
No such hypothesis is needed for (14a)--(14c).

More generally, for ANY dominant generically finite rational map
G=(G_1,...,G_n) of C^n with nonzero Jacobian determinant, the field
generated by G and ALL its first partial derivatives is C(x_1,...,x_n).
Indeed the inverse Jacobian expresses every dx_i as a differential form
over that field, and the same differential-trace argument puts every
x_i in it. This assertion is field-theoretic and does not require a
constant Jacobian or source polynomiality. It supplies no integrality.

For ANY p!=0 define the finite source-derivative algebra

    J_delta=C[delta^j(P),delta^j(Q): 0<=j<=M],
    M=max(deg_x(P),deg_x(Q)), E_delta=Frac(J_delta).

Higher derivatives vanish because delta=p(y)*partial_x. This algebra
contains B and is delta-stable, so (7) applies:

    [L:E_delta] divides gcd(d,deg p+1).                          (13)

In particular for EVERY bare source translation in EVERY coordinate
frame,

    C(P,Q,partial_x P,partial_x Q,...,
      partial_x^M P,partial_x^M Q)=C(x,y).                       (14)

If gcd(d,deg p+1)=1, the same full-field conclusion holds for that
replica. If the generic deck group is trivial it holds for EVERY
monomial replica, regardless of k, by Section5.

Equation(14a) strengthens (14) to the FIRST derivatives in one FIXED
source direction. Higher derivatives supply no additional function-field
generation in the bare-translation case.
It is different from target Weyl/Nambu derivatives and from the full
jet-space map K times an affine identity retained in the jet audit.
It says the fixed-direction derivative map is BIRATIONAL onto its
image. It does not say its raw image is normal, that it is injective
at every point, or that the first two outputs have a polynomial inverse.

For a simple NON-KELLER control, the Ga-equivariant map

    (x,y)->(x,y^2-1,y*(y^2-1))

is birational with differential rank two everywhere, yet identifies
(x,1) and (x,-1) for every x. Its image algebra is nonnormal. No
constant-Jacobian pair is claimed inside it. This guards the exact
field-versus-global-algebra distinction in using (14).

## 7. What this changes in the main programme

We read the complete retained KV0, KX0, LO0, MQ0 and PF53 proofs.
The local slice, rational-curve and coordinate-primitive arguments
are reproduced at the needed scope. No new external theorem is adopted.
Primary literature on LND local slices and bracket width was screened;
it supplied no theorem making these field generators integral, making
the raw algebra normal, or supplying a Darboux pair on an ML1 surface.
No primary download/cache or predecessor mathematical replay is needed.

The new preregistered checker checks degree/ramification arithmetic,
slice/volume transcriptions, finite derivative controls, the nonnormal
birational control, source seals and scope. The unbounded arguments
above establish the conclusions; finite examples do not classify maps.

The field obstruction is now stronger than just failure of descent to
F: bare translations cannot descend to ANY proper intermediate field,
and already their two first target velocities generate the full field;
and replica degree restricts every other additive-equivariant factor.
Nevertheless the source-derivative algebra can still fail to supply
the ACTUAL full finite normalization or its trace-dual lattice. For
MP0's common source translation, J_delta lies inside A=B[Lambda];
full fraction field therefore does not make A finite over B or equal R.
Ramified omitted primes and marked-energy poles remain unresolved.
Full HC4/JC2 stays OPEN. The published quartic PLANE exclusion is
retained; independent external review remains unclaimed.
