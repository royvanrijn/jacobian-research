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
assumption on their generated field (Section 6). Section6.1 also proves that the transverse coordinate is integral over
this one-column jet algebra. In every generic source direction those
two first velocities already have the whole source as finite normalization;
the FULL first-jet algebra has that property in every dimension.
Section6.2 gives the exact RAW equality R=S[delta_v P,delta_v Q] for
generic v, with S the original full finite normalization of the target.
These are reconstructions over the jet algebra or S, not integrality over
the target ring B or inversion. Nonconstant replicas must not be divided
out in an unknown intermediate field.

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

### 6.1. Regular differentials give integral coordinate recovery

Here is a useful normal-domain lemma. Let A be a normal finitely generated
complex domain, with fraction field M, and f in M. If df belongs to the
image of Omega_(A/C) inside Omega_(M/C), then f belongs to A. At a
smooth general point of any prime divisor choose a transverse parameter
u. A pole f=u^-q*g, q>=1 and g a unit, would give

    df=u^(-q-1)*(-q*g*du+u*dg).                              (15)

Its transverse coefficient modulo u is nonzero in characteristic zero,
so df is not regular there. Thus f has no height-one pole. The normal
intersection-of-DVRs criterion gives f in A. This does not claim that
every regular differential is exact; it treats a specified rational
primitive, and is false in positive characteristic.

For the fixed plane first-column algebra put

    J_x=C[P,Q,P_x,Q_x], A_x=normalization of J_x in L.

Section6 gives Frac(J_x)=L. Normality of R and J_x subset R give A_x
subset R. Equation(14b) is a REGULAR differential expression over A_x,
so the lemma proves

    y belongs to A_x; equivalently y is INTEGRAL over J_x.     (16)

This is stronger than its earlier field membership. It does not make
the first velocities integral over B. It also does not put the other
coordinate x into A_x.

If P or Q has a nonzero CONSTANT leading coefficient as a polynomial
in x over C[y], then x is integral over B[y]. Equation(16) and normality
put both x,y in A_x. Hence

    A_x=R and R is FINITE over J_x                           (16a)

under that source-finiteness condition. This condition holds in EVERY
generic constant source direction: if P has ordinary degree D and
highest homogeneous part P_D, a direction v with P_D(v)!=0 makes the
x^D coefficient exactly that nonzero constant, in any complementary
linear source coordinate. The excluded projective directions are at
most D. A translation or an arbitrary preselected polynomial source
frame does not alter this argument within that frame. Genericity is
only sufficient; some excluded directions may satisfy monicity too.

Thus for every generic direction the FOUR-function first-column map
(P,Q,partial_v P,partial_v Q) is finite and birational onto its closed
image, and its normalization is the ENTIRE affine-plane source. Finiteness
here uses the standard normalization theorem stated below. This result
requires only the TWO FIRST velocities, not all higher source jets.
The field statement(14a) keeps its stronger ALL-direction scope.

There is a global geometric refinement. The source map X->Spec(A_x)
is birational and quasi-finite, since its composite with the target map
is the original quasi-finite Keller map. Zariski--Main therefore makes
it an open immersion. Over C(y), one of P,Q has positive x-degree.
Its leading coefficient is a unit in that field, so x is integral over
C(y)[P,Q] and hence over (A_x)_(C[y] minus zero). Localization preserves
normality, and x is in L; therefore

    (A_x) tensor_(C[y]) C(y)=C(y)[x].                         (17)

Consequently some nonzero h(y) satisfies

    A_x[h^-1]=C[x,y,h^-1].                                   (18)

Thus the normalized first-column image is an affine A1-fibration over
the ACTUAL transverse source coordinate y; all its added boundary lies
in finitely many fibres of y. This statement keeps the marked source
open and does not turn that image into a finite normalization over B.

A source replica extends to this entire normalized image. Choose finite
C-algebra generators g_i of A_x. By(18), their derivatives partial_x g_i
lie in A_x[h^-1]. A common power h^N clears these finitely many
denominators, so

    D=h(y)^N*partial_x preserves A_x.

It is locally nilpotent because A_x subset C[x,y] and it kills h(y).
Its invariant ring is exactly C[y]. This supplies an actual additive
action on Spec(A_x), preserving the whole source open. The power N can
depend on the map and direction. Bare partial_x need not preserve A_x,
and neither this action nor its invariant coordinate is asserted to
descend to B or F.

For the FULL first-jet algebra the same finite-normalization conclusion
holds in EVERY source dimension n, without a direction restriction. Let G be any complex polynomial Keller map and set

    J_full=C[G_i,partial_j G_i : 1<=i,j<=n],
    R=C[x_1,...,x_n], L=Frac(R).

The general field result in Section6 gives Frac(J_full)=L. The inverse
Jacobian is the adjugate divided by its nonzero CONSTANT determinant,
so its entries belong to J_full. Hence every dx_i is a regular
J_full-linear combination of dG_1,...,dG_n. Apply(15) in the normalization
A_full: it contains every x_i. Conversely A_full subset R because R is
normal and contains J_full. Therefore

    A_full=R, R is FINITE over J_full.                        (19)

The finiteness is the finite-normalization theorem for varieties over
a field: [Stacks, Lemma33.27.1](https://stacks.math.columbia.edu/tag/0BXQ)
applies to the finitely generated complex domain J_full. Its exact
statement and proof were inspected; normality of a raw jet image was
not assumed. Equivalently, the full first-jet map is finite and
birational onto its closed image, and the original affine source is its
normalization. This does not force that raw image to be normal or its
projection to the first n target outputs to be finite.

The constant-Jacobian restriction matters. For the NON-KELLER polynomial
map G=(x^2*y^2,y), its full first-jet algebra is
C[y,x*y^2,x^2*y]. Along x=t^-1,y=t^2, all its generators have nonnegative
t-order but x has a pole. The t-adic valuation therefore proves x is
not integral over that algebra, although its fraction field is C(x,y).
Its Jacobian 2*x*y^2 is nonconstant. This is a failure of a more general
polynomial-map assertion, not a plane Keller counterexample.

No new conclusion about finite landing over B follows from(16)--(19),
including the generic two-first-velocity conclusion(16a).
In particular (19) is compatible with the retained noninvertible Keller
maps in dimension at least three. The new integral extension has base
J_full, while the unresolved finite-cover and trace questions have base
B. Those two bases cannot be interchanged.

### 6.2. Two generic velocities recover the raw source over the finite normalization

Let S be the FULL finite normalization of B in L. For any direction v
satisfying the sufficient monicity condition in(16a), put

    a=delta_v P, b=delta_v Q.

Then there is an exact equality of ACTUAL embedded source algebras,

    R=S[a,b]=B[Lambda,a,b].                                  (20)

The second equality uses the plane canonical trace-dual embedding
S subset Lambda subset R, already proved by the finite normalization,
trace and unramified source valuation arguments in KU0/1. No differential
IC-socle or semisimple source-layer theorem is needed. Any finite
B-generating family of Lambda suffices in(20).

Here is the general open-immersion lemma proving the assertion. Suppose
A subset C subset R are domains and Spec(R)->Spec(A) is an open
immersion. If R is finite over C, then C=R. Indeed lying-over makes
Spec(R)->Spec(C) surjective. Every point of Spec(C) thus maps into the
specified open Spec(R) in Spec(A). The map Spec(C)->Spec(A) consequently
factors through that open, giving a ring homomorphism R->C. Its composite
with the given inclusion C->R is the identity: the corresponding source
map becomes the identity after composition with the open immersion,
which is a monomorphism. Hence that inclusion is surjective, proving
C=R. No normality of C was assumed, and the factorization treats all
points, not only prime divisors or the smooth locus.

Apply this with A=S and C=S[a,b]. Equation(16a) says R is finite over
J_v=B[a,b]. A finite J_v-module generating set of R also generates it
as a C-module because J_v subset C subset R. Thus R is finite over C.
KN0 supplies the canonical OPEN immersion X=Spec(R)->N=Spec(S), so
the lemma gives the first equality in(20). Since

    S[a,b] subset B[Lambda,a,b] subset R,

the second equality follows immediately. Equivalently, one can apply
the lemma to the actual canonical etale completion U=Spec(B[Lambda])
using KU0/1. The finite-normalization argument already suffices.

There is also an unrestricted-dimensional version. For ANY complex
polynomial Keller map G in n variables, with full finite target
normalization S, equation(19) and the same open immersion give

    C[x_1,...,x_n]=S[partial_j G_i : 1<=i,j<=n].               (21)

In particular, for an actual four-variable unit-Hessian potential h,
apply (21) directly to G=grad(h). Symmetry leaves at most TEN generators:

    C[x_1,x_2,x_3,x_4]=S[h_(ij) : 1<=i<=j<=4].              (22)

Here S is the full finite normalization for the gradient map itself.
No Tang classification is needed for (22). Integrality of these Hessian
entries over the gradient target remains unproved.

Thus raw source reconstruction over the ORIGINAL full finite normalization
requires only first source derivatives, and in the plane only two first
velocities in a generic direction. These statements keep S itself in
the generators; they do not identify the raw first-jet image with its
normalization or supply a new polynomial inverse.

The open-immersion condition in the lemma is essential. The finite
birational normalization

    C[u,v]/(v^2-u^2(u+1)) subset C[y],
    u=y^2-1, v=y(y^2-1),

identifies y=1 and y=-1. It is finite because y^2=u+1, yet the two
rings differ. Taking a product with A1 gives the nonnormal source-image
control in Section6. Such a normalization map is not an open immersion,
so it cannot be substituted for KN0 in the proof of(20).

The remaining landing condition is now particularly exact. For a generic
v as above,

    a,b in S  iff R=S  iff K has a polynomial inverse.        (22)

The forward direction follows from(20); R=S makes K finite etale, and
characteristic-zero affine-space covering triviality gives the inverse.
The converse is immediate. Similarly a,b in B[Lambda] is equivalent
to R=B[Lambda], namely absence of omitted UNRAMIFIED primes; it can
leave ramified omitted primes outside that completion. These are two
different landing conditions. Generic first-velocity integrality over
B, full trace-dual perfection and JC2 remain unproved.

## 7. What this changes in the main programme

We read the complete retained KV0, KX0, LO0, MQ0 and PF53 proofs.
The local slice, rational-curve and coordinate-primitive arguments
are reproduced at the needed scope. Section6.1 additionally inspects the standard finite-normalization input
for varieties. Primary literature on LND local slices and bracket width was screened;
it supplied no theorem making these field generators integral, making
the raw algebra normal, or supplying a Darboux pair on an ML1 surface.
The inspected finite-normalization statement and current application
limits are retained in
[sources/first-jet-normalization-context.json](sources/first-jet-normalization-context.json).
No unrelated predecessor mathematical checker was replayed.

The retained checker checks degree/ramification arithmetic,
slice/volume transcriptions, finite derivative controls, the nonnormal
birational control, source seals and scope. The unbounded arguments
above establish the conclusions; finite examples do not classify maps.
The open-immersion corollary in Section6.2 was independently checked
internally; it uses finite lying-over and open factorization, not a
finite symbolic census or the unverified differential KQ0 package.

The field obstruction is stronger than just failure of descent to
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
