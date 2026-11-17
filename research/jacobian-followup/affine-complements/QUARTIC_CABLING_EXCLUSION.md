# The marked Wright surface has no finite quartic Keller map

14 November 2026. PF2RULED35 closes the finite degree-four map problem
on the actual marked Wright surface, at every pole height.

The proof uses the full marked end pair, the residual-cover longitude,
and a finite classification of degree-four cabling behavior. One cabling
step can survive, but it reduces the number of cabling stages while
preserving the same marked lens-space condition. The other possible
transition goes to alternating monodromy and is excluded separately.
No coefficient-degree or tangency-order cutoff is used.

This excludes the compatible ruled quartic realization, not an arbitrary
quartic Keller normalization. Existence of a compatible completion for an
arbitrary actual normalization is not proved. Higher finite field degrees
on the Wright ring are also not excluded by this quartic argument.

## 1. The actual end pair and the longitude

Keep the actual ring and finite-map hypotheses from
[PF2WRIGHT1](WRIGHT_MARKED_MODEL.md), [PF2RULED33](RESIDUAL_DOUBLE_COVER.md)
and [PF2RULED34](TORUS_INFINITY_OBSTRUCTION.md):

    A=C[Y,XY,X^2Y,X^3Y+X] inside C[X,Y],
    U=A2_(X,Y)=N minus R, N=Spec A,
    pi=(P,Q):N->A2 finite, [C(N):C(P,Q)]=4,
    J_(X,Y)(P,Q)=c != 0.

The whole branch curve is C=pi(R), with R=A1, e_R=2, f_R=1,
finite birational normalization gamma:A1->C, and global monodromy S4.

### The end knot is a core, not just a nonzero homology class

The actual two-chart transition is

    u=1/X, v=X^3 Y+X;   X=1/u, Y=u^3 v-u^2.

These charts define an affine-line torsor over P1, with linear
transition v=X^3Y, namely the line bundle O(-3). They give the entire
marked N. The second plane chart is the complement of the retained
line X=0; on its overlap with U the displayed inverse is exact.
Its u=0 fiber is precisely the actual R, with no extra points deleted.

An explicit smooth section is

    Y_s=-X*conjugate(X)^3/(1+|X|^6),
    v_s=X/(1+|X|^6)
       =u^2*conjugate(u)^3/(1+|u|^6).

Translation by that section identifies the torsor smoothly with O(-3),
preserving the entire fiber R over infinity. The fiber metric
|Y-Y_s|^2*(1+|X|^2)^3 is global and proper over the compact base.
Its circle-bundle end pair is

    M=L(3,1) up to the irrelevant orientation convention,
    K_R=a circle fiber,  M minus int nu(K_R)=S1 x D2.       (1)

In particular H1(M)=Z/3 and K_R generates it. This also recovers
the earlier 4mu=mu calculation from E.B=4 and B^2=3.

For the target radial cross-section M it is enough to compare the
stable groups, rather than assume that any two four-dimensional
collars are identical. Cofinal product neighborhoods of the pair
give pi1(M)=C3, pi1(M minus K_R)=Z and the same nonzero marked class.
These are neighborhoods in N with R removed, not the usual end
neighborhoods of U alone.

Elliptization, in Perelman Section1.4, recognizes a closed oriented
M of group C3 as a cyclic spherical space form, hence a lens space.
It is irreducible. Since K_R is homologically essential, every
sphere in its exterior bounds a ball disjoint from K_R; the exterior
is irreducible too. Its torus boundary cannot inject into its cyclic
fundamental group. The loop theorem and compression then identify
that exterior with a solid torus. Thus (1) holds for the actual
radial cross-section as well. This recognition input is imported
and pinned explicitly.

Its restriction is a simple four-sheeted branched cover f:M->S3
of the boundary knot K, with ramification knot K_R. The
infinity-group epimorphism retains full S4 monodromy.

### Its preferred longitude has identity monodromy

PF2RULED33 gives the complete residual quadratic cover

    w^2=product_(m_i odd) (t-a_i)(t-b_i).

Its squarefree degree is even; the split case has two affine-line
components instead. In either case its two points at infinity are
unramified over the boundary normalization. Thus transport around the
end of C fixes both residual sheets.

Let mu,lambda be the meridian and preferred zero-linking longitude
of K. The meridian maps to a transposition tau. The longitude
commutes with tau and fixes the two residual letters. Hence its image
is either 1 or tau. A preferred knot longitude is in the commutator
subgroup, whose image lies in A4. Therefore

    rho(mu)=tau=(12),  rho(lambda)=1.                       (2)

Changing the push-off framing only multiplies by a power of tau,
which fixes the residual letters; it does not affect this deduction.
In particular a normal-bundle framing has not been mistaken for the
preferred longitude.

### What must now be excluded

We prove that no iterated torus knot admits a degree-four cover with

    (S) monodromy S4, meridian (12), longitude 1,
        M a lens space with H1(M)=Z/3,
        ramification-knot exterior a solid torus.           (3)

Neumann--Norbury Section4 describes every one-place-at-infinity
polynomial curve knot by a finite iterated-torus cabling sequence.
This standard fact is imported with its actual one-place hypothesis.
Neither a positive local-algebraic-knot inequality nor transversality
to the projective line at infinity is assumed.

## 2. Topological and group facts used in the cabling reduction

Write K as a cable of a nontrivial companion J, with winding p>=2
and signed coprime second parameter q. The regular Seifert fiber of
the cable space has peripheral slopes

    h = p*lambda_J + q*mu_J
      = lambda_K + p*q*mu_K,                              (4)

with harmless simultaneous orientation changes allowed. The inner
pq framing, the outer section/fiber bases and the change of cabling
notation are recorded in the pinned Ito source. Our p is the winding
number; Ito calls that number q. Only these slope facts are imported.

Let V=nu(J). In the branched preimage of V, there is one component
B containing the ramification knot. Let m be its covering degree.
All other components cover V unbranched and are solid tori.
Let H be the image of the companion group, and write

    sigma=rho(mu_J), l=rho(lambda_J).

A knot group is normally generated by its meridian, so sigma normally
generates H. Also l belongs to H' and commutes with sigma.

The following literal subgroup facts are elementary and checked using
all thirty subgroups of S4:

* a knot-group image normally generated by one transposition is S_j
  on one orbit of j letters, j=2,3,4, fixing the other letters;
* a subgroup normally generated by a three-cycle is C3 on three
  letters or A4 on four letters;
* a subgroup normally generated by a double transposition is its C2;
* in A4 a three-cycle has centralizer C3 and A4'=V4. Its longitude
  therefore always has identity image.

If the fiber image is identity, let alpha be the exceptional-value
monodromy on B and let t be the branch meridian on that orbit.
The pattern image is generated by alpha,t. Up to inversions that
do not change any cycle types below,

    alpha^p=1, delta=alpha*t,
    sigma=delta^p, l=delta^(-q).                           (5)

These are peripheral statements in the fiber/section basis; in
particular the p-th power cannot be replaced by delta itself.
The base Euler characteristic of B is

    chi(base B)=c(alpha)+c(t)-m,                           (6)

and its outer boundary count is c(delta). An alpha cycle of length a
gives a cone of order p/a. Order-one entries are omitted below.

If h has nontrivial image, it is central in the pattern image. A
central permutation on a transitive orbit has no fixed point unless
it is the identity. This restricts the branched orbit to the support
of the branch permutation in every such case used below.

### Cutting, compression, and filling

A nontrivial knot exterior has incompressible boundary; every finite
unbranched cover has the same property. An orientable Seifert manifold
with nonempty boundary also has incompressible boundary unless it is
a solid torus. Annular bases, including T2 x I, are retained.

When all sides of the cutting tori are boundary-incompressible,
the graph-of-groups injections give Z2 in the fundamental group.
That is impossible in a lens space or a solid torus. A cycle in the
dual graph also gives a free Z quotient of first homology.

We repeatedly use this compression consequence. Suppose a torus
separates off a piece containing a homologically essential knot in an
irreducible lens space, and a compressing disk exists on the opposite
side. That opposite side is a solid torus. Compressing the torus
gives a sphere; the alternative ball containing the knot would make
its homology class zero. This argument excludes the ball alternative
rather than silently discarding it. The loop theorem and compressible
torus facts are imported from Hatcher.

Filling a boundary of a Seifert piece along a nonfiber slope adds a
cone of order its distance from the fiber. Filling along the fiber
kills that fiber. If the marked knot is regular, it is then null in
homology. For a disk with at least two genuine cones, fiber filling
has a noncyclic free-product quotient. These elementary presentation
facts exclude the degenerating filling separately.

## 3. Exclude the alternating intermediate type

First exclude the auxiliary type

    (A) transitive A4 monodromy on four letters,
        meridian a three-cycle, longitude 1,
        the same marked lens-space/core conditions as in (3).

The longitude condition follows from the subgroup facts above.
The branch fiber has one index-three knot and one unramified sheet.

For a torus knot T(a,b), its central fiber dies in A4. Finite first
homology forces a sphere base. Cycle counting gives c(alpha)+c(beta)=4.
Both nonidentity even permutations have two cycles. Generating A4 and
coprimality of a,b force one three-cycle and one double transposition.
Write a=3u, b=2v, exchanging them if necessary. The lifted cone orders
are u,3u,v,v, including order-one entries. Hence

    |H1(M)|=(a^2*b^2/12)*4/(a*b)=a*b/3=2u*v,               (7)

which is even and cannot be three. The unknot is excluded by its
cyclic complement group.

Now let the target knot be a genuine cable. By (4) the pattern fiber
image is t^(pq), with t a three-cycle.

If 3 does not divide pq, the branched orbit has exactly three
letters and the pattern group there is C3. The fourth letter is
unbranched. The companion image must be A4, since otherwise the
whole action would fix that fourth letter. Both companion longitude
and its peripheral action are therefore as above.

The branched piece B has a disk base with cones p and 3. Its
ramification knot is the order-three exceptional fiber. The one
unbranched cap leaves a companion piece Y with one boundary.
B is not a solid torus. Incompressibility and the compression
consequence force Y to be a solid torus. Fiber filling would give
the noncyclic quotient C_p * C_3. A nonfiber lens-space filling must
add an order-one cone, since three genuine cone orders have a
noncyclic triangle-group quotient. Thus the closed base has orders
p,3. Since gcd(p,3)=1, its orbifold abelianization is zero and the
regular fiber generates H1(M). But that regular fiber is three
times the ramification knot, hence zero in Z/3. Contradiction.

Suppose instead that 3 divides pq, so h maps to 1.

If m=3, alpha is 1 or a three-cycle. If 3 divides p, then
sigma=delta^p=1, forcing the companion image to be trivial; the
whole group could not be A4. Otherwise alpha=1, 3 divides q, and
B has a disk base with three cones p,p,p. The companion image is A4
and the remaining unbranched cap again leaves a one-boundary Y.
Compression forces Y to be a solid torus. Fiber filling kills the
marked regular fiber; any nonfiber filling still has at least
three genuine cone orders and cannot give a lens space.

If m=4, the pattern image is A4, alpha has two cycles, and (6)
gives an annulus base with two outer boundary tori. There are no
unbranched caps. B and all companion covers have incompressible
boundary, contradicting finite fundamental group.

This exhausts type (A) for every nontrivial cable, independently of
the companion's cabling depth. Type (A) is therefore excluded on
all iterated torus knots.

## 4. The simple symmetric type: every outer cable case

Return to type (S), with branch transposition t. Equation(4) gives
rho(h)=t^(pq). Thus m is 2,3 or 4; for m>=3 the fiber image must
be 1. In the latter case (6) reduces to chi(base B)=c(alpha)-1.

### Four-letter branched orbit

Here all four sheets lie in B, so there are no unbranched caps.

If alpha is a four-cycle, B has an annular base. Its boundary and
all companion-cover boundaries are incompressible, impossible.

Otherwise alpha has two cycles joined by t. Their lengths are
either 1,3 or 2,2. B has a disk base with the corresponding orders
p,p/3 or p/2,p/2. If two genuine cones remain, B and the companion
cover are boundary-incompressible, impossible. Thus the only
possibilities are:

* p=3 and cycles 1,3. B is a solid torus with one cone of order
  three. Its marked regular fiber winds three times around its
  core and is zero in H1(M)=Z/3, impossible.
* p=2 and cycles 2,2. The companion meridian is delta^2, a double
  transposition. Its normal-generation subgroup cannot act
  transitively on four letters, while the one boundary of B and
  connectedness of M require a connected companion cover.
  Equivalently (5) makes its longitude a four-cycle, incompatible
  with that companion image.

Thus m=4 is excluded.

### Three-letter branched orbit with alpha a transposition

The pattern group is S3. Here p is even; delta is a three-cycle.
If 3 divides p the companion meridian is identity, impossible
for the full S4 action. Otherwise its image is a three-cycle.
The companion image must be A4 on all four letters; the C3
alternative fixes the fourth letter globally. Its longitude is 1.

There is one unbranched cap. After it is attached, the companion
piece Y has one boundary. B has disk base with orders p/2,p.

For p>=4, compression forces Y to be a solid torus. Fiber filling
makes the marked regular fiber null. Nonfiber filling of B minus
that regular fiber gives a disk with at least the two genuine
cones p/2,p, so its fundamental group is not cyclic. This
contradicts the solid-torus marked-knot exterior.

For p=2, B itself is a solid torus with one cone of order two.
Its ramification knot is a regular fiber of winding two in B.
The exterior B minus that knot is a cable space. If Y had
incompressible boundary, gluing it would inject Z2 into the
marked-knot exterior. Thus Y compresses and is a solid torus.

The companion boundary lattice on this three-sheeted component is

    m0=3*mu_J, l0=lambda_J.

The branched map B->V sends l0 to a generator of H1(V), so its
map H1(B)=Z -> H1(V)=Z is an isomorphism. It sends m0 to zero.
The primitive slope m0 is therefore exactly the meridian of B,
not a different surgery slope.

Replacing B by the usual order-three branched filling changes
neither the closed manifold nor its core knot. It gives a cover
of J of type (A): A4 monodromy, index-three meridian, the same
lens space, and ramification core. Type (A) is already excluded.
This closes the p=2 transition.

### Three-letter branched orbit with alpha a three-cycle

Now 3 divides p and delta is a transposition. Even p gives
identity companion meridian and is impossible. Thus p is odd,
the companion meridian is a transposition, and its longitude is 1.
For the latter assertion, in the only possible larger group S4
a nontrivial longitude would be a double transposition; since p
is odd, (4) would leave a residual transposition in rho(h)=1.
The smaller S2 and S3 images have trivial longitude directly.

The companion image can be S4 on four letters, or S3 containing
the pair and the letter outside B. Smaller support would fix
the fourth letter globally.

In the S4 case, the companion cover is connected, has three
boundary components, and receives the one unbranched cap.
The remaining two boundaries attach to the two boundaries of B.
The dual graph has a cycle, so first homology is infinite.

In the S3 case, one sheet of the companion is a copy of the
original exterior X_J. The other three sheets, after their
unramified cap, give a one-boundary piece Y. The branched piece
B has annulus base with one cone of order p/3. Its two boundaries
attach to X_J and Y.

If Y were boundary-incompressible, all the gluing tori would
inject. Thus Y compresses and is a solid torus, using the
homologically essential marked knot on the opposite side.
Fiber filling would make that regular knot null. Any nonfiber
filling of Y leaves, after removing the marked knot, an annular
Seifert piece (possibly T2 x I) with incompressible boundaries.
Its gluing to X_J injects Z2 into the marked-knot exterior.
That exterior is a solid torus, a contradiction.

Thus m=3 is excluded in both alpha cases.

### Two-letter branched orbit: the one possible removal step

In this orbit all branch meridians are the same transposition.
A meridian disk of V meets the cable p times, so

    sigma=t^p.

Even p would make the companion image trivial and leave the
whole action nontransitive. Hence p is odd. A companion image
normally generated by t is S_j on one orbit. Its longitude is
1 except possibly a double transposition in S4. Equation(4),
with p odd and rho(h)=t^(pq), rules that alternative out.
Thus the companion image is S4, its longitude is 1, and the
other two pattern sheets are separate unbranched caps.

After those caps, the companion piece Y has one boundary.
For even q the fiber image is 1; B has disk base with two
cones p,p and the ramification knot is regular. Compression
forces Y to be a solid torus. Fiber filling kills the marked
class; nonfiber filling cannot make the knot exterior a solid
torus in the presence of those two genuine cones. So q is odd.

For odd q the fiber degree is two. B has disk base with cones
p and 2, and the ramification knot is the order-two exceptional
fiber. Again Y is a solid torus. Fiber filling cannot occur:
twice the marked knot is nonzero in Z/3. Nonfiber filling must
have distance one from the regular fiber, since removing the
order-two fiber leaves a solid torus only if no further genuine
cone is added.

Fill the missing companion ramification torus in the usual way.
On its boundary the primitive cover lattice is

    m0=2*mu_J, l0=lambda_J.

The two cones of B are coprime, so H1(B)=Z, with no torsion.
The map B->V sends l0 to the generator of H1(V); consequently
it induces an isomorphism of these groups and sends m0 to zero.
Its boundary homology map is the same as that of the usual
ramification solid torus W_J, whose meridian is m0.

Therefore the replacement

    M_K=Y union B  --->  M_J=Y union W_J

preserves first homology Z/3. Both Y and W_J are solid tori,
so M_J is again a lens space, and its ramification knot is its
Heegaard core. It is again a simple S4 cover, with longitude 1.

This is a statement about the actual gluing lattice. It does
not identify the two closed manifolds merely because their
homology groups agree. A possibly different lens space is
allowed by the induction.

One can also track the marked class. The source regular fiber
maps twice around the target fiber, so in this basis it is

    h_B= +/- q*m0 + 2p*l0.

Twice the ramification class equals h_B. Since H1(B) is free,

    [K_R]= +/- p*[l0] in H1(B).

Thus 3 does not divide p; the removed winding is coprime to six.
The companion ramification core is still homologically essential.

We have reduced a hypothetical type (S) cover of K to a type
(S) cover of its strictly simpler companion J, retaining every
condition in (3).

## 5. Finite descent closes the marked quartic problem

Every one-place-at-infinity knot here has finitely many cabling
stages. The only surviving step in Section4 removes one stage.
A transition to type (A) is impossible by Section3. Repeated
removal must reach a torus knot or the unknot, both excluded
by PF2RULED34. This proves nonexistence of type (S).

Sections1 and2 show that every actual finite quartic Keller
map on the marked Wright surface would supply type (S).
Therefore no such P,Q exist, at any pole height.

Together with the earlier period and recognition reductions,
there is no compatible ruled quartic realization in the
specified boundary class. For a hypothetical arbitrary quartic
normalization, this obstructs a compatible completion; it does
not prove that every normalization has one or exclude either
full quartic packet without it.

## 6. Evidence and precise assurance

The primary [source manifest](sources/quartic-cabling/manifest.json)
pins Neumann--Norbury's iterated-cabling description, Hatcher's
compression and incompressibility results, Ito's cable slopes and
Perelman's finite-group end recognition.
The earlier pinned Scott and Leidy--Maxim inputs remain explicit
dependencies. Their full proofs are imported, not independently
re-proved here.

The [checker](verify_quartic_cabling.py) and
[receipt](quartic_cabling_receipt.json) verify the two-chart
transition, finite S4 subgroup facts, all group-theoretic
fiber-trivial residue cases modulo12, nontrivial-fiber cases,
integral homology and the descent's strict depth decrease.
The modulo12 coverage is exhaustive only for the permutation
and peripheral residue layer; it is not a formal verification
of the topological gluing proof. Unbounded winding parameters,
slopes and cabling depth are covered by the written arguments.

This is a written proof supported by exact controls. It has
not received independent whole-proof implementation, formal
verification or external review. No novelty priority is claimed.

The previous coefficient fronts remain retained algebraic
experiments, not solved unit-ideal systems. Their realization
question is now excluded by this global theorem; no assertion
is made that every truncated necessary system becomes empty.
The unrestricted higher-field-degree Wright problem and the
arbitrary-normalization quartic problem remain open.
