# Torus-curvette ends: coprime poles and prime boundary orders

14 November 2026. PF2RULED42 excludes the whole coprime regular-coordinate
sector in field degree four, uniformly in all Laurent coefficients.
Section7 extends the exclusion to all prime boundary orders and,
more generally, to every reference branch with one characteristic pair.
Its geometric end calculation and explicit monomial controls are
independent of the field degree.

Let N be an existing smooth affine surface with marked source
U=A2_(x,y), sole complement R=A1 and div(dx wedge dy)=R.
Suppose y is a regular source coordinate, shifted to boundary value zero,

    v_R(y)=r>=2, v_R(x)=-q, gcd(q,r)=1.

PF2RULED40 gives q>=r-1. Define

    m=q-r+2>=1,  n=qr+m.                                (1)

For every full Laurent orbit of that theorem, the marked end is

    n-surgery on T(q,r), with the omitted knot its surgery dual;
    H1=Z/n, omitted knot primitive;
    closed Seifert orders q,r,m, |e|=1/(qr)+1/m;
    omitted-knot exterior X_(T(q,r)).                    (2)

Both the closed source and the marked exterior are atoroidal.
All lower Laurent coefficients are retained: they move the free
blowup centers, but not the graph, puncture or filling slopes.

A stronger quartic theorem excludes simple full-S4 covers of finite
iterated torus knots when the source is irreducible and atoroidal,
has cyclic odd H1 of order n>3, and has primitive ramification with
atoroidal exterior and identity preferred longitude. A lens-space
hypothesis is no longer needed. The odd order-three exception has
an actual topological control and is not silently included.

Thus no coprime model above admits a finite quartic Keller extension.
The sole order-three pair (q,r)=(1,2) was already excluded on Wright.
In all finite degrees, the marked end still excludes target torus/unknot
infinity and prime boundary-coordinate degrees.

For any remaining actual quartic normalization with exactly one omitted
prime of index two and a regular source coordinate, the necessary region is

    r>=2, q>r, gcd(q,r)>1;
    q>=r+ell(r), ell(r)=smallest prime divisor of r.       (3)

In particular q>=r+2; if r is prime, q>=2r. Section7 further requires composite r and at least two characteristic
pairs in the remaining reference branch. The remaining noncoprime sector,
absence of a regular coordinate, extra primes, higher field degrees
and full objective remain open. These are written all-parameter proofs
with exact supporting controls, not formal or external verification.

## 1. The entire actual Laurent chart determines the blowup contact

Use the full root charts of PF2RULED40, not a generic valuation model:

    y=z^r, x=z^(-q)*(H(z)+z^m*t),
    H(z)=sum_(j=r-1)^q c_j*z^(q-j), c_q!=0.             (4)

Choose integers A,B with -q*A+r*B=1 and set toric coordinates

    s=x^A*y^B=z*(H+z^m*t)^A,
    a=x^r*y^q=(H+z^m*t)^r.
    x=a^B/s^q, y=s^r/a^A.                              (5)

These are regular coordinates near the open torus of the exceptional
divisor with ray (q,r): s has valuation one and a has nonzero residue
c_q^r. Negative powers of a or H cause no problem there, since both
are units along the whole omitted line.

For t=0, s_0=z*H(z)^A has nonzero linear coefficient. Invert it
formally and put psi(s)=H(z(s))^r. Let phi(s) be its truncation
through degree m-1. The m free blowups along the smooth reference
germ a=psi(s) have last coordinate

    T=(a-phi(s))/s^m.                                  (6)

The contact order is exactly m, with the actual whole R present.
Indeed changing t gives

    s-s_0=O(z^(m+1)),
    a-H(z)^r=r*c_q^(r-1)*z^m*t+O(z^(m+1)).

Consequently, on R the coordinate in (6) is

    T=psi_m+r*c_q^(r-1-A*m)*t.                          (7)

Its slope is nonzero, so every finite t on each root omitted line
is included. The root components are permuted freely, and the
expressions s,a,T are invariant; they give the same whole omitted
A1 after the quotient.

The map to the last blowup chart is regular and etale along that
entire A1, by (7) and the simple transverse parameter s. It is
birational and hence an open immersion near it. On the whole source
plane it is the identity. Thus the deleted-boundary compactification
below recovers the actual marked affine N, rather than just its
formal neighborhood. No affine converse for arbitrary Laurent data
is being asserted.

The volume confirms the number of free blowups without using a
bounded coefficient calculation:

    dx wedge dy=-s^(r-q-1)*a^(B-A-1) ds wedge da
               =-s*a^(B-A-1) ds wedge dT.               (8)

The toric divisor has volume order r-q-1=1-m; after m free
blowups the last divisor has order one, exactly as required.

## 2. The full punctured graph and surgery coefficient

Start with P1_x times P1_y, with original deleted curves
D_x={x=infinity} and D_y={y=infinity}. Refine the corner
(x=infinity,y=0) until the ray (q,r) appears. The usual adjacent-ray
sum procedure terminates by the Euclidean algorithm because gcd(q,r)=1.
At every step the target coordinates in the current unimodular cone
are positive coprime integers and their sum decreases after the choice
of the next cone. The target divisor E is the last toric exceptional
and has self-intersection minus one.

The exceptional curves form a linear chain. Blow up m times at the
successive free reference points in Section1. The first step decreases
E's self-intersection by one; the later steps form a chain with
m-1 deleted negative-two curves and final retained Rbar(-1).
Rbar is transverse to the last deleted curve, or to E if m=1.
Keep this puncture throughout.

D_y is a zero leaf and D_x has just one remaining attached branch.
The retained Neumann Proposition2.1 R6 removes this pair, with
genus zero, one attachment, no deleted-fiber boundary on the removed
pair and no extra summand. The remaining graph is the actual corner
exceptional chain, with E's weight decreased by one and the
negative-two free chain attached at E.

Here is a uniform way to identify its arms. Let C_0 be the positive
matrix of the corner exceptional chain before the free blowups.
It is unimodular: it is the exceptional matrix of point blowups of
a smooth surface. If its rays are (a_i,b_i), the principal divisors
of w=1/x and y give

    C_0*(a_i)_i=e_left, C_0*(b_i)_i=e_right.

For a tridiagonal chain with off-diagonal minus one, its inverse
entry between a vertex and an endpoint is the determinant of the
opposite arm. At E the ray is (q,r), so the two arm determinants
are q,r and (C_0^-1)_(E,E)=qr. This proves the arm orders without
a list of small continued fractions.

The closed original chain is S3: it blows down to the empty smooth
corner plumbing. A regular fiber at E is the link of its curvette,
parametrized at the original corner by w~z^q,y=z^r. It is T(q,r).
Equivalently its punctured graph is the disk-base Seifert exterior
with orders q,r and regular fiber

    h=qr*mu+lambda,                                    (9)

in the meridian and preferred longitude of that torus knot.

After the first free blowup, the regular-fiber meridian changes
from mu to mu+h. The remaining negative-two arm has fraction
m/(m-1), with an empty arm for m=1. Propagating the actual meridian
from the punctured endpoint gives

    mu_R=m*(mu+h)-(m-1)*h
        =m*mu+h=(qr+m)*mu+lambda.                      (10)

Thus it is integer n-surgery with the surgery dual marked, not only
a closed Seifert manifold of the correct determinant. Descriptions
are up to a consistent simultaneous orientation/mirror choice.

The Euler number changes from 1/(qr) to 1/(qr)+1/m. The closed
determinant is therefore n=qr+m. Independently, deleting the marked
vertex has determinant n-1 (replace m by m-1, or delete E at m=1).
The marked meridian has order n, and H1 is cyclic. Along the free
chain the central fiber satisfies

    h=m*[k] in H1(M),                                  (11)

where k is the omitted knot, with one consistent sign. Since 1<=m<n,
that central class is nonzero.

All coefficients of H only alter the finite jet phi and hence the
locations of free centers. None changes this graph or its marking.

## 3. Atoroidality and the actual target-radial section

The source in (2) has Seifert base S2 with at most three genuine
cones and nonzero Euler number. It is irreducible, by the retained
Hatcher criterion; the exceptional reducible cases have zero Euler
number or incompatible first homology. An incompressible vertical
torus would enclose at most one cone on one side and would compress.
An incompressible horizontal torus is impossible, because Hatcher
Proposition2.2(b) requires Euler number zero for a closed horizontal
surface. Hence the closed source is atoroidal. The complete relevant
criterion and proof, printed pages26-27, were visually checked.

The torus-knot exterior is irreducible, boundary-incompressible and
atoroidal when q,r>1, by the horizontal/vertical surface argument
retained in PF2RULED41. Its peripheral structure, including mu_R,
has been computed above.

For an actual finite Keller map, cofinal product neighborhoods of
the end pair identify the closed radial-section group, exterior
group, peripheral subgroup and omitted meridian with this model.
There is a new issue here: the closed model need not be a lens
space or have finite fundamental group.

It is not necessary to classify closed Seifert manifolds from
their abstract groups. By (11), its regular-fiber central element
is nontrivial. A nontrivial free product of nontrivial groups has
trivial center. Kneser--Milnor prime decomposition, followed by
the Poincare theorem for any simply connected prime summands,
therefore makes the actual closed radial section prime. Its
finite first homology excludes S1 times S2, so it is irreducible.
Only this removal of homotopy-sphere summands uses retained
Perelman; no general group recognition of a small Seifert space
is assumed.

Primitivity and n>1 make the actual omitted knot nonzero in H1.
As before, this makes its exterior irreducible: the ball bounded
by a sphere cannot contain the knot. Its group is the noncyclic
torus-knot group, so a boundary compression cannot occur.
The peripheral isomorphism to the model exterior is realized by
a map of pairs to that aspherical, sufficiently large target.
The full Waldhausen1968 Theorem6.1 argument of PF2RULED41 applies:
its exceptional product case has the wrong boundary count, and
the resulting covering has degree one. Thus the exteriors are
homeomorphic with mu_R preserved. Extending across that filling
also identifies the actual closed radial pair with the model.
This transfers both atoroidality assertions rigorously.

The only pair with a solid-torus exterior in our range is (q,r)=(1,2).
It is the already recognized Wright/lens-core case.

## 4. Quartic obstruction with an atoroidal closed source

There is no connected simple degree-four cover of a finite iterated
torus knot with all of the following properties:

* full S4 monodromy and identity preferred longitude;
* closed source M irreducible and atoroidal;
* H1(M)=Z/n with n>3 odd;
* unique ramification knot k primitive, with atoroidal exterior.

This replaces the lens hypothesis in PF2RULED41. We reuse its exact
thirty-subgroup/peripheral layer and check every change of hypothesis.
A single incompressible embedded torus now contradicts closed
atoroidality; it need not inject into a finite fundamental group.

### 4.1 Regular fibers with repeated cone factors cannot be primitive

For a closed Seifert filling whose original disk piece has cone
orders A,B with g=gcd(A,B)>1, quotient its homology by the regular
fiber h. The orbifold abelianization still surjects onto C_g:
send the first two cone generators to 1,-1 and every further cone
generator to zero. This works regardless of the third filling cone.
If h=[k] were primitive in cyclic H1, that quotient would be zero.
Fiber filling instead kills h and contradicts its nonzero class.

This replaces the part of the old proof that used a lens filling
to discard a third cone.

### 4.2 Four-letter and three-letter transposition patterns

All boundary-incompressible four-letter patterns glue to unbranched
companion covers along incompressible tori. Such a torus is forbidden
in the closed atoroidal M. The remaining solid pattern is p=3 with
a regular marked fiber and a genuine order-three cable exterior.
Its gluing to the incompressible companion cover gives a
non-boundary-parallel incompressible torus in E(k), also impossible.
The p=2,(2,2) alternative is still ruled out by the exact peripheral
subgroup conditions.

For the three-letter transposition pattern p>=4, the disk piece has
cones p/2,p. Its incompressibility and closed atoroidality force the
partially filled companion Y to compress, hence to be a solid torus.
The compression-ball alternative is excluded by nonzero [k] on
the opposite side, as in the retained proof. Section4.1 now excludes
the primitive regular knot, without requiring the resulting M
to be a lens space.

For p=2 the disk piece is solid with one order-two cone. If Y
were boundary-incompressible, E(k) would contain the genuine cable
torus. Thus Y is solid. Replacing this piece by the exact standard
index-three branch solid torus, with m0=3mu_J,l0=lambda_J, produces
a lens/core alternating A_n cover of the companion. Its H1 is the
same Z/n. PF2RULED39 excludes it. That replacement is a lens space
because its two pieces are solid, not because the original M
was assumed to be one.

### 4.3 Two-letter patterns

The winding p is odd. For even cable slope a the disk piece has
cones p,p and the marked knot is regular. Closed atoroidality
and compression force Y solid. Section4.1 excludes primitivity,
including every possible nonfiber filling; fiber filling kills [k].

For odd a the disk cones are p,2 and k is the order-two fiber.
Again Y is solid. Fiber filling kills 2[k], nonzero at odd n.
For a nonfiber filling, replace this disk piece by the standard
simple-branch solid torus. Its actual boundary homology map agrees
on m0=2mu_J,l0=lambda_J: H1 of the disk piece is Z since p,2 are
coprime, l0 maps to the generator of H1(V), and m0 maps to zero.
Mayer--Vietoris preserves H1=Z/n for any gluing slope.
The replacement source is Y union a solid torus, hence a lens
space with a ramification core. It is a full-S4 identity-longitude
cover of the simpler companion, excluded by PF2RULED39.

Unlike the former lens-source argument, no distance-one hypothesis
on the original filling is needed for this replacement.

### 4.4 The final three-letter three-cycle pattern

The full-S4 companion case has a dual graph cycle and infinite H1.
In the S3 companion case there is an actual unbranched copy of
X_J and a one-boundary partially filled three-sheeted cover Y.
The branched piece has annular base with a possible cone p/3.

If Y were incompressible, all gluing tori would inject, contrary
to closed atoroidality. Thus Y is solid. Its filling is not along
a fiber, which would kill [k]. After filling, the branched side
must itself be a solid torus: otherwise its gluing to the
boundary-incompressible X_J gives an incompressible torus in M.
If this solid torus had a genuine cone, deleting the regular
marked knot would give a genuine cable torus in E(k).
Therefore p=3 and the filling distance from the fiber is one,
with E(k)=X_J.

Atoroidality makes the iterated companion a torus knot.
The unchanged three-sheeted torus argument of PF2RULED41 then
forces it to be a trefoil or its mirror: its closed three-fold
simple cover has genus-zero Seifert base and finite H1, and the
solid Y makes its ramification a lens/core knot; the torus lemma
forces H1 zero and the remaining disk cones A,2 force A=1.

Use the exact slopes already proved there. With epsilon=1 or -1
for the trefoil mirror choice,

    meridian(Y)=epsilon*m0+l0,
    m0=-mu, l0=lambda+(a/2)*mu,
    distance(meridian(Y),h_B)=|a/2-3epsilon|.

Distance one makes a/2=3epsilon+delta, delta in {1,-1}.
Thus M is filling of that trefoil exterior along

    (2epsilon+delta)*mu+lambda.                        (12)

It has first homology of order one or three, contradicting n>3.
We do not call these fillings impossible in an arbitrary closed
atoroidal source. This is exactly the new order threshold.

The torus base case is excluded by the all-degree torus lemma
because [k] would be zero; the unknot has cyclic monodromy.
All finite iterated-cabling possibilities have now been covered.

### 4.5 The order-three exception is an actual topological cover

The n>3 threshold cannot be dropped from this topological statement.
Keep the explicit permutations of PF2RULED41, but take cable (3,8)
over the positive trefoil. They still define a simple full-S4 cover,
with identity longitude. Its companion filling in (12) is 3mu+lambda.
The closed source has Seifert orders 2,3,3 and H1=Z/3; its Euler
number is nonzero, so it is irreducible and atoroidal.
The pattern filling distance is |4-3|=1, so the ramification
exterior is the trefoil exterior and is atoroidal. Its class
8mu=2mu modulo three is primitive.

This is an actual topological branched cover, not a finite Keller
map or an assertion of polynomial infinity realization. It prevents
a false extension of the new theorem to n=3. The older lens-source
theorem still applies there because this particular source is not
a lens space.

## 5. Explicit affine exact-volume controls for every coprime pair

The monomial whole family can be constructed, not merely treated
as formal chart data. Fix coprime q,r with r>=2,q>=r-1, and put

    m=q-r+2, k=ceil(q/m), f(v)=v^r-1,
    A=v*f(v)/z^m, B=v*f(v)^k/z^q,
    D_(q,r)=C[z,v,A,B] inside C[z,z^-1,v].               (13)

Here A,B name ring generators, distinct from the Bezout integers
of Section1. If r=2 then q=m,k=1 and A=B, which causes no problem.

For a root ci of f let g_i=f/(v-ci). The principal opens at f and
at all v*g_i cover Spec D_(q,r), by sum_i v*g_i=r*(f+1).
They recover whole plane opens:

* at f, x=B/f^k gives v=z^q*x and A=z^(q-m)*x*f;
* at v*g_i, t_i=A/(v*g_i) gives v=ci+z^m*t_i and
  B=z^(km-q)*v*g_i(v)^k*t_i^k.

Both exponents q-m=r-2 and km-q are nonnegative. Thus these are
actual polynomial chart rings. The source plane A2_(x,z) is covered
by D(z) and D(f); at z=0, f=-1. The r omitted root lines are
entire A1's, with charts

    x=ci*z^(-q)+z^(2-r)*t_i, y=z^r.                   (14)

Let mu_r act by z->zeta*z, v->zeta^q*v, A->zeta^(q-m)*A
and B->B. Coprimality makes its action free and transitive on
the omitted lines. Along the retained source zero line it is the
reflection (x,z)->(x,zeta*z). Its affine quotient N_(q,r) is smooth,
has marked source A2_(x,y), sole omitted A1, and the required
orders and source volume divisor R.

These controls even have an explicit regular volume primitive.
The commuting torus action z->lambda*z, v->v descends to N_(q,r)
and acts on its source by x->lambda^(-q)*x,y->lambda^r*y.
Its Euler vector field contracts the source volume to

    theta=(-q*x dy-r*y dx)/(r-q),  dtheta=dx wedge dy.   (15)

The denominator is nonzero for coprime q,r in this range.
The action is regular globally, so the contraction is regular.
On a root omitted chart its numerator is exactly

    -r*m*z*t_i dz-r*z^2 dt_i,

also verifying the cancellation of the apparent source-coordinate
poles directly. No target equivariance is assumed.

Thus affine geometry and exact volume alone really allow every
coprime pair in this range. Their finite quartic targets are
excluded by the marked-end argument, not by a false geometric
nonexistence statement. An affine converse for general mixed
Laurent data remains unproved.

## 6. Apply to actual finite Keller maps and retain the real gap

If a coprime model has an actual finite Keller target, the marked
end has cyclic H1 and primitive omitted knot. In every field
degree the torus-infinity lemma excludes target torus/unknot
infinity. The retained coordinate-degree/shear arguments exclude
coprime boundary-normalization degrees and every prime boundary
coordinate degree; all such degrees are composite and at least four.

For field degree four, PF2RULED36 gives intrinsic identity longitude,
and the boundary knot is a finite iterated torus knot. The actual
radial pair has both atoroidality properties by Section3.

The integer n=qr+q-r+2 is odd: coprimality excludes the only even
parity case q,r both even. Apart from (q,r)=(1,2), one has q,r>=2,
hence n>=7. Section4 therefore excludes every such quartic map.
The exceptional Wright pair is already excluded by PF2RULED35/41.

Now combine PF2RULED39-41. Every remaining regular-coordinate
quartic normalization has r>=2,q>r; this result adds gcd(q,r)>1.
Since gcd(q,r)=gcd(q-r,r), the positive difference q-r is at least
the smallest prime divisor ell(r). This proves (3), including
q>=r+2 and q>=2r for prime r.

The full Laurent support, common coefficient-r residue, monic norm
H_C^r*H_0, horizontal conductor and actual finite algebra remain
mandatory. The remaining noncoprime sector can have further
characteristic branches and essential gluing tori; it must not be
identified with the atoroidal end proved here. No entire noncoprime
sector, arbitrary omitted-prime packet, no-regular-coordinate case
or higher finite Wright degree is excluded by this theorem.

All universal compactification, surgery, group-recognition and
cabling arguments are written with explicit primary imports.
The checker supplies exact arithmetic and graph controls and pins
the sources, not a formal verification or independent full reproof.

## 7. Integral-power removal excludes all prime boundary orders

This extends the coprime leading-pair result to the entire
one-characteristic-pair sector in the remaining q>r region.

Keep the normalized root Laurent series h(z) of PF2RULED40 and the
reference germ at the original corner

    w=1/h(z), y=z^r, w=1/x.

Its parametrization is primitive. Indeed C((z))/C((z^r)) is cyclic,
and the stabilizer of h has order gcd(r,support(h))=1. Thus
C((z^r,h(z)))=C((z)). In particular its characteristic gcds
eventually reach one.

If q is not divisible by r, its first nonintegral exponent is q.
It has one characteristic pair precisely when gcd(q,r)=1,
the case already excluded above. If that gcd is greater than one,
primitivity requires at least one additional characteristic pair.

Suppose now r divides q. Put

    j=max {s:c_s!=0 and r does not divide s}.

Such j exists by the support-gcd condition, and r-1<=j<q.
Separate the r-divisible powers of h from the other powers.
The reciprocal of its divisible part is a series in y. Its
first nondivisible correction comes from the leading remaining
term c_j*z^(-j) and has order

    beta=2q-j,
    w=phi(y)+C*z^beta+O(z^(beta+1)), C=-c_j/c_q^2!=0.   (16)

Here phi is the polynomial consisting of the integral-power
terms preceding beta. Higher powers of the nondivisible part
have strictly larger order since q>j, so this leading term
cannot cancel. The local coordinate change u=w-phi(y) has
Jacobian one and leaves y unchanged.

If gcd(r,j)=1, then gcd(r,beta)=1. The reference branch has
exactly the one characteristic pair (r,beta), so its link is
T(r,beta). The whole actual omitted chart has variation

    w=1/[h(z)+z^(2-r)*t],
    delta w=-c_q^(-2)*z^(2q-r+2)*t+higher terms.

Relative to the leading u~C*z^beta, its order is

    k=2q-r-beta+2=j-r+2>=1.                            (17)

The same fixed-transverse-coordinate argument as Section1
therefore gives exactly k free blowups and a nonzero affine
coordinate on the entire omitted line. Equivalently the last
toric divisor has ordinary canonical discrepancy r+beta-1;
the original source volume w^-2 dw wedge dy has order
r+beta-1-2q there, and k free blowups make that order one.

There is a boundary detail to keep. The original deleted D_x
is w=0, hence u=-phi(y); it detaches from the reference at
its integral contact q/r, strictly before the beta characteristic
point. It meets one earlier exceptional transversely. Its only
other deleted neighbor is the zero leaf D_y at infinity.
The same one-attachment R6 removes this pair with no extra
summand, even when the attachment vertex is internal to the
exceptional chain. The remaining local exceptional chain is
the usual unimodular resolution in (u,y), and its two arms
have orders r,beta. No local coordinate change has been
mistaken for a source polynomial automorphism.

Consequently the full marked end is

    N=r*beta+k surgery on T(r,beta),
    closed Seifert orders r,beta,k,
    H1=Z/N with primitive dual omitted knot.            (18)

As before both the closed source and the marked exterior
are atoroidal, and the peripheral comparison identifies the
actual radial pair. The coefficient is odd, since

    N=(r-1)*(beta-1)+2q+1,                              (19)

and coprimality makes the first product even. Also N>=7.
The closed-atoroidal quartic theorem excludes every such model.
All integral and lower Laurent terms in (16) have been retained;
they locate the free centers without changing this end pair.

If r is prime, a q not divisible by r is already coprime to r,
while for q divisible by r the index j in (16) is automatically
coprime to r. These two cases exhaust all prime r. Together
with the previous simple-zero and borderline exclusions,
no actual quartic regular-coordinate normalization can have
prime boundary vanishing order.

Thus the current necessary region is more sharply

    r composite (in particular r>=4),
    q>r, gcd(q,r)>1, q>=r+ell(r),
    reference germ has at least two characteristic pairs. (20)

If r divides q, its largest nondivisible support index j must
also satisfy gcd(r,j)>1; otherwise (18) excludes it.
The remaining characteristic gcd must drop further at a later
index by primitivity. The extended translation vector field of any
remaining quartic model vanishes along R to order at least two: its
free r=2 sector and its r=3 sector are now entirely excluded. This identifies the genuine additional
branching that the next marked-end argument must retain,
rather than creating a deeper bounded coefficient search.

The full single-characteristic-pair end calculation, including
(18), is all-degree geometry. Its primitive marked class excludes
target torus/unknot infinity and prime boundary-coordinate degrees
in every possible finite Keller degree. The finite-map exclusion
here remains quartic. No whole multi-characteristic sector,
extra-prime packet or higher finite Wright degree is closed.
