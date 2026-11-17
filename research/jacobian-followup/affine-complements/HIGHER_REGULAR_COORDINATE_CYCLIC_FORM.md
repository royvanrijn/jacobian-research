# Higher regular-coordinate vanishing: cyclic charts and an all-degree diagonal gap

14 November 2026. PF2RULED40 addresses the regular-coordinate sector left
after PF2RULED39. It is a necessary global normal form and an exact
volume criterion for an existing smooth affine normalization. It does
not assert an affine converse for arbitrary chart data.

Let N be smooth affine, with marked open U=A2_(x,y), sole complement
R=A1, and div_N(dx wedge dy)=R. Suppose the source coordinate y is
regular and, after subtracting its boundary value,

    v_R(y)=r>=2,   v_R(x)=-q<0.                         (1)

There is no bound on r, q, polynomial degree or finite field degree.
Then

    q>=r-1.                                            (2)

Normalize y=z^r. The entire cyclic cover has r+1 plane charts: its
retained source plane and one plane for each of its r omitted lines.
With z on every chart and suitable fiber coordinates t_i, i mod r,

    X_src=A2_(x,z),
    X_i=A2_(z,t_i),
    x=z^(2-r)*t_i+h_i(z),                              (3)
    h_i(z)=sum_(j=r-1)^q c_j*zeta^(i*j)*z^(-j),
    c_q!=0, zeta a primitive r-th root.

(The source chart X_src is separately named; indices i on boundary
charts run from 0 to r-1.) Every pairwise overlap is the same
Gm_z times A1. The exponents in the nonzero support satisfy

    gcd(r,{j:c_j!=0})=1.                               (4)

The source volume is exact on N if and only if

    c_r=0,                                            (5)

where an absent coefficient is zero. In particular q=r is impossible
for an actual finite Keller extension, in every finite field degree.

The source derivation partial_x extends to a locally nilpotent
derivation on N. It is free for r=2; for r>=3 its vector field has
zero divisor exactly (r-2)R. Extension of this action to the source
normalization does not imply invariance of the target algebra.

The claims are written all-r proofs with exact supporting controls.
They are not formal verification, independent implementation or
external review, and have no novelty-priority assertion.

## 1. Retained affine lines cannot acquire endpoints

Every source coordinate line is closed in N and is disjoint from R.
Indeed, the normalization of its closure is a normal affine curve
containing A1 as an open. Any added point would have to be its one
point at infinity in the smooth completion P1. Adding that point
makes P1, which is not affine. Normalization is finite for these
complex affine curves, so its affineness is assured.

Thus L={y=0} is closed and disjoint from R, and

    div_N(y)=L+rR.                                    (6)

The pole q in (1) really is positive. If x were regular too, normality
would give O(N)=C[x,y] intersect O_(N,eta_R)=C[x,y], and the given
source plane open would already be all of N. The minimum pole over
all source coordinates equals q by PF2RULED38 Section1. Also its only
regular source coordinates are affine functions of y. In particular
r is unchanged by choosing another regular source coordinate and
subtracting its boundary value.

## 2. Normalize the finite regular-coordinate root

Let rho:N'->N be the normalization in C(x,z), y=z^r. It is finite,
connected and cyclic of degree r. Above L it is locally the smooth
map (x,z)->(x,z^r), with ramification divisor (r-1)L'. Above R it is
etale: locally y=s^r*a with a a unit, so the normalization is obtained
by adjoining a root of a and putting z=s*a^(1/r).

The restriction to R is a finite etale cover of A1 over C. It splits
into r copies R_i of A1, each mapping isomorphically to R, and the
deck group permutes them transitively. This also follows from the
absence of nontrivial finite etale covers of the complex affine line.
There are no other possible singular points: the cover is etale
away from L and is smooth at L. Therefore N' is smooth affine and

    rho^(-1)(U)=A2_(x,z),  div(z)=L'+sum_i R_i,
    div(rho^*(dx wedge dy))=(r-1)L'+sum_i R_i.          (7)

Here the root is taken at the finite boundary value of a regular
coordinate. It is different from the root of a pole coordinate in
PF2RULED37; the locations of the special fibers are different.

Remove L' and all R_j with j!=i, and call the resulting open X_i.
It is affine: these are effective Cartier divisors on the smooth
affine N', and the complement of an effective Cartier divisor in an
affine scheme is affine. The map z:X_i->A1 is faithfully flat, with
generic fiber A1 and every geometric fiber a reduced A1. Its zero
fiber is precisely R_i; every nonzero fiber is a source line.

Miyanishi TIFR60, ChapterI Theorem4.1.1 (printed page24 / PDF page40)
makes it an affine-line bundle. Both its linear bundle and its torsor
class on A1 vanish, so X_i=A2_(z,t_i). This is the same explicitly
pinned primary import used in PF2RULED38; see
[sources/coordinate-pole-sector/manifest.json](sources/coordinate-pole-sector/manifest.json).
No formal local neighborhood is being substituted for the whole chart.

Together with the retained source X_src these charts cover all of N'.
For any two distinct charts, the overlap is exactly z!=0 and is
Gm times A1. On the source/boundary overlap an automorphism of the
polynomial algebra over C[z,z^-1] gives

    x=a_i*z^k_i*t_i+h_i(z), a_i!=0, h_i Laurent.

The pulled-back volume equals

    r*z^(r-1)*dx wedge dz
      =-r*a_i*z^(r-1+k_i)*dz wedge dt_i.

Its order along R_i is one by (7). Thus k_i=2-r. Rescale t_i
to set a_i=1, obtaining (3).

## 3. Gauge, poles and separation

Replacing t_i by t_i+p_i(z), with p_i polynomial, removes exactly
the terms of h_i of exponent at least 2-r. There is a unique
representative containing only powers z^(-j), j>=r-1. These changes
do not alter the source coordinate x.

This representative cannot be zero. Otherwise the source curve
x=beta, z!=0 in X_i has t_i=beta*z^(r-2), which extends across zero
to a finite point of R_i. But the source affine line x=beta is closed
in the affine N', by Section1. This contradicts its acquiring an
endpoint on R_i. Thus the representative has a largest pole q>=r-1.
Its leading coefficient is a nonzero constant, and its pole strictly
dominates the t_i term of order r-2. Hence this q is exactly
-v_(R_i)(x)=-v_R(x), proving (2), not merely a chart bound.

Label R_i so a generator g, given on the source by g(z)=zeta*z,
g(x)=x, sends R_i to R_(i+1). Uniqueness of the normalized Laurent
representative gives

    h_(i+1)(zeta*z)=h_i(z),
    g^*(t_(i+1))=zeta^(r-2)*t_i.                       (8)

This yields the displayed formula for all h_i from h_0=h. No
additional polynomial translation survives the chosen gauge.

The h_i must all be different. If h_i=h_j with i!=j, the overlap
identifies the two planes by t_j=t_i for z!=0. The two zero lines
are still distinct and their equal finite parameters give a limit
of the diagonal outside the diagonal, violating separatedness of N'.
Equivalently the two maps Spec C[[z]] into these charts agree on
Spec C((z)) and differ at their origins.

For the orbit of h under z->zeta*z, the stabilizer has order
gcd(r,{j:c_j!=0}). This is elementary character theory in the cyclic
group: g^a fixes h exactly when r divides a*j for every support j.
Separatedness therefore gives (4).

These are necessary conditions on existing affine N. No converse
affineness assertion follows here just from unequal Laurent parts
or from separatedness of the abstract gluing. Outside the explicit borderline family of Section6, Laurent controls
in the receipt are not claimed geometric realizations. None is a
Keller realization.

## 4. One common residue is the complete volume test

Write omega=rho^*(dx wedge dy). Polynomial primitives on the charts are

    theta_src=r*z^(r-1)*x dz,   theta_i=r*z*t_i dz,
    theta_src-theta_i=r*z^(r-1)*h_i(z) dz.               (9)

Only j=r contributes a residue, equal to r*c_r*zeta^(i*r)=r*c_r.
It is the same for every omitted chart.

If omega is exact, restrict a global primitive to each plane.
Every closed polynomial one-form on a plane is polynomial-exact.
Thus the difference in (9) must be an exact regular one-form on
Gm times A1, and has zero residue. This proves necessity of c_r=0.

Conversely, if c_r=0, set

    K_i(z)=r*sum_(j!=r) c_j*zeta^(i*j)*z^(r-j)/(r-j).

Then dK_i=theta_src-theta_i on the common overlap. These functions
and their differences form a Cech 1-cocycle for O_(N') on the affine
plane cover. Since N' is affine, H1(N',O)=0. Choose regular chart
functions f_src,f_i with f_src-f_i=K_i. The forms theta_src-df_src and
theta_i-df_i agree and give a global regular primitive of omega.
This argument uses affineness of the existing normalization, rather
than asserting it for arbitrary gluing data.

Exactness on N implies exactness on N'. For the reverse, average
the primitive over the finite deck group. An invariant regular
one-form descends regularly: this is etale descent away from L',
and near L' the action is (x,z)->(x,zeta*z). The coefficient of dx
is a polynomial/power-series in z^r, while the invariant coefficient
of dz has the form z^(r-1)*b(x,z^r); hence it is b*dy/r. The same
calculation works in etale-local regular coordinates. Therefore the
averaged primitive descends across the branch line as well.
This proves the exact equivalence (5).

For any actual finite extension pi:N->A2_(P,Q) of a polynomial Keller
map, dP wedge dQ=c*dx wedge dy with c!=0, and P dQ/c is a global
regular primitive on N. Thus c_r=0 is necessary in every field degree.
If q=r, c_r=c_q!=0, a contradiction. The full diagonal q=r is excluded,
without coefficient bounds or an assumption about boundary collisions.
Passing the residue test does not construct a finite extension.

## 5. Translation extends, but target equivariance does not follow

On the source take delta=partial_x, with z fixed. In a boundary chart,

    delta=z^(r-2)*partial_(t_i).                       (10)

Since r>=2, this is regular on every chart. It is locally nilpotent
on each chart, and commutes with the deck action in (8). It descends
to N and restricts to partial_x on U. Its derivation on O(N) is
locally nilpotent also directly from O(N) subset C[x,y].

For r=2 it has no zero anywhere, so the induced Ga action is free
(over characteristic zero a nontrivial stabilizer in Ga is Ga).
For r>=3 its vector field has exact zero divisor (r-2)R: it has
that order on every etale boundary chart and has no zero on U.
This is an intrinsic direction up to nonzero scaling, since regular
source-coordinate uniqueness forces any complementary coordinate
to be a*x+p(y). The fixed divisor records a structural distinction
between r=2 and r>=3.

Nothing proves delta(C[P,Q]) subset C[P,Q]. The action on the
normalization must not be pushed to the target without this extra
invariance. This is a possible finite-algebra gate, not an exclusion.

## 6. Explicit affine controls for every borderline q=r-1

The entire borderline family has actual smooth affine geometric
realizations. This is a converse for this one family only, not for
arbitrary Laurent data in (3). Fix r>=2, put q=r-1, and let

    f(v)=v^r-1,   w_j=v*f(v)^j/z^j,  1<=j<=q,
    B_r=C[z,v,w_1,...,w_q] inside C[z,z^-1,v].          (11)

Write ci for each root of f and g_i(v)=f(v)/(v-ci).
The principal opens at f and at the v*g_i cover Spec B_r: these
polynomials in v have no common zero. More explicitly,
sum_i v*g_i(v)=r*v^r=r*(f+1), a unit Bezout identity.
Their rings are explicit:

* at f, x=w_q/f^q gives v=z^q*x and
  w_j=z^(q-j)*x*f^j, so the open is D(f) in A2_(z,x);
* at v*g_i, t_i=w_1/(v*g_i) gives v=ci+z*t_i and
  w_j=v*g_i(v)^j*t_i^j, so the open is D(v*g_i) in A2_(z,t_i).

All identities hold in the displayed subring of a field; no extraneous
component has been introduced by only writing equations. These opens
prove smoothness and recover the entire charts. At z!=0 the ring is
C[z,z^-1,v]=C[z,z^-1,x]. The source charts D(z) and D(f) glue to
the entire plane A2_(x,z), since f=-1 along its z=0 line. Their
complement consists of r disjoint affine lines z=0,v=ci. Every whole
boundary plane maps in by v=ci+z*t_i and the displayed w_j; near
z=0 it is the principal open above, and off zero it is the source
overlap. Thus its chart really is the entire plane, with transition

    x=z^(1-q)*t_i+ci*z^(-q)
      =z^(2-r)*t_i+ci*z^(1-r).                         (12)

Now let the cyclic group act by

    z -> zeta*z, v -> zeta^q*v,
    w_j -> zeta^(q-j)*w_j.

It preserves B_r, fixes x on the source and permutes the r omitted
lines freely and transitively. The quotient N_r=Spec(B_r^(mu_r))
is affine. It is smooth: near the retained zero line the action is
the reflection (x,z)->(x,zeta*z), and near an omitted line its
stabilizer is trivial. Away from these sets the action is free.

The marked source quotient is exactly A2_(x,y), y=z^r. Its sole
complement is R=A1, with v_R(y)=r and v_R(x)=-(r-1).
The chart volume calculation gives div_N_r(dx wedge dy)=R.
Moreover x dy is already a regular global primitive: on a root chart

    rho^*(x dy)=r*(z*t_i+ci) dz.

It descends regularly as in Section4. Thus every r>=2 gives an
actual smooth affine exact-volume control, and r>=3 gives q>=2.
By PF2RULED38 none of these latter marked surfaces has a hidden
source coordinate of pole one. No finite Keller target is supplied.
This family rules out a geometric claim that exact volume, the
pole invariant, the root gluing or the extended Ga action alone
excludes every higher regular-coordinate order.

## 7. The actual finite algebra has a repeated special norm fiber

Now additionally assume pi is the actual finite normalization of a
Keller map, with exactly one omitted prime of generic index two.
PF2OMIT7 gives residue degree one on R. The all-degree Euler/Chau
argument retained in PF2RULED37 Section5 makes C=pi(R) nonimmersive.

For every beta, L_beta={y=beta} is closed in N and is disjoint from R.
Its image D_beta is finite with birational immersive normalization L_beta:
factor the finite map through the normalization of D_beta, which is A1.
The resulting polynomial A1->A1 has no finite critical point, since
the original Keller map is etale on L_beta; it is affine linear.
Thus D_beta differs from the nonimmersive C.

With H_C,H_beta reduced irreducible equations, push forward the
principal divisors. By (6), up to nonzero constants,

    Norm(y)=H_C^r*H_0,
    Norm(y-beta)=H_beta for beta!=0.                    (13)

Every nonzero-parameter norm has reduced irreducible divisor, so
the norm tower proves C(N)=C(P,Q)(y), just as in PF2RULED38.
The monic minimal polynomial chi(T) therefore has degree d and

    chi(0)=constant*H_C^r*H_0,
    chi(beta)=constant_beta*H_beta for beta!=0.         (14)

The special norm fiber has multiplicity r along C, whereas every
nonzero member is reduced irreducible. This is a polynomial family
of plane curves, not a linear pencil or a target polynomial fibration.

The primary Gwozdziewicz one-line injectivity theorem, pinned in the
same coordinate-pole manifest, forces every D_beta to self-collide.
Consequently C[P,Q,y] is not normal; its normalization is O(N), and
its conductor support has a component dominating the y-line. This
extends that PF2RULED38 finite-algebra condition to all r>=2; it does
not turn generation of the function field by y into generation of
the normal ring.

The composite N'->A2 has degree rd and full C_r wreath S_d monodromy.
Indeed its entire branch locus is C union D_0; a meridian at C is r
disjoint transpositions and a meridian at D_0 is one r-cycle. The
intermediate cyclic cover gives d blocks of r sheets. Conjugating
that single-block r-cycle by the transitive S_d quotient generates
all r-cycles independently on the blocks, giving the full wreath
product, by the PF2RULED37 Section5 argument. Its source Jacobian is
c*r*z^(r-1), so this lifted cover is not a Keller map for r>=2.
Removing the other boundary charts loses finiteness over the target;
one cannot apply a one-boundary finite-map exclusion to X_i.

## 8. Consequences and the remaining gate

For an actual quartic normalization with exactly one omitted prime
of index two, PF2OMIT7 gives the smoothness and boundary line used
above. If it has a regular source coordinate, PF2RULED39 first forces
r>=2 and PF2RULED36/PF2RULED38 force intrinsic minimum pole q>=2.
The present result sharpens the necessary region to

    r>=2, q>=2, q>=r-1, q!=r,
    c_r=0, gcd(r,{j:c_j!=0})=1.                        (15)

The inequality, residue condition, support condition and Ga extension
hold independently of finite field degree under the stated smooth
one-line and volume hypotheses. The quartic lower bound q>=2 uses
the separate quartic theorem and is not asserted here in all degrees.

When q<r, the only possible value is q=r-1; then the normalized
Laurent part is a single monomial, with coprime support automatically.
Section6 realizes this entire family by actual smooth affine models
(up to rescaling its nonzero coefficient), all with exact volume. When r=2,
at least one odd exponent must occur to separate the two boundary charts,
and c_2=0.
The data (r,q,h)=(3,2,z^-2) are realized by Section6;
(2,3,z^-3) and (3,4,z^-4) remain formal controls against false
stronger local exclusions here. None is a Keller construction.

The full finite-quartic problem in the r>=2 sector remains open,
as do the no-regular-coordinate case, extra omitted primes, all
higher finite Wright degrees and the complete research goal.
A useful next step is to use the cyclic chart structure to determine
the marked end pair or constrain the actual finite target subalgebra.
Do not restart a bounded source-automorphism or deeper coefficient
search merely because those global steps remain open.
