# Seifert framing cancels every simple critical-pair twist

PF2RULED126, 16 November 2026. This proves a componentwise peripheral
theorem for EVERY genuine finite intermediate of an actual complex plane
Keller normalization. Involutory branch inertia forces its preferred
longitude to fix EVERY moved letter. For binary affine intermediates this
excludes the order-four and mixed peripheral profiles left admissible by
PF125: the longitude has order at most TWO and every critical shear is
ZERO. It does not exclude the entire binary affine monodromy family or
close the full finite algebra. The shared HC4 empty-fibre gate remains OPEN.

## 1. Scope of the componentwise theorem

Retain the ORIGINAL finite normalization N, its WHOLE source A2 and ALL
omitted primes. Let N->M->A2 be ANY nontrivial genuine sheet-block
intermediate. Choose a full branch component C whose generic monodromy
q is an INVOLUTION. Other branch components and their indices are
unrestricted. Let K be C's oriented knot at infinity, mu a positive
meridian and lambda the preferred Seifert longitude.

PF123 gives disjoint critical A1s in M and cyclic quotient-axis germs
with n<=e-1. Over this C every critical index is two; hence M is SMOOTH
along its critical curves E_1,...,E_r above C. PF119's LOCAL simple-axis
calculation applies there. The map E_j->C^nu is locally unramified: each
chosen image branch has consecutive orders (m,m+1), so its source
parameter is already its normalization parameter. Since E_j=A1 and the
finite image normalization is A1, that map is an isomorphism. This also
works when M has singularities elsewhere.

Write h(t)=(X(t),Y(t)) for C's polynomial normalization, set

    g(t)=gcd(X'(t),Y'(t)),  M_C=deg(g),
    G(t)=h'(t)/g(t).

Scale g to be monic. The components of G are coprime; choose a polynomial
vector w(t) by Bezout so det(G(t),w(t))=1. Thus w is a nowhere-zero
normal direction to the normalized Gauss tangent, including at cusps.
The finite nonimmersion parameters have branch types (m,m+1), with
ord(g)=m-1. There is no restriction on their number or on the contact
or number of different branches over a finite value in this theorem.

The conclusion is

    rho(lambda)(x)=x for EVERY x in supp(q).             (1)

Preserving each transposition pair would allow its two letters to swap;
equation(1) rules out those swaps as well.

## 2. The holomorphic normal framing has linking 2 delta + M_C

Choose a sufficiently large target ball containing every finite singular
value of C and all intersections with other branch components. The
normalization of its part of C is a compact disk. Along its boundary,
push h(t) in direction w(t) by a small nonzero complex epsilon and move
the resulting curve radially to the sphere. Call that framed parallel
nu_w. It is a push-off of K, not a new choice of preferred longitude.

The linking n=lk(K,nu_w) equals the intersection number of the two
normalized disk chains h and h+epsilon*w inside the ball. This is the
ordinary boundary-linking/intersection identity, with complex orientations.
All intersections occur near finite singular values when epsilon is
sufficiently small. There are none near an immersed single branch,
because w is transverse to its Gauss tangent. For different local
branches i,j, the two ORDERED intersections contribute twice their
intersection multiplicity I_ij. This includes arbitrary tangency and
multiple values, not only ordinary nodes.

For a single branch of type (m,m+1), choose its local parameter so
x=t^m, y=c*t^(m+1)+higher terms, c!=0. A transverse normal displacement
can be deformed locally to a constant vertical one without changing the
intersection number. The equality x(s)=x(t) gives s=zeta*t, zeta^m=1.
The case zeta=1 gives no intersection for nonzero displacement. Each
of the m-1 other roots gives

    ord_t(y(t)-y(zeta*t))=m+1,

because zeta^(m+1)=zeta!=1. The contribution is therefore m²-1, with
ALL higher analytic terms retained. A smooth branch m=1 contributes0.
Since 2 delta_branch=m(m-1), summing gives

    n=2*delta_finite(C)+M_C,       n=M_C modulo2.         (2)

Equivalently one may use the usual polar-intersection formula, but the
root calculation above supplies the special case actually needed here.
The disk-chain argument uses a fixed large ball and small displacement;
it does not assume that an arbitrary polynomial perturbation preserves
the entire link at infinity uniformly in its parameter.

## 3. The same M_C is the ramified-pair monodromy

Identify each E_j with A1_t. The differential has rank one everywhere
on E_j, even at a nonimmersive image parameter. Its kernel L_j and its
target cokernel are line bundles on A1, hence trivial. The intrinsic
second derivative along L_j is a regular section

    a_j:L_j^2 -> coker(d pi_M)|E_j.

Use w to trivialize the target cokernel and any global nowhere-zero
section to trivialize L_j. At a fold a_j is a unit. At a branch point
with multiplicity m, the PF119 local form is pi_M=(u,F(u,z)),

    E_j={u=phi(z)}, F_z=V(u,z)*(u-phi(z)), V(0,0)!=0,
    ord(phi)=m.

The kernel is generated by partial_z; modulo the image of d pi_M its
second derivative is F_zz|E_j=-V(phi(z),z)*phi'(z). Thus a_j has a zero
of EXACT order m-1. These are all its zeros. Changing the chosen global
trivializations multiplies a_j by a nonzero constant: units on A1 are
constant. Consequently a_j/g has no zero or pole and is constant.

On the large normalization circle the two nearby sheets over the
normal displacement are the two fold roots, governed by

    a_j(t)*s²=epsilon up to nonzero analytic factors.

The small-root pair can be continued by the holomorphic Morse lemma
with parameter on this circle. Its exchange parity is the winding of
a_j, namely M_C modulo2. No ramified pair is joined to a different
pair: E_j->C^nu has degree one. Thus rho(nu_w) restricts to q^M_C on
the ENTIRE moved set of q. Other local cover factors do not change
this statement, even where other target branch components intersect C.

With the positive meridian convention nu_w=mu^n*lambda in peripheral
homology. Since these peripheral elements commute and q²=1, equations
(2) and the pair parity give rho(lambda)|supp(q)=1. This proves(1).
Using an uncorrected holomorphic framing would instead leave the cusp
pair twist; it is essential to pass to the ACTUAL Seifert framing.

## 4. Full filling data and the binary affine exclusion

For every component in Section1 the critical tori have

    (e,f,shear)=(2,1,0).

All remaining tori have e=1, hence shear0, but their longitudinal degrees
need not be one. In PF101's normal filling convention EVERY shear is
therefore zero. Whenever the cores are rationally null-homologous,

    F=F^T,                 F*e=0.                       (3)

No integral-null-core assumption or intermediate coverage theorem is
silently added. The original PF111 integral reconstruction still needs
its ORIGINAL whole-plane/all-omitted-knots coverage.

Now specialize to a binary affine intermediate of degree b=2^a, a>=2.
Use PF125's coordinates q(x)=x+l(x)v, H=ker(l), S={l=1}. An affine map
fixing S POINTWISE is uniquely

    T_u(x)=x+(l(x)+1)u,          u in H.                 (4)

Indeed its affine difference from the identity vanishes on the affine
hyperplane S, so equals (l+1) times a constant vector; invertibility
forces that vector into H. Hence

    rho(lambda)=T_u,   T_u²=1.                          (5)

For u=0 every residual letter is fixed. For u!=0 the residual action on
H is fixed-point-free translation, hence consists entirely of two-cycles.
The ONLY actual filling profiles left are:

| Case | Critical tori | Unramified tori | Total tori |
|---|---|---|---|
| lambda=1 | b/4 copies (2,1,0) | b/2 copies (1,1,0) | 3b/4 |
| lambda=T_u!=1 | b/4 copies (2,1,0) | b/4 copies (1,2,0) | b/2 |

Thus PF125's order-four profile and its mixed fixed/two-cycle residual
profile are EXCLUDED for ALL actual binary affine intermediates, not
merely for a bounded dimension. Its affine order-four countercontrol
remains correct algebra; it lacks this Seifert-framing compatibility.
The previous PF60/PF125 degree-four longitude theorem is credited.
Here it also follows immediately: for b=4 a nontrivial T_u is a single
transposition on H, but a preferred longitude has even sign. Therefore
u=0. This does not assume a degree-four intermediate to be a plane.

## 5. Binary residual covers are elementary abelian two-covers

Here C is the FULL branch of the binary affine intermediate. Take the
normalization R of the GENERICALLY UNRAMIFIED part of
M times_(A2) C^nu, of total rank b/2 over A1. This discards the doubled
critical sections before normalizing; it does not discard source primes
from the original finite normalization N. The residual cover may be
disconnected. Let v_o be the number of double values with ODD contact.

At either normalization preimage of such a double value the other fold
gives z²=t^k*unit, so residual monodromy is a fixed-point-free involution.
In the Gauss-normal framing it fixes S pointwise: the selected folds
extend there. Thus it is a translation of H by some NONZERO u_i in H.
Even contact gives identity. At the optional ordinary cusp, the local
factorization

    z³-3t²*z+2t³=(z-t)²*(z+2t)

makes the residual sheet a section, so residual monodromy is trivial.
There is no other finite branching: the full intermediate branch is C.
Conjugating to a common normalization base point preserves translations
(they are normal in the affine group of H).

Consequently the residual monodromy is the translation subgroup

    W=span_F2{u_1,...,u_(2 v_o)} subset H,
    u_infinity=sum_i u_i.                              (6)

Each connected residual component is a GALOIS cover of A1 with group W,
has degree w=|W|, and there are (b/2)/w such components. It is NOT
assumed that the two columns belonging to one double value are equal.
The preferred residual longitude is translation by u_infinity; the
meridian framing correction acts trivially on H.

When v_o=0, W=0, R splits into b/2 copies of A1 and lambda=1 in EVERY
binary degree. When v_o>0, set epsilon_infinity=1 if u_infinity!=0 and0
otherwise. Every compactified residual component has

    genus=1-w+(w/4)*(2 v_o+epsilon_infinity),
    infinity_points=w/(2^epsilon_infinity).             (7)

This is ordinary Riemann--Hurwitz: all 2 v_o finite branch points have
inertia two, and infinity has inertia two exactly in the nonzero-sum
case. A negative or nonintegral genus forbids a proposed branch-label
matrix. In particular for EXACTLY ONE odd-contact double value:

* equal columns: w=2, infinity splits, genus0 with two points removed;
* unequal columns: w=4, infinity ramifies, genus0 with two points removed.

Thus EVERY residual normalization component in that case is C*, even
when the original source is imprimitive or has higher visible indices.
These are curve-cover consequences, not a proof that the intermediate
surface is a plane or that the whole monodromy family is excluded.

## 6. A common symplectic form forces trivial longitude in every degree

Suppose the global linear part H_lin of the binary affine group preserves
a nondegenerate alternating bilinear form B on V. No classification of
transvection groups is assumed. This hypothesis includes the binary
symplectic and orthogonal linear subgroups in their defining modules.

For the chosen q, its linear part is I+v*l, so symplecticity forces
l=B(v,-). Every odd-contact column u_i in (6) gives an affine
transvection T_(u_i) in the actual group. Its linear part I+u_i*l must
also preserve B, so l=B(u_i,-) and nondegeneracy forces u_i=v.
For completeness, preservation of B applied to x in ker(l) and any y
with l(y)=1 gives B(x,u_i)=0. Thus the nonzero functional B(u_i,-)
has the same kernel as l, and over F2 they are equal. This proves the
assertion without importing a group classification.

There are exactly 2 v_o columns. Hence their sum vanishes and

    rho(lambda)=1 in EVERY binary affine degree
    whose linear part preserves such a common B.        (8)

If v_o>0, W=<v>, so the residual normalization consists of b/4 copies
of the SAME double cover of A1, with squareclass

    product_(odd contacts j) (t-a_j)*(t-b_j).

The 2 v_o branch parameters are distinct. Each component has genus
v_o-1 and two points at infinity. For v_o=0 retain the b/2 split A1
components. This extends the credited degree-four quadratic-residual
parity argument to all these linear subgroups and all even dimensions.
It excludes their nontrivial preferred-longitude profiles; it does NOT
exclude these entire affine monodromy groups or all general linear parts.
Even with lambda=1, PF92's integral-null-core hypotheses are not automatic;
use the full PF101 peripheral calculation unless they are established.

## 7. Comparison, checks and the remaining goal

The column-span/sum description gives a finite linear-algebra interface
to the full peripheral recursion: it decides the remaining longitude
profile and component genera without expanding a new knot diagram.
The rational linking ENTRIES still need actual linking/chain data.
The pinned authors-code comparison in PF92/PF101 is retained on its
odd-dihedral overlap, not rerun or relabeled as a binary affine cover.

The checker supplies cheap exact controls for the polynomial Gauss frame,
intrinsic Hessian and local intersection counts, framing parity, affine
translations and branch-label Riemann--Hurwitz. Written parameter proofs
above establish the all-degree statements; finite samples do not.
No new external classification or theorem import is used. The ordinary
linking/intersection identity and finite-curve-cover conventions are the
same ones retained in PF101. A polar-formula literature check was a
comparison only; the needed local calculation is proved explicitly.

Natural A/S and the remaining binary affine GLOBAL compatibility,
the shared HC4 EMPTY ORIGINAL full-sheet even-radial-zero profile, and
the full finite algebra remain OPEN. The all-original-prime condition,
shared trace/marked-energy priorities, published quartic PLANE retirement
and sealed predecessors are retained. No historical priority, independent
replay, formal verification, fresh PDF visual audit or external review
is claimed. The full research goal remains ACTIVE.
