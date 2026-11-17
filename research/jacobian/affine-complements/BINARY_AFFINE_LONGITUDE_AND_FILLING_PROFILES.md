# Binary affine longitudes have order at most four

PF2RULED125, 16 November 2026. This is an all-degree adaptation of the
user's knot-theory route to the ACTUAL intermediate covers retained by
PF122--124. It gives necessary peripheral conditions, not a new excluded
monodromy family. The global finite algebra and the shared HC4 empty-fibre
obstruction remain OPEN. No intermediate is assumed to be a plane.

## 1. Geometric input and the longitude condition

Retain the ORIGINAL finite normalization N of a complex plane Keller map,
the WHOLE source A2, and ALL omitted primes. Let N->M->A2 be any genuine
binary affine sheet-block quotient, of degree b=2^a, a>=2. PF122--123 give
one full branch component C, smooth M, and r=b/4 disjoint critical A1s.
Each critical curve maps ISOMORPHICALLY to C^nu=A1. The image is immersive
or has one ordinary cusp; double values have two smooth branches with
arbitrary finite contact. The original map may be imprimitive and its
visible indices may be greater than two.

The polynomial curve C has one place at infinity. Let K be its oriented
knot at infinity, mu its positive meridian and lambda its preferred
Seifert longitude. Restrict the actual intermediate monodromy to this
knot exterior. Its image need not be all the global affine group Q.
Nevertheless q=rho(mu) is the retained affine transvection. Translate
its fixed hyperplane to put

    q(x)=x+l(x)v,   0!=l in V*, 0!=v in ker(l), V=F2^a.

Write H=ker(l) and S={l=1}. Every critical A1 has just one end and maps
with longitudinal degree ONE. Therefore the joint mu,lambda orbit above
each critical core has size TWO: lambda preserves EACH individual pair
{x,x+v} in S. This is stronger than centralizing q or preserving S.
It does not require lambda to fix the two letters separately.

## 2. Parameter proof of the affine normal form

Every affine g commuting with q and preserving each pair in S is uniquely

    g(x)=x+(l(x)+1)u+m(x)v,
    u in H, m in V*, m(v)=0,                            (1)

where m=0 is allowed. Conversely every such formula is invertible and
satisfies both conditions. To prove necessity write g(x)=Ax+u. Commuting
with q gives l(u)=0, Av=v and lA=l. The pair condition implies
(A-I)x+u in <v> for l(x)=1; taking differences gives
(A-I)H subset <v>. Consequently A-I+u*l has image in <v>, hence equals
v*m for a unique linear m. Av=v gives m(v)=0. Conversely the linear
part I+u*l+v*m has nilpotent deviation of cube zero, so is invertible;
the conditions follow directly.

Set q_c(x)=x+(l(x)+1)v, the complementary affine transvection. A direct
composition, valid for EVERY a and every parameter in (1), gives

    g^2=q_c^(m(u)),       g^4=1.                         (2)

Thus every actual binary affine LONGITUDE HAS ORDER DIVIDING FOUR.
If it has order four its square is exactly q_c, fixes S pointwise and
is fixed-point-free on H. This does not assert that q_c belongs to Q
unless this order-four longitude actually occurs. Generic higher-degree
longitude is NOT proved trivial. The previously obtained empty-fibre
radial involution is a different peripheral element.

## 3. Complete peripheral orbit and shear profiles

On H, g is x->x+u+m(x)v. If m|H=0 this is a translation. Otherwise:
if u in <v> it has exactly half the points fixed; if u outside <v> it
has no fixed point. Equation(2) then determines whether the latter has
only two-cycles or only four-cycles. These exhaust the possibilities:

| Longitude on the b/2 unramified sheets H | Unramified boundary tori (e,f,shear) | Total tori, including the r critical ones |
|---|---|---|
| 1^(b/2) | b/2 copies (1,1,0) | 3b/4 |
| 2^(b/4) | b/4 copies (1,2,0) | b/2 |
| 4^(b/8) | b/8 copies (1,4,0) | 3b/8 |
| 1^(b/4) 2^(b/8) | b/4 copies (1,1,0), b/8 copies (1,2,0) | 5b/8 |

The last two rows require a>=3. The critical tori ALWAYS have
(e,f,shear)=(2,1,s_j), s_j in {0,1}. In (1) a pair represented by x in S
has s_j=m(x), well-defined because m(v)=0. On S/<v> this is an affine
linear function. Thus the number of s_j=1 is precisely one of

    0, r, r/2.                                          (3)

The half case requires a>=3 and occurs exactly when m|H is nonzero.
In particular the four-cycle and mixed rows have exactly r/2 nonzero
critical shears. These are restrictions on the ACTUAL preferred framing.
They are not arbitrary choices of torus basis. The total sheet count
is sum(e_i*f_i)=b in every row; the unramified longitudinal degrees can
be two or four even though all critical longitudinal degrees are one.

Use the PF101 normalized stabilizer lattice

    (e_i,0), (s_i,f_i), 0<=s_i<e_i.

Its normal filling and lifted longitude are
m_i=lift(mu^e_i), L_i=lift(mu^s_i lambda^f_i). Whenever the cores are
rationally null-homologous, their rational linking matrix obeys

    F=F^T,       F*e=s.                                 (4)

This supplies a uniform input filter for the existing full-peripheral
recursions, including longitude JOINING. It does not reconstruct the
entries of F from an orbit profile. The PF111 complete-coverage integral
reconstruction concerns the ORIGINAL whole-plane normalization. It is
NOT automatically transferred to M: the intermediate critical knots
are not proved to generate its entire end homology. Retain integral
torsion defects unless the required coverage is separately established.

## 4. Degree-four intermediates have trivial preferred longitude

When a=2 the global affine group Q is S4: its irreducible linear part is
generated by transvections, so is S3. Besides its ramified pair, mu has
two fixed sheets. Their generic residual algebra along C^nu has rank two.
Its finite normal cover of A1 can branch only as follows.

At a double value of contact k, restricting the OTHER fold to either
selected branch gives z^2=t^k*unit. Odd k causes simple quadratic
monodromy at BOTH distinct normalization parameters; even k causes none.
At an ordinary cusp the cubic local factor pulls back along its
normalization to

    z^3-3t^2*z+2t^3=(z-t)^2*(z+2t).                     (5)

Its residual sheet is a section, so the cusp adds NO residual ramification
over the normalization. The fourth local sheet is etale. At every other
finite normalization parameter the residual cover is unramified.

Consequently its quadratic squareclass has an even number of finite
odd-order zeros: one pair for each odd contact. A squareclass in C(t)
with precisely these odd valuations can be represented by their
squarefree polynomial product, whose degree is even. The residual cover
is therefore unramified and split at infinity. This includes the globally
split residual algebra. The actual longitude fixes each of the two
mu-fixed letters. Its remaining possible action is 1 or mu. A preferred
longitude lies in the commutator subgroup of the knot group; its image
has even sign in S4, whereas mu is odd. Hence

    rho(lambda)=1,  b=4, BOTH cusp profiles.              (6)

This extends PF60's rank-two residual argument to actual degree-four
INTERMEDIATES, including their ONE-CUSP case. PF60's original immersive
quartic result is credited, not claimed new. The proof does not invoke
the published quartic PLANE inversion theorem on a nonplane intermediate.
In particular the degree-four one-cusp L(p,1) survivor of PF124 must
have three peripheral tori, (2,1,0),(1,1,0),(1,1,0), and F*(2,1,1)=0.
This is a necessary condition, not an exclusion of that survivor.

## 5. Exact controls, comparison and remaining scope

The checker verifies formula(1), the square and the four orbit profiles
on small explicit affine spaces, including an order-four countercontrol.
Its written parameter proof, not those finite controls, establishes
the all-degree theorem. It also checks the cusp factorization and the
quadratic parity argument. No large representation search is needed.

The user's dihedral comparison remains the pinned, unmodified authors'
code comparison in PF92/PF101. Cahn et al., *Dihedral Linking Invariants*,
https://arxiv.org/html/2112.14790v2, treats odd-dihedral covers; the authors'
https://github.com/patriciacahn/DihedralLinkingInvariants also includes
pseudo-branch linking for signature calculations. Neither those inputs
nor a dihedral representation are substituted for the actual binary
affine cover. No new run of that sealed comparison is asserted here.

The retained shared-pencil countermodels and radial finite-lattice result
do not exclude the full-sheet empty-fibre profile. Natural A/S actions,
binary affine global compatibility, the specific missing-point obstruction,
and the full finite algebra remain OPEN. No new family exclusion,
historical priority, independent replay, formal verification, fresh PDF
visual audit or external review is claimed.
