# Complementary corner plumbing and the boundary JSJ endpoint

15 November 2026. PF2RULED77 gives a degree-free source-end reduction
for an exact-volume smooth plane extension with TWO entire omitted
lines, one fold R and one unramified D, in complementary source pencils.
The R-knot exterior has a path of Seifert pieces with its marked torus
boundary at an END. Its canonical boundary JSJ piece consequently has
at most ONE adjacent JSJ torus.

Every non-trefoil S3/S3 actual-binomial two-stage QUARTIC target has
TWO adjacent JSJ tori at that boundary. Stable full peripheral systems,
with Waldhausen's hypotheses checked, distinguish these exteriors.
Thus the ENTIRE non-trefoil S3/S3 representation sector is excluded
in complementary sources as well as the common sources closed by PF76.
This does not exclude other representations of the same target keys.
Trefoil S3/S3 complementary sources and the other passing complementary
sectors remain open. No full quartic or higher-degree exclusion is claimed.

## 1. Hypotheses and the actual complementary corner

Let N be smooth affine over C with marked U=A2_(x,f),
N minus U=R disjoint-union D, both ENTIRE affine lines, and

    div_N(dx wedge df)=R, with dx wedge df algebraically exact.

Assume the intrinsic pencils are complementary, as in PF50:
f is regular and constant at a finite value on R, and x is regular
and constant at a finite value on D. The other coordinate has a pole
there. Their WHOLE charts, and the affine component-open statements
needed for PF43/46, are retained hypotheses. There is no field-degree,
Laurent-support or pole-order bound in the source theorem. It is not
a claim for arbitrary higher-degree normalizations or any number of
omitted lines. The actual complementary quartic packet satisfies them.

Choose the complementary coordinate pair, shift f|R to0, and write
rho=v_R(f), q=-v_R(x). PF50's exact-volume and negative-definite
pair inequality give the COMPLETE ramified alternatives

    rho=1: q=2, x=c*f^-2+A+f*T;
    rho>=2: q=rho-1, f=z^rho,
             x=c*z^-q+z^(2-rho)*T.                         (1)

c is nonzero. Polynomial terms can be absorbed into the whole chart;
for rho=1 the simultaneous constant A is removed by translating x,
which only translates the finite value x|D. No noncomplementary
coordinate change or deletion of D from the finite map is used.

Both coordinates extend as morphisms N->P1 by their entire charts.
Resolve (x,f):N->P1 times P1 outside N, as in PF67. Prune boundary
blowups not on R or D ancestry in reverse order; their contractions
remain outside N. R lies over a finite point of H_x={x=infinity},
D over a finite point of H_f={f=infinity}, and neither lies over
the intersection H_x meet H_f. The original corner edge is untouched.
Cut its plumbing torus, keeping the whole boundary. This separates
the R-side B_R from the D-side Y_D. The end is their ACTUAL gluing.

Let h_x,h_f be the positive normal circle fibres of the original
sections. On the R-side boundary write d=h_f for the base curve,
using the ordinary positive-edge plumbing exchange. The pair(h_x,d)
is a fixed integral boundary basis. On the D side the fibre is h_f
and its corresponding base curve is h_x. These actual framings
must be kept when converting a rational filling slope.

## 2. Minimal higher-order R: a solid torus and a unit cable

For rho>=2 set q=rho-1 and n=q^2+q+1=rho^2-rho+1. The local
curvette (u=1/x,f) has orders(q,q+1). Its exact Euclidean extraction
has first multiplicity q, followed by q+1 unit multiplicities. After
the first blowup the coordinate f/u=c*z+z^2*T has order1. There is
one free b=1 creation, q-1 crossings with the old b=0 component,
and the final free b=2 creation R. This proves the pattern for
UNRESTRICTED q; it is not inferred from sampled blowup matrices.

Remove R and retain the boundary puncture toward H_f. The R-side
weighted graph, in order from that puncture, is

    [-1, -(q+1), -2 repeated q times].                     (2)

The omitted R arrow is on the FIRST -2 vertex, the third vertex of
this chain. All off-path boundary blowups have already been undone.
A linear plumbing with one end torus is a solid torus, by successive
model solid-torus/annulus absorptions. The relevant continued fractions
have nonzero denominators. In particular

    [q+1,2,...,2]=n/(q+1),
    [1,q+1,2,...,2]=q^2/n.                               (3)

This is a geometric linear-plumbing argument, not a solid-torus
recognition from cyclic homology.

Let g be the far-end fibre/core generator. Successive chain relations
give h_x=n*g, d=q^2*g. The primitive solid-torus meridian and a fixed
longitude are therefore

    mu=q^2*h_x-n*d,
    ell=(1-q)*h_x+q*d, [ell]=g,
    det((mu,ell) in (h_x,d))=1.                           (4)

The R arrow is the fibre P of that first -2 vertex. The first two
punctured chain vertices are torus products, so transport to the outer
torus gives an ACTUAL curve, not only its class in the solid-torus group:

    P=q*h_x-(q+1)*d=mu+q*ell.                             (5)

Thus the R knot is the unit-meridian cable of the filling core with
winding q. For q=1 it is the core and drilling it leaves T2 times I;
for q>=2 its exterior in B_R is an annulus Seifert piece with one
genuine cone q. On the D-side boundary the filling meridian is

    q^2*h_x-n*h_f.                                       (6)

The entire higher-rho source end, with its R knot retained, is now an
explicit filling of Y_D with this marked solid-torus pattern. We make
no assertion that an arbitrary proposed Y_D/filling is affine or finite.

## 3. The simple R exception is also a boundary-end piece

For rho=1,q=2 the exact chart in (1) has five unit centres. The
R-side after removing R has all weights-2, with edges

    H_x -- E_2 -- E_3 -- E_4 -- arrow R,
             |
             E_1,

and the remaining torus puncture at H_x toward H_f. E_1 is a leaf.
This follows directly from the two initial section blowups, then
three free blowups; the last creates R. B_R is Seifert over a disk
with cone orders2,3. The R knot is the core of the order3 leg.
Drilling it therefore gives an annulus piece with one genuine cone2,
between the R boundary and Y_D. No oriented named-knot identification
or solid-torus assumption is made for this exceptional B_R.

## 4. Arbitrary D support still gives a rooted path

Write the unramified entire chart, allowing every characteristic gcd
chain and every coefficient permitted by exactness,

    x-lambda=w^sigma,
    f=h(w)+w^(1-sigma)*V,
    h=sum_(j>=sigma+1)^p c_j*w^-j,
    gcd(sigma,support h)=1, c_p nonzero.                  (7)

For sigma=1, exactness removes c_1 and p>=2. At the point on H_f
use u=1/f and x-lambda. The reference branch (u=1/h(w),w^sigma)
is a SINGLE irreducible plane branch. Its characteristic resolution
has a linear chain of rupture vertices: each Euclidean block has
only one strict-branch continuation, with its other exceptional
directions terminating in linear legs. Subsequent characteristic
blocks continue that same chain. They do not create two continuations.

The original smooth axis u=0 separates at the initial integral or
first nonintegral contact, BEFORE any later characteristic block.
If p is divisible by sigma, it separates before the first nonintegral
exponent; otherwise it separates in the first characteristic block.
Thus the boundary puncture at H_f lies at the FIRST end of the
rupture chain, possibly behind a linear collar. This is why the root
location is known, not just the unrooted shape of a plane-branch graph.

As in PF43, let j_* be the last support index decreasing the characteristic
gcd and beta_*=2p-j_*. The first chart-varying term in u occurs at
2p+1-sigma. Hence the number of final free contacts is

    kappa_D=j_*+1-sigma>=2.                               (8)

This replaces PF43's index-two free-tail count by the actual index-one
count. Removing the last D leaf leaves a linear negative-two tail;
it becomes the final exceptional order kappa_D rather than a second
rupture continuation. With sigma>=2 the resulting Y_D is a path of
annulus pieces with genuine cones, terminated by a disk piece with
two genuine cones. If sigma=1 there are no characteristic blocks:
the same direct integral-contact extraction gives a single disk piece
with cone orders p,p, p>=2. The corner boundary is at the first end
in either case. No length or support bound has been used.

All these pieces are irreducible with incompressible torus boundary.
Annulus pieces with genuine cones have negative orbifold Euler number.
The terminal disk piece has two cone orders at least2; the (2,2)
case is the twisted I-bundle over the Klein bottle and is INCLUDED,
not incorrectly treated as having a unique fibre. Its boundary is
incompressible and it is irreducible as well. See Hatcher
Propositions1.12-1.13 and the explicitly recorded exceptions in
Lemmas1.14-1.15. Product collars can be removed.

## 5. The source boundary JSJ invariant

Glue Section2's annulus-cone-q piece, its q=1 collar, or Section3's
annulus-cone-2 piece to Y_D. This is an ACTUAL path decomposition
of the R-knot exterior into irreducible boundary-incompressible
Seifert pieces. The R-boundary component is at the first end.
The whole exterior is irreducible and boundary-irreducible, by the
innermost-circle sphere/disk argument and boundary injection.

Some adjacent fibrations may agree and merge. Product collars or
the terminal Klein-bottle alternative can affect which tori remain;
no unique-fibre assertion is used for that exception. Starting from
this Seifert decomposition, choose an inclusion-minimal subcollection
of its incompressible tori giving Seifert or atoroidal pieces.
Hatcher Theorem1.9 identifies that collection as the canonical JSJ
collection up to isotopy. Deleting or merging edges of a path leaves
a path with the marked boundary at an end. Consequently

    the JSJ piece containing the R boundary has at most
    ONE adjacent canonical JSJ torus.                     (9)

The TWO-omitted-line hypothesis is essential here. Several unramified
clusters on H_f can create extra continuations; no arbitrary-count
version of (9) is asserted. The local solid-torus calculation (2)-(6)
can still be reused when its corner hypotheses hold.

## 6. The non-trefoil S3/S3 target has the boundary in the middle

In every PF63 S3/S3 representation write e=2b,o=3t,r=3s, with
s,t positive odd, gcd(2b,3t)=1, and (b,t)!=(1,1). Keep all actual
fixed-sheet fillings as in PF73/75. Remove ONLY the simple ramification
core, leaving the other meridian-cycle fillings in place. The exterior
has three actual Seifert pieces:

    E_1: the positive T(e,o) exterior,
    P_R: base with THREE boundary circles and one cone s,
    E': base disk with cones b,2b,t.

The order1 cone s=1 is omitted, so P_R is then pair-of-pants times S1,
not T2 times I. Its R boundary is the third boundary. Both companions
have at least two genuine cones, negative orbifold Euler number and
incompressible boundary at every non-trefoil first pair (PF76).
The core drilled from the cable is a regular fibre since its true
lifted branch meridian meets h once; no extra exceptional order is
introduced at that boundary.

On the fixed edge, the fibres are h=k*M+r*L and H=eo*M+L.
On the central degree2 edge use m=M^2,L and divide those M coefficients
by2. Thus the two fibre-intersection numbers are

    r(e+o)-1 and (r(e+o)-1)/2, BOTH nonzero.                (10)

All three pieces have negative base-orbifold Euler number. None is
a product/Klein exception. Hatcher fibre uniqueness therefore prevents
merging across either companion torus. Those two tori form the minimal
JSJ decomposition by Theorem1.9. The piece P_R containing the marked
boundary consequently has EXACTLY TWO adjacent JSJ tori.
This holds for every integer lift and representation in the sector,
not only r=3 or s>=3, and keeps both fixed-sheet companions.

## 7. Peripheral recognition, not recognition from homology

The actual full N and proper embedded R have cofinal collared end pairs
from the compactification and, if a finite map exists, from target-radial
sections. Interlacing their cofinal product neighbourhoods gives the
stable exterior-group isomorphism and full torus peripheral system,
including the normal R meridian. This is PF44's retained end-pair
argument applied to the FULL two-prime N. Neither component open is
treated as a finite cover after deleting the other prime.

Both exteriors above are irreducible and boundary-irreducible. The
target exterior is aspherical (a tree of injectively glued aspherical
Seifert pieces) and sufficiently large: it contains the two genuine
incompressible companion tori. Realize the peripheral group isomorphism
by a map of pairs to that aspherical target, with a homeomorphism on
the one boundary torus. Waldhausen Theorem6.1 applies. Its exceptional
product over a CLOSED orientable surface has two boundary components,
whereas the source has one. The resulting covering has degree one
because pi1 is an isomorphism, and its homotopy preserves the boundary
map. We obtain a boundary-preserving homeomorphism of exteriors.

Such a homeomorphism preserves the JSJ piece containing that boundary.
But (9) gives at most one adjacent torus, whereas (10) gives exactly
two. Contradiction. This excludes ALL non-trefoil S3/S3 complementary
sources at all target exponents/windings and all source orders/poles.
PF76 already excludes their common sources, so the ENTIRE non-trefoil
S3/S3 representation sector is closed in BOTH source pencils.

## Evidence and remaining work

[The checker](verify_complementary_corner_jsj.py) pins the retained
source/classification/peripheral claims, audits exact solid-torus
framings and cable slope identities, verifies arbitrary-support source
cluster controls (including terminal Klein cases), and checks the
complete S3/S3 incidence and target fibre-intersection formulas.
[The receipt](complementary_corner_jsj_receipt.json) and
[manifest](sources/complementary-corner-jsj/manifest.json) retain
complete current Hatcher/Waldhausen page audits. The unrestricted
Euclidean/characteristic and recognition arguments, not the bounded
controls, prove the new claims. Imported topology is not independently
reproved; no novelty/priority or external-review claim is made.

Trefoil S3/S3 has a solid-torus main companion; that edge can disappear
and (10) does not force two boundary-adjacent JSJ tori. It is NOT
declared excluded in the complementary pencil. Other surviving
complementary sectors, other/non-binomial keys, longer sequences,
higher degree and the full goal remain OPEN. No finite Keller map,
full quartic exclusion, formal verification or independent full
reproof is asserted. PF63 polynomial immersions still exist.
