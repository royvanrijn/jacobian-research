# Radial moment laws and a spherical covering obstruction

16 November 2026. Internal written all-degree proofs, without novelty,
independent replay, specialist review or formal verification claims. See the
[checker](verify_radial_covering_controls.py),
[receipt](radial_covering_verification.json) and
[source audit](radial_covering_source_audit.json).

The full MZ questions on SL_2/2O and SL_2/2I remain OPEN. The theorem here
excludes an entire moment-law realization mechanism, including arbitrary
nonlinear regular coordinate choices of every degree. It does not exclude
all four-term circuits or all phase-invariant moment tables.

## 1. The newly inspected primary construction and its boundary

Zongmin Liu's [Coherent-State Circuit Cancellation, version2.2](https://zenodo.org/records/21689760)
was previously known here only through its public record. Its complete
14-page PDF is now retained and read. Sections1-3 and7 are the relevant
construction; rendered pages4,5,11 were inspected to preserve conjugations
that plain text extraction loses. In particular Section7 uses z_1 bar(z_3)
and z_2 bar(z_3), not z_1 z_3 and z_2 z_3.

The paper's factorial radial specialization assumes ACTUAL conjugate
coordinates A,bar(A),B,bar(B) with

    L(A^r bar(A)^s B^u bar(B)^v)
      = delta_(r,s) delta_(u,v) rho_(r+u) r! u!.             (1)

For F=(a+d bar(B))(bB+c A bar(A)), binomial balance gives

    L(F^m)=rho_m m! (db+ac)^m.                              (2)

If db+ac=0 and abcd!=0, the paper's fixed markers satisfy

    L(bar(B)^k F^m)
      =rho_m m! a^k d^(m-k) b^m binomial(m-1,k-1),          (3)

for 1<=k<=m; the output is zero for k>m. The full proofs of these
identities were read. Their argument is a finite binomial sum at arbitrary
m, not an inference from the ancillary checks. The raw source of the
radial checker is retained and inspected, but its larger parameter grid
is not rerun. Small local controls check transcription only.

The coherent-state theorem supplies (1) on a higher-rank group using a
highest-weight vector with an adjacent unit root string. The projective
construction supplies another (1) on a sphere using three complex coordinate
planes. Neither construction asserts descent to arbitrary noncentral finite
right quotients of SU(2). Such a quotient is a homogeneous space, not a
compact-group quotient. Its highest-weight line need not descend.

The paper's connected compact-group classification overlaps
[LONGCOMPACT20261114](COMPACT_GROUPS.md). Its virtually toral remainder,
left open in Section8, is already closed by the retained
[component-trace theorem](DISCONNECTED_COMPACT_GROUPS.md). No classification
or old frontier is reopened. The paper's sparse Gaussian minimality
certificate is not replayed or imported here. The source metadata says
version2.2, publication20July2026, created29July and modified30July; these
different fields are preserved rather than conflated into one date.

## 2. Full factorial moments force unitary radial symmetry

Let M be a compact connected smooth d-dimensional manifold without boundary,
carrying a probability measure mu with continuous strictly positive density
in smooth charts. Let q>=2 and Z=(Z_1,...,Z_q):M->C^q be C^1. Suppose,
for EVERY multiindex alpha,beta,

    integral_M Z^alpha bar(Z)^beta dmu
       =delta_(alpha,beta) rho_(|alpha|) alpha!,            (4)

where rho_0=1. The functions bar(Z_i) in this statement are actual complex
conjugates. The parameters rho_m are initially arbitrary. Positivity and
compactness are properties of the actual pushforward nu=Z_*mu, not extra
assumptions on an abstract table.

For complex vectors u,v, binomial expansion of (4) gives

    integral (u.Z)^m overline((v.Z)^l) dmu
       =delta_(m,l) rho_m m! (u.bar(v))^m.                  (5)

Consequently every unitary transformation of Z has the same polynomial
moments as Z. Polarization of (5) covers all mixed monomials. Both measures
have compact support, and real polynomial approximation on their compact
union makes equality of all these moments equality of measures. Thus nu
is U(q)-invariant. This step would be invalid if only holomorphic moments,
or only the pure circuit equations, were known.

The support of nu is exactly Z(M): the image is compact, and the inverse
image of a neighborhood of any image point is a nonempty open set with
positive mu-mass. Hence Z(M) itself is U(q)-invariant.

If d<=2q-1 and ||Z|| is nonconstant, its range contains a nondegenerate
interval of positive radii. Unitary invariance then makes Z(M) contain an
open annulus in R^(2q). This is impossible: the image of a compact C^1
d-manifold in R^(2q), with d<2q, has zero 2q-dimensional volume. For
completeness cover finitely many compact coordinate charts by O(epsilon^-d)
cubes; their Lipschitz images have diameters O(epsilon), giving total
2q-volume O(epsilon^(2q-d)), which tends to zero.

Therefore ||Z||=R is constant. If Z is not identically zero, R>0 and
the image is the full sphere S_R^(2q-1). If d<2q-1 the same covering
estimate in sphere charts rules this out too. At d=2q-1, nu is the unique
unitary-invariant probability on that sphere, namely normalized area.
Also, summing (4) over |alpha|=m yields

    integral ||Z||^(2m) dmu =rho_m (q)_m,
    rho_m=R^(2m)/(q)_m,  (q)_m=(q+m-1)!/(q-1)!.             (6)

So at the critical dimension the radial parameter sequence is forced to be
a single-sphere sequence; it cannot be an arbitrary mixture of radii.

## 3. Smooth positive-volume preservation in equal dimensions forces a cover

Here is the needed local fact with its proof. Let F:M^d->N^d be a C^1
map between compact connected smooth manifolds without boundary. Give M
a continuous strictly positive density and N a smooth strictly positive
density, both normalized to probability. If F pushes the first measure
to the second, F is a finite covering map.

Surjectivity follows from compactness and full support. Suppose DF has
rank r<d at p. In smooth coordinates, for a radius-epsilon ball about p,

    F(p+h)=F(p)+Lh+o(epsilon), uniformly for ||h||<=epsilon.

After fixed orthogonal coordinate changes its image is in a rectangular
box with r side lengths O(epsilon) and d-r side lengths o(epsilon).
Its target measure is o(epsilon^d), because the target density is locally
bounded above. The source ball has measure at least c epsilon^d.
But pushforward equality implies

    measure_M(ball) <= measure_N(F(ball)),                  (7)

contradiction. Closed small balls have compact measurable images, so there
is no measurability issue in (7). Thus DF is invertible everywhere. The
inverse-function theorem makes F a local diffeomorphism; compactness makes
it proper, and a proper surjective local diffeomorphism is a finite cover.
This argument uses actual positive measures; cancellation cannot repair a
critical point. It is not a statement about arbitrary signed or holomorphic
period functionals, or about maps between unequal dimensions.

**GRADIALCRITDIM20261116.** Under (4), a nonzero C^1 coordinate tuple
requires d>=2q-1. At d=2q-1 it has constant positive norm, (6) holds, and

    Z:M->S_R^(2q-1)

is a finite covering map. Since q>=2 the target sphere is simply connected,
so the cover from connected M has one sheet. In particular M is
C^1-diffeomorphic to a sphere and Z is a measure-preserving diffeomorphism.
The q>=2 condition is essential to this last assertion: a circle permits
nontrivial finite covers. No analogous conclusion is asserted for
d>2q-1, where the paper's higher-dimensional constructions operate.

## 4. Consequence for every nontrivial finite SU(2) quotient

**GFINITESU2RADIALNO20261116.** Let Gamma be ANY nontrivial finite subgroup
of SU(2), and let M=SU(2)/Gamma with normalized Haar quotient measure.
There are no nonzero C^1 complex functions A,B on M satisfying the ENTIRE
law (1), for any sequence rho_m. In particular this is impossible for
regular functions on SL_2/Gamma, in every matrix-entry degree.

Indeed the right action is free, so M is a compact connected smooth
three-manifold with positive smooth Haar density. Its universal cover is
SU(2)=S^3 and pi_1(M)=Gamma. Apply section3 with q=2. A nonzero pair would
make M diffeomorphic to S^3, contradicting Gamma!=1.

This excludes every two-coordinate factorial-radial realization of the
paper's cancellation mechanism on BOTH actual exceptional quotients,
Gamma=2O and2I. It also excludes it on central and cyclic quotients; it
uses no degree restriction, no equivariance, no selected generators, and
no hypothesis that the proposed coordinates are linear or holomorphic.
It is stronger than checking descent of one highest-weight line.

The trivial quotient is an essential control. On SU(2)=S^3 the two
fundamental complex coordinates have (1) with rho_m=1/(m+1)!, and the
radial map is the identity sphere parametrization. The paper's SO(3)
rank-one circuit is also consistent: it uses a different three-real-
coordinate sphere law, not (1) for two complex coordinates on SO(3).

There is no conflict with the retained algebraic
[Gaussian radial bridge](GAUSSIAN_RADIAL_TRANSFER_AND_EXCEPTIONAL_UNIFORMITY.md).
That bridge homogenizes arbitrary quotient phases in FOUR ambient Gaussian
variables and uses a radial factorial factor plus finite averaging. It does
not construct two regular complex coordinate functions on the quotient
with the full law (1). Its ordinary ideals, uniformity and all protected
records remain unchanged.

## 5. Exact limitation and next major gate

Section1 of the paper permits arbitrary phase-invariant tables mu_(r,u).
For its circuit, pure nullity requires only

    sum_(j=0)^m (-1)^j binomial(m,j)^2 mu_(j,m-j)=0.         (8)

Equations (8) alone do NOT give the full isotropic moment identity (1),
and the preceding proof does not exclude an actual quotient realization
of a more general such table. Nor does pure nullity of a general phase
imply that it comes from a two-coordinate circuit. Thus full exceptional
MZ remains OPEN, with all non-torus-invariant, mixed-isotypic and
image-interior cases retained. The actual complement-logarithm obstruction
is still the general target. An arbitrary measurable map, abstract
moment functional, or finite matching prefix does not meet (4).

The checker confirms low exact unitary-moment transformations, the sphere
norm sequence and a projective law with nonconstant radius moments, the
paper's small pure/mixed circuit coefficients, and a local critical-map
mass scaling. These are transcription and boundary controls only. The
universal dimensional and covering conclusions are the written proofs of
sections2-4; no finite sample proves them. No large source grid, new sparse
minimum, general period producer or external theorem solving either MZ
case is claimed. External review/outreach remains deferred.
