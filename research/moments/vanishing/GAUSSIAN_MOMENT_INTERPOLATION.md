# Gaussian and sphere moment interpolation from pure polynomial shears

Let k be an algebraically closed field of characteristic zero. Use split
coordinates Z,W,Y and the normalized Gaussian functional

    L_G(Z^a W^b Y^c)=0 unless a=b and c is even;
    L_G(Z^a W^a Y^(2c))=a!(2c-1)!!.

The normalized sphere functional L_S for q=2ZW+Y^2=1 divides the
Gaussian value of each homogeneous monomial of even spatial degree 2r
by (2r+1)!!. These are two different functionals. Each assertion below
holds for either one separately; it does not prescribe them simultaneously.
Coordinates may have complex coefficients. No positive real probability
distribution with arbitrary prescribed moments is asserted.

For d>=3 let V_d=k[Z,W,Y]_(degree<=d), let m=3d+1, and set

    F_j(P)=L(P^j), M_m(P)=(F_1(P),...,F_m(P)).                  (I1)

**Theorem MomentInterpolation.** For either functional and every integer
d>=3, M_(3d+1):V_d->k^(3d+1) is surjective on k-points. Every prescribed
list of 3d+1 scalar moments is attained by a polynomial of spatial
degree at most d in three variables.

For data in any characteristic-zero field K, the polynomial can be chosen
over a finite algebraic extension of K. This is existence, not a rational
section, a unique interpolant or a bound on the coefficient extension.

## 1. Classical mechanism: openness and contracting weights

**Lemma WeightedOpenImage.** Let V be a finite-dimensional vector space over
k, and let f_j:V->k be homogeneous polynomial functions of positive degrees
w_j, j=1,...,m. If there is v0 with f(v0)=0 and rank(df(v0))=m, then
f=(f_1,...,f_m) is surjective.

Choose m linear directions on which the differential is invertible.
The restriction to their affine span through v0 is a polynomial map
g:A^m->A^m with invertible Jacobian at v0. After inverting its Jacobian
determinant, its graph presentation is a standard smooth algebra of
relative dimension zero: the equations are g_j(u)-y_j and their u-Jacobian
is invertible. This is exactly the Jacobian presentation in
[Stacks, Definition 10.137.5 and Lemma 10.137.6](https://stacks.math.columbia.edu/tag/00T1).
A [smooth morphism is open, Stacks Lemma 29.35.10](https://stacks.math.columbia.edu/tag/056G).
Hence the image contains a Zariski open neighborhood U of 0.

For a prescribed b=(b_1,...,b_m), the morphism

    alpha_b:A^1->A^m, alpha_b(t)=(t^w_1 b_1,...,t^w_m b_m)

takes 0 to 0. Its inverse image of U is open and contains 0. Since k is
infinite, it contains some t!=0. Choose v with f(v)=alpha_b(t).
Homogeneity gives f(t^(-1)v)=b, as required. Nonempty finite-type fibers
over an algebraically closed field have k-points.

Over C this has a short analytic proof: the inverse function theorem
gives a neighborhood of 0 in the image, positive weighted dilations send
any target into it, and scalar rescaling of v recovers the original target.
The algebraic proof also covers characteristic-zero fields not embedded
in C. Openness alone would not give surjectivity without the contracting
positive weights.

This is a classical openness-and-scaling argument. The specific input
requiring calculation is a pure polynomial with an invertible moment
Jacobian on the stated original spatial degree space.

## 2. The exact Jacobian in original markers

Put

    P*=W+ZY-Z^3/2,
    phi(Z)=Z, phi(Y)=Y-Z^2, phi(W)=W+ZY-Z^3/2.

The [original shear calculation](MIXED_CUTOFF_JETS_AND_INFINITE_TAILS.md)
proves that phi preserves both functionals, that P* is pure, and that

    Phi_B(t)=L(B exp(tP*))=L(phi^(-1)(B) exp(tW)).              (I2)

Its proof is self-contained: the locally nilpotent derivation
Z(-Z partial_Y+Y partial_W) preserves q, has zero divergence, and has
zero Gaussian expectation on every polynomial. Its homogeneous degree
shift also gives sphere invariance.

For n=0,...,3d choose the marker B_n by writing n=3a+r, 0<=r<=2:

    B_(3a)=W^a, B_(3a+1)=ZW^a, B_(3a+2)=YW^a.                (I3)

Only indices <=3d are used. Thus every B_n has original degree <=d,
and the markers are distinct. The inverse shear sends B_n to a
polynomial whose unique highest weighted term contributing to its period
is (-1/2)^a Z^n. Every other term has strictly lower period degree.
One can see this directly by giving Z,Y,W weights 1,2,3 before the
substitution, and noting that an unshifted Y or W strictly lowers the
exponent of t in L(Z^u W^v Y^h exp(tW)). Consequently Phi_(B_n) has
degree n and leading coefficient

    c_n=(-1/2)^floor(n/3)                         for L_G;
    c_n=(-1/2)^floor(n/3)/(2n+1)!!               for L_S.    (I4)

The differential of the j-th actual moment is

    dF_j(P*)[B_n]=j L(B_n P*^(j-1))
                =j! [t^(j-1)] Phi_(B_n)(t).                 (I5)

In row order j=1,...,m and column order n=0,...,m-1 this matrix is
upper triangular with diagonal (n+1)!c_n. In particular its determinant is

    prod_(j=1)^(3d+1) j! * (-1/2)^(d(3d-1)/2)               for L_G, (I6)

and the same expression divided by prod_(n=0)^(3d) (2n+1)!!
for L_S. Each is nonzero in characteristic zero. Since F_j(P*)=0 for
all j>=1, Lemma WeightedOpenImage applies with weights w_j=j and proves
Theorem MomentInterpolation.

The determinant and the all-order shear argument concern original
markers. No change to an uncontrolled marker filtration is made.
The general d statement uses the same fixed cubic base phase; it does
not require a classification of pure phases of degree d.


## 2.1. Higher shears and a quadratic interpolation bound

The cubic is the first member of a polynomially Gaussian-preserving family.
For r>=1 set

    P_r=W+Z^r Y-Z^(2r+1)/2,
    delta_r=Z^r(-Z partial_Y+Y partial_W).                      (I9)

The derivation fixes Z, kills q, is locally nilpotent and has zero divergence.
Gaussian integration by parts gives L_G(delta_r A)=0. On a homogeneous
input its output has one homogeneous degree, so the same Gaussian-to-sphere
normalization gives L_S(delta_r A)=0. Thus exp(delta_r) preserves each
functional and sends W to P_r. In particular P_r is pure.

For B=Z^i Y^j W^k, the inverse shear substitutes

    Y -> Y+Z^(r+1), W -> W-Z^rY-Z^(2r+1)/2.

Its period has leading degree

    n=i+(r+1)j+(2r+1)k

and leading coefficient (-1/2)^k for Gaussian, or
(-1/2)^k/(2n+1)!! for sphere. To check uniqueness, replacing one chosen
top Y term by an unshifted Y lowers Z-degree by r+1; replacing one top
W term by -Z^rY lowers it by r+1; replacing it by W lowers the
Gaussian period degree by 2r+2. Every nonzero non-top contraction has
strictly lower period degree. Hence distinct leading degrees give the
same triangular Jacobian construction as (I5).

**Lemma MarkerWeightInterval.** If d>=2r+1, the least nonnegative integer
not expressible as i+(r+1)j+(2r+1)k with i,j,k>=0 and i+j+k<=d is

    M(r,d)=(2r+1)(d-r+1)+r.                                    (I10)

Here is a direct bounded additive-basis proof. For fixed total degree
ell=i+j+k the weights are exactly ell+r*s, 0<=s<=2ell:
write s=2k+j, take k=floor(s/2), j=s mod2 and i=ell-j-k.
For a fixed residue a modulo r, the possible ell are a,a+r,...
up to ell_a<=d. Their weight intervals, all in that residue class,
join without gaps, since the next lower endpoint ell+r is at most
the previous upper endpoint (2r+1)ell plus r. The union in this residue
class is therefore all integers from a to (2r+1)ell_a, in steps of r.
The r largest degrees ell_a are precisely d-r+1,...,d. The first missing
integer across all residue classes is their smallest upper endpoint
plus r, giving (I10). This also proves that M(r,d) itself is absent;
no assumption about a generic additive basis is used.

**Theorem HigherShearInterpolation.** For every d>=3 and
1<=r<=floor((d-1)/2), the first M(r,d) moments of either functional
can be prescribed arbitrarily over k by a polynomial of degree<=d.

Choose one original marker at every weight 0,...,M(r,d)-1 using
Lemma MarkerWeightInterval. The moment Jacobian at P_r is triangular
with nonzero diagonal, and the base phase has degree 2r+1<=d.
Lemma WeightedOpenImage proves the claim.

The bound is increasing in r throughout its allowed interval. The next
argument improves it when r>=3, by using subleading coefficients of two
additional original markers.

## 2.2. Projected Jacobians give longer interpolation

**Theorem ProjectedShearInterpolation.** Let r>=3 and d>=2r+1. For either
functional separately, the first

    M^+(r,d)=M(r,d)+2r                                      (I12)

moments can be prescribed arbitrarily over k by a polynomial of degree<=d.

All Jacobian columns will remain original markers of degree<=d. Only the
time coefficients relevant to these first moments are projected: this
does not replace the original spatial degree filtration or assert a
mixed-moment cutoff at the projected time degree.

For an original marker B=Z^a W^b Y^c, expand P_r^s using i,j,k copies
of its three terms. Nonzero Wick contributions require

    i+j+k=s, (r+1)(j+2k)=s+b-a, c+j even.

Thus every nonzero coefficient of Phi_B has exponent congruent to

    a-b+(r+1)c modulo L=2(r+1).

Its leading degree n=a+(r+1)c+(2r+1)b has this same residue. Its period
therefore has possible exponents n,n-L,n-2L,... only.

The two needed subleading coefficients can be calculated without an
inverse coordinate change. Circular contraction and the ordinary
Gaussian source formula give the formal identity

    L_G(exp(tP_r+beta W+gamma Y))
      =exp(gamma^2/2+gamma*t*(t+beta)^r
                         -t*beta*(t+beta)^(2r)/2).

Differentiating in beta and, for the second marker, once in gamma gives

    Phi_(W^b)=b! sum_(h=0)^b (-1/2)^h/h!
                   binom(2rh,b-h)t^(Lh-b),

    Phi_(YW^b)=b! sum_(h=0)^b (-1/2)^h/h!
                   binom(2rh+r,b-h)t^(Lh+r+1-b).

Terms with negative exponents have zero binomial coefficients. Taking
h=b-1 yields the Gaussian coefficients

    [t^((2r+1)b-L)] Phi_(W^b)
       =2rb(b-1)(-1/2)^(b-1),

    [t^((2r+1)b-r-1)] Phi_(YW^b)
       =rb(2b-1)(-1/2)^(b-1).                              (I13)

To transfer precisely these coefficients to the sphere, fix s and a,b,c
in the original multinomial expansion. Set H=(s+b-a)/(r+1). Every
balanced contribution has the same ordinary spatial degree

    a+b+c+s+rH=2b+2s-H+c.

The sphere coefficient is therefore the Gaussian coefficient divided
by (2b+2s-H+c+1)!!, independently of the summation indices. At each
subleading coefficient in (I13), this spatial degree is 2s+2, so the
divisor is (2s+3)!!. At a leading coefficient it is instead (2s+1)!!.
All these divisors are nonzero in characteristic zero; the support and
nonvanishing used below hold for both functionals separately.

Put M0=M(r,d) and truncate periods above T=M0+2r-1. The leading-weight
construction supplies original markers with distinct period degrees
0,...,M0-1. Append the following 2r original markers in projected degree
order:

| Original marker | Highest degree after projection to degree<=T |
|---|---:|
| YW^(d-r+2) | M0 |
| Z^(j-1)YW^(d-r+1), 1<=j<r | M0+j |
| W^(d-r+3) | M0+r |
| Z^(j-1)W^(d-r+2), 1<=j<r | M0+r+j |

The first special marker has actual leading degree M0+L, and the second
has actual leading degree M0+r+L. Both exceed T. Since L>2r, their next
possible degrees are exactly M0 and M0+r, respectively, with the
nonzero coefficients (I13); their W exponents are greater than one.
The other appended markers have their displayed ordinary leading period
degrees. All original marker degrees are at most d: the special ones
have degree d-r+3<=d, and the other two families have degree at most d.
The markers are distinct, including from those of leading weight<M0.

The projected periods have all distinct degrees 0,...,M0+2r-1 with
nonzero leading coefficients. Their coefficient matrix is triangular
with nonzero diagonal. By (I5), multiplying row s by (s+1)! gives
the actual Jacobian of the first M0+2r moments at the pure phase P_r.
Lemma WeightedOpenImage proves Theorem ProjectedShearInterpolation.
The discarded higher period coefficients are not needed by that finite
Jacobian; they have not been claimed to vanish.

For the best proved count, M(r,d) and M^+(r,d) are increasing in r on
their allowed intervals. Taking the largest available r gives

    M_d=10,13,22,27  for d=3,4,5,6, respectively;
    M_d=(d^2+6d-3)/2  for odd d>=7;
    M_d=(d^2+6d-10)/2  for even d>=8.                        (I11)

The first values at d=7,8,9 are 44,51,66. Thus any universal geometric
purity prefix for degree<=d in three variables is at least M_d, by
prescribing the first M_d-1 moments to be zero and the M_d-th to be one.
These are quadratic lower bounds. Optimality of the interpolation
lengths or of the resulting purity lower bounds is not asserted.
The bounded additive-basis and projected-Jacobian calculations are
elementary connections; no historical novelty is inferred from them.

## 3. Genuine points, normalized witnesses and reverse implications

**Corollary NormalizedPrefixWitness.** For either functional, every d>=3
admits a polynomial of degree <=d with

    F_1=...=F_(3d)=0, F_(3d+1)=1.                            (I7)

This is Theorem MomentInterpolation with target (0,...,0,1).
The same conclusion, with 3d+1 replaced by M(r,d) or M^+(r,d), follows
from Theorem HigherShearInterpolation or ProjectedShearInterpolation.
All directions in (I8) below also hold at those larger indices.
For K=Q the coefficients can be chosen in Qbar. For d=3 this is an
actual algebraic cubic with its first nine moments zero and its tenth
equal to one. It is not merely a tangent vector, a formal jet or a
finite computation of a conjecturally pure polynomial.

Let R be the polynomial coefficient ring of V_d, and let
J_(m-1)=(F_1,...,F_(m-1)) in R. The following statements are equivalent:

    (a) There is P with F_1(P)=...=F_(m-1)(P)=0, F_m(P)!=0.
    (b) There is such P with F_m(P)=1.
    (c) F_m is not in radical(J_(m-1)).                     (I8)

For (a)=>(b), choose an m-th root lambda of F_m(P)^(-1) and replace
P by lambda P. The reverse implication is immediate. The equivalence
(a)<=>(c) is precisely the
[Hilbert Nullstellensatz, Stacks tag 00FS](https://stacks.math.columbia.edu/tag/00FS):
vanishing at every common zero is radical membership. These directions
are recorded separately; surjectivity is a stronger statement than any
one normalized witness.

Thus 3d+1 is a lower bound for any universal geometric purity prefix
on V_d. It is not asserted optimal. In the cubic Gaussian case the
separate [orbit and height argument](GAUSSIAN_CUBIC_PREFIX_GEOMETRY.md)
gives a stronger lower bound, conditional only on its explicitly stated
classification input. Interpolation itself does not use that classification.

**Corollary MomentAlgebraicIndependence.** F_1,...,F_(3d+1) are
algebraically independent in R. Indeed a polynomial relation would vanish
at every point of k^(3d+1) by surjectivity and hence would be the zero
polynomial. The same proof applies to the first M(r,d) or M^+(r,d)
moments from the higher-shear theorems, including the optimized M_d. The reverse
implication is not asserted: algebraic independence alone only gives
dominance, not surjectivity.

## 4. Exact checks and the remaining research question

The standalone [rational verifier](verify_gaussian_moment_interpolation.py)
expands P_r directly, uses the monomial contraction rule, and checks the
Jacobian matrices for both functionals at d=3,...,6 with r=1,
and at (r,d)=(2,5),(2,6),(3,7). It additionally checks the projected
bases at (r,d)=(3,7),(3,8),(4,9), with 44,51,66 moments. It checks
the full bounded marker-weight sets and their first missing weights,
and compares the balanced one-dimensional contraction sum with the
literal two-dimensional multinomial expansion on separate small inputs.
It does not call the inverse-shear implementation. The complete matrices,
triangular diagonals, exact determinants, special subleading sphere
normalizations, finite pure-moment controls and corrupted-matrix rejection
are checked. The retained
[receipt](gaussian_moment_interpolation_verification.json) gives the actual
run. The proof for every d and the open-image/scaling step are written
arguments, not inferences from those finite matrices.

For cubic Gaussian moments the exact universal purity prefix remains
open. The orbit-and-height argument proves that even fifteen vanishing
Gaussian moments do not force purity. Determining the full local moment ideal at P*, and then
the global coefficient scheme, requires additional equations and
geometric information.

All proof and computation assurance here is internal. No external review,
formal verification, independent whole-proof replay, constructive coefficient
algorithm, least prefix or historical novelty is asserted.
