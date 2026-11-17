# Three linear factors, the nodal conic gate, and correlated probes

Current status,15 November2026: the [full ternary cubic theorem](TERNARY_CUBIC_GMC_CLOSURE.md)
closes the degree<=3 GMC decision as an internal written proof. Any
open cubic frontiers below are historical; the earlier results retain
their stated scopes. Independent review and numerical global pure-only
gates remain separate.


14 November 2026. Exact algebra and written all-order Wick proofs. No
external review, proof-assistant verification or priority is claimed.
The global arbitrary-cubic Gaussian question remains open. The subsequent
[homogeneous nodal Newton theorem](NODAL_NEWTON_EXCLUSION.md) closes the
homogeneous nodal subclass and completes all homogeneous ternary cubics;
the general nonhomogeneous conic growth problem below remains open. The
subsequent [Gaussian-phase critical profiles](NODAL_GAUSSIAN_CRITICAL_PROFILES.md)
locate all seven critical points for irreducible leading form, exhaust all
escaping branches and prove stable width on the compact conic chart; they
do not themselves supply a uniform global period splitting or center-line
decay. The subsequent [center-line theorem](NODAL_CENTER_LINE_DECAY.md)
proves normalized center-line cut decay and excludes every zero-critical-line
member and the larger closed-half-plane critical-value class.

## Complete classification for three centered linear factors

Let X be a standard real Gaussian vector in any dimension, and let A,B,C
be complex linear forms with no constants. Put P=A*B*C and use ordinary
complex bilinear covariances

    S_ij=E(L_i*L_j), (L_1,L_2,L_3)=(A,B,C).

These are not Hermitian covariances. Make a graph on the three labeled
factors with an edge i--j exactly when S_ij!=0, i!=j. Repeated or zero
factors are permitted. After a permutation, the following table is complete.

| Nonzero off-diagonal edges | Pure-null condition | Sufficient prefix to force the condition |
|---|---|---:|
| None | At least one of S_AA,S_BB,S_CC is zero | 2 |
| Only A--B | S_CC=0 | 4 |
| A--B and A--C | S_BB=S_CC=0 | 6 |
| All three | Impossible | 8 |

Here pure-null means E(P^m)=0 for every m>=1. Every listed null profile
satisfies, for every polynomial Q,

    E(Q*P^m)=0 whenever m>deg Q.                    (1)

Thus GMC holds for every product of three centered linear Gaussian forms,
in any dimension. Vanishing of the first eight pure moments suffices for
the classification and (1). Both the universal prefix eight and the
uniform mixed coefficient one are sharp. The four graph-specific prefixes
2,4,6,8 are sharp as uniform bounds over realizable covariance matrices;
three real Gaussian coordinates suffice for their sharpness examples.
No affine constants in A,B,C, no added lower terms, and no sums of such
products are covered by this statement.

In particular a homogeneous ternary cubic which factors completely into
linear forms is GMC-safe. If its three factors are linearly independent,
its Gram matrix is nonsingular, so it is not pure-null at all: some pure
moment among the first eight is nonzero. This excludes the homogeneous
projective triangle locus, including frames unavailable to the previous
W-linear theorem. A cubic with the same leading triangle and arbitrary
lower terms is still outside this result.

## Exact covariance proof, including every zero-edge face

All odd moments vanish because P is homogeneous of odd degree. Wick's
rule computes the even moments as universal polynomials in the six
symmetric covariance entries. Rescaling each form multiplies P by one
nonzero scalar and preserves all the asserted vanishing properties.

With three edges, normalize all off-diagonal covariances to one. This is
possible over C by solving lambda_i*lambda_j*S_ij=1, using one nonzero
square root. Write the normalized diagonal entries as a,b,c and their
elementary symmetric functions as s=a+b+c, q=ab+ac+bc, r=abc. The second
moment is

    M2=r+2s+8.

On M2=0 the fourth moment is

    M4=-36*(7s^2+88s-20q+288).

Hence M2=M4=0 forces

    r=-2s-8, q=(7s^2+88s+288)/20.                 (2)

Substituting these values into the next two even moments gives

    M6=-1620*F(s),  M8=381024*G(s),
    F=107s^3+2196s^2+11280s+17152,
    G=337s^4+12736s^3+109968s^2+359424s+406272.    (3)

The retained [integer Bezout certificate](three_factor_bezout.json) gives
polynomials U,V with

    U*F+V*G=28,653,363,200 !=0.                   (4)

So M2,M4,M6,M8 cannot all vanish. This proof does not assume the covariance
matrix nonsingular or the diagonal entries nonzero.

With two edges normalize S_AB=S_AC=1, S_BC=0. If both leaf variances
S_BB,S_CC vanish, all moments vanish by the partner count below. Otherwise
exchange the leaves if necessary and use the remaining rescaling to put
S_BB=1. Write S_AA=a, S_CC=c. M2=0 forces c!=0 and

    a=-2*(1+c)/c.

On this locus

    M4=-36*(7c^2-6c+7),
    M6=32400*(c+1)*(2c^2-17c+2).                  (5)

A second integer Bezout identity in the same certificate has constant2140
and shows these cannot both vanish. No leaf-variance-zero chart has been
dropped: c=0 is already incompatible with M2=0 after normalizing the
other leaf variance.

With one edge normalize S_AB=1 and put S_AA=a,S_BB=b,S_CC=c. C is
uncorrelated with both A and B, so Wick's rule factors its contributions.
The first two even moments are

    M2=c*(ab+2),
    M4=3c^2*(9a^2b^2+72ab+24).

If c!=0 and M2=0, M4=-252c^2!=0. If c=0, C can pair only with an
external multiplier and all positive pure moments vanish. With no edges,
M2=S_AA*S_BB*S_CC, and its zero set is exactly the listed null profile.
This exhausts all covariance graphs and proves finite determination.

For sufficiency and the mixed bound, use Wick partners in an individual
monomial of Q of degree d. In the two-edge null profile, the 2m leaf
factors B,C in P^m have zero covariance with each other and with themselves.
They must pair with the m factors A or the at most d external factors.
Thus a nonzero contraction requires 2m<=m+d. In the one-edge and zero-edge
null profiles, a null factor orthogonal to the other factors occurs m
times and must pair entirely with Q, again requiring m<=d. This is an
all-order counting proof, not an extrapolation from the first eight checks.

The mixed coefficient one is sharp: for a circular pair E(ZW)=1, take
P=Z^2*W and Q=W^d. All pure moments vanish, but

    E(Q*P^d)=E(Z^(2d)*W^(2d))=(2d)! != 0.

## The finite Gram moment scheme: radical and characteristic boundary

In the polynomial ring k[a,b,c,d,e,f], let a,b,c be the three diagonal
covariances, d=S_AB,e=S_AC,f=S_BC, and let I8=(M2,M4,M6,M8). The exact
set classification is the radical equality

    sqrt(I8)=J=(a*f,b*e,c*d,a*b*c,d*e*f).           (6a)

The squarefree monomial ideal J has the six coordinate-prime components

    (a,d,e), (b,d,f), (c,e,f),
    (b,c,f), (a,c,e), (a,b,d).                     (6b)

These are respectively the three isolated null-row and three null-leaf
star profiles. Membership in their intersection gives precisely the five
minimal monomial generators in (6a); the checker exhausts all64 possible
squarefree supports. The graph case proof supplies the reverse radical
inclusion, while the partner count makes every moment vanish on all six
components.

The same polynomial proof works over every algebraically closed field of
characteristic outside

    {2,3,5,7,11,23,79,107}.                       (6c)

All moment polynomials are first formed over Z by Wick pairing counts.
The only nonzero integers inverted in the chart proof divide products of
the displayed Bezout constants and the explicit scalars2,3,20,36,32400,
1620,381024,252; their prime divisors are exactly the list (6c). This is
an integral certificate argument, not promotion from selected modular tests.
No assertion that the theorem fails in those excluded characteristics is
made. Probability and sharpness statements here remain over C.

The scheme is not reduced: a*f belongs to sqrt(I8) but not I8, because
it has Gram coefficient degree two whereas all nonzero moment generators
have degree at least three. Thus literal ideal containment cannot replace
the radical.

### Generic local thickness in characteristic zero

The [exact function-field calculation](three_factor_thickness.json) now
computes the generic local algebra on every component, not just its
support. Factor permutations give two types:

| Generic component | Local length | Maximal-ideal nilpotence index | Hilbert function |
|---|---:|---:|---|
| Isolated null row, e.g. (a,d,e) | 24 | 9 | 1,2,3,4,4,4,3,2,1 |
| Null leaves, e.g. (b,c,f) | 6 | 4 | 1,2,2,1 |

For the row component normalize b=c=1, f=z, where z is a transcendental
parameter. The second moment eliminates

    a=-2*(d^2+e^2+4z*d*e)/(1+2z^2).

Over Q(z), the remaining fourth and sixth moments are a homogeneous
binary regular sequence of degrees4 and6; their computed standard basis
contains pure powers d^4 and e^9 among its leading monomials. Their
quotient has24 standard monomials and largest surviving degree8. The
eighth moment reduces to zero in this ideal. Thus the complete-intersection
Hilbert function above and nilpotence index9 are exact at the generic
point, rather than a conclusion from a numeric parameter sample.

For the star component all generic a,d,e are nonzero. Rescaling the
three forms normalizes a=d=e=1, leaving transverse variables b,c,f at0.
The exact rational Groebner basis of the four moments has6 standard
monomials. Reduction of every degree-four monomial is zero, while at
least one degree-three monomial survives. The ranks of powers of the
maximal ideal are6,5,3,1,0, proving the second row. The two covariance
normalizations are invertible after a field extension at the generic
point, so these are generic local lengths of the original six-variable
scheme, not arbitrary transverse samples.

Let I_inf be the ideal of all positive pure moments. These same generic
local algebras hold for I_inf. On the row normalization, the reduced
moment M_(2m) is homogeneous of degree2m in d,e by the number of occurrences
of the A factor. All degrees at least9 vanish in the quotient, so all
later even moments vanish. On the star normalization, P^(2m) has4m leaf
factors but only2m center factors. At least m leaf-to-leaf Wick pairs
are required, so M_(2m) lies in (b,c,f)^m; the fourth power already
vanishes. Odd moments are identically zero. This proves persistence at
all orders without computing a long numerical prefix.

There are three linear components of each type. Consequently both the
finite and infinite Gram moment schemes have top-dimensional degree

    3*24+3*6=90.

These lengths, Hilbert functions and indices are asserted in characteristic
zero; (6c) certifies the radical identity, not automatically the same
thickness in every remaining characteristic. The subsequent [global moment-ideal theorem](THREE_FACTOR_GLOBAL_MOMENT_IDEAL.md)
now proves I_inf=I10 strictly larger than I8:8 detects the set but10
generates the full ideal. The [colon geometry refinement](THREE_FACTOR_COLON_GEOMETRY.md) now
completes the difference-class annihilator and exact global indices14/13.
The [constructive primary theorem](THREE_FACTOR_PRIMARY_DECOMPOSITION.md)
now supplies the full universal I10 decomposition and intrinsic origin
torsion77; constrained coefficient-pullback geometry remains open. Generic agreement
alone did not establish the global answer.

## An exact seven-zero moment prefix, and sharpness on every graph

Choose any complex root sigma of the cubic F in (3). Set q,r by (2), and
choose a,b,c to be the roots of

    z^3-sigma*z^2+q*z-r.

The symmetric matrix with these diagonal entries and all off-diagonal
entries one has determinant

    r-sigma+2=-3*(sigma+2)!=0,

because F(-2)=2520. Every nonsingular complex symmetric matrix is L*L^T
for an invertible complex matrix L; the rows of L therefore realize A,B,C
as independent forms in three standard real Gaussian variables. Their
product has

    M1=...=M7=0, M8=381024*G(sigma)!=0.             (6)

The final nonvanishing follows from (4). The construction is exact over
an algebraic extension; it uses no approximate root or floating matrix
factorization. This proves that eight cannot be replaced by seven in the
three-factor theorem. It also supplies a lower bound eight for the least uniform initial pure-
moment gate covering all ternary cubics. Such a finite gate exists
nonconstructively: the universal coefficient moment ideal is finitely
generated, and an initial segment containing those generators defines
all-order pure vanishing. No explicit global upper bound is obtained here. It does not improve the separate lower bound for smooth
leading cubics: this witness has singular triangle leading form.

For the two-edge sharpness choose c a root of 7c^2-6c+7, b=1 and
a=-2*(1+c)/c. The Gram determinant is -3*(1+c)!=0, and (5) gives five
initial zero moments and nonzero sixth. For one edge use a=-2,b=c=1:
M1=M2=M3=0 and M4=-252, with determinant -3. For no edges the identity
covariance has zero first and nonzero second moment. All are realizable
in three real Gaussian coordinates by complex linear forms.

## Exact nonisotropic-node reduction and affine critical catalogue

Write P=F3+F2+F1+F0 in three Gaussian coordinates. Suppose [v] is an
ordinary projective node of F3, with F2(v)=0 and v nonisotropic for the
Gaussian bilinear form. A complex orthogonal coordinate T along v gives

    P=T*H(x,y)+G(x,y),
    deg H<=2, deg G<=3, det Hess H2!=0.            (7)

Indeed grad F3(v)=0 kills the T^3,T^2*x,T^2*y coefficients; F2(v)=0
kills T^2; the nodal Hessian rank gives the nondegenerate quadratic H2.
Conversely (7) has this leading node on the T-axis. Constant or linear H
would have lower nodal Hessian rank and is not silently included here.

By an affine invertible change of x,y, put H=u*v+h. This change generally
is not orthogonal: retain the complete transformed Gaussian quadratic
q(u,v), including its linear and constant terms, and its constant volume
Jacobian. It cannot be treated as a new independent standard Gaussian pair.
The critical-point classification of P itself is nevertheless purely
algebraic and is unaffected by this change.

For h!=0, H=0 is the smooth irreducible conic v=-h/u, u!=0. Put

    R(u)=G(u,-h/u).

It is a Laurent polynomial with exponents between -3 and3. The gradient
equations are exactly

    H=0,  u*G_u-v*G_v=0,  T=-G_v/u,
    u*G_u-v*G_v=u*R'(u) on H=0.                  (8)

If R is nonconstant, multiplying u*R' by u^3 gives a nonzero polynomial
of degree at most six, so there are at most six affine critical points,
counted with their multiplicities. If R is constant r, primeness of H
and degree division give

    G=r+H*L, deg L<=1,
    P=H*(T+L)+r.                                 (9)

Its entire conic H=0,T=-L is critical. These are exactly the positive-
dimensional affine-critical cases for h!=0. Equation (9) has affine
factors and possibly r!=0, so the centered three-linear-factor theorem
does not close it.

For h=0, H=uv has two branches and a singular center. On v=0,u!=0,
critical points satisfy

    d/du G(u,0)=0, T=-G_v(u,0)/u.                 (10)

The derivative has degree at most two. It vanishes identically exactly
when G is constant on this axis; then G=r+v*L2, deg L2<=2, and

    P=v*(u*T+L2)+r.                              (11)

The other axis is symmetric. At the center u=v=0 there is a further
condition independent of these axis derivatives away from zero:

    G_u(0,0)=G_v(0,0)=0                          (12)

makes the entire T-axis critical. Otherwise there is no center critical
point. Thus the affine critical set is finite precisely when both axis
restrictions are nonconstant and (12) fails; in that case it has at most
four points, counted with multiplicity. This includes an empty critical
set. The example G=u^3+v^3 has both restrictions nonconstant but still
has the whole center T-axis critical. Ignoring this branch gives an
incorrect isolated-critical-point test.

## The conic period chart: exact identity, unresolved width control

Formal conditional integration in the original orthogonal T coordinate gives

    I(t)=E_xy exp(t*G+t^2*H^2/2).                 (13)

This is a coefficientwise Wick identity, not a convergent real quartic
integral assertion. On the chart u!=0, set

    s=t*(uv+h), v=(s/t-h)/u,
    du dv=(1/t)*(du/u) ds.

The resulting phase and measure, up to a constant volume factor, are

    Phi=s^2/2+t*G(u,(s/t-h)/u)-q(u,(s/t-h)/u),
    (1/t)*exp(Phi)*(du/u) ds.                    (14)

The u-Laurent support is within [-3,3] and the s-degree is at most three.
The chart complement u=0 and critical conic fiber s=th must still be
accounted for in the original rapid-decay homology. No assertion of a
complete cycle basis or zero-action decay follows from this coordinate
identity alone.

In particular the existing one-variable logarithmic lemma cannot simply
be tensored with ds. For a concrete negative control take

    Phi0=u+1/u-s^2/(2t^2).

On the product of a circle about u=0 and the descending Gaussian s-line,
its period normalized by 1/t is a nonzero constant: the s-integral is
sqrt(2*pi)*t and the circle residue is

    2*pi*i*sum_(k>=0) 1/(k!)^2 !=0.

There is an exact algebraic check of the same cancellation. For
D_s=partial_s-s/t^2,

    [s^2]=t^2,
    [-1/t+partial_t Phi0]=[-1/t+s^2/t^3]=0.

The growing Gaussian width cancels the prefactor. This control is not
asserted to arise from (7); it refutes the unrestricted two-variable
extension of the old logarithmic estimate. A proof for (14) must use
its special phase structure to control that width and all infinity nodes.

## October correlated-Gaussian revision and its exact data boundary

The latest [Assaad revision, 2608.14475v2, 12 October 2026](https://arxiv.org/html/2608.14475v2)
changes the title to *Correlated Gaussian Moments Determine Polynomials up
to Orthogonal Equivalence*. The [scoped audit](correlated_gaussian_source_audit.json)
records the actual sections read and the distinct hypotheses. Its real
orbit reconstruction varies replica covariances. Its local origin Wick
identity is unconditional for degree bound m>=3; global leading-Jacobian
rank and intrinsic period projection have additional hypotheses. It is
neither an external review of our internal proofs nor a theorem that pure
moments provide all correlated probes.

Here is an explicit obstruction to that last inference, using the already
classified cubic Gaussian shear

    P=W+Z*T-Z^3/2, E(ZW)=E(T^2)=1.

Circular translation puts Z=t and gives

    E exp(tP)=exp(-t^4/2)*E exp(t^2*T)=1,          (15)

so all positive pure moments vanish and the earlier shear theorem makes
P GMC-safe. For four replicas let the real Gaussian covariance between
copies be the star matrix S, with diagonal one, S_1j=rho for j=2,3,4,
and all leaf-to-leaf covariances zero. Direct Wick expansion gives

    E(P(X_1)*P(X_2)*P(X_3)*P(X_4))=-3*rho^3.      (16)

At rho=1/4 this is -3/64 and S is positive definite (its leading
principal minors are 1,15/16,7/8,13/16). Thus all pure moments zero does
not force even these positive-definite correlated probes to vanish.
The example uses complex coefficients as GMC does; it does not contradict
the paper's real-coefficient orbit separation theorem. Perfect correlation
instead gives the ordinary fourth moment zero. Extra covariance data
must be supplied or derived, not inferred from the pure series.

## Reproduction and the surviving global gate

Run `.venv/bin/python research/moments/vanishing/verify_nodal_product_reduction.py`.
The [receipt](nodal_product_verification.json) compares recursive Wick
moments with independent pair-count enumeration in all six covariance
parameters through order eight, verifies both integer Bezout identities,
all graph faces, the radical coordinate decomposition and certified characteristic boundary, the exact generic local lengths and nilpotence indices, exact sharpness formulas, the universal conic tangent
identity and transformed phase support, four exact gradient-ideal examples,
the Gaussian-width negative control and the four-replica obstruction.
The all-order statements use the adjacent written partner counts and
polynomial division; no finite check is promoted by extrapolation.

The [cubic problem](../../MATH_STATUS.json) remains open. The subsequent [full reducible theorem](REDUCIBLE_CUBIC_CLASSIFICATION.md)
closes every homogeneous reducible ternary cubic, including tangent and
transverse line-conic forms. A nonhomogeneous polynomial with the same
leading factors is not settled. The next nodal work must
separate the finite-affine-critical cases of (8)/(10), the product loci
(9)/(11), and the center line (12). The latter includes homogeneous
irreducible nodal forms. An analytic exclusion needs uniform control at
all infinity nodes and the omitted conic chart; a generic isolated/Morse
calculation cannot close the exceptional loci. The subsequent
[sharp cusp theorem](CUSPIDAL_CUBIC_SHARP_THRESHOLD.md) excludes the entire
homogeneous irreducible cuspidal class by order10. Nonhomogeneous cuspidal,
tangent and repeated leading components beyond the stated subclasses
also remain. The preceding [smooth-leading exclusion](COERCIVE_GRADIENT_MOMENTS.md)
and [complete W-linear theorem](RANK_ONE_WT2_CLASSIFICATION.md) retain their
scoped bounds; none of8,188,888 is an arbitrary-cubic upper bound.
