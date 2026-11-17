# All-degree structure and barriers

15 November 2026. Current qualitative research handoff, following the user's
explicit priority for general results across degrees. This synthesis does
not replace the retained cubic proof or change its least-prefix intervals.
The elementary propositions below have all-index written proofs; the small
[checker](verify_all_degree_structure.py) checks transcription and examples.
It does not certify the universal proofs. No novelty claim is made.

## 1. Keep the quantifiers explicit

Use independent standard real Gaussians and circular coordinates
Z_i=(X_i+iY_i)/sqrt(2), W_i=(X_i-iY_i)/sqrt(2), with real spectators T_j.
Then

    E(Z^a W^b T^c)=1_(a=b) product_i(a_i!) product_j E(T_j^c_j),

where the last factors vanish for odd exponents and equal (c_j-1)!!
for even exponents, including (-1)!!=1. Exponential expectations in this
note mean formal generating series; analytic integrability is not assumed.

There are different questions:

| Quantifier | Answer or retained frontier |
|---|---|
| Fixed dimension n and degree bound D: finite test for pure nullity? | Yes, some K(n,D) exists, even for ordinary moment ideals (section 2). |
| One finite test length for all degrees, already in dimension two? | Impossible (section 3). |
| One finite test length for homogeneous ternary forms of arbitrary degree? | Impossible (section 4). |
| One coefficient C, independent of deg P, in m>C deg Q for every GMC-safe P? | Impossible already for triangular shears and linear Q (section 5). |
| A complete support-only test for vanishing for every coefficient choice? | Yes, strict angular separation (section 6). |
| Does the same test classify each fixed coefficient choice? | No. Cancellation permits both safe shears and GMC failures on balanced supports. |
| Is GMC true for arbitrary degrees and dimensions? | No; the retained three-variable quartic witness already excludes that target. |

## 2. Fixed-degree finite determination is qualitative

Fix n,D and let R=Q[c_a: |a|<=D] parameterize all polynomials P of degree
at most D. Each M_m=E(P^m) belongs to R. Set I_K=(M_1,...,M_K) and
I_infinity=(M_m:m>=1). The polynomial ring R is Noetherian, so this
ascending chain stabilizes: I_infinity=I_K for some K=K(n,D).
Equivalently, choose finitely many generators of I_infinity; each uses
only finitely many M_m, so their largest index is sufficient. This gives
ordinary ideal equality, and therefore also equality of the corresponding
zero sets after complex specialization. No GMC assumption is used.

The ring theorem is [Stacks, section 10.31, lemmas 10.31.1 and
10.31.3](https://stacks.math.columbia.edu/tag/00FM). This argument gives no
numerical bound, effective stopping test, or algorithm for recognizing the
last generator. In particular I_K=I_(K+1) alone is insufficient: an ascending
chain may pause and later grow. The retained effective cubic radical-search
argument has additional finite-target input; it is a different assertion.

## 3. Arbitrarily delayed first moment in two variables

For every integer d>=2 put P_d=Z+W^d. Expanding its mth power, a Wick
contraction survives precisely when the number r of W^d factors and the
number m-r of Z factors satisfy m-r=dr. Consequently

    E(P_d^m)=0                         if (d+1) does not divide m,
    E(P_d^((d+1)r))=((d+1)r)!/r!      for r>=0,
    E exp(s P_d)=exp(s^(d+1)).

Indeed the unique surviving term contributes binomial(m,r)*(dr)!=m!/r!.
The first nonzero positive moment is d+1. For any proposed finite K choose
d>=max(2,K): the first K moments vanish, but the (d+1)st does not. Thus no
degree-independent pure test exists, even in dimension two. Adding unused
independent coordinates gives the same obstruction in every larger dimension.
This family is not pure-null and is not a counterexample to GMC(2).

## 4. The same obstruction for homogeneous ternary forms

For d>=2 take the homogeneous degree-d polynomial

    H_d=Z^d+W^(d-1) T.

If there are a factors of the first summand and b of the second, circular
balance requires da=(d-1)b. Since d and d-1 are coprime,
a=(d-1)k, b=dk, m=(2d-1)k. Real parity additionally requires dk even.
There is one surviving occupation, so

    E(H_d^m)=binomial(m,dk)*(d(d-1)k)!*(dk-1)!!

when m=(2d-1)k and dk is even; it is zero otherwise. Every displayed
nonzero value is a positive integer. Hence the first nonzero index is

    nu(d)=2d-1    for even d,
          4d-2    for odd d.

This tends to infinity. Restricting to homogeneous forms in exactly three
variables therefore does not restore a degree-independent finite test.
These are growing-degree lower bounds, not new lower bounds for fixed cubics.

## 5. Safe all-degree shears require a degree-dependent mixed cutoff

The existing [affine-W classification and sharp shear
theorem](CUBIC_AFFINE_W_CLASSIFICATION.md), claims G3AFFINEW20261114 and
G3SHEARCUTOFF20261114, prove for arbitrary polynomial u and a!=0 that

    S=aW+u(Z)T-Zu(Z)^2/(2a)

is pure-null and GMC-safe. Its full source identity is

    E exp(sS+alpha Z+beta W+gamma T)
      =exp(alpha(as+beta)+gamma^2/2+gamma s u(as+beta)
            -beta s u(as+beta)^2/(2a)).

Specialize a=1,u(Z)=Z^k, k>=1. For S_k=W+Z^k T-Z^(2k+1)/2 this gives

    E exp(sS_k)=1,
    E(T exp(sS_k))=s^(k+1),
    E(W exp(sS_k))=-s^(2k+1)/2.

Thus the fixed degree-one multiplier Q=W has nonzero mixed moment at
m=2k+1=deg S_k. No coefficient C independent of deg P can make
m>C deg Q sufficient across all GMC-safe P. This uses a positive class,
not a family of GMC failures. The existing sharp bound m>(2k+1)deg Q for
each individual shear remains valid. This is a consequence of the retained
source theorem, not a new shear classification.

## 6. A complete support-universal criterion, in every degree and dimension

Fix a finite nonempty monomial set A in (Z_1,W_1,...,Z_r,W_r,T_1,...,T_s).
For a in A define its angular charge chi(a)=(a_Zi-a_Wi)_i, and put
P_c=sum_(a in A)c_a x^a, treating the c_a as independent coefficients.
The following are equivalent:

1. Every positive moment M_m(c) is the zero coefficient polynomial.
2. Every complex coefficient choice is pure-null.
3. Zero is outside conv{chi(a):a in A}.
4. There is an integral vector eta with eta.chi(a)>=delta>0 for all a.

For 3=>4 use strict separation of a rational finite polytope and clear
denominators. For 4=>1 every monomial of P_c^m has positive eta-charge,
whereas a Wick contraction requires charge zero. Assertions 1 and 2 are
equivalent because a polynomial vanishing at every complex point is zero.

For the converse, if zero belongs to the hull, rational linear feasibility
gives integers k_a>=0, not all zero, with sum k_a chi(a)=0. Double all
k_a to make every real spectator exponent even. This occupation contributes
to M_m(c), m=sum 2k_a, the coefficient

    m!/product_a((2k_a)!)*E(x^(sum_a 2k_a a)) > 0

of product_a c_a^(2k_a). No other occupation has the same coefficient
monomial. Thus M_m is not identically zero. This proves necessity without
assuming positivity of the specialized coefficients.

More precisely, let mu(A) be the least positive sum k_a over nonnegative
integer occupations satisfying circular balance and even real exponents.
If no such occupation exists, mu=infinity. The expansion proves that every
M_m with m<mu is identically zero and M_mu is nonzero when mu is finite.
Outside the proper hypersurface M_mu=0 in the coefficient torus, the first
nonzero moment is exactly mu. This is a generic statement, not a necessity
for each specialized polynomial: complex coefficients can cancel different
occupation monomials after specialization.

For separated support and any fixed Q let L be the least eta-charge of a
monomial in Q. Every monomial of QP^m has charge at least L+m delta. Hence

    E(QP^m)=0 whenever m>-L/delta.

In particular m>(max_i |eta_i|/delta)deg Q suffices. Spectator variables
carry charge zero. For Q=0 the conclusion is immediate. This gives an
all-degree certificate with an explicit cutoff; it does not claim that all
pure-null polynomials admit this certificate. The separation mechanism is
already present in the retained one-sided proofs and Long's theorem.

## 7. General theorems already available, and exact transfer boundaries

| Retained input | General scope to use | Boundary to retain |
|---|---|---|
| [Long's factorially weighted multitorus theorem](../../foundations/extended-geometry/FACTORIALLY_WEIGHTED_MULTITORUS_THEOREM.md), FWT1 | Arbitrary angular rank; eventual pure vanishing is equivalent to strict support separation, with eventual mixed vanishing. | One radial variable with moment n!, not a product of independent radial factorials. The manuscript states it was not yet human verified or formalized; local audit is separate assurance. |
| [Prime-separating radial extension](../../foundations/extended-geometry/PRIME_SEPARATING_RADIAL_MOMENTS.md), PSMT1 | General normalized prime-separating rational radial sequences, with arbitrary angular rank. | Still one radial variable; its exact arithmetic hypotheses are essential. |
| [Restricted higher-dimensional supports](../../foundations/extended-geometry/HIGHER_DIMENSIONAL_GMC_RESTRICTED_SUPPORTS.md) | Internal multiradial proof when the balanced radial polytope has an attained coordinatewise least vector. | One scalar minimizer or radial rank one without monotonicity is insufficient. Odd spectators require parity input. This synthesis does not upgrade the prior note's assurance or register an unrestricted theorem. |
| [Affine-W classification](CUBIC_AFFINE_W_CLASSIFICATION.md), G3AFFINEW20261114 | Arbitrary polynomial h,v,u,r in P=W h(Z)+v(Z)T^2+u(Z)T+r(Z); classifies every pure-null member and decides safe versus failing branches. | A polynomial architecture, not every ternary polynomial. It already has a degree-D finite gate 5D-3 for D>=2. |
| [Quadratic-W logarithmic periods](QUADRATIC_W_LOG_PERIODS.md), GWNONSING20261114 | Nondegenerate quadratic h(Z,T), arbitrary-degree G in P=W h+G; internal all-degree exclusion. | Independent review of the analytic period/seed bridge remains outstanding; finite exact checks do not establish it. |
| [Dvorsky GVC failure and completion obstruction](../../foundations/extended-geometry/DVORSKY_ONE_PAIR_SCHUR_OBSTRUCTION.md), DVG1/DVG4 | Explicit five-variable order-three failure; the local obstruction excludes every polynomial degree in the specified one-pair harmonic completion route. | This is not an order-two Laplacian failure in five variables. The all-dimensional SIC=>GVC=>JC implications do not preserve dimension. |
| [Two-pair SIC failure](../../foundations/extended-geometry/TWO_PAIR_IMAGE_MATHIEU_COUNTEREXAMPLE.md), SIC2C4 | Sharp two-pair image-Mathieu failure and its retained all-degree extensions. | Its nonseparable coefficient matrix cannot be substituted into a separable/rank-one GVC dictionary. Independent replay of SIC2C4 remains retained, not re-executed here. |

Selected primary bodies were revisited on 15 November: [Long's manuscript](https://raw.githubusercontent.com/octonion/mathematics/main/gmc/gmc2_stronger_arbitrary_torus.tex),
especially its one-radial hypothesis and separation/prime-isolation proof,
and [Dvorsky v1](https://arxiv.org/html/2608.07338v1), sections 0--3.
This is a scoped source refresh, not a claim of exhaustive new literature
coverage, external review, or verification of producer code.

## 8. Research priority and completion boundaries

The next work should seek a structural classification with explicit safe
and failing branches, a hypothesis that controls incomparable radial
degrees, or a valid transfer preserving the required operator order and
separability. First audit reusable universal lemmas and their exact scopes,
especially the analytic period/seed bridge. A new calculation is justified
when it decides such a lemma or unlocks a genuinely new degree, dimension,
or architecture. Further cubic cutoff reductions are secondary.

The first broad task is to determine when balanced-support coefficient
cancellation can persist to every order and when it is conjugate to a
Gaussian shear, rather than identifying generic nonvanishing again. The
support criterion closes the coefficient-universal and generic questions;
it deliberately leaves this exceptional locus open. A statement that every
balanced-support pure-null point is safe is already false by the affine-W
linear branch and Long's quartic witness.

The second broad task is to establish a reusable, independently auditable
criterion for transferring a nonzero formal Gaussian period to the exact
twisted quotient/analytic cycle used in all-degree exclusions. Check its
hypotheses at the family level before optimizing any rank or pole count.

The third is to formulate a valid lower-order transfer, or an obstruction
covering an entire proposed architecture, retaining differential order,
dimension, harmonicity and rank-one tensor hypotheses. Do not retry the
already excluded one-pair completion or the nonseparable SIC substitution.

The [current cubic numerical result](RANK_ONE_TWO_BRANCH_GLOBAL94.md) remains
general/odd/homogeneous sufficient prefix 94, with least intervals
general/odd [20,94] and homogeneous [14,94]. Both lower witnesses are
complete. Its optimality and independent analytic review remain separate.
The qualitative cubic question stays internally closed; the overall
all-degree research goal and the new structural frontier remain open.
