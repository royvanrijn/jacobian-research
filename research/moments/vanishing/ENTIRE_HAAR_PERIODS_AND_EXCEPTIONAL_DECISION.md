# Entire Haar periods and effective exceptional-model tests

15 November 2026. Internal written deductions connecting the repository's
[holonomic source theorem](HOLONOMIC_SOURCES_AND_CUTOFF_LOCI.md),
[finite-family decision mechanism](ALGEBRAIC_SPECIALIZATION_AND_FINITE_FAMILY_DECISION.md)
and [Haar shears](HAAR_SHEARS_AND_GENERAL_NULLCONE_CLASSIFICATION.md).
The universal algorithms are written, not implemented or executed.
Independent specialist review and formal verification remain outstanding.
No novelty claim is made. The [source audit](haar_period_source_audit.json)
and [finite checker](verify_haar_period_decision.py) delimit the evidence.

The global Mathieu--Zhao questions for SL_2/2O and SL_2/2I remain open.
This note supplies a structural criterion for failure, a terminating
**theoretical** decision in each fixed finite coefficient family, and a
complete classification of a concrete exceptional family. It does not
report a computation deciding either whole model or its first full degree.

**Later uniform transfer.** The [Gaussian radial bridge](GAUSSIAN_RADIAL_TRANSFER_AND_EXCEPTIONAL_UNIFORMITY.md)
now supplies a common linear matrix-entry cutoff for all safe pure-null
phases at each fixed degree on both exceptional quotients, and proves their
pure-safe loci closed. Determinant homogenization, radial integration and
finite averaging preserve every ordinary cutoff coefficient ideal. Gaussian
nilpotent lifting and formal-safe completion therefore apply as well. These
are internal consequences of the retained Gaussian theorem, not computations
of the common constants or decisions of either global MZ question.

**Later executed result.** The [octahedral lowest-type computation](OCTAHEDRAL_LOWEST_TYPE_AND_SHARP_PREFIXES.md)
now closes the full degree<=11 space on SL_2/2O,with sharp pure prefix7
and ordinary ideal prefix10. The general source solver remains unimplemented;
the icosahedral first space and both all-degree model questions remain open.

## 1. Compactness removes the Gaussian finite-pole obstruction

Let h be integration against a finite measure on a compact space, and let
A be a complex algebra of continuous functions. Assume h(1)=1. For fixed
P,Q in A put

    F_Q(t)=h(Q exp(tP))=sum_(m>=0) h(Q P^m) t^m/m!.

The bound |F_Q(t)|<=||h|| ||Q||_infinity exp(|t| ||P||_infinity)
shows that F_Q is entire. Therefore F_Q is rational if and only if it is
polynomial. More strongly, it is algebraic over C(t) if and only if it
is polynomial: an algebraic equation with polynomial coefficients gives
a polynomial growth bound outside a large disk by the elementary root
bound; Cauchy estimates then make the entire function a polynomial.

Define the **actual period image**

    V_P=sum_(Q in A) C[t] F_Q(t) subset O(C).

It is a C[t]<partial_t>-module because partial_t F_Q=F_(PQ).
This is an image with the specified compact integration cycle, not the
abstract algebraic integration module itself.

Suppose P is pure-null, so F_1=1. Then C[t] is contained in V_P, and

    every fixed Q has h(Q P^m)=0 eventually
      iff V_P=C[t]
      iff dim_(C(t))(C(t) tensor_(C[t]) V_P)=1.             (1)

Polynomiality is exactly eventual coefficient vanishing. If the generic
rank is one, every F_Q is a rational multiple of 1; entireness makes it
polynomial. The converse is immediate. Finiteness of rank is not needed
for this argument; it holds in the holonomic examples below.

Moreover V_P/C[t] has no C[t]-torsion: r(t)F(t) polynomial, with nonzero
r and F entire, forces F polynomial. Every bad marker thus gives a
**transcendental entire** period and actual generic rank greater than one.
It cannot be explained by the Gaussian image K[t,(1-t)^(-1)] and its
rank-zero pole defect. This does not exclude singularities of an abstract
differential equation or assert anything about arbitrary quotients.

For any reductive homogeneous pair with its compact integration model,
(1) applies to Reynolds. It does not identify the actual period rank with
the generic rank of abstract cohomology. Haar-shear examples have actual
rank one despite being outside the original acting-group nullcone.

## 2. An explicit holonomic joint source on SL_2

Let B=K[x,y,z,w]/(xw-yz-1), over a specified effective finitely generated
characteristic-zero field K. Adjoin i if needed. Haar moments are exactly

    h(x^a y^b z^c w^d)=0                 unless a=d and b=c,
    h(x^a y^b z^b w^a)=(-1)^b a! b!/(a+b+1)!.             (2)

All polynomial moments are therefore computable in K. For any polynomial
representative P define

    F(t,s)=h exp(tP+s_x x+s_y y+s_z z+s_w w).

Write D_i=partial_(s_i), and use substitution derivations

    E=x partial_z+y partial_w,
    F_-=z partial_x+w partial_y,
    H=x partial_x+y partial_y-z partial_z-w partial_w.

The following five polynomial operators annihilate the joint source:

    D_x D_w-D_y D_z-1,
    partial_t-P(D),
    s_z D_x+s_w D_y+t (EP)(D),
    s_x D_z+s_y D_w+t (F_-P)(D),
    s_x D_x+s_y D_y-s_z D_z-s_w D_w+t (HP)(D).             (3)

The first two use det=1 and differentiation. The last three follow from
h(Ea)=h(F_-a)=h(Ha)=0; their plus signs matter. They hold analytically
for all complex t,s and coefficientwise over K. Another representative
modulo det-1 gives the same source.

These equations are a supplied holonomic system. Before Fourier
transformation in the matrix variables, start with the cyclic module
with relations det-1,E,F_-,H. On det=1 the three tangent fields span the
tangent bundle by left translation. Its characteristic variety is
contained in the conormal bundle, of dimension four. The module is
nonzero, and Bernstein's inequality makes it holonomic. External product
with O_t and twisting by exp(tP) preserve holonomicity, adding the relation
partial_t-P. Algebraic Fourier transformation x_i->D_i,
partial_(x_i)->-s_i preserves holonomicity and gives (3). The tangent
relations acquire harmless minus signs; their scalar divergences vanish,
including 1+1-1-1=0 for H.

The actual compact cycle can be written

    g=[[a+ib,c+id],[-c+id,a-ib]], a^2+b^2+c^2+d^2=1.

Its surface measure is a constant multiple of delta(r^2-1). Tangent
rotations and r^2-1 give the same conormal construction after this linear
change of variables. This does not equate SU_2 with a real determinant
hypersurface in four independent real matrix entries.

The classical inputs are holonomic Fourier transformation, integration,
and the integration b-function with a lift; see
[Oaku 2018, Theorems 5.1 and 5.4, Algorithm 5.6](https://www.lab.twcu.ac.jp/oaku/MSJSI2015oakuLecture.pdf),
printed pp.308--314. The exact cached source was already retained by the
Gaussian work. Compactness justifies coefficientwise integration by parts.

## 3. Pointwise decision, invariant markers and linear cutoffs

Let Gamma be an explicitly presented finite subgroup of SU_2, allowing
an algebraic conjugate embedding by transport, and let P in A=B^Gamma.
The marked condition on A is equivalent to that on B, since

    h(QP^m)=h((E_Gamma Q)P^m),
    E_Gamma Q(g)=|Gamma|^(-1) sum_gamma Q(g gamma).          (4)

Averaging does not increase matrix-entry degree. A bad ambient monomial
thus gives an invariant fixed marker of at most the same degree. Its
average is nonzero because its marked series is bad. No power or norm
of a noninvariant failure witness is assumed pure-null.

Apply the retained general Taylor-source theorem to (3). In principle it
computes a nonzero polynomial b and a finite normal-order lift

    (b(E_s)+R)F=0, E_s=sum_i s_i D_i,
    R=sum_l c_l t^(k_l) partial_t^(j_l) s^(u_l) D^(v_l),
    r_l=|u_l|-|v_l|>=1.                                   (5)

For ordinary Taylor coefficients F_alpha=[s^alpha]F,

    b(|alpha|) F_alpha
      =-sum_(l:alpha>=u_l) c_l
          (alpha-u_l+v_l)!/(alpha-u_l)!
          t^(k_l) partial_t^(j_l) F_(alpha-u_l+v_l).       (6)

All right-side source degrees are lower. F_0 has constant coefficient one,
so b(0)=0. Let B_0 be the largest nonnegative integer root of b. All marked
series are polynomial exactly when those with |alpha|<=B_0 are polynomial.
Each has a computable scalar annihilator; the retained origin/infinity
indicial tests decide polynomiality from a finite exact jet. The same
zero test separately decides F_0=1. This is a terminating theoretical
pointwise test for whether a supplied exactly represented P is an MZ
counterexample. It uses the full source, not the single pure annihilator.

On failure, (4) returns an invariant fixed marker and a finite differential
and jet certificate for an infinite nonzero tail. On safety, exact marked
polynomials can be computed at every requested bounded marker degree.

If P is pure-null and passes, F_0=1. Put

    C_rec=max(0,max_l (k_l-j_l)/r_l),
    C=max(C_rec,max_(1<=|alpha|<=B_0,F_alpha!=0)
                       deg_t(F_alpha)/|alpha|).          (7)

Then (6) proves deg_t F_alpha<=C|alpha| by induction, hence

    h(QP^m)=0 whenever m>C deg_matrix_entries(Q).          (8)

C is a sufficient computable rational number depending on P, not a
uniform coefficient produced by this pointwise construction. The later
radial bridge gives a separate fixed-degree uniform coefficient for all
pure-null safe points. A non-pure all-marker-safe P
similarly has an affine bound. Generic certificates must not be specialized
blindly. No general exceptional b-function computation is reported here.

## 4. Finite coefficient families are theoretically decidable

Fix Gamma and a finite polynomial coefficient family of invariant functions,
with exact algebraic/rational equations and inequations. This includes
the full invariant space of any specified matrix-entry degree d: a basis
and invariance constraints follow by reducing modulo det-1 and averaging
a finite spanning set.

The following theoretical construction terminates.

1. Compute the full ordinary pure-moment ideal by ideal-membership zero
   tests. For a proposed ideal J, compute primary components and inject each
   primary quotient into a finite algebra over a normalization field, as in
   the retained primary oracle. Formula (2) computes every coefficient.
   Multiplication by P in the finite algebra gives a commuting matrix-valued
   exponential connection; nilpotent coefficients alternatively expand into
   finitely many parameter derivatives of scalar exponential connections.
   All component series are holonomic, so the zero oracle decides whether
   every pure moment belongs to J. If not, enumeration returns one outside
   J. Strict ordinary-ideal ascent terminates, retaining nilpotents.
   Inequations use inverse-variable charts vq-1=0.

2. Take the finitely many minimal primes of the pure-null ideal. Run the
   all-marker test on each effective generic field.

3. If every generic component passes, the family has no failure. For each
   fixed Q its generic tail coefficients vanish in the component domain,
   so vanish at every point, including certificate-denominator exceptions.
   This specializes coefficient identities, not generic certificates.

4. If a generic component fails, apply the retained bounded-space algebraic
   specialization lemma to its selected bad marked series. The norm of its
   infinity indicial polynomial through a finite Noether-normalization model
   has bounded coefficient degree in the normalization variables, independent
   of the candidate integer degree. A Kronecker/Eisenstein specialization
   preserves every nonzero integer evaluation simultaneously, a nonzero
   coefficient above the generic degree bound, and all required denominators.
   No new polynomial degree is possible. Pure identities remain zero, and
   (4) yields a number-field invariant failure marker. The detailed finite
   algebra and simultaneous-specialization proofs are in the linked notes;
   no new analytic specialization assertion is used.

Thus each fixed bounded-degree model question is decidable in principle,
with an algebraic witness on failure. The later lowest-type calculation
closes the full octahedral degree<=11 space; the general solver remains
unimplemented. Any complex failure can be replaced by an algebraic
one with the same P-degree bound, and an invariant marker of at most the
selected ambient monomial degree.

Enumerating d and applying these terminating decisions semidecides failure
on each whole model. No finite list of safe degrees proves all-degree
safety; the enumeration need not terminate for a safe model.

## 5. Entire nonpolynomial control: the retained normalizer circuit

For P_N,U from the [spectral note](EXCEPTIONAL_SPECTRAL_BARRIER_AND_POWER_AVERAGES.md),

    h(P_N^m)=0, h(U P_N^m)=c_(2m)>0 (m>=1),
    c_k=4^k(k!)^2/(2k+1)!.

Let f=F_U=sum_(m>=1)c_(2m)t^m/m! and G=1+f. Then

    [t^m]G=(1/2)_m/((3/4)_m(5/4)_m).

With theta=t partial_t, put

    L=(theta-1/4)(theta+1/4)-t(theta+1/2).

The coefficient ratio proves LG=-1/16, Lf=t/2. Therefore Mf=0 for
M=(t partial_t-1) composed L, namely

    M=1/16+(-t/16-3t^2/2) partial_t
          +(2t^2-t^3) partial_t^2+t^3 partial_t^3.         (9)

The infinity indicial polynomial is -v(v+1/2), bounding polynomial
solutions by degree zero. The origin indicial polynomial is
(v-1)(v-1/4)(v+1/4). The nonzero first coefficient 8/15 is a finite
certificate that f is nonpolynomial. Compactness makes it entire and
transcendental. Its actual period image contains 1 and f and has generic
rank at least two, rather than a rank-zero Gaussian pole defect.
This is not an exceptional invariant marker; the failed averaging attempts
remain as recorded in the spectral note.

## 6. A complete two-axis family on both exceptional models

Use the real invariant binary forms from the Haar-shear note, with n=8
for 2O and n=12 for 2I. Write p=f_n(x,y), q=f_n(z,w),
Phi_t=exp(tqF_-), P_t=Phi_t(p). The full three-parameter family is

    R_(a,b,t)=a P_t+b q=Phi_t(a p+b q).                    (10)

For every k>=1, moment preservation and weight balance give

    h(R^(2k))=binom(2k,k) C_k(ab)^k,
    h(R^(2k+1))=0, C_k=h((pq)^k)>0.                       (11)

Consequently its full ordinary pure-moment ideal is

    I_infinity=(ab) in C[a,b,t].                          (12)

The second moment generates it; the first is identically zero. Pure-nullity
is exactly ab=0, and every such member is eventually marked-null by the
automorphism transfer. The q-axis has sharp coefficient 1/n. The P_t-axis
has coefficient (n+1)/n for every t, sharp for t!=0; at t=0 it drops to
1/n. The zero function has coefficient zero. Sharpness follows from

    h(p^s P_t^((n+1)s))=t^(ns) C_((n+1)s), deg(p^s)=ns.

Changing s changes the marker and does not create an infinite fixed tail.
The second moments are 16ab/15 for 2O and 25ab/546 for 2I. At b=0,
a!=0,t!=0, the family retains the marked-safe non-nullcone points of the
Haar-shear theorem. This nodal pure-null parameter scheme contains a
whole non-nullcone safe branch, and no exceptional MZ failure.

## 7. Next gate

Seek actual period-rank excess on a pure-null component of a right-invariant
coefficient family. A bad marker can be certified by one scalar equation
and a finite exact jet. Abstract cohomology rank and non-nullcone separators
do not suffice. Family (10) is now closed and should not be reopened by
varying its three parameters. The full octahedral space through degree11
is now closed by the linked lowest-type computation; the first icosahedral
space, larger mixed-type families and all-degree MZ classification remain
open. The radial bridge closes fixed-degree safe-locus uniformity, while
structural classification of that locus remains a separate question.
The finite-family construction is a justified route for future
implementation, not evidence of its numerical outcome.
