# Every additional polynomial image meets the complete gap21 first image

14 November 2026. PF2OMIT17 closes the disjoint-image gate left by
PF2OMIT16. All distinct additional omitted image supports meet C,
without a degree or high-contact hypothesis on those images. This
gives connected nonproperness support and forces singleton values in
the shared-image and index-three rows. It does not prove properness
or identify a global Keller normalization.

## 1. The remaining two-fibre obstruction

Retain exactly the complete birational polynomial degree(6,10),gap21
first-image family of PF2R7FAM1. Use the polynomial chart L,d in C,
or the polynomially scaled Laurent chart K,z=x^3 nonzero. In the latter
chart p is replaced by x^6*p(t/x), and q by x^10*q(t/x). This is an
invertible diagonal target change for every allowed x. Both new
polynomials have coefficients polynomial in z and are monic of degrees
6,10. Thus every original Laurent point is covered; z=0 is not added.

Put h=H_C=Res_t(p-U,V-q), monic of V-degree6. The written proof of
[PF2OMIT16](DISJOINT_IMAGE_GATE.md) gives the following universal facts:

* h-lambda is reduced irreducible with one infinity branch of delta24;
* its generic affine fibre has genus12 and Euler characteristic-23;
* its total affine critical Milnor number is24, at every family point;
* a distinct disjoint polynomially normalized image would be another
  degree10 fibre h=lambda, lambda!=0;
* both that fibre and C would have twelve nodes, exhausting the Milnor
  budget, with precisely the two critical values0 and lambda.

Let T be multiplication by h on its24-dimensional Jacobi algebra

    J=C[U,V]/(h_U,h_V).

At a critical point its local length is the Milnor number and h acts
by its critical value plus a nilpotent. A disjoint polynomial image
therefore requires

    tr(T^j)=12*lambda^j, j=1,2,3,
    M2=12tr(T^2)-tr(T)^2=0,
    M3=144tr(T^3)-tr(T)^3=0.                       (T1)

Vanishing of either moment alone is not sufficient over C. We exclude
their simultaneous vanishing, rather than making such an inference.

## 2. Traces extend across every basis denominator

The computations below initially use a basis over C(d) or C(z), and
its coefficients have parameter denominators. Dropping their zeros
would not close the whole family. Here is the required extension
argument, using the independently established constant Milnor budget.

Fix any allowed parameter b0. Its gradient has finitely many isolated
zeros, with total length24. Choose small disjoint Milnor neighbourhoods
around these points whose boundaries contain no gradient zero. Local
conservation of isolated complex intersection multiplicity says that
the nearby gradient zeros inside these neighbourhoods still have total
length24. The universal budget24 accounts for every zero in every
nearby parameter fibre. Consequently no additional critical points
can occur outside these neighbourhoods or escape to infinity.

The local critical scheme is therefore finite over a neighbourhood
of b0, and its fibre length is constant24. Over the smooth
one-dimensional parameter base it is flat: for a finite module over
the local discrete valuation ring, any torsion would increase the
special fibre dimension above the generic rank. Equivalently one may
use the analytic local complete intersection and conservation theorem.
Thus its multiplication traces are analytic at b0 and coincide with
the rational traces on the dense open where the chosen basis exists.

This holds at every d in C and every nonzero z. Hence any polynomial
formulas for the generic traces are the actual traces at all these
parameters, including every zero of a certificate denominator. This
argument does not require choosing a fresh basis at each such zero.
The constant Milnor length is essential; an arbitrary rational-basis
matrix would not have this specialization property.

## 3. Two coprime trace moments on each complete chart

The exact trace formulas tr(T),tr(T^2),tr(T^3) have parameter degrees
30,60,90 on each of L and scaled K. M2,M3 are polynomials of degrees
60,90, with nonzero rational scalar denominators. Their primitive
integer numerator polynomials are recorded in
[critical_trace_ledger.json](critical_trace_ledger.json). For **each**
chart these two polynomials are coprime.

The replay verifies a literal Bezout identity modulo101. Both leading
coefficients remain nonzero modulo101, and their degrees are preserved.
Thus a common rational primitive factor could not disappear or lose
degree on reduction. The modular identity proves that the integer
Sylvester resultant is nonzero, hence there is no common complex root.
This is a fixed univariate finite-dimensional inference with checked
degrees, not a generic modular unit ideal in an unbounded affine chart.

By Section2 the computed moments are the actual moments throughout
the allowed parameter domains. They cannot vanish simultaneously.
Condition(T1) is impossible at every L point and every nonzero K point,
so **no polynomially normalized irreducible curve can be disjoint
from any complete gap21 first image C**.

This excludes every possible image degree: different degrees were
already excluded from disjointness by the constant-restriction lemma;
the only remaining degree10 case is now excluded as well.
Equivalently, h=0 is the unique rational fibre of h: every other fibre
still has one infinity place, so genus zero would give a polynomial
A1 normalization and contradict the theorem. Thus all nonzero fibres
have projective normalization genus at least one.

## 4. What the exact certificates prove

Discovery used Singular over Q(d) and Q(z). The default
[replay](verify_critical_trace.py) does not invoke Singular or trust its
standard-basis or trace output. It reads the frozen
[L certificate](critical_trace_L.json.gz) and
[K certificate](critical_trace_K.json.gz), each containing sparse rational
functions with denominators depending only on the family parameter.

For each chart it reconstructs h from the canonical parametrization.
It then checks, by direct denominator-cleared polynomial multiplication:

1. Each of six displayed polynomials G_i is an explicit combination
   of h_U,h_V, using the frozen two-row lift matrix.
2. Their leading U,V monomials leave exactly24 standard monomials.
   Reduction by these leading terms spans a quotient of dimension at
   most24. Since the independent geometric Milnor dimension is24 and
   (G_i) is contained in the gradient ideal, the ideals are equal and
   these24 monomials form a basis. No reverse solver membership or
   unverified Groebner claim is needed.
3. Every column of the displayed24-by24 matrix satisfies the literal
   identity h*b_j-sum_i T_ij*b_i=sum_k Q_kj*G_k.
4. A common scalar-parameter denominator D converts T to a polynomial
   matrix A. Direct products verify tr(A^j)=D^j*tr(T^j), j=1,2,3.
5. The two moment formulas, their nonzero scalar relation to the
   primitive integer polynomials, leading-degree preservation, and
   both modular Bezout identities hold exactly.

The arithmetic uses sparse rational multivariate polynomials and
univariate polynomial matrices. Negative controls mutate one matrix
entry and one Bezout coefficient. Their identities must fail.

The rank24 geometry and flat extension in Sections1--2 remain written
proof arguments. The certificates are not a formalization or external
review. Earlier symbolic Python reductions reached terminal caps;
their partial dimension output is not a proof premise. The native
discovery and frozen cofactor identities replace that incomplete route.

## 5. Complete consequences for omission, images and exceptional lengths

Apply the theorem to each distinct additional image in PF2OMIT13's
four oriented residual rows. Their normalizations are A1, so all meet
C. Therefore the full nonproperness support is connected in **every**
row and grouping. This gives no contradiction from a theorem requiring
a smooth nonproperness set: our support is singular.

Shared-image row A and cyclic-index-three row C now necessarily have
singleton source values, without any degree, contact or coefficient-open
restriction on B. Every point of C intersect B is such a value:

| Grouping | Original omitted local lengths | Retained fibre | Radial omitted indices |
|---|---|---|---|
| A, H=C or H=B | 2+2+1 | one reduced point | 2+2+1 |
| C | 2+3 | one reduced point | 2+3 |

All original omitted residue degrees remain one. In A the length-one
factor is unramified and is a distinct omitted point even over the
shared image. Its deletion is essential to obtaining retained length1.
The inverse coefficient/index bounds of PF2OMIT10 apply at every one
of these necessarily nonempty exceptional loci.

If deg B=5k, the intersection polynomial H_B(p,q) is now always
nonconstant, of degree in <6,9,10> intersect[6,30k]. Its degree is the
total affine intersection multiplicity, not the number of distinct
target values. The stronger9k odd bound and the k1 lower bound10 remain
valid. Actual tangencies and values still require actual B equations.

In rows B,N and in A with three distinct images, a C/B intersection
does not by itself force a singleton. The cusp incidence and the
unramified image branch count remain necessary. In particular this
does not overturn the high-contact degree(3,5) no-singleton conclusions
of PF2R7B5N1 and PF2R7B5B1.

Actual additional equations, all exceptional values and their normal
fibre factors, infinity attachments, the S6 representation and the
identification of the original source A2 remain missing. Even in the
two groupings with mandatory singletons, existence of a singleton is
not a properness proof. The full research goal stays open.
