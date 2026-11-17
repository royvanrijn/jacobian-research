# Generic-field certificates and numerical Gaussian family bounds

16 November 2026. Internal deduction and executable general producer.
Independent, formal and external assurance remain false. This implements
the numerical generic-component gate identified by the
[family theorem](GENERAL_GAUSSIAN_FAMILY_LOCUS_AND_NORMALIZATION.md).
It does not classify all safe phases or compute least cutoffs.

The [producer](gaussian_generic_cutoff_producer.py),
[relative Weyl backend](gaussian_generic_field_core.m2),
[focused checker](verify_gaussian_generic_cutoff_producer.py),
[receipt](gaussian_generic_cutoff_producer_verification.json) and
[source audit](gaussian_generic_cutoff_producer_source_audit.json) retain
the certificates and implementation boundary. Earlier sealed artifacts
are unchanged.

## 1. Exact generic fields and original operators

Let B=QQ[c_1,...,c_r], let p be a proper prime, and set K=Frac(B/p).
The phase P(c,X) is a polynomial over B, and the symmetric nonsingular
Gaussian precision H has rational entries. Its covariance is H^-1.
The relative exponential connection is

    I=(D_i+(HX)_i-t P_i, D_t-P) in K<X,t,D_X,D_t>.

The parameters are constants for this relative Weyl algebra. This is
distinct from the full JOINT source, including parameter derivatives,
used by the family solver to find safe components. We do not replace
joint constructibility by a claim about arbitrary central parameters.

The direct backend follows the full integration presentation in
[Oaku's integration paper](https://www.lab.twcu.ac.jp/oaku/Oaku2013.pdf),
physical pp14-15, before its cyclic extraction. It uses the public
Macaulay2 `gbw` and `bFunction` operations, Fourier transforms, and the
maximal integer root B0 of the integration polynomial b. It retains all
positive original amplitudes X^alpha with |alpha|<=B0 and the full
projected relation matrix before pruning. No time localization is used.

The native `Dintegration` route over K reached a QQ-only b-function
factor formatter and failed. A Singular parameter/minimal-polynomial
attempt rejected a multivariable minimal polynomial and continued with
an unintended transcendental parameter; it is retained as a FAILURE.
The implemented route constructs the finite presentation directly and
does not modify installed libraries or interpret either failure as a
successful generic computation.

Macaulay2 matrices represent left-module maps: entrywise ring products
cannot be used to guess coefficient orientation. The retained controls
check that exported lift entries are LEFT coefficients. Python replays
their Weyl products in that order.

## 2. Integer roots over a quotient function field

For f(z)=sum a_i z^i in K[z], write a_i=n_i/d_i with nonzero d_i in B/p.
Multiply by the common nonzero denominator and reduce every coefficient
to its normal form modulo p. Expand the result in parameter standard
monomials:

    (product d_i) f(z)=sum_e c^e q_e(z), q_e in QQ[z].

For an integer m, f(m)=0 in K iff every q_e(m)=0, since standard monomials
are linearly independent over QQ. The integer roots are therefore exactly
the integer roots of gcd_e q_e. Rational factorization finds them; only
linear factors with integer root contribute. Nonlinear factors are not
assumed to split over QQ. The same operation computes integration roots
and scalar indicial roots. Denominator products, coefficient normal forms,
rational content polynomial and integer roots are retained.

These are GENERIC integer roots. Root sets may change at special points.
No assertion of constancy of integer roots under specialization is used.
Primality is checked by Macaulay2; the finite Python replayer assumes the
supplied quotient is a domain and does not independently prove primality.

## 3. Finite certificates select the actual formal periods

The proof in the sealed
[pointwise note](WHOLE_GAUSSIAN_SOURCE_HOM_AND_EFFECTIVE_CUTOFFS.md)
works over K with precision H. This implementation retains and replays:

1. An original LEFT Weyl identity L=sum W_i I_i=b(Theta)+R, where
   Theta=-sum D_i X_i and every remainder monomial X^u t^k D_X^v D_t^j
   has degree drop |v|-|u|>=1.
2. Polynomial IBP primitives for every relation column, expressed using
   the original connection partial_i-(HX)_i+t P_i.
3. A nonzero scalar annihilator A_alpha for EVERY base amplitude, its
   LEFT relation witnesses, and its directly checked original IBP identity.
4. Both scalar indicial polynomials, using falling factorials (z)_b:
   the maximum a-b bounds polynomial degree by D_alpha; the minimum a-b
   bounds the initial order of any nonzero formal solution by R_alpha.
5. Exact original Wick coefficients E[X^alpha P^k]/k! through
   max(D_alpha,R_alpha), using rational covariance H^-1.

The implementation can avoid a polynomial-Hom linear solve. Any polynomial
actual period has degree <=D_alpha, so it MUST equal the actual truncation
through D_alpha. Check that this candidate solves A_alpha and that its
jet agrees with the actual series through R_alpha. If both checks pass,
origin-jet injectivity proves equality with the actual formal series.
If either fails, no polynomial actual period exists. Empty nonnegative
root sets mean bound -1 and the zero candidate, with the same argument.

The recurrence gives, for every original monomial marker,

    b(|alpha|) F_alpha
      =-sum_R c (-1)^|v| (alpha+u)_v t^k D_t^j F_(alpha+u-v).

It covers all larger amplitudes by induction. Thus all base periods are
polynomial iff every original marker period is polynomial. This finite
decision needs the retained original identities, roots and actual jets;
it does NOT need an independently proved completeness of the CAS relation
matrix or a computation of the actual-period image kernel.

For a safe generic phase set

    C=max(0,max_R (k-j)/(|v|-|u|)),
    E=max(0,max_base (deg F_alpha-C|alpha|)).

The recurrence proves E[Q P^m]=0 whenever m>C deg(Q)+E. The actual unit
period is tested for equality to1 before declaring generic purity. When
pure, the finite positive-degree base ratios give a sufficient pure
linear bound. Independently C+E is a valid conditional pure linear slope:
positive marker degrees absorb E, and pure nullity handles degree zero.

All numerical bounds here are sufficient; no minimality is claimed.

## 4. Bounds extend through every degeneration of a safe component

Suppose p is a safe component returned by the exact joint family solver.
Every forbidden original moment E[X^alpha P^m] is an element of B/p.
If its image in K vanishes, it is already zero in B/p. Consequently the
generic C,E hold at EVERY field point of the component, including singular
points and points where coefficient denominators of the generic operators
vanish. The operators themselves specialize only where their retained
denominators are nonzero. These are separate assertions.

Take maxima of C and E over the finitely many safe prime closures. This
computes a common affine bound on all safe field fibers of ANY polynomial
parameter family. There is no parameter dimension, spatial dimension,
phase degree or parameter degree cutoff in the interface. Native CAS
computation may exhaust the default30second cap and return UNKNOWN; zero
disables the cap. Errors fail closed. No efficiency bound is asserted.

In particular the complete bounded coefficient space is constructed by
`gaussian_family_producer.universal_family(n,d)`. Passing it to the new
family-bound producer is a general algorithm for numerical constants for
each n,d. This does not assert that large universal spaces were executed.

## 5. Quantitative ordinary torsion over presented coefficient rings

Let A=B/I0, with I0 kept ORDINARY. Let J be the computed radical ideal of
safe field parameters in B, and let N=Gamma_(JA)(A) be the exact polynomial
Hom normalization ideal from the sealed family theorem. The family output
retains a positive L with N*(JA)^L=0. The new Python layer separately
checks this finite ordinary annihilation without radicalizing I0. Exact
maximality of N and complete safe-component coverage remain contracts of
the trusted joint family solver, rather than consequences of that finite
annihilation check.

Let D be the spatial degree of the ORIGINAL P in B[X], before any generic
prime reduction, and r the number of ambient parameters. Write
mu_(alpha,m)=E[X^alpha P^m] in B. A parameter derivative of order k is a
sum of terms

    (m)_j E[X^alpha (product_(l=1..j) partial^(nu_l)P) P^(m-j)],
    0<=j<=k, sum_l |nu_l|=k, |nu_l|>=1.

The formula follows by induction, differentiating the power or one of
the existing derivative factors. Each marker in the product has spatial
degree <=|alpha|+Dj. If all safe field fibers satisfy the bound C,E,
all parameter derivatives of order <k vanish on V(J) whenever

    m>C|alpha|+E+(CD+1)(k-1).

The previously imported Zariski--Nagata theorem in the REGULAR ambient
ring B identifies this differential vanishing with J^(k), the symbolic
power. The retained Ein--Lazarsfeld--Smith containment, using height <=r,
gives J^(rL) subset J^L. Therefore, for L>=2,

    mu_(alpha,m) in J^L
      if m>C|alpha|+E+(CD+1)(rL-1).                      (1)

For L=1, ordinary vanishing itself gives mu in J under m>C|alpha|+E;
no symbolic-power bound is needed. These imports are applied in B,
NEVER in the possibly singular/nonreduced quotient A.

Multiplying (1) by ANY normalization scalar a in N proves

    a E[Q P^m]=0 for m>C deg(Q)+E_ring,
    E_ring=E+(CD+1)(rL-1) for L>=2; E_ring=E for L=1.

This is a numerical normalization bound over the entire presented ring.
It applies to a=1 precisely when the family solver certifies normalization1.
Residue-field purity alone does not imply ring purity: for P=eX over
QQ[e]/(e^3), the actual unit period is1+e^2 t^2/2. Its retained output
has C=0,E_ring=2, L=3, proving safety without asserting purity.

For the full affine cutoff ideal
K_N=(mu_(alpha,m):m>N(|alpha|+1)), the sufficient forward cofinal index is
ceil(max(C,E+(CD+1)(rL-1))) for L>=2 and ceil(max(C,E)) for L=1.
This is quantitative ordinary cofinality BEFORE a pure quotient.
The algorithm still does not compute every ideal K_N or its least index.

## 6. Executed scope, credit and remaining work

Seven generic controls include a transcendental unsafe line, a safe
isotropic cone with denominators vanishing at its vertex, a split cubic,
an elliptic function field with both coefficients used, a quartic with
two amplitudes, a zero-dimensional prime and a number-field prime.
Four presented-family controls retain ordinary thick torsion, a nonpure
nilpotent safe phase, the split cubic and the full universal affine family.
All original identities and jets replay. Root-content controls and
deliberately damaged certificates check the new boundaries.

Oaku integration, scalar indicial bounds, Gaussian IBP, classical
constructibility/local-cohomology methods and differential/symbolic powers
are credited to their sources. This is a general implementation and
combination of retained results, with no priority claim. The
[Oaku1609.04676 packet](OAKU_GRAPH_CERTIFICATES_AND_ACTUAL_POLYNOMIAL_DECISION.md)
keeps its original operator identities, torsion witnesses and specialization
thresholds. The
[classical-credit boundary](BERNSTEIN_CLASSICAL_CREDIT_AND_FORMAL_EXTENSIONS.md)
credits Sabbah/Maisonobe for the convergent rational-slope mechanism;
arbitrary nonconvergent formal extension and all-ring descent remain
separate candidate extension scopes with priority unassessed.

Full ordinary cutoff-ideal producers, least-bound extraction and structural
safe normal forms remain general obligations. Earlier constructibility,
field uniformity and formal/all-ring theorems stay closed at their recorded
scopes. The broader research goal remains active; the historical capped
cusp family remains UNKNOWN and was not restarted.
