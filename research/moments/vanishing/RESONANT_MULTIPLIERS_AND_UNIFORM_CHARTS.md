# Resonant multiplier degrees and uniform certificate charts

15 November 2026. Internal structural refinement of the
[b-function cutoff recurrence](LINEAR_MIXED_CUTOFF_FROM_B_FUNCTION.md).
No novelty, independent review or formal verification is claimed.
The [checker](verify_resonant_multiplier_charts.py) replays exact recurrence,
residue and parameter-family controls. It is not a general b-function or
coefficient-stratification solver. Uniform bounds on the whole universal
degree-D coefficient space remain unproved.

## 1. Only integer-root degrees can introduce a new bad multiplier

Retain a valid polynomial Weyl identity

    (b(theta)+R) exp(Psi)=0,
    theta=-sum_i partial_i X_i,
    Psi=tP-sum_i X_i^2/2,

with every normal-order term c*t^k*Dt^j*X^u*partial_X^v of R satisfying
r=|v|-|u|>=1. The preceding proof gives the module identity and its Gaussian
image recurrence

    b(d) F_alpha
      =-sum c*(-1)^|v|*(alpha+u)_v*t^k*Dt^j F_(alpha+u-v),
    d=|alpha|.                                                (1)

Every child has smaller amplitude degree. Let

    S={d>=0 integer:b(d)=0}.

The earlier proof shows 0 is in S. All the following statements apply to
any valid certificate, including a nonminimal one.

**Resonance theorem.** The images of monomial amplitudes whose degrees
belong to S generate the entire polynomial integration module over the
polynomial time Weyl algebra. Hence the eventual mixed property is
equivalent to polynomiality of F_alpha only for |alpha| in S.

For the proof, use induction on degree. At a nonroot degree, (1) expresses
every new class in smaller-degree classes. At a root degree retain the
class as a generator. This also works in the full module before evaluating
the Gaussian functional. No rational division in t occurs.

Suppose P is pure-null and has a nonpolynomial mixed series. Let e be the
least degree of a multiplier giving one. A monomial of degree at most e
must already give a nonpolynomial series, since finite linear combinations
of polynomials are polynomial. All smaller-degree monomial images are
polynomial. Equation (1) then forces

    e>=1 and b(e)=0.                                          (2)

In particular the least failure-multiplier degree must be a positive
integer root of the minimal integration b-function, not merely less than
its largest root.

To make the last statement precise, the polynomials b admitting a lift of
the displayed form are an ideal in C[s]: multiply a certificate by any
polynomial in theta, whose weight is zero. Their monic generator is the
minimal integration b-function for the fixed cyclic generator and weight.
It divides every certificate polynomial, and the argument applies to its
own lift. No equality between arbitrary nonminimal root lists is asserted.

A useful sufficient criterion follows immediately:

    If P is pure-null and a valid b has no positive integer root,
    then P is GMC-safe.                                      (3)

This is an all-degree, all-dimension class criterion. It is not a necessary
condition, even if the b-function is minimal; section 3 gives a safe
example whose full module requires positive resonance.

The cutoff extraction can also use only resonances. Put
C_rec=max(0,max_terms (k-j)/r). If all resonant images are polynomial,
their degrees set the affine offset

    A=max(0,max_(|alpha| in S,F_alpha!=0)
                         (deg F_alpha-C_rec |alpha|)).

For a pure-null polynomial take instead

    C=max(C_rec,max_(|alpha| in S,|alpha|>0,F_alpha!=0)
                                    deg F_alpha/|alpha|).

Induction at nonroot degrees and the given bounds at root degrees prove
the same affine and homogeneous cutoffs as before. Filling every degree
below max S is unnecessary.

## 2. Exact controls include a deliberately gapped root list

For P=W in a circular Gaussian pair, the preceding note retains

    b_0(s)=s,   R_0=-2 partial_Z partial_W+t partial_Z.

Multiplying this annihilator on the left by theta-3 gives

    b_3(s)=s(s-3),
    R_3=2Z partial_Z^2 partial_W+2W partial_Z partial_W^2
        +10 partial_Z partial_W-tZ partial_Z^2
        -tW partial_Z partial_W-5t partial_Z.

This certificate has root degrees {0,3}, with gaps at one and two.
The checker starts only from monomials of degrees zero and three and
reconstructs every other tested amplitude by (1), then compares the result
to the exact all-order formula for P=W. The general induction above, not
the finite replay, proves the infinite statement. The extracted slope is
one; no minimality of b_3 is claimed.

## 3. Positive resonance can belong to a safe kernel direction

Take the pure-null, GMC-safe polynomial

    P=Z^2 W,       E(ZW)=1.

Every monomial has positive angular weight. In fact

    F_(Z^a W^b)= (2b-a)!/(b-a)! * t^(b-a) if b>=a,
                 0 otherwise.                               (4)

This follows by matching the circular exponents in
Z^a W^b (Z^2 W)^m. It proves safety with sharp coefficient one.
Thus the actual Gaussian image is the polynomial module.

Nevertheless the full polynomial integration module H is not generated
by [1]. Here is an algebraic certificate, without an analytic cycle.
Let g(Z,t)=Z-tZ^2. For A=sum_b A_b(Z,t)W^b define a rational time function

    L_infinity([A])
      =Res_(Z=1/t) sum_b b!*A_b(Z,t)/g(Z,t)^(b+1) dZ.         (5)

For the phase Psi=-Wg, its total derivative in W is partial_W-g.
In (5), the partial_W contribution cancels the -g contribution by shifting
b; the leftover polynomial term has residue zero. The Z total derivative
is partial_Z-Wg_Z; its image is the residue of an exact Z derivative.
Finally nabla_t=partial_t+Z^2 W agrees with differentiating (5), since
g_t=-Z^2. Differentiation of a moving-pole residue commutes with the
parameter derivative: changing the center contributes an exact Z
derivative with zero residue. Consequently (5) is a well-defined
D_t-linear map from H to C(t).

It has

    L_infinity([1])=-1,      L_infinity([Z])=-1/t.              (6)

Every polynomial Weyl operator applied to -1 gives a polynomial in t.
Therefore [Z] cannot lie in D_t[1]. If a valid b had only zero as a
nonnegative integer root, the module induction in section 1 would make
D_t[1]=H. This contradiction proves that every valid integration b for
this safe P, including the minimal b, has a positive integer root.

This is the required boundary to (3). Resonance locates degrees at which
a new class can appear; it does not say the actual Gaussian functional
is nonzero on that class or that its image is nonpolynomial. The residue
map (5) is not the Gaussian map: it has a pole at t=0 and selects the
other algebraic branch. No numerical minimal b-function is computed.

## 4. What data suffice to prove uniformity on a coefficient chart

Consider a polynomial coefficient family over a locally closed rational
chart U. Its coordinate ring can be represented by a polynomial ring
modulo equations, with one inverse variable for an excluded denominator.
Suppose the following finite, specialization-valid data are supplied:

1. A single nonzero scalar b in Q[s] and a lift b(theta)+R valid over the
   chart coordinate ring, hence at every point, with regular parameter
   coefficients and a fixed finite
   negative-weight normal-order support for R.
2. For each positive root degree d in S and each monomial X^alpha of that
   degree, a uniform integer H_alpha>=0 such that any polynomial series
   F_alpha at any point of U has degree at most H_alpha.

For example condition 2 follows from scalar annihilators whose
specializations remain nonzero and whose infinity indicial polynomials
have no nonnegative integer roots above H_alpha at any point. A generic
annihilator or generic root bound alone is insufficient.

**Uniform chart theorem.** Under these supplied hypotheses the pure-null
GMC-safe locus in U is an effectively computable relatively closed
algebraic set. Every member of it has the common sufficient coefficient

    C_U=max(C_rec,max_positive_resonance H_alpha/|alpha|).     (7)

If there are no positive resonance degrees, the second maximum is absent.
Vanishing coefficients of R after specialization can only improve the
support-based upper bound C_rec.

Indeed impose all positive pure moments zero and, for every listed alpha,

    E(X^alpha P^m)=0 for all m>H_alpha.                       (8)

If P is safe, hypothesis 2 forces (8). Conversely (8) makes the resonant
series polynomial, while pure nullity handles the zero-degree root.
Section 1 then proves safety and (7). Conditions (8) and pure nullity
are polynomial coefficient equations, so they define a closed set in U.

Their full ordinary ideal is computable by the preceding
[primary-thickening oracle](PRIMARY_THICKENING_MOMENT_IDEAL_ALGORITHM.md).
For clarity, the needed marked-tail extension is as follows. For a finite
list of seeds Q_j and starting indices N_j, the series

    E(Q_j P^(N_j) exp(tP))

contains precisely the moments E(Q_j P^m), m>=N_j. In each finite algebra
of a full primary decomposition, its amplitude is a seed in the same
finite-rank exponential connection. Scalar holonomic annihilators and
origin indicial tests decide whether the entire algebra-valued series
is zero. This decides membership of every tail coefficient in a supplied
ordinary ideal. Adjoining a returned coefficient outside the ideal
strictly increases it, so ACC terminates and computes the ideal of all
the finitely many tails. Include Q=1,N=1 for the pure equations and
Q=X^alpha,N=H_alpha+1 for (8). Denominator inverses are represented by
their polynomial relations. This retains primary thickness without
claiming an intrinsic scheme structure for the abstract safety property.

If a finite such chart cover of the universal fixed-n,D coefficient
space is constructed, the union of the chart loci is constructible and
max_U C_U is a uniform fixed-n,D safe cutoff. This is a sufficient finite
certificate route. Such a cover has not been proved or constructed for
the full universal Gaussian family. The theorem does not assert its
existence, implement a general stratification algorithm, or close global
Gaussian constructibility or uniform bounds.

## 5. A genuine uniform chart in every degree

Let P=p(Z)=sum_(k=0)^D c_k Z^k in a circular pair. This whole family has
the fixed certificate b(s)=s with

    R=-2 partial_Z partial_W
      -t sum_(k=1)^D k*c_k*(-1)^k partial_W^k.                 (9)

For Psi=t p(Z)-ZW, the amplitude of theta-2 partial_Z partial_W
is t Zp'(Z); the remaining terms cancel it. Thus (9) is an exact identity
over the universal coefficient ring, including all zero coefficients.

There are no positive resonances. The pure series is exp(c_0 t), and the
pure-null/safe locus is precisely c_0=0. Its ordinary positive-pure ideal
is (c_0), since every moment is c_0^m. The whole pure-null chart has the
uniform cutoff coefficient one.

For a nonzero polynomial with c_0=0 let k_0 be its least positive nonzero
exponent. Specializing the support in (9) gives C=1/k_0. This is sharp:
F_(W^k_0)=k_0!*c_(k_0)*t, so any smaller coefficient fails.
For P=0 use C=0. These elementary examples validate the all-degree chart
mechanism; they are not presented as a new classification of one-sided
supports. The checker verifies a symbolic degree-four certificate as a
transcription control, while (9) proves every degree.

## 6. Why a known Bernstein-Sato stratification theorem does not close the gap

[Leykin, arXiv:math/0003155v1, 24 March 2000](https://arxiv.org/pdf/math/0003155v1)
was read in the introduction and sections 4--6, including lemma 4.2,
algorithms 6.2/6.5 and theorem 6.4 with its proof. It treats the classical
Bernstein-Sato polynomial of f^s, gives specialization control for Weyl
Groebner computations and obtains constructible constant-b strata.
The later journal title is *Constructibility of the Set of Polynomials
with a Fixed Bernstein-Sato Polynomial: an Algorithmic Approach* (2001).
The retained body is the 2000 arXiv version.

This is not a theorem supplying the two chart hypotheses above for the
Gaussian exponential integration module. The modules and b-functions
differ, and a bound on polynomial degrees of resonant amplitude images
is an additional requirement. No theorem transferring those hypotheses
is imported; no producer implementation was run.

The scalar holonomic family

    (1+t) f'-a f=0,      f(0)=1

shows why generic specialization and pointwise decidability alone cannot
supply the missing uniform degree bound. Its series (1+t)^a is polynomial
exactly for a in Z>=0, with degree a. The degree identity follows from the
infinity indicial polynomial z-a. Thus its polynomiality locus is an
infinite, proper countable subset of A^1_C, which is not constructible:
a constructible subset of A^1_C is finite or cofinite.
There is no uniform polynomial degree bound on that locus.

This is an obstruction in general holonomic parameter families. It is not
a Gaussian-moment family and does not disprove uniformity or constructibility
for Gaussian polynomials. It prevents a false inference from the available
generic algorithms.

The next general task is to establish the chart hypotheses for Gaussian
exponential modules, or identify a genuine Gaussian obstruction. Root-degree
resonance and the safe kernel example must be retained in that investigation.
The broad research goal remains active; fixed-cubic optimization stays
secondary. [Source audit](resonant_multiplier_source_audit.json).

