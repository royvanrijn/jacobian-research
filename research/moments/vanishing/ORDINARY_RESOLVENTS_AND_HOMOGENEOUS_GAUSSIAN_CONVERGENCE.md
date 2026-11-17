# Ordinary resolvents, exponential periods and homogeneous Gaussian convergence

16 November 2026. Sourced clarification and consequences of RETAINED
registered theorems; no new full-closure claim or duplicate theorem entry.
Internal written arguments, without independent replay, formal verification
or external review. Both full exceptional Mathieu--Zhao cases remain OPEN.
External review remains deferred and every earlier goal/result is preserved.

The user's generating-function distinction is essential. The new weighted
Gaussian theorem gives algebraic EXPONENTIAL periods in its architecture.
The full Mathieu algebraic criterion concerns ORDINARY ACTUAL marked
resolvents. Algebraicity is not preserved by a general factorial/Borel
conversion. The previous weighted-symmetry proof remains valid because it
obtains POLYNOMIALITY of the actual exponential periods by the fixed-linear-
marker criterion; it does not substitute algebraicity for polynomiality.

There is an additional valid CONDITIONAL conversion for the ACTUAL canonical
homogeneous four-Gaussian lifts of degree D>2. It is already a consequence
of GMCHOMBOREL20261116 and MVCOMPACTREGULAR20261116: convergence of a marked
Gaussian exponential series forces that marked tail to terminate. It uses
homogeneity, the exact radial factors and compact regular-singular periods.
None of these facts proves convergence or algebraicity for every marked
period of the surviving exceptional phases. The new weighted theorem does
not supply that missing coverage.

## 1. Three series and their exact coefficients

Fix an original phase f on SL2/Gamma and ONE fixed original marker q.
Write a_m=h(q f^m), and let P,Q be their canonical determinant-homogeneous
lifts of prescribed even degrees D=2d,E=2e in FOUR Gaussian variables.
The exact radial bridge gives

    M_m=E_G[Q P^m]=2^(dm+e)(dm+e+1)! a_m.             (1)

Keep the following objects separate:

| Object | Series | Coefficient normalization |
| --- | --- | --- |
| Ordinary Haar series | R_q(z)=sum a_m z^m | no factorial |
| Haar exponential period | E_q(t)=sum a_m t^m/m! | divide by m! |
| Gaussian exponential period | G_Q(t)=sum M_m t^m/m! | multiply by (1), then divide by m! |

The actual ordinary Cauchy resolvent is

    C_q(s)=sum a_m s^(-m-1)=R_q(1/s)/s.               (2)

Haar exponential periods are entire by compactness. Gaussian exponential
periods in the last row need not converge at all. Being homogeneous of
ordinary degree D is a DIFFERENT hypothesis from having Gaussian
covariant weight2 in a supplied coordinate weighting.

## 2. The actual full Mathieu criterion uses the ordinary object

The retained
[ordinary algebraic-splitting theorem](ALGEBRAIC_RESOLVENT_SPLITTING_AND_ZERO_VALUE_LOGARITHMS.md)
(MVALGORDSAFE20261116) says, under purity and its established exact
actual-source/time-zero hypotheses,

    full safety iff EVERY original C_q(s) is algebraic. (3)

Its trace/principal-part retraction acts on algebraic ORDINARY Laurent
functions in s. It is not a retraction supplied by algebraicity of G_Q(t),
or by a generic Gaussian time connection with finite monodromy.

The retained sufficient monodromy condition is on the COMPLEMENT/CAUCHY
Gauss--Manin connection for

    U_s=X minus {f=s}, forms q Omega/(s-f),

at PHASE VALUE s=0: semisimplicity on generalized eigenvalue1 suffices.
The exact local criterion is absence of integer-exponent logarithms in
ALL actual original marked C_q there. Generic-fibre monodromy or a time-
localized exponential connection cannot replace this connection.

Thus the genuine whole-quotient obligation remains to establish (3)'s
premise, this sufficient complement condition, the exact intrinsic
splitting criterion, or an actual pure failure. Purity makes the UNIT
exponential period polynomial automatically; algebraicity of that ONE
period provides no all-marker conclusion.

## 3. An actual weighted Gaussian countercontrol

The completed [algebraic Gaussian construction](ALGEBRAIC_FIXED_POINT_LIFTS_AND_UNIVERSAL_GAUSSIAN_IMAGES.md)
already supplies an exact countercontrol. Take its FIVE-variable version
with H(Z)=1+Z. In split centered Gaussian coordinates Z,W,Y,V,T with
E[ZW]=E[YV]=E[T^2]=1, its phase is

    P=W(1+Z)+V(1+Y)
        -(1/2)(1+Y)(4+6Y+4Y^2+Y^3)T^2.             (4)

It is PURE at every positive order. Its fixed ORIGINAL linear marker Z
has the actual exponential period

    G_Z(t)=g(t)=t/(1-t),
    E_G[Z P^m]=m!, m>=1.                            (5)

ALL marked exponential periods in this construction are algebraic; for
this member they are rational functions of t. Yet the ordinary marked
series sum_(m>=1) m! z^m has radius ZERO and is transcendental as a formal
series. An algebraic formal power series over C(z) has positive radius,
so it cannot be this factorial series.

The phase (4) has the supplied weights0 on Z,Y,2 on W,V,1 on T. It has
linear terms and is NOT homogeneous of one ordinary degree. It therefore
neither violates the no-linear-term corollary nor the homogeneous
convergence theorem below. It is not an exceptional invariant phase.
This is an actual all-order Gaussian control imported from the completed
construction, not an abstract module or a newly inferred finite example.

## 4. A conditional conversion that is valid for the actual homogeneous lift

Let D>2, and retain ONE fixed original q and its homogeneous lift Q.
For these actual lifts the following are equivalent:

    G_Q has positive convergence radius,
    G_Q is algebraic,
    G_Q is polynomial,
    R_q is polynomial,
    a_m=0 for every sufficiently large m.             (6)

The nontrivial implication is convergence=>termination. If G_Q has positive
radius, Cauchy's coefficient estimate and (1) give constants C,A>0 with

    |a_m| <= C A^m m! / [2^(dm+e)(dm+e+1)!].          (7)

Here d=D/2>1. Stirling's formula makes the mth root on the right tend to0.
Hence R_q is ENTIRE. But it is an ACTUAL compact algebraic resolvent.
MVCOMPACTREGULAR20261116 puts it in a regular-singular system at every
point of P1 including infinity. An entire branch has only power-logarithmic
growth at infinity and is a polynomial by Cauchy estimates. This proves
termination. The same conclusion is the D>2,s=1 case of the retained
[homogeneous Borel trichotomy](COMPLEX_HOMOGENEOUS_GROWTH_AND_BOREL_TRICHOTOMY.md).

Algebraic formal series over C(t) converge, so algebraicity implies
convergence. Termination makes both series polynomial; the nonzero radial
factors make their finite tails equivalent. No purity premise is needed
for this FIXED-marker statement. Its proof does not claim that Borel
conversion preserves algebraic functions.

Quantifying (6) over EVERY original q shows: if ALL actual marked Gaussian
exponential series of the canonical homogeneous lift converge (in
particular, if they are all algebraic), then ALL original ordinary
resolvents are finite Laurent polynomials and the original phase is safe.
One can alternatively use every Gaussian marker and GMCHOMBOREL directly;
Reynolds and the radial theorem retain the entire original-marker packet.
The marker is fixed before taking m to infinity; it does not change with m.

This is a conditional bridge, already supported by retained theorems.
The premise for the surviving exceptional phases is UNPROVED. It is not
supplied by the unit series, by symmetry of linear markers alone outside
the completed weighted architecture, or by an algebraic solution of a
larger presentation that has not been identified with the actual periods.

The threshold D>2 matters. At D=2 and e=0, a_m=1 gives the actual radial
Gaussian series

    sum 2^m(m+1)! t^m/m! = (1-2t)^-2,

algebraic and nonpolynomial. The corresponding original phase is the
constant1, not pure. Higher-degree homogeneity cannot be dropped either,
as the actual inhomogeneous countercontrol (4)-(5) demonstrates.

## 5. Why this does not supply the missing exceptional coverage

The [weighted symmetry/coverage proof](WEIGHTED_GAUSSIAN_SYMMETRY_AND_EXCEPTIONAL_COVERAGE.md)
was checked at its actual logical steps. It uses the fixed-linear-
marker POLYNOMIALITY criterion, then exact safety transfer. Its independent
coverage theorem says that a nonzero canonical homogeneous lift of even
D>2 admits the supplied LINEAR weighting iff some left SL2 basis gives a
single weight |k|=D orD-2. These cases are already torus safe. Surviving
mixed phases are outside that architecture. Both the symmetry argument
and (6) therefore give no new surviving case from the weighted Gaussian
theorem as presently stated.

For the full goal, the next meaningful advance must control ALL actual
ordinary marked C_q, the correct complement monodromy, or the intrinsic
boundary extension on the entire pure locus. A proof of actual marked
Gaussian convergence for the entire surviving homogeneous lift family
would also suffice BY (6), but none is established. This is an alternative
expression of the same all-marker obligation, not a new decision procedure
or a replacement success criterion.

Both full exceptional cases remain OPEN. Retain the capacity/support
bounds, all automorphism/primitive constraints, completed family transfer,
fixed original bad-marker requirement and every earlier result/goal.
External review remains deferred. No new degree/prefix campaign, Gaussian
backend replay or external theorem import is introduced.
