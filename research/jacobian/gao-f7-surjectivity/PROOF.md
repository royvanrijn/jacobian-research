# Gao's F7 is surjective over C

14 November 2026. Registered claim **GAOF7SURJ1**.

For the specific F7 of [Gao, arXiv:2608.00222v2](https://arxiv.org/html/2608.00222v2),
every complex target has a preimage. This settles the empty-fibre question
left open in our F7 follow-up. It does not settle the three-point question,
all multiplicity strata, F5, or the planar Jacobian conjecture.

The conclusion concerns this explicit five-dimensional map. Surjectivity
in other Keller families has already been reported, for example in the
[threshold notes and their linked Corollary S](https://github.com/dasjoms/jacobian-conjecture-counterexample-exploration/blob/main/jacobian_threshold_deck.md).
We claim neither the first surjective Keller counterexample nor established
novelty priority. No theorem from that other family is used in this proof.

## Prior inputs and the exact target

Gao supplies F7, its triangular fibre reconstruction, its curve-type
simple-root criterion and generic degree twelve. The earlier
[GAOF7BND1 proof](../gao-f7-fibres/PROOF.md) supplies the complete boundary
count: a target (0,A,B,D,Z) has one point if A=0 and two if A!=0.
Thus only the chart c!=0 remains for surjectivity.

On that chart set Y=(cA,c^2 B,c^3 D,c^4 Z), and define

    b(t)=Y4-3Y2*t^2+2Y1*t^3-t^4,
    a(t)=Y3-2Y2*t+Y1*t^2-t^3-t^2*b(t),
    R(t)=a(t)^2+a(t)+t*b(t)+t^3+Y1*t-Y2.

R is monic of degree twelve. Precisely its simple roots lift to actual
source points, by the prior criterion gamma=R' modulo R and the unique
stage inverse. We prove that R has at least one simple root for every
Y in C4. All variables and factors below are over C.

## Every possible empty fibre reduces to two factorizations

Complete the square, writing F=a+1/2:

    R=F^2+L,
    L=t*b(t)+t^3+Y1*t-Y2-1/4.

F is monic of degree six, whereas L has degree five with leading
coefficient -1. Therefore R is never a square: if R=H^2 with H monic
of degree six, then (H-F)(H+F)=L. The second factor has degree six;
the first is either zero or has nonnegative degree. Neither possibility
gives the nonzero degree-five polynomial L.

Suppose now that R has no simple root. Every root multiplicity is at
least two, and at least one is odd. The number of odd multiplicities is
even, at most four, and hence either two or four.

If there are two, form Q2 as the monic product of these two distinct
roots. Each odd multiplicity is at least three. Subtracting three from
both makes every remaining multiplicity nonnegative and even, of total
degree six. Thus

    R=Q2^3*T3^2,                                      (I)

with Q2 monic quadratic and squarefree, T3 monic cubic. T3 may have
repeated roots or share roots with Q2; neither possibility is omitted.

If there are four odd multiplicities, all four must be three, since
their sum already equals twelve. Thus

    R=Q4^3                                           (II)

for a monic quartic Q4. These alternatives exhaust the ten nonsquare
empty-fibre partitions. No genericity or distinctness assumption is made
beyond the necessarily squarefree Q2 in (I).

## Excluding a quartic cube

Let Q4=t^4+u*t^3+v*t^2+w*t+z, and let c_j=[t^j]Q4^3. Comparing
coefficients of R and Q4^3 in degrees 11,10,9 forces

    Y1=-3u/4,
    Y2=u^2/8+v/2,
    w=u^3/24-u*v/2-2/3.

Degrees eight and six then force

    Y4=(6Y1+9Y2^2-c8)/2,
    Y3=(c6-14Y1*Y2+6Y2*Y4-2)/2.

All denominators are nonzero constants. The remaining coefficients are
degrees 0,1,2,3,4,5,7. Substitution gives seven polynomial equations in
(u,v,z). Their primitive integer versions f_i are retained verbatim in
[quartic_cube-input.json](quartic_cube-input.json). The checker independently
rederives each coefficient equation and verifies that its ratio to the
retained f_i is a nonzero rational constant. It also checks all eliminated
coefficients vanish identically, so no additional equation is hidden.

The [certificate](certificates.json) gives explicit rational polynomials
lambda_i with

    sum_i lambda_i*f_i = 1.

This is a literal polynomial identity over Q, checked by FLINT
multiplication and addition. The exported seven multipliers contain 291
terms in total. It excludes every quartic cube, including nongeneric or
non-squarefree quartics. A Groebner-basis decision is not trusted in place
of this independently replayed identity.

## Centering the quadratic-cube/cubic-square case

Write Q2=t^2+u*t+v, d=u^2-4v. Since Q2 is squarefree, d!=0. Set
s=t+u/2, so Q2=s^2-d/4, and express the arbitrary monic cubic as

    T3=s^3+alpha*s^2+beta*s+kappa.

The monic polynomial part of sqrt(Q2^3*T3^2) is the following degree-six
polynomial, written using s=t+u/2:

    F=s^6+alpha*s^5+(beta-3d/8)*s^4+(kappa-3alpha*d/8)*s^3
      +(-3beta*d/8+3d^2/128)*s^2
      +(-3kappa*d/8+3alpha*d^2/128)*s
      +3beta*d^2/128+d^3/1024.

The checker verifies directly that Q2^3*T3^2-F^2 has degree at most five.
This F must equal a(t)+1/2 if (I) holds. Indeed two monic degree-six
polynomials whose squares differ by degree at most five must coincide:
factor their difference of squares as in the preceding square argument.

Put f_j=[t^j]F, and L=Q2^3*T3^2-F^2. The actual F7 formulas force

    f3=-1, 3f1+2f4=0,
    L5=-1, L4=-f5, L3=1-f4,
    L2=0, L1=-f5-f2, L0=-f4/3-1/4.                  (A)

These are exact coefficient identities, not estimates of a truncated
infinity expansion. To check the last three terms, recover
Y1=-f5/2, Y2=f4/3, Y4=Y1-f2, Y3=f0-1/2 in the defining expression for L.

The first three L equations simplify sharply. Define

    tau1=3kappa*d^2/128+alpha*d^3/1024,
    tau2=beta*d^3/1024+3d^4/32768,
    tau3=kappa*d^3/1024+3alpha*d^4/32768.

Equating the t^5,t^4,t^3 coefficients in (A), with the shift s=t+u/2,
is triangularly equivalent to

    tau1=-1/2, tau2=-u/4, tau3=1/2-u^2/8.

The checker verifies the three universal triangular identities before any
substitution. The linear coefficient determinant in (alpha,beta,kappa)
is a nonzero rational multiple of d^9, so the solution is unique for d!=0:

    alpha=2048*(24-6u^2+d)/(5d^4),
    beta=-256u/d^3-3d/32,
    kappa=-64/(3d^2)-alpha*d/24.

This is the decisive reduction: the original five factor coefficients
have become just u,d, with the sole nonconstant denominator d. No case
with squarefree Q2 has been discarded.

## The remaining five equations force d=0

Insert these three expressions in the remaining conditions (A):

    f3+1=0,
    3f1+2f4=0,
    L2=0,
    L1+f5+f2=0,
    L0+f4/3+1/4=0.

Their primitive integer numerators g1,...,g5 are retained in
[centered-quadratic-input.json](centered-quadratic-input.json). Clearing
denominators loses no solution on d!=0: every denominator is a nonzero
constant times a power of d, checked explicitly. The five degrees are
7,9,9,10,11, with 8,15,10,16,22 terms respectively.

The [second explicit certificate](certificates.json) supplies rational
polynomials mu_i such that

    sum_i mu_i*g_i = d.

The five multipliers have 130 terms in total. Independent FLINT arithmetic
replays this exact identity. Hence every common zero of the necessary
equations has d=0, contrary to the squarefreeness of Q2. This excludes (I)
globally. The fact that clearing denominators creates extraneous points
on d=0 is harmless: the argument excludes the entire open d!=0.

## Conclusion and assurance

Neither possible factorization of an empty-fibre eliminant exists.
Therefore every R(t;Y) has a simple root. That root gives an actual source
point above every c!=0 target. Every c=0 target was already counted by
GAOF7BND1. This proves that F7:C5->C5 is surjective.

Gao's degree-twelve generic fibre still proves noninjectivity, and the
previously determined nonproper locus remains present. The result does
not assert properness, injectivity, the earlier sufficient uniform gcd
bound, or any planar-Jacobian conclusion. Together with the occurrence
certificates in GAOF7CARD1, only the global cardinality three remains
undecided; the full stratum classification remains open.

The direct five-variable system for (I) timed out at thirty seconds.
Its input, empty log and run receipt are retained. The centered two-variable
system completed in about a hundredth of a second; generation timings
are observations, not proof. Likewise the quartic cube was initially a
unit-ideal decision and was promoted only after exporting and replaying
the full identity. Every run had a thirty-second cap, 1536MB address-space
cap and 1MB output cap; default verification reruns none of them.

From the repository root:

```sh
timeout 60s .venv/bin/python research/jacobian/gao-f7-surjectivity/verify.py
```

SymPy independently reconstructs the coefficient equations; FLINT replays
the Singular-produced identities. The written factorization exhaustion,
uniqueness of F, necessary-equation transport and lifting criterion explain
why the identities prove the theorem. This is not formal verification or
an independent end-to-end proof review. No external review, author endorsement,
publication or novelty priority is claimed. MATH_STATUS.json is authoritative.
