# Exact finite identity certificates at arbitrary singular initial points

Internal lemma packet, 20 November 2026. The statements concern formal
power series over characteristic zero. No regular-singular or ordinary-
point hypothesis is imposed. The application to sphere periods requires
a VERIFIED annihilator; matching finitely many coefficients is not a
substitute for that premise.

## 1. The lowest-valuation indicial polynomial

Let M=sum_(j=0)^r p_j(z)D^j be a NONZERO differential operator, with
p_j inQ[z], D=d/dz. Omit the zero p_j. Set

    v_j=ord_0 p_j, c_j=[z^(v_j)]p_j,
    s0=min_j(v_j-j),
    I_M(N)=sum_(v_j-j=s0) c_j*(N)_j,                 (I1)

where (N)_j=N(N-1)...(N-j+1) and(N)_0=1.

**Lemma SingularFormalInitialCutoff.** The polynomial I_M is nonzero.
If a nonzero formal series H with M(H)=0 has order N>=0, then I_M(N)=0.
Consequently, if

    B(M)=max({0} union {N>=0 integer:I_M(N)=0}),       (I2)

then M(H)=0 together with [z^n]H=0 for0<=n<=B(M) proves H=0.

Proof. The largest participating j contributes the unique highest
polynomial degreej in(I1), with nonzero coefficient c_j. Thus I_M is
nonzero. Write H=h_N z^N+... with h_N!=0. The coefficient of
z^(N+s0) in M(H), viewed in the Laurent-series field if necessary,
is exactly h_N I_M(N). Every later H-term and every later p_j-term
has strictly larger total exponent. Falling factorials already handle
j>N, where the derivative of z^N vanishes. Therefore I_M(N)=0.
The finite initial-segment conclusion follows immediately. QED.

The finite set in(I2) is computable exactly. Clear rational denominators
in I_M. If it is constant, there are no roots. Otherwise the elementary
Cauchy bound1+max_(i<degI)|I_i/I_top| bounds every complex root's
absolute value; check the finitely many nonnegative integers up to its
floor by exact rational evaluation. This computes B(M) without a root
finder, algebraic-number approximation or numerical tolerance.

The bound depends on the ACTUAL operator coefficients. For example
M=zD-N has I_M(s)=s-N and nonzero solution z^N. Fixed differential
order and fixed coefficient degrees alone cannot justify replacing
B(M) by a constant independent of the coefficient N. This example does
not rule out separately proved bounds in a phase-degree parameter.

## 2. Homogenizing the comparison with the constant series

**Lemma ConstantComparisonHomogenization.** Suppose a known polynomial-
coefficient operator L!=0 annihilates G inQ[[z]], G(0)=1. To test G=1,
put H=G-1 and q=-L(1), which is just minus L's coefficient of D^0.
If q=0 take M=L. Otherwise put

    M=(qD-q') composedwith L.                        (I3)

Then M!=0 and M(H)=0. The coefficient-specific test of Lemma(I1)
therefore proves G=1 exactly by checking its moments1 throughB(M).
If B(M)=0 the known constant term suffices.

Indeed L(H)=q, and(qD-q')(q)=0. The leading differential coefficient
of(I3) is q times the leading coefficient of L, hence nonzero. If
L initially has coefficients inQ(z), clear a common polynomial
denominator before defining q. No singularity created by this process
invalidates Lemma(I1). In coefficient form,

    M_j=q*p_(j-1)+q*p_j'-q'*p_j,                     (I4)

with absent p_j interpreted as zero.

## 3. Exact operator fixtures

The checker includes delayed solutions z^N to zD-N, repeated-root and
irregular examples, and two actual sphere-period operators. Each
example keeps annihilator validity distinct from initial coefficient
verification.

For p=t, normalized sphere moments are0 at odd m and1/(m+1) at evenm.
Thus G(z)=atanh(z)/z. The exact formal derivative identity

    (zD+1)G=1/(1-z²)

is checked termwise and proves the homogeneous operator

    L=z(1-z²)D²+(2-4z²)D-2z

annihilates G. Taking M=(zD-1)L, which differs from(I3) by a nonzero
scalar, gives

    I_M(N)=N(N-2)(N+1), B(M)=2.

The first moment is0 but the second is1/3, so this verified finite
certificate detects nonpurity at the singular point z0.

For the genuine three-height phase

    p=t/x*(1+x²-t4)=t*x+2t*y*(1+t²),

constant-term extraction gives G(z)=J(4z²), where

    J(s)=integral_0^1 [1-s(t²-t6)]^(-1/2) dt.

The companion period packet supplies a verified order-four L for J,
with coefficients in ascending differential order

    c0=-12(4s²-84s+7),
    c1=-6(232s³-3612s²-14s-147),
    c2=-12s(232s³-3052s²-94s+2121),
    c3=-8s²(152s³-1764s²-378s+5103),
    c4=-16s³(2s-21)(4s²-27).

With q=-c0, formula(I3) gives

    I_M(N)=-10584*N(N-1)(2N-1)(6N-7)(6N+1).

Its nonnegative INTEGER roots are exactly0,1. Therefore one J-coefficient,
or moments through order2 of G, decides comparison with1 using this
annihilator. J'(0)=2/21, so Lambda(p²)=8/21 and the phase is not pure.
The annihilator's validity is an imported EXACT certificate from the
companion period proof; the indicial calculation alone cannot establish it.

## 4. Exact recognition of polynomial solutions and fixed-marker tails

**Lemma FixedOperatorPolynomialityDecision.** Let L!=0 inQ[z]<D> be a
VERIFIED annihilator of a formal series A(z). Omit zero coefficients
p_j and put

    s_infinity=max_j(deg(p_j)-j),
    I_infinity(N)=sum_(deg(p_j)-j=s_infinity)
                         lc(p_j)*(N)_j.             (I5)

This polynomial is nonzero by the largest-j argument. Every NONZERO
polynomial solution has degree N among its nonnegative integer roots:
the coefficient of z^(N+s_infinity) in L(A) is its leading coefficient
times I_infinity(N). All lower A-terms have strictly smaller exponents.
Falling factorials also cover j>N.

If the root set is empty, the only polynomial solution is0, so test
A=0 by Lemma(I1). Otherwise let Nmax be its largest root and compute

    P(z)=sum_(m=0)^Nmax [z^m]A*z^m.

Then A is a polynomial if and only if A=P. This last equality is
an exact zero test: put q=-L(P); use M=L if q0 and otherwise
M=(qD-q')L. The same proof as(I3) gives M(A-P)=0 and M!=0. Lemma(I1)
therefore supplies a finite certificate. In the empty-root case P0,
the constant coefficient of A-P is not automatically zero and must
also be checked. This is a coefficient-SPECIFIC polynomiality test,
not a uniform degree bound in a family of operators.

In the following n-pair identities, Lambda denotes normalized average on
the complex unit sphere; at n2 it agrees with the Hopf functional above.
For any fixed balanced original marker G of bidegree(e,e), the weighted
sphere series A_G(z)=sum Lambda(g*p^m)z^m has a rational-period
annihilator by the companion global proof. Its coefficients eventually
vanish exactly when A_G is polynomial. The normalized weighted period
must use the UNWEIGHTED volume denominator I_1(0)!=0, not I_G(0): a
marker such as x may have mean0. The factorial identity

    E_n(GF^m)=(dm+e+n-1)!/(n-1)! * Lambda(g*p^m)

transfers eventual vanishing for fixed bihomogeneous F,G.

ANY fixed original polynomial marker can be handled while preserving
cancellations. Decompose G into bidegrees(u,v). Terms with u>v give zero
under E_n after multiplication by F^m. For each delta=v-u>=0 group
G_delta=sum_(v-u=delta)G_(u,v). Its output has homogeneous z_out-degree
delta. For every multi-index beta with |beta|=delta,

    E_n(xi^beta G_delta F^m)
       =partial_z_out^beta E_n(G_delta F^m)
       =beta!*[z_out^beta]E_n(G_delta F^m).            (I6)

For this fixed delta,beta put J_v=xi^beta G_(v-delta,v), a balanced
marker of degree(v,v). Let e0 be the least occurring v and let A_v(z)
be the angular weighted series of J_v. With theta=zD_z, form

    B_beta(z)=sum_v product_(k=0)^(v-e0-1)
                       (d*theta+e0+n+k) A_v(z).       (I7)

The empty product is1. At orderm the multiplier is exactly
(dm+v+n-1)!/(dm+e0+n-1)!. Consequently(I6) equals

    [(dm+e0+n-1)!/(n-1)!]*[z^m]B_beta(z).

The common factor is nonzero, so its presence does not affect tail
vanishing. Every same-output-degree cancellation is retained in(I7).
Different delta give different homogeneous output degrees and cannot
cancel. There are finitely many delta,beta, so eventual vanishing of
E_n(GF^m) is equivalent to polynomiality of ALL these B_beta.

The Euler operators can be applied to the corresponding rational tube
integrands and summed BEFORE computing an annihilator. The result is
still one rational period; no additional D-finite closure theorem is
required. This decides every FIXED polynomial marker, including mixed
bidegrees. It does not decide safety quantified over EVERY possible G.
For example with n2,F1,G=R-2, the separate weighted series of R and-2
are nonpolynomial but their aggregate is identically zero. Testing the
pieces separately would give the wrong answer.

Two actual weighted sphere fixtures make the zero-mean case explicit.
For p=y,g=x, torus weight gives A_G=z/3, so the tail is zero although
A_G(0)=0. The annihilator zD-1 gives Nmax1 and exact comparison with z/3.
For p=t,g=t,

    A_G=(atanh(z)/z-1)/z=z/3+z³/5+... .

Let L_t be the verified operator of Section3. The homogeneous operator
M_A=(zD-1)(L_t composedwith multiplicationbyz) annihilates A_G, because
L_t(1+zA_G)=0 and L_t(1)=-2z. Its infinity polynomial is

    I_infinity=-(N+1)(N+2)(N+3),

so only the zero polynomial could solve it. Its origin indicial is
(N-1)(N+1)(N+2); the coefficient [z]A_G=1/3 certifies that A_G is NOT
polynomial. Thus this fixed marker survives infinitely often.

The same algebraic statements hold over any characteristic-zero field.
For a finitely generated effective coefficient field, integer-root
enumeration reduces to rational polynomial candidates by clearing
parameter denominators and taking one nonzero standard-monomial
coordinate, as detailed in the companion global proof.

## 5. Application boundary

For rational-coefficient sphere polynomials, an effective Picard--Fuchs
or regular-certificate creative-telescoping construction provides L.
Lemmas(I1)--(I3) then make identity with1 decidable by finitely many
EXACT moments, even when z0 is irregular or the leading differential
coefficient vanishes there. The cutoff is specific to the computed
operator. No uniform phase-degree bound is asserted in this lemma
packet, and no arbitrary coefficient guessing is treated as a proof.

Application proof and certificate: [global sphere recognition](GLOBAL_SPHERE_PURITY_RECOGNITION.md). Exact support: [checker](verify_formal_indicial_purity_certificates.py) and [receipt](formal_indicial_purity_certificates_verification.json). The unrestricted formal lemmas are proved above; finite controls and internal replay do not establish external or proof-assistant verification.
