# Exact coalescing-binomial saddle and amplitude bridge

## Scope

This note proves local algebraic and analytic identities for the ACTUAL consecutive factorial-moment pair. It does not prove all-gap Strong Factorial, all lower exponents at gap two, or a uniform coalescing-saddle tail. The missing global and uniform estimates are stated explicitly in Section6. The local identities are not a transfer from arbitrary positive measures or Gamma convolutions.

Let N>=2 be an integer, let a>0 be real and let d>=1 be an INTEGER; put D=a+d and g(s,w)=s^a(w+s^d). For integer a,d these are the actual binomial phases. On any local logarithm near s!=0,g!=0 put Phi=log(g)-s. The actual normalized pair on the common N scale is I_N(w)=int_0^infinity exp(NPhi)ds and J_N(w)=int_0^infinity exp(NPhi)/g ds. When a is an integer, all deformation questions refer to the original entire polynomial integrands, not to an unrestricted logarithmic contour. For real noninteger a, only the stated local-branch geometry and the positive-axis identities are asserted; the origin is a branch point.

## 1. Rational stationary map and its critical-value circle

**Lemma CoalescingSaddleCriticalValueCircle.** The stationary equation is w=Omega(s)=s^d(D-s)/(s-a), and

Omega'(s)=-d*s^(d-1)*Q(s)/(s-a)^2,
Q(s)=s^2-(D+a-1)s+aD.

Assume a>(d-1)^2/4. Put h=(d-1)/2, m=a+h and v^2=a-h^2>0. The two nontrivial critical points are s_plus=m+iv and s_minus=m-iv. They satisfy

|s|^2=aD, |s-a|^2=a, |D-s|^2=D.

Consequently their critical values have the exact common modulus

|Omega(s_plus)|^2=|Omega(s_minus)|^2=a^(d-1)D^(d+1).

Proof. Differentiate the displayed rational map and factor its numerator. The discriminant of Q is(d-1)^2-4a. The three modulus identities are direct quadratic identities in a,h. Substitution in |s|^(2d)|D-s|^2/|s-a|^2 gives the result. Neither the circle nor the two points are an approximate saddle calculation.

The critical values themselves can be computed without complex square-root approximations. In the quadratic algebra Q(a)[s]/(Q),

1/(s-a)=(D-1-s)/a,
Omega(s)=s^d[D(d-1)-d*s]/a.

Let L_0=2,L_1=D+a-1,L_k=(D+a-1)L_(k-1)-aD L_(k-2). The trace of the critical values is [D(d-1)L_d-dL_(d+1)]/a and their product is a^(d-1)D^(d+1). Thus they are the roots of an explicit rational quadratic.

For positive INTEGERS a,d the critical values are algebraic integers, not merely roots of a rational quadratic. Reduce the numerator s^d[D(d-1)-d*s] modulo a and Q: s^2=(d-1)s modulo a, so both reduced numerator coefficients vanish modulo a. Division by a leaves an integer linear polynomial in the integral element s. Hence its trace is integral, its norm is the integer displayed above, and in the nonreal regime its exact minimal polynomial is an irreducible imaginary quadratic.

For gap d2 this specializes to

wcrit=(1-4a)*s+3a(a+2),
trace(wcrit)=-2a^2+10a+1,
norm(wcrit)=a(a+2)^3,
trace^2-4norm=(1-4a)^3.

For d>=2 their arguments never collide on the real axis inside this regime. Along the upper critical point, the derivative in a is

d/da log Omega(s_plus(a))=d/s_plus(a).

The stationary derivative term vanishes because Omega'(s_plus)=0. Its imaginary part is negative. For d>1 an explicit continuous argument, with every arctangent in(0,pi/2), is

d*atan(v/m)-atan(v/(h+1))-atan(v/h).

This expression is strictly greater than-pi, starts at0 at the discriminant boundary and tends to-pi as a increases to infinity. Its derivative is Im(d/s_plus)<0. Thus it stays strictly between-pi and0, without any unproved winding assumption. Hence one critical value is strictly in the lower half-plane and the other in the upper half-plane. For d1 the explicit values are1-a plus/minus2i*sqrt(a), which give the same conclusion.

## 2. Simple Airy turning points and the quartic boundary

At a stationary point, g=d*s^D/(s-a), so it is nonzero at the two points in Section1. Also -Phi_ss=Q(s)/s^2 along the stationary curve. Differentiating along that curve at a critical point, where Omega'=0, gives

Phi_sss(s,wcrit)=-Q'(s)/s^2=-(2s-D-a+1)/s^2.

It is nonzero in the strict coalescing regime, and

Phi_sw(s,wcrit)=-(s-a)^2/(d*s^(d+1))!=0.

These are the local simple-turning-point/nontrivial-unfolding hypotheses for a cubic CFU coordinate. This identifies the local type; it does not identify the global integration cycle.

The equality boundary is a different singularity. If d>1, a=h^2,d=2h+1 with h>0, then D=(h+1)^2, s0=h(h+1), and

w0=h^(2h)(h+1)^(2h+2)>0.

There Phi_s=Phi_ss=Phi_sss=0, while Phi_ssss=-2/s0^2!=0. Indeed Q has a double zero, Omega' and Omega'' vanish there, and a second differentiation along the stationary curve gives -Phi_ssss=(Q/s^2)''=2/s0^2. Thus an eventual theorem crossing this boundary would need a quartic/Pearcey treatment, not just two simple Airy turning points.

## 3. The actual two amplitudes remain independent at confluence

**Lemma ActualPairAiryJetNondegeneracy.** At every stationary point, (1/g)_s=-1/g. In any local cubic coordinate

Phi(s(v,w),w)=Phi0(w)+v^3/3-zeta(w)*v,

centered at a simple turning point v0=zeta0=0, put b=s_v(0,wcrit). Then b^3=2/Phi_sss. The two transformed amplitudes for the actual pair are

p(v)=s_v(v), q(v)=s_v(v)/g(s(v),w).

Their zero/first-jet determinant is

p(0)q'(0)-p'(0)q(0)=-b^3/g(scrit)=-2/[g(scrit)*Phi_sss(scrit)]!=0.

Proof. Phi_s=0 means g_s/g=1, which proves the first assertion. Differentiate the cubic normal form three times at its double stationary point to get Phi_sss*b^3=2. Applying the product rule to q gives the determinant. This exact algebra shows that the limiting Airy/Airy-derivative coefficient matrix is invertible; the unnormalized difference1/g(s_plus)-1/g(s_minus) vanishing at confluence is not itself an obstruction.

This lemma is local. In particular it does NOT bound its inverse or error terms uniformly when a,d grow.

## 4. Exact radial transversality of the Airy parameter

**Lemma RadialAiryParameterTransversality.** In the same cubic coordinate,

zeta_w(wcrit)*wcrit=b/d.

Indeed at the double saddle Phi_sw=-zeta_w/b. Multiply by the explicit stationary value Omega(s)=s^d(D-s)/(s-a), and use(s-a)(D-s)=s, which is equivalent to Q(s)=0. This gives the formula.

For d>=2 the radial direction is uniformly transverse to the three standard rotated Airy zero rays. Here is a purely algebraic angle bound. The exact identity

m^2-(d^2-1)v^2=[a-d(d-1)/2]^2>=0

implies theta=arg(s_plus)<=arcsin(1/d)<=pi/6. Moreover

Phi_sss(s_plus)=-2iv/s_plus^2,
arg(Phi_sss)=-pi/2-2theta.

Every cube root b of2/Phi_sss is therefore separated from any ray whose cube is negative real by angle at least pi/18. Those rays are the zero rays of the standard three rotated Ai functions. This is an angle statement conditional on identifying the actual Airy cycle; it does not assert that every Airy linear combination has zeros on those rays.

## 5. The exact phase shift and the uniformity obstruction

The weight1/g need not be treated as a slowly varying amplitude. For n=N-1, substitution s=(n/N)r gives the exact identity

J_N(w)=(n/N)^(Dn+1)*I_n((N/n)^d*w).

This changes BOTH actual order and radial parameter. It is not a fixed-Gamma-shape replacement.

**Corollary ActualConsecutivePhaseShiftEquivalence.**
For each fixed positive integers a,d and integer N>=2, put n=N-1. Then
J_N(w)=0 if and only if I_n((N/n)^d*w)=0, and the original actual pair
I_N,J_N has a common zero w if and only if I_N(w)=I_n((N/n)^d*w)=0.
The scalar (n/N)^(Dn+1) is nonzero. For a nonzero binomial with AB!=0,
the actual factorial moments are B^N*N^(DN+1)*I_N(w) and
B^n*N^(Dn+1)*J_N(w), with A/B=N^d*w; both scalar factors are nonzero.
The monomial cases have nonzero factorial moments. Thus the equivalence
concerns the two ACTUAL coupled shapes, not a fixed-shape replacement.
This identity does not prove either global noncollision assertion.


For FIXED a,d as N tends to infinity near a simple turning point, the radial change is asymptotic to d*wcrit/N. By Section4 its leading effect on the normalized Airy argument is asymptotic to b*N^(-1/3). No such joint-parameter asymptotic is asserted without separate derivative control. It is transverse to the selected Airy zero rays, with the angle bound above. This suggests using the two separately shifted Airy representations of I_N and I_(N-1), rather than a nonuniform Taylor expansion of1/g.

The uniformity problem is real: in the natural large-a coordinate t=(s-a)/sqrt(a), the log derivative of1/g is -sqrt(a) at a stationary point. Across an Airy scale t=O(N^(-1/3)) its variation can be exp(O(sqrt(a)*N^(-1/3))). Thus fixed-amplitude CFU estimates alone do not give uniform error constants at N comparable to a/d. The exact phase-shift identity avoids assuming that false bound, but does not itself supply the needed shifted-root estimate.

## 6. Precise remaining obligation, with gap two as first test

A proof in the coalescing regime must still establish: a root-localizing contour for every relevant w; its decomposition into the actual right-sector Airy cycles with controlled endpoints/other saddles; a uniform error bound for the TWO orders at the radially shifted arguments; and a quantitative Airy-vector or adjacent-root separation estimate that dominates those errors. The local jet determinant alone does not prove any of these global claims.

The sharp initial target is fixed gap d2 for every lower exponent. The [retained general lower-exponent theorem](GAMMA_PRODUCT_ROOT_MODULUS.tex) already covers a>=2(2N-3), hence N<=a/4+3/2. Least exponents1,2 are already closed at every gap, and the retained mixed-adjoint theorem closes a3,d2 for every positive Gamma shape. The retained [lower-four](LOWER_FOUR_QUADRATIC_SECTOR_SFC.tex) and [lower-five](LOWER_FIVE_QUADRATIC_SECTOR_SFC.tex) quadratic-sector theorems also settle their respective families at every positive real Gamma shape. The genuinely remaining ordinary family is therefore

d2, a>=6, N>(a+6)/4.

A rigorous uniform shifted-Airy separation theorem in this domain would close an entire exponent-gap family, not merely another finite set of supports. No claim that it has been proved is made here.

## 7. Actual classical antecedents and scope barriers

[DLMF2.4(v)](https://dlmf.nist.gov/2.4#v), equations2.4.17-21, gives the cubic CFU change and explains the contour-dependent Airy/Scorer alternatives. [DLMF36.12(ii)](https://dlmf.nist.gov/36.12#ii), equations36.12.11-12, gives the two-amplitude Airy/Ai-prime uniform approximation. These provide the local normal-form framework, not the explicit global two-order bounds in Section6. The Airy no-common-zero fact follows directly from y''=z*y and uniqueness: y=y'=0 at one finite point forces the solution to be identically zero.

Dunster--Gil--Ruiz-Antolin--Segura, [ETNA65(2026)140-155](sources/dunster_gil_ruiz_segura_etna65_2026_140_155.pdf), treat ordinary generalized Bessel polynomials (_2F0). Section2 condition(2.4) keeps(a_Bessel-2)/(n+1/2) in a bounded interval. Our gap-two original polynomial is _3F0. The retained cubic-adjoint polynomial is an ADDITIVE convolution of two _2F0 factors; its first factor has a_Bessel/n approximately a-1, unbounded here. The published theorem therefore does not directly imply our uniform actual-pair claim.

Wolfs, [arXiv2401.08312v1](sources/wolfs-2024-multiple-orthogonal-2401.08312v1.pdf), Theorem4.6 matches the generalized hypergeometric class, but its Section4 works on the unit circle rather than an AT system and does not establish the required general complex-zero separation.

The older retained [Gamma-convolution obstruction](GAMMA_CONVOLUTION_APPELL_OBSTRUCTION.tex) already gives positive Gamma convolutions whose ordinary and normalized-cumulant Hankel matrices are BOTH strictly totally positive while consecutive Appell polynomials share a nonreal quadratic. Thus a generic positivity or total-positivity transfer would be false. Any proof must use the single Gamma-power density and the exact aN+1/a(N-1)+1 coupling.


## Exact evidence

Run from a checkout:

    python3 research/factorial/image-rigidity/verify_coalescing_binomial_saddle_bridge.py --repo-root . --output /tmp/coalescing_saddle_review.json

The checker and [receipt](coalescing_binomial_saddle_bridge_verification.json)
support the explicit finite controls. The [shared source ledger](coalescing_binomial_local_source_ledger.json)
pins both written accounts, both checkers and their retained/context sources.
The uniform assertions use the written arguments. All evidence and review are
internal; no external human review, formal verification or historical priority
is asserted. Retained20November2026.
