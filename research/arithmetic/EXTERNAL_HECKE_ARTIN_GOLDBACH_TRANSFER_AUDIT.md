# Hecke, Artin and Goldbach transfers

The source is the [solved-conjectures catalogue](https://github.com/1337isnot1337/solved-conjectures/tree/758066055ef49c2c66cb1087918ae62ad75e7ac1)
and its [research snapshot](https://github.com/1337isnot1337/solved-conjectures/tree/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/rh-transfer).
Here the finite class-field and character-projection statements are
unconditional. The window criterion uses the pinned abstract Landau
theorem. All deductions from a zero-free hypothesis retain that
hypothesis.

## Hecke windows including angular order and real norm twists

The [finite-order source](https://github.com/1337isnot1337/solved-conjectures/blob/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/rh-transfer/hecke-bootstrap.md)
is an application of
[Deduction.landau_window](https://github.com/1337isnot1337/riemann-hypothesis-solved/blob/9d5316ec73f5df89c4d95703ed106837b8a4060e/extracted/lean/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauTheorems/landauWindow.lean)
to [PoleData](https://github.com/1337isnot1337/riemann-hypothesis-solved/blob/9d5316ec73f5df89c4d95703ed106837b8a4060e/extracted/lean/src/Development/AutomorphicArithmetic/DedekindLFunctions/LandauZeroFreeCriterion.lean).
Both declarations were inspected at the given archive pin. The argument
extends to finite families of unitary Hecke characters, including the
Eisenstein angular family in the [local transfer audit](EISENSTEIN_HECKE_TRANSFER_AUDIT.md#uniform-prime-ideal-windows-conditional-on-the-78-theorem).

For such a character \(\chi\), write \(D_\chi=L'/L\) and set
\[
T_\chi(x)=\sum_{N\mathfrak a\le x}\Lambda_K(\mathfrak a)\chi(\mathfrak a)
            \left(\frac{x}{N\mathfrak a}-1\right).
\]
Use the finite Euler-product \(L\)-function. If it has no pole, put
\(G_\chi=T_\chi+xD_\chi(1)\). If it has the pole \(\kappa=1+i\tau\ne1\)
of a nonzero real norm twist of the principal character, put
\[
G_\chi=T_\chi+xD_\chi(1)-\frac{x^\kappa}{\kappa(\kappa-1)}.
\]
If it has a pole at \(1\), write
\(D_\chi(s)=-1/(s-1)+\gamma_\chi+O(s-1)\) and use
\(G_\chi=T_\chi-x\log x+(1+\gamma_\chi)x\).
The Mellin transform of \(T_\chi\) for \(\Re s>1\) is
\(-D_\chi(s)/(s(s-1))\), obtained by integrating each summand.
The corrections remove the artificial pole at \(1\) and the genuine
logarithmic-derivative pole at \(\kappa\), or both terms at the double
pole when \(\kappa=1\). These claims use the usual analytic
continuation and nonvanishing on \(\Re s=1\) for unitary Hecke
\(L\)-functions, after shifting the principal pole for norm twists.
The trivial and nontrivial archimedean orders do not change the
coefficient bound \(\left|\sum_{N\mathfrak a=n}
\Lambda_K(\mathfrak a)\chi(\mathfrak a)\right|\le [K:\mathbb Q]\Lambda(n)\).
Finite Euler-factor deletions add zeros only on \(\Re s=0\).

Take positive weights \(a_i\) and put
\(g(u)=\sum_i a_i\Re G_{\chi_i}(e^u)\) for \(u\ge0\), extending by
zero below zero. The conjugate-transformed expressions use
\(\overline\chi_i\) and the conjugate pole corrections. Their average
has residue
\[
-\frac{\sum_i a_i(m_{\chi_i}(\rho)+m_{\overline\chi_i}(\rho))}
       {2\rho(\rho-1)}
\]
at any zero \(\rho\) in \(\Re\rho>0\). This is nonzero. For a real zero
it is positive because \(0<\rho<1\). The coefficient bound gives
\(|g(u)|\ll e^{2u}\). These verify the analyticity, isolated-pole,
residue, growth and Laplace hypotheses of the abstract theorem.
If for some \(0<\theta<1\), \(\lambda>1\) the weighted real sum is at
most \(C x^\theta\) throughout arbitrarily large \([Y,Y^\lambda]\),
use \(U=\log Y\) and majorant \(|C|e^{\theta u}\). Its weighted
integral over \([U,\lambda U]\) is \(O(U)\), hence at most \(e^{\eta U}\)
for each fixed \(\eta>0\) at large \(U\). The theorem excludes zeros
with real part greater than \(\theta\). Taking \(\theta\downarrow1/2\)
and applying the functional equation yields the corresponding critical
line statement if those windows are supplied.

This extension admits the actual infinity types, even the principal
norm twists. It supplies no windows. The external archival
\(x^{1/2+9\varepsilon}\) bound consumes a moving cyclotomic
\(\operatorname{IsChiField}\) host and the all-prime finite Selmer-length
bound \(\operatorname{VVle}\) at \(t=\lfloor x^{401}\rfloor\) and
\(B=\exp(x^\varepsilon/2)\). Those hypotheses are absent for the
Eisenstein angular family. Even a smooth prime-ideal GRH estimate
would control a single character sum. The central shell at
\(H\le N^{2/3}\) needs joint cancellation of actual quaternionic
degree components, and the sharp energy and switching nodes require
additive correlations with their local main terms. In the switched
range \(Q\asymp N^{3/4+o(1)}\), even a GRH-scale individual character
prime error \(N^{1/2+o(1)}\) stays above the \(N^{1/4+o(1)}\) progression
scale. This is an interface bound for the direct estimate. It leaves
the central range of
[ARITHMETIC-EISENSTEIN-POWER-WINDOWS](EISENSTEIN_MOVING_WINDOWS.md),
[GAP-ARITH-FBS-SHARP-FOUR-FORM-ENERGY](PRESCRIBED_PRIME_BINARY_VARIANCE.md)
and [GAP-ARITH-M12-CRITICAL-SIEVE-REMAINDER](ASYMMETRIC_CRITICAL_SWITCHING.md)
open.

## Fixed CM splitting conditions on class-prime pools

The [CM source](https://github.com/1337isnot1337/solved-conjectures/blob/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/rh-transfer/cm-character-descent.md)
preserves odd anticyclotomic character order. The same argument gives
a useful sharpened norm statement. If \(K\) is imaginary quadratic and
\(F/K\) is CM, then
\[
\mathrm{Cl}(K)^2\subseteq N_{F/K}\mathrm{Cl}(F).
\]
Let \(H_K\) be the Hilbert class field and set \(M=F\cap H_K\).
Artin reciprocity identifies the annihilator of the norm image with
the characters of \(\operatorname{Gal}(M/K)\). In fact the pullback
corresponds to character restriction, whose image on \(G_F\) is
\(\operatorname{Gal}(H_K/M)\). Conjugation on
\(\operatorname{Gal}(H_K/K)\) is inversion by the ideal norm identity.
The intersection \(M\) is CM since it is a nonreal subfield of \(F\),
and is Galois over \(\mathbb Q\) since inversion preserves the
corresponding subgroup. Complex conjugation of a Galois CM field is
central: every embedding intertwines the CM involution, and in a common
normal closure complex conjugation thus commutes with the Galois
action. Therefore inversion acts trivially on
\(\operatorname{Gal}(M/K)\), so this quotient has exponent at most two.
This proves the square inclusion. Applying the same argument to a
finite cyclic anticyclotomic character field of order \(r\) shows that
the pullback can lose a factor dividing \(\gcd(2,r)\). Equality occurs
for the quadratic character of \(K(\sqrt2)/K\) with
\(K=\mathbb Q(\sqrt{-3})\), after pullback to \(K(\sqrt2)\).

For application, suppose \(F/\mathbb Q\) is a **fixed Galois CM field**
containing \(K\), and \(A\) is any odd quotient of the ideal class
group of \(K\). Let \(L/K\) be its class field. Then \(L\cap F=K\)
by the same inversion argument, now applied to a quotient of \(A\).
In the fiber product
\(\operatorname{Gal}(LF/\mathbb Q)\) over
\(\operatorname{Gal}(K/\mathbb Q)\), the element \((a,1)\) for
\(a\in A\) has conjugacy class \(\{(a,1),(a^{-1},1)\}\), collapsing to
one element when \(a=1\). Chebotarev therefore gives rational primes
which split completely in \(F\) and have one prime over them of class
\(a\) in \(A\), with density
\[
\frac{|\{a,a^{-1}\}|}{[F:\mathbb Q]|A|}.
\]
For a real-rational congruence instance take
\(F=K\,\mathbb Q(\zeta_M)\). All the selected primes can satisfy
\(p\equiv1\pmod M\). This extends
[EXT-ODD-IDEAL-CLASS-PRIME-POOLS](PRIMITIVE_NORM_ODD_CLASS_GROUPS.md#existence-of-the-finitely-many-class-primes):
when the full class group has odd order, the fixed reserved primes of
every absorber block can additionally split in \(F\). The original
pair admissibility and counting proof then uses these new finite
pools with exact exponent one. Only the reserved primes have the
extra splitting condition. There is no restriction here on the
Frobenius of the residual cofactor primes.

The analytic descent in the external source is distinct. For abelian
\(F/K\), Artin induction gives
\(L_F(s,\chi\circ N)=\prod_{\psi\in\widehat{\operatorname{Gal}(F/K)}}
L_K(s,\chi\psi)\). Every factor is holomorphic away from a possible
pole at \(1\) of the principal factor. Thus the hypothesized strip
\(H(\delta)\) on cyclotomic hosts descends to finite-order characters
over abelian base fields, with that pole allowed. Odd-order
preservation is unnecessary for zero descent away from \(1\) when the
host hypothesis includes trivial characters, but is essential for
the odd class-pool argument above.

## Artin’s weighted family and its residual correlation

The [uniform splitting](https://github.com/1337isnot1337/solved-conjectures/blob/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/rh-transfer/artin-uniform-splitting.md)
and [weighted average](https://github.com/1337isnot1337/solved-conjectures/blob/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/rh-transfer/artin-average.md)
arguments concern a fixed integer base \(a\) with \(|a|>1\), neither\na square nor \(-1\), and a growing prime \(q\).
Their separate input \(H(\delta)\) says that all finite-order Hecke
\(L\)-functions over cyclotomic fields containing \(\mu_{12}\) have
no zeros in \(\Re s>1-\delta\), allowing the pole at one. The
[OpenAI source](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/Primitive-roots-for-every-admissible-integer-base-October-4-2026/build/sections/01-introduction.tex)
claims \(\delta=10^{-6}\). The preceding theta estimate for that
claim remains outside this audit.

For \(K_q=\mathbb Q(\mu_q,a^{1/q})\), the degree is \(q(q-1)\)
for all large prime \(q\). A fixed valuation of \(a\) at a prime
\(\ell\mid a\) stays prime to \(q\) in \(\mathbb Q(\mu_q)\) for
large \(q\), since \(\ell\ne q\) is unramified there. The valuation
of a \(q\)-th root has denominator \(q\), which forces relative
degree \(q\). The relative discriminant divides the discriminant of
\(X^q-a\), so \(\log|D_{K_q}|\ll_a q(q-1)\log q\). Adjoin
\(\mu_{12}\): the enlarged field is cyclic Kummer over the cyclotomic
host, and abelian over \(K_q\). The two zeta factorizations descend
\(H(\delta)\) to \(\zeta_{K_q}\). In the
[source smooth explicit formula](https://github.com/openai/math/blob/4eaee87548250b2682affda7d05d947c05554bae/preprints/Primitive-roots-for-every-admissible-integer-base-October-4-2026/build/sections/07-splitting.tex),
divide the resulting smoothed prime-ideal estimate by
\(q(q-1)\log x\). Each split rational prime in \((x,2x]\) contributes
at least this amount before division. The bound is
\(O_a(x/(q(q-1)\log x)+x^{1-\delta}\log(2q)/\log x)\).
Summing this over \(y<q\le x^\eta\), \(\eta<\delta\), saves the small
tail. For \(q>\sqrt x(\log x)^A\), each candidate prime divides
one of \(a^m-1\) at \(m\ll\sqrt x(\log x)^{-A}\). The log of their
product is \(O_a(x(\log x)^{-2A})\) and each distinct counted prime
contributes \(\log x\). The latter bound is unconditional.

Write \(N_q\) for the split-prime count on \((x,2x]\) and \(A_q\)
for the count with \(p\equiv1\pmod q\). In the affine Galois group
\(C_q\rtimes\mathbb F_q^\times\), the augmentation character of the
permutation action takes values \(q-1,-1,0\) at the identity, a
nonidentity translation, and a nontranslation. Hence its prime sum
\(T_q\) has exact normalization \(qN_q=A_q+T_q\). In
\(I=(x^\eta,\sqrt x/(\log x)^B]\), Brun–Titchmarsh gives
\(\sum_{q\in I}A_q/q=o(x/\log x)\), so
\(\sum|T_q|/q=o(x/\log x)\) is equivalent to a negligible middle
splitting tail. The squared-character identity gives diagonal mass
\(\sum_{q\in I}((q-2)N_q+A_q/q)\) in \(\sum|T_q|^2/q\).
Each \(p-1\) has \(O_\eta(1)\) eligible prime divisors so this
diagonal is \(O_\eta(\sqrt{x}(\log x)^{-B}\pi(2x))\), below the
sufficient \(o(x^2/(\log x)^2)\) second-moment scale. The off-diagonal
remains. The Artin factor has multiplicity \(q-1\) in
\(\zeta_{K_q}=\zeta_{\mathbb Q(\mu_q)}L_q^{q-1}\), so normalization
by the degree gives weight \(1/q\) on its zeros. Dirichlet GRH
controls the cyclotomic term in the external signed zero-sum
reduction and does not supply the nonlinear family saving.

The diagonal mechanism is an exact template for a *family whose
actual Frobenius indicators are those of this affine group*. None
of the three additive DAG targets above has that identification.
A pair or four-form sum would also have to retain the joint local
main term for its actual coefficients.

## Exact Goldbach projection and the shrinking Hardy test

The [Goldbach branch note](https://github.com/1337isnot1337/solved-conjectures/blob/dab1ac8d1ecea8dc34c38932fcd2a4b3f56ce0ef/research/rh-transfer/goldbach-spectral-bridge.md)
identifies the norm gap in
[Li–Liu](https://mrc.sdu.edu.cn/__local/6/14/8E/28753731EE31D7635BAF60732C2_0C644604_2CE958.pdf),
§§2.4–2.5 and Theorems 6.1–6.2. The unconditional finite projection
applies to general real weights. For \(w\) supported in
\([N/3,2N/3]\), and \(N>2\), character orthogonality gives
\[
\sum_{\substack{n+m=N\\(n,N)=1}}w(n)w(m)
=\sum_{(n,N)=1}w(n)^2-\frac{2}{\varphi(N)}
    \sum_{\chi(-1)=-1}\left|\sum_n\chi(n)w(n)\right|^2 .
\]
Indeed the odd-character kernel is the difference of the unit-class
equal and opposite indicators. On the support, those congruences
are respectively \(m=n\) and \(m+n=N\). For nonunits, partition
pairs by \(d=(n,N)=(m,N)\) and apply the same identity to
\(w(d\,\cdot)\) at modulus \(N/d\). The case \(N/d=2\) has only its
diagonal term, and the opposite and equal indicators coincide.
This gives an exact projection on the centered semigroup pair
count, including all valuations. The logarithmic-derivative
arithmetic kernel and its beta model still use von Mangoldt
weights. Replacing those weights requires a new residual identity
and estimate.

There is an independent width cost in the Hardy norm. For the
original von Mangoldt operator \(T_N\) under \(\mathrm{GRH}(N)\),
take a fixed nonzero real even \(b\in C_c^\infty(-1,1)\), with
\(0\le b\le1\), and set
\(\phi_\Delta(u)=b((u-\tfrac12)/\Delta)\) with
\(\Delta=H/N\to0\). Proposition 7.4 applies for each such test.
The norm formula (7.15) gives
\[
\|v_{N,\Delta}\|^2
=\frac1{2\pi}\int
\frac{(Nu)^{2/\log N}|\phi_\Delta(u)/2-u\phi_\Delta'(u)|^2}
 {1+(Nu)^{-2(1/8-1/\log N)}}\,du
\sim \frac{e^2}{8\pi\Delta}\int |b'(y)|^2\,dy .
\]
For the asymptotic substitute \(u=\tfrac12+\Delta y\). On the
support the quotient of the displayed weights tends uniformly
to \(e^2\), and the derivative-square term gives the leading mass.
The cross term is \(O(1)\) and the undifferentiated term is
\(O(\Delta)\). Meanwhile \(c_{\phi_\Delta}=\Delta\int b^2\)
exactly. Thus a norm-only route to relative asymptotics at width
\(H\) needs a bound at the scale
\[
\|T_N\|=o(\Delta^2\mathfrak S(N)).
\]
Even the paper’s proposed Hypothesis R scale
\(\|T_N\|\le N^{-1/2}(\log N)^{3/2}\sqrt{\mathfrak S(N)}\)
would only certify fixed powers \(H=N^\theta\) with
\(\theta>3/4\) by this norm bound, since \(\mathfrak S(N)\) is bounded
below and \(O(\log\log N)\). The existing Eisenstein central bound
already reaches \(\theta>2/3\) for a different arithmetic weight.
A direct estimate of the one moving quadratic form, which can be
smaller than its operator norm bound, would still require new
arithmetic cancellation.

The spectral zero equivalence is algebraic once the operator is
given. The factor from eigenvalue \(\mu\ne0\) has zeros
\(q=\tfrac12\pm\tfrac12\sqrt{1-\delta_N^2/\mu^2}\). This verifies
the \(\|T_N\|\le\delta_N\) scale for Hypothesis R and the necessity
of an extra quantitative assumption in the spectral deduction.
Li–Liu’s softer Hypothesis S asks for
\(\|T_N\|\le\mathfrak S(N)/(2h(N)\sqrt{1-\eta_N})\) with
\(h(N)\to\infty\) and \(\eta_N=h(N)^2(\log N)^3/(N\mathfrak S(N))\).
The source's GRH-only \(O((\log N)^2)\) norm estimate does not
supply it. A later [Yilin Liu preprint](https://hexagonmath.org/pdf/2611.00788v1)
claims, under \(\mathrm{GRH}(N)\) at many moduli, a lower bound
exceeding the Hypothesis R scale by arbitrary fixed powers of
\(\log N\). The preprint's test compares \(N-m\) with \(N+m\)
after cancellation of the continuous model. Its full arithmetic
construction is outside this audit and is unused in the width
calculation.

The external algebraic and analytic interfaces were checked at
the displayed pins. Neither the imported Lean archive nor the
full independent analytic manuscripts are replayed or certified
here. The new class-pool restriction, norm-cokernel refinement,
finite projection with gcd strata, and Hardy width calculation have
the direct proofs above.
