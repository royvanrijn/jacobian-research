# Eisenstein Hecke transfer audit

The full unweighted variance theorem for \(C_q\) extends from
\(\theta>9/10\) to every fixed \(\theta>5/6\), with its original mean-square
error and sharp four-form energy constant. The proof uses classical
Dirichlet convexity and zero density. Winarski's attributed \(7/8\)
half-plane gives a conditional prime-ideal estimate uniform in angular
order and finite modulus. The central quaternionic window threshold and
the critical \(C_{12}+C_{12}\) switching estimate require separate inputs.
[Revised variance proof](PRESCRIBED_PRIME_POWER_WINDOW_VARIANCE.md).

## Sources and verification scope

- Winarski, **2610.00716v1**, Theorem 1.2 and §§2, 4–8. The stated boundary
  is strict. The family includes all integral angular orders compatible
  with the six units, primitive and imprimitive finite parts, and real
  unitary norm twists. The reported verification concerns the repaired
  import closure at commit 79830bf6a346a6742fde1ace977a8bb805e638e7.
  [Versioned manuscript](https://hexagonmath.org/pdf/2610.00716v1),
  [frozen release](https://github.com/TheRealDubb/angular-hecke-7-8-xenoproof/releases/tag/v7-8-2026-10-30).
- Winarski, **2610.00458v1**, Theorem 1.2, Proposition 4.1 and
  Proposition 6.3. The reported Lean evidence is attached to the frozen
  release at commit 67bb5d2df315fe5ebe7b6ffe0eb4529ffef225de.
  [Versioned manuscript](https://hexagonmath.org/pdf/2610.00458v1),
  [frozen release](https://github.com/TheRealDubb/angular-hecke-15-16-xenoproof/releases/tag/v15-16-2026-10-18).
- Bhattacharjee–Bhattacharya, **2610.00256v3**, Theorems 1.1, 1.5 and 2.3,
  Proposition 3.7, Corollary 3.8 and Appendix B. Theorem 1.5 uses the
  OpenAI \(11/12\) finite-order input reproduced in Appendix B. The authors
  report a separate Lean check using OpenAI's \(7/8\) finite-order source.
  [Versioned manuscript](https://hexagonmath.org/pdf/2610.00256v3),
  [source and reported verification archive](https://doi.org/10.5281/zenodo.23672710).

The arguments below check the interfaces against the manuscript statements.
Formal verification relies on the authors' reports. The three preprints
remain attributed sources, and the unconditional variance proof uses only
classical inputs.

## Full-count variance and completed rows

For \(f=\mathbf 1_{C_q}\), write \(r_{q,H}(N)\) for the ordered binary
count with \(|a-b|\le H\), \(H=X^\theta\). Let
\(\mu_X=c_q(\log X)^{1/\varphi(\operatorname{lcm}(2,q))-1}\) and let
\(S_q(N)\) be the exact local product. The required estimate is a centered
additive second moment of \(r_{q,H}(N)-H\mu_X^2S_q(N)\). The equivalent
sharp upper energy is
\[
 \sum_{X<N\le2X}\sum_{\substack{|2a-N|\le H\\ |2b-N|\le H}}
 f(a)f(N-a)f(b)f(N-b)
 \le (K_q+o(1))XH^2\mu_X^4 .
\]
Here \(K_q\) is the arithmetic mean square of \(S_q\).
At \(q=12\) the energy scale is \(XH^2/(\log X)^3\); a bound with an
unspecified leading constant cannot give the variance. The local Fourier
route can use the uniform centered estimate (P12), of size
\(o(H^3\mu_X^4)\) per interval of length comparable to \(H\), with enough
room for the final interval averaging.
[Exact reduction Eq](PRESCRIBED_PRIME_BINARY_VARIANCE.md#the-exact-power-window-variance-reduction),
[local estimate and averaging](PRESCRIBED_PRIME_POWER_WINDOW_VARIANCE.md#comparing-the-local-binary-convolutions).

Winarski's Proposition 4.1 fixes angular type \(-n\), a base finite
conductor \(q_0\), a smooth profile \(W\) compactly supported on
\((0,\infty)\), and \(\delta,\varepsilon>0\). For \(X_0,K_0\ge1\), a
nonzero mask \(m\), a squarefree label \(f_0\), and row ideals
\(K_0/2\le Nj\le K_0\), it controls the sum of squared completed rows by
\[
 \ll_{n,q_0,W,\delta,\varepsilon}
 (N\mathfrak q)^\varepsilon(K_0N\mathfrak q)^{10\delta+\varepsilon}
 \left(K_0+\frac{K_0^2N\mathfrak q}{X_0}\right),\qquad
 \mathfrak q=(me_{q_0})f_0 .
\]
Here \(e_{q_0}\) records the fixed excluded primes. The row coefficient is
the actual sextic twist \(c_\chi(a)R_{(a)}(m_{q_0}^6 f^4 z)\).
Completion uses normalized cubic Gauss coefficients, a cube variable and
a theta reflection. Partitioning by powerful and masked parts yields
injective squarefree fibers away from the bad primes. The quadratic sieve
uses a common column vector after separating the reflected kernel. The
finite ray orbit and the state margins control the moving mask and the
two scales. In power notation \(K_0=X_0^\kappa\),
\(N\mathfrak q=X_0^\rho\), the two terms have exponents \(\kappa\) and
\(2\kappa+\rho-1\) before the displayed losses. The second term is at
most the first exactly when \(K_0N\mathfrak q\le X_0\). Proposition 6.3
instead bounds a fixed-character, fixed-profile smoothed Möbius sum with
excluded primes and totient damping. Its exponent \(15/16\) comes from
the sixth-power sampling budget \((11+5\sigma)/12\) as
\(\sigma\downarrow 1/20\).
[Winarski, §§3–6, especially (4.2), (4.4) and (6.5)](https://hexagonmath.org/pdf/2610.00458v1).

Two tests rule out direct substitutions.

1. The function \(\mathfrak a\mapsto f(N\mathfrak a)\) cannot be a finite
   linear combination of fixed admitted angular characters, even off a
   finite set of masked primes. Choose a rational prime
   \(p\equiv-1\bmod\operatorname{lcm}(q,3)\) outside the masks and
   conductors. For \(\mathfrak a=(p^r)\), all angular factors are one.
   If \(r\) is a common multiple of the orders of the finite residue
   characters, all the character values are one, exactly as on the unit
   ideal. The semigroup values are respectively zero and one since
   \(N(p^r)=p^{2r}\).
2. Allowing arbitrary bounded row-dependent columns invalidates the
   quadratic sieve bound. Choose disjoint families of prime ideals with
   norms at most \(Y\), each of size \(\gg Y/\log Y\), over distinct
   rational primes \(1\bmod3\), avoiding \(2,3\). If the row matrix is
   \(Q_j(g(I))\), set \(b_{j,I}=\overline{Q_j(g(I))}\). Coprimality makes
   every summand in each row one. The squared energy is
   \(\gg Y^3/(\log Y)^3\), whereas the common-vector bound with
   \(O(Y)\) squared column mass is \(O_\varepsilon(Y^{2+2\varepsilon})\).
   This fails for \(\varepsilon<1/2\).

A valid completed-row transfer would require an exact reflected-fiber
identity for the centered additive problem, or a decomposition into
common coefficient vectors with norm and multiplicity control. The
decomposition must retain the singular-series contribution, treat prime
powers, and pay the moving mask and smoothing costs at the target \(H\).
The Fourier separation of a smooth coupled kernel is available once the
arithmetic identity is supplied.

## A classical improvement and the type-zero specialization

The original variance proof used \(T=t^\tau\), a phase partition at
\(t^\gamma\), and the conditions \(\tau<1/5\),
\(1-\tau<\gamma<2\theta-1\). In the revised proof, write the twisted
Euler product as \(G_\chi(s)\prod_{\psi\bmod L}L(s,\chi\psi)^{1/\varphi(L)}\).
Log-loss Dirichlet convexity costs the exponent \((1-v)/2\) in height on
each slit. For \(v\le4/5\), Montgomery's zero-density exponent
\(3(1-v)/(2-v)\) is at most \((5/2)(1-v)\). For \(v\ge4/5\) Jutila's
\((2+\varepsilon)(1-v)\) with \(\varepsilon=1/4\) has the same upper
bound. The moduli are polylogarithmic, so their factors cost only log
powers. At abscissa \(v\), the combined height and Perron cost relative to
interval length is \(t^{-(1-3\tau)(1-v)}\), apart from log powers.
The left edge at \(\sigma_0\le3/5\) costs at most
\(t^{-2/5+(6/5)\tau}\), and truncation costs \(t^{1-\tau}\) before logs.
Thus \(1-\tau<\gamma<2\theta-1\) with \(\tau<1/3\) is possible for
\(\theta>5/6\). In the minor arcs, the revised \(X^{1/6}\) divisor
cutoffs and \(X^{5/6}\) largest-prime split make \(X^{5/6-\theta}\) the
first Type II error term. The Type I index is at most \(X^{1/2}\).
These changes extend (P1) and (P14) with their original log savings.
[Full proof](PRESCRIBED_PRIME_POWER_WINDOW_VARIANCE.md),
[Jutila, Theorem 1, (1.8)](https://journals.msp.org/mscand/article/download/1943/1942/1974).

There is a genuine finite-order interface with Winarski's \(7/8\)
theorem. For a Dirichlet character \(\xi\bmod d\), \(\xi\circ N\) is
admitted with type zero and modulus \(d\mathcal O_K\). The split, inert
and ramified Euler factors give
\[
 L_K(s,\xi\circ N)=L(s,\xi)L(s,\xi\chi_{-3}).
\]
Here \(\chi_{-3}\) is the quadratic character of \(K\), with products
interpreted as imprimitive at the displayed moduli. There is no pole
of either factor in \(7/8<\Re s<1\), and the principal pullback treats
\(\zeta\). Thus the attributed result would give a common upper bound
\(v\le7/8\) for the slits. With only Montgomery's density bound and
log-loss convexity the resulting budget is
\(\tau(8/3+1/2)<1\), which gives \(\theta>16/19\) after the minor-arc
retuning. The classical \(\theta>5/6\) range is stronger. In the revised
classical proof, a fixed strip would instead give a power saving in the
short-interval character error, since \(1-v\ge1/8\). That gain does not
lower the \(\tau<1/3\) edge budget or the minor-arc bound.
[Winarski, Definition 1.1 and Theorem 1.2](https://hexagonmath.org/pdf/2610.00716v1),
[Bhattacharjee–Bhattacharya, p. 3 for the factorization](https://hexagonmath.org/pdf/2610.00256v3).

## Uniform prime-ideal windows conditional on the 7/8 theorem

The \(7/8\) proof closes each finite angular cap with a common boundary.
The probe estimates may have character-dependent thresholds and retain
a presentation-modulus-squared factor in the detector. The quantitative
prime consequence below combines the zero-free conclusion with the
ordinary degree-one Hecke functional equation. Its constants come from
this separate argument.
[Winarski, Theorem 1.2, (2.3), (5.1), Proposition 5.1](https://hexagonmath.org/pdf/2610.00716v1).

Assume the stated half-plane with \(\beta=7/8\). Let \(\chi\) have angular
order \(k\) and presentation modulus \(\mathfrak m\), and let
\(\mathbf1_{\mathrm{pr}}\) indicate the principal character. For any
fixed \(0<\varepsilon<1-\beta\), the uniform bound on
\(\beta+\varepsilon\le\sigma\le2\) is
\[
 \left|\frac{L'}{L}(\sigma+it,\chi)
             +\frac{\mathbf1_{\mathrm{pr}}}{\sigma+it-1}\right|
 \ll_\varepsilon \log\!\big(2N\mathfrak m(3+|k|+|t|)\big).
\]
In the principal case the expression is continued across \(1\). To
prove it, the primitive functional equation with complex gamma factor
\(\Gamma(s+|k|/2)\), Stirling on a fixed left line, and the bounded Euler
product at \(2\) give polynomial growth in \(N\mathfrak m(3+|k|+|t|)^2\)
on fixed strips by Phragmén–Lindelöf. The finitely removed Euler factors
also have polynomial growth. Regularize the principal factor by
multiplying by \((s-1)/(s+1)\). On the disk centered at \(2+it\) with
radius \(2-\beta-\varepsilon/2\), the regularized function has a
holomorphic logarithm. Its real part has the stated logarithmic upper
bound, and the Euler product bounds the real part at the center on
both sides. Borel–Carathéodory on an intermediate disk and Cauchy's
derivative estimate on radius \(2-\beta-\varepsilon\) give the bound.
The regularization adds only a bounded \(1/(s+1)\) term. A real norm
twist replaces \(t\) by \(t-\tau\) and shifts the principal pole to
\(1+i\tau\).

Write \(\Lambda_K(\mathfrak p^j)=\log N\mathfrak p\), with zero on other
nonzero ideals. For fixed \(w\in C_c^\infty(0,1)\), \(X\ge2\) and
\(0<H\le X\), put \(\Delta=H/X\). The profile
\(w((u-1)/\Delta)\) has Mellin transform bounded by
\(O_{w,A}(\Delta(1+|t|\Delta)^{-A})\) on the fixed shift strip, by
changing to the logarithmic variable and integrating by parts.
Shifting \(-L'/L\) to \(\beta+\varepsilon\) gives
\[
 \sum_{\mathfrak a}\Lambda_K(\mathfrak a)\chi(\mathfrak a)
     w\!\left(\frac{N\mathfrak a-X}{H}\right)
 =\mathbf1_{\mathrm{pr}}H\int w
  +O_{w,\varepsilon}\!\left(
    X^{\beta+\varepsilon}
    \log\!\big(2N\mathfrak m(3+|k|+X/H)\big)\right).
\]
The pole at \(1\) supplies the main term. The error follows by integrating
the Mellin bound against the displayed logarithmic derivative bound.
The total contribution from higher prime powers is
\(O(\sqrt X\log^2 X)\), by counting ideals of norm at most \(\sqrt X\)
and summing over possible exponents. The same estimate therefore
holds for log-weighted prime ideals alone.
For polynomially bounded order and finite modulus this gives a
power-saving relative error for \(H\ge X^{\beta+\delta}\), with
\(0<\varepsilon<\delta\). In a fixed ray modulus, a smooth angular bump
of aperture \(\eta\) on the unit-quotient circle has Fourier coefficients
\(O_A(\eta(1+\eta|k|)^{-A})\). Its principal mean is of size \(H\eta\).
Summing the admissible modes costs a bounded absolute Fourier mass and a
logarithmic height factor with \(\eta^{-1}\) inside the logarithm. The
sufficient power budget is \(H\eta\ge X^{\beta+\delta}\) when
\(\eta^{-1}\) is polynomially bounded.

EisPower and EisLocationWindows need quaternionic matrix coefficients
at a fixed shell index. At finite cofactor levels the estimate must
handle the two component norms and the fixed-sum relation. A
one-variable Hecke prime-ideal estimate supplies no bound for those
matrix coefficients or their joint degree sum. The moving bump uses
absolute Legendre mass \(O(\sqrt{M/H})\) and the shell error has size
\(n^{1/2+o(1)}\sqrt{M/H}\) against a main of size \(H\). This gives
\(3\theta>1+\gamma\) at \(H=n^\theta\), \(M\asymp n^\gamma\), before
choosing the finite-level savings. At the center this is \(3\theta>2\).
The positive convolution example in the repository gives an empty
central band even with the existing single-index operator envelope and
scalar coset means. Moving below this threshold requires information
such as joint degree cancellation or relations between actual Hecke
indices.
[Moving-window proof](EISENSTEIN_MOVING_WINDOWS.md),
[single-index operator barrier](EISENSTEIN_SINGLE_INDEX_OPERATOR_BARRIER.md).

## Critical switching and Pólya–Vinogradov

The required uniform absolute estimate is
\(\sum_{d<D,\ d\mid P_7(z)}|r_d(N)|\ll X_N(\log N)^{-A}\) at
\(D=N^{1/2}(\log N)^{-B}\), \(z=N^{1/2-\epsilon}\), where
\(r_d(N)=\sum_{d\mid b}U(N-b)V(b)-X_Ng_N(d)\),
\(U=\tau^2\mathbf1_{C_{12}}\), and \(V=\tau\mathbf1_{A_3}\).
The signed alternative keeps the actual lower-sieve weights and requires
a lower remainder bound \(-o(X_N\mathcal V_N(z))\), uniformly in \(N\).
The positive lower main has size \(\sqrt\epsilon\,X_N\mathcal V_N(z)\);
the bad survivors cost \(O(\epsilon^{3/2}X_N\mathcal V_N(z))\).
The large-prime part of the lower weight and the complementary branch
both enter the remainder.
[Critical criterion and signed form](ASYMMETRIC_CRITICAL_SWITCHING.md#5-the-critical-level-implication),
[exact switched branches](ASYMMETRIC_WELL_FACTORABLE_SWITCH.md).

The transfer theorem gives a prefix bound \(M(\chi)\ll_\beta\sqrt Q
\log\log Q\) for a primitive character modulo \(Q\ge3\) under a fixed
zero-free rectangle for every twist by a character of modulus at most
\(\lceil\log Q\rceil^{16}\), to height \(5Q^2\). It applies via partial
summation to an isolated character factor with a bounded-variation
weight. The numerical bound saves a logarithmic factor over classical
Pólya–Vinogradov. Corollary 3.8 also gives a one-variable prime sum
\(O_\beta(x^{1-\eta_\beta}\log(Q+2))\) for a nonprincipal character under
its zero-free rectangle to height \(4x^2+2\), with
\(\eta_\beta=(1-\beta-6/2000)/4\) at \(\beta=11/12\) or \(7/8\).
Even at \(7/8\) this yields \(\eta_\beta=61/2000\). Absolute summation
over nonprincipal characters in progression orthogonality keeps this
same error scale, with no \(1/\varphi(Q)\) saving. At \(x\asymp N\),
\(Q=N^{3/4+o(1)}\), the available exponent \(1939/2000\) exceeds the
expected progression length exponent \(1/4\). Theorem 2.3 and the
exceptional-character result concern a character average. Control of
the target-dependent switched weights against those exceptional
characters would require an additional estimate.
[Bhattacharjee–Bhattacharya, Theorems 1.1, 1.2, 2.3 and Corollary 3.8](https://hexagonmath.org/pdf/2610.00256v3).

After the exact compression \(c=k^2\ell\), even the \(k=1\) branch can
take \(c\) near \(D\) and the smaller companion divisor \(e\) near
\((N/c)^{1/2}\). The progression modulus is then
\(Q=ce\asymp (Nc)^{1/2}=N^{3/4+o(1)}\). Reduced classes can have twists
of conductor of order \(Q\), and nonunit classes require gcd strata.
The Pólya–Vinogradov bound at that scale is \(N^{3/8+o(1)}\), so it is
weaker than trivial counting on an isolated free factor of length
\(N^{1/4}\) in a balanced four-factor box. This leaves the divisor
weights, complementary factor, and sum over coefficients untreated.
The exact sign-transfer example makes the issue concrete. The weight
\(W=\chi_4 V\) has prefixes \(O_\nu(x^{1/2+\nu})\), but its binary
correlation at \(c=1\) is \(X_N+O(N^{1-\eta})\), because all actual
pairs force \(\chi_4(t)=1\). At general \(c\) the sign changes the
local main by \(\chi_4(c)\). The residual target is dispersion for
the centered, signed correlation with that main and all target
valuations present.
[Combined moduli, (27)–(28)](ASYMMETRIC_WELL_FACTORABLE_SWITCH.md#7-why-the-critical-application-still-lies-beyond-these-inputs),
[exact resonance, (1)–(8)](ASYMMETRIC_SIGN_TRANSFER.md#two-nonprincipal-factors-with-an-additive-resonance).
