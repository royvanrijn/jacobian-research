# Exact sign transfer and the remaining nonperiodic coefficient

Let \(C=C_{12}\), \(A=A_3\), \(U=\tau^2\mathbf1_C\), and
\(V=\tau\mathbf1_A\), as in
[AsymPointwise](ASYMMETRIC_GAUSSIAN_EISENSTEIN_SUMS.md), and let
\(N\equiv2\pmod{12}\). Write \(\chi_4=\chi_{-4}\), lifted
literally to modulus 12 when an \(L\)-function is used. Thus
\(\chi_4(p)=-1\) at primes \(p\equiv7\pmod{12}\).
Every integer with positive \(V\)-weight is a unit modulo 6.
Write \(C_c(N)\), \(\mathcal M_c(N)\), and \(X_N\) as in
[AsymAuxiliaryAverage](ASYMMETRIC_AUXILIARY_AVERAGE.md), and let
\(\chi_3=\chi_{-3}^{(6)}\), \(\chi_{12}=\chi_3\chi_4\) on units.
The function \(g_{E,6}=(1*\chi_{-3})\mathbf1_{(\cdot,6)=1}\)
is the two-character core from AsymPointwise.

## Proposition AsymCharacterResonance

Fix \(1<z<D\) and \(K\ge1\), and use the actual lower
half-dimensional weights
\[
 \lambda_d^-=
 \mu(d)\mathbf1_{d\mid P_7(z)}\mathbf1_{d\in\mathcal D_D^-}.
\]
In the finite-kernel compression let \(c=k^2\ell\), where
\(k,\ell\) are squarefree, coprime, and supported on primes 7 modulo
12, with \(k\le K\). Then \(d=k\ell=\operatorname{rad}(c)\), and
\[
 \mu(k)2^{\omega(\ell)}\lambda_{k\ell}^-
 =2^{\omega(\ell)}\mu(\ell)
      \mathbf1_{\operatorname{rad}(c)\mid P_7(z)}\mathbf1_{\operatorname{rad}(c)\in\mathcal D_D^-}
 =\chi_4(c)\Psi(c),\tag{1}
\]
where
\[
 \Psi(c)=2^{\omega(\ell)}
   \mathbf1_{k\le K}\mathbf1_{\operatorname{rad}(c)\mid P_7(z)}
   \mathbf1_{\operatorname{rad}(c)\in\mathcal D_D^-}
\]
is nonnegative and vanishes unless \(c\) is cube-free and supported on primes 7 modulo 12.
The common-square cutoff \(k\le K\) is included in \(\Psi\).
Any further allowed coefficient restriction may be included there too.

On every pair \(a+ct=N\) with \(U(a)V(t)>0\),
\(a\equiv ct\equiv1\pmod4\). Consequently
\[
 \chi_4(c)=\chi_4(t),\qquad
 C_c^W(N):=\sum_{a+ct=N}U(a)\chi_4(t)V(t)
             =\chi_4(c)C_c(N).\tag{2}
\]
Therefore the ACTUAL compressed linear form equals
\[
 \sum_c\chi_4(c)\Psi(c)C_c(N)
       =\sum_{a+ct=N}U(a)\Psi(c)W(t),\qquad W=\chi_4V.\tag{3}
\]
Its remainder, however, is
\[
 \sum_{a+ct=N}U(a)\Psi(c)W(t)
       -\sum_c\chi_4(c)\Psi(c)\mathcal M_c(N).\tag{4}
\]
The second, signed main-model term is retained. It is not zero by any
orthogonality used in(2).

## Two nonprincipal factors, with an additive resonance

The Dirichlet series of \(W\) is
\[
 \sum_{n\ge1}\frac{W(n)}{n^s}
 =E_V(s)L(s,\chi_4^{(12)})L(s,\chi_{12}),\qquad
 E_V(s)=\prod_{p\ge5,p\equiv2\ (3)}(1-p^{-2s}).\tag{5}
\]
Both displayed characters are nonprincipal. The correction is unchanged
by twisting because its indices are squares. A two-variable hyperbola
and bounded partial sums of these fixed characters, followed by the
absolutely convergent square kernel, give
\[
 \sum_{n\le X}W(n)\ll_\nu X^{1/2+\nu}\qquad(\nu>0).\tag{6}
\]
For the kernel step,
\(\sum_{v\le X}|e(v)|v^{-1/2}
 \le X^\nu\sum_v|e(v)|v^{-1/2-\nu}\), and the latter norm is finite.

Nevertheless (2) has a nonzero pointwise main. In particular
\[
 \sum_{a+t=N}U(a)W(t)=C_1(N)=X_N+O(N^{1-\eta}),\tag{7}
\]
by AsymPointwise, for all sufficiently large admissible targets.
This is an exact counterexample to inferring additive cancellation
from the two nonprincipal factors or the prefix bound(6).
For every fixed coefficient \(c\) supported on primes 7 modulo 12,
the main of \(C_c^W\) is \(\chi_4(c)\mathcal M_c\), as also follows
from [AsymAuxiliaryAverage](ASYMMETRIC_AUXILIARY_AVERAGE.md). Thus the sign
transfer preserves, rather than removes, the target-dependent main.

The resonance also appears in the divisor switch. In the base gcd
stratum, \(g\in C_{12}\). After the normalized \(G_g\) expansion, its
kernel index \(v\) is 1 modulo 4. This support assertion uses that
restriction on \(g\); it is not asserted for arbitrary shifts in \(A_3\). For
\(u\equiv1\pmod6\),
\[
 \chi_4(u)g_{E,6}(u)
 =2\sum_{q\mid u,\ q<\sqrt u}\chi_{12}(q)\chi_4(u/q)
   +\mathbf1_{u=\square}\chi_{12}(\sqrt u)\chi_4(\sqrt u).
\]
On every actual \(a+cvqr=N\) with \(U(a)>0\),
\[
 \chi_{12}(q)\chi_4(r)=\chi_4(c)\chi_3(q).\tag{8}
\]
Indeed \(\chi_4(cvqr)=1\), \(\chi_4(v)=1\), and
\(\chi_{12}=\chi_4\chi_3\) on units.
Thus freezing the complementary variable's congruence class restores
exactly the original \(\chi_3(q)\) phase and the \(\chi_4(c)\) outer sign. It does
not turn the large-c family into an unweighted inner modulus sum.
This is a failure of that proposed simplification, not a proof that
no further dispersion or switching argument can exploit(3).

## The nonnegative coefficient has an exact convolution identity

Define
\[
 M_7(n)=\mu(n)^2\mathbf1_{p\mid n\Rightarrow p\equiv7\ (12)},
 \qquad b=M_7*M_7,
\]
where * is Dirichlet convolution. Then
\[
 b(c)=\begin{cases}
 2^{\omega(\ell)},&c=k^2\ell,\ k,\ell\text{ squarefree, coprime and p7-supported},\\
 0,&\text{otherwise}.
 \end{cases}\tag{9}
\]
At a prime 7 modulo 12 its local series is \((1+z)^2\); all other
local factors are 1. In (1), \(\Psi\) is precisely \(b(c)\) with
the indicated Rosser, prime, and common-square cutoff restrictions.
Without these restrictions there is an exact signed convolution identity
with nonnegative output:
\[
 b*W=\tau\mathbf1_C.\tag{10}
\]
At primes 7 modulo 12 the two factors cancel as
\((1+z)^2(1+z)^{-2}=1\). At primes 1 modulo 12 only \(W\)
remains, with factor \((1-z)^{-2}\). At every other prime both
sides have local factor 1. This proves (10) coefficient by coefficient.
Thus omitting the selector and extending the entire convolution does
not create a simpler zero-main sequence: it reconstructs the desired
modulus-12 second summand.

## Proposition AsymRosserAutomaticRange

For squarefree \(d=p_1\cdots p_r\), with \(p_1>\cdots>p_r\),
the standard full lower half-dimensional support is
\[
 \mathcal D_D^-=\{d:
    p_1\cdots p_{2j-1}p_{2j}^2<D\text{ for every }2j\le r\}.
 \tag{11}
\]
For every even index \(2j\),
\[
 p_1\cdots p_{2j-1}p_{2j}^2
 =\left(\prod_{i\le2j}p_i\right)p_{2j}
 \le d^{1+1/(2j)}\le d^{3/2}.
\]
Consequently
\[
 d<D^{2/3},\quad d\mid P_7(z)
       \quad\Longrightarrow\quad\lambda_d^-=\mu(d)=\chi_4(d).
 \tag{12}
\]
This checks ALL support conditions, not just the first pair.
The exponent \(2/3\) is sharp for a uniform automatic-support statement.
Choose distinct primes \(p>q\), both 7 modulo 12 and of comparable
large size, and put \(D=pq^2/2\). Then \(d=pq\notin\mathcal D_D^-\),
while \(\log d/\log D\to2/3\). For any fixed \(\gamma>2/3\),
eventually \(d<D^\gamma\). These primes are also below the critical
\(z\), whose order exceeds \(D^{1/3}\).

The actual selector is not determined by a fixed residue class.
For example, take \(D=1000\) and \(31<z\le1000\). Then
\[
 d_1=19\cdot7=133,\quad19\cdot7^2=931<D,
 \qquad d_2=31\cdot7=217,\quad31\cdot7^2=1519>D.
\]
Both integers are 1 modulo 12, but only \(d_1\) belongs to
\(\mathcal D_D^-\). The \(k=1\) slice of (9) therefore has selected
weights 4 and 0 at the same residue. This issue occurs beyond the sharp automatic zone.

The primary source for (11) is H. Iwaniec, *The half dimensional sieve*,
Acta Arithmetica **29** (1976), 69–95, printed p.75, equations
(3.5)–(3.6): [DOI](https://doi.org/10.4064/aa-29-1-69-95).
The definition is first written with a finite pairing cutoff \(R\);
choose \(R\) large enough to cover every possible prime-factor count
for \(d<D\). Its displayed conditions then give (11).
The primary spread on PDF page 4 was inspected in the same scan
identified in AsymSemilinear, SHA256
`ca6e84aafb32ee62dcd5ce76679c65652bf677ccd4e36d200576024fac1eca67`.

## Theorem AsymPeriodicReplacementObstruction

For either \(f=b\) or
\[
 f(n)=b_{sf}(n):=2^{\omega(n)}M_7(n),
\]
and every fixed integer \(m\ge1\), there is \(c_f>0\) such that
\[
 \inf_{h_X\ m\text{-periodic}}
       \sum_{X<n\le2X}|f(n)-h_X(n)|^2\sim c_fX.\tag{13}
\]
The coefficients of \(h_X\) may be complex, unbounded, and depend on \(X\).

Both choices are nonnegative and multiplicative, with \(f(p)=2\)
exactly at primes 7 modulo 12 and zero at all other primes. Their
prime-square and higher-power values are bounded, and \(f\le\tau\).
The one-variable Frobenian mean estimate therefore gives
\[
 \sum_{X<n\le2X}f(n)\ll X/\sqrt{\log X},\qquad
 \sum_{X<n\le2X}f(n)^2\sim c_fX,\quad c_f>0.\tag{14}
\]
For the second estimate, \(f(p)^2=4\) on a prime set of density
\(1/4\), so its prime mean is 1. The higher-prime-power series
converge absolutely beyond \(1/2\). Both functions meet the cited
lemma's divisor growth and Frobenian prime hypotheses. Only the order
of the first mean is used.

Write \(n_r=\#\{X<n\le2X:n\equiv r\pmod m\}\) and
\(S_r=\sum_{X<n\le2X,\ n\equiv r\ (m)}f(n)\). Orthogonal projection
onto the functions of residue class gives
\[
 \inf_{h_X}\sum_{X<n\le2X}|f(n)-h_X(n)|^2
       =\sum_{X<n\le2X}f(n)^2-\sum_{r\ (m)}S_r^2/n_r.
\]
For fixed \(m\), \(n_r\sim X/m\). Positivity yields
\[
 0\le\sum_r S_r^2/n_r\ll_m X^{-1}\left(\sum_r S_r\right)^2
                   \ll_m X/\log X=o(X).
\]
Equations (14) prove the sharp asymptotic (13). No AP mean theorem
for \(f\) is needed.

This obstruction applies even before the Rosser selector becomes
nontrivial: on the \(k=1\) slice, if
\(2X<\min(D^{2/3},z)\), that selector is identically 1, yet the
coefficient is \(b_{sf}\) and has (13). Thus neither the prime-factor
support nor the factor \(2^{\omega(\cdot)}\) can be replaced by a fixed-periodic
coefficient with a negligible error on its natural mean-square scale.
It prevents invoking FixedPeriodicBoxes through THAT replacement.
It does not rule out estimates that retain these nonperiodic weights,
or establish failure of any every-target representation conjecture.

The mean input in (14) is
[Loughran–Matthiesen, Lemma 2.8](https://arxiv.org/html/1904.12845v2#S2),
with the [cached primary PDF](sources/loughran-matthiesen-1904.12845v2.pdf).
No growing-modulus or short-interval mean is asserted. The remaining
critical signed remainder is open.

## Dependencies and finite controls

The exact sign and convolution identities are elementary. The nonzero
additive main in (7) uses AsymPointwise. The lower-weight support is
Iwaniec's displayed definition, also retained in
[AsymSemilinear](ASYMMETRIC_SEMILINEAR_SIEVE.md). The fixed-coefficient
main comparison uses AsymAuxiliaryAverage. The periodic projection
argument uses only the one-variable means from the cited Lemma 2.8.

The [standard-library checker](eisenstein-restricted-prime-sums/verify_asymmetric_sign_transfer.py)
and its [receipt](eisenstein-restricted-prime-sums/asymmetric_sign_transfer_checks.json)
retain exact finite convolution, binary sign, Rosser-support and
least-squares projection controls. The receipt binds this note and its
primary and canonical dependencies. These controls do not establish the
asymptotic means, a dispersion estimate, or the critical remainder.
The mathematical review recorded here is internal; external specialist
review, formal verification and priority are not established. No new
finite binary coverage scan is used.
