# Binary sums with prime factors congruent to three modulo four

This is the \(3\bmod4\) specialization of
[Theorem PrimeCosetSelfSums](PRIME_COSET_SELF_SUMS.md#the-pointwise-result).
The arithmetic formulas below make its negative character branch and
leading constant explicit; the analytic comparison is proved once in
the general theorem.

It uses exactly inputs ZF and G of the
[existing input audit](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md): the
fixed all-conductor zero-free half-plane and physical weighted-dilation
root-moment theorem. This is an internally reviewed application of
those inputs, without external human review or formal verification.

## 1. Pointwise count and local criterion

Put
\[
 C^-_4=\{n\ge1:p\mid n\Longrightarrow p\equiv3\pmod4\},\qquad
 R_-(N)=\#\{(a,b)\in(C^-_4)^2:a+b=N\}.
\]
We include 1 and count ordered numerical pairs. Define
\[
 A_-=\frac12\prod_{p\equiv1(4)}
                  \left(1-\frac2{p(p-1)}\right)>0,\qquad
 J_-(N)=\prod_{\substack{p\mid N\\p\equiv1(4)}}\frac{p-1}{p-2}.
                                                               \tag{M.1}
\]

**Theorem PrimeThreeModuloFourCorrelation.** Under the same ZF and G
inputs as the subgroup correlation theorem, uniformly over sufficiently
large even integers \(N\),
\[
 R_-(N)=A_-J_-(N)\frac N{\log N}
       \left(1+O((\log N)^{-1/16384})\right).                    \tag{M.2}
\]
Odd targets have zero count identically. Evenness is also exactly the
condition for representability in every finite quotient. Thus every
sufficiently large locally admissible target has a representation.

The error is uniform in all target prime valuations. There is no
numerical cutoff or exact small-exception classification. No assertion
for prescribed residue states of the individual summands, or for all
non-subgroup prime prescriptions, is part of this theorem.

Here the generated group is
\(H=\langle3\bmod4\rangle=(\mathbb Z/4)^\times\), but the prime density
is \(\rho=1/2\). In particular \(\rho\ne1/[(\mathbb Z/4)^\times:H]\).
The subgroup density-index relation must not be used.

## 2. The marginal and the negative secondary branch

Write \(f=\mathbf1_{C^-_4}\), \(\chi=\chi_{-4}\), and
\[
 F(s)=\prod_{p\equiv3(4)}(1-p^{-s})^{-1},\qquad
 P_3(s)=F(2s),\quad P_3=P_3(1).
\]
Euler factors give
\[
 F(s)^2=(1-2^{-s})\frac{\zeta(s)}{L(s,\chi)}P_3(s),\qquad
 F_\chi(s):=\sum_n\chi(n)f(n)n^{-s}=\frac{F(2s)}{F(s)}.
                                                               \tag{M.3}
\]
Near one, with positive determinations on the real line to its right,
\[
 F(s)=(s-1)^{-1/2}h_0(s),\qquad
 F_\chi(s)=(s-1)^{1/2}h_\chi(s),
\]
where
\[
 h_0(1)^2=\frac{2P_3}{\pi},\qquad
 \frac{h_\chi(1)}{h_0(1)}=\frac\pi2,\qquad
 c_-=\frac{h_0(1)}{\Gamma(1/2)}=\frac{\sqrt{2P_3}}\pi.
\]
The fixed-function Selberg--Delange theorem therefore gives
\[
 \sum_{n\le X}f(n)\sim c_-\frac X{\sqrt{\log X}},\qquad
 \sum_{n\le X}\chi(n)f(n)
 =-\frac{\pi c_-}{4}\frac X{(\log X)^{3/2}}
       +O\left(\frac X{(\log X)^{5/2}}\right).                  \tag{M.4}
\]
The corrected author edition of de la Bretèche--Tenenbaum, Theorem 1.2,
p. 3, permits \(-1/2\); (1.10) gives
\(\Gamma(-1/2)^{-1}=-1/(2\sqrt\pi)\). Both odd states have leading
density \(c_-/2\), but arbitrary logarithmic centering requires the
secondary branch as well.

## 3. Two-branch rough density

Set
\[
 L=\log X,\quad y=\exp(L^{2/5}),\quad S=\{p:p\le y\},\quad
 V=\prod_{p\le y}(1-1/p).
\]
Write \(n=A v\), with \(A=n_S\) and \((v,\prod_{p\le y}p)=1\).
Removing the allowed small Euler factors gives
\[
 F_R(s)=(s-1)^{-1/2}h_{0,R}(s),\qquad
 F_{\chi,R}(s)=(s-1)^{1/2}h_{\chi,R}(s).
\]
Exactly,
\[
 F_{\chi,R}(s)=\frac{F_R(2s)}{F_R(s)},\qquad
 \frac{h_{\chi,R}(1)}{h_{0,R}(1)}
   =\frac\pi2\prod_{\substack{p\le y\\p\equiv3(4)}}\frac{1+1/p}{1-1/p}
   \sim e^\gamma\log y.                                      \tag{M.5}
\]
In particular \(h_{0,R}(1)\asymp(\log y)^{-1/2}\) and
\(h_{\chi,R}(1)\asymp(\log y)^{1/2}\).
A disk of radius \(c/\log y\) gives, for each fixed \(j\),
\[
 |h_{\sigma,R}^{(j)}(1)|/j!
       \ll_j h_{\sigma,R}(1)(C\log y)^j,\qquad \sigma=0,\chi.
                                                               \tag{M.6}
\]
This is the Cauchy estimate in the general coset model.

For a fixed sufficiently deep \(K\), put
\[
 \begin{split}
 M_0(t)&=\sum_{j=0}^K\frac{h_{0,R}^{(j)}(1)}{j!}
             \frac{(\log t)^{-1/2-j}}{\Gamma(1/2-j)},\\
 M_\chi(t)&=\sum_{j=0}^K\frac{h_{\chi,R}^{(j)}(1)}{j!}
             \frac{(\log t)^{-3/2-j}}{\Gamma(-1/2-j)}.
 \end{split}                                                  \tag{M.7}
\]
Uniformly on every fixed power bulk \(X^a\le t\le X^C\),
\[
 M_0(t)>0,\qquad |M_\chi(t)|\ll L^{-3/5}M_0(t).                 \tag{M.8}
\]
Use the general model's common smooth extension, setting \(M_\chi=0\)
below a smaller power bulk and extending \(M_0\) positively with
\(M_0\ll h_{0,R}(1)L^{-1/2}\), while preserving (M.8). Define
\[
 D(n)=\frac{f(A)}V\{M_0(v)+\chi(v)M_\chi(v)\},\qquad
 D_0(n)=\frac{f(A)}V M_0(v),\qquad \xi(n)=f(n)-D(n).             \tag{M.9}
\]

**Lemma MinusFourInvariantModel.** The model is nonnegative,
\(D(n)\ll L^{-3/10}\), and for every allowed \(p\le y\),
\[
 f(pn)=f(n),\quad D(pn)=D(n),\quad \xi(pn)=\xi(n)               \tag{M.10}
\]
exactly, including repeated factors. It has arbitrary logarithmic
rough centering against bounded periodic tests of period at most
a fixed power of \(L\), with fixed polynomial logarithmic smooth
costs, after \(K\) is chosen sufficiently large.

This is the specialization of the model and Lemma
`lem:coset-centering` in
[cosets.tex](../../unpublished/pointwise-binary-splitting-semigroups/cosets.tex),
with \(m=4,K=\{1\},a=3\). Invariance follows because multiplying by
an allowed small prime changes \(A\) and leaves \(v\) unchanged.
The character is evaluated at \(v\). The size bound follows from
\(h_{0,R}(1)\asymp L^{-1/5}\), \(V\asymp L^{-2/5}\).
For an inflated character \(\psi\), the specialized factorization is
\[
 F_\psi(s)^2=\frac{L(s,\psi)}{L(s,\psi\chi)}
       \prod_{p\equiv3(4)}
                (1-\psi(p)^2p^{-2s})^{-1}.                    \tag{M.11}
\]
Only \(1,\chi\) have singular branches, with exponents \(1/2,-1/2\).
Lemma `lem:coset-centering` proves the moving-\(y\) estimate for every
polylogarithmic period and retained record \(A\le e^{\sqrt L}\),
using the full two-branch model. Its analytic-log bound controls the
reciprocal factor in (M.11); its Hankel remainder is
\(O_K(L^{-3(K+1)/5})\) relative to the principal density.
The model depth is chosen after the requested accuracy. A scalar
model would leave the negative branch uncanceled.

## 4. Checks for the shared comparison theorem

The general coset proof checks all four hypotheses of
[the shared comparison theorem](SPARSE_ARITHMETIC_MODEL_COMPARISON.md).
Here its data are \(m=4,c=3,A=B=1\),
\(\rho=1/2,H=(\mathbb Z/4)^\times,d=2,g_S=1\), and
\(\mathcal S_N=NJ_-(N)/\log N\).
The specialized independent checks are:

* (C1)--(C2): (M.8)--(M.10) give nonnegative bounded models and exact
  invariance for all allowed small primes.
* (C3): Lemma MinusFourInvariantModel retains both character means,
  with the depth chosen after the precision.
* (C4): the forbidden primes are exactly \(1\bmod4\). Section
  “Common banks and the independent sieve check” of `cosets.tex`
  gives denominator orders \(\log w\), \((\log y\log w)^{1/2}\),
  and \(\log y\) for the three pair types. The factors
  \(D\ll L^{-3/10}\) restore density \(L^{-1}\). Its exact bank
  law, interval factor and CRT error are unchanged here.

The bound \(1\le J_-(N)\ll(\log\log(3N))^2\) puts this sparse
scale in the theorem's permitted range. None of these checks assumes
a binary cancellation or residue-energy estimate.

**Lemma MinusFourCosetFrame.** The common progression frame applies to
\(3\pmod4\) with unchanged kernel and reproduction bounds: equivalent
characters differ on \(t\) labels by \((-1)^t\), and both frame
products have the same slot type. This is the shared slot-phase check.

The comparison conclusion now gives, with \(\delta=1/16384\),
\[
 \sum_{a+b=N}\bigl(f(a)f(b)-D(a)D(b)\bigr)
       \ll \frac{NJ_-(N)}{L^{1+\delta}},
\]
and the same absolute discrepancy for interval restrictions.
The target-grid energy and reconstruction are conclusions of the shared
proof. Its saving is \(1/(8192\varphi(4))=1/16384\); the secondary
gap costs only \(3g_S/10=3/10\).

## 5. Leading model and arithmetic factors


On the retained bulk and records,
\[
 |D-D_0|\ll L^{-3/5}D_0,\qquad
 |D*D-D_0*D_0|\ll L^{-3/5}(D_0*D_0).                          \tag{M.16}
\]
The general main-term calculation gives finite factor two on even
targets: odd pairs occupy half the line and each odd marginal has
normalization \(1/2\). At \(p\equiv1(4)\), the normalized factor is
\[
 \frac{1-(2-\mathbf1_{p\mid N})/p}{(1-1/p)^2};
\]
allowed good primes have factor one. Thus the model coefficient is
\(2c_-^2 C_1J_-(N)\), where
\(C_1=\prod_{p\equiv1(4)}(1-(p-1)^{-2})\).
The identity
\(\prod_{p\ {\rm odd}}(1-p^{-2})^{-1}=\pi^2/8\) gives
\[
 2c_-^2C_1
 =\frac{4P_3 C_1}{\pi^2}
 =\frac12\prod_{p\equiv1(4)}\left(1-\frac2{p(p-1)}\right)=A_-.
\]
For this specialization \(\Gamma_S(N)=2\) on even targets and zero
on odd targets. The finite-image lemma `lem:coset-finite-images`
and Theorem `thm:coset` therefore give exactly the stated local
criterion and (M.2). Their main-term proof evaluates the leading
model by the fundamental sieve lemma; (M.16) and the positive-density
variation \(O(L^{-1/2})\) are smaller than the comparison remainder.

## 6. Historical scope and limits

Dietmann--Elsholtz, *Sums of two numbers having only prime factors
congruent to one modulo four*, Banach Center Publications 118 (2019),
177--179, proves the plus-class existence theorem. Its final remark,
printed p. 179, explicitly poses the prime-\(3\bmod4\) sum problem
on \(N\equiv2\bmod4\). The
[published primary PDF](sources/dietmann-elsholtz-prime-factors-2019-published.pdf)
is retained. The problem is an established historical question.
That 2019 statement does not determine its current literature frontier.

The Gaussian/Greaves argument cannot simply exchange prime classes:
\(5=1^2+2^2\) is a primitive norm divisible by a \(1\bmod4\) prime.
This application and the [general coset theorem](PRIME_COSET_SELF_SUMS.md)
reuse the source graph and zero-free theorems. No historical priority,
external review or exact onset is asserted. Arbitrary residue sets
and summand-state filters still require independent checks.
