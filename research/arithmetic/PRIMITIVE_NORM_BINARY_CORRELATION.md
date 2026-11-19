# Pointwise primitive-norm correlations: a conditional dilation-graph proof

**Status.** The theorem below assumes a uniform fixed zero-free half-plane
for Dirichlet \(L\)-functions and imports the precise weighted-dilation
root-moment theorem identified below. It is a conditional proof of the
pointwise correlations, with those inputs exposed. It does not establish
the unconditional research objective or independently certify the recent
manuscripts supplying the inputs.

The proof handles the Gaussian case, every fixed imaginary quadratic
order via its local shadow, and arbitrary fixed abelian splitting
semigroups, including \(C_{12}+C_{12}\). Its additional ingredients are
proved in [the supporting arguments](PRIMITIVE_NORM_BINARY_CORRELATION_SUPPORT.md):
a projected rough density, a sparse bank envelope, fixed-progression
ideal kernels, a large-shift graph transfer, and a relative target-grid
energy estimate. No binary correlation theorem is assumed as an input.

## 1. Inputs and the pointwise semigroup theorem

We use two explicit inputs.

**ZF.** There is a fixed \(\sigma_0<1\) such that every Dirichlet
\(L\)-function is nonzero in \(\Re s>\sigma_0\), apart from its
possible principal pole at one. The bound is uniform in conductor and
height. GRH implies this hypothesis. The recent OpenAI manuscript
*The Quasi-Riemann Hypothesis* claims it in Theorem 1.1, with
\(\sigma_0=7/8\); its earlier Theorem 3.1 claims \(11/12\), which would
also suffice. Actual statements and the scope of the proof inspected
are recorded in the source ledger. That claim is not silently promoted
to an independently verified theorem here.

**Graph input G.** Import Theorem 3.5, the absolute row/column bound of
Lemma 3.4, and the uniform residue-replacement Lemma 3.6 in
*Weighted dilation graphs, smooth shifted primes and totient fibers*,
24 September 2026. Only the \(k=1\) physical root-moment theorem is
imported. Its arbitrary finite prime-group hypotheses are retained.
The fixed-AP ideal construction and large-\(k\) pairing needed here
are proved in the support, rather than presumed printed in that source.
The graph input is independent of the desired binary asymptotic.

Fix an even modulus \(L_0\), subgroups \(H_1,H_2\) of
\(G_0=(\mathbb Z/L_0)^\times\), and
\[
 \rho_i=|H_i|/|G_0|,\qquad
 \mathcal C_i=\langle p:p\nmid L_0,\ p\bmod L_0\in H_i\rangle,
 \qquad f_i=\mathbf1_{\mathcal C_i}.
\]
Let \(c_i>0\) be their Gamma-normalized one-variable constants:
\[
 F_i(s)=\sum_{n\ge1}f_i(n)n^{-s},\qquad
 c_i=\frac1{\Gamma(\rho_i)}
              \lim_{s\downarrow1}(s-1)^{\rho_i}F_i(s).
 \tag{B.1}
\]
These constants have explicit convergent Euler-product expressions from
\[
 F_i(s)=E_i(s)\prod_{\eta\in H_i^\perp}L(s,\eta)^{\rho_i}.
\]
Here characters are literal modulo \(L_0\). At a prime whose coset has
order \(r>1\) in \(G_0/H_i\), the correction is
\((1-p^{-rs})^{1/r}\); it is absolutely convergent to the right of
one half. The omitted bad Euler factors are carried by the literal
principal character. The positive real determination at \(s=1+\)
fixes every constant.

Write \(\mathcal E\) for the good primes forbidden on both sides:
\[
 \mathcal E=\{p\nmid L_0:p\bmod L_0\notin H_1,\
                                      p\bmod L_0\notin H_2\},
 \quad
 C_{\mathcal E}=\prod_{p\in\mathcal E}
                           \left(1-\frac1{(p-1)^2}\right),
 \quad
 J_{\mathcal E}(N)=\prod_{\substack{p\mid N\\p\in\mathcal E}}
                                      \frac{p-1}{p-2}.
 \tag{B.2}
\]
The prime two is excluded by \(L_0\), so all these factors are positive.

Fix coprime positive \(L_0\)-supported coefficients \(A,B\), and
fixed subsets \(U\subset H_1,V\subset H_2\). Define
\[
 \Gamma(N)=\frac{L_0}{|H_1||H_2|}
   \#\{(u,v)\in U\times V:Au+Bv=N\bmod L_0\}.
 \tag{B.3}
\]

**Theorem ConditionalBinarySemigroups.** Under ZF and graph input G,
for every such fixed data and every target with \(\Gamma(N)>0\),
\[
 \begin{split}
 &\#\{Ax+By=N:x\in\mathcal C_1,\ y\in\mathcal C_2,\
                     x\bmod L_0\in U,\ y\bmod L_0\in V\}\\
 &\qquad=
 \big(c_1c_2\Gamma(N)C_{\mathcal E}+o(1)\big)
 \frac{N J_{\mathcal E}(N)}
             {AB(\log N)^{\,2-\rho_1-\rho_2}} .
 \end{split}
 \tag{B.4}
\]
The \(o(1)\) is uniform in the target and its permitted valuations;
all displayed data are fixed. A zero finite density gives an
identically zero count. In particular no nesting between \(H_1,H_2\)
is required.

The proof below establishes the limit. A quantitative logarithmic
error requires keeping the indicated cutoff seminorm and parameter
orders through the support; no numerical onset is asserted.
If a subgroup is the full group, its centered model is exact and
that side needs no cancellation argument.

## 2. The complete proof and its error order

Put \(\mathcal L=\log N\),
\(y=\exp(\mathcal L^{2/5})\), and
\(\mathcal S=\{p:p\le y\}\). All bank families and density models
remain tied to this outer scale throughout dilations.

For the \(i\)-th rough core,
\[
 F_{i,\mathcal S}(s)=
 F_i(s)\prod_{\substack{p\le y\\p\ {\rm allowed\ on\ }i}}(1-p^{-s})
                    =(s-1)^{-\rho_i}G_{i,\mathcal S}(s),
\]
take a fixed, sufficiently deep Hankel expansion
\[
 M_i(t)=\sum_{j=0}^{K}
 \frac{G_{i,\mathcal S}^{(j)}(1)}{j!}
           \frac{(\log t)^{\rho_i-j-1}}{\Gamma(\rho_i-j)}.
 \tag{B.5}
\]
Extend it smoothly, boundedly and nonnegatively outside the relevant
rough bulk. For \(n=a v\), \(a=n_{\mathcal S}\), define
\[
 D_{i,U}(n)=
 [G_0:H_i]\,\mathbf1_{n\bmod L_0\in U}\,
            f_i(a)\frac{M_i(v)}{\prod_{p\le y}(1-p^{-1})},
 \qquad
 \xi_{i,U}=f_i\mathbf1_U-D_{i,U}.
 \tag{B.6}
\]
The same formula uses \(V\) on the second side.

The character projection is necessary. All twists in \(H_i^\perp\)
have the same principal branch, and a fixed residue filter produces
further fixed-period combinations of them. Formula (B.6) cancels
their full rough means. Keeping only a scalar density would leave
an uncancelled character contribution of main order.

The density and centering module proves the following, with no binary
correlation assumption. On every bulk interval the model is bounded
by one for sufficiently large \(N\), and all centered principal means
at logarithmic conductors/frequencies have arbitrarily strong
logarithmic savings after choosing \(K\). Its proof uses the classical
zero-free region and Siegel's theorem only at those logarithmic
parameters. In particular \(K\) is fixed after the desired accuracy;
it does not grow with \(N\).

Choose small and big prime banks only among \(p=1\bmod L_0\), with
ranges
\[
 [\exp(\mathcal L^{1/8}),\exp(\mathcal L^{1/4})],
 \qquad
 [\exp(\mathcal L^{1/3}),\exp(\mathcal L^{2/5})],
\]
respectively. Choose their initial reciprocal masses in a fixed
interior interval, and delete primes dividing \(N\). There are still
\(c_{L_0}\log\log N\) banks of each kind with masses bounded above
and below. The deleted harmonic mass is smaller than every power
of \(1/\log N\). Every original, model and centered function is
exactly invariant under multiplication by any remaining bank prime.

For clarity, here is the reduction without a hidden deweighting step.

1. A direct Selberg sieve on the reflected line gives a bank-record
   envelope of order
   \(N J_{\mathcal E}(N)(\log N)^{\rho_1+\rho_2-2}\).
   It applies to each of \(f_1f_2,D_1f_2,f_1D_2,D_1D_2\).
   The joint bank valuation law is the ordinary two-root law, since
   bank primes avoid \(NAB\). Its independence across banks and a
   Chernoff bound make the contribution with too few successful banks
   a smaller logarithmic order. Complete finite-record CRT counts
   and the endpoint budget are proved in the support.

2. Inclusion-exclusion on the successful banks reduces each centered
   correlation to a polynomial-logarithmic family with exactly one
   prime from each active bank on both sides. The total coefficient
   cost is a fixed power of \(\log N\). No factors \(\tau(n)\) remain:
   our actual indicators and models are already dilation invariant.

3. The AP ideal-kernel lemma repairs the interval-prime graph
   construction. Its character frame keeps one representative from
   each class
   \[
       \eta\sim\eta'
       \quad\Longleftrightarrow\quad
       \operatorname{cond}(\eta\overline{\eta'})\mid L_0 .
   \]
   These are exactly the characters equal on all bank products.
   The frame is positive, has norm at most two, and reproduces every
   small-denominator rational phase with the required error.
   The unmodified frame would count, for example, both \(1\) and
   \(\chi_{-4}\) as distinct constants on Gaussian banks. All precisions
   and coefficient exponents are fixed before the eventual mark length.

4. Graph input G now applies uniformly to the growing shift \(k=N\).
   After its CRT residue replacement, multiply every root residue by
   \(k^{-1}\) at each bank prime. This identifies each closed-walk
   expectation exactly with the \(k=1\) expectation. The replacement
   errors are uniform in offset values. In the endpoint pairing,
   partition intervals also by residue modulo \(k\); the operator
   preserves these classes, so the Gram-block dimension is
   \(O(D_{\max})\), rather than \(O(kD_{\max})\).
   This removes the large-shift norm loss.

5. Reflection causes no new sign assumption. Use the even extensions
   of the numerical indicators, supported on negative first and
   positive second endpoints. A contributing edge has
   \(\Psi(|n|/(N D))\ne0\), hence \(D\ll |n|/N\) and
   \(|n+ND|\ll|n|\). Insert a fixed-size second-endpoint cutoff before
   Mellin inversion. It equals one on all contributing edges and
   remains present after inversion. No change to shared-label removal
   or an additional \(D\)-localization is required.

6. The graph comparison has two independent free big-prime products.
   The elementary bilinear estimate makes their multiplier negligible
   off the preimages of small-denominator arcs. On those arcs, the
   frequencies lie on grids of denominator \(Nr\), \(r\le\mathcal L^C\),
   with physical width \(\mathcal L^C/Y\), where
   \(Y\asymp N D_1D_2\). The new target-grid energy lemma proves, for
   every prescribed \(T>0\),
   \[
    \sum_{j\bmod q}\int_{|\beta|\le\mathcal L^C/Y}
       |\widehat F(j/q+\beta)|^2\,d\beta
                \ll_T Y\mathcal L^{-T},
       \qquad q=Nr,\quad r\le\mathcal L^C.
    \tag{B.7}
   \]
   Here \(F\) is the actual centered, exactly banked coefficient from
   the comparison, not a model without its arithmetic filters.
   This is where ZF is indispensable in the present proof.

To explain the last estimate, separate a unique small-bank prime
\(p\asymp P\). All other bank weights remain in an arbitrary bounded
coefficient \(\gamma(m)\). For nonexceptional characters modulo any
actual target gcd profile, ZF gives
\[
 \sum_{\substack{p\in{\cal P}\\p\nmid q}}
          w(p/P)\eta(p)p^{it}
     \ll P^{\sigma_1}\mathcal L^C,
 \qquad \sigma_0<\sigma_1<1,
 \tag{B.8}
\]
uniformly in conductors up to \(Y^{O(1)}\) and \(|t|\le\mathcal L^C\),
with the fixed exponent required in the application.
The actual Mellin frequencies here are logarithmic: the factor \(N\)
in the target grid cancels the outer scale in its arc width.

The remaining polynomial has length \(M=Y/(dP)\), modulus \(q/d\), and
\[
          M/(q/d)=Y/(qP)\gg\mathcal L^T
\]
for every fixed \(T\). Character orthogonality bounds its mean square,
with no multiplicativity needed for \(\gamma\). Thus the nonexceptional
energy has an exponential-in-\(\mathcal L^{1/8}\) saving.
Exceptional characters are precisely the finitely many restrictions
from the fixed modulus. They are canceled by the full model (B.5)--(B.6).
The exact identity summing all \((n,q)=d\) profiles, and
\(\sum_{d\mid q}1/d\ll\log\log q\), retain arbitrary target valuations.
Large smooth records, target prime divisors above \(y\), and selected
prime squares are negligible in residue energy. The support proves
these facts before discarding them.

Cauchy--Schwarz with (B.7), followed by Parseval on the other endpoint,
makes each major comparison smaller than any requested log power.
The minor comparison has that saving by the multiplier bound.
Pay the finite comparison and reconstruction coefficients in order,
and obtain
\[
 \sum_{Ax+By=N}\Phi(Ax/N)\,\xi_{1,U}(x)w_2(y)
       =o\!\left(NJ_{\mathcal E}(N)
                    (\log N)^{\rho_1+\rho_2-2}\right)
 \tag{B.9}
\]
for \(w_2=f_2\mathbf1_V\) or \(D_{2,V}\).
Fixed bad-supported coefficients are embedded by
\(\widetilde f(n)=\mathbf1_{A\mid n}f(|n|/A)\).
Their Fourier sums are evaluated at \(A\alpha\); the target grid
denominator changes by \((q,A)\), a fixed factor, so (B.7) still applies.
This avoids incorrectly treating the embedded values as unit residues.

Swap the sides in (B.9) and use the exact identity
\[
 f_1\mathbf1_U\,f_2\mathbf1_V-D_{1,U}D_{2,V}
       =\xi_{1,U}f_2\mathbf1_V+D_{1,U}\xi_{2,V}.
\]
It remains to evaluate model times model. The fundamental lemma is
asymptotically exact at \(y=\exp(\mathcal L^{2/5})\) and any small
power level in \(N\). Its finite root density is (B.3).
A prime forbidden on only one side has relative factor one.
A prime forbidden on both sides contributes
\[
 \frac{1-(2-\mathbf1_{p\mid N})/p}{(1-1/p)^2};
\]
the product is \(C_{\mathcal E}J_{\mathcal E}(N)\).
The Gamma-normalized marginal constants give \(c_1c_2\).
The real line has length \(N/(AB)\); its bulk integral is
\(\int_0^1\Phi(t)\,dt\). There is no extra beta factor.
An upper sieve bounds the endpoint strips by \(O(\varepsilon)\)
of the full scale; let \(N\to\infty\), then \(\varepsilon\to0\).
This proves (B.4).

The parameter order is noncircular: choose the bank reconstruction
accuracy, then the graph error target, then its finite comparison
parameters, then the target-energy precision and Hankel depth.
For bounded indicators/models, exact-list endpoint moments are
\(O(X)\), independent of the eventual fixed mark length. Increasing
the model depth changes its threshold but not this bound.

## 3. Primitive Gaussian norms

Let
\[
 A_G=\frac12\prod_{p=3\bmod4}
                          \left(1-\frac2{p(p-1)}\right),
 \qquad
 J_G(N)=\prod_{\substack{p\mid N\\p=3\bmod4}}\frac{p-1}{p-2}.
\]
For \(P_G=\{x^2+y^2>0:\gcd(x,y)=1\}\), one has
\(P_G=C_4\sqcup2C_4\).

**Corollary ConditionalGaussianLeading.** Under the same two inputs,
\[
 \#\{a+b=N:a,b\in P_G\}\sim
 A_G J_G(N)\frac N{\log N}
 \begin{cases}
 1,&N\bmod8\in\{2,3,6,7\},\\
 1/2,&N\bmod8=4 .
 \end{cases}
 \tag{B.10}
\]
The count is identically zero on \(0,1,5\bmod8\).
For the odd core, the main is
\(A_GJ_G(N)N/\log N\) on \(N=2\bmod4\).

Indeed (B.1) gives \(c_{C_4}=1/(4K_{\rm LR})\), and
\(4c_{C_4}^2C_{\rm in}=A_G\).
The odd-plus-even contribution uses the two fixed coefficients \(1,2\);
the both-even contribution is the odd-core sum at \(N/2\).
These three disjoint cases give (B.10).

## 4. Every imaginary quadratic order

For a fixed order \(\mathcal O\) use the notation of
[the all-orders theorem](PRIMITIVE_NORM_ALL_ORDERS.md).
For \(p\in S=\{p:p\mid2\operatorname{disc}(\mathcal O)\}\), let
\(d_p(N)\) be the Haar density of the primitive local VALUE convolution.
Let \(\mathcal C\) be its good split core,
\(\alpha=\prod_{p\in S}(1-1/p)\), and let \(c_{\mathcal C}\),
\(C_{\rm in}\), \(J_{\mathcal O}\) have their
[previous normalization](PRIMITIVE_NORM_LEADING_CONSTANT.md).

**Corollary ConditionalAllOrderLeading.** Under ZF and graph input G,
uniformly through the locally admissible targets,
\[
 R_{P_{\mathcal O}}(N)\sim
 \frac{4c_{\mathcal C}^2}{\alpha^2}
 \left(\prod_{p\in S}d_p(N)\right)
 C_{\rm in}J_{\mathcal O}(N)\frac N{\log N}.
 \tag{B.11}
\]
Each bad density is an explicitly computable rational periodic
function; the remaining product converges absolutely.

Theorem (B.4) supplies exactly BinaryCoreCorrelation in the earlier
leading-constant note. Its proved uniform conductor-prefix tail and
Haar stratum identity sum the fixed-prefix asymptotics. Its lower-order
ring-class comparison transfers the result from the local shadow to
the primitive principal norm values. No uniqueness of a genus attached
to a conductor-divisible integer is assumed.

## 5. The \(C_{12}+C_{12}\) test case

Let \(C_{12}\) be supported on primes \(1\bmod12\), put
\[
 c_{12}=\frac1{\Gamma(1/4)}
      \lim_{s\downarrow1}(s-1)^{1/4}
                      \prod_{p=1\bmod12}(1-p^{-s})^{-1},
\]
and let \(\mathcal E_{12}=\{p:p\nmid6,\ p\ne1\bmod12\}\).

**Corollary ConditionalC12Leading.** Under the two explicit inputs,
\[
 \#\{a+b=N:a,b\in C_{12}\}
 \sim12c_{12}^2
 \prod_{p\in\mathcal E_{12}}\left(1-\frac1{(p-1)^2}\right)
 \prod_{\substack{p\mid N\\p\in\mathcal E_{12}}}\frac{p-1}{p-2}
             \frac N{(\log N)^{3/2}},
 \qquad N=2\bmod12.
 \tag{B.12}
\]
Other residue classes have count zero.

This is (B.4) with both subgroups equal to the identity subgroup modulo
12: \(\rho_1=\rho_2=1/4\), \(\Gamma=12\).
For a more explicit Euler normalization,
\[
 c_{12}=\frac1{\Gamma(1/4)}
 \left[
 \frac23 L(1,\chi_{-4})L(1,\chi_{-3})L(1,\chi_{12})
 \prod_{p\in\mathcal E_{12}}(1-p^{-2})^2
 \right]^{1/4}.
\]
The characters in this last formula are primitive. The factor \(2/3\)
comes from converting the four literal characters modulo 12 to them.

This conditional route uses new graph and all-conductor information.
It does not claim that the older square-root sieve barrier was removed
by its existing positive-restoration argument.

## 6. Assurance and remaining unconditional work

The source-level derivations, arithmetic normalization and parameter
interfaces above have internal written review. The infinite graph
input and the all-conductor zero-free input are distinct external
dependencies. The source ledger pins their precise editions and the
limits of inspection. Neither finite computation nor agreement among
agents independently proves those imported results.

In particular the unconditional binary-correlation gap remains open in
the repository. To discharge it by this route, the two imported deep
claims must be accepted with a complete applicable proof, or the
all-conductor prime estimate must be replaced by a proved sufficient
averaged estimate. The current classical short-progression theorem
alone gives too large an absolute variance error for this sparse scale.
No unconditional completion, historical priority or external review
is asserted.
