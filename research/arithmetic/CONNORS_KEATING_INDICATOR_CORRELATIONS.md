# Connors and Keating indicator correlations

**Scope.** Indicator correlations at fixed and logarithmically growing
shifts, Poisson variance at the mean-spacing scale, and fixed gap
frequencies under the imported zero-free and weighted-dilation theorems
of the [pointwise binary theory](PRIMITIVE_NORM_BINARY_CORRELATION.md).
Each numerical value is counted once. Sections 1--6 treat fixed shifts.
Section 7 proves uniformity and the moment limits. Section 8 treats
fixed consecutive gaps and binary patterns with two prescribed ones.

## 1. Statement and imported inputs

Let \(b(n)\) be one if the positive integer \(n\) is a sum of two
integer squares, and zero otherwise; set \(b(n)=0\) for \(n\le0\).
Zero is permitted as either coordinate. Let
\[
 K=\frac1{\sqrt2}\prod_{p\equiv3\ (4)}(1-p^{-2})^{-1/2}
\]
be the Landau--Ramanujan constant. For \(h\ne0\), put
\[
 W_2(h)=
 \begin{cases}
 1,&2\nmid h,\\
 2-3\,2^{-v_2(h)},&2\mid h,
 \end{cases}
 \qquad
 C(h)=\frac{W_2(h)}2
 \prod_{\substack{p\mid h\\p\equiv3\ (4)}}
       \frac{1-p^{-v_p(h)-1}}{1-p^{-1}}.                 \tag{CK.1}
\]
Valuations and divisibility of negative shifts refer to \(|h|\).

**Theorem CKFixedShift.** Under inputs ZF and G below, for every
fixed nonzero integer \(h\),
\[
 \sum_{1\le n\le x}b(n)b(n+h)
       \sim C(h)\frac{x}{\log x}.                       \tag{CK.2}
\]
For every \(\Phi\in C_c^\infty((0,\infty))\), the corresponding
smooth statement is
\[
 \sum_n\Phi(n/x)b(n)b(n+h)
 =C(h)\frac{x}{\log x}\int_0^\infty\Phi(t)\,dt
                         +o_{h,\Phi}(x/\log x).         \tag{CK.3}
\]
In particular \(C(1)=1/2\), \(C(2)=1/4\), \(C(3)=2/3\),
and \(C(4)=5/8\).

The hypotheses in the modular statement are:

- **ZF:** one fixed \(\sigma_0<1\) is a zero-free right
  half-plane for every Dirichlet \(L\)-function, with its principal
  pole allowed and with no conductor or height restriction.
- **G:** the physical root-moment theorem, Theorem 3.5, and
  Lemmas 3.4 and 3.6 of *Weighted dilation graphs, smooth shifted
  primes and totient fibers*, 24 September 2026, with the finite-group,
  ideal-norm, coefficient and complexity hypotheses retained.

These are the existing imported nodes
`HYP-BINARY-DIRICHLET-FIXED-ZEROFREE` and
`HYP-BINARY-WEIGHTED-DILATION-ROOT-MOMENT`.
Their exact primary editions and internal source-validation scope are
in [the input audit](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md).
The recorded ZF import uses only its Part I, Theorem 3.1,
\(\sigma_0=11/12\); the stronger Part II is unnecessary.
The conditional formulation displays the dependence even when those
sources are accepted as proved inputs. This note does not replace
their proofs or supply external human review of them.

Walsh's [v2 revision, 12 November 2026](https://arxiv.org/abs/2610.40439v2),
Theorem 1.1, printed p.2, proves fixed-shift integer-order divisor
correlations. Section 1.2.4, pp.5--6, explicitly proposes the
sums-of-two-squares indicator application. That direction is credited
to Walsh. The argument here supplies its indicator model, sparse
envelope, local constants and completed fixed-shift application.
It does not substitute a fractional order into his theorem.

Walsh also remarks on p.2 that his divisor theorem can be recovered
using classical estimates. His written prime-sum/high-frequency proof,
Lemma 5.1 and Section 6.4, uses the recent zero-free input. That remark
is not a completed classical replacement for the indicator proof here.
Its ZF dependence is retained; the remaining full-height estimate is
separate from the moderate-height classical range recorded in the
[shared comparison note](SPARSE_ARITHMETIC_MODEL_COMPARISON.md).

## 2. The indicator density and both principal branches

Write \(L=\log x\), \(y=\exp(L^{2/5})\),
\(\mathcal S=\{p:p\le y\}\),
\(V=\prod_{p\le y}(1-p^{-1})\), and \(E=\exp(\sqrt L)\).
All models and banks stay tied to this outer scale under dilation.
For \(n>0\), write \(n=a v\), where \(a=n_{\mathcal S}\)
and \((v,\prod_{p\le y}p)=1\).

The exact Euler factors are
\[
 F_p(s)=
 \begin{cases}
 (1-p^{-s})^{-1},&p=2\text{ or }p\equiv1\ (4),\\
 (1-p^{-2s})^{-1},&p\equiv3\ (4).
 \end{cases}
\]
Thus
\[
 F(s)=\sum b(n)n^{-s}
 =\left\{\zeta(s)L(s,\chi_4)(1-2^{-s})^{-1}
           \prod_{p\equiv3\ (4)}(1-p^{-2s})^{-1}\right\}^{1/2}.
                                                               \tag{CK.4}
\]
The positive real value for \(s>1\) fixes the branch.
The rough series is
\[
 F_R(s)=F(s)\prod_{p\le y}F_p(s)^{-1}
                    =(s-1)^{-1/2}H_R(s),
 \qquad
 \frac{H_R(1)}{\Gamma(1/2)}
                    =K\prod_{p\le y}F_p(1)^{-1}.        \tag{CK.5}
\]
For a sufficiently large fixed integer \(J\), use the density
\[
 M_J(t)=\sum_{j=0}^J\frac{H_R^{(j)}(1)}{j!}
                 \frac{(\log t)^{-1/2-j}}{\Gamma(1/2-j)}.
                                                               \tag{CK.6}
\]
On fixed power ranges in \(x\), the derivative bound
\(H_R^{(j)}(1)/j!\ll_j H_R(1)(C\log y)^j\) makes the
successive relative terms \(O_j(L^{-3j/5})\). Extend the density
smoothly and nonnegatively below those ranges, with the same upper
bound as its leading term. Define
\[
 D(n)=b(a)(1+\chi_4(v))\frac{M_J(v)}V,
 \qquad \xi(n)=b(n)-D(n).                                \tag{CK.7}
\]
On relevant scales \(D\ll L^{-3/10}\), hence \(|\xi|\ll1\).
For every split prime \(p\le y\),
\[
 b(pn)=b(n),\qquad D(pn)=D(n),\qquad \xi(pn)=\xi(n)       \tag{CK.8}
\]
exactly, including prime powers.

There are two principal branches. On odd representable integers
\(\chi_4(n)=1\). Consequently both the principal twist and the
\(\chi_4\) twist of the rough indicator coincide with the rough
indicator itself. A scalar model misses the second branch.

Here is the one-variable estimate required for centering. For
\(q\le L^B\), refine to \(\operatorname{lcm}(q,4)\).
The square of a twisted rough series has factors
\(L(s,\eta)L(s,\eta\chi_4)\), the removed small-prime factors,
and an Euler correction absolutely convergent for \(\Re s>1/2\).
The inert correction has the inverse sign to that for the odd
primitive Gaussian core; its absolute convergence and fixed-order
derivative bounds are the same. Its two principal branches occur
exactly when \(\eta=1\) or \(\chi_4\), inflated to the ambient
modulus. Every prime dividing that modulus is in \(\mathcal S\).

Apply the contour/Hankel argument of Module D, Section 2, and
Module C of [the binary support](PRIMITIVE_NORM_BINARY_CORRELATION_SUPPORT.md).
Explicitly, at height \(\exp(L^{1/10})\), the classical zero-free
region, with Siegel's theorem for a possible real character of
polylogarithmic conductor, permits the line
\(\Re s=1-c/\log y\). Its saving is
\(\exp(-cL^{3/5})\). At either principal branch, Taylor expansion
to depth \(J\) gives (CK.6), with relative remainder
\(O_J(L^{-3(J+1)/5})\). For any prescribed log power, choose
\(J\) after the conductor, periodic-weight and smooth-derivative
budgets. The unweighted rough population has principal density
\(V\); multiplying its model by \(1+\chi_4(v)\) gives exactly
the two matching principal means.

Therefore \(\xi\) has arbitrary logarithmic saving in the centered
low-frequency, fixed-small-prime-record progression sums required
in Walsh's Lemma 5.2 and Section 6.3. The record coefficients cost
only fixed powers of \(L\), since
\(\sum_{a\ \mathcal S\text{-smooth}}d_k(a)/a=L^{O_k(1)}\).
The discarded records \(a>E\) have total count
\(\ll x^C L^{O(1)}\exp(-cL^{1/10})\) on \([1,x^C]\),
by the tilted Euler moment used in Module D, (W8).
This is a one-variable assertion and assumes no correlation.

In the actual banked low-frequency expression, write \(n=av\)
before applying this assertion. Every bank weight, including the
Mellin-deformed exact-record weight \(W_t^u\), satisfies
\(W_t^u(av)=W_t^u(a)\): all bank primes are in \(\mathcal S\).
The coefficient \(b(a)W_t^u(a)\), bounded by a fixed power of
\(L\), is outside the inner rough sum. Its remaining arithmetic
test is periodic, such as \(e(cav/q)\), together with the stated
smooth and Mellin factors. Thus the centering argument does not
assert cancellation against an arbitrary bounded rough-variable
weight.

## 3. An envelope at the indicator correlation scale

Choose the small and big banks from split primes \(p\equiv1\pmod4\),
using the ranges \([\exp(L^{1/8}),\exp(L^{1/4})]\) and
\([\exp(L^{1/3}),\exp(L^{2/5})]\). Delete primes dividing \(h\).
There are \(\asymp\log L\) disjoint banks of each kind, with
reciprocal masses in a fixed interval such as \([1/2,2]\).
The prime progression has fixed positive density, so it supplies
these banks with the required fixed parameters.

For a bank prime, the law of the exact valuation record of
\((n,n+h)\) is
\[
 \mu_p(0,0)=1-2/p,\quad
 \mu_p(e,0)=\mu_p(0,e)=(1-1/p)p^{-e}\ (e\ge1),\quad
 \mu_p(e,f)=0\ (e,f\ge1).                                \tag{CK.9}
\]
Let \(\mu\) be the product over the banks. For every nonnegative
record function \(F_0\le L^B\), any interval \(I\subset[x,2x]\),
and \(w_1,w_2\in\{b,D\}\),
\[
 \sum_{n\in I} F_0(n_{\mathcal P},(n+h)_{\mathcal P})
                     w_1(n)w_2(n+h)
 \ll_h \frac{|I|}{L}\mathbb E_\mu F_0
       +xL^{O(1)}e^{-cL^{1/10}}+O(x^{1/2}).              \tag{CK.10}
\]

To prove this, first discard records larger than \(E\) by the
preceding tilted moment. Retained records \(r,s\) are coprime;
their divisibilities define one progression with \(|I|/(rs)+O(1)\)
terms in the chosen interval. The CRT endpoint bound is independent
of both interval endpoints.
Exactness at bank primes contributes precisely (CK.9): two roots
outside \(rs\), one remaining root at a prime of \(rs\).

At every inert prime, a necessary condition for \(b(m)=1\) is
\(v_p(m)\ne1\). Its forbidden classes modulo \(p^2\) number
\(p-1\). For \(p\nmid h\), the pair of forms forbids
\(2(p-1)\) classes. All bank-record slopes are units at inert
primes. Drop the finitely many inert primes dividing \(h\) for
an upper bound. The resulting sieve dimension is one, and the
prime-square terms differ from the usual two-root product only
by a convergent \(1+O(p^{-2})\) correction.

The upper sieve also permits these prime-square atoms directly.
If their forbidden probability is \(\rho_p\), put
\(g(p)=\rho_p/(1-\rho_p)\). For the Selberg denominator
\(G(z)=\sum_{d\le z}\mu^2(d)g(d)\), the normalized subset law
of the Euler product has
\(\mathbb E\log d=\sum_{p\le w}\rho_p\log p\ll\log w\).
Choosing a sufficiently small fixed power \(w\) of \(z\)
and applying Markov gives \(G(z)\gg\prod_{p\le w}(1+g(p))\).
The same argument allows the different cutoffs on model and actual
sides. Optimal Selberg weights and the CRT root counts give the
\(z^{4+o(1)}\) endpoint budget used next.

Apply the Selberg upper sieve with prime-square atoms, level a small
fixed power of \(x\), and \(w=x^\epsilon\) for a sufficiently
small fixed \(\epsilon>0\). The same elementary CRT error budget
as Module D, Section 5, is \(O(z^{4+\epsilon})\), before the
\(E^2=x^{o(1)}\) record sum. Choose, for example, \(z=x^{1/100}\),
then \(w=z^{1/16}\); the error is power smaller than \(x^{1/2}\).
For \(b b\), the nonbank denominator has order \(\log w\).
For \(D b\), it has order
\((\log y\log w)^{1/2}\); the bound \(D\ll L^{-3/10}\)
restores the factor \(L^{-1}\). For \(D D\), it has order
\(\log y\), and \(D^2\ll L^{-3/5}\) again restores
\(L^{-1}\). This proves (CK.10) for all four pairs.

This is the positive record envelope used by the shared comparison
theorem below. Its bank selection and reconstruction are part of
that theorem, rather than additional arithmetic estimates here.

## 4. Applying the fixed-shift comparison theorem

Use the fixed-shift case of the
[arithmetic-model comparison theorem](SPARSE_ARITHMETIC_MODEL_COMPARISON.md),
with \(m=4\), common prime class \(c=1\), \(f_1=f_2=b\), and
\(\mathcal S_x=x/L\). We verify its independent hypotheses.

* **(C1), bounded models.** The positive extension after (CK.6)
  gives \(0\le D\ll L^{-3/10}\), hence \(D\le1\) eventually
  for every fixed depth. The bound and the earlier endpoint moment
  exponents do not increase with the final depth.
* **(C2), a common invariant prime class.** Equation (CK.8) proves
  exact invariance at every selected split prime, including powers.
  Bank primes dividing \(h\) are deleted; for a fixed shift they
  are eventually below the bank ranges in any event.
* **(C3), all recordwise rough means.** Section 2 proves these
  against arbitrary permitted periodic tests, with both principal
  branches retained and the density depth chosen after the required
  accuracy. Every structured bank coefficient is a function of the
  small-prime record and is outside the inner rough-variable sum.
  No arbitrary bounded rough-variable test is being assumed.
* **(C4), a positive sparse record envelope.** The interval version
  (CK.10) proves this for every pair in \(\{b,D\}^2\), by explicit
  prime-square forbidden atoms and the uniform CRT endpoint bound.
  It uses only an upper sieve and does not assume cancellation.

There is no unit-support hypothesis at two. In a gcd profile for
\(Q=\operatorname{lcm}(q,4)\), the residual modulus \(Q/d\) may
lose factors of four. The shared proof therefore uses a primitive
character \(\eta\) with conductor dividing both four and \(Q/d\),
and the record phase \(\eta(r/d)\). It does not evaluate a literal
modulo-four character on a possibly even record after the gcd split.
The periodic rough means above cover every character that remains.

The fixed-shift clause is essential. Its major windows satisfy
\(|\beta|Y\asymp xL^B/|h|\), rather than a polynomial-logarithmic
bound. The shared proof treats these windows by the existing
additive-to-Dirichlet conversion and low/high Mellin-frequency
argument: (C3) supplies the low range, and the unique small-bank
prime, the same zero-free prime estimate and the ordinary Dirichlet
mean-value theorem supply the high range. The two free big products
make the cofactor length exceed the relevant height. Thus this
application does not invoke the reflected target-grid estimate at
a frequency range outside its statement.

The comparison theorem gives, for the centered term and either
\(w\in\{b,D\}\),
\[
 \sum_n\Phi(n/x)\xi(n)w(n+h)=o_h(x/L).                 \tag{CK.11}
\]
For a general fixed smooth bulk cutoff use partial summation and
a finite dyadic partition; the theorem's uniform interval error
has the same saving. Interchanging endpoints and replacing \(h\)
by \(-h\) gives the other orientation. The exact identity
\[
 b(n)b(n+h)-D(n)D(n+h)
       =\xi(n)b(n+h)+D(n)\xi(n+h)                      \tag{CK.13}
\]
reduces the smooth asymptotic to the arithmetic model product below.
There is no separate graph or energy hypothesis in these checks.

## 5. Evaluation of the model product

Define
\[
 \nu_2(n)=\mathbf1_{\operatorname{oddpart}(n)\equiv1\ (4)},
 \qquad
 b_{3,y}(n)=\mathbf1_{v_p(n)\text{ even for all }p\le y,
                                                    \ p\equiv3\ (4)}.
\]
After deleting \(n_{\mathcal S}>E\) and
\((n+h)_{\mathcal S}>E\), (CK.5)--(CK.7) give uniformly on
each fixed positive bulk interval
\[
 D(n)=\lambda_y\nu_2(n)b_{3,y}(n)(1+O_J(L^{-1/2})),
 \qquad
 \lambda_y=\frac{2K}{\sqrt L}
                     \prod_{\substack{p\le y\\p\equiv3\ (4)}}(1+1/p).
                                                               \tag{CK.14}
\]
Indeed \(\prod_{p\le y}F_p(1)^{-1}/V\) equals the displayed
inert-prime product. If \(b(a)=1\), the odd part of \(a\)
is \(1\pmod4\), so the rough character projection is exactly
\(2\nu_2(n)\). The discarded record error is superlogarithmic.

Let \(\delta_p(h)\) denote the Haar probability that both local
conditions hold at \(p\). An asymptotically exact sieve gives
\[
 \sum_n\Phi(n/x)\nu_2(n)\nu_2(n+h)b_{3,y}(n)b_{3,y}(n+h)
 =x\delta_2(h)\prod_{\substack{p\le y\\p\equiv3\ (4)}}\delta_p(h)
              \int\Phi+o_h(x/\log y).                   \tag{CK.15}
\]
The accuracy and its uniformity in \(y\) can be seen directly as
follows. Write the inert parts as \(d^2,e^2\), with \(d,e\)
supported on inert primes at most \(y\). Discard \(d>E\) or
\(e>E\), and dyadic factors exceeding \(E\). A union bound for
square divisors and geometric dyadic tails costs
\(O_h(x/E+\sqrt x)\). For retained prefixes, the divisibilities
and conditions \(n\equiv2^r\pmod{2^{r+2}}\),
\(n+h\equiv2^s\pmod{2^{s+2}}\), specify either an empty
stratum or one progression modulo
\(Q=\operatorname{lcm}(d^2,e^2,2^{r+2},2^{s+2})=x^{o(1)}\).
The residual two forms must be units at inert primes at most \(y\).
Each prime excludes zero, one, or two roots; if it excludes the
entire residue field the stratum is empty.

Use upper and lower fundamental-lemma weights of level
\(x^\epsilon\). The dimension is at most one, and
\(\log(x^\epsilon)/\log y\asymp L^{3/5}\to\infty\).
Thus the relative sieve error is superlogarithmic, and the total
CRT error over \(E^2(\log E)^2\) strata is
\(x^{O(\epsilon)+o(1)}=o(x/L^A)\) for every fixed \(A\),
after choosing a sufficiently small fixed \(\epsilon\).
The sum of prefix main masses is bounded: compatible inert
prefixes have \((d,e)^2\mid h\), and their masses are
\(\ll_h1/\operatorname{lcm}(d^2,e^2)\), with summable
dyadic factors. Summing the disjoint valuation profiles yields
exactly the local product in (CK.15). The omitted Haar tails are
\(O_h(E^{-1})\). This proves (CK.15) with the required absolute
error, without a moving-modulus fractional-order progression theorem.

The local densities are
\[
 \delta_p(h)=\frac p{p+1}(1-p^{-v_p(h)-1})
       \quad(p\equiv3\ (4)),\qquad
 \delta_2(h)=\frac{W_2(h)}4.                              \tag{CK.16}
\]
For an inert prime, partition \(v_p(n)<v_p(h)\), equality,
and \(v_p(n)>v_p(h)\), retaining even valuations; their geometric
sums give (CK.16). At two, if \(h\) is odd, one endpoint is odd.
The odd endpoint must be \(1\pmod4\); summing the permitted
valuation classes of the even endpoint gives \(1/4\).
If \(v_2(h)=1\), odd endpoints cannot both be \(1\pmod4\),
and the even endpoints give \(\delta_2(h)=\delta_2(h/2)/2=1/8\).
If \(4\mid h\), the odd contribution is \(1/4\), and the
even contribution is \(\delta_2(h/2)/2\). This recurrence gives
\(1/2-3/2^{v_2(h)+2}\) and proves the dyadic part of (CK.16).

By (CK.14)--(CK.16), the model coefficient is
\[
 \frac{4K^2\delta_2(h)}L
 \prod_{\substack{p\le y\\p\equiv3\ (4)}}
                    (1+1/p)^2\delta_p(h).
\]
For \(p\nmid h\), its local factor is \(1-p^{-2}\).
For \(p\mid h\), the ratio to that factor is
\((1-p^{-v_p(h)-1})/(1-p^{-1})\). Since
\(K^2\prod_{p\equiv3\ (4)}(1-p^{-2})=1/2\), the
coefficient tends to \(C(h)/L\). The missing Euler tail is
absolutely convergent. Together with (CK.13), this proves (CK.3).

Finally, nonnegative smooth upper and lower approximations on
dyadic intervals remove the smooth weight. The indicator pair
upper bound from (CK.10), used on ordinary dyadic scales with no
record weight, bounds the total contribution below \(\epsilon x\)
by \(O_h(\epsilon x/\log x)+o(x/\log x)\). Let the smoothing
error and then \(\epsilon\) tend to zero. The finitely many
endpoints with \(n+h\le0\) do not affect the limit.
This proves (CK.2).

## 6. Verification scope

The application uses the shared arithmetic-model comparison
theorem through the independent checks in Sections 2--4. The
rough model, prime-square upper sieve and local-factor evaluation
are the arithmetic arguments retained here. Its source dependence on
ZF and G remains explicit. In particular neither the existing
positive-sum theorem nor Walsh's integer-divisor theorem is used
as if it already asserted (CK.2).

The [exact finite checker](verify_connors_keating_local.py) checks local residue
densities and their valuation tails.
The command `python3 research/arithmetic/verify_connors_keating_local.py`
passed 640 residue-pair checks, 320 disjoint valuation-partition checks
and ten singleton checks for shifts 1 through 64 and primes
2, 3, 7, 11 and 19, using exact rational arithmetic. These computations
check the arithmetic factors; they do not verify
the infinite graph or zero-free estimates. External human review,
formal verification, an effective numerical onset, and historical
priority are separate from the written deduction.

Internal application review separately examined the local-model sieve
and the analytic comparison. The written proof records the resulting
mark-fiber, periodic-centering, cofactor-length and sharp-bank smoothing
conditions explicitly. No independent verification of the imported
theorems is inferred from that application review.

## 7. Logarithmically growing shifts and the second moment

The same ZF and G inputs give the following uniform form of the
Connors--Keating asymptotic. Use \(T_h=C(h)\) from (CK.1), and let
\[
 C_h(X)=\sum_{X<n\le 2X}b(n)b(n+h),\qquad
 L=\log X,\qquad \delta=\frac1{16384}.
\]

**Theorem CKLogarithmicShifts.** For every fixed \(A>0\), uniformly
for integers \(1\le |h|\le L^A\),
\[
 C_h(X)=T_h\frac X L\bigl(1+O_A(L^{-\delta})\bigr).
 \tag{CK.17}
\]
In this estimate \(X\) can be real. The proof uses the same model
\(\mathcal D\) as in Sections 2--5, at a fixed depth chosen for the
entire shift range.

### 7.1. Uniformity in the graph and frequency estimates

We use the interval comparison in Section 2 of the
[shared comparison note](SPARSE_ARITHMETIC_MODEL_COMPARISON.md),
with (C1)--(C3) as stated there and with (C4) uniform in the shift,
at a possibly shift-dependent scale \(\mathcal S_{X,h}\). Here is
the extra argument for the shift range. Set the common bank class
to \(1\bmod 4\).

Every bank prime exceeds \(L^A\) for large \(X\), so all the banks
are common to this range. In the [graph reduction](../../results/pointwise-binary-splitting-semigroups/graph.tex),
the physical displacement is \(h\) times an edge product.
Definition 3.3, Theorem 3.5 and Corollary 3.11 of G cover this
uniformly with the single exponent \(C_2=A\). The progression-bank
ideal is estimated uniformly in its real frequency and has no
dependence on \(h\). The shared-label reconstruction uses that
\(h\) is a unit at the bank primes and that both endpoint positions
are \(O(X^{1+o(1)})\). Those bounds are uniform here. On the support
of a dilated bulk cutoff the two endpoints remain comparable,
since their displacement divided by the dilated scale is
\(O(|h|/X)\).

In the [Fourier reduction](../../results/pointwise-binary-splitting-semigroups/prime-estimates.tex),
put \(U_*=U_1U_2\) for the product of the two dyadic scales and
\(Y=XU_*\). Approximate \(h\theta\bmod 1\) with Dirichlet
denominator ceiling \(\lfloor U_*/L^B\rfloor\), with \(B\)
chosen after the fixed comparison costs. If its denominator
\(r\) exceeds \(L^B\), the integer bilinear estimate gives
the prescribed minor-arc saving by taking \(B\) large. For
\(r\le L^B\), the lifts to \(\theta\bmod 1\) have reduced
denominators at most \(|h|r\), and the total number of arcs is
\[
 O\left(|h|\sum_{r\le L^B}\varphi(r)\right)=O(L^{A+2B}).
\]
Each lifted arc is inside a window of radius \(2L^B/U_*\)
about its center, since lifting divides the radius by \(|h|\).
Use this common radius for every shift.

The [wide-arc energy proof](../../results/pointwise-binary-splitting-semigroups/fixed-shift-comparison.tex)
now runs with window parameter \(W=U_*L^{-B}\) and
\(q\le L^{A+B}\). Its additive-to-Dirichlet inequality uses the
fixed radius constant \(2\). The tail cutoff remains
\(T_0=(Y/W)L^{B_1}=XL^{B+B_1}\) for the chosen fixed \(B_1\).
At low frequency the periodic test on the rough variable has
period dividing \(q\), so (C3) applies with one enlarged fixed
exponent. At high frequency, splitting the small-bank prime and
cofactor into residue classes costs at most \(q^2\) pairs. The
smooth prime estimate uses modulus \(\operatorname{lcm}(q,4)\)
and polynomial height, uniformly, and the unaltered cofactor ratio
satisfies
\[
 \frac{T_0}{Y/P}=\frac{P L^{B+B_1}}{U_1U_2}
 \le L^{B+B_1}\exp(O(L^{1/4})-2L^{1/3})=o(1).
\]
These are all the modulus and height costs in the energy estimate.
For every fixed \(D>0\), its bound \(O(YL^{-D})\) on each
window follows by choosing the low-frequency model depth after
the residue, smoothing and Mellin budgets. Cauchy--Schwarz and
Parseval then absorb \(O(L^{A+2B})\) windows, as well as the
fixed comparison costs, by choosing \(D\) before that depth.
The same bounds apply with the centered endpoint in either slot,
using the opposite sign of the displacement for the reversed pairing.

The explicit reconstruction in
[assembly.tex](../../results/pointwise-binary-splitting-semigroups/assembly.tex)
now uses its existing rare-bank probability \(O(L^{-2\delta})\)
and smooth endpoint transitions of width \(XL^{-\delta}\).
Both costs use the uniform (C4) bound at \(\mathcal S_{X,h}\).
The scale is allowed between \(XL^{-\kappa_0}\) and
\(XL^{\kappa_0}\) for a fixed \(\kappa_0\), so the earlier
arbitrary logarithmic savings absorb the reconstruction costs
with one accuracy for all shifts. The resulting interval comparison
error is \(O_A(\mathcal S_{X,h}L^{-\delta})\).

### 7.2. The sieve scale and the model tails

Use the envelope of Section 3 with
\[
 J(h)=\prod_{\substack{p\mid h\\p\equiv3\ (4)}}\frac{p-1}{p-2},
 \qquad \mathcal S_{X,h}=\frac{XJ(h)}{L}.
\]
The inequalities \(1\le J(h)\le |h|\) put this scale in the range
required above. To get the displayed factor \(J(h)\), when an inert
prime divides \(h\) and both endpoints are being sieved at that prime,
keep the necessary condition \(v_p(n)\ne1\) at one endpoint.
At other inert primes with both tests, the excluded classes are
disjoint. With \(\rho_p=(p-1)/p^2\), the cost of using one test
instead of two is therefore at most
\[
 \frac{1-\rho_p}{1-2\rho_p}\le\frac{p-1}{p-2}.
\]
At primes with just one active test keep that test, even if the prime
divides \(h\). The nonbank product without these corrections
has the same three orders as in Section 3 by Mertens' theorem
in the class \(3\bmod4\). They are \(1/\log w\),
\(1/\sqrt{\log y\log w}\) and \(1/\log y\) for two actual
endpoints, one model endpoint and two model endpoints,
respectively. Each model factor is \(O(L^{-3/10})\).
The exact bank law is unchanged since \(p\nmid h\) at every bank
prime. The Selberg denominator proof in Section 3 is uniform
since all forbidden fractions are at most \(2/p\), and the
CRT endpoint bound only counts roots. With level \(z=X^{1/100}\),
the \(O(z^{4+o(1)})\) error per retained record, summed over at
most \(E^2\) records, is power smaller than \(X^{1/2}\).
For omitted records the one-variable Rankin bound applies on
\([1,3X]\) to either endpoint. This gives (C4) with constants
independent of \(h\).

For \(v=v_p(h)\ge1\), the factor of \(T_h\) at an odd inert
prime is \(S_p(v)=1+1/p+\cdots+1/p^v\). Uniformly in \(v\),
\[
 1\le \frac{(p-1)/(p-2)}{S_p(v)}
 \le\frac{(p-1)/(p-2)}{1+1/p}
 =1+\frac{2}{(p-2)(p+1)}.
\]
The upper bounds have a convergent product. The dyadic factor
\(W_2(h)/2\) lies in \([1/4,1]\). Consequently
\(T_h\asymp J(h)\) with absolute constants.

For the model product, the full small-prime records exceeding
\(E=\exp(\sqrt L)\) have total integer count
\(O(XL^{O(1)}\exp(-cL^{1/10}))\), uniformly by the same
one-variable bound on \([1,3X]\). On the other records
(CK.14) holds with its uniform relative \(O(L^{-1/2})\).
For (CK.15) use the half-open interval and the inert square
prefixes \(d^2,e^2\), retaining \(d,e\le E\) and dyadic powers
at most \(E\). Either endpoint's square-prefix tail is bounded
by
\[
 \sum_{E<d\le\sqrt{3X}}\left(\frac X{d^2}+1\right)
       \ll X/E+\sqrt X.
\]
A dyadic tail consists of multiples of the first power of two
above \(E\) and costs \(O(X/E+1)\). The corresponding Haar
tails have mass \(O(E^{-1})\) by the same divisibility union
bound, since translation by \(h\) preserves Haar measure.

In each retained compatible prefix progression \(n=a+Qt\),
the normalized endpoints \(n/d^2\) and \((n+h)/e^2\) are
integral linear forms in \(t\). At an inert prime at most \(y\),
the two unit conditions exclude at most two roots in \(t\), or
the cell is empty if an endpoint is identically zero there.
This holds also at primes dividing \(Q\). The sieve dimension
bound follows uniformly from the product of
\((1-2/p)^{-1}\) over primes \(p\equiv3\bmod4\), using
Mertens' theorem in that class. With fundamental-lemma level
\(D_0=X^{1/100}\), the relative error is
\(O(\exp(-cL^{3/5}))\), by the
[fundamental lemma](https://books.google.com/books?id=GJSKAwAAQBAJ)
of Friedlander and Iwaniec, Lemma 6.11.
The CRT remainder at a squarefree divisor \(m\) is
\(O(2^{\omega(m)})\) per cell. Over all
\(O(E^2(1+\log E)^2)\) cells, the endpoint error is at most
\(X^{1/50+o(1)}\).

The main coefficient for each cell is its exact Haar mass,
because restriction to the progression fixes the prefix and
the residual root densities multiply. The retained valuation
cells are disjoint, hence the sum of their Haar masses is at
most one. Adding the omitted Haar tails gives
\(\delta_2(h)\prod_{p\le y,\ p\equiv3\ (4)}\delta_p(h)\),
with all the errors above smaller than every fixed logarithmic
power relative to \(X\) in the count. Every prime divisor of
\(h\) lies below \(y\), so the only missing Euler factors in
(CK.14)--(CK.16) are \(1-p^{-2}\) at \(p>y\). Their tail
has relative error \(O(y^{-1})\). This proves
\[
 \sum_{X<n\le2X}\mathcal D(n)\mathcal D(n+h)
    =T_h\frac X L\bigl(1+O_A(L^{-1/2})\bigr).
\]
The uniform comparison gives an error of size
\(O_A(XJ(h)L^{-1-\delta})\).
Using \(J(h)\ll T_h\) proves (CK.17).

### 7.3. The mean, second factorial moment and variance

**Corollary CKPoissonVariance.** Under ZF and G, for every fixed
\(\lambda>0\), let \(H=\lfloor\lambda\sqrt L/K\rfloor\) and
let \(t\) be uniform on \([X,2X]\cap\mathbb Z\). Then
\(Z_X(t)=\sum_{j=1}^H b(t+j)\) satisfies
\[
 \mathbb E Z_X\longrightarrow\lambda,\qquad
 \mathbb E[Z_X(Z_X-1)]\longrightarrow\lambda^2,\qquad
 \operatorname{Var}(Z_X)\longrightarrow\lambda.
 \tag{CK.18}
\]

Write \(N_X=|[X,2X]\cap\mathbb Z|=X+O(1)\) and
\(B(x)=\sum_{n\le x}b(n)\). Translating the \(j\)-th inner interval
changes the count from \(B(2X)-B(X)\) by \(O(j+1)\). Landau's
asymptotic gives
\[
 \mathbb E Z_X=
 \frac{H}{N_X}\bigl(B(2X)-B(X)\bigr)+O(H^2/X)
 =\frac{HK}{\sqrt L}(1+o(1))+O(H^2/X).
\]
For a pair \(j,j+h\), the same translation changes the inner
pair count from \(C_h(X)\) by \(O(j+1)\). Accounting for
both orders gives
\[
 \mathbb E[Z_X(Z_X-1)]
 =\frac{2}{N_X}\sum_{1\le h<H}(H-h)C_h(X)+O(H^3/X).
\]
In our normalization Smilansky's
[Theorem 1](https://arxiv.org/abs/1208.2419v2) gives
\[
 2\sum_{1\le h<H}(H-h)T_h=K^2H^2+O_\varepsilon(H^{1+\varepsilon}).
\]
In particular (CK.17) and positivity imply
\[
 \sum_{1\le h<H}(H-h)\left(C_h(X)-T_h\frac X L\right)
 \ll_\lambda XL^{-\delta}.
\]
The second factorial moment is therefore \(K^2H^2/L+o(1)\).
The variance follows by adding the first moment and subtracting
its square.

For the original sampling range \(1\le t\le X\), one can use
\(\sum_{1\le n\le X}b(n)b(n+h)\) in the same proof. To see that
it satisfies (CK.17) with the same main term, use dyads
\((Y,2Y]\) with \(Y=X/2^j\ge \sqrt X\) and bound the remaining
\(O(\sqrt X)\) terms trivially. For large \(X\) all these dyads
have \(|h|\le(\log Y)^{A+1}\).
Their comparison remainders sum to \(O_A(T_h XL^{-1-\delta})\).
The main terms sum using
\(\sum_jY/\log Y=X/L+O(X/L^2)\), which follows from
\(\left|(\log Y)^{-1}-L^{-1}\right|\ll j/L^2\) and
\(\sum_{j\ge1}j2^{-j}<\infty\). The unused mass at the lower
endpoint is absorbed since \(T_h\ge1/4\).
The argument permits \(H+1\) terms with the same limiting
\(HK/\sqrt L\). Taking an inclusive interval gives
\(\mathbb E Z_X^2\to\lambda+\lambda^2\) in the sampling of
Gorodetsky--Rodgers [Conjecture 1, \(k=2\)](https://arxiv.org/abs/1810.06002v4).

## 8. Fixed consecutive gaps and patterns with two ones

**Corollary CKFixedPatterns.** Assume ZF and G. Let \(u<v\) be
fixed integers and let \(\mathcal Z\) be a fixed finite subset of
\(\mathbb Z\setminus\{u,v\}\). Then
\[
 P_{u,v,\mathcal Z}(X)
 :=\sum_{1\le n\le X} b(n+u)b(n+v)
                   \prod_{z\in\mathcal Z}(1-b(n+z))
 \sim T_{v-u}\frac{X}{\log X}.
 \tag{CK.19}
\]
In particular, if \(s_1<s_2<\cdots\) are the positive sums of two
squares, then for every fixed \(h\ge1\),
\[
 \#\{j:s_j\le X,\ s_{j+1}-s_j=h\}
 \sim T_h\frac{X}{\log X}.
 \tag{CK.20}
\]

The Selberg sieve bound of [Nowak, Section 2](https://doi.org/10.2298/PIM0591071N),
also recorded by
[Freiberg--Kurlberg--Rosenzweig, Section 2.1](https://arxiv.org/abs/1701.01157v1),
gives, for three distinct fixed integer offsets \(u,v,z\),
\[
 \sum_{1\le n\le X}b(n+u)b(n+v)b(n+z)
 \ll_{u,v,z}\frac{X}{(\log X)^{3/2}}.
\]
Indeed, take the normal extension \(\mathbb Q(i)/\mathbb Q\) of
degree two and the three monic linear forms with offsets \(u,v,z\).
Their pairwise determinants are nonzero. The finite exceptional-prime
product depends only on the offsets. In the Gaussian field the
integral-ideal norm indicator is exactly \(b\).

Let \(P_{u,v,\varnothing}(X)\) be the pair count. Each excluded pair
has \(b(n+z)=1\) for some \(z\in\mathcal Z\). The union bound gives
\[
 0\le P_{u,v,\varnothing}(X)-P_{u,v,\mathcal Z}(X)
 \le\sum_{z\in\mathcal Z}\ \sum_{1\le n\le X}
                      b(n+u)b(n+v)b(n+z)
 \ll_{u,v,\mathcal Z}\frac{X}{(\log X)^{3/2}}.
\]
Translating \(n\) by \(u\) changes at most \(O_u(1)\) endpoint
terms when applying (CK.2) at shift \(v-u\). Since \(T_{v-u}>0\),
the pair asymptotic gives (CK.19). Taking \(u=0\), \(v=h\) and
\(\mathcal Z=\{1,\ldots,h-1\}\) gives exactly the count in (CK.20).
For \(h=1\) the interior is empty.

These are frequencies for each fixed gap size. Bounds on gap size
as \(s_j\) grows are the subject of
[Erdős Problem 222](https://www.erdosproblems.com/222).
