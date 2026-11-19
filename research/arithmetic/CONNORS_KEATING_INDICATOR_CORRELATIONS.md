# Connors and Keating indicator correlations

**Scope.** A fixed-shift indicator-correlation argument using the same
imported zero-free and weighted-dilation theorems as the
[pointwise binary theory](PRIMITIVE_NORM_BINARY_CORRELATION.md).
The precise inputs are stated below. The argument concerns values,
each counted once, and does not use a representation-weighted
correlation in place of an indicator correlation. Growing shifts and
higher correlations are outside its scope.

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

Walsh's [12 November revision](https://arxiv.org/html/2610.40439v2),
Theorem 1.1, proves integer-order divisor correlations. Its Section
1.2.4 identifies the indicator problem as further work. The proof
below supplies an indicator density, an envelope at scale
\(x/\log x\), and the fixed-shift application of the graph comparison.
It does not substitute a fractional order into that theorem.

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
record function \(F_0\le L^B\) and \(w_1,w_2\in\{b,D\}\),
\[
 \sum_{n\asymp x} F_0(n_{\mathcal P},(n+h)_{\mathcal P})
                     w_1(n)w_2(n+h)
 \ll_h \frac{x}{L}\mathbb E_\mu F_0
       +xL^{O(1)}e^{-cL^{1/10}}+O(x^{1/2}).              \tag{CK.10}
\]

To prove this, first discard records larger than \(E\) by the
preceding tilted moment. Retained records \(r,s\) are coprime;
their divisibilities define one progression of length \(O(x/(rs))\).
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

In the product law, a bank has exactly one distinct prime at each
endpoint with probability
\(V_g^2e^{-2V_g}+o(1)\), bounded below by a positive constant.
Chernoff and (CK.10) make the contribution with too few successful
small or big banks \(O_h(xL^{-1-\eta})\), for some \(\eta>0\).
Selected bank prime squares cost less than every log power.
The exact inclusion-exclusion and normalization in Module D,
Section 6, then reduce centered correlations to banked correlations,
with a total coefficient cost \(L^{C_0}\). There is no divisor
deweighting: (CK.8) already gives exact invariance.

## 4. The fixed-shift graph comparison

For clarity, the comparison needed after reconstruction is
\[
 \sum_n\Phi(n/x)\xi(n)w(n+h)=o_h(x/L),
                         \qquad w\in\{b,D\}.             \tag{CK.11}
\]
Its reflected counterpart follows by interchanging the endpoints
and replacing \(h\) by \(-h\).

Here are the graph and frequency interfaces. Apply the construction
of Walsh, Proposition 6.1 and Section 7, to the exactly banked
bounded invariant coefficients in (CK.8). The replacement for the
unrestricted-prime ideal kernel is Module I of the binary support,
with fixed modulus four. Characters differing by \(\chi_4\) have
the same restriction to these banks. Keep one representative per
restriction class, rather than counting the two constants twice.
The conditional Gram frame is positive with norm at most two;
the signed probe construction has polynomial-logarithmic mass
independent of the eventual fixed mark length.

The imported root moment G applies to this ideal kernel. Every
bank prime avoids \(h\), so the uniform residue replacement and
the unit change of variable by \(h^{-1}\) identify its physical
root moments with the \(k=1\) input. Here \(h\) is fixed; no
growing-shift theorem, reflection of a positive-sum line, or
target-sized denominator assertion is being used. On each physical
dyadic interval of length \(O(X)\), the absolute endpoint moments
of fixed order are \(O_r(X)\), with log exponent independent of the
eventual fixed mark length. For a fixed ordered prime list of
product \(P=x^{o(1)}\), counting its multiples gives \(O(X/P)\);
summing the normalized reciprocal list masses gives this bound.
The numerical first-endpoint factor is a function of the integer
position alone, constant on the mark fiber, as required by the
root-moment endpoint consequence. Indicator and model factors are
bounded and do not enlarge these moments.
The coefficients retain the big-label restriction of G: dependence
on freed small labels is excluded before forming comparisons.
The permutation averages, shared-label deletion, collision budget
and endpoint pairing are those of Module I and Walsh Section 7.5.
They depend on the record and graph structure and on the endpoint
moment bound, not on integer divisor orders.

Consequently, the banked fixed correlation is a signed polynomial
family of comparisons with two independent free big-prime products,
up to \(O_A(xL^{-A})\). The elementary bilinear multiplier estimate
in Walsh Section 6.1 is uniform in their bounded big-label weights;
it gives arbitrary logarithmic saving away from small-denominator
arcs. For fixed \(h\), the resulting major denominators are
\(O_h(L^C)\). Walsh Sections 6.2--6.4 reduce their contribution
to low-frequency centered rough means and high-frequency prime
sums, with only polynomial-logarithmic coefficient losses.

The low-frequency means are supplied by Section 2, including both
principal branches. For the high-frequency part extract the unique
small-bank prime. Exact invariance makes its cofactor coefficient
independent of that prime; its other record restrictions remain in
the bounded cofactor sequence. On a smooth dyadic interval, ZF
gives the estimate (P1), Module P of the binary support:
\[
 \sum_p\eta(p)w_0(p/P)p^{it}
 =\mathbf1_{\eta\text{ principal}}
          \int_0^\infty w_0(u/P)u^{it}\frac{du}{\log u}
       +O(P^{\sigma_1}L^{C_1}),                          \tag{CK.12}
\]
uniformly for \(q\le L^C\), \(|t|\le x^C\), and fixed
\(\sigma_0<\sigma_1<1\). The proof shifts the smooth
von-Mangoldt Mellin integral; the logarithmic derivative costs
\(O(\log(q(|t|+3)))\), so high height costs only powers of \(L\).

Expand the split restriction using \((1+\chi_4(p))/2\).
Both possible principal integrals in (CK.12) decay faster than any
fixed power of \(|t|\) after integration by parts. Above a large
polylogarithmic threshold they are negligible. The error, divided
by the prime mass, is
\(P^{\sigma_1-1}L^{O(1)}\), smaller than every negative log
power because \(P\ge\exp(L^{1/8})\).
As in Walsh Section 6.4, first smooth sharp bank endpoints at
width \(P L^{-K_0}\). Their normalized boundary mass is
\(O(L^{-K_0}+P^{-1})\), uniformly in \(t\), and the cofactor
mean square pays only a fixed log power. Choose \(K_0\) after
that loss. The sharp partial-summation estimate (PS2), whose error
grows linearly in \(|t|\), is not used at these heights.

The cofactor mean square has the required length budget. With
free big products \(H_1,H_2\), physical length
\(Y\asymp xH_1H_2\), small-prime size \(P\), and height
\(T_0=xL^C\), the cofactor length is \(R\asymp Y/P\).
Then
\(T_0/R\ll P L^C/(H_1H_2)=o(1)\), because the small-prime
logarithmic size is at most \(L^{1/4}\), whereas each big
product has logarithmic size at least \(L^{1/3}\).
The Dirichlet mean-value theorem therefore bounds its normalized
energy by a fixed power of \(L\). This pays the character,
tuple and smoothing losses before the prime saving is applied.

Choose the requested saving after the reconstruction and comparison
coefficient exponents, the density depth after the low-frequency
budget, and the fixed graph mark length after the ideal precision.
The bounded endpoint-moment exponent does not increase with that
mark length. Each comparison now has enough saving to pay every
coefficient. Adding the rare-bank error proves (CK.11).

The exact identity
\[
 b(n)b(n+h)-D(n)D(n+h)
       =\xi(n)b(n+h)+D(n)\xi(n+h)                        \tag{CK.13}
\]
therefore reduces (CK.3) to evaluation of the model product.

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

The new application uses the previously proved projected-density,
fixed-progression ideal-kernel and smooth-prime modules with the
specific changes stated in Sections 2--4. Its source dependence on
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
