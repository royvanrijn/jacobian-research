# Averaged triples at the mean spacing of sums of two squares

Let \(b(n)\) indicate positive sums of two squares, with zero coordinates
allowed, and let \(K\) denote the Landau–Ramanujan constant. Fix
\(\lambda>0\), put \(L=\log X\) and \(H=\lfloor\lambda\sqrt L/K\rfloor\),
and choose \(t\) uniformly among the \(N_X\) integers in \((X,2X]\).
Set
\[
 Z=\sum_{j=1}^H b(t+j),\qquad
 W=\sum_{j=1}^H D_X(t+j),\qquad
 \xi_X(n)=b(n)-D_X(n).
\]
Here \(D_X\) is the nonnegative two-character Hankel model (CK.7)
at the *same outer scale \(X\)*, with a sufficiently deep fixed
truncation. The inputs ZF and G are the Dirichlet zero-free and
weighted-dilation inputs of the logarithmic-shift pair theorem.
[The indicator model and logarithmic-shift theorem](CONNORS_KEATING_INDICATOR_CORRELATIONS.md)
· [Exact analytic inputs](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md)

**Proposition (averaged mixed triples and the centered remainder).**
Under ZF and G, the expected sum over ordered distinct triples of
shifts converges to \(\lambda^3\) for each of the three products with at least
one model factor, namely \(bbD_X\), \(bD_XD_X\) and \(D_XD_XD_X\).
Consequently
\[
 \mathbb E(Z)_3
  =\lambda^3+
   \frac{6}{N_X}\sum_{1\le u<v<H}(H-v)
      \sum_{X<n\le2X}\xi_X(n)\xi_X(n+u)\xi_X(n+v)
   +o_\lambda(1),                                             \tag{1}
\]
where \((Z)_3=Z(Z-1)(Z-2)\).
Also \(\mathbb E(Z)_3\ll_\lambda1\) and
\[
 \mathbb E\sum_{\substack{i,j,k\le H\\\mathrm{distinct}}}
            |\xi_X(t+i)\xi_X(t+j)\xi_X(t+k)|
       =\mathbb E(Z)_3+7\lambda^3+o_\lambda(1).                \tag{2}
\]
The positive sieve estimate giving the boundedness of
\(\mathbb E(Z)_3\) is independent of ZF and G. The signed sum in
(1) remains the open step.

Freiberg, Kurlberg and Rosenzweig already prove the mean-one theorem
for the tuple singular series in Proposition 1.3. Their Theorem 1.2
uses averaged tuple asymptotics through all fixed orders to obtain
Poisson laws. The centered remainder in (1) concerns the arithmetic
indicator and its binary comparison model at order three.
[Freiberg–Kurlberg–Rosenzweig, Theorem 1.2 and Proposition 1.3](https://arxiv.org/abs/1701.01157)

## A uniform upper sieve for the window

The fixed-tuple prime-square sieve extends uniformly to distinct shifts
\(a_1,\ldots,a_k\in[1,H]\) for any fixed \(k\). This uses the same
valuation-one obstruction as the upper sieve for fixed norm patterns.
[Prime-square upper envelope for fixed norm patterns](ERDOS_675_TRANSLATION_PATTERNS.md#8-a-higher-correlation-upper-bound-with-no-existence-consequence)

Write \(y=\exp(L^{2/5})\), and let \(q\) of the \(k\) factors \(w_i\)
be actual indicators \(b\), with the others equal to \(D_X\).
The model is \(O(L^{-3/10})\) and its support requires even valuations
at all inert primes at most \(y\). Each actual indicator requires even
valuations at all inert primes. Put
\[
 \Delta(\mathbf a)=\prod_{i<j}|a_i-a_j|,\qquad
 J_k(\Delta)=\prod_{p\mid\Delta}(1+2k/p).
\]
Use the empty product for \(k=1\). Then, uniformly in the tuple,
\[
 \sum_{X<t\le2X}\prod_{i=1}^k w_i(t+a_i)
       \ll_k X L^{-k/2} J_k(\Delta) + X^{1/40}.               \tag{3}
\]

To see this, drop the finitely many primes at most \(2k\) and the primes
dividing \(\Delta\). At an inert prime \(p\le y\) every slot excludes
valuation one. At \(y<p\le z\), where \(z=X^{1/(100s)}\) for a
sufficiently large fixed \(s=s(k)\), keep the tests on the \(q\)
actual slots. The excluded classes for distinct active slots
are disjoint because their shifts are distinct modulo the
retained \(p\). Thus the removed proportion is
\(j_p(p-1)/p^2\) with \(j_p=k\) or \(q\), and the number of excluded
classes modulo \(p^2\) is at most \(kp\).

Apply the fundamental lemma at level \(D_0=X^{1/100}\), with sieve
dimension at most \(k\). The dimension estimate follows from
\(j_p(p-1)/p^2\le k/p\) and the prime Mertens bound, uniformly even
after dropping the collision primes. The retained proportions at
\(p>2k\) are bounded below. The CRT remainder for a squarefree atom
\(d\) is \(O(k^{\omega(d)}d)\), independent of the sizes of the
shifts. Richert's remainder weight \(3^{\omega(d)}\) thus gives
total error \(O_k(D_0^2(\log D_0)^{O_k(1)})\), which is absorbed
by the displayed \(X^{1/40}\). The same argument on a subinterval
has its length in the main-term bound and the same endpoint error.
[Richert, Lectures on Sieve Methods, Theorem 11.7](sources/richert-sieve-methods-1976.pdf)

For the main term, a retained factor differs from
\((1-1/p)^{j_p}\) by \(1+O_k(p^{-2})\). Restoring the missing
comparison factors at \(p\mid\Delta\), \(p>2k\), costs at most
\((1-k/p)^{-1}\le1+2k/p\) at each such prime. Mertens' theorem
for primes \(3\bmod4\) gives the unweighted product bound
\[
 \ll_k J_k(\Delta)(\log y)^{-(k-q)/2}(\log z)^{-q/2}.
\]
The \(k-q\) model weights supply \(L^{-3(k-q)/10}\), giving (3)
since \(\log y=L^{2/5}\) and \(\log z\asymp_k L\).

The collision loss in (3) has bounded mean over the ordered
distinct tuples. For \(k\ge2\) put \(m=\binom{k}{2}\) and
\(j(d)=\prod_{p\mid d}(1+2k/p)\). We have
\[
 J_k(\Delta)\le\prod_{i<j}j(|a_i-a_j|).
\]
The expansion of \(j(d)^m\) is \(\sum_{r\mid d}g(r)\), where
\(g\ge0\) is supported on squarefree \(r\), is multiplicative there,
and \(g(p)=(1+2k/p)^m-1\ll_k 1/p\). Consequently
\[
 \sum_{d\le H}j(d)^m
    \le H\prod_p(1+g(p)/p)\ll_k H.
\]
For each pair of slots, the sum of \(j(|a_i-a_j|)^m\) over distinct
tuples is at most \(H^{k-2}\sum_{d<H}2(H-d)j(d)^m\ll_k H^k\).
Hölder now gives \(\sum_{\mathbf a\ {\rm distinct}}J_k(\Delta)
\ll_k H^k\). Combining this with (3) proves that every fixed-order
mixed sum over distinct shifts has bounded expectation. In
particular \(\mathbb E Z^4+\mathbb E W^4\ll_\lambda1\), since
repeated shifts in the fourth powers reduce to at most four
distinct slots using \(0\le b,D_X\le1\).

## Model concentration and the exact signed target

The mean of \(W\) tends to \(\lambda\). Indeed, summing the
one-variable centering estimate (C3) on a smooth dyadic cutoff
gives \(\sum_{X<n\le2X}(D_X(n)-b(n))=o(X/\sqrt L)\).
For this use, approximate the endpoints by smooth transitions
of length \(XL^{-B}\) with fixed \(B>1\). Their derivatives
are polylogarithmic at scale \(X\), and boundedness controls
the transition errors. The sum of reciprocal retained records
costs only a power of \(L\), and the large-record tail is
superlogarithmically small. Choose the model depth after
these costs. Landau's mean for \(b\) and translation of the
\(H\) endpoints now give \(\mathbb E W\to\lambda\).

For the second moment of \(W\), the uniform model-pair
evaluation in Section 7.2 and Smilansky's weighted mean in
Section 7.3 of the indicator note give the contribution
\(\lambda^2+o(1)\) from unequal shifts. The diagonal is at
most \(M_X\mathbb E W=o(1)\), where \(M_X=\sup_{X/3<n\le3X}D_X(n)
\ll L^{-3/10}\). Thus
\[
 \|W-\lambda\|_{L^2}\longrightarrow0.                         \tag{4}
\]
The actual first and second factorial moments of \(Z\) are
the ones proved in Section 7.3, namely
\(\mathbb E Z\to\lambda\) and \(\mathbb E(Z)_2\to\lambda^2\).
[Uniform model pairs and the two Poisson statistics](CONNORS_KEATING_INDICATOR_CORRELATIONS.md#72-the-sieve-scale-and-the-model-tails)

After removing the collisions involving model slots, the three mixed
sums in the proposition are respectively \((Z)_2W\), \(ZW^2\),
and \(W^3\) in expectation, up to \(o(1)\).
For \(bbD_X\), removal of the two collisions with its model slot
costs at most \(2M_X(Z)_2\). For \(bD_XD_X\), all three collisions
cost at most \(3M_X ZW\), and for \(D_XD_XD_X\) at most \(3M_XW^2\).
These bounds tend to zero in expectation. Replacing \(W\) by
\(\lambda\) in the three expressions costs \(o(1)\) by (4) and
Cauchy–Schwarz with the fourth-moment bounds. For example,
the error in \(\mathbb E ZW^2\) is at most
\(\|W-\lambda\|_2\|Z(W+\lambda)\|_2\), and the second norm is
bounded by Hölder. Hence all three mixed limits are \(\lambda^3\).

Expanding the ordered distinct triple in \(b-D_X\) now gives
\(\mathbb E(Z)_3-\lambda^3+o(1)\). Sort the three shifts and
translate the first to zero. The multiplicity for gaps
\(u,v\) is \(H-v\), and the three slots have six orders.
Changing the endpoints in the inner \(n\)-sum costs at most
\(O(H)\) per tuple since \(|\xi_X|\le1\), so the total normalized
error is \(O(H^4/X)\). This proves (1).

For (2), use the identity \(|b-D_X|=b+D_X-2bD_X\).
Replacing each absolute value by \(b+D_X\) costs at most
\(O(M_X\mathbb E[Z(Z+W)^2])=o(1)\) in the ordered triple sum.
The bound follows by telescoping the three factors and using
\(bD_X\le M_X b\). In the resulting positive expansion
the all-\(b\) term is \(\mathbb E(Z)_3\). The seven other
terms each tend to \(\lambda^3\). Thus averaging shifts
does not make the absolute centered triple negligible.
The upper sieve gives a bounded third factorial moment, and
(1) requires a signed estimate at that bounded scale.

For uniform \(t\) in \([1,X]\cap\mathbb Z\), the same conclusions
hold with the inner sums in (1) over \(1\le n\le X\).
One can discard \(t<X/L^{10}\) in all moments through order
four at cost \(O(H^4/L^{10})\). On the retained range the
model-pair CRT argument still applies with the same outer model.
For records at most \(E=\exp(\sqrt L)\) the rough logarithm is
\(L+O(\sqrt L)\), and the CRT endpoint errors are unchanged.
The one-variable centering is summed over bulk dyads
\(Y\gg X/L^{10}\), all inside the allowed \([X^{1/2},X^2]\).
The actual first and second factorial moments on \([1,X]\)
are already covered in Section 7.3 of the indicator note.

## The third slot can be restricted to raw frequencies of order \(1/H\)

Write \(T_X\) for the signed sum in (1) before multiplication by
\(6/N_X\). Choose \(0\le\Phi\le1\) in \(C_c^\infty((1/2,3))\) equal to one on
\([3/4,9/4]\), and put \(f_X(n)=\Phi(n/X)\xi_X(n)\), extended by zero.
Use \(e(t)=\exp(2\pi i t)\) and
\(\widehat f_X(\alpha)=\sum_{n\in\mathbb Z} f_X(n)e(-n\alpha)\).
Define
\[
 P_{u,X}(\alpha)=\sum_{X<n\le2X}\xi_X(n)\xi_X(n+u)e(n\alpha),
 \qquad
 A_{u,H}(\alpha)=\sum_{u<v<H}(H-v)e(v\alpha).
\]
There is an even nonnegative Schwartz function \(\phi\) of integral one
with \(\operatorname{supp}\widehat\phi\subset[-1,1]\). Under the same pair inputs, for \(2<R\le H\) we have
\[
 \left|T_X-
 \int_{-1/2}^{1/2}\widehat f_X(\alpha)\widehat\phi(R\alpha)
       \sum_{u=1}^{H-2}A_{u,H}(\alpha)P_{u,X}(\alpha)\,d\alpha
 \right|
 \ll_\lambda X\bigl(R/H+L^{-\delta}\bigr),
 \qquad \delta=1/16384.                                      \tag{5}
\]

To construct \(\phi\), take the squared modulus of the inverse
Fourier transform of a nonzero even real smooth bump supported in
\([-1/2,1/2]\) and normalize the integral. Poisson summation gives
\(\sum_j K_R(j)=1\) for \(K_R(j)=R^{-1}\phi(j/R)\) and \(R>2\).
The periodized multiplier is the discrete Fourier transform of
\(K_R\). Consequently the integral is the triple sum after
replacing the third slot by \(\sum_j K_R(j)f_X(n+v+j)\).

For the error estimate, set \(w_u(v)=\mathbf1_{v>u}(H-v)_+\).
For \(|j|\le H^2\), \(\Phi((n+v+j)/X)=1\) throughout the original
sum. After shifting \(v\) the change is given by
\(w_u(v-j)-w_u(v)\) against the original centered product. The
sieve proof of (3) for distinct shifts gives
\[
 \sum_{X<n\le2X}|\xi_X(n)\xi_X(n+u)\xi_X(n+v)|
 \ll X L^{-3/2}J_3(|uv(v-u)|)+X^{1/40}                         \tag{6}
\]
uniformly in this range. The same CRT proof applies to negative
offsets and to offsets \(O(H^2)\), since the shifted arguments
stay in the bulk and the root counts only use their differences.

Let \(d(m)=\prod_{p\mid m}(1+6/p)\) for \(m\ge1\).
The divisor expansion used after (3) gives
\(\sum_{m\le B}d(m)^k\ll_k B\) for every fixed \(k\ge1\).
Since \(J_3(|uv(v-u)|)\le d(u)d(|v|)d(|v-u|)\), Hölder gives a
total \(O(B^2)\) over \(0<u\le B,\ |v|\le B,\ v\notin\{0,u\}\).
On the strip \(0<|v-u|\le s\le H\), with \(0<u<H\) and
\(v\ne0\), the total is \(O(Hs)\). Indeed, fix \(v-u=a\),
bound \(\sum_u d(u)d(|u+a|)\ll H\) by Cauchy–Schwarz, and sum
\(d(|a|)\) over the strip.

For \(1\le |j|\le H\), the weight change is at most \(|j|\)
away from the jump at \(v=u\). The additional change is at
most \(H\) on \(|v-u|\le |j|\). The support lies in
\(|v|\le2H\). Summing (6) thus costs \(O_\lambda(X|j|/H)\)
on distinct shifts. For \(H<|j|\le H^2\), use the weight bound
\(H\) and the \(O(B^2)\) estimate with \(B=H+|j|\). This
costs \(O_\lambda(X(1+|j|/H)^2)\). The \(X^{1/40}\) errors
are absorbed by these bounds.

The terms \(v=0,u\) created by smoothing must retain their signs.
Since \(\xi_X^2=b+O(M_X(b+D_X))\), the single-centered pair
comparison in Section 7.1 and the positive bound (3) give
\[
 \left|\sum_{X<n\le2X}\xi_X(n)^2\xi_X(n+u)\right|
 +\left|\sum_{X<n\le2X}\xi_X(n)\xi_X(n+u)^2\right|
 \ll X L^{-1-\delta}J_2(u).
\]
Here \(M_X\ll L^{-3/10}\), and the shift factor \(J(u)\) in
the pair comparison is at most \(J_2(u)\).
The weight at either collision is at most \(H\), and
\(\sum_{u<H}J_2(u)\ll H\), giving \(O_\lambda(X L^{-\delta})\).
Finally, for \(|j|>H^2\) use boundedness and the trivial
\(O(XH^3)\) bound. The Schwartz tail contributes \(O(XH^{-2})\).
Averaging the distinct-shift bounds against \(K_R\) costs
\(O_\lambda(XR/H)\) by its first and second absolute moments.
These estimates prove (5).
[Single-centered logarithmic-shift comparison](CONNORS_KEATING_INDICATOR_CORRELATIONS.md#71-uniformity-in-the-graph-and-frequency-estimates)

For fixed \(A>1\), put \(R=H/A\), and denote the integral
in (5) by \(I_{X,A}\). Its support in \(\alpha\) is
\([-A/H,A/H]\) for large \(X\). The signed target is equivalent to
\[
 \lim_{A\to\infty}\limsup_{X\to\infty}\frac{|I_{X,A}|}{X}=0.    \tag{7}
\]
For any fixed \(B>0\), after choosing the fixed model depth
sufficiently large, one can omit from this integral the arcs
\(|\alpha-a/q|\le L^B/X\) modulo one, \(q\le L^B\).
Indeed, apply (C3) to each smooth record \(r\le E\) using the
periodic factor \(e(-arv/q)\) in the rough variable \(v\) and
the cutoff \(\Phi(t)e(-(\alpha-a/q)Xt)\) at scale \(X\).
Summing \(1/r\) costs \(O(\log y)\), and the record tail is
superlogarithmically small. This gives
\(\widehat f_X(\alpha)\ll_{B,D}XL^{-D}\) on the arcs for any
prescribed fixed \(D\) at sufficient depth. The arcs have combined
length \(O(L^{3B}/X)\) and
\(\sum_u|A_{u,H}P_{u,X}|\ll XH^3\), so taking \(D>3B+2\)
gives \(o(X)\). Equation (7) still requires joint control of the
centered Fourier transform and the additively twisted pair sums
on the surviving frequencies.
[Recordwise centering (C3)](SPARSE_ARITHMETIC_MODEL_COMPARISON.md#1-independently-checkable-hypotheses)

The separate square-mean bounds have the natural scale of the
remaining integral. For fixed \(A\), take
\(h=\lfloor H/(4A)\rfloor\). On \(|\alpha|\le A/H\),
\(|\sum_{j=1}^h e(j\alpha)|\gg h\). Parseval and (3) give
\[
 \int_{-A/H}^{A/H}|\widehat f_X(\alpha)|^2\,d\alpha
  \ll h^{-2}\sum_n\left|\sum_{j=1}^h f_X(n+j)\right|^2
  \ll_{\lambda,A} X/H^2.
\]
In the last step the \(h\) diagonal sums cost \(O(Xh/\sqrt L)\).
The distinct-shift sums cost \(O(Xh^2/L)\) by the bounded mean
of \(J_2\). The sieve applies throughout the support of \(\Phi\)
by the same interval version of (3).

To test a separate \(L^2\) estimate on the pair factor, replace
\(\xi_X(n)\) by independent centered Bernoulli variables of
success probability \(p\asymp 1/H\). For their variance
\(\sigma^2=p(1-p)\), independence gives exactly
\[
 \mathbb E\left|\sum_u A_{u,H}(\alpha)P_{u,X}(\alpha)\right|^2
       =N_X\sigma^4\sum_u|A_{u,H}(\alpha)|^2.
\]
The third-slot energy has the corresponding exact integrand
\(\mathbb E|\widehat f_X(\alpha)|^2=\sigma^2\sum_n|\Phi(n/X)|^2\),
so its integral over a band of length \(\asymp H^{-1}\) is
\(\asymp X/H^2\). Distinct pairs leave at least one singly
occurring centered variable. On \(|\alpha|\le 1/(12H)\) and \(u\le H/2\), the
real part of \(A_{u,H}(\alpha)\) is \(\gg H^2\). The pair
energy integrated even over this narrower band is therefore
\(\asymp XH^2\) in expectation. Omitting the preceding
polylogarithmic major arcs preserves this order. In the same
independent model distinct sorted triples are orthogonal. Summing
the squared weights costs \(O(XH^4)\) before multiplication by
\(\sigma^6\), so the variance of \(T_X\) is \(O(XH)\).
Thus the separate square-mean scales already multiply to
order \(X\) for an independent model in which the normalized
signed triple tends to zero.

## The reduction at every fixed order

Let \(s_X=\sqrt L/K\). Take disjoint bounded half-open intervals
\(I_i=(a_i,b_i]\) for \(1\le i\le m\), with
\(\lambda_i=b_i-a_i>0\). Let \(Z_i\) and \(W_i\) sum \(b(t+j)\) and
\(D_X(t+j)\), respectively, over integer \(j\) with \(j/s_X\in I_i\).
Use the same uniform dyadic \(t\) and the same outer scale in all
slots. For a fixed vector \(\mathbf r=(r_1,\ldots,r_m)\) of
nonnegative integers let
\[
 F_{\mathbf r}(X)=\mathbb E\prod_i(Z_i)_{r_i},\qquad
 C_{\mathbf r}(X)=\mathbb E\sum_{\mathbf j}^{\!*}
                              \prod_{i=1}^m\prod_{\ell=1}^{r_i}
                                      \xi_X(t+j_{i,\ell}).
\]
The star denotes ordered distinct shifts in each interval,
with \(j_{i,\ell}/s_X\in I_i\). Disjointness handles distinctness
across intervals. Set \(F_{\mathbf0}=C_{\mathbf0}=1\).

Under the pair inputs, for every such fixed vector,
\[
 \begin{split}
 C_{\mathbf r}(X)
  &=\sum_{\mathbf q\le\mathbf r}
       \binom{\mathbf r}{\mathbf q}
       (-\boldsymbol\lambda)^{\mathbf r-\mathbf q}
       F_{\mathbf q}(X)+o(1),\\
 F_{\mathbf r}(X)
  &=\sum_{\mathbf q\le\mathbf r}
       \binom{\mathbf r}{\mathbf q}
       \boldsymbol\lambda^{\mathbf r-\mathbf q}
       C_{\mathbf q}(X)+o(1).                                 \tag{8}
 \end{split}
\]
The powers and binomial coefficients are products over the
intervals. In each term of the expansion with \(q_i\) actual
slots in interval \(i\), the sum differs in expectation from
\(\prod_i (Z_i)_{q_i}W_i^{r_i-q_i}\) by \(o(1)\). Any excluded
collision involves a model slot. Its weight is bounded by \(M_X\)
times the product with that slot deleted. Dropping other
distinctness restrictions only enlarges this positive bound.
The summed error is \(O(M_X)\) by (3) and the bounded fixed
ordinary moments of the counts on a containing window.

For each \(i\), (4) gives \(W_i\to\lambda_i\) in \(L^2\).
Translating the starting point by \(O(s_X)\) changes the
uniform sampling endpoints by \(O(s_X)\), and rounding an
endpoint changes the model count by \(O(M_X)\).
Replacing model counts one at a time by \(\lambda_i\) in
\(\prod_i (Z_i)_{q_i} W_i^{r_i-q_i}\) costs \(o(1)\) in
expectation by Cauchy–Schwarz and (3) through order \(2\sum_i r_i\).
This gives the first line of (8). Finite binomial inversion gives
the second.

The signed sums of total orders one and two tend to zero.
Order one follows from the two first moments. For order two,
the single-centered pair estimate applied separately in each slot
bounds the distinct-shift error by \(O(L^{-\delta})\) in expectation
after averaging the factor \(J_2\) over the containing window.
Translating the first shift back to zero introduces only
\(O(s_X^3/X)\) in the normalized sum.

Thus for any fixed collection of intervals the vector of counts
converges to independent Poisson variables of means \(\lambda_i\)
if and only if \(C_{\mathbf r}(X)\to0\) for every nonzero
\(\mathbf r\). For sufficiency, (8) gives all mixed factorial
moments. Bounded higher moments from (3) yield uniform integrability
on every subsequential limit, and the joint factorial moments
identify that limit. In detail, for nonnegative \(z_i\), the
binomial theorem and monotone convergence give
\(\mathbb E\prod_i(1+z_i)^{Z_i'}=\exp(\sum_i\lambda_i z_i)\)
for any such limit \((Z_i')\), identifying its joint distribution.
Necessity follows from the same uniform integrability and the
first line of (8). Taking all collections of disjoint intervals
would give the finite-dimensional count laws of the unit-rate
Poisson process. For one interval, the signed order \(r\ge2\)
sum takes the span-weighted form
\[
 C_r(X)=\frac{r!}{N_X}
   \sum_{0<u_2<\cdots<u_r<H}(H-u_r)
   \sum_{X<n\le2X}\prod_{\ell=1}^r \xi_X(n+u_\ell)+o(1),
 \qquad u_1=0,                                                \tag{9}
\]
by sorting and translating as in (1).

Each prescribed gap tuple in (9) contributes only
\(O(HL^{-r/2})=o(1)\) for fixed \(r\ge2\), by (3).
For the sums-of-two-squares part of Erdős #675, the corresponding
sufficient arithmetic result would be a positive asymptotic for the
actual indicators of the fixed affine forms \(1+(M_N/a)U\),
\(a\in A_N\), for every \(N\), with the other positions excluded
by the fixed progression in (T.6)–(T.7). The local densities for
these affine forms are positive. The averaged window criterion
(8) leaves that fixed-pattern estimate open.
[Fixed norm-pattern reduction and local densities](ERDOS_675_TRANSLATION_PATTERNS.md#5-an-exact-norm-polynomial-reformulation)

## Tests of direct use of the two-endpoint graph

The existing graph uses physical edge positions \(n,n'=n+ud\),
where \(d\) is the edge's dilation product. Carrying a third
raw slot of gap \(v\) adds a factor at \(n+vd\). Equivalently,
on a raw edge with shared product \(d\) and \(n=dw\), the third
raw value is evaluated at \(w+v\). Averaging over \(v\) leaves
the raw weight
\[
 A_{u,H}(w)=\sum_{u<v<H}(H-v)b(w+v).
\]
This is a function of \(n/d\) on the physical state.
The graph reduction needs a mark-independent first function,
and its raw identity uses permutation-symmetric endpoint
lists. The displayed averaged weight still detects which
prime was selected for the shared product.
[Graph-to-arithmetic reduction, restrictions 1 and 3](../../results/pointwise-binary-splitting-semigroups/graph.tex)
· [Earlier fixed-third-endpoint control](ERDOS_675_TRANSLATION_PATTERNS.md#8-a-higher-correlation-upper-bound-with-no-existence-consequence)

Here are two finite checks of the exact raw identity. Each uses one
bank of primes \(1\bmod4\), two endpoint marks, one shared mark,
\(u=1\), and \(H=4\). All raw pairs are represented values, and each
physical endpoint has exactly two bank prime divisors.

| Bank interval | \(w\) | \(d\) | \(n=dw\) | \(n'=d(w+1)\) | \(2b(w+2)+b(w+3)\) |
|---|---:|---:|---:|---:|---:|
| \([13,233]\) | 193 | 233 | 44969 | 45202 | 1 |
| \([13,233]\) | 233 | 193 | 44969 | 45162 | 0 |
| \([17,233]\) | 136 | 233 | 31688 | 31921 | 0 |
| \([17,233]\) | 232 | 137 | 31784 | 31921 | 2 |

For the first bank the second raw endpoints are
\(194=2\cdot97\) and \(234=2\cdot3^2\cdot13\).
The third-slot value \(196=14^2\) differs from \(236=4\cdot59\);
the other two possible third-slot values are \(3\bmod4\).
For the second bank the first raw endpoints are
\(136=2^3\cdot17\) and \(232=2^3\cdot29\).
The third-slot values at gap two are \(138=2\cdot3\cdot23\)
and \(234\), and at gap three they are \(3\bmod4\).
The first two rows keep the source integer fixed while exchanging
its selected mark. The last two do the same at the target.
They test exact label invariance for the direct insertion with
averaged, tapered third slots. Their scope is the exact raw identity for the direct insertion.

There is also a quantitative limitation to a term-by-term
Fourier insertion into the physical pairing. For a third function
\(c\) supported on a physical interval of \(O(Y)\) integers, Fourier
inversion attaches phases \(e(\theta n)e(\theta v d)\).
Since \(e(\theta v d)=e((\theta v/u)(n'-n))\), the second phase is a
unitary conjugation by an integer-position phase. The first phase
is mark independent. The pair estimate can therefore be used
uniformly for each frequency before the arithmetic reconstruction.
Integrating this uniform bound by the triangle inequality costs
\(\|\widehat c\|_1\). The generic bound supplied by Plancherel is
\(O(Y^{1/2})\) for \(|c|\le1\), whereas the graph gives arbitrary
*fixed logarithmic* savings. A term-by-term average over the
polylogarithmically many \(v\)'s has the same loss on the normalized
average. This route would need an additional joint frequency estimate
or signed control after grouping by the third position.

The unresolved estimate can be stated exactly as
\[
 \sum_{1\le u<v<H}(H-v)
      \sum_{X<n\le2X}\xi_X(n)\xi_X(n+u)\xi_X(n+v)=o_\lambda(X).
\]
It is sufficient and necessary for the desired third factorial limit
given (1). It asks for a single averaged cancellation, with neither
a pointwise triple asymptotic nor a uniform relative tuple error.
A proof at this order alone would leave the all-orders Poisson
distribution and the prescribed configurations in Erdős 675 open.
