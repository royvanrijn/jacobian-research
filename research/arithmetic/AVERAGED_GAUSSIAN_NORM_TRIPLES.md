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

## Independent Bernoulli control for the joint tests

Let \(\eta_n=B_n-p\), where \(B_n\) are independent Bernoulli
variables with \(p\asymp H^{-1}\), and \(\sigma^2=p(1-p)\).
Use the same finite interval of length \(N_X\), Fourier cutoff
\(\Phi\), and quadratic factor \(\sum_u A_{u,H}P_{u,X}\) as in
(5), replacing every \(\xi_X\) by \(\eta\). The expectation
of their Fourier product is zero at each frequency. For every
indexed term one has a product of three centered variables with
a singleton index, even when the third index coincides with
one of the first two.

The entire smoothed joint integral in (5) also has variance
\(O(XH)\), uniformly for \(2<R\le H\). In fact the assertion holds
if the third slot is convolved with any nonnegative kernel of
mass one. The new coefficient of the third offset \(v\) is
\(\widetilde w_u(v)=(w_u*K_R)(v)\), with
\(\|\widetilde w_u\|_\infty\le H\) and
\(\sum_v\widetilde w_u(v)^2\le \sum_v w_u(v)^2\ll H^3\).
The cutoff \(\Phi\) only decreases the bound. For distinct
indices, group by unordered triples. At most six terms yield
the same triple, and products from different triples are
orthogonal. This portion has variance \(O(\sigma^6 XH^4)=O(XH)\).
For the two collision types use
\(\eta_n^2=(1-2p)\eta_n+\sigma^2\). The quadratic coefficients
are \(O(H)\) at \(O(XH)\) pairs, and the linear coefficients
are \(O(\sigma^2 H^2)\) at \(O(X+H)\) sites. Their variances
are \(O(\sigma^4 XH^3)\) and
\(O(\sigma^6 XH^4)\), respectively. Thus the joint integral
divided by \(X\) tends to zero in \(L^2\), at the same
parameters for which the separate energies in (5) admit an
order-\(X\) Cauchy–Schwarz bound.

For the bounded Palm test, each \(\eta_{t+j}\) is independent
of \(Z-B_{t+j}\). The expectation of the Palm summand is zero
for every \(0\le z\le1\), with \(0^0=1\). The entire window
summand has bounded second moment since its modulus is at
most \(Z+Hp\). Disjoint windows are independent, so averaging
over \(N_X\) starting sites has variance \(O(H/X)\), uniformly
in \(z\).

## The clean anchored graph has a parity obstruction

Fix \(0<u<v\) below the smallest bank prime. Write \(a\) for
the physical third position of an edge of dilation \(d\).
The endpoints are \(a-vd\) and \(a-(v-u)d\). If \(S(a)\) is
the sum of signed weights of the edges with third position
\(a\), then squaring after this grouping pairs
\[
 a-vd_1,\quad a-(v-u)d_1,\quad a-vd_2,\quad a-(v-u)d_2 .
\]
For the following test the endpoint functions are constant and
the archimedean phase is zero. The ideal law samples independent
bank labels with their reciprocal-prime weights and samples the
residues of \(a\) independently and uniformly at the bank
primes. A clean pattern has its selected primes distinct in
each bank, counting a shared label once. Require each endpoint
to have exactly its marked bank-prime divisors.

Put \(s_g=J-f_g\) for the number of shared labels in bank \(g\),
where \(f_g\) probes have been freed. A shared prime forces
\(a=0\bmod p\). An endpoint-only prime forces one prescribed
nonzero residue for \(a\), since it does not divide \(d\).
A freed auxiliary prime forces \(a\ne0\bmod p\), since it
divides \(d\) and is absent from both endpoint lists.
An unused prime leaves \(p-2\) residues for \(a\), exactly
one of which is zero. These four local rules give the exact
parity expectation in a bank, conditional on distinct labels
and the endpoint constraints,
\[
 \mathbb E\bigl[(-1)^{\omega_g(a)}
              \mid\text{labels and endpoint constraints}\bigr]
   =(-1)^{s_g}\prod_{p\text{ unused in }g}\left(1-\frac{2}{p-2}\right).
                                                               \tag{10}
\]
In particular the auxiliary prime supplies no replacement
for the parity of the removed shared prime.

For the large-bank signed family of the binary proof, let
\(\mathcal B_1,\mathcal B_2\) be its two blocks, and let
\(P_1=m|\mathcal B_1|\) be the probe count in the first
block. Test the normalized signed anchor density \(S(a)\)
against \(\chi(a)=(-1)^{\sum_{g\in\mathcal B_1}\omega_g(a)}\).
Normalize by \(c_{\rm loc}=\prod_{g,p\in\mathcal P_g}(1-2/p)\).
To see the normalization, a clean pattern in bank \(g\)
has \(c_g=M+t_g\) distinct marked endpoint primes,
with \(M=J+1\) and \(t_g=1+f_g\).
The exact-count probability equals \(c_{\rm loc}\) times
the product of \(1/p\) at these marked primes, up to a
relative \(O((\log L)/p_{\min})\) for fixed \(J,m\).
Indeed, removing each marked prime from the avoidance product
costs a factor \((1-2/p)^{-1}\), and an auxiliary prime
changes \(1-2/p\) to \(1-1/p\). The weight \(V_g^{-c_g}\)
converts the product of \(1/p\) into the independent label law.
Collisions have probability \(O((\log L)^2/p_{\min})\) under
that law. The total supremum of the signed coefficients is a
fixed power of \(L\), so these errors are superlogarithmically
small.

By (10) every pattern carries, up to the same accuracy,
the factor \(C(-1)^{|A|}\) when its first-block freed set
is \(A\), where
\[
 C=(-1)^{J|\mathcal B_1|}
     \prod_{g\in\mathcal B_1}\prod_{p\in\mathcal P_g}
           \left(1-\frac{2}{p-2}\right).
\]
At zero phase the comparison kernels reproduce constants
to arbitrary logarithmic precision in mean. This follows
from the principal-character case of the frame reproduction
estimate and the bound for the discarded log cells. The
ordinary raw term contributes \(C\), while the comparison
over nonempty \(A,B\) contributes \(-C(-1)^{|B|}\) to
the parity pairing. Summing in \(B\) and then in \(A\) gives
\[
 \mathbb E_a[S(a)\chi(a)]
      = C\,2^{P_1}+O_D(L^{-D}).                              \tag{11}
\]
Here the frame accuracy is chosen to absorb the subset count.
Since \(|C|=\exp(-2\sum_{g\in\mathcal B_1}V_g+o(1))\), this
pairing is bounded below for \(m\ge6\) when \(V_g\le2\).
In fact its modulus grows if \(|\mathcal B_1|\to\infty\).
By Cauchy–Schwarz \(\mathbb E_a|S(a)|^2\ge
|\mathbb E_a[S(a)\chi(a)]|^2\). Thus the original signed
family has no small squared norm in this independent-residue
anchor geometry, even with constant endpoints.
[Binary frame reproduction and signed subsets](../../results/pointwise-binary-splitting-semigroups/graph.tex)
· [Walsh, Section 7.2](https://arxiv.org/html/2610.40439v2#S7.SS2)

There is an exact finite example with four distinct endpoint
positions. Use \(u=1,v=2,J=m=1\) and one bank in each block,
\[
 (13,17,29,37,41,53),\qquad (61,73,89,97,101,109).
\]
Take \(d_1=13\cdot61=793\), \(d_2=17\cdot73=1241\), and
\(a=33704162568656802976\). In each listed bank the residues
of \(a\) are \(0,0,2d_1,d_1,2d_2,d_2\), taken modulo the
respective primes. The endpoint prime sets are
\[
\begin{array}{c|l}
 a-2d_1 &13,29,61,89\\
 a-d_1  &13,37,61,97\\
 a-2d_2 &17,41,73,101\\
 a-d_2  &17,53,73,109 .
\end{array}
\]
Both raw edges therefore contribute at the same anchor.
A clean comparison frees the sole shared probe in each bank.
It uses four endpoint-only primes and one auxiliary prime
per bank, none of which can divide \(a\). At most one bank
prime then remains available to divide \(a\). The displayed
anchor is divisible by two in each bank, so there are no clean
comparison edges at that anchor. The squared sum has a strictly
positive cross term from the two distinct raw edges.

## Cancellation against the local anchored graph weights

The parity obstruction belongs to the residue part of the
anchored operator. The centered arithmetic indicator cancels
against the entire clean local residue weight, including this
parity test, as follows.

**Lemma (recordwise centering for moduli supported on small primes).**
Assume ZF. Let \(Q\le\exp(L^{1/2})\) have all prime divisors
at most \(y=\exp(L^{2/5})\). Let \(X^{1/2}\le Y\le X^2\),
\(r\le E=\exp(\sqrt L)\) be a \(y\)-smooth record, and
let \(g_r\) be any complex function of modulus at most one
on residue classes modulo \(Q\). Write \(m_r\) for the
mean of \(|g_r|\) on the reduced classes. For smooth cutoffs
\(\Psi\) with the support and derivative budgets of (C3),
and sufficiently deep fixed model depending on the desired
\(D\),
\[
 \left|\sum_{v\ {\rm rough}}\xi_X(rv)\Psi(rv/Y)g_r(v)\right|
 \ll \frac Yr\bigl(m_r L^{-D}+\exp(-cL^{3/5})\bigr).           \tag{12}
\]
Here rough means that all prime factors exceed \(y\).

Expand \(g_r\) into characters modulo \(\operatorname{lcm}(4,Q)\).
Each coefficient has modulus at most \(m_r\), and their total
absolute mass is at most \(\sqrt{4Q}\). The principal and
\(\chi_4\) coefficients are handled by (C3) at period four.
For each other character the two \(L\)-factors of the rough
indicator series, \(L(s,\eta)\) and \(L(s,\eta\chi_4)\),
are holomorphic and zero-free in the strip
\(\Re s\ge1-1/\log y\) for large \(X\) by ZF. The Euler
correction is holomorphic there, and the removed Euler
factors cost \(L^{O(1)}\), since \(p^{1/\log y}\le e\)
for \(p\le y\). With \(\sigma=1-1/\log y\), the bound
\[
 |L(\sigma+it,\eta)|\ll (\log y)\bigl(Q(2+|t|)\bigr)^{1-\sigma}
\]
follows by summing up to \(4Q(2+|t|)\) and then using the
periodic partial-sum bound \(4Q\). The same bound applies
to \(\eta\chi_4\) with adjusted constants. Mellin inversion
can consequently use the full vertical line \(\Re s=\sigma\),
with the cutoff transform absorbing the height by fixed
derivatives. The saving \((Y/r)^{-1/\log y}\) gives
\(\exp(-cL^{3/5})\) even after the factor \(\sqrt{4Q}\).
For the model term use the rough Dirichlet series of
\(1\) with twists \(\eta,\eta\chi_4\). Put its smooth
density \(M_J(v)/V\) in the Mellin cutoff. Its fixed rescaled
derivatives cost fixed powers of \(L\) by (CK.6).
This proves (12).
[Model and recordwise centering](CONNORS_KEATING_INDICATOR_CORRELATIONS.md#2-the-indicator-density-and-both-principal-branches)
· [Exact ZF hypothesis](PRIMITIVE_NORM_BINARY_IMPORTED_INPUTS.md#zf-a-fixed-all-conductor-zero-free-half-plane)

For the graph application, let \(G(a)\) be formed from any
fixed-length clean patterns above, with their exact endpoint
bank constraints and the binary list normalizations. Auxiliary
draws use the bank laws. Allow signed label coefficients whose
total supremum over patterns is \(L^{O(1)}\).
For every fixed \(D\), under ZF at sufficient model depth,
\[
 \sum_a \xi_X(a)\Psi(a/Y)G(a)\ll_D Y L^{-D}.                   \tag{13}
\]
Dividing by \(c_{\rm loc}\) is allowed by increasing the
requested precision. One can also include \(\chi(a)\) in (13).

For fixed labels let \(P\) be the product of the forced
endpoint primes. There is one prescribed residue of \(a\)
at each \(p\mid P\). At unused and auxiliary primes there
are at most two forbidden residues. Use inclusion–exclusion
cut off before degree \(k=\lceil C_1\log L\rceil\).
For a parity test in \(\mathcal B_1\), at a
nonforced prime \(p\) write the *combined* local factor as
\(1+h_p(a)\), where \(h_p\) is supported on at most three
residues and has modulus at most two. The product is of
modulus at most one. If \(b\) such factors have \(h_p(a)\ne0\),
the truncation before degree \(k\) is exact for \(b<k\);
for \(b\ge k\) the error is at most
\(1+\sum_{j<k}2^j\binom bj\le 3^k\binom bk\).
The CRT main bound for the unweighted truncation error of
each forced list is therefore at most
\[
 \ll (Y/P)\,3^k\,\frac{(\sum_{g,p}3/p)^k}{k!}.
\]
This also covers the test without parity. Parity at forced
primes is a constant for each list. The endpoint normalization
\(V_g^{-c_g}\) turns \(1/P\) into independent draw weights,
so the error has arbitrary logarithmic savings for large
fixed \(C_1\). In all CRT counts the endpoint \(O(1)\) error
can be summed trivially over all selections and labels.
Their number is \(\exp(O(L^{2/5}\log L))=X^{o(1)}\).

Each retained term prescribes residues at a squarefree
modulus \(Q\) with \(O(\log L)\) bank primes, so (12) applies
on every record. Its principal term contains \(m_r\).
Normalize the weights \(1/r\) over all smooth records by
\(Z_y=\prod_{p\le y}(1-1/p)^{-1}\ll\log y\). Under this law
\(p\mid r\) independently with probability \(1/p\).
If \(p\mid r\) the tested residue of \(rv\bmod p\) is zero.
If \(p\nmid r\) it is uniform among the units when calculating
\(m_r\). Thus
\[
 \sum_{r\le E,\ r\ {\rm smooth}}\frac{m_r}{r}
 \le \frac{Z_y}{P}\prod_{p\ {\rm selected}}\frac{3}{p}.
\]
Summing the retained terms with their coefficients costs at
most \(\prod_{g,p}(1+6/p)=L^{O(1)}\) after the \(1/P\)
normalization. The nonprincipal error in (12) pays even for
the trivial number of labels and selections.

To discard \(a_{\mathcal S}>E\) before inclusion–exclusion,
choose a smooth divisor \(r'\mid a\) with \(E<r'\le Ey\)
by multiplying prime factors one at a time. The CRT bound
with the forced residues is \(O(Y/\operatorname{lcm}(r',P)+1)\).
Put \(\eta=1/\log y\). The main sum over such divisors is
bounded above by
\[
 \frac YP E^{-\eta}
   \prod_{p\le y}(1-p^{-1+\eta})^{-1}
   \prod_{p\mid P}\bigl(1+(p-1)p^{-1+\eta}\bigr).
\]
The Euler factors cost a fixed power of \(L\), since
\(p^\eta\le e\) and \(\omega(P)=O(\log L)\). The factor
\(E^{-\eta}\) is superlogarithmically small. The \(O(1)\)
terms total \(X^{o(1)}\) even before using the label weights.
The proof of (13) is uniform under any common translate of
the anchor residue conditions. It uses the cardinalities of the
allowed sets at the forced and excluded primes, and the
large-record CRT bound also permits translated residues. In
particular take \(G_0(a)=G(a)\Phi_0(a/Y)\) for a bounded smooth
compactly supported cutoff, and \(f(a)=\xi_X(a)\Psi(a/Y)\).
For \(K_R\) as in (5), first apply (13) to each translate and
then average with \(\sum_j K_R(j)=1\). Parseval gives the joint
Fourier bound
\[
 \left|\int_{-1/2}^{1/2}\widehat f(\alpha)\widehat\phi(R\alpha)
                     \widehat G_0(-\alpha)\,d\alpha\right|
        \ll_D YL^{-D} \qquad(2<R\le H).                      \tag{13a}
\]
After translating the anchor, the product of the two cutoffs
has uniform support and derivative bounds at scale \(Y\).
The bound holds for every allowed gap in the mean-spacing
window. It controls the local part of the quadratic window
statistic. The centered arithmetic values at the two
endpoints still vary along the admissible residue classes
and are the remaining source of joint correlations.

## Bounded Palm tests

For \(R_X(z)=\mathbb E\sum_j \xi_X(t+j)z^{Z-b(t+j)}\)
and \(G_X(z)=\mathbb E z^Z\), the model concentration gives
\[
 R_X(z)=G_X'(z)-\lambda G_X(z)+o(1)\qquad (0<z<1).             \tag{14}
\]
Indeed the sum with the \(b(t+j)\) factors is \(G_X'(z)\).
Replacing the \(D_X\)-weighted test by \(Wz^Z\) costs
\(O_z(M_X\mathbb EZ)\), and \(W\to\lambda\) in \(L^2\).
The counts are tight by their bounded first moments.
Along a weak limit, both \(z^Z\) and \(Zz^{Z-1}\) are
bounded tests at each fixed \(0<z<1\). Vanishing of \(R_X(z)\)
throughout that interval therefore forces
\(G'(z)=\lambda G(z)\) and \(G(1-)=1\), identifying the
Poisson law. This is the generating-function specialization
of the reduced Palm identity.
[Chen and Xia, Stein's method and Palm theory](https://doi.org/10.1214/009117904000000027)

One genuinely bounded part of the Palm test cancels for all
\(z\). For an integer \(n\), let \(r(n)\) be its full
\(y\)-smooth record. For any modulus \(Q\) as in (12)
and arbitrary functions \(g_j\) of modulus at most one,
\[
 \left|\mathbb E\sum_{j=1}^H \xi_X(t+j)
           g_j\bigl(r(t+j),(t+j)\bmod Q\bigr)\right|
       \ll_D L^{-D}.                                         \tag{15}
\]
Apply (12) with \(v\mapsto g_j(r,rv\bmod Q)\) and sum over
records, using \(\sum_{r\ {\rm smooth}}1/r\ll\log y\).
The large-record tail follows from the argument for (13)
with \(P=1\). Smooth the interval endpoints with transition
width \(XL^{-B}\) for a sufficiently large fixed \(B\),
using boundedness on the transitions. The factor \(H\)
is absorbed by choosing the accuracy before the depth.
In particular take \(g_j\) to be the conditional mean of
\(z^{Z-b(t+j)}\) given this anchor record and residue.
Equation (15) cancels that part of the actual bounded Palm
test, uniformly in \(z\). It also allows any bounded
window test computed from residues modulo \(Q\) alone.
For example the window can use valuation-one exclusions
at all primes \(p\le c\sqrt L\), with a sufficiently small
fixed \(c>0\), since its modulus is a product of the
squares of those primes.

For the full Palm discrepancy the pair estimates show
uniformly for \(0\le z\le1\) that
\[
 |R_X(z)|\le o(1)+O_\lambda((1-z)^2).                          \tag{16}
\]
Use \(z^j=1-j(1-z)+O(j(j-1)(1-z)^2)\). The constant
term is centered by the first moments. The linear term
is the sum of single-centered pairs, controlled by the
logarithmic-shift comparison. The positive remainder
costs \(O((1-z)^2)\) by the mixed triple upper sieve (3).
In particular the unresolved error at fixed \(z\) already
has the correct second-order onset at \(z=1\).

Exact pair independence itself admits a discrepancy of this
order. In a window of \(H\) exchangeable Bernoulli sites,
choose \(Z=0,1,2\) with respective probabilities
\((H-1)/(2H),1/H,(H-1)/(2H)\) and then choose the occupied
sites uniformly given \(Z\). Every site has mean \(1/H\)
and every pair has joint occupancy probability \(1/H^2\).
With \(D_j=1/H\) the Palm discrepancy tends to
\(z-(1+z^2)/2=-(1-z)^2/2\). This example isolates the
dependence beyond pairs that must be controlled in (14).
The conditional fluctuation of the Palm test left after
(15) has yet to be estimated for the arithmetic indicator.

## Anchor-preserving exchanges

Let \(m\) be the common anchor. The dilations \(D_1,D_2\)
are squarefree products of bank primes with exactly \(J\)
prime factors per bank. Write \(D_i=g d_i\), where
\(g=(D_1,D_2)\), and write \(m=g d_1d_2 r\).
For endpoint offsets \(a_i,b_i\) in \([1,H]\), stripping
the dilation leaves
\[
 (d_2r-a_1,\ d_2r-b_1),\qquad
 (d_1r-a_2,\ d_1r-b_2).                                      \tag{17}
\]
An exchange only reallocates the exclusive primes in the same
bank and leaves the common primes in both dilations. It
preserves the union of the dilation primes and hence their
least common multiple.

For any bank prime \(p\mid m\), the bound \(p>H\) gives
\(p\mid m-cD_i\) if and only if \(p\mid D_i\), for every
\(1\le c\le H\). So the exchange preserves the number of
forced prime divisors per bank at each endpoint.
For \(D_1\ne D_2\), the equality \(aD_1=bD_2\) with positive
\(a,b\le H\) would give \(d_1\mid b\) and \(d_2\mid a\).
Each nontrivial \(d_i\) has a prime divisor exceeding \(H\),
so such an equality is impossible. The argument includes
\(d_i=1\) on either side individually.

The exchange \(61\leftrightarrow73\) in the finite example
following (11) keeps the same anchor. Its complete endpoint
bank-prime records are
\[
 \begin{array}{c|l|l}
  D& m-2D &m-D\\\hline
  949&13,73&13,73\\
  1037&17,61&17,61 .
 \end{array}
\]
Thus it drops the extra prime in both banks at all four
endpoints. Requiring exactly \(J+1=2\) bank primes at each
endpoint rejects both exchanged dilations. Exchanging only
in the first bank produces the same two dilations in the
opposite order. Exchanging in both banks merely interchanges
the two original edges. In this example the unordered
admissible pair has no genuine exchange.

### The ideal on a fixed overlap profile

In bank \(b\) put \(k_b=J-|D_1\cap D_2|_b\), where the
intersection means the common prime set. Fix that common set
and the union of the \(2k_b\) exclusive primes. The allocations
are \(k_b\)-subsets of those \(2k_b\) primes, assigned to
the first edge, with the complement assigned to the second.
The product of the reciprocal-prime label weights is constant
on this fibre. The forced-anchor CRT weight is also constant
there, since each union prime is imposed once. Thus the
label law, before imposing the extra endpoint marks, is the
product of uniform laws on these finite sets. This is the
Bernoulli–Laplace exchange walk on each bank.
[Diaconis and Shahshahani, Time to Reach Stationarity in the Bernoulli–Laplace Diffusion Model](https://doi.org/10.1137/0518016)

Let \(S_b\) exchange one uniform prime on each side and take
\(P_b=(I+S_b)/2\) for \(k_b>0\). For \(k_b=0\) take the
identity on the singleton. Put \(P=\bigotimes_b P_b\) and
let \(\Pi\) average over the fibre. For a \(j\)-subset \(B\)
of the exclusive primes write \(F_B(S)=\mathbf1_{B\subset S}\).
Counting the cases \(B\subset S\), \(|B\setminus S|=1\),
and \(|B\setminus S|\ge2\) gives
\[
 S_b F_B
 =\left(1-\frac{j(2k_b-j+1)}{k_b^2}\right)F_B
   +\frac{k_b-j+1}{k_b^2}\sum_{q\in B}F_{B\setminus\{q\}}.
\]
The spaces generated by the indicators of degree at most \(j\)
form a filtration of the entire function space. Self-adjointness
therefore gives the displayed coefficient as the eigenvalue
on each new orthogonal degree. Every nonconstant lazy
eigenvalue lies in \([0,1-1/k_b]\). It follows that
\[
          \|P^t-\Pi\|\le (1-1/J)^t.                         \tag{18}
\]
If all \(k_b=0\), the left-hand side is zero. Taking
\(t\gg_A J\log L\) gives arbitrary logarithmic savings for
this difference of two orbit averages.

The raw-plus-exchanges operator has a different norm. On any
nonempty fibre, give the raw identity coefficient one and use
any signed comparisons supported on genuinely changed numerical
allocations. At zero phase the resulting matrix has diagonal
entries all one and hence operator norm at least one.
The same conclusion holds on restriction to any nonempty set
of admissible allocations. Repeated swaps that return to the
original allocation can change this diagonal, but a small
norm would require their diagonal coefficient to be
\(-1+o(1)\) on every state. Those returns leave all four
numerical endpoint positions unchanged. Estimate (18) gives
mixing after subtracting the orbit average. It does not
give the small norm of a raw-plus-genuine-exchanges family.

### The extra endpoint marks isolate local exchanges

The extra marks can be analyzed in the independent residue
model, conditional on the anchor-dividing prime sets.
Let \(s\asymp\log L\) be the number of banks, with
\(v_-\le V_b=\sum_{p\in b}1/p\le v_+\). Condition on
sets \(A_b\) of bank primes dividing \(m\), containing the
fixed union of the dilation primes. Assume the total size
of these sets is bounded by a power of \(L\), as it is for
physical anchors \(m\le X^{O(1)}\).
At \(p\notin A_b\) the residue of \(m\) is uniform over the
nonzero residues, independently across primes.

Test any fixed collection of distinct numerical offsets
\(cD_i\), with \(1\le c\le H\). Outside the primes dividing
one of their nonzero differences, their residues at a
non-anchor prime are distinct and nonzero. The indicators
for extra marks at that prime have a mutually exclusive
law with marginal \(1/(p-1)\). The total variation distance
to independent Bernoulli indicators of parameter \(1/p\)
is \(O_\ell(p^{-2})\) for \(\ell\) tested offsets. At exceptional
primes the distance is \(O_\ell(p^{-1})\). Their number is
\(O_\ell(\log(2HD_{\max})/\log p_{\min})\), where \(D_{\max}\)
is the largest tested dilation and \(p_{\min}\) is the
smallest bank prime. Summing the distances and the
\(O_\ell(\sum_{p\in\bigcup_b A_b}1/p)\) cost of omitting the
anchor primes gives a superlogarithmically small error.
Coupling each remaining Bernoulli variable with a Poisson
variable of mean \(1/p\) costs another \(O_\ell(\sum p^{-2})\).
The extra mark counts therefore have, to this accuracy,
independent Poisson laws of means \(V_b\) over banks and
distinct endpoint positions.

Put \(\pi_b=V_b e^{-V_b}\). When \(D_1\ne D_2\), the
probability of four admissible endpoints is
\(\prod_b\pi_b^4(1+o(1))\). A genuine exchange that changes the
unordered dilation pair gives four additional distinct
endpoints. Indeed \(D_1D_2\) is fixed, so sharing a dilation
would force the unordered pair to agree. Thus the
conditional probability that the new configuration is
admissible, given the original one, is also
\(\prod_b\pi_b^4(1+o(1))\). The absolute error in the coupling
is superlogarithmically small, while these probabilities
are bounded below by fixed inverse powers of \(L\), so
the relative assertions are uniform. Since \(\pi_b\le1/e\),
a union bound over \(O_J(\log L)\) one-prime exchanges shows
that with probability \(1-o(1)\) an admissible state has
no admissible genuine one-prime exchange. The same holds
for allocations reachable by any fixed number of swaps
before testing admissibility. This estimate is for the
independent residue model, including the endpoint masks.
It does not exclude a family with many simultaneous swaps
and different signed weights.

For identical dilations all \(k_b=0\). If the two ordered
offset pairs agree, the extra-mark test then involves only
two distinct endpoints. More generally it involves two,
three or four according to the size of the union of the
offset pairs. The corresponding factor is
\(\prod_b\pi_b^{\,\ell}\) with that value of \(\ell\), rather
than the four-endpoint factor in the distinct-dilation sector.

One can also quantify overlap under the normalized ideal
law that counts anchored pairs. Before extra marks, the
mass of profile \(k\) in bank \(b\) is
\[
 w_b(k)=\frac{V_b^{J+k}}{(J-k)!(k!)^2}(1+o(1)),
              \qquad 0\le k\le J.                           \tag{19}
\]
The exact mass is the elementary symmetric sum of degree
\(J+k\) in the \(1/p\), times
\(\binom{J+k}{J-k}\binom{2k}{k}\) for the common and ordered
exclusive sets. Its approximation excludes repeated
draws of probability \(O_J(1/p_{\min})\).
For fixed \(k_0\ge1\), the successive mass ratios
\(V_b(J-k)/(k+1)^2\) give probability \(O_{k_0}(J^{-1})\)
of \(k<k_0\), uniformly over banks for large \(J\).
Independence across banks then bounds the probability
that at least half are in this range by
\(2^s(O_{k_0}(J^{-1}))^{s/2}\). The same estimate with
an additional fixed exponential factor in \(s\) applies
after imposing admissibility. Indeed, conditional on the
union being forced, at unused primes the residues are
uniform and the preceding Poisson argument gives the
four-endpoint factor for distinct dilations. The diagonal
may have the larger two-endpoint factor, a ratio of at
most \((\inf_b\pi_b)^{-2s}\). Taking the fixed \(J\) large
gives arbitrary logarithmic relative savings for the
sector with too few exclusive primes in most banks,
including the identical-dilation sector.

This last estimate concerns relative mass. The unnormalized
phase-zero square can have a large diagonal. For instance,
with equal gap pairs and the binary raw ordered-list weight
\(V_b^{-(J+2)}\), each unordered \(J\)-set has \(J!\) orders
and the two extra labels are unique when the endpoints have
exactly \(J+1\) marks. The identical-dilation contribution
to the local residue mean of the square is therefore, up
to \(1+o(1)\),
\[
 \prod_b (J!)^2 V_b^{-2(J+2)}
           \frac{V_b^J}{J!}\pi_b^2
       = \prod_b J!\,V_b^{-J-2}e^{-2V_b}.
\]
This is the model with constant endpoint functions and no
archimedean localization. It prevents using relative overlap
rarity as an absolute root-moment bound in that model.

## One actual endpoint with a logarithmic-period mask

Under ZF and G, fix \(A,B>0\), let \(Q\le L^B\), and
let \(g_0\) be a bounded \(Q\)-periodic function of modulus
at most one. For \(1\le |h|\le L^A\), the existing
single-centered comparisons extend to this mask. In
particular, uniformly on dyadic subintervals \(I\),
\[
 \sum_{n\in I}\xi_X(n)\xi_X(n+h)g_0(n)
 \ll_{A,B}\frac{XJ_{\rm CK}(h)}{L^{1+\delta}},
                                                               \tag{20}
\]
where \(J_{\rm CK}(h)=\prod_{p\mid h,\ p\equiv3(4)}(p-1)/(p-2)\)
and one may use the exponent \(\delta\) of Section 7.1.
The model depth is fixed after \(A,B\).

Transport the mask along the edge product. For the residue
class \(a\bmod Q\), insert the physical edge condition
\(\mathbf1_{n\equiv aD\bmod Q}\), where \(D\) is that edge's
full dilation product and \(n\) is its first physical
position. The following expansion has coefficient mass one.
\[
 Q^{-1}\sum_{t\bmod Q}e(tn/Q)e(-taD/Q).
\]
On an edge \(n'-n=hD\), each second phase is a unitary
integer-position conjugation and the first phase preserves
mark independence. Thus the same physical comparison holds
with arbitrary fixed logarithmic accuracy. All bank primes
exceed \(Q\) for large \(X\). When \(D=D_R D_C\),
\(n=D_R w\), the condition becomes \(w\equiv aD_C\bmod Q\).
The raw member has \(D_C=1\), giving the requested mask.

For the other members put \(q=Q/d\), \(d=\gcd(a,Q)\).
Separating the condition on \(w,D_C\) gives
\[
 \mathbf1_{w\equiv aD_C\bmod Q} =
 \frac{\mathbf1_{d\mid w}}{\varphi(q)}
  \sum_{\psi\bmod q}
       \psi(w/d)\overline{\psi(a/d)}\,\overline{\psi(D_C)} .
\]
Use the trivial character of modulus one when \(q=1\).
The character of \(D_C\) splits across the two auxiliary
products. The factor
\(\mathbf1_{d\mid w}\psi(w/d)\) in the centered endpoint is
a bounded periodic test of period \(Q\). It leaves the
minor-arc multiplier estimate intact and adds at most
\(L^B\) to the periodic budget in the wide-arc energy
argument. In small-bank prime extraction this factor
contributes \(\psi(p)\) at a prime \(p>Q\), which is
constant after splitting the prime into its residue classes.
The same centered record means and high-frequency
prime bounds apply with the enlarged polynomial-logarithmic
moduli. Thus each separated comparison still has arbitrary
fixed logarithmic savings. The multiplicative expansion
has coefficient mass one, and summing over residue classes
for \(g_0\) costs at most \(Q\).

Apply the well-marked reconstruction to the original bounded
mask before expanding it into residues. The discarded bank
records and endpoint strips are bounded by the unmasked
positive envelope and incur no \(Q\)-loss. In the retained
part the \(Q\)-loss and all periodicity costs are absorbed
by the earlier arbitrary savings. This gives the
single-centered bound with \(b\) or \(D_X\) in the other
slot. Taking their difference gives (20).
[Physical comparison](../../results/pointwise-binary-splitting-semigroups/graph.tex)
· [Uniform shift energy and reconstruction](CONNORS_KEATING_INDICATOR_CORRELATIONS.md#71-uniformity-in-the-graph-and-frequency-estimates)

Take fixed smooth cutoffs \(\Phi_0,\Phi_1\) supported in
the interior of the dyadic interval, and \(|h|=O(H)\).
Set \(f(n)=\Phi_0(n/X)\xi_X(n)\) and
\(F_h(a)=\Phi_1(a/X)g_0(a)\xi_X(a-h)\). For the
Fourier multiplier of (5) with \(2<R\le H\), one obtains
\[
\begin{split}
 \int_{-1/2}^{1/2}\widehat f(\alpha)\widehat\phi(R\alpha)
                         \widehat F_h(-\alpha)\,d\alpha
 &=K_R(h)\sum_a \Phi_0((a-h)/X)\Phi_1(a/X)g_0(a)
                                      \xi_X(a-h)^2\\
 &\quad+O_{B,\lambda}(XL^{-1-\delta/2}).                \tag{21}
\end{split}
\]
Write the smoothed first slot as
\(\sum_j K_R(j)f(a+j)\) and use (20) on the
noncolliding translates \(|j|\le H^2\), with the translated
cutoffs absorbed by partial summation. For these shifts
\(J_{\rm CK}(h+j)\ll_\epsilon L^\epsilon\) for every fixed
\(\epsilon>0\), from its divisor product. The remaining
Schwartz tail is absorbed by boundedness. The coincident
translate \(j=-h\) gives the displayed term, since \(K_R\)
is even. In the independent sparse Bernoulli model the
same diagonal has expectation \(\sigma^2 K_R(h)\) times
the sum of its smooth periodic weights, and the off-diagonal
expectation is exactly zero. The diagonal cannot be included
in the cancellation claim.

For a single dilation \(D\) made of split bank primes,
one can apply (21) after passing from \(m=Dw\) to \(w\).
The invariance in (CK.8) turns \(\xi_X(m)\) and
\(\xi_X(m-cD)\) into \(\xi_X(w)\) and \(\xi_X(w-c)\).
The periodic mask here is in the quotient \(w\),
and the original physical scale is \(XD\). For two
distinct dilations, stripping \(D_1\) leaves the additional
condition \(d_2\mid w\) at the common anchor, even before
accounting for the extra endpoint masks. This divisor is
\(y\)-smooth and its divisibility condition is compatible
with one-variable recordwise centering. However, the
existing pair estimate pays the full dyadic scale \(X\)
and supplies no relative error on the progression of
density \(1/d_2\). Physical transport of one residue class still has coefficient
mass one when the modulus is coprime to the bank labels.
The wide-arc estimates and the relative reconstruction in
(20) only give the stated logarithmic modulus range. Extending (13a) to one actual endpoint
with those constraints requires a relative two-point
comparison under the anchor conditions.

## The second-order Palm error

Let \(\Delta_X(z)=G_X'(z)-\lambda G_X(z)\). The bounded
fourth moment from (3), Taylor's formula at \(z=1\),
and the known first two factorial limits give
\[
 \Delta_X(1-\epsilon)=
       \frac{\epsilon^2}{2}\left(\mathbb E(Z)_3-\lambda^3\right)
          +O_\lambda(\epsilon^3)+o_X(1)
     \qquad (0<\epsilon<1).                                  \tag{22}
\]
The derivative remainder uses \(\mathbb E(Z)_4\), and
the remainder for \(G_X\) uses \(\mathbb E(Z)_3\).
Thus
\(\lim_{\epsilon\downarrow0}\limsup_{X\to\infty}
 |\Delta_X(1-\epsilon)|/\epsilon^2=0\)
would already imply the centered triple by (1). The same
criterion can use the Palm fluctuation after subtracting
the conditional mean in (15), since its other terms vanish
in \(X\) at every fixed \(\epsilon\).

## Relative estimates on the moving progression

Let \(Q\le\exp(L^{9/20})\) be a product of split primes at
most \(y=\exp(L^{2/5})\), with repetitions allowed, and put
\(Y=X/Q\). Suppose \(0<|h|\le L^A\) and \((h,Q)=1\).
This range includes the exclusive dilation products in
(17), whose logarithms are \(O_J(L^{2/5}\log L)\).
Use the same outer model at scale \(X\) throughout.

On \(X<n\le2X,\ Q\mid n\), the model product satisfies
\[
 \sum_{Q\mid n}D_X(n)D_X(n+h)
    =\frac{Y T_h}{L}\bigl(1+O_A(L^{-1/2})\bigr).              \tag{23}
\]
Under ZF, with a sufficiently deep fixed model, the mixed
centered terms satisfy, for every prescribed \(D>0\),
\[
 \left|\sum_{Q\mid n}\xi_X(n)D_X(n+h)\right|
 +\left|\sum_{Q\mid n}D_X(n)\xi_X(n+h)\right|
       \ll_{A,D} Y L^{-D}.                                    \tag{24}
\]
Thus at these moduli the only main-scale replacement still
needed for the two-centered product is the two-actual
asymptotic on the progression.

For (23), use (CK.14)–(CK.16) and the uniform valuation
sieve of Section 7.2. All the moduli in its inert and dyadic
tests are coprime to \(Q\). The substitution \(n=Qu\)
is a bijection modulo each such test modulus. It replaces
the interval length by \(Y\), keeps the same Haar factor,
and leaves the \(O(1)\) count error per progression unchanged.
For the discarded smooth-record tail at \(n\), the full
record equals \(Q\) times the full record of \(u\), so the
omitted record at \(u\) exceeds \(E/Q\). Here
\(\log(E/Q)=(1+o(1))\sqrt L\). For \(n+h=Qu+h\),
all its divisors are coprime to \(Q\). If its smooth
record exceeds \(E\), it has a smooth divisor in \((E,Ey]\),
and each such divisor prescribes a single residue of \(u\).
Rankin's argument on those counts bounds both omitted tails
by \(Y L^{O(1)}\exp(-cL^{1/10})+X^{o(1)}\).
The same sieve truncation and Euler tails as in Section 7.2
give (23).
[Model product and the uniform valuation sieve](CONNORS_KEATING_INDICATOR_CORRELATIONS.md#72-the-sieve-scale-and-the-model-tails)

For (24), the proof of (12) also works for moduli
\(R\le\exp(C\sqrt L)\) with any fixed \(C\), provided
their prime divisors are at most \(y\). The change costs
at most \(\exp(O_C(\sqrt L))\) in the character sum and
\(\exp(O_C(L^{1/10}))\) in the bound for the \(L\)-factors.
Both are absorbed by the contour saving
\(\exp(-cL^{3/5})\). The principal term retains the
factor \(m_r\). This factor measures the proportion of
permitted units after the record is fixed.

First write the summands as
\(\xi_X(u)D_X(Qu+h)\) using split-prime invariance. For the
model endpoint fix its full smooth record \(s\le E\).
Its rough variable is \(v=(Qu+h)/s\).
The factor \(1+\chi_4(v)\) on rough \(v\) is
\(2\mathbf1_{v\equiv1\bmod4}\). Separate this residue
condition, the divisibility by \(s\), and the smooth density
\(M_J(v)/V\). The last factor can be included in the
rescaled cutoff at \(Y\). The terms in (CK.6) have
successive relative sizes \(O_J(L^{-3j/5})\), and fixed
rescaled differentiation preserves the logarithmic bound,
with its exponent independent of the depth. The previous
tail argument justifies
working only in this range.

To impose roughness, expand the exclusion of
\(p\mid v\), \(p\le y\), up to degree
\(k=\lceil C_1\log L\rceil\) by Bonferroni, taking \(C_1\)
large after the desired precision. For a squarefree
selection \(d\) the added condition is \(sd\mid Qu+h\).
It is empty if a prime of \(sd\) divides \(Q\), and
otherwise it gives one progression of \(u\) modulo \(sd\).
Combining with \(v\equiv1\bmod4\) gives a compatible
progression of modulus \(R\) between \(sd\) and \(4sd\).
For the Bonferroni error use the unweighted \(k\)-fold
divisibility count, bounded by \(Y/(sd)+O(1)\) per term.
Summing uses
\[
 \sum_{s\ {\rm smooth}}\frac1s\ll\log y,\qquad
 \sum_{d\ {\rm squarefree},\,\omega(d)=k,\,p\mid d\Rightarrow p\le y}
          \frac1d \le \frac{(\sum_{p\le y}1/p)^k}{k!}.
\]
The main error is arbitrarily small on the scale \(Y\),
by increasing \(C_1\). The \(O(1)\) errors have
\(\exp(O(\sqrt L)+O(L^{2/5}\log L))=X^{o(1)}\) total count.
The same estimates discard the centered endpoint's
records above \(E\) before applying the truncated sum.

For each retained progression \(a\bmod R\) apply (12)
on the centered endpoint's smooth record \(r\).
Summing the principal bound over \(r\) gives
\[
 \sum_{r\le E,\ r\ {\rm smooth}}\frac{m_r}{r}
       \le \frac1R\prod_{p\le y}(1-1/p)^{-1}.
\]
This also holds for prime powers in \(R\). Under the
normalized weights \(1/r\), the smooth valuations are
independent geometric variables. At each prime dividing
\(R\), multiplying a variable of that valuation by a
uniform unit modulo the specified prime power gives
uniform Haar residues after averaging the valuation.
The progression therefore has density \(1/R\) in this
calculation. Summing over \(s,d\) now pays only powers of
\(L\) in the principal error. The nonprincipal error
\(\exp(-cL^{3/5})\) pays even for the trivial number of
records and terms. This proves the forward bound in (24).

In the reverse orientation use \(D_X(u)\xi_X(Qu+h)\).
The model record \(s\) and rough divisor \(d\) now
give a progression of modulus \(Qsd\), up to the factor
of four, in the centered variable \(Qu+h\). Apply (12)
at scale \(X\) in that variable. The principal computation
has the extra factor \(1/Q\), thus gives \(Y\) again.
The centered record tail is relative to \(Y\) by the
same single-residue divisor count in \(u\), since a
divisor of \(Qu+h\) is coprime to \(Q\). All remaining
bounds are unchanged. To use sharp endpoints, smooth
the interval in \(u\) over transitions of length
\(YL^{-C_2}\) and use boundedness there, with \(C_2\)
chosen after the desired error. This proves (24).
[Recordwise centering with smooth moduli](#cancellation-against-the-local-anchored-graph-weights)

There is also a relative two-actual estimate throughout the
logarithmic range. For every fixed \(B>0\), under ZF and G,
uniformly over the same data with \(Q\le L^B\),
\[
 \sum_{\substack{X<n\le2X\\Q\mid n}}
      b(n)b(n+h)
     = \frac{Y T_h}{L}\bigl(1+O_{A,B}(L^{-\delta})\bigr),
       \qquad \delta=1/16384.                                \tag{25}
\]
In (20), the retained bank reconstruction and arithmetic
comparison already save arbitrary fixed powers of \(L\).
Request enough precision there to pay the factor \(Q\).
For the rare bank records and interval transitions, apply
the positive envelope on the progression itself. All bank
primes are coprime to \(Q\) in this range, and all inert
primes are coprime to \(Q\) because it is split. For a
retained bank record the CRT moduli for the exact split
counts and the inert sieve are therefore unchanged in
density under \(n=Qu\). The progression length in the
envelope is \(Y\) and the count error per term remains
\(O(1)\). If a bank record exceeds \(E\), the divisor
argument above still applies, including to the record of
\(Qu\) at the chosen bank primes. Thus the same relative
rare-bank and transition savings hold at scale \(YJ(h)/L\).
The single-centered comparisons then give the relative
version of (20). Adding the model product (23) yields
(25), reducing the exponent if needed to that of the
original logarithmic-shift theorem.
[Positive bank envelope and reconstruction](CONNORS_KEATING_INDICATOR_CORRELATIONS.md#3-an-envelope-at-the-indicator-correlation-scale)

For the independent centered Bernoulli benchmark, a sum
over \(Q\mid n\) of the off-diagonal pair products has mean
zero and variance \((Y+O(1))p^2(1-p)^2\).
Distinct pairs share at most one site and hence have
zero covariance. At \(p\asymp L^{-1/2}\) this has variance
\(O(Y/L)\), so the progression density belongs in the
scale of the pair sum itself.

### Graph localization and the target norm

For a physical dyad of length comparable to \(U\),
the root-moment estimate itself localizes to an
individual progression modulo a bank-coprime
\(Q\le\exp(L^{9/20})\). In the root-replacement lemma
fix a path and its active congruences as usual. Intersect
each CRT cell with the root progression. The normalized
independent residue probabilities are unchanged, and
the count error is now \(O(Q/U)\) in each cell. The
active modulus and the number of Bonferroni terms are
\(\exp(o(L))\). With \(U\) in the bulk range of G the
relative CRT counting error is \(\exp(-cL+o(L))\), and
the independent-residue truncation error is unchanged.
The original ideal path estimate then
gives the same logarithmic saving with root count
\(O(U/Q)\) in place of \(O(U)\).

This root restriction alone leaves a different scale at
the second endpoint. For a mark-independent first vector
on the progression, the two unweighted exact-list norms
are \(O((U/Q)^{1/2})\) and \(O(U^{1/2})\). The block
Hölder argument therefore only gives
\(O(U Q^{-1/2}L^{-E_0})\) with a prescribed fixed
\(E_0\), even using the localized root trace. The
arithmetic target at scale \(U/Q\) needs a second
progression factor or further cancellation.

A direct target mask that accepts every raw edge has a
quantifiable mark cost. Fix the target list and its big-bank
omission choices in the usual phase-zero comparison, with
one omitted prime per small bank, and let \(F\) be the product
of its small-bank marks. Write \(D_{\rm big}\) for the fixed
big shared product and \(W\) for the product of the small
omitted primes. At a root \(n\equiv0\bmod Q\), the target
residue must satisfy
\[
              n'\equiv hD_{\rm big}F/W\pmod Q .
\]
Suppose \(Q\) exceeds every such \(W\), as in the regime
\(Q\ge\exp(L^{1/3})\) with the indicated small banks.
Unique factorization gives \(\mathcal D=(J+1)^{s_{\rm small}}\)
distinct choices of \(W\), and all factors including \(h\)
are units modulo \(Q\). Hence the \(\mathcal D\) target
residues are distinct. They give \(\mathcal D\) equal
main CRT densities when counting marked lists over integer
positions. Keeping the union of these residues as a
target cutoff costs \(\sqrt{\mathcal D}\) in the target
norm at the progression scale. In the existing rare-transfer
proof the \(J\)-dependent edge saving is
\(L^{O(1)}\mathcal D^{-1/2}\), from the same small-bank
omission choices. Thus increasing \(J\) in that bound
does not pay for this target cutoff. The argument concerns
the direct union mask and the existing rare-transfer bound.
It leaves open a comparison that jointly uses the omitted
labels and the target residue.
[Weighted dilation graphs, Lemmas 3.6 and 3.9 and Corollary 3.11](https://raw.githubusercontent.com/openai/math/main/preprints/Weighted-Dilation-Graphs-Smooth-Shifted-Primes-and-Totient-Fibers-September-24-2026/paper.pdf)

## Average absolute cancellation over product moduli

Put \(L=\log X\) and retain the earlier model \(\xi_X=b-D_X\).
For a fixed nonzero \(|h|\le L^A\) write
\[
 E_Q(h;I)=\sum_{\substack{n\in I\\Q\mid n}}
                       \xi_X(n)\xi_X(n+h),\qquad
 F(n)=\sum_{Q\mid n}Q\nu(Q)
\]
for subintervals \(I\subset(X,2X]\). Let \(\nu\) be a
probability law on split moduli \(Q\le\exp(L^{9/20})\)
supported on primes at most \(y\). All the moduli are coprime
to the primes of the fresh comparison banks, and the fresh
banks are available in both ranges with the usual prime
distribution and with \(\asymp\log L\) banks of masses
\(1+o(1)\) in each range. Assume
\[
                   \mathbb E_{\nu\times\nu}(Q_1,Q_2)
                      \ll L^{C_0}                            \tag{26}
\]
for a fixed \(C_0\). These conditions hold for the product
law below when its original and fresh banks are reserved
disjointly.

The independent sparse Bernoulli benchmark supports even
absolute averaging. For \(\eta_n=B_n-p\) with independent
Bernoulli \(B_n\), the progression pair sum has mean zero
and variance \((|I|/Q+O(1))p^2(1-p)^2\). Distinct pairs
share at most one variable and so are orthogonal. Thus
the expected \(\nu\)-average of \(Q\) times its absolute
value is \(O(p\sqrt{XQ_{\max}})\), where \(Q_{\max}\)
bounds the moduli. At \(p\asymp L^{-1/2}\) this is
\(o(X/L)\).

**Proposition (averaged absolute cancellation over split moduli).**
Under ZF and G, at sufficiently deep fixed model depth,
there is a fixed \(\eta>0\) such that
\[
 \sup_{|c_Q|\le1}
       \left|\sum_Q\nu(Q)c_Q Q E_Q(h;I)\right|
  =\sum_Q\nu(Q)Q|E_Q(h;I)|
  \ll_{A,C_0} XJ_{\rm CK}(h)L^{-1-\eta}.                     \tag{27}
\]
In particular this is \(o(X/L)\), uniformly for logarithmic
shifts, since \(J_{\rm CK}(h)=L^{o(1)}\) in that range.
The distribution may vary with \(X,h\), subject to the
stated uniform bounds. The sign \(c_Q\) may depend on the
entire progression sum.

The argument also gives the following version with smooth
records. Let \(r_0(n)\) be the part of the \(y\)-smooth
record outside all fresh comparison banks. For arbitrary
complex functions \(U_1,U_2\) of this record, with
\(|U_1(r_0(n))|\le F(n)\) and \(|U_2(r_0(n))|\le1\),
one has
\[
 \sum_{n\in I}U_1(r_0(n))U_2(r_0(n+h))
                    \xi_X(n)\xi_X(n+h)
   \ll_{A,C_0} XJ_{\rm CK}(h)L^{-1-\eta}                      \tag{28}
\]
in absolute value. This includes \(c_Q\) depending on the
exterior smooth record of the first endpoint when they are
inserted in its sum over divisors.

To prove this, first note the exact counting identity with
an endpoint error
\[
 \sum_{X<n\le2X}F(n)^2
  = X\mathbb E_{\nu\times\nu}(Q_1,Q_2)
      +O\bigl((\mathbb E_\nu Q)^2\bigr).
\]
Indeed, the pair of divisibility conditions has density
\(1/[Q_1,Q_2]\), and \(Q_1Q_2/[Q_1,Q_2]=(Q_1,Q_2)\).
The error is \(X^{o(1)}\) by the size bound on the moduli.
Choose \(T=L^K\) with a fixed \(K\) large after \(C_0\).
Truncating \(U_1(r_0(n))\) on \(F(n)\le T\) costs at most
\(T^{-1}\sum F(n)^2\), which can be made \(O(XL^{-D})\)
for any requested fixed \(D\). The same global truncation
is invariant under multiplication by every fresh bank
prime.

The upper envelope saves at the original correlation
scale even with the untruncated weight. For each fixed
\(Q\), apply the valuation sieve of Sections 3 and 7.2
of the indicator note to \(n=Qu\), with its fresh-bank
record function \(F_0\). The exact record law is unchanged
because the fresh bank primes and the inert primes are
coprime to \(Q\). A record at a fresh endpoint exceeding
\(E=\exp(\sqrt L)\) has an appropriate smooth divisor
between \(E\) and \(Ey\), coprime to \(Q\). Counting the
residue classes in \(u\) and applying Rankin as in the
relative progression proof bounds this tail by
\((X/Q)L^{O(1)}\exp(-cL^{1/10})+X^{o(1)}\).
The remaining Selberg and CRT count errors are
\(O(X^{1/2+o(1)})\), uniformly in \(Q\).
Multiply by \(Q\nu(Q)\) and sum. For \(0\le F_0\le L^{C_1}\)
and \(w_i\in\{b,D_X\}\) the resulting envelope is
\[
 \sum_{n\in I} F(n)F_0(\text{fresh records})
                             w_1(n)w_2(n+h)
 \ll \frac{|I|J_{\rm CK}(h)}L\mathbb E_\mu F_0
       + XL^{O(1)}e^{-cL^{1/10}}+X^{1/2+o(1)}.               \tag{29}
\]
In particular the main bound uses \(\sum_Q\nu(Q)=1\).
It carries no factor \(T\).

Use (29) with the rare-bank indicator and with the small
endpoint transition strips in the reconstruction. The
Chernoff saving for unsuccessful banks is unchanged.
For the retained good-bank contribution only, insert the
bounded endpoint factors
\[
  T^{-1}U_1(r_0(n))\mathbf1_{F(n)\le T}\,\xi_X(n),
       \qquad
  U_2(r_0(n+h))\,\xi_X(n+h).
\]
Their invariance under the fresh bank primes is exact.
Recordwise centering for the first factor follows from
(C3), since its additional multiplier is constant on
each fixed full smooth record. In the Mellin energy
step that selects a small-bank prime \(p\), removing
\(p\) leaves both the record multiplier and its cutoff
unchanged. Thus the graph pairing, the signed
comparisons and the high- and low-frequency energy
bounds for the retained portion apply with arbitrary
fixed logarithmic accuracy. Pay the factor \(T\) by
increasing that accuracy before choosing the model depth.
The rare-bank and transition contributions were already
bounded using (29), without this rescaling. This proves
(28). Taking \(U_1(r_0(n))=\sum_{Q\mid n}c_Q Q\nu(Q)\),
\(U_2=1\), and interchanging the finite sums gives (27).
[Record envelope and bank reconstruction](CONNORS_KEATING_INDICATOR_CORRELATIONS.md#3-an-envelope-at-the-indicator-correlation-scale)
· [Rare banks and precision order](../../results/pointwise-binary-splitting-semigroups/assembly.tex)

An interval depending on the modulus is also permitted
at a smaller saving. Partition the dyad into
\(\lceil L^{\eta/2}\rceil\) intervals of equal length.
Bound the full subinterval contributions by applying
(27) to each grid interval and summing. At the two
boundary intervals use the individual relative positive
envelope just proved, with \(F_0=1\). Thus
\[
 \sum_Q\nu(Q)Q\,\sup_{I\subset(X,2X]}|E_Q(h;I)|
       \ll XJ_{\rm CK}(h)L^{-1-\eta/2}.                       \tag{30}
\]
The same bound permits modulus-dependent smooth weights
whose suprema and total variations are bounded by a fixed
constant, by partial summation. This accommodates differing
archimedean cutoffs in the individual progression sums.

### The product law and the exchange marginal

Suppose there are \(O(\log L)\) original disjoint banks,
and take at most \(m_0\) independent primes in each,
where \(m_0\) is fixed and the reciprocal-prime laws
are \(\mu_b(p)=1/(V_b p)\) with \(V_b\) bounded below.
Repetitions can either be kept with their multiplicities
or excluded by conditioning on distinctness. The counts
per bank may first be chosen in a joint profile, provided
the label sampling conditional on that profile has the
specified independent laws across banks.

For two samples of a single bank, match equal primes
between their two lists, matching each occurrence at most
once. There is a maximal matching that contributes their
gcd, including all repeated values. Therefore the gcd is
at most the sum of the weights over all partial injective
matchings, where a matched pair has weight
\(p\mathbf1_{p=p'}\). The expected weight of one matched
pair is \(\sum_p p\mu_b(p)^2=1/V_b\), and disjoint
matched pairs are independent. The bank contribution to
the gcd moment is bounded by
\[
            \sum_{j=0}^{m_0}
                 \binom{m_0}{j}^2j! V_b^{-j}.
\]
If distinctness is enforced the bound changes by at most
a bounded factor per bank, since the bank primes tend
uniformly to infinity. Taking the product over the
banks proves (26). The estimate is uniform over the
two profiles, so conditioning on correlated bank counts
before drawing the labels is allowed. It also holds for
a marginal dominated by \(L^{C'}\) times such a law, with
an adjusted \(C_0\).

In the independent residue model of (19), fix an overlap
profile and let \(U_b\) be the full forced union of
\(J+k_b\) primes in bank \(b\). Its weight is proportional
to \(\prod_{p\in U_b}1/p\). Once the union is fixed,
the allocations among the common and exclusive sets are
uniform. The marginal density of a specified \(k_b\)-set
\(S\) of second-edge exclusive primes is proportional to
\[
        \prod_{p\in S}\frac1p\;
          e_J\bigl((1/p)_{p\in\mathcal P_b\setminus S}\bigr),
\]
where \(e_J\) is the elementary symmetric polynomial.
The last factor is
\(V_b^J/J!\,(1+O_{J,k_b}(1/p_{\min}))\) uniformly in \(S\).
Conditional on a distinct-dilation profile, the extra
endpoint mark probability is the common factor
\(\prod_b (V_b e^{-V_b})^4\), up to superlogarithmically
small relative error, by the residue estimate following
(18). Thus the ideal conditioned exchange marginal has
the requisite domination by the product of distinct-draw
laws. The diagonal profile has \(Q=1\) and causes no
problem for this normalized marginal.

A stationary exchange walk on the allocations leaves
this uniform partition law unchanged. For a genuine
exchange proposed by a label-only rule, conditioning on
both configurations passing the bank-mark tests weights
the proposal by \(\prod_b(V_b e^{-V_b})^8(1+o(1))\)
in the independent residue model. All eight tested
endpoint positions are distinct, and the error is
uniform by the residue argument after (18).
Consequently, on profiles where the proposal is a
stationary exchange conditioned to change the unordered
pair, this further bank-mask conditioning preserves
the same marginal domination. The probability of a
genuine move under the stationary walk is constant on
the fibre by its allocation symmetry. These statements
concern the bank-mark tests. The actual Gaussian norm
values can impose further conditioning on the law.

There is room for simultaneous bank choices in the usual
ranges. For split primes the reciprocal mass in
\([\exp(L^a),\exp(L^b)]\) is
\(((b-a)/2+o(1))\log L\). Divide the small range at
\(3/16\) and the big range at \(11/30\), in their
exponents. Each half has enough mass for the
canonical \(\lfloor(\log L)/64\rfloor\) banks of
mass \(1+o(1)\). Use one half for the original banks
and the other half for the new comparison. The prime
length and scale separation in Section 7.1 persists.
If the original canonical banks have already been chosen
as disjoint intervals in the full ranges, their complement
still has enough harmonic mass for the same count.
Fresh banks formed greedily there can be unions of
\(O(\log L)\) intervals. The interval prime estimates
then cost an additional logarithmic factor. In the
smooth prime sums, smoothing all those boundaries at
a sufficiently small fixed logarithmic width preserves
the derivative budgets, and the total transition mass
is bounded by the sum of the individual transition
masses. In either construction the fresh bank law
and the uniform comparison estimates remain available.

A bounded signed coefficient on the original labels may
be averaged conditional on \(Q\), yielding a coefficient
of modulus at most one in (27). In particular a restriction
to a label event can be kept with its original mass, even
if it is rare. The maximal bound (30) permits the interval
or a uniformly bounded-variation weight to depend on
all the labels by taking the supremum before averaging.

### Scope in the four-endpoint geometry

The conclusion controls the averaged pair sums in (27),
including varying cutoffs in (30). It uses that a
divisibility condition \(Q\mid n\) is attached to one
of the two centered arguments, or to the other argument
after translating and reversing the shift. For a nonzero
residue \(a\bmod Q\), the condition \(Q\mid n-a\) is
generally not invariant under multiplying \(n\) by a fresh
prime. A sum over such conditions needs a separate
transport or comparison estimate.

Also, after stripping the two dilations in (17), setting
\(n=d_2r\) makes the second pair depend on \(r=n/d_2\).
If this second pair is left as actual arithmetic factors,
its value varies within the first progression sum.
Choosing \(c_Q\) after evaluating a whole progression
sum, as in (27), does not allow this point-dependent
multiplier. Allowing arbitrary point-dependent signs would
fail even in the independent Bernoulli benchmark. Choosing
the sign of \(\eta_n\eta_{n+h}\) at each point gives
expected sum \(4p^2(1-p)^2\sum_{n\in I}F(n)\), of order
\(X/L\) in the full window. Formula (28) covers separable
weights on the exterior smooth records of the two
tested endpoints. Connecting the actual four-endpoint
weights, which see rough cofactors, to the averages
proved here is a remaining arithmetic step.

## External banks at growing root-moment order

The affine residue weight has the same second-moment bound
as the divisor weight. For fixed \(Q_1,Q_2\), its joint
CRT density is \(1/[Q_1,Q_2]\) if the two residues agree
modulo \((Q_1,Q_2)\), and zero otherwise. Taking absolute
coefficients therefore gives
\[
 X^{-1}\sum_{X<n\le2X}|F_c(n)|^2
       \le \mathbb E (Q_1,Q_2)+O(X^{-1}Q_{\max}^2).           \tag{31}
\]
The remaining issue is the placement of the weight on the
physical edge. If \(n'=n+hD\), the physical third position
is \(n+vD\). When \((h,Q)=1\), testing its divisibility by
\(Q\) is equivalent to \(Q\mid (h-v)n+vn'\).

There is an independently marked CRT replacement for
these edge weights at the order used in the root moment.
Fix a path of \(N=O(L^{0.52})\) edges, all fresh-bank
labels and dilation choices, and the external modulus
\(Q_i\) at each edge, where \(Q_i\le\exp(L^{9/20})\).
Require the primes of the external moduli to be disjoint
from the fresh banks. Write \(s_i\) for the offset of the
third position from the integer root \(m\). The external
conditions are
\[
                m\equiv r_{Q_i}-s_i\pmod {Q_i}\quad(1\le i\le N).
\]
Their joint factor under independent external residues is
zero for incompatible congruences and \(1/[Q_1,\ldots,Q_N]\)
otherwise. Add these congruences to the active
congruences in the proof of Lemma 3.6 of G. They are
coprime to all fresh-bank conditions. The product of the
external moduli is at most \(\exp(O(NL^{9/20}))=X^{o(1)}\).
The same holds for the product after adding the fresh active
modulus and every term of the Bonferroni truncation.
The relative CRT counting error is therefore
\(\exp(-cL+o(L))\) on a bulk root dyad, and the
independent-residue truncation error for the fresh banks
is unchanged. This proves the same pathwise relative
replacement for each compatible external label choice.
After expanding a positive external weight, the
relative error may be summed termwise. With signed
coefficients it is an error relative to the corresponding
absolute path sum. This puts the affine weights into the
residue model but provides no signed norm estimate for
that model.
[Weighted dilation graphs, Lemma 3.6](https://raw.githubusercontent.com/openai/math/main/preprints/Weighted-Dilation-Graphs-Smooth-Shifted-Primes-and-Totient-Fibers-September-24-2026/paper.pdf)

For external product laws with at most a fixed number of
labels per bank, this replacement error can be summed over
the signed root moment without losing precision. Indeed
drop the external compatibility indicator for an upper
bound. In one bank with \(k=O(N)\) external draws,
\(\prod_i Q_i\) cancels the \(1/p\) of every occurrence in
its sampling weight. The CRT denominator still pays
\(1/p\) for each distinct prime, since the least common
multiple contains every used prime. Summing equal-label
partitions bounds the bank factor by \(V^{-k}T_k(V)\),
also for products with multiplicities. Here
\(T_k(V)\le(V+k)^k\), as proved below. Distinct-draw
conditioning costs at most
\(\exp(O(N(\log L)/P))\) over \(O(\log L)\) banks.
It follows that the independent external expectation
is at most \(\exp(O(N(\log L)^2))\), uniformly in the
fixed fresh path and even in the residue assignments.
The same conclusion holds with an extra fixed
logarithmic domination of each sampled law.

In the fresh-bank CRT proof the Bonferroni degree may
be increased to \(C N(\log L)^2\). The sum of excluded
prime densities is still \(O(N\log L)\), so the
truncation error then has arbitrarily prescribed
precision \(\exp(-D N(\log L)^2)\), relative to the
fresh independent mass. The modulus sizes and the
number of truncated terms remain \(\exp(o(L))\).
This precision absorbs the displayed external mass
and the \(L^{O(N)}\) absolute fresh-path mass from G.
It gives a replacement error smaller than any
\(L^{-D N}\) also after summing signed paths. The
result is a signed moment in the enlarged residue
space, whose norm is a separate problem.


### Repeated anchors and exceptional external primes

Consider one external bank of reciprocal mass \(V\asymp1\)
and minimum prime \(P\). Group the list of \(K\) offsets
by equality, writing \(s\) for the distinct offsets and
\(k_s\) for their multiplicities. Write \(T_j(V)\) for the
\(j\)-th moment of a Poisson variable of mean \(V\), with
\(T_0(V)=1\), and \(T_j(V;z)\) for the moment after a
shift by \(z\). Let \(\mathcal E\) be the bank primes
dividing a nonzero difference of two offsets. Put
\(V_0=\sum_{p\notin\mathcal E}1/p\) over the bank, and in
the independent residue model put
\[
 Z_s=\sum_{p\in\mathcal E}\mathbf1_{m\equiv-s\bmod p}.
\]
For \(K^2/P=o(1)\) and \(V_0\gg1\), the following
relative formula holds:
\[
 \mathbb E_{\rm CRT}\prod_s W(m+s)^{k_s}
  =\bigl(1+O(K^2/P)\bigr)V^{-K}
             \mathbb E_{\mathcal E}\prod_s T_{k_s}(V_0;Z_s).
                                                               \tag{32}
\]
For good primes, a label can occur in slots belonging to
only one anchor. Partition the good slots at each anchor
into equal-label blocks, then assign distinct primes to
all the blocks over all the anchors. If there are \(b\)
blocks, their reciprocal-prime sum is
\(V_0^b(1+O(b^2/P))\), by the union bound for a collision
among \(b\) independent harmonic draws. Summing the
partitions gives independent Poisson moments. Condition
on the exceptional residues and expand the shifted
powers by the binomial theorem. All terms are
nonnegative and the relative error is uniform in their
degrees, which proves (32).

Here is one criterion for replacing the exceptional part
by the contribution of equality of integer anchors.
Let \(b_{p,t}\) count the *distinct* offsets of residue
\(t\bmod p\) and put \(M=\max_{p\in\mathcal E,t} b_{p,t}\).
The Poisson identity
\(T_j(V_0;z)=V_0^{-z}\mathbb E[(N_0)_z N_0^j]\)
for integer \(z\ge0\), where \(N_0\) is Poisson, gives
\[
 \frac{T_j(V_0;z)}{T_j(V_0)}
    \le \left(\frac{V_0+j+z}{V_0}\right)^z .
\]
Indeed \(T_{j+1}=V_0(T_j+T_j')\le(V_0+j)T_j\) by the
partition polynomial. With \(B\ge1+(K+|\mathcal E|)/V_0\),
independence over exceptional primes consequently gives
\[
 1\le
 \frac{\mathbb E_{\mathcal E}\prod_s T_{k_s}(V_0;Z_s)}
      {\prod_s T_{k_s}(V_0)}
 \le \exp\!\left(
    \sum_{p\in\mathcal E}\frac1p
          \sum_{t\bmod p}(B^{b_{p,t}}-1)\right).              \tag{33}
\]
In the relevant root range the number of offsets and the
logarithms of their nonzero differences are bounded by
fixed powers of \(L\), and \(|\mathcal E|\) is then also
polynomial in \(L\). The excluded prime mass is
\(O(|\mathcal E|/P)\). In particular, if
\(M\log B\le\tfrac12\log P\), the error in replacing
\(V_0\) by \(V\) and dropping the exceptional part is
superlogarithmically small. The resulting moment is
\(V^{-K}\prod_s T_{k_s}(V)\). This includes arbitrary
multiplicities of equal integer anchors. For a fixed
number of all distinct offsets it gives \(1+o(1)\).

A uniform replacement by this last product at growing
order fails without a restriction on congruence
collisions. Fix an integer \(t>v_+\), where \(V\le v_+\),
and choose \(t\) external primes with product \(P_*\).
Take \(K=\lfloor L^{0.52}\rfloor\) distinct offsets
\(s_i=iP_*\). Their logarithmic size is \(O(L^{2/5})\)
for \(t\) fixed in the bank ranges. The probability that
\(P_*\mid m\) is \(P_*^{-1}\), and on this event every
\(W(m+s_i)\ge t/V\). The joint moment is at least
\((t/V)^K/P_*\), which diverges exponentially in
\(L^{0.52}\). The proposed product for distinct
integer offsets is one. This counterexample concerns
uniform pathwise residue independence. A signed
root-moment bound would have to control the paths
with a large congruence collision count in the sum
over dilation labels.

## Separating the exclusive label laws

Fix a shared set \(g\) and the counts of the exclusive
labels in each bank. In the CRT model, before fixing the
union, the reciprocal-prime law on the two exclusive
sets is the product of their separate distinct-draw
laws outside \(g\), conditioned on their being disjoint.
This follows from the factor \(1/(g d_1d_2)\) for a
forced anchor and the constant profile multiplicities.
Since the maximum harmonic atom is \(O(1/P)\), the
probability of a conflict under the product law is
\(O_J((\log L)/P)\). The disjoint conditioning has the
same total variation cost. The four-extra-mark CRT
probability after averaging the quotient residue is
\(\prod_b\pi_b^4(1+O_D(L^{-D}))\) for distinct
dilations, uniformly in the labels. This follows from
the independent residue calculation after (18), or
from its quotient version with roots
\(r\equiv c\overline{d_i}\bmod p\) at unused primes.
Thus conditioning on these bank tests after
averaging the quotient preserves factorization in
total variation to superlogarithmic precision.

Conditioning on the *union* of exclusive primes gives
a different result. One set is then the complement of
the other in each bank. On an allocation fibre of size
\(M_{\rm alloc}\), the product of the two uniform
marginals assigns mass only \(1/M_{\rm alloc}\) to
complementary pairs. Its total variation distance from
the joint law is \(1-1/M_{\rm alloc}\). In particular
a fibre with two exclusive primes in one bank already
has distance \(1/2\). Factorization applies when the
shared set and the overlap profile are fixed and the
union is averaged.

There is also separation of the raw endpoint masks
with the actual endpoint values retained. Take
\(m=g d_1d_2r\) with \(r\) in a common range of length
\(\asymp R\), and dilation products in comparable
dyadic boxes. For the arm with dilation \(D_i=g d_i\),
strip it from the endpoints, giving
\(d_jr-a_i,d_jr-b_i\), \(j\ne i\).
The indicators and their centered models are
unchanged by this strip since all removed primes
are split and at most \(y\). A prime of \(D_j\setminus D_i\) cannot divide the
stripped values because its size exceeds the
offsets. The only dependence on \(D_i\) in the
extra-prime count is the exclusion of primes of
\(D_i\) from the stripped counts. For \(p\mid D_i\),
\(p\nmid d_j\), and the event that \(p\) divides
either stripped value has at most two residues for
\(r\bmod p\). Over the quotient interval, dropping
all these exclusions changes the normalized sum
against bounded endpoint tests by
\[
                    O_J\bigl((\log L)(P^{-1}+R^{-1})\bigr).
                                                               \tag{34}
\]
After this replacement the first mask is a function
of \(d_j,r\) only, and the second is a function of
\(d_i,r\) only. The normalized error remains
superlogarithmically small after bounded or fixed
logarithmic endpoint weights.

For example, the product of the two raw harmonic
weights at the first endpoint of each arm becomes
\[
 \prod_{i=1}^2
  \frac{\Phi_i((m-a_iD_i)/(XD_i))}{m-a_iD_i}
 =
 \frac1{g^2d_1d_2}
       \prod_{i=1}^2
       \frac{\Phi_i((d_jr-a_i)/X)}{d_jr-a_i},
       \qquad j\ne i .
\]
The factor \(1/(d_1d_2)\) supplies the reciprocal
label law exactly. Comparability windows and
dyadic restrictions on the two arms can be
included as separate weights, without conditioning
on a dyadic cell of possibly tiny mass. With the
shared labels and profile fixed, (34) and the
disjointness error give a separated product of
the two dilation averages in this normalized
raw harmonic expression, up to superlogarithmic
error. This separation does not make either
average independent of \(r\), and it does not
supply the sign cancellation of the centered
average.

## The off-diagonal four-form variance

Take distinct \(1\le a,b\le H\) and
\(d,e\le\exp(L^{9/20})\) composed of split bank
primes in \((H,y]\). Suppose the four arguments \(dr-a,dr-b,er-a,er-b\) are
comparable with \(X\) on a range of length \(\asymp R\),
where \(R=X^{1-o(1)}\). Positive cutoffs and bounded
separated masks may be retained. If \(d\ne e\),
all six pair determinants are nonzero. For
example \(bd=ae\) would force the reduced coprime
parts of \(d,e\) to divide \(a,b\), which are
smaller than any of their prime divisors. An
inert prime at which roots can coincide divides
\[
        ab(a-b)(d-e)(bd-ae)(ad-be).
\]
The valuation-one sieve of (3) applies directly
to these four forms. Every slope is a unit at
inert primes, each retained prime-square atom
has the same root count as for four shifts,
and the interval CRT errors have the same
bound. Thus the absolute four-factor average
is at most \(O(L^{-2}J_4(\Delta))+X^{1/40}/R\),
where \(\Delta\) is the displayed product and
\(J_4\) may be restricted to inert primes.

For the independent product law of the exclusive
labels, one can average the variable part of
this sieve factor. Assume at least one prime is
sampled in each list and that the bank laws are
the reciprocal-prime laws on progression
intervals, with the bounded exclusions used
above. Write \(J_0=J_4(ab(a-b))\). For inert
\(p\nmid ab(a-b)\), after fixing \(e\) and all
but one label of \(d\), the three variable
determinants forbid at most three unit residues
of that label modulo \(p\). For every squarefree
\(q\le L^B\) of these primes, Siegel–Walfisz and
partial summation bound the corresponding
divisibility probability by
\[
            \ll 3^{\omega(q)}/\varphi(q)+O_D(L^{-D}),
\]
at arbitrarily fixed precision. This is uniform
under the choice of the fixed labels. The upper
bound applies to dyadic cutoffs left as
unnormalized nonnegative restrictions by
discarding them.

To use this only at logarithmic moduli, omit
\(p>L^2\) from the variable factor at cost
\(1+o(1)\), since \(\log|\Delta|=O(L)\).
The normalized subset expansion of
\(\prod_{p\le L^2,\,p\mid\Delta}(1+8/p)\) has
expected log divisor \(O(\log L)\). Markov's inequality
shows that a sufficiently large fixed \(B\)
captures at least half of its mass in divisors
\(q\le L^B\). The expected truncated expansion
away from \(ab(a-b)\) is bounded by the
convergent product
\(\prod_p(1+24/(p(p-1)))\), plus a negligible
Siegel–Walfisz remainder. Thus the averaged
off-diagonal positive bound is \(O(L^{-2}J_0)\).
For the identical-dilation term use the two-form
bound, since \(|\xi_X|\le1\) in the bulk. Its
contribution is at most
\(O(L^{-1}J_4(a-b)\sum_d\nu(d)^2)+X^{1/40}/R\).
For this diffuse product law \(\max_d\nu(d)\ll_J
1/P\). The same assertion with a logarithmic
factor holds under logarithmic domination of
the label law. Thus this diagonal is negligible.
Consequently the positive bound for the proposed
\(\mathbb E_r|K(r)|^2\) has scale
\(L^{-2}J_0\), with \(J_0=L^{o(1)}\). The same
bound applies on the other arm. A signed fixed
power saving beyond that scale would suffice
in the Cauchy–Schwarz step.

In the independent sparse Bernoulli benchmark the
distinct-dilation terms in this variance have expectation
zero for each \(r\) in the bulk. The four arguments are
distinct. If arguments from different dilations agreed,
then \(|d-e|r\le H\), impossible for this range.
The expected variance is exactly
\(p^2(1-p)^2\sum_d\nu(d)^2\) without cutoffs, for
Bernoulli parameter \(p\asymp L^{-1/2}\). Thus the
proposed fixed logarithmic improvement is consistent
with the benchmark at every comparable dyadic scale.


One signed part of the four-form expansion also
cancels. Under ZF, uniformly for the fixed
dilations and offsets above, a product of one
\(\xi_X\) and three \(D_X\) factors at these
four arguments, with a smooth cutoff in \(r/R\),
has sum \(O_D(RL^{-D})\) for any prescribed
fixed \(D\), with sufficiently deep fixed
model. This estimate even allows \(d=e\)
when the four factors are read with
multiplicity.

For the proof, fix the exact smooth records
\(s_1,s_2,s_3\le E\) of the three model
arguments and include the smooth factors
\(M_J/V\) in the cutoff. The odd rough part must
be \(1\bmod4\). Each model form \(l_i(r)\)
then imposes \(l_i(r)\equiv s_i\bmod4s_i\).
A prime of its slope cannot divide its
nonzero offset, so compatible conditions from
these prefixes give one progression modulo
\(M=\operatorname{lcm}(4s_1,4s_2,4s_3)\).
In the progression variable, each rough-unit
condition at \(p\le y\) is deterministic or
excludes one residue modulo \(p\). Discard a
prefix with a deterministic failure. At most
three residues are excluded at each remaining
prime. Use Bonferroni to degree \(C\log L\).
The error is bounded by \(R/M\) times a
\(k\)-fold sum of products \(3/p\), plus the
subpower CRT error. Choose \(C\) large after
the required precision. The prefix main masses
are summable with a logarithmic cost because
\[
 \sum_{s_1,s_2,s_3\ {\rm smooth}}\frac1{\operatorname{lcm}(s_1,s_2,s_3)}
      =\prod_{p\le y}\sum_{j\ge0}
               \frac{(j+1)^3-j^3}{p^j}
      \ll L^{O(1)}.
\]
Each retained term fixes \(r\) in a residue
class modulo \(Mq\), with \(q\) a product of
\(O(\log L)\) small primes and with at most
\(3^{\omega(q)}\) classes. In the centered
coordinate, whose slope is \(d\) or \(e\),
the required modulus is \(dMq\) or \(eMq\).
It is smooth and at most \(\exp(C'\sqrt L)\),
so the extension of (12) used in (24) applies.
Its principal bound recovers the exact
residue density, hence costs \(R/(Mq)\) per
class at arbitrary logarithmic accuracy.
The exponentially small nonprincipal error
pays the subpower progression count and
the slope. Smooth-record tails are discarded
by the divisor argument of (24), since any
smooth divisor of a form has its slope
invertible. Derivatives of the three density
factors are fixed logarithmic powers as in
(24). This gives the signed bound.

The retained four-actual or four-centered
off-diagonal product does not fit this
recordwise argument. Applying the binary
comparison to one arm would put the second
arm's two actual values inside its endpoint
test. Multiplying by a fresh prime in the
first arm changes the other affine forms,
so the available invariance and recordwise
centering do not cover that test. The
variance saving required beyond the positive
bound remains open.

## Signed affine moments for a bounded number of external labels

Separate the absolute CRT replacement from the signed
moment in the independent residue model. A useful
strengthening of the fresh-bank moment bound is uniformity
under *different* prescribed additive phases on each edge.
For a fixed \(N=2\lceil L^{0.52}\rceil\), multiply the
\(i\)-th oriented edge of a length-\(N\) fresh path by
\(e(\theta_i D_i)\), where \(\theta_i\in\mathbb R\) are
fixed independently of the fresh labels. After replacing
the root by fresh residues and keeping the integer return
condition, the signed path sum still has arbitrarily
prescribed bound \(L^{-E N}\), with the parameter order
of Theorem 3.5 of G.

This follows inside the proof of G. In the clean edge
estimate of Lemma 3.10, condition on the shared small
product \(d'\). The new phase is \(e(\theta_i d'D_B)\),
and just changes the real ideal frequency on that edge.
For a reversed edge it can be absorbed into the frequency
of the adjoint in the same way. The ideal estimate is
uniform in every real frequency. Absolute row and column
bounds, including the rare memory transfers, are
unchanged. The equality-rank argument restoring distinct
births still only loses the signed bound on edges affected
by birth-equality multipliers. All other edges keep their
small norm with the new frequency. Finally the integer
return still has its unit-mass additive Fourier integral.
None of these estimates requires the edge frequencies to
agree.
[Weighted dilation graphs, Sections 3.6–3.8](https://raw.githubusercontent.com/openai/math/main/preprints/Weighted-Dilation-Graphs-Smooth-Shifted-Primes-and-Totient-Fibers-September-24-2026/paper.pdf)

**Proposition (signed affine root moments with bounded external labels).**
Use fresh banks and the signed family of G, with a compatible
phase range as in the binary comparison. On the \(i\)-th
edge insert an external weight of the form
\[
 \sum_{Q_i}c_{i,Q_i}Q_i\nu_i(Q_i)\,
       \mathbf1_{\,m+\lambda_i(Q_i;\mathbf D)\equiv r_{i,Q_i}
                                                          \bmod Q_i}.
                                                               \tag{35}
\]
Here \(m\) is the integer root, \(|c_{i,Q_i}|\le1\),
and, conditional on \(Q_i\), \(\lambda_i\) is an integer
linear combination of the displacements \(D_j\), with
coefficients independent of the fresh labels. For a
physical edge weight each such expression is of course
given by the current edge and position. The external moduli are at most
\(\exp(L^{9/20})\), are coprime to the fresh banks,
and each contains at most a fixed total number of labels
counted with multiplicity, drawn from harmonic prime
banks as in (26). Mixtures of bounded-count profiles and fixed
logarithmic domination are allowed. The laws and
residues can vary by edge. Then the resulting signed
root moment retains the bound \(L^{-E N}\) for every
requested fixed \(E\), after increasing the graph
precision. In particular, it holds for the independently
weighted affine edge condition
\(\mathbf1_{n_i\equiv a_{Q_i} D_i\bmod Q_i}\) and for
a third-position test, with the corresponding conjugate
weights on reversed edges.

To prove this, first apply the affine CRT replacement
above at arbitrary precision. Freeze the \(Q_i\).
Averaging the independent external residues gives the
factor \(1/[Q_1,\ldots,Q_N]\) times their compatibility
indicator. For each prime power in the least common
multiple, choose a reference incidence of its maximal
exponent. Compatibility at every other incidence
requires equality with the reference modulo the prime
power at that incidence. Thus each compatibility condition
has a finite additive Fourier expansion of coefficient
mass one. Its argument is a difference of two affine
expressions in the \(D_j\). After collecting the Fourier
variables, the resulting factor on the fresh path is a
constant phase times \(\prod_j e(\theta_jD_j)\).
It has the signed bound just proved for every choice of
the Fourier variables. No edge loses its small norm from
this Fourier insertion.

The absolute cost of the remaining external factors is
bounded by
\[
   \mathbb E\frac{\prod_i Q_i}{[Q_1,\ldots,Q_N]}
       \le (1+o(1))\prod_b V_b^{-k_b}T_{k_b}(V_b) .
                                                               \tag{36}
\]
Here the count profile is fixed, with \(k_b\) the total
number of draws in bank \(b\) across the \(N\) edges.
The factor \(1+o(1)\) allows conditioning the draws to
be distinct within each modulus. For random profiles
apply the bound before averaging. The proof is the partition bound following
(31), including the negligible cost of distinct-draw
conditioning. With \(\sum_b k_b=O(N)\),
\(T_{k_b}(V_b)\le(V_b+k_b)^{k_b}\) shows that the cost
is \(L^{O(N)}\), which is paid by the arbitrary precision
of the fresh signed bound. More generally the same proof
works with any external weight family for which the
absolute external cost on the left of (36) is
\(L^{O(N)}\) and the CRT sizes in (31) remain valid.

For the many-bank law with bounded labels *per bank*, a
useful part of the moment still has this bound. Partition
the external slots into equal-prime blocks within each
bank and define its rank as the total number of slots
minus the number of blocks. For any fixed \(C\), the
contribution of partitions of rank at most \(CN\) has
arbitrarily strong signed logarithmic savings. Indeed a
rank-\(r\) partition can be specified by \(r\) edges of a
forest on the \(O(N\log L)\) slot addresses. There are
at most \((O(N\log L))^{2r}\) choices. Its harmonic
weight is at most \(\prod_b V_b^{-r_b}\), where \(r_b\)
is the rank in that bank, by dropping distinctness
between blocks. Thus the total cost is \(L^{O_C(N)}\).
Apply the same compatibility phases for each partition.
This argument includes repeated numerical anchors and
highly occupied individual congruence classes if the
total external rank is in this range.

For a common external test of divisibility at the third
position, there is also a positive restriction on the
high-rank sector. In a fixed path, call a block crossing
if its labelled prime must divide the difference of two
distinct numerical anchors. A crossing block has
external harmonic mass at most
\(O(L^{C_2}/P_{\min})\), by counting prime divisors of a
nonzero difference. The total number of slots is
\(O(N\log L)\), so even the unrestricted absolute
partition sum before the crossing penalties is at most
\(\exp(O(N(\log L)^2))\). Consequently paths and
partitions with at least
\(C_3N(\log L)^2/\log P_{\min}\) crossing blocks are
negligible by an absolute bound, with \(C_3\) chosen
sufficiently large. Repeats at identical numerical
anchors have no such penalty. This restriction uses a
common residue at a fixed anchor, and does not apply
to arbitrary \(r_Q\).

## A many-bank counterexample to the untruncated root bound

For the usual product of external bank marks, the
growth of (36) is reflected in the physical signed
operator, rather than only in its upper bound. Work
at zero phase, \(\ell=1\), and fixed nonzero \(h\).
Use disjoint external banks, \(s_{\rm ext}\asymp\log L\)
of them, and let the edge multiplier be the
untruncated third-position weight
\[
 F(n+vD)=\prod_{b\ {\rm external}}
         \frac{1}{V_b}\sum_{p\in b}\mathbf1_{p\mid n+vD}.
\]
This is one harmonic draw in each external bank. It is
covered by the CRT replacement of the preceding section,
including when the root has an unusually large external
record.

Choose \(J\) shared fresh primes in each bank and fix
their product \(D_*\). In each fresh bank choose one
additional source prime and one additional target prime.
Restrict roots \(n\) to the local event \(\mathcal C\)
that the complete fresh-bank prime lists of \(n\) and
\(n+hD_*\) are exactly the chosen lists. The prescribed
divisibilities occupy compatible residues. Unused fresh
primes exclude at most two residues, so the independent
fresh density of \(\mathcal C\) is at least
\(\exp(-O_J(L^{2/5}\log L))\) in the bank ranges of
the indicator argument.

On this event the signed matrix entry from any fixed
compatible ordering at \(n\) to a compatible ordering at
\(n+hD_*\) receives a positive raw contribution. No
comparison can reach that integer endpoint from \(n\).
Its dilation would have to equal \(D_*\) and it would
then use an auxiliary prime dividing both endpoints,
which is forbidden because their mark lists are complete.
The same argument applies to all source permutations.
The endpoint symmetrizations and the row and state
weights in G cost at most \(L^{C_J}\) in this lower bound.
Thus with \(u_n\) the mark-independent root vector
and \(B_F=A_F A_F^*\) for the untruncated weighted
operator, the spectral Jensen inequality gives
\[
  \langle u_n,B_F^{N/2}u_n\rangle
         \ge L^{-C_JN} F(n+vD_*)^N
                \qquad(n\in\mathcal C).                     \tag{37}
\]
Indeed \(\|A_F^*u_n\|^2\) already has the square of a
noncancelling raw contribution, and
\(\|u_n\|^2\) is between fixed powers of \(L\) on
\(\mathcal C\). The argument also applies to the
operator compressed to exact fresh counts.

The positive expansion of \(F(n+vD_*)^N\) has
\(O(N\log L)\) external label slots. By the same
uniform CRT argument as in (31), its average on
\(\mathcal C\) is its independent-residue expectation
times \(1+o(1)\), including the fresh density.
This follows term by term since the external
congruences all have the same shifted root, and
the total modulus, even after the fresh
Bonferroni restrictions, has logarithm \(o(L)\).
In each external bank its moment is
\(V_b^{-N} T_N(V_b)(1+o(1))\), uniformly with
error \(O(N^2/P_{\min})\) from distinct block
labels as in (32). Taking just the Poisson atom
at \(\lfloor N/\log N\rfloor\) gives
\[
 \log(V_b^{-N}T_N(V_b))
        \ge N\log N-N\log\log N-O(N).
\]
The product over \(s_{\rm ext}\asymp\log L\)
has logarithm \(\gg N(\log L)^2\). This dominates
both the fresh event cost
\(O_J(L^{2/5}\log L)\) and the normalization loss
\(O_J(N\log L)\) in (37). In particular the
averaged root moment in (37) is at least
\(\exp(cN(\log L)^2)\) for large \(X\).
The obstruction already uses the common
divisibility residue at the third position.
It neither challenges the vertex-weight argument
with truncation in (27) nor rules out a different
affine operator or a truncation compatible with
its signed graph proof. The four-centered
off-diagonal saving is also separate from this
lower bound.
