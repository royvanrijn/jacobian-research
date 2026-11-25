# A computable cutoff in the bounded shared quartic window

Keep the coefficient phase \(F\), field \(w\), and notation of the
[shared quartic comparison](MULTIVARIATE_SHARED_MOVING_QUARTIC_WINDOW.md).
All unequal rows are fully shared ordinary rows. Let \(q_*\) be the
real critical center \(F'''(q_*)=0\) near the maximum of
\(\rho(q)=q(1-q)G''(q)\), put \(w_*=F'(q_*)\),
\(a_2=F''(q_*)/2,\ g=F''''(q_*)/24\), and
\(\ell=(-4g)^{-1/4}\). For order \(j\), use
\[
 \lambda_j=j^{-1/4},\qquad
 T_j=-2a_2\ell^2/\lambda_j^2,\qquad
 H_j=\ell(w_*-w_j)/\lambda_j^3 .
\]

**Theorem.** Given rational \(0<W_-\le W_+\) and rational \(B_T,B_H\ge1\),
one can compute positive \(\delta_0\) and an integer \(N_0\).
Whenever \(W_-\le W\le W_+\), \(|\max\rho-1|\le\delta_0\),
and \(|T_j|\le B_T,\ |H_j|\le B_H\) for \(j=N-1,N\ge N_0\),
the two moments have no common physical root. The cutoff is
uniform in the profile, including the total-degree difference.

## Effective center and localization on the real lattice

Put \(r_i=d_i/(2m_i)\) with \(m_i=(a_i+b_i)/2\) and \(d_i=b_i-a_i\).
There is an explicit rational \(R<1\) bounding every \(|r_i|\).
If the larger exponent in a row is at least twice its lower
exponent, the absolute gap \(d\) bounds the lower exponent and the
row contributes at least \(2d/3\) to \(W\). Since the lower exponent
is an integer at least one, the larger exponent is at most
\(1+\max(1,3W_+/2)\) times the lower one. For
\(B=1+\max(1,3W_+/2)\) take \(R=(B-1)/(B+1)\).
The row formula for \(\rho''\) therefore gives
\[
       \rho''(q)\le-c,\qquad
       c=\frac{2(1-R^2)W_-}{(1+R)^3}>0 .
\]
Let \(q_m\) maximize \(\rho\) and \(\delta=\rho(q_m)-1\).
If \(|\delta|\le1/2\), then \(\rho(q_m)\ge1/2\) and
\(\rho(q)\le W_+ q(1-q)/(1-R)\) on the real interval. Thus
\(q_m,1-q_m\ge2\eta\) for
\(\eta=\min(1/8,(1-R)/(4W_+))\).

On the complex disk of radius \(2b\) about \(q_m\), where
\(b=\min(\eta/4,(1-R)/16)\), use
\[
 F''(q)=\sum_i\frac{d_i^2/m_i}{1+r_i(2q-1)}
                -\frac1q-\frac1{1-q}.
\]
The denominators in the row sum have modulus at least
\((1-R)/2\) and the endpoint distances are at least \(\eta\).
Cauchy's bound with radius \(b\) gives explicit upper bounds for
all derivatives of \(F''\) on the smaller disk. At \(q_m\) the
identity \(F''=(\rho-1)/(q(1-q))\) gives
\[
 |F'''(q_m)|\le2|\delta|/\eta^2,\qquad
 F''''(q_m)\le-4c+4|\delta|/\eta^3 .
\]
Choose an effective radius \(b_1<b/4\) so the Cauchy bound on
\(|F^{(5)}|\) times \(b_1\) is at most \(c\), then decrease the
allowed \(|\delta|\) so that \(F''''\le-2c\) on
\(q_m\pm b_1\) and \(2|\delta|/\eta^2<c b_1\).
The intermediate value theorem and monotonicity give a unique \(q_*\)
on this interval, with \(|q_*-q_m|\le |\delta|/(c\eta^2)\).
The bounds on \(F''''\) give computable positive lower and upper
bounds for \(-g\), hence also for \(\ell\).

There are effective \(r>0\) and \(K\ge1\) such that on
\(|h|\le2r\) one has the analytic expansion
\[
 F(q_*+h)-w_j(q_*+h)-F(q_*)+w_jq_*
       =(w_*-w_j)h+a_2h^2+gh^4+h^5 E(h),\qquad |E(h)|\le K .
 \tag{1}
\]
For instance first take \(2r<b_1/2\) and bound the tail of the
Taylor series of \(F''\) by its Cauchy bounds, then shrink \(r\)
until the quartic term dominates the remainder by a factor four.
Shrink it further using the Cauchy bound on \(E'\) so the
derivative of the remainder is bounded by \((-g)|h|^3\)
on the real interval. All these choices depend only on \(W_\pm\).
Decrease \(\delta_0\) to ensure \(|q_*-q_m|\le r/4\) and
\(\delta_0\le cr^2/64\).
For real \(|q-q_*|\ge r/2\),
\(\rho(q)-1\le\delta_0-c(q-q_m)^2/2<0\),
so the real phase is concave on each outer interval.
At \(h=\pm r/2\), the two bounds on \(T_j,H_j\) imply
\(a_2=O(\lambda_j^2)\) and \(w_*-w_j=O(\lambda_j^3)\) with
computable constants. Equation (1) gives a strict negative
real-phase gap there and strict opposite derivative signs, for
all \(j\) above a computable threshold. Concavity on the outer
intervals extends the gap through both endpoints \(q=0,1\).

Normalize each exact summand by the constant factorial coefficient.
At the retained grid it is \(j^{-1/2} A_j(q)e^{j\Phi_j(q)}\),
where \(\Phi_j=F-w_jq\) and \(A_j\) is the Stirling amplitude
from the shared quartic comparison. On the entire integer grid
the summand has modulus at most \(C e^{j\Re\Phi_j(q)}\) with a
computable common \(C\). For this bound use the entropy inequality
for the binomial coefficient. The shared rows' remaining
half-log factors are bounded by \(W_+/(1-R)\), since
\(\sum |r_i|\le W_+/2\), and their Stirling remainders have
a uniform sum using \(\sum 1/m_i\le W_+\).
This gives a computable exponential bound for all omitted summands.

## A linked expansion through fourth order

Put \(\lambda=\lambda_j,\ h_j=\lambda^3/\ell\), and rescale
\(q=q_*+\ell\lambda y\). In the retained interval \(|q-q_*|\le r\),
divide the summands by \(j^{-1/2}A_j(q_*)e^{j\Phi_j(q_*)}\)
and multiply the sum by \(h_j\). This nonvanishing normalization
leaves the positions of zeros unchanged. At both retained grid
endpoints \(|q-q_*|\ge r/2\) for \(j\ge2/r\).
The global gap shows that the omitted sum and both trapezoid
endpoint terms have magnitude bounded by \(C_1 j^{3/4}e^{-\gamma j}\)
with computable positive \(C_1,\gamma\). This is \(O(\lambda^5)\)
with a computable constant, using \(e^{-\gamma j}\le2/(\gamma j)^2\).

Write \(A_j=A_0 e^{R_j}\) on the retained interval, using the
relative real Stirling expansion. The derivative bounds from
[Binet's remainder](MULTIVARIATE_EFFECTIVE_BOUNDED_PRIVATE_CUBIC.md#uniform-real-lattice-expansion)
give \(|R_j'|+|R_j''|\le C_2/j\) there. The endpoint arguments of
the binomial and the row Gamma factors are uniformly positive
multiples of \(j\) and \(jm_i\), respectively. In the row sum
one uses \(\sum_i |d_i|/m_i^2\le W_+\) and
\(\sum_i d_i^2/m_i^3\le W_+\). For the binomial terms
\(q,1-q\ge\eta\). It follows that
\[
 |R_j(q)-R_j(q_*)|\le C_2|q-q_*|/j
                         \le C_2\ell\lambda^5 |y|.
\]
The exact normalization at \(q_*\) thus removes the
order-\(\lambda^4\) constant amplitude correction.

In the rescaled exponent, put
\[
       \Psi_0(y)=-y^4/4-T_j y^2/2+H_j y .
\]
For \(0\le e\le\lambda\) and \(|\ell\lambda y|\le r\), the real
part of \(\Psi_0(y)+e\ell^5y^5 E(\ell e y)\) is bounded by
\(-y^4/8+B_T y^2/2+B_H|y|\), by the choice of \(r\).
Cauchy's bounds on \(E\) and on \(A_0(q)/A_0(q_*)\) control their
derivatives on this interval. For \(A_0\), use its explicit
square-root factors and the same bounded sum of relative logarithms
as above, with the log branches positive on the real interval.
Its real value at \(q_*\) has a uniform positive lower bound.

The second derivative of the exact normalized rescaled summand
is bounded by a fixed polynomial in \(|y|\) times the last
quartic envelope. This follows by differentiating its exponent,
its \(A_0\) factor, and \(R_j(q)-R_j(q_*)\) twice using the
displayed bounds. The real trapezoid error is at most
\(h_j^2/8\) times the integral of this derivative bound, hence
\(O(\lambda^6)\) with a computable constant. Replacing
\(e^{R_j(q)-R_j(q_*)}\) by one in the integral costs \(O(\lambda^5)\)
by the preceding difference estimate. Now Taylor expand
\[
 \frac{A_0(q_*+\ell e y)}{A_0(q_*)}
       \exp(e\ell^5 y^5 E(\ell e y))
\]
through degree four in \(e\) at zero. The product and chain rules
and the same quartic majorant bound the fifth-order error by
\(\lambda^5\) times an integrable polynomial times the majorant.
The four correction coefficients are polynomials in \(y\)
whose coefficients depend on the profile but agree at the two orders.
Extending the coefficient integrals to the full real line
costs \(O(\lambda^5)\): on the added tails
\(|y|\ge r/(2\ell_{\max}\lambda)\), use
\(1\le(2\ell_{\max}\lambda |y|/r)^5\).

For the normalized exact moments \({\cal Q}_j\), this gives
\[
 {\cal Q}_j=Z(T_j,H_j)+\sum_{p=1}^4\lambda_j^p C_p(T_j,H_j)
                         +O(\lambda_j^5),                 \tag{3}
\]
where \(Z(T,H)=\int_{\mathbb R}e^{-y^4/4-Ty^2/2+Hy}\,dy\).
The error and bounds for all needed derivatives of \(Z,C_p\)
are computable on any bounded \(T,H\) rectangle. Multiplication
by powers of \(y\) and \(-y^2/2\) controls the derivatives by
the same moment bounds. To compute any such bound, cut at a
sufficiently large integer \(Y\) so that the quadratic and
linear terms beyond \(Y\) sum to at most \(y^4/16\).
Integrate with interval bounds on \([-Y,Y]\) and bound the
tails by polynomial moments of \(e^{-y^4/16}\).

## Compact residuals and termination

Let \(n=N-1,\ L=\log(1+1/n),\ \lambda=n^{-1/4}\), and use the
same physical \(t\) at the two orders. Then
\[
 T_N=e^{L/2}T_n,\qquad H_N=e^{3L/4}H_n+v,\qquad
 \lambda_N=e^{-L/4}\lambda,\qquad v=\ell\,\Delta\,N^{3/4}L\in\mathbb R.
\]
For \(\Delta\ne0\) the integer degree difference gives
\(|v|\ge(\ell_{\min}/2)\lambda\), since \(L\ge1/N\).
The bounded fields also imply \(|v|\le 2B_H+1\) for large enough \(n\).
Put \(T=T_n,\ H=H_n\).

The first model residual uses the continuous divided difference
\[
 B(T,H,v)=\int_0^1 Z_H(T,H+sv)\,ds .
\]
If the exact polynomials both vanish, (3) implies
\(|Z(T,H)|\le C_3\lambda\). For \(\Delta\ne0\), divide the exact
consecutive difference by \(v\). The two model fields differ from
\(T,H+v\) by \(O(L)\), and the differences of the correction terms
are bounded by a constant times \(\lambda^p(|v|+L)\) for each \(p\).
The two remainders cost \(O(\lambda^5/|v|)\).
Thus \(|B(T,H,v)|\le C_4\lambda\) for a computable \(C_4\).

For \(\Delta=0\) divide the difference by \(L/4\).
The model term equals \(2T Z_T(T,H)+3H Z_H(T,H)+O(L)\)
by two differentiations along the field path.
Each correction now differs by \(O(L\lambda^p)\). The
remainders have size \(O(\lambda^5/L)=O(\lambda)\).
This gives \(|A(T,H)|\le C_5\lambda\) for
\(A=2T Z_T+3H Z_H\), again with computable \(C_5\).

Both residual sums \(|Z|+|B|\) and \(|Z|+|A|\) have positive
infima on the corresponding compact rectangles.
The [quartic model theorems](MULTIVARIATE_DIFFUSE_MOVING_CRITICAL_WINDOW.md#the-bounded-moving-critical-field)
give simple, purely imaginary zeros for real \(T\) and
\(3y_k(T)-2T y_k'(T)>0\) at every positive zero height.
For real \(v\ne0\) two zeros of \(Z\) therefore cannot differ
by \(v\), and simplicity treats \(v=0\).
At a model zero \(H_k(T)\), the other residual is
\(A=Z_H(3H_k-2T H_k')\), which is nonzero.
Negative heights follow by conjugation and there is no real zero.

Compute positive lower bounds for both sums by subdividing the
rectangles. Use real quadrature with the quartic tail bounds at
the centers and the uniform derivative bounds across each cell.
For \(B\) use its integral representation for bounding the
variation, and the divided difference at a center with nonzero
\(v\). As quadrature and mesh tolerances tend to zero, the
strict positive infima prove termination. Taking \(n\) large
enough that the bounds \(C_3\lambda+C_4\lambda\) and
\(C_3\lambda+C_5\lambda\) fall below these lower bounds gives the
claimed numerical cutoff algorithm.
