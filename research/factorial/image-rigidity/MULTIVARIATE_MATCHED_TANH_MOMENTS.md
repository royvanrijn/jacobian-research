# Matched tanh moments at consecutive orders

For a fully shared ordinary exponent profile, write
\(\ell_i=\min(a_i,b_i)\) and \(d_i=|a_i-b_i|\) on its unequal rows.
Match the Gauss factors indexed by \(1\le h\le d_i\) at orders \(n\)
and \(n+1\), for an integer \(n\ge2\). On such a pair put
\[
 x=\frac{d_i}{n\ell_i+h},\qquad
 x'=\frac{d_i}{(n+1)\ell_i+h},\qquad a=\frac{h}{d_i}.
\]
Write \(s=\log(1+(n-1)x)/2,\ s'=\log(1+nx')/2\) and
\(T=\tanh s,\ T'=\tanh s'\).

**Lemma.** For \(R=n^2/(n^2-1)\), set
\[
 \eta=\frac{t}{2(n+1)+(n^2-2)t}
\]
if all lower-order reciprocals satisfy \(x\le t\).
Without a reciprocal cap one may set \(\eta=1/(n^2-2)\).
For every matched factor,
\[
 \frac{T'}{T}
 =R\left(1+\frac{(2a-1)x}{2(n+1)+(n^2-2a)x}\right)
 \in [R(1-\eta),R(1+\eta)].
\]
For every \(p>0\), write \(M_p=\sum_i\epsilon_i T_i^p\),
\(M_p'=\sum_i\epsilon_i(T_i')^p\), and \(V_p=\sum_i T_i^p\),
where the sums are over matched factors and the same signs
\(\epsilon_i\in\{-1,1\}\) are used at both orders. Then
\[
 \left|M_p'-C_pM_p\right|\le E_p V_p,\qquad
 C_p=\frac{R^p}{2}\bigl((1+\eta)^p+(1-\eta)^p\bigr),\qquad
 E_p=\frac{R^p}{2}\bigl((1+\eta)^p-(1-\eta)^p\bigr).
\]

The reciprocal identity \(x'=nx/(n+1-ax)\) gives
\[
 \frac{T'}T=\frac{n^2}{n-1}
       \frac{2+(n-1)x}{2(n+1)+(n^2-2a)x}.
\]
Dividing by \(R\) and subtracting one gives the displayed
fraction. Its absolute value is at most
\(x/[2(n+1)+(n^2-2)x]\), since \(0\le a\le1\).
The last expression increases with \(x\) and has limit
\(1/(n^2-2)\). Raising the positive endpoint bounds to \(p\)
gives \(|(T'/T)^p-C_p|\le E_p\). Multiply by
\(\epsilon_iT_i^p\) and sum using the triangle inequality.

## Bounds within one order

Set \(U=\sum_i s_i\), and put \(P=M_2/2\). For \(q\ge2\) Jensen's
inequality with weights \(T_i/\sum_jT_j\) gives
\[
 V_q\ge
 \frac{\big(\sum_i T_i^2\big)^{q-1}}{\big(\sum_i T_i\big)^{q-2}}
 \ge \frac{(2|P|)^{q-1}}{U^{q-2}} .
\]
The empty product profile has all moments zero and can be treated
directly. Otherwise \(U>0\). Given an upper bound \(T_*\) with
\(T_*^2\le1-1/q\), the upper bound
\[
 V_q\le \frac{T_*^q}{\operatorname{artanh}T_*}\,U
\]
follows by differentiating \(T^q/\operatorname{artanh}T\).
The derivative has the sign of
\(q(1-T^2)\operatorname{artanh}T-T\), which is nonnegative
because \(\operatorname{artanh}T\ge T\). For an integer \(l\ge2\)
and any \(T_*<1\), the residual bound
\[
 \left|\frac{M_{2l}}{2l}
       -\frac{T_*^{2l-2}}{2l}\,P\right|
 \le \frac{T_*^{2l-2}}{4l}\,V_2
\]
follows by centering each \(T_i^{2l-2}\) at \(T_*^{2l-2}/2\).
These inequalities sum the matched Gauss factors with arbitrary
orientations and do not depend on the number of shared rows.

The odd unsigned and even signed moments are the coefficients
of the reflected mean-defect and skew Taylor vectors,
respectively. Their formulas and geometric tails are derived in the
[private linear comparison](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md#grouping-the-shared-logarithmic-corrections).
The [fully shared terminal-nine proof](MULTIVARIATE_TERMINAL_NINE_SHARED_MATCHED_LIGHT.md#matched-skew-parameter)
uses \(p=2\) and \(\eta=1/206\) at \(n=8\), and uses the
within-order bounds for each subsequent vector.

## Exact tanh transport and moment errors at every order

In the \(T\) coordinate the elementary transfer takes the form
\[
 T'=\frac{n^2T}{n^2-1+(1-2a)T},\qquad
 T'-T=\frac{T(1+(2a-1)T)}{n^2-1+(1-2a)T}>0.
\]
In particular every unsigned \(V_p\) strictly increases for a nonempty
profile. The midpoint comparison also has a version with one higher
moment on the error side. If \(T_i\le T_*\) for all factors, then
\[
 \left|M_p'-R^pM_p\right|
       \le \frac{p n^{2p}}{(n^2-1-T_*)^{p+1}}\,V_{p+1}.
 \tag{4}
\]
Write \(T'=RT(1-\xi)^{-1}\) with
\(\xi=(2a-1)T/(n^2-1)\). The mean value bound on
\((1-\xi)^{-p}\) is
\(p|\xi|(1-T_* /(n^2-1))^{-p-1}\) for its deviation from one.
Multiply by \(R^pT^p\), use \(|2a-1|\le1\), and sum.
For fixed \(p\), (4) retains an \(n^{-2}\) factor and gains one
power of \(T_i\) in the error.

For fully shared ordinary rows put \(S=\sum_i d_i^2/(\ell_i+d_i)\).
Since \(\ell_i\ge1\), the row curvature is at least
\((d_i/\ell_i)^2/(1+d_i/\ell_i)\ge d_i/\ell_i-1\).
Each of the \(d_i\) factors of this row has
\(T_{i,h}\le d_i/(2\ell_i+d_i)\) at either order. Thus one
can use \(T_*=(S+1)/(S+3)\) in (4) and
\(V_{p+1}\le T_*^p S\), because
\(\sum_h T_{i,h}\le d_i^2/(2\ell_i+d_i)\) on each row.

Signed skew can reverse across an adjacent pair at arbitrarily
small total curvature. Fix \(n\ge2\), put \(K=4n+3\), choose an
odd integer \(A\ge3\) and let \(B=KA+4\). Use one increasing
unit row with smaller exponent \((A-1)/2\) and \(K^2\)
decreasing unit rows with smaller exponent \((B-1)/2\).
At order \(j\) the respective factors are
\(T_A=(j-1)/(jA+1)\) and \(T_B=(j-1)/(jB+1)\). Hence
\[
 \frac{T_A}{T_B}-K=\frac{4j+1-K}{jA+1},
\]
and \(T_A^2-K^2T_B^2\) is negative at \(j=n\) and positive
at \(j=n+1\). The total curvature tends to zero as \(A\to\infty\).
Any proof using the signed second moment must allow this reversal.

## A failure of independent errors in the near-flat angular comparison

An all-order extension of the terminal-nine grouped model needs to
carry cancellations through the Cayley transform. In fact the
quadratic exponential model with independent coefficient errors admits
false shared roots at every fixed even lower order \(n=2m\ge28\),
even with its matched moments fixed at exact ordinary values.

Take the two mirrored unit-row pairs
\((\ell,\ell+1),(\ell+1,\ell)\) and
\((2\ell,2\ell+1),(2\ell+1,2\ell)\), and let \(\ell\to\infty\).
Denote the normalized polynomials by \(P_n,P_{n+1}\), and let
\(\delta=1-u\) be the lower outer deficit. Each of the two different
lower-order factors \(T\) is \(\Theta(\delta)\). At \(k=m-1,m,m+1\)
the normalized log coefficient relative to its Gaussian reference is
\(\sigma=\tau+E_{\rm true}\). Here \(\tau\) retains the cubic and
quintic mean-defect terms and \(E_{\rm true}\) comprises the higher
ones. In this mirrored profile the signed terms vanish. To check
their signs, set \(b_j=2j/(n-1)-1\) and \(k_0=\min(k,n-k)\). The
log coefficient correction from one mirrored pair is
\[
  2\sum_{j=0}^{k_0-1}
      \left(\operatorname{artanh}(b_jT)-b_j\operatorname{artanh}T\right).
\]
This follows by pairing the normalized Gauss logarithms at \(j\)
and \(n-1-j\). For \(j<k_0\) one has \(-1\le b_j<0\), so
every summand has nonnegative coefficients in its odd-power series.
At all three indices the cubic coefficient is strictly positive.
Consequently \(0<\tau=\Theta(\delta^3)\) and \(E_{\rm true}\ge0\).

With the retained quantities set to their exact values, use the
quadratic exponential model
\[
 e^\sigma \in 1+\tau+\tfrac12\tau^2+[-W,W],\qquad
 W=\frac{L^3}{6(1-L)}+E(1+L),\qquad
 L\ge\tau+E,\quad E\ge E_{\rm true},\quad L<1 .
\]
These are the coefficientwise error intervals from the quadratic
remainder and omitted logarithm estimate. Each such central
interval has room of order at least \(\delta^{12}\) on both sides of
the true coefficient. For the upper side use
\[
 e^\sigma-1-\tau-\tfrac12\tau^2
   \le e^L-1-L-\tfrac12 L^2+E(1+L).
\]
The difference between \(L^3/[6(1-L)]\) and the exponential tail on
the right is at least \(L^4/8\), by comparing the series termwise.
For the lower side \(e^\sigma-1-\tau-\tau^2/2\ge0\), and
\(W\ge L^3/6\). Multiplication by the positive Gaussian reference
coefficients preserves these orders of slack.

By the [linked Hermite limit](MULTIVARIATE_FIXED_TERMINAL_HERMITE_SEPARATION.md#exact-cancellation-in-the-cayley-transform)
the upper polynomial has a root \(t=e^{i\theta}\) with
\(\pi-\theta=\Theta(\sqrt\delta)\) for all sufficiently large \(\ell\).
Indeed a simple nonzero imaginary Hermite root gives a unique nearby
Cayley root, and reciprocal conjugation of the mirrored profile
fixes that root. In the same bounded Cayley scale the lower polynomial is
\(O(\delta^{m})\). Its reciprocity makes \(G=(-t)^{-m}P_n(t)\) real.
If \(G\le0\), add \(-G(-z)^m\) to \(P_n(z)\). If \(G>0\), add
\[
       \frac{G}{2+2\cos\theta}(-z)^{m-1}(1+z)^2.
\]
The first increment evaluated at \(t\), after division by \((-t)^m\),
is \(-G\). The second has the same property because
\((1+t)^2/(-t)=-(2+2\cos\theta)\).
Both increments are nonnegative for negative real \(z\), and their
coefficient modifications have size \(O(\delta^{m-1})=o(\delta^{12})\).
They fit the three independent remainder intervals. The modifications
leave the outer coefficients and all retained model parameters fixed.
They also preserve coefficient positivity and binomial-normalized
strict log convexity, since the true logarithmic second differences
there are positive of order \(\delta\).
The independent coefficient mean-defect caps allow this perturbation
as well. For each of these three indices, the mean-defect series for one
factor divided by \(\operatorname{artanh} T\) is \(c_k T^2+O(T^4)\)
with \(c_k>0\). The actual defect per unit outer loss is a weighted
average of these ratios at the two tanh values. The ratio of those
values tends to two and the two weights tend to \(2/3,1/3\).
Thus the actual defect is separated from the cap that uses the largest
tanh value by \(\Theta(\delta^3)\). The perturbation is smaller than
this margin. The modified pair has the exact common root \(t\).
The construction concerns the enlarged error model and leaves the
exact ordinary angular question open.

## What a paired angular estimate must control

An antisymmetric logarithmic perturbation of a reciprocal reference
has zero first angular derivative at every simple unit-circle root.
Indeed, let \(P_\epsilon(z)=\sum_k b_k e^{\epsilon h_k}z^k\)
with \(b_{n-k}=b_k\) and \(h_{n-k}=-h_k\), all real. For real
\(\epsilon\), the nearby branch from a simple unit root satisfies
\(z(-\epsilon)=1/\overline{z(\epsilon)}\), by the coefficient
reflection and uniqueness of the branch. Its argument is therefore
even in \(\epsilon\). This puts the first signed-skew contribution
to the angle at quadratic order. The symmetric mean defect also
contributes to the angle. A uniform argument needs adjacent-order
bounds on the two angular changes, including the transformed
remainders, compared with the gaps of the linked reference pair.

Grouping by degree before taking absolute bounds preserves the
near-flat cancellation. With \(X=\sum x\) and \(\lambda=\sqrt X\),
the order-\(X^r\) log coefficient is a polynomial in the index
of degree at most \(r+1\) by the Gauss power sums. After exponentiation
the order-\(X^p\) coefficient has degree at most \(2p\). In the
order-\(j\) Cayley transform a factor of \(\lambda^h\) adds degree at
most \(h\). The alternating binomial sum kills the entire group
when \(2p+h<j\). At \(2p+h=j\) only the highest index degrees
survive, giving the linked Hermite anchors from the quadratic log
term. The moment bounds can be applied to these groups while
preserving both facts.

The [paired angular expansion](MULTIVARIATE_MATCHED_SECOND_ORDER_ANGLES.md)
makes this comparison through the second and fourth orders in the
near-flat Cayley scale, uniformly in the degree for the respective
truncated models. The skew correction starts quadratically and the
leading fourth-order coefficient depends on the order through \(jR_j\).
The exact angular remainder at the linked Hermite gap scale remains open.
