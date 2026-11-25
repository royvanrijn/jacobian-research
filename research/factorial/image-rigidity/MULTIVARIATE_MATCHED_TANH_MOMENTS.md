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
