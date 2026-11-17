# The exact exceptional Laurent correction

The exceptional distinguished coefficient remains invariant after an
explicit quadratic correction. A zero final head-minus-one coefficient
is sufficient but unnecessary for this scalar Laurent statement.
Existence of the shared heads and their transport from actual operators
remains a separate obligation.

## Hypotheses and exact formula

Use the fixed polynomial chain and shared-head hypotheses of
[the conditional Laurent lemma](CONDITIONAL_LAURENT_COEFFICIENT_INVARIANCE.md):
\[
 H_0(X,Y)=Y,\qquad
 H_{j+1}(X,Y)=H_j(X,Y)^p-\epsilon_j^pX^{q_j}.
\]
Here \(I\ge1\), \(p\ge2\) is even, all \(q_j\) are nonnegative integers divisible by 2,
\(q_{j+1}<pq_j\), and \(\epsilon_j\in k^\times\) are fixed over a
characteristic-zero field. We restrict to its exceptional arithmetic
\(d_0=2,l_0=1\), with final gap \(N\), where \(N>1\).

Explicitly put
\[
 D=2N+1,\quad p_N=pD/2,\quad q_{N,j}=Dq_j/2\quad(j<I),
 \quad P_N=t^{p_N},
\]
\[
 L_{N,k}=(p-1)\sum_{j<k}q_{N,j},\quad
 r_N=L_{N,I}-p_N,\quad q_{N,I}=r_N+N+1.
\]
Require the shared heads
\[
 H_j(P_N,Y_N)=\epsilon_jt^{q_{N,j}}+\text{lower}\quad(j<I),
 \quad\operatorname{ord}H_I(P_N,Y_N)=q_{N,I},
\]
and \(pq_{N,I-1}-q_{N,I}=N\).
These are the exceptional order hypotheses; they are not conclusions
about an arbitrary operator construction. Set
\[
 A_N=[t^{q_{N,I}}]H_I(P_N,Y_N),\quad
 B_N=[t^{q_{N,I}-1}]H_I(P_N,Y_N),\quad
 C_N=[t^{r_N}]H_I(P_N,Y_N),\quad
 \tau_N=[t^{-p_N}]Y_N.
\]
In particular \(A_N\ne0\), while \(B_N\) and \(\tau_N\) are unrestricted.

Construct the fixed common head root as follows. Write
\[
 W_{I-1}(z)=\epsilon_{I-1}z^{q_{I-1}/2},\qquad
 W_j(z)=\epsilon_jz^{q_j/2}
 \left(1+\frac{W_{j+1}(z)}
 {\epsilon_j^p z^{pq_j/2}}\right)^{1/p},\quad W=W_0.
\]
Then \(H_I(t^{p_N},W(t^D))=0\). Put
\[
 \lambda_j=p^j\prod_{h<j}\epsilon_h^{p-1},\qquad
 \epsilon=\epsilon_{I-1}.
\]

**Corrected invariance theorem.** Under these hypotheses, for every
\(N>1\),
\[
 \boxed{\quad C_N-\lambda_I\tau_N
       -\frac{p-1}{p\epsilon^p}A_NB_N
       =-\lambda_I[z^{-p/2}]W.\quad}                         \tag{1}
\]
The right side depends only on the fixed chain data. In particular,
under the trace hypothesis \(\tau_N=0\), the combination
\(C_N-(p-1)A_NB_N/(p\epsilon^p)\) is independent of \(N\).
No missing-head assumption \(B_N=0\) is needed.

## Proof and the quadratic coefficient

Put \(E=Y_N-W(t^D)\). The shared-head factorization in the conditional
lemma gives
\[
 \operatorname{ord}E=\gamma_N=-p_N+N+1,\qquad
 J_{N,k}=\partial_YH_k(P_N,W(t^D))
                  =\lambda_kt^{L_{N,k}}+\text{gap at least }D.
\]
All coefficient exponents of these derivatives are in \(D\mathbb Z\).
Let \(a=[t^{\gamma_N}]E\) and \(b=[t^{\gamma_N-1}]E\).

Expand the final polynomial at its common root:
\[
 H_I(P_N,W(t^D)+E)=J_{N,I}E+K_{N,I,2}E^2+\cdots .
\]
The earlier-prefix quadratic contributions have an additional gap
of at least \(D\); the only top quadratic contribution comes from
the last power. Its leading coefficient and exponent are
\[
 \kappa=\binom p2\epsilon^{p-2}\lambda_{I-1}^2,\qquad
 \operatorname{ord}K_{N,I,2}=L_{N,I}+p_N-D.                    \tag{2}
\]
The exponent follows as well from
\(q_{N,I-1}=L_{N,I-1}-p_N+D\), which is equivalent to the final
\(N\)-gap and order equality.

The linear term contributes
\[
 [t^{r_N}]J_{N,I}E
 =\lambda_I[t^{-p_N}]E
 =\lambda_I\tau_N-\lambda_I[z^{-p/2}]W.
\]
Any lower derivative term would require an exponent of \(E\) above
\(\gamma_N\), since \(D>N+1\).

The quadratic product has top possible order \(r_N+1\). At \(r_N\),
only the leading coefficient in (2) can contribute: its lower
coefficients have gap at least \(D\). The required coefficient of
\(E^2\) is
\[
                      [t^{2\gamma_N-1}]E^2=2ab.
\]
All terms of degree at least three have order at most
\(q_{N,I}-2N<r_N\), because \(N>1\).
At the final head and head-minus-one, all nonlinear terms are too
low, so
\[
                              A_N=\lambda_Ia,\qquad
                              B_N=\lambda_Ib.
\]
Finally
\[
 \frac{2\kappa}{\lambda_I^2}
 =\frac{p-1}{p\epsilon^p}.
\]
Combining these exact contributions proves (1). For an earlier prefix \(k<I\), put \(\ell_k=L_{N,k}/D\).
The corresponding exact linear formula, including the unrestricted
trace coefficient, is
\[
 [t^{L_{N,k}-p_N}]H_k(P_N,Y_N)-\lambda_k\tau_N
 =[z^{\ell_k-p/2}]W_k-\lambda_k[z^{-p/2}]W.
\]
Its nonlinear terms are already strictly below the distinguished order.
With \(\tau_N=0\), this is precisely the earlier restricted formula.

## Nonzero-head controls and application boundary

[The independent Fraction checker](check_exceptional_laurent_correction.py)
keeps the parameters of \(A_N,B_N,\tau_N\) and the last shared head
\(\epsilon\) indeterminate. It verifies twelve finite families at
\(N=3,5,9\), with \(p=2,4\) and first head coefficient 1 or 2.
The \(p=2\) controls have a nonzero fixed baseline
\(-\lambda_I[z^{-p/2}]W\). Cross-multiplication verifies (1) without
numerically specializing the nonzero last head coefficient.
The [receipt](exceptional_laurent_correction_controls.json) pins the
proof and arithmetic sources. These controls support the explicit
formula; the all-\(N\) assertion is proved above.

For the earlier \(p=4,\epsilon=1\) resonant control,
\(A_N=16\) and \(B_N=16v\), so the correction is
\((3/4)16(16v)=192v\), exactly its previously retained varying
distinguished coefficient.

This repairs the scalar implication at the exceptional step of
Zheglov's Lemma 6.3, printed pages 70--71 of the
[pinned primary PDF](sources/zheglov-2410.06959v5.pdf),
arXiv:2410.06959v5. It imports no assertion that a simultaneous
Weyl pair supplies the shared heads, trace, commuting-curve model,
or required conjugations. For actual transport, a monomial-face
gap can independently give \(B_N=0\), even with unequal edge centroids.
That issue and the first-component PBW calculation are treated in
[the separate transport audit](PBW_TRANSPORT_AND_FINAL_COEFFICIENT_COMPARISON.md).
There is no actual-pair construction or whole-theorem conclusion here.
