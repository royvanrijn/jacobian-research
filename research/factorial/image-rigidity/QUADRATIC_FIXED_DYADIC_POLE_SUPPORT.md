# Fixed pole support at every dyadic precision for quadratic moment series

Let
\[
 f=\sum_i(\alpha_iX_i^2/2+d_iX_i)+\sum_{i<j}\beta_{ij}X_iX_j
\]
have zero constant and integral natural coefficients at a
place above two. Write \(E(T)=\sum_{n\ge0}\Lambda_m(f^n)T^n/n!\),
and \(U_i(T)=1-\alpha_iT-d_i^2T^2\).
The valuation ring is \(\mathcal O\) with \(v(2)=1\).
All statements allow ramified fields.

**Theorem.** At every integer precision \(k\ge1\), the
series \(E\) is congruent modulo \(2^k\mathcal O[[T]]\)
to an integral rational series with denominator dividing
both
\[
        \prod_i U_i^{\,2k-1}
       \qquad\text{and}\qquad
        \left(\prod_i U_i\right)L^{\,4k-4},
                   \quad L=\operatorname{lcm}_i U_i.        \tag{1}
\]
The least common multiple is normalized to have constant
coefficient one. Thus the set of possible poles is fixed
at all precisions, and the excess denominator growth is
controlled by \(L\) even when many vertices share the
same factors. Modulo \(4\), the
[explicit quadratic correction](QUADRATIC_MATCHING_SECOND_DYADIC_LAYER.md)
gives sharper leading pole orders and identifies their
graph contributions.

## Derivative weights at one vertex

Write \(U=1-\alpha T-d^2T^2\),
\(F=(1+dT)/U\), and \(Q=J_0\), using the one-variable
\(J_r\) of the
[second matching layer](QUADRATIC_MATCHING_SECOND_DYADIC_LAYER.md#the-single-vertex-correction).
The exact single-vertex recurrence gives
\[
 Q=F+2\mathcal L Q,\qquad
      \mathcal L=\frac{\alpha T^2}{U}\frac{d}{dT}.
\]
Hence \(Q\equiv\sum_{a=0}^{M-1}2^a\mathcal L^a F
\bmod 2^M\mathcal O[[T]]\), for every \(M\ge1\).
The numerator of \(\mathcal L^a F\) is integral
with denominator \(U^{2a+1}\), by differentiation.

Let \(D_T=T\partial_T+\alpha\partial_\alpha\).
Before evaluating the vertex coefficients, the exact
operator identity is
\[
                 r!J_r=\prod_{j=1}^r(j+D_T)Q .             \tag{2}
\]
In fact \(D_T\) multiplies the \(T^{a+b}\alpha^a d^b\)
coefficient by \(2a+b\). Equation (2) multiplies it by
\((2a+b+1)_r\), precisely the factorial weight of the
insertion \(X^r\). This uses polynomial derivatives on
the universal coefficients, so it preserves all the
coefficientwise congruences. Each step in (2) can raise
the denominator exponent by one only at the cost of
another factor of two, since
\[
                        D_T U=-2(\alpha T+d^2T^2).
\]
Thus applying the operator in (2) to the \(a\)-th term
from \(Q\) gives a finite sum of integral rational terms
of the form
\[
                      \frac{2^{a+b}P_{a,b}(T)}{U^{2a+1+b}},
                          \qquad 0\le b\le r .              \tag{3}
\]
One obtains this directly by keeping the factors of two
when a derivative lands on the denominator. Numerator
differentiations preserve the already extracted powers.

## The global precision and the common denominators

For a cross count array \((j_e)\), set \(r_i=\sum_{e\ni i}j_e\)
and \(W=\prod_i r_i!/\prod_e j_e!\).
The exact cross terms are
\[
      W T^{\sum j_e}\prod_e \beta_e^{j_e}\prod_i J_{r_i,i}.
\]
The multiplier \(W\) is integral. Put \(w=v_2(W)\) and
\(S=\sum_e v_2(j_e!)\). The
[carry identity](QUADRATIC_MATCHING_SECOND_DYADIC_LAYER.md#the-graph-terms-of-valuation-one)
gives \(S\le w\). In addition \(v_2(r_i!)\le w\) at
every vertex: cancel the incident edge factorials in
\(W\) against the vertex multinomial at \(i\); the
remaining factors outside \(r_i!\) are integers.
Thus if \(w<k\), all \(j_e\le2k-1\), and there are
only finitely many contributory cross arrays for a
fixed precision. Arrays with
\(w\ge k\) contribute zero modulo \(2^k\).

At a vertex use (2) with the truncation of \(Q_i\)
at \(M_i=k-w+v_2(r_i!)\). The division by \(r_i!\)
then leaves precision \(k-w\). Because the resulting
approximation to \(J_{r_i,i}\) has integral formal
coefficients, these approximations can be multiplied
and weighted by \(W\) to approximate \(E\) modulo \(2^k\).
In (3) for a vertex we can first discard every term
with \(a+b\ge M_i\), at that same precision after
division by \(r_i!\).
The remaining denominator exponents satisfy
\[
                  2a+1+b\le2M_i-1\le2k-1,
\]
proving the first bound in (1).

For the second bound, distribute the rational summands
rather than first clearing the vertex denominators.
A product in (3), after division by every \(r_i!\)
and multiplication by \(W\), has scalar valuation at
least \(\sum_i(a_i+b_i)-S\). The cross powers may
increase the valuation. Thus it suffices to retain
\(\sum_i(a_i+b_i)<k+S\). Its additional denominator powers
beyond \(\prod U_i\) have total exponent at most
\[
           \sum_i(2a_i+b_i)\le2(k-1+S)\le4k-4.
\]
Each extra \(U_i\) divides \(L\), proving the second
bound in (1). For a given cross array one may use the
sharper bound \(2(k-1+S)\).

To see integrality of these common denominators, in a
splitting extension write each \(U_i\) as a product
of factors \(1-\mu T\). The \(\mu\) are integral by the
reversed monic polynomials. Products taking maximum
multiplicity at each reciprocal root give \(L\)
and the normalized common factors over the original
field by Galois invariance, with integral coefficients.
The two denominators in (1) have constant one, and
their normalized greatest common divisor is integral
by the same argument. Multiplying by a common
denominator gives an integral numerator: the resulting
polynomial is integral as a formal series, since the
rational approximation is congruent to the integral
series \(E\). This also justifies regrouping individual
terms of negative scalar valuation in the second
argument.

For example, if all nonconstant vertex denominators
divide a fixed quadratic, the second bound adds at most
\(8k-8\) degrees beyond \(\deg\prod_i U_i\), irrespective
of the dimension or the number of cross coefficients.
For a chosen precision one can use the resulting pole
multiplicities in a coefficient recurrence or in a
unit-pole removal test. The theorem controls the
rational approximation denominators; a zero window
still depends on its numerator and the relevant
coefficient precision.
