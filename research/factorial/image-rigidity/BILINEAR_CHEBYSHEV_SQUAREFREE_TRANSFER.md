# A Chebyshev derivative test for bilinear three-windows

Let \(f=Bxy+dx+ey\) with \(Bde\ne0\), and let \(P_i=\Lambda_2(f^i)/i!\). The following assertions hold over every characteristic-zero field.

**Theorem.** If \(e=-d\), then \(f\) satisfies the three-window at every positive starting index, with no valuation restriction on \(B,d\).

More generally, adjoin a square root of \(de\), scale the phase by it, and set
\[
                  x=\frac{d+e}{2\sqrt{de}},\qquad
                  b=\frac{B}{2\sqrt{de}}.
\]
With \(U_n\) and \(T_n\) the Chebyshev polynomials of the second and first kind, respectively, define
\[
                  G_n(x,b)=\sum_{i=0}^n U_n^{(i)}(x)b^i .
\]
A three-zero block ending at \(k=n+1\) is equivalent to \(T_{n+1}(x)=0\) and a multiple zero of \(G_n(x,\cdot)\) at \(b\ne0\). This gives a finite list of axis ratios at each terminal and a single-polynomial squarefreeness target for each ratio.

## An exact consecutive identity

Normalize \(de=1\) by the indicated phase scaling. The [bilinear cross expansion](ZERO_CONSTANT_BINARY_QUADRATIC_MATCHING_AND_BOUNDARY_LIFTS.md#a-valuation-one-lift-for-the-bilinear-three-term-phase) gives
\[
                   \sum_{n\ge0} P_n T^n
                     =\sum_{i\ge0}\frac{i!(2b)^iT^i}
                                           {(1-2xT+T^2)^{i+1}}.
\]
Since \(\sum_{n\ge0}U_n(x)T^n=(1-2xT+T^2)^{-1}\), differentiation in \(x\) proves \(P_n=G_n(x,b)\).

Write \(\mathcal R_b H=\sum_i b^i H^{(i)}(x)\) for a polynomial \(H\). The identities
\[
 \mathcal R_b(xH)=x\mathcal R_b H+b(b\partial_b+1)\mathcal R_b H,\qquad
 \mathcal R_b(T_{n+1})=T_{n+1}+(n+1)bG_n
\]
come from the product rule and \(T'_{n+1}=(n+1)U_n\). Applying them to \(U_{n-1}=xU_n-T_{n+1}\) gives
\[
           G_{n-1}=(x-nb)G_n+b^2\partial_b G_n-T_{n+1}(x).    \tag{1}
\]
The two Chebyshev identities also follow directly by substituting \(x=(z+z^{-1})/2\), using \(T_j=(z^j+z^{-j})/2\) and \(U_j=\sum_{\ell=0}^j z^{j-2\ell}\), and differentiating where the substitution is finite and unramified. They are polynomial identities and extend to all \(x\).

The [exact bilinear boundary recurrence](ZERO_CONSTANT_BINARY_QUADRATIC_MATCHING_AND_BOUNDARY_LIFTS.md#a-valuation-one-lift-for-the-bilinear-three-term-phase) gives
\[
                     P_{n+1}=2(n+1)bP_n+P_{n-1}+2T_{n+1}(x).
\]
Hence three zeros force \(T_{n+1}(x)=0\). Equation (1) then makes the two earlier zeros equivalent to \(G_n=\partial_bG_n=0\) at the same \(b\ne0\). Conversely these conditions give all three zero moments. The polynomial \(G_n\) has degree \(n\) in \(b\), with leading coefficient \(2^n n!\), so this is the ordinary repeated-root condition.

## The opposite-axis ratio

When \(e=-d\), the normalized \(x\) is zero. The boundary condition \(T_{n+1}(0)=0\) forces \(n=2m\) even. Here \(m\ge1\) because the block starts at a positive index. In the cross expansion the axis denominator is \(1+T^2\), so only even \(i\) contribute at index \(2m\). Their coefficients are
\[
             G_{2m}(0,b)=
             \sum_{r=0}^m(-1)^r\frac{(2m-r)!}{r!}(2b)^{2m-2r}
               =(2b)^{2m} B_{m,m+1}\bigl(-(2b)^{-2}\bigr),     \tag{2}
\]
where
\[
 B_{m,\alpha}(z)=\sum_{j=0}^{m}\binom mj(\alpha)_j j!\,z^{m-j}.
\]
The coefficient identification uses \((m+1)_j=(m+j)!/m!\). For \(m=1\) the polynomial \(B_{1,2}\) is linear. For \(m\ge2\) the [Robin resolvent squarefreeness theorem](MULTIVARIATE_ONE_PRIVATE_ROBIN_RESOLVENT.md), \(B_{m,\alpha}\) is squarefree for \(\alpha>m-1\), which includes \(\alpha=m+1\). The change \(b\mapsto -(2b)^{-2}\) has nonzero derivative at every \(b\ne0\), so (2) has no multiple zero there. This excludes the three-zero block.

All coefficient and squarefreeness statements used in (2) are over \(\mathbb Q\). The rational coprimality identity for \(B_{m,m+1}\) and its derivative applies after the indicated algebraic scalar extension, so the exclusion holds over every characteristic-zero field.

## Exact checks

The [finite replay](verify_bilinear_chebyshev_squarefree_transfer.py) checks 150 instances of the coefficient identity (1) and 860 Robin coefficients in (2) by exact arithmetic. The uniform squarefreeness input is the cited resolvent theorem.
