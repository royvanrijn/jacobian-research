# Further consequences of the octahedral conditional laws

The notation and exact normalized laws come from [the octahedral compact extension](MATHIEU_CLASSIFICATION_FIELD_AND_COMPACT_EXTENSIONS.md#the-two-octahedral-arrangements-and-a-uniform-conditional-coordinate) and [the stripped Hopf calculation](HOPF_DIRICHLET_LAW_AND_QUINTIC_SUBALGEBRA_OBSTRUCTION.md).

## Trace transfer with regular angular markers

Write \(p,q\) for the octahedral reflection coordinates and \(J\) for the orientation coordinate. In the two-arrangement law of (H14)–(H15), put \(Z=1-X-V\) and
\[
 B_2=Z-1=-(X+V),\qquad C_2=-X.
\]
Then \((-C_2,C_2-B_2,1+B_2)=(X,V,Z)\). The two independent coordinates have full support on the open triangle in the real model and so are algebraically independent over \(\mathbb Q\). This independence persists over any characteristic-zero field by extending the rational coordinate algebra. The assignment \(B\mapsto B_2,C\mapsto C_2\) is a trace-preserving isomorphism from the stripped Hopf polynomial algebra \(K[B,C]\) onto the rational octahedral height algebra \(K[X,V]\): the moment formula (4) and the two-arrangement Dirichlet law agree on all monomials. The corresponding conclusion including polynomial markers goes in both directions. In particular, every normalized phase
\[
 (X+V)X+\alpha X^2-\beta(X+V)-\gamma X+\delta
\]
has a nonzero pure moment at some order from one through five.

To project regular angular markers, take the recurrence \(C_n(Z)\) from the two-arrangement law. With \(r=1-3p\), its conditional moments are
\[
 \mathbf E[r^n\mid X,V]=\left(\frac{1-X}{2}\right)^n C_n(Z).
\]
For nonnegative integers \(i,j\), use the exact identity
\[
 p^i q^j=\frac{X^j}{9^j3^{i+j}}(1-r)^{i+j}
\]
to obtain
\[
 P_{ij}(X,V):=\mathbf E[p^i q^j\mid X,V]
 =\frac{X^j}{9^j3^{i+j}}
    \sum_{n=0}^{i+j}(-1)^n\binom{i+j}{n}
        \left(\frac{1-X}{2}\right)^n C_n(Z)\in\mathbb Q[X,V].
\]
Every reflection-invariant polynomial marker has a polynomial conditional trace by linearity. A regular rotation-invariant marker has the form \(H_0(p,q)+JH_1(p,q)\). Inversion fixes \(X,V,p,q\) and reverses \(J\), so the odd summand has zero conditional trace. For coefficients of a central Laurent marker, apply the same projection to each angular coefficient and retain the central powers. Consequently all mixed moments of any phase in \(K[t,t^{-1},X,V]\) against these regular angular Laurent markers agree with those against their projected height-algebra markers. The constant term in \(t\) commutes with the angular projection because the expansions are finite in \(t\). Thus eventual vanishing against every height-algebra Laurent marker implies eventual vanishing for every regular reflection- or rotation-invariant angular Laurent marker. The direct transfer does not require regularity of the rational phase itself. The reverse equality of mixed moments uses relative markers on the two isomorphic Dirichlet algebras. No surjectivity of the projection from regular angular markers is asserted.

There is also a weighted function-field description of the projection. At fixed \(X,V\) the two sheets exchange \(S\) with \(S^{-1}\) and fix \(B\). The two weights computed in the source are \(S/(S+1)\) and \(1/(S+1)\). Thus for any bounded rational angular function its conditional value is
\[
 \frac{S f(S,B)+f(S^{-1},B)}{S+1},
\]
the field trace weighted by \(S/(S+1)\). The polynomial projection of regular markers is genuinely smaller than the whole height-marker algebra: each displayed \(P_{ij}(X,V)\) specializes at \(X=1\) to a constant independent of \(V\), so \(V\) itself is outside its image. This makes the independent-marker restriction in the isomorphism explicit.

The analogous result for the earlier law (H12) follows by putting \(U=(1-v)R\). Conditional on \(v,U\), the two values of \(p-A(v)\) are equally likely with square \((1-v)^2(9-v)U/576\). Expand \(p^i q^j=9^{-j}v^j[A(v)+(p-A(v))]^{i+j}\). The odd terms vanish and the even terms are polynomial in \(v,U\). The same orientation and scalar Laurent arguments apply. All transfer identities hold over characteristic zero by the rational formulas for the resulting moments. The use of real conditional expectations provides those identities for rational polynomial coefficients, followed by scalar extension.

## An exact three-moment criterion on the affine Dirichlet sector

Let \(D_1,D_2,D_3\) be the three barycentric coordinates in any of the three Dirichlet algebras above, and \(L=\sum_i\lambda_iD_i\). This includes every affine phase because \(\sum_iD_i=1\). The Dirichlet moment formula gives the formal series
\[
 \sum_{m\ge0}\frac{(3/2)_m}{m!}\Lambda(L^m)z^m
   =\prod_{i=1}^3 (1-\lambda_i z)^{-1/2}.
\]
If the first three positive moments vanish, the series is \(1+O(z^4)\). The inverse square of this series is the polynomial \(\prod_i(1-\lambda_i z)\). It must be \(1\), so the three elementary symmetric polynomials in the \(\lambda_i\) vanish and every \(\lambda_i=0\). In fields containing a primitive cube root of unity \(\zeta\), the choice \((1,\zeta,\zeta^2)\) has first two moments zero and third nonzero. Thus three is the sharp bound uniform in characteristic-zero coefficient fields. This is a finite criterion just for the affine barycentric sector.

## Exact scalar marker formulas for the full-field rational endpoint algebra

Fix \(N\ge1\), let \(\mathcal A=K[E,W,W^{-1},B]\) in (H16) and put \(g(U)=G(U)^N\). The elements \(E,W,B\) are algebraically independent because their field is \(K(t,S,B)\), with inverse formulas in (H14)–(H16), and \(S,B\) are two independent angular coordinates. Let \(\nu\) denote the normalized product measure of \(S=(\sqrt S)^2\) and \(B\) from (H14). Its Laurent polynomial moments are rational. To see this for \(S^j\), substitute \(u=\sqrt2\sqrt S\), which is uniform on \([1,2]\). Then \(S^j=(u^2/2)^j\), and the antiderivative of \(u^{2j}\) has denominator \(2j+1\ne0\) for integral \(j\). The \(B\) moments are arcsine moments.

Let \(r\ge1,\ j\in\mathbb Z,\ \ell\ge0,\ k\in\mathbb Z\). The endpoint substitution of (H16) gives a prefactor \(t^{k-r}T^{r-1}\), with \(T=1+t\), and replaces the variable integration interval by the original one with error of order at least \(Nr+1\). In particular, for \(k\ge1\),
\[
 \Lambda(t^k E^r W^j B^\ell)=
 \binom{r-1}{k-1}\,\nu\!\left(g(S)^r S^j B^\ell\right),
\]
where the binomial is zero if \(k>r\). For \(k=0\) the coefficient is zero because \(T^{r-1}\) has degree \(r-1\). For \(k=-d<0\), if \((N-1)r\ge d\), even the endpoint error contributes zero and the original interval multiplies the missing coefficient of \(T^{r-1}\) at order \(r+d\). Thus these negative-marked moments vanish under the displayed order condition. The same endpoint calculation works over every characteristic-zero field by rational coefficients.

For \(H\in E\mathcal A\) and \(k\ge1\), the positive-marker formula can be written as
\[
 \Lambda(t^k H)
  =\nu\!\left[
      \frac{g(S)^k}{(k-1)!}
      \left.
      \frac{\partial^{k-1}}{\partial E^{k-1}}\frac{H(E,W,B)}{E}
      \right|_{E=g(S),\,W=S}
    \right].
\]
In particular the fixed marker \(t\) reads evaluation at \(E=g(S)\), \(W=S\).

Over any formally real field, purity in \(\mathcal A\) holds exactly for \(F\in E\mathcal A\). Indeed, \(\Lambda(F^m)=\nu(F(0,S,B)^m)\) since all higher \(E\) powers have zero unmarked trace. For any finite Laurent polynomial support, the Gram matrix for \(\nu\) is rational and positive definite: a nonzero Laurent polynomial in two independent real variables cannot vanish on the open rectangle where the measure has positive density. Positive definiteness over any ordered coefficient field follows from the same rational positive principal minors. Thus the second moment already forces \(F(0,S,B)=0\). Centered \(B-\frac12\) shows that the first moment alone does not decide purity in general.

Moreover, for a pure \(F\in\mathcal A\) over a formally real field, eventual marked vanishing against every element of \(\mathcal A[t]\) is equivalent to \(F\in E(E-g(W))\mathcal A\), for every \(N\ge1\). If \(F(E=g(W))\ne0\), the marker \(t\) has nonzero moment at every positive even order by the Gram argument and the evaluation identity. If the evaluation vanishes, write \(F=E F_1\); the domain property and \(g(W)\ne0\) imply \(F_1(g(W),W,B)=0\), so division by the monic polynomial \(E-g(W)\) gives the indicated ideal. For \(F\) in that ideal and fixed \(Q\in\mathcal A\), the derivative above applied to \(H=F^m Q\) vanishes for \(m\ge k\) because \(H/E\) is divisible by \((E-g(W))^m\). The case \(k=0\) has zero trace since \(H\in E\mathcal A\). Finite linear combinations give all the specified polynomial scalar markers.

For \(N\ge2\), this equivalence also holds against every marker in \(\mathcal A[t,t^{-1}]\). A fixed marker \(t^{-d}Q\) is covered once \((N-1)m\ge d\) because every term of \(F^m Q\) has \(E\)-exponent at least \(m\). The reverse uses the same fixed marker \(t\). This extension concerns the rational endpoint algebra and the stated scalar marker ring. The established regular-intersection result continues to govern descent to the original regular compact product.

## The original regular angular markers and the exceptional case \(N=1\)

The original reflection-invariant regular Laurent markers \(K[t,t^{-1},p,q]\) belong to \(\mathcal A[t,t^{-1},T^{-1}]\). Indeed, \(S=W/T^2\), and (H14) writes \(p,q\) as Laurent polynomials in \(S\) and polynomials in \(B\). Multipliers can therefore be expanded as finite sums of \(t^kT^aQ(E,W,B)\), with \(k,a\) integers and \(Q\in\mathcal A\).

If \(F=E(E-g(W))F_1\), the change of variables \(W=s^2\) makes \(E=(T/t)g(s^2)\) and \(E-g(W)=g(s^2)/t\). In a monomial whose remaining factor \(F_1^mQ\) has \(E\)-exponent \(j\ge0\), the scalar prefactor after substitution is
\[
 t^{k-2m-j}T^{a+m+j-1},
\]
and the integrand contains \(G(s^2)^{N(2m+j)}\). For \(m\) large enough that \(a+m-1\ge0\) and \(m\ge k+a\), the required coefficient of \(T^{a+m+j-1}\) is at order \(2m+j-k\), strictly above its degree. Corrections from the moving endpoints, after multiplication by the prefactor, have order at least \(k+(N-1)(2m+j)+1\). This is positive for \(m\) large when \(N\ge2\). This proves eventual vanishing for the entire enlarged marker ring, hence in particular for the original regular angular markers. Odd orientation markers vanish by inversion.

For \(N=1\), the same sufficiency holds after requiring one further factor \(G(W)\) in \(F\). The extra factor in \(F^m\) adds \(m\) to the endpoint vanishing order without increasing the pole order in \(t\).

Conversely for \(N=1\), assume already \(F\in E(E-G(W))\mathcal A\) and suppose its restriction \(F(E,w_i,B)\) to one endpoint \(w_i\in\{1/2,2\}\) is nonzero. The fixed regular multiplier \(t^{-1}p\) selects \(w_i=1/2\), whereas \(t^{-1}(p-1/4)\) selects \(w_i=2\). To check this, set \(s_i=\sqrt{w_i}\) and write \(s=s_i(1+tu)\) on a moving endpoint segment, with \(0\le u\le1\). The factor \(t^{-1}\) is canceled by \(ds=s_i t\,du\). At \(t=0\) we have \(W=w_i\), \(T=1\), and \(E=2w_iG'(w_i)u\). The angular coordinate \(p\) takes the constant value \(1/4\) at the lower endpoint and \(0\) at the upper endpoint, by (H14). The fixed-interval contribution vanishes for sufficiently large \(m\) by the prefactor calculation in the preceding paragraph, with the fixed powers of \(T\) coming from the multiplier. The other moving endpoint contribution vanishes by the chosen multiplier. The selected endpoint contributes a nonzero rational constant times
\[
 \int_0^1\int_0^1\!
      F\bigl(2w_iG'(w_i)u,w_i,B\bigr)^m\,d\nu_{\rm arc}(B)\,du .
\]
The multiplier values, normalized endpoint factors \(s_i/(\sqrt2-1/\sqrt2)\), and all moments in this expression are rational. The factor \(2w_iG'(w_i)\) is nonzero at both simple endpoints. The same positive rational Gram argument over any formally real field shows that the displayed integral is nonzero at every positive even order when the restriction is nonzero. This proves necessity of vanishing on both endpoint divisors, hence of divisibility by \(G(W)\). If \(F=E(E-G(W))F_1\), reduction at either \(w_i\) sends \(E(E-G(W))\) to \(E^2\), so both linear factors \(W-w_i\) divide \(F_1\) as well.

Together with the scalar markers \(1,t\), these arguments classify every safe \(F\in\mathcal A\) over a formally real field for the original regular angular Laurent markers. The ideal is \(E(E-g(W))\mathcal A\) for \(N\ge2\), and \(E(E-G(W))G(W)\mathcal A\) for \(N=1\). The same classifications hold for the larger \(\mathcal A[t,t^{-1},T^{-1}]\) marker ring, and for \(\mathcal A[t,t^{-1}]\) in the \(N=1\) case as well: for the last assertion replace \(p\) or \(p-1/4\) in the necessity proof by a linear polynomial in \(W\) vanishing at the other endpoint. The phases in these assertions are rational endpoint phases. Regular phases in the endpoint algebra still satisfy the separate regular-intersection theorem.
