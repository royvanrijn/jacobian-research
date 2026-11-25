# Polynomial carriers and the third support layer

The [second-layer minimal-repair example](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#correct-heights-and-repair-through-the-second-support-layer) gives balanced radical caps, correct repair heights and simultaneous divisibility by the square of the support radical, with a nonzero additive error. There is a sharp third-layer obstruction for polynomial divisibility identities. For individual integer evaluations the obstruction leaves square-full values of a fixed polynomial factor. A fixed-field curve estimate gives a power saving whenever such a factor has degree at most four.

## A bound on three simultaneous contact orders

Let \(A,B,C\) be pairwise coprime nonconstant polynomials over a field of characteristic zero, in one or more variables. Write \(E=A+B-C\) and \(r_A=\deg\operatorname{rad}(A)\), using total degree and the polynomial radical. For each irreducible \(P\mid A\), put \(m_{A,P}=\min(v_P(A),v_P(E))\), with \(v_P(0)=+\infty\). Then
\[
 \sum_{P\mid A}(m_{A,P}-1)_+\,\deg P\le r_B+r_C-1,
                                                               \tag{PC1}
\]
and there are corresponding inequalities with the letters interchanged. In particular
\[
 \sum_{J=A,B,C}\ \sum_{P\mid J}(m_{J,P}-1)_+\,\deg P
                 \le 2(r_A+r_B+r_C)-3.                         \tag{PC2}
\]

Choose a constant-coefficient directional derivative \(\partial\) for which the derivatives of the three pairwise ratios are nonzero. Such directions exist because each nonconstant rational function in characteristic zero has a nonzero partial derivative, so finitely many proper linear subspaces of directions can be avoided. The polynomial
\[
 W_{B,C}=\operatorname{rad}(B)\operatorname{rad}(C)
                   \left(\frac{\partial B}{B}
                                  -\frac{\partial C}{C}\right)
\]
is nonzero of degree at most \(r_B+r_C-1\). The degree bound follows by expanding the logarithmic derivatives over the irreducible factors, so every fraction in parentheses has a single simple denominator. At a factor \(P\mid A\), \(B/C\) is a unit in the localization at \(P\) and
\[
               B/C-1=(E-A)/C
\]
has order at least \(m_{A,P}\). A derivative lowers this order by at most one, and \(\operatorname{rad}(B)\operatorname{rad}(C)\) is a unit there. Thus \(v_P(W_{B,C})\ge(m_{A,P}-1)_+\), proving (PC1). At \(P\mid B\) use \(A/C-1\), and at \(P\mid C\) use \(A/B+1\). Dividing by these respective unit ratios gives the same logarithmic derivative argument and proves (PC2).

Suppose every irreducible factor of each of \(A,B,C\) has multiplicity at least three in its own polynomial. If \(E=0\), then every \(m_{J,P}\ge3\), which (PC2) excludes. If \(E\ne0\), at least one factor \(P\mid ABC\) has \(v_P(E)\le2\). Consequently there is no simultaneous polynomial divisibility identity \(\operatorname{rad}(ABC)^3\mid E\) under these hypotheses. This polynomial threshold is sharp. For coprime Pythagorean polynomial forms \(\ell,m,n\), take \(A=(\ell^2)^e\), \(B=(m^2)^e\), \(C=(n^2)^e\) at odd \(e\ge3\). The identity \(\ell^2+m^2=n^2\) gives \((\ell mn)^2\mid E\) when the three forms are squarefree, as in the balanced second-layer construction.

## A reduction for integer evaluations

Use fixed \(A,B,C\in\mathbb Z[t]\) under the preceding hypotheses and take a primitive integral irreducible \(P\mid ABC\) with \(j=v_P(E)\le2\). Here \(E\ne0\) automatically. Write \(E=P^j H\), so \(H\) is integral and \(\operatorname{Res}(P,H)\ne0\). Add the prime divisors of this resultant to any fixed prime pool \(S\) to obtain \(S'\). For any integer \(t\) with \(A(t)B(t)C(t)E(t)\ne0\) satisfying
\[
          \operatorname{rad}_{\notin S}
                (A(t)B(t)C(t))^{k}\mid E(t)
                     \qquad(k\ge3),                            \tag{PC3}
\]
every prime \(p\notin S'\) dividing \(P(t)\) satisfies
\[
                        j\,v_p(P(t))\ge k.                    \tag{PC4}
\]
Indeed \(p\) divides the evaluated product \(A(t)B(t)C(t)\), and it does not divide \(H(t)\). Thus for \(j=0\) all prime factors of \(P(t)\) belong to \(S'\), and for \(j=1,2\) the part of \(P(t)\) outside \(S'\) is \(\lceil k/j\rceil\)-full. At the third layer it is in particular square-full. The exact equalities \(E(t)=0\) in this fixed family amount to at most \(\deg E\) integer parameters. Statement (PC4) concerns local acceptance with \(E(t)\ne0\).

## Square-full values outside a fixed prime pool

Let \(P\in\mathbb Z[t]\) be fixed and separable of degree \(d\ge2\), and \(S'\) fixed. The number of \(1\le t\le T\) with \(P(t)\ne0\) and square-full part outside \(S'\) in the sense of (PC4) is
\[
                     O_{P,S',\varepsilon}(T^{d/5+\varepsilon}).
                                                               \tag{PC5}
\]
Here the condition means that every positive valuation outside \(S'\) is at least two. For \(d=1\) the corresponding bound is \(O_{P,S',\varepsilon}(T^{1/2+\varepsilon})\). For a fixed nonzero sign write
\[
                      |P(t)|=b u^2 v^3
\]
with \(b\) a divisor of the product of primes in \(S'\) and \(u,v\) positive integers. Every positive exponent outside \(S'\) is either even, or odd at least three. At the primes in \(S'\), put a first power into \(b\) exactly when the exponent is odd. This proves the decomposition. For \(d=1\) the elementary square-full count follows by summing \((T/v^3)^{1/2}\) over \(v\). For \(d\ge2\), at least one of \(u,v\) is \(O_P(T^{d/5})\). Each value of \(v\) gives a square twist of \(P(t)\) and each \(u\) a cube twist, with nonzero integer coefficient bounded by \(T^{O_P(1)}\). It suffices to give a subpower bound for each such twist.

For \(q=2\) when \(d\ge3\), and \(q=3\) when \(d\ge2\), consider
\[
                         P(t)=a y^q,\qquad y\in\mathbb Z,
                    \qquad 0<|a|\le T^{O_P(1)}.
\]
Take once and for all a number field \(K\) containing the roots of \(P\) and the cube roots of unity. Include in a fixed finite set of places the denominators of the roots, the nonunit places of the leading coefficient and the places where two roots have equal reduction. Then enlarge by the places over \(q a\) to obtain \(\Sigma\). Away from \(\Sigma\), the factors \(t-\alpha\) are integral and at most one has positive valuation. Its valuation is divisible by \(q\) by the twist equation. Choose \(q=2\) with three distinct roots, or \(q=3\) with two distinct roots, and let \(g(t)\) be the product of the chosen linear factors. The class of \(g(t)\in K^\times/K^{\times q}\) has valuations divisible by \(q\) outside \(\Sigma\).

There are at most \(q^{|\Sigma|+O_K(1)}\) such classes. The map to the \(q\)-torsion of the \(\Sigma\)-class group sends an element to the class of the fractional ideal given by its outside valuations divided by \(q\). Its kernel consists of \(\Sigma\)-units modulo \(q\)-th powers, and the \(\Sigma\)-class group is a quotient of the fixed class group. This yields the bound by the unit theorem. Since \(|\Sigma|=O_P(1+\omega(a))\), the number of classes is \(T^{o(1)}\) uniformly in the allowed twists.

For each nonempty class choose an actual parameter \(t_0\) and set \(h=g(t_0)\ne0\). Both \(h\) and every \(g(t)\) have polynomial height in \(T\), and all counted parameters in that class satisfy \(h w^q=g(t)\) for some \(w\in K\).

* At \(q=2\), take \(x=t/h\), \(z=w/h\). The model is \(z^2=\prod_{i=1}^3(x-\alpha_i/h)\), a nonsingular elliptic curve with three nonzero \(K\)-rational points of order two.
* At \(q=3\), take \(x=4hw\), \(z=4h(2t-\alpha_1-\alpha_2)\). The model is
  \[
                z^2=x^3+\bigl(4h(\alpha_1-\alpha_2)\bigr)^2 .
  \]
  The point with \(x=0\), \(z=4h(\alpha_1-\alpha_2)\) has the horizontal tangent of intersection multiplicity three and is a nonzero point of order three.

In each case the coordinate change recovers \(t\) from the elliptic point, and the model and point heights are polynomial in \(T\). Dujella's fixed-number-field bound gives \(T^{o(1)}\) points per class, since a point of prime order two or three is rational over that fixed field. [Dujella, Theorem 1.1 and Corollary 1.3](https://arxiv.org/html/2312.03655v1#S1)

For the square twist when \(d=2\), write \(P(t)=A_2 t^2+B_2 t+C_2\) and \(\Delta=B_2^2-4A_2C_2\ne0\). The equation \(P(t)=a y^2\) gives
\[
                      (2A_2t+B_2)^2-4A_2a y^2=\Delta .
\]
If \(A_2a\) is a square, factor the left hand side over the integers and use the divisors of \(\Delta\). Otherwise remove the square part of \(A_2a\) and view the left hand side as the norm of an integral element in a real or imaginary quadratic field. Its ideal has norm \(|\Delta|\), so there are \(O_\Delta(1)\) possible ideals uniformly over those fields. For any ideal, all possible generators differ by units. In an imaginary field there are boundedly many units. In a real field the logarithmic distance between absolute values of distinct units up to sign is bounded below by an absolute positive constant, as follows from their integral traces and norms \(\pm1\). The two embeddings of the counted generator have polynomial size because \(t,a,y\) do, and inverse polynomial lower bounds by the nonzero fixed norm. This gives \(O_P(1+\log T)\) generators per ideal. Thus also this square twist has a subpower parameter count. Summing over the \(O_P(T^{d/5})\) smaller variables and the finitely many choices of \(b\) and sign proves (PC5).

The required twist bound also follows directly from Evertse–Silverman's
Theorem 1, in the square case for \(d\ge3\) and in the cube case for
\(d\ge2\). For \(P(t)=a y^q\), apply it over \(\mathbb Q\) to
\(a^{q-1}P\) and the places dividing \(a\operatorname{disc}(P)\)
together with the infinite place. The discriminant is a unit away
from those places. One may use the fixed splitting field of \(P\)
for the required two or three roots. At fixed \(P,q\) the bound is
\(\exp(O_{P,q}(1+\omega(a)))\), even without a parameter height cutoff.
This proves the twist input to (PC5) by classical means, with
the quadratic square branch covered by the elementary norm argument
above.
[Evertse and Silverman, Theorem 1](https://ir.cwi.nl/pub/1788/1788D.pdf#page=2)

Wongcharoenbhorn–Meemark's Theorem 1.2 gives exponent about \(0.4769\)
for admissible quadratics and obtains \(2/5\) conditional on an
all-height Mordell curve conjecture. Plymesser's Appendix B gives
the fixed quadratic exponent \(2/5\) by classical estimates,
with coefficient dependence, and Section 10 uses the same split
for \(1+A t^3\) at exponent \(3/5\) with polynomial dependence on
\(A\). The form used here for arbitrary fixed separable factors has
a fixed exceptional prime pool. The Dujella argument supplies an
alternative route at bounded parameter height.
[Wongcharoenbhorn and Meemark](https://doi.org/10.1017/S0004972725000255)
· [Plymesser, Corollary 1.4, Section 10 and Appendix B](https://hexagonmath.org/pdf/2611.00699v1)

## Consequences for third-layer tests

For any fixed one-parameter carrier as above, a factor \(P\) with \(v_P(E)=0\) bounds the number of locally admitted \(t\le T\) by \(O(T^\varepsilon)\), since there are only polylogarithmically many \(S'\)-supported nonzero values of polynomial size. If a factor with \(v_P(E)\le2\) has degree \(1,2,3\) or \(4\), there is an unconditional saving in \(T\) by (PC5), with respective upper exponents \(1/2,2/5,3/5,4/5\) up to \(\varepsilon\).

For higher depths \(k\), another bound for \(j>0\) is \(O_{A,B,C,S,k,\varepsilon}(T^{d/\lceil k/j\rceil+\varepsilon})\). To see this, write a \(D\)-full value outside \(S'\) as \(\prod_{i=D}^{2D-1} a_i^i\) with the last \(D-1\) roots squarefree, by taking valuations modulo \(D\). Once the latter roots are fixed, the first root costs at most the \(1/D\)-power of the remaining height. The sums over the other roots converge since \(i/D>1\). The part on \(S'\) costs only a polylogarithmic factor. Each value allows at most \(d\) parameters. In particular, at fixed \(k>2d\) even a factor of arbitrary fixed degree \(d\) gives a strict saving in \(T\), because \(j\le2\). A test at depth \(k\asymp 1/\lambda\) on all three supports therefore saves in any fixed carrier with bounded \(d\) at sufficiently small fixed \(\lambda\). Such a simultaneous test is a separate condition from raising the repair order on the target alone.

These bounds are on a fixed carrier's integer parameters. For example, a family indexed by \(t\le T\) with repaired height scale \(Y\asymp T^M\) inherits the same saving relative to \(Y^{1/M}\). The hypotheses of the least-repair equivalence alone do not provide a fixed polynomial carrier, a factor of degree at most four, or that parameter-to-height comparison. A general counting argument must also address sparse families that have no such parametrization.
