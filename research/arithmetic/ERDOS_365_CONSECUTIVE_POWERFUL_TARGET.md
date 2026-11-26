# Erdős #365: a subquarter bound for consecutive powerful numbers

Let \(N(X)=\#\{n\le X:n,n+1\text{ are powerful}\}\), where all
positive prime valuations of a powerful integer are at least two. We prove
\[
 N(X)\ll_\varepsilon X^{12127/49920+\varepsilon}
 \qquad (X\ge2,\ \varepsilon>0).                         \tag{1}
\]
The exponent is approximately \(0.242929\). The proof combines
Plymesser's published preprint estimate for reduced determinant windows
with a new use of a fixed quadratic norm fiber and the first moment of
quadratic class group \(3\)-torsion. The polylogarithmic bound in
[Erdős #365](https://www.erdosproblems.com/365) remains open.

[Reuss, Theorem 4](sources/reuss-consecutive-1212.3150v2.pdf) gives
the exponent \(29/100\).
[Plymesser, Theorem 1.1 and Remark 9.1](sources/plymesser-powerful-products-hexagon2611.00602-v1.0.1.pdf)
give \(1/4\) and identify the binding corner of that proof. There one
cube kernel is bounded, the other has size \(X^{1/4}\), and its
square variable has size \(X^{1/8}\). The fiber estimate below handles
both orientations of the corner. The optimization then uses the same
estimate throughout the box domain.

## Fixed fibers from a trace equation

**Proposition (Uniform fixed fiber).**

Let \(d,v\ge1\), with \(d\) squarefree, and let
\(\sigma=\pm1\). For \(K=\mathbb Q(\sqrt{-\sigma d})\) quadratic write
\(h_3(K)=|\operatorname{Cl}(\mathcal O_K)[3]|\). Use \(h_3=1\) in
the split case \(d=1,\sigma=-1\). The number of positive integer
solutions \((u,e)\) of
\[
 d^3u^2+\sigma=v^2e^3
\]
is \(O(2^{\omega(v)}h_3(K))\) with an absolute constant. In particular it
is \(O_\varepsilon(d^{1/2+\varepsilon}v^\varepsilon)\).

In the quadratic case put \(\alpha=1+du\sqrt{-\sigma d}\). Its trace
is \(2\) and its norm is \(\sigma v^2e^3\). Factor its nonzero principal
ideal as \((\alpha)=\mathfrak D\mathfrak A^3\) with cube free integral
ideal \(\mathfrak D\) and integral ideal \(\mathfrak A\). A prime above
an odd rational prime cannot divide both \(\alpha\) and its conjugate,
by the trace identity. Hence an odd prime dividing \((\alpha)\) lies
above a split prime \(p\), and at most one of the two primes above
\(p\) occurs. Its exponent is \(2v_p(v)+3v_p(e)\) by the norm identity.
The residue modulo \(3\) is fixed by \(v\), with at most two possible
choices of prime when this residue is nonzero. There are at most nine
choices of residual exponents at the primes above \(2\). Thus there
are at most \(9\cdot2^{\omega(v)}\) possibilities for \(\mathfrak D\).

For each \(\mathfrak D\), the class of \(\mathfrak A\) solves
\([\mathfrak A]^3=[\mathfrak D]^{-1}\). Each nonempty fiber of this map
has \(h_3(K)\) elements. Fix an integral ideal \(\mathfrak J\) for
one such class and a generator \(\gamma\) of
\(\mathfrak D\mathfrak J^3\). Then \(\mathfrak A=\beta\mathfrak J\) for
some \(\beta\in\mathfrak J^{-1}\), and
\(\alpha=\eta\gamma\beta^3\) for a unit \(\eta\).
An arbitrary cube in the unit can be absorbed into \(\beta\), so
Dirichlet's unit theorem leaves at most nine choices for \(\eta\)
(one free generator at most and a cyclic group of roots of unity).

Write \(L(r,s)=r w_1+s w_2\) in a \(\mathbb Z\)-basis of
\(\mathfrak J^{-1}\). The trace equation becomes
\[
 F(r,s)=2,\qquad F=\operatorname{Tr}_{K/\mathbb Q}(\eta\gamma L^3).
\]
This binary cubic has integral coefficients because
\(\gamma(\mathfrak J^{-1})^3=\mathfrak D\subseteq\mathcal O_K\).
The two linear forms \(L,\overline L\) are independent over \(K\):
otherwise \(w_1/w_2\) would be rational. Over \(K\) the cubic is an
invertible linear change of a diagonal binary cubic with two nonzero
coefficients. It has distinct projective roots in characteristic zero.
Its content divides \(2\) if there is a solution, and every solution
has \(\gcd(r,s)=1\), since its cube divides \(2\).
[Stewart, the consequence of Theorem 1 on p. 796](https://uwaterloo.ca/pure-mathematics/sites/default/files/uploads/documents/s0894-0347-1991-1119199-x_0.pdf)
therefore bounds the number of pairs by an absolute constant after
dividing out the content. The theorem allows reducible forms with
nonzero discriminant. Each pair determines \(\alpha\), hence \(u\)
and at most one positive \(e\).

In the split case \(u>1\) and the two nonzero integers \(1+u,1-u\)
have sum \(2\) and product \(-v^2e^3\). Factor each as \(a_j b_j^3\)
with \(a_j>0\) cube free and \(b_j\in\mathbb Z\), retaining the sign
in \(b_j\). The same odd prime argument and the at most nine choices
at \(2\) give \(O(2^{\omega(v)})\) coefficient pairs. The equation
\(a_1b_1^3+a_2b_2^3=2\) again meets Stewart's hypotheses. This covers
both signs and the exceptional square radicand. Finally
\(h_3(K)\le h_K\ll_\varepsilon d^{1/2+\varepsilon}\) by Minkowski's bound
and counting ideals of norm up to \(O(\sqrt d)\) with the divisor
bound. Also \(2^{\omega(v)}\ll_\varepsilon v^\varepsilon\).

A related general input is
[Evertse–Silverman, Theorem 1(a), pp. 237–238](sources/evertse-silverman-superelliptic-counts-1986.pdf).
Take the base field \(\mathbb Q\), \(f(T)=(d^3T^2+\sigma)/v^2\),
and \(S\) consisting of infinity and primes dividing \(2dv\).
Its discriminant \(-4\sigma d^3/v^4\) is an \(S\)-unit and both
roots lie in \(K\) (in \(\mathbb Q\) in the split case).
That theorem gives \(O_\varepsilon((dv)^\varepsilon h_3(K))\),
which also suffices for (1). The trace proof keeps the explicit
\(2^{\omega(v)}\) dependence and avoids including primes dividing \(d\)
in the branch budget.

## Averaging over the cube kernel

For each sign the fundamental discriminant of \(K\) has absolute
value at most \(4d\), and the squarefree radicands parametrize distinct
quadratic fields. The first moment in
[Bhargava–Shankar–Tsimerman, Theorem 2](sources/bst-2013-davenport-heilbronn.pdf)
therefore gives \(\sum_{d\le D}h_3(K)\ll D\) over the squarefree
kernels, including the split convention. For any \(D,V\ge1\), the
number of solutions with \(d\le D,\ v\le V\) is consequently
\(O_\varepsilon(DV^{1+\varepsilon})\) for either sign. For fixed
\(D\) and \(V=X^{1/8}\) this gives the bounded-kernel
\(O_{D,\varepsilon}(X^{1/8+\varepsilon})\) bound.

## Global bound from four estimates

Every positive powerful number has a unique representation \(d^3u^2\)
with \(d\) squarefree and \(u>0\). Use \(e^3v^2-d^3u^2=1\) and divide
\(n=d^3u^2\) into dyadic ranges of size comparable to \(x\). Group the
kernels at dyadic scales \(D,E\), and apply Plymesser's Lemma 2.3 to use
\(O(1)\) exact windows per group, with comparable \(D,E,x\) and
\(U=x^{1/2}D^{-3/2},\ V=x^{1/2}E^{-3/2}\).
To use open endpoints in all four variables, round \(3d/4,3e/4,n/4\)
down to the grid \(\{2^{j/8}:j\in\mathbb Z\}\) to obtain \(D,E,x\).
Then \(d/D,e/E\) are strictly between \(1\) and \(2\), and
\(u/U=(D/d)^{3/2}(n/x)^{1/2}\) lies between \(1\) and \(2\)
because \((3/4)2^{-1/8}<D/d\le3/4\) and \(4\le n/x<4\cdot2^{1/8}\).
The same bounds with the extra factor \(((n+1)/n)^{1/2}\) give
\(1<v/V<2\) for \(n\ge2\). This grid still uses \(O(1)\)
choices per dyadic group.
Write the kernel exponents as \(a,b\), and set
\[
 s=a+b,\quad t=1-3s/2,\quad m=9st/8,\qquad
 \bar u=(1-3a)/2,\quad \bar v=(1-3b)/2.
\]
A nonempty window has exponents within \(O(1/\log x)\) of the square
\([0,1/3]^2\). Clipping to the square changes every expression below
by \(O(1/\log x)\), which costs only a constant factor. Assume
\(a\le b\) by exchanging the pairs and using the corresponding sign
in the fiber. The determinant input can always be applied to the
original orientation with shift \(1\) since all of its exponents
are symmetric. Each counted solution has \(\gcd(d,e)=\gcd(u,v)=1\),
by the shift equation.

The available box count exponents, up to \(\varepsilon\), are

1. \(s\) and \(t\) by Plymesser's Lemma 4.1 with shift \(1\). The
   two estimates fix the kernel pair and the square variable pair
   respectively. The first is the uniform quadratic norm count in a
   polynomial height range. The second uses Stewart on the binary cubic
   with right-hand side \(1\), whose coefficients are coprime since
   \(\gcd(u,v)=1\).
2. \(a+\bar v\) by the averaged fixed fiber, with exponent \(O(\varepsilon)\)
   absorbed by applying that estimate with smaller \(\varepsilon\).
3. \(\max(m,B)\) by Plymesser's Theorem 6.3 when the actual
   \(s,t\ge 1/10\), where
   \[
     B=m+\tfrac3{10}\max((s-m)/2,a-m,b-m)
          +\tfrac15\max((t-m)/2,\bar u-m,\bar v-m).
   \]

The third input uses Reuss's **Lemma 7**, with \(k=3,l=2,h=1\),
and \(M=\lceil x^{(1+\delta)m}\rceil\) for a fixed sufficiently small
\(\delta>0\). In the region where it is needed below,
\(s,t\ge 6/25-O(1/\log x)\), so \(\log M\asymp\log x\),
\(M\le x^{1/2}\) and \(x\ge 3M\) for large \(x\).
The construction in Sections 5–6 uses the two reduced lattices and
primitive pairs with gcd set \(\{1\}\). The fixed ratios use Lemma 4.1.
The bidegree-one cases use Proposition 5.3 with bounded bidegree and
polynomial coefficient height. The remaining components have bidegree
\((\alpha,\beta)\) with \(\alpha,\beta\ge2\). Lemma 5.5 embeds their
actual integral coordinates in dimension 12 on an affine curve of
degree \(2\alpha+3\beta\ge10\). This is the input for
[Castryck–Cluckers–Dittmann–Nguyen, Theorem 3 and Remark 4.1.5](https://arxiv.org/pdf/1904.13109v1).
The coefficient height bound uses Gelfond's inequality as in
Plymesser's Theorem 3.5. The first minimum tails of Lemmas 6.1–6.2
are at this same \(M\), and the two moment exponents after
Cauchy–Schwarz are \(3/5,2/5<2\). This application uses the
reduced count of Theorem 6.3 and Lemma 4.1. It needs neither the
separate Reuss box bound in Theorem 3.1 nor the anisotropic estimate
of Theorem 3.4 or the optimization in Proposition 7.1(ii).

### Exact optimization

Let \(\theta=12127/49920\). If any of \(s,t,a+\bar v\) is at most
\(\theta\), use that estimate. Otherwise set \(\Delta=b-a\ge0\).
The fiber inequality gives
\[
 a+\bar v=\tfrac12-\tfrac14s-\tfrac54\Delta>\theta,\qquad
 \Delta<(2-s-4\theta)/5 .
\]
Here \(6/25<\theta\) and \(s,t>\theta\), so
\(6/25\le s\le38/75\). On this interval \(m\ge(26/25-s)/5\):
the difference is concave, and it is positive at both endpoints
by substitution. Thus \(\Delta\le m\).
The first maximum in \(B\) is \((s-m)/2\) and the second is
\((t-m)/2+\tfrac12\max(0,3\Delta/2-m)\). Since \(3s+2t=2\),
\[
 B=\max\left(\tfrac1{10}+\tfrac34m,\;
              \tfrac1{10}+\tfrac{13}{20}m+\tfrac3{20}\Delta\right).
\]
The first term is at most \(77/320<\theta\), as \(m\le3/16\).
The second is at most
\[
 \tfrac4{25}-\tfrac3{25}\theta+\tfrac{561}{800}s-\tfrac{351}{320}s^2
 \ \le\ \tfrac4{25}-\tfrac3{25}\theta+
             \frac{(561/800)^2}{4(351/320)}=\theta .
\]
The last equality defines the displayed rational exponent. Also
\(m\le3/16<\theta\). There are \(O((\log X)^3)\) windows over
all dyadic sizes and kernel scales. Their count gives (1) by
taking a smaller exponent loss in each window and absorbing the
logarithms.

The fixed fiber and optimization are written proofs. The determinant
estimate is an imported theorem whose stated windows and proof
interfaces above were checked. See also the
[primary source audit](PLYMESSER_2026_POWERFUL_COUNTS_SOURCE_SCOPE.md).
The argument has no independent mathematical review or formal replay.
