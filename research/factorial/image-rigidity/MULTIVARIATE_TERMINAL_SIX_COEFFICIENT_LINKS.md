# Coefficient constraints at terminal six

At terminal six, the order-five and order-six angular envelopes must
keep track of how both coefficient arrays concentrate in the same
exponent rows. We derive constraints on that concentration, and
exhibit a shared root in the envelope that keeps the within-order
constraints and discards the cross-order ones. We also show that the
within-order constraints suffice near the flat limit.

Use the endpoint normalization from the
[terminal-five proof](MULTIVARIATE_TERMINAL_FIVE_ANGULAR_SEPARATION.md).
The binomial multipliers are
\[
 (1,u e^\alpha,v e^\chi,v e^{-\chi},u e^{-\alpha},1),\qquad
 (1,w e^\beta,y e^\eta,z,y e^{-\eta},w e^{-\beta},1).
\]
In one Gauss factor let \(x=1/T\), \(g_j=\log(1+jx)\), and
\(L_k=\sum_{j=0}^{k-1}g_j\). The normalized log multiplier at order
\(n\) and index \(k\) is \(L_k-(k/n)L_n\). The corresponding
\(x,x'\) at orders five and six satisfy \(5x/6\le x'\le x\).

## Within-order bounds

The order-five geometric bounds and skew link are as in the
terminal-five proof. For order six,
\[
\begin{gathered}
 u^{5/4}\le w\le\sqrt{\frac{24u^2}{25-u^2}},\qquad
 w^{8/5}\le y\le w\sqrt{\frac{1+4w^2}{4+w^2}},\\
 w^{1/5}\le z/y\le\sqrt{\frac{2+3w^2}{3+2w^2}}.
\end{gathered}                                                     \tag{1}
\]
The local squared outer mean is \(1/(1+5x)\). The squared
ratios \(y/w\) and \(z/y\) locally are
\((1+x)/(1+4x)\) and \((1+2x)/(1+3x)\).
They give (1) by the power comparison and product inequality for
the fractional linear functions in the terminal-five proof.
For the cross-order outer bounds, the ratio \(5x'/(4x)\)
lies in \([25/24,5/4]\). This gives the lower power bound and
the upper fractional linear bound by the same product argument.

The signed local skew and symmetric defect satisfy
\[
 4\beta-5\eta=\tfrac12(\Delta^4g_0+\Delta^4g_1),\qquad
 J=-\log z+\tfrac32\log y-\tfrac35\log w=\tfrac1{20}\Delta^5g_0.
                                                                    \tag{2}
\]
These identities follow by substituting \(L_k-(k/6)L_6\),
averaging reflected entries for the means, and subtracting them
for the skews. Reversing the row changes the sign of the first
expression and leaves the second unchanged. The fourth and fifth
derivatives of \(g(t)=\log(1+xt)\) give
\[
 |4\beta-5\eta|\le6\bigl((w^{-2}-1)/5\bigr)^4,\qquad
 0\le J\le\tfrac65\bigl((w^{-2}-1)/5\bigr)^5.                         \tag{3}
\]
To sum the single-factor bounds, use
\(\sum x_i\le(\prod_i(1+5x_i)-1)/5\) and
\(\sum x_i^j\le(\sum x_i)^j\).

Both bounds also admit single-factor concentration caps. For
\(q=1+5x\) let \(D=(q+4)^2(2q+3)^2(3q+2)^2(4q+1)^2\).
The second derivatives in \(\log q\) of the two nonnegative local
magnitudes \(|4\beta-5\eta|\) and \(J\) are
\[
 \frac{24q(q-1)^2(101q^4+471q^3+731q^2+471q+101)}{D},
 \qquad
 \frac{300q(q-1)^3(q+1)(q^2+3q+1)}{D}.
\]
Both functions vanish at zero in the \(\log q\) variable and are
convex there. Superadditivity therefore bounds their sums by
the respective single-factor values at \(q=w^{-2}\).

## An exact collision in the separate-order envelope

Set \(u=3/40\), \(v=u^{3/2}\), and
\(w=\sqrt{24u^2/(25-u^2)}\). Use zero skews and vary the
sextic means from
\[
 (y_0,z_0)=(w^{8/5},w^{9/5})\quad\text{to}\quad
 y_1=w\sqrt{(1+4w^2)/(4+w^2)},\qquad
 z_1=y_1\sqrt{(2+3w^2)/(3+2w^2)}.
\]
Both endpoints satisfy (1), the concentration cap on \(J\),
and strict log convexity. For the first endpoint \(J=0\),
and for the second it equals its single-factor cap. Linear
interpolation of \((\log y,\log z)\) preserves all these
conditions, including strict log convexity.

For symmetric arrays the quotient polynomials in \(X=t+t^{-1}\),
after removing \(t+1\) from the quintic, are
\[
 B(X)=X^2+(5u-1)X+10v-5u-1,\qquad
 C(X)=X^3+6wX^2+(15y-3)X+20z-12w .
\]
The larger root \(X_*\) of \(B\) lies in \((1.438,1.439)\).
At the two sextic endpoints,
\[
 C_{y_0,z_0}(X_*)\in(-.801,-.789),\qquad
 C_{y_1,z_1}(X_*)\in(.094,.106).
\]
For direct rational enclosures use \(w\in(.07349,.07350)\),
\(v\in(.02053959,.02053960)\),
\(y_0\in(.01534,.01535)\), \(z_0\in(.00910,.00911)\),
\(y_1\in(.03711,.03712)\), and \(z_1\in(.03037,.03038)\).
These follow by squaring or taking fifth powers. The two signs
give a common complex unit-circle root at an intermediate point.
This is a collision in the relaxed coefficient envelope.

## Cross-order transfers

Define \(G_a(s)=as-\operatorname{artanh}(a\tanh s)\) for
\(0<a<1\) and \(s\ge0\). If \(S=-\log u,\ W=-\log w\), set
\[
 V=\log(v/u^{3/2}),\quad
 Y=\log(y/w^{8/5}),\quad Z=\log(z/(y w^{1/5})).
\]
For the Gauss factors write their contributions to \(S,W\) as
\(s_i,s'_i\). The three defects are sums of \(G_{1/2}(s_i)\),
\(G_{3/5}(s_i')\), and \(G_{1/5}(s_i')\), respectively.
Then
\[
 \tfrac{128}{125}V\le Y\le2V,\qquad
 \tfrac25 V\le Z\le V.                                                \tag{4}
\]
For the lower inequalities use \(s_i'\ge s_i\) and compare the
derivatives \(G_a'(s)=a(1-a^2)\tanh^2 s/(1-a^2\tanh^2 s)\).
For the upper bounds let \(t=\tanh s_i\). The relation
\[
 \tanh s_i'=\frac{\lambda t}{1+(\lambda-1)t},
 \qquad 25/24\le\lambda\le5/4
\]
allows the upper comparison at \(\lambda=5/4\).
The ratio of the derivative in \(t\) of
\(G_a(s_i')\) to that of \(G_{1/2}(s_i)\) is
\(2N/D_a\) for \(a=3/5\) and \(N/D_a\) for \(a=1/5\), where
\[
 N=(1+t)(1-t^2/4),\qquad
 D_a=((1+t/4)^2-a^2(5t/4)^2)(1+3t/2).
\]
In the two cases \(D_a-N\) equals
\(t(2-t)(1+t)/2\) and \(t(t+2)^2/4\), respectively.
Integrate from zero to obtain (4). In the collision above \(V=0\)
and \(Y>0\), so it violates this transfer.

Write \(D_5=G_{1/2}(S)-V\). The concentration deficits also obey
\[
 Y\le G_{3/5}(W)-\tfrac{128}{125}D_5,\qquad
 Z\le G_{1/5}(W)-\tfrac5{16}D_5.                                    \tag{5}
\]
In general, for \(s'_i\ge s_i\ge0\), the concentration deficit
of \(G_b\) at the primed list is at least the concentration
deficit of \(G_a\) at the unprimed list times
\[
 \frac{b(1-b^2)}{a(1-a^2)}
       \min\left(1,\left(\frac{1-a^2}{1-b^2}\right)^2\right).
\]
Indeed, for a convex \(G\) vanishing at zero,
\(G(\sum s_i)-\sum G(s_i)\) is increasing in each \(s_i\).
The displayed constant is a lower bound on \(G_b''/G_a''\),
as one sees from
\(G_a''(s)=2a(1-a^2)t(1-t^2)/(1-a^2t^2)^2\),
where \(t=\tanh s\). Subtract the corresponding multiple of
\(G_a\) and use the superadditivity of the resulting convex
function. The mean ratios for adjacent indices \(k-1,k\)
at order \(n\) use \(a=(n+1-2k)/(n-1)\) when \(0<a<1\).

One can also require
\[
 W-\tfrac12\log(1+\tfrac{25}{24}(e^{2S}-1))
 \ge \tfrac3{200}D_5.                                                 \tag{6}
\]
Let \(f(s)=\tfrac12\log(1+\tfrac{25}{24}(e^{2s}-1))\).
Then \(s_i'\ge f(s_i)\), and \(-f''(s)\ge(2/25)e^{-2s}\),
while \(G_{1/2}''(s)\le(16/3)e^{-2s}\).
The convex function \(-f-(3/200)G_{1/2}\) is superadditive,
giving (6).

The skew magnitudes also obey a same-order concentration
constraint. Write \(B_a(s)\) for the magnitude of the
single-factor outer skew at the order under consideration.
These functions are increasing and convex in \(s\), from
their positive-coefficient sums of \(g_j''\) in the variable
\(\log(1+(n-1)x)=2s\). If \(m=\max_i s_i\), \(T=\sum_i s_i\)
and \(H=\sum_iG_a(s_i)\), then
\[
 H\le T\,G_a(m)/m,\qquad
 |\text{outer skew}|\ge \max(0,B_a(m)-B_a(T-m)).
\]
The first inequality uses monotonicity of \(G_a(s)/s\).
The second follows by summing the magnitudes of every
other factor and applying superadditivity. In particular,
equality \(H=G_a(T)\) forces one nonzero factor and its
associated skew magnitude. Sampling with (4)–(5) alone
still shows overlapping angular bands near \(u=.2\),
where an artificially zero sextic skew is paired with a
saturated sextic geometric bound. The same-order constraint
excludes that choice. This diagnostic does not certify
angular separation on the remaining domain.

## Near-flat angular separation

Let \(u=1-\varepsilon\), \(w=1-r\varepsilon\), and
rescale \(t=-1+\sqrt{\varepsilon}\,\xi\). By (1), any limiting
\(r\) belongs to \([25/24,5/4]\). The terminal-five skew link
gives the order-five limit
\(\xi^5+5\xi^3+(15/4)\xi\) after division by \(\varepsilon^{5/2}\).
At order six the geometric bounds give
\(y=w^{8/5}+O(\varepsilon^3)\),
\(z=w^{9/5}+O(\varepsilon^3)\), and (2)–(3) give
\[
 z-w^{9/5}=\tfrac32(y-w^{8/5})+O(\varepsilon^4).
\]
The combined central perturbation therefore has the form
\(15(y-w^{8/5})t^2(1+t)^2+O(\varepsilon^4)\) on bounded
\(\xi\)-sets. The linear skew terms combine using
\((t-t^5)+2(t^2-t^4)=t(1-t)(1+t)^3\). The discrepancy in
their coefficients is \(O(\varepsilon^3)\), since
\(4\beta-5\eta=O(\varepsilon^4)\) and \(w-y=O(\varepsilon)\).
Each individual skew is \(O(\varepsilon^2)\) uniformly, since its
constant and linear terms in each local \(x_i\) cancel and
\(\sum_i x_i^2\le(\sum_i x_i)^2=O(\varepsilon^2)\).
The remaining even skew errors are \(O(\varepsilon^4)\).
Consequently the order-six limit after division by
\(\varepsilon^3\) is
\[
 \xi^6+6r\xi^4+\tfrac{36}{5}r^2\xi^2+\tfrac{24}{25}r^3.
\]
Here the limiting symmetric coefficients follow by expanding
\(w,w^{8/5},w^{9/5}\) through degree three in \(\varepsilon\).
The positive imaginary root magnitudes squared are \(r\) times
the three roots of
\(X^3-6X^2+(36/5)X-24/25\), isolated respectively in
\((.152,.153),(1.427,1.428),(4.420,4.421)\) by signs at the
rational endpoints. These ranges interlace the order-five
values \((5-\sqrt{10})/2,(5+\sqrt{10})/2\) throughout the
allowed \(r\)-interval. Coefficient convergence on bounded
sets and simplicity of the limiting roots then give uniform
angular separation for \(u\) sufficiently close to one.
The [complete terminal-six proof](MULTIVARIATE_TERMINAL_SIX_ALL_PROFILES.md)
assembles the subsequent linked certificates for the intermediate range.
