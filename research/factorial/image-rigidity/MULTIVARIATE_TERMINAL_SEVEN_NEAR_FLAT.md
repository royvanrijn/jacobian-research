# Shared terminal-seven separation when the outer mean is at least .95

Let \(a,b\) be distinct ordinary exponent columns whose unequal rows
are all shared. Let \(u\) be the order-six outer geometric mean after
normalizing the constant and leading coefficients to one and removing
the binomial factors. If \(u\ge19/20\), the order-six and order-seven
polynomials have disjoint root arguments up to independent positive
rescalings. The dimension and exponents are arbitrary.

The [exact checker](verify_multivariate_terminal_seven_near_flat.py)
uses outward dyadic interval arithmetic on the Cayley variable. It
can be replayed from the repository root with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_terminal_seven_near_flat.py
~~~

## Scaling and the order ratio

Put \(p=1-u=\lambda^2\), \(q=1-w=rp\), where \(w\) is the
order-seven outer mean. Use \(0\le\lambda\le9/40\),
\(513/500\le r\le129/125\), and \(q<.055\). Every actual profile
in the stated interval has \(\lambda\) in the specified range.
The [linked outer-mean comparison](MULTIVARIATE_LINKED_GAMMA_ANGLE_COMPACTIFICATION.md)
gives \(w\le u(1+(1-u^2)/35)^{-1/2}\), whence
\[
 r\ge1+\frac{u(1+u)}{70(1+(1-u^2)/35)^{3/2}}.
\]
A Gauss reciprocal \(x=d/(6\ell+h)\) is at most
\((u^{-2}-1)/5<.022\). The corresponding order-seven reciprocal
\(x'=d/(7\ell+h)\) satisfies
\[
 \frac{6x'}{5x}=\frac{36}{5(7-xh/d)}\le\frac{36}{5(7-.022)} .
\]
The power comparison for the normalized outer means and the tangent
inequality for \(u^\gamma\), \(\gamma\ge1\), give the upper bound on
\(r\). These endpoint inequalities are checked exactly.

Use \(t=(\lambda X-1)/(\lambda X+1)\) and multiply the normalized
order-\(j\) polynomial by \((1+\lambda X)^j/(2^j\lambda^j)\).
The order-six enclosure is already proved in the
[near-flat terminal-six extension](MULTIVARIATE_TERMINAL_SIX_NEAR_FLAT_EXTENSION.md).
Here use its scaled sextic formula with its ratio argument equal to
one, since \(p=1-u\) is now the sextic deficit. Its within-order
Cayley coefficient bounds require only the sextic deficit cap and
hold with this substitution.

## Linked order-seven coefficients

Write the normalized septic multipliers as
\[
 (1,we^\alpha,ye^\beta,ze^\gamma,ze^{-\gamma},
                         ye^{-\beta},we^{-\alpha},1).
\]
In a single order-seven factor put \(g_h=\log(1+hx')\) and
\(s=g_6/2\). Its normalized log coefficient at \(k\) is
\(\sum_{h=0}^{k-1}g_h-(k/7)\sum_{h=0}^{6}g_h\). Subtracting
the reflected coefficients gives the signed links
\[
 L=\alpha-\beta=\frac{2\Delta^4g_0+3\Delta^4g_1+2\Delta^4g_2}{14},
 \qquad K=\alpha-3\beta+5\gamma=\frac{\Delta^6g_0}{14}.
\]
Thus \(|L|\le3(x')^4\) and \(|K|\le(60/7)(x')^6\)
for one factor, by integrating the derivatives of
\(\log(1+x't)\) over unit cubes. The local outer skew is at most
\((5/2)(x')^2\): pair the logarithms at \(i,6-i\), using
\(0\le g_i+g_{6-i}-g_6\le i(6-i)(x')^2\).
For the local \(\gamma\), use the three centered second differences
of steps one, two and three, whose absolute bounds sum to
\(14(x')^2\), and divide by fourteen. Reversing the orientation
negates all the signed links.

The sum of the septic reciprocals is at most
\((w^{-2}-1)/6\le.365q\), by expanding their product.
Summing the single-factor estimates gives
\[
 |\alpha|\le.334q^2,\quad |\beta|\le.335q^2,\quad
 |\gamma|\le.134q^2,\qquad |L|\le.054q^4,\quad |K|\le.021q^6 .
\]
All displayed decimals here and below are exact rationals.

For the means let \(G_a(s)=as-\operatorname{artanh}(a\tanh s)\).
Averaging reflected log coefficients gives
\[
 H=\log(y/w^{5/3})=\sum G_{2/3}(s),\qquad
 \log(z/w^2)=\tfrac95 H+J,\quad
 J=\sum\big(G_{1/3}(s)-\tfrac45G_{2/3}(s)\big).
\]
The derivative
\(G'_a(s)=a(1-a^2)\tanh^2(s)/(1-a^2\tanh^2(s))\)
gives \(0\le-J\le\sum s^5/25\). For \(H\) the same formula,
with \(\tanh s\le s\), yields
\[
 H\le \frac{10}{81}
       \frac{(-\log w)^3}{1-(4/9)(-\log w)^2}.
\]
Use \(-\log w\le q/(1-q)\) and \(e^h-1\le h/(1-h)\) for \(0\le h<1\).
It follows that \(y-w^{5/3}=Dq^3\) with \(0\le D\le.15\).
For the central coefficient expand \(w^{-1}y^{9/5}\) at
\(y=w^{5/3}\), and then multiply by \(e^J\). The first derivative
at this point is \((9/5)w^{1/3}\); the Taylor error is at most
\((18/25) w^{-4/3}(y-w^{5/3})^2\). The exponential reduces the
value by at most \(w^{-1}|J|\) since \(y\le1\).
Therefore
\[
 z=w^2+\tfrac95Dq^3+E_{\rm geom}q^4,\qquad
              -.1\le E_{\rm geom}\le.0001.
\]
The checker verifies the rational caps used in this calculation,
including \((1-w^{1/3})/q\le.35\) and
\((18/25)w^{-4/3}\le.8\).

## Stable Cayley expansion

Set \(a=w\cosh\alpha\), \(b=y\cosh\beta\), \(c=z\cosh\gamma\)
and \(A=w\sinh\alpha\), \(B=y\sinh\beta\), \(C=z\sinh\gamma\).
Before division by \(2^7\lambda^7\), the Cayley coefficients in
\(Z=\lambda X\) are obtained by expanding
\(\sum_k \binom7k m_k(Z-1)^k(Z+1)^{7-k}\), where \(m_k\) are
the normalized multipliers. Their values are
\[
\begin{array}{c|r}
7&70c+14a+42b+2\\
6&70A+126B+70C\\
5&-210c+126a+42b+42\\
4&70A-210B-210C\\
3&210c-70a-210b+70\\
2&-126A+42B+210C\\
1&-70c-70a+126b+14\\
0&-14A+42B-70C .
\end{array}
\]
For the same-parity coefficients substitute
\[
\begin{aligned}
 a&=1-q+A_0q^4,&0\le A_0&\le .056,\\
 b&=1-\tfrac53q+\tfrac59q^2+(R+D)q^3+B_0q^4,&0\le B_0&\le .057,\\
 c&=1-2q+q^2+\tfrac95Dq^3+Eq^4,&-.1\le E&\le .0091 ,
\end{aligned}
\]
where \(5/81\le R\le5/81+.022q\) is the cubic Taylor term
and tail of \(w^{5/3}\). The hyperbolic corrections follow from
\(\cosh h-1\le.501h^2\) on \(|h|\le.01\). The checker
performs the cancellations in these formulas before dividing
by the powers of \(\lambda\).

For the opposite-parity coefficients use
\(\alpha=(5\gamma+3L-K)/2\),
\(\beta=(5\gamma+L-K)/2\). Retain the common parameter
\(\gamma/q^2\) in the four degrees. In degree two its contribution
in the unscaled numerator is
\(-105((w-y)+2(w-z))\gamma\); in degree zero it is
\(-35(w-3y+2z)\gamma\). Use the stable expansions
\[
 \frac{w-y}{q}=\tfrac23-\tfrac59q-(R+D)q^2,\qquad
 \frac{w-z}{q}=1-q-\tfrac95Dq^2-E_{\rm geom}q^3,
\]
\[
 \frac{w-3y+2z}{q^2}
     =\tfrac13+(-3R+\tfrac35D)q+2E_{\rm geom}q^2 .
\]
The last expression cancels both the constant and linear terms
of the original numerator. The remaining bounds in degrees
\(6,4,2,0\) are respectively
\[
 .0003 r^2\lambda^3,\quad .00002r^2\lambda,\quad
 .0045r^3\lambda,\quad .0004r^4\lambda .
\]
They use the bounds for \(L,K\), \(0\le w-y\le 2q/3\), and
\(|\sinh h-h|\le |h|^3/5\) on \(|h|\le.01\).
For example the degree-zero \(L\)-contribution has numerator
\(-21(w-y)L\), and the \(K\)-contribution has numerator
\((7w-21y)K\). The scalar checks cover all four remainders.

## Rays and sector counts

Leading-term dominance bounds all scaled roots by \(|X|<2\).
Since \(\lambda\le9/40\), all physical roots then lie in the
open left half-plane. For each positive ray parameter \(b\),
the physical argument is \(\pi-2\arctan(\lambda b)\). A point
\(X=x+iY\) on it satisfies
\[
 (Y-b)(1+\lambda^2 bY)=-b\lambda^2x^2,\qquad 0<Y\le b .
\]
The checker takes three interval updates from \(0\le Y\le b\)
and tests signed projections
\(\Im(P_j(X)\overline X^{\,h})\) for \(0\le h\le j\).
It suffices to use \(x\ge0\), since reflection
\(X\mapsto-\overline X\) reverses the opposite-parity
coefficients and their joint enclosure is symmetric.
The five ray parameters are affine in \(p\), with the rational
coefficients in the checker. The exact replay used
\(4045,415,1961,541,175\) boxes for the sextic and
\(7599,1401,3191,847,657\) for the septic.

To reach the anchor, multiply all matched Gauss reciprocals
by a common parameter tending to zero. Their cross-order
factor ratios are fixed, so the lower and upper bounds on
\(r\) remain valid. The scaled limits as \(\lambda\to0\)
are the imaginary Hermite forms
\[
 X^6+\tfrac32X^4+\tfrac9{20}X^2+\tfrac3{200},\qquad
 X^7+\tfrac74rX^5+\tfrac{35}{48}r^2X^3+\tfrac{35}{576}r^3X .
\]
Exact rational brackets for the positive squared imaginary
coordinates show that the successive upper sectors contain
orders \(7,6,7,6,7,6\) throughout the permitted range of \(r\).
The checker verifies the brackets by evaluating the two
cubics in the squared coordinate, dividing that coordinate
by \(r\) for order seven. At the anchor the ray equations
reduce to \(Y=b\), so the verified exclusions preserve the
sector counts through the deformation.

The final joined sector has one septic root, which is real
by conjugacy, and two sextic roots. The actual sextic is
strictly positive on the real line, since it is the product
exponential integral of the sixth power of a nonzero real
binomial in two distinct monomials. Its two roots are nonreal
and consequently have distinct arguments from the septic root.
This proves disjointness, and the rational Bézout identity
gives the assertion over every characteristic-zero field.
