# Consecutive factorial moments with one private linear variable

Let \(K\) be a field of characteristic zero, and let
\(\mathcal L(X^\nu)=\prod_i\nu_i!\). For two exponent vectors
\(\alpha,\beta\) put
\[
 D_\alpha=\sum_{\beta_i=0}\alpha_i,\qquad
 D_\beta=\sum_{\alpha_i=0}\beta_i.
\]
Thus \(D_\alpha\) measures the degree of \(X^\alpha\) on variables
absent from \(X^\beta\). Assume \(D_\alpha=1\) and \(D_\beta\ge2\),
possibly after interchanging the vectors.

**Theorem.** For every \(n\ge360\) and every nonzero
\(F=A X^\alpha+B X^\beta\), at least one of \(\mathcal L(F^n)\)
and \(\mathcal L(F^{n+1})\) is nonzero. The terminal threshold
\(N=n+1\ge361\) is independent of the number of variables, the
private and shared exponents, and the coefficient field.
The conclusion holds for every \(n\ge1\) when \(D_\beta\ge3\).
For \(D_\beta=2\) it also holds for every \(n\ge1\) if
\[
 S_{\mathrm{sh}}=\sum_{i:\alpha_i\beta_i>0}
        \frac{(\alpha_i-\beta_i)^2}{\max(\alpha_i,\beta_i)}
 \ \ge\ 4\log 2 .
\]

The argument uses angular bounds for both ends of the normalized
moment polynomial. The two ends have different coefficient scales
because exactly one private degree equals one. The
[nondivisibility proof](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md)
also covers terminals 2 through 4 and proves the all-order result
when both private degrees are at least two. For \((D_\alpha,D_\beta)=(1,2)\) with \(S_{\mathrm{sh}}<4\log 2\),
the threshold above leaves a finite range of terminals, with unbounded
exponents there.

## 1. Factorial coefficient bounds

For \(0\le k\le n\), put
\[
 C_{n,k}=\prod_i ((n-k)\alpha_i+k\beta_i)!,\qquad
 T_n(z)=\sum_{k=0}^n c_{n,k}z^k,\qquad
 c_{n,k}=\binom nk\frac{C_{n,k}}{C_{n,0}^{1-k/n}C_{n,n}^{k/n}}.
\]
Both endpoint coefficients of \(T_n\) are one. The moment polynomial
\(\sum_k\binom nk C_{n,k}t^k\), divided by \(C_{n,0}\), equals
\(T_n(t(C_{n,n}/C_{n,0})^{1/n})\).

Write \(f_{n,j}=j!/(n!)^{j/n}\). For each private exponent \(e\ge1\)
the sequence \((ej)!/(j!)^e\) is log-convex. Its consecutive ratio
\(\prod_{\ell=1}^e(e-(e-\ell)/(j+1))\) increases with \(j\).
The sequence starts at one, so its value at \(j\le n\) is at most
its value at \(n\) to the power \(j/n\). The normalized contribution
of the row is therefore at most \(f_{n,j}^e\), with \(j=n-k\) or \(k\)
according to which variable is private. Shared rows contribute at
most one, by log-convexity of \(j!\) along the arithmetic progression
in that row. Since \(f_{n,k}\le1\) and
\(\binom nk f_{n,n-k} f_{n,k}=1\), the private-degree assumption gives
\[
 0<c_{n,k}\le \binom nk f_{n,n-k}f_{n,k}^{\,2}=f_{n,k}.
 \tag{1}
\]

Put \(h=\lfloor n/2\rfloor\), split the coefficients at \(h\), and write
\[
 L_n(z)=\sum_{k=0}^{h}c_{n,k}z^k,\qquad
 U_n(w)=\sum_{j=0}^{n-h-1}c_{n,n-j}w^j.
\]
Thus \(T_n(z)=L_n(z)+z^n U_n(1/z)\) for \(z\ne0\). Both half-polynomials
have constant coefficient one and degree at most \(h\).

## 2. Bounds on the half-polynomials

We first prove the following bounds for any \(n\ge360\) and any
\(r>0\) with \(|\log r|\le3/n\). Let \(\ell_j,u_j\) be the respective
coefficients of \(L_n(rz),U_n(r^{-1}z)\), with zero coefficients outside
their degree ranges. Then
\[
 \begin{array}{lll}
 \ell_j\le (39/5000)(31/200)^{j-1}&(1\le j\le20),&
 \ell_j\le(13/25)^j\quad(j>20),\\
 u_j\le(19/50)^j&(1\le j\le9),&
 u_j\le(51/100)^j\quad(j>9).
 \end{array}
 \tag{2}
\]

The factorial bounds used below are
\[
 m\log m-m+1\le\log(m!)\le m\log m-m+1+\tfrac12\log m .
 \tag{3}
\]
The lower bound integrates \(\log x\) on the intervals immediately
before each integer. The upper bound follows because the trapezoidal
sum of the concave function \(\log x\) on \([1,m]\) is at most its
integral. We use \(e<11/4\), \(e^{-1}<46/125\), and \(\log360<6\).
The first follows by summing through \(1/3!\) and bounding the rest
by \((1/24)/(1-1/5)\). For the second, the sum through \(1/6!\)
exceeds \(125/46\). The third follows from \(e>27/10\) and
\((27/10)^6>360\).

For the left coefficients, use (1) and \(r^j\le e^{3j/n}\).
The ratio of \(e^{3j/n}f_{n,j}\) to the preceding term equals
\(j e^{3/n}/(n!)^{1/n}\), which by (3) and
\(e^x\le(1-x)^{-1}\) for \(0\le x<1\) is bounded by
\(11j/(4(n-3))\). For \(j=1\) this is below \(39/5000\), and for
\(j\le20\) it is below \(31/200\).
For any \(1\le j\le h\), monotonicity of \((j!)^{1/j}\) and (3) give
\[
 \bigl(e^{3j/n}f_{n,j}\bigr)^{1/j}
 \le \frac12\exp\!\left(\frac{1+\tfrac12\log h}{h}+\frac3n\right)
 \le \frac{180}{349}<\frac{13}{25}.
\]
Here \(h\ge180\), the function \((1+\tfrac12\log x)/x\) decreases on
\(x\ge1\), and \(\log180<6\).

For the right coefficients, use (1) at \(k=n-j\) and
\(r^{-j}\le e^{3j/n}\). Then
\[
 \begin{aligned}
 \frac1j\log(e^{3j/n}f_{n,n-j})
 &=\frac3n+\log\frac{(n!)^{1/n}}{n}
   +\frac1j\sum_{i=0}^{j-1}-\log(1-i/n)\\
 &\le -1+\frac7{360}
   +\frac1j\sum_{i=0}^{j-1}-\log(1-i/n).
 \end{aligned}
 \tag{4}
\]
We used that \((4+\tfrac12\log n)/n\) decreases and \(\log360<6\).
For \(1\le j\le9\), the sum divided by \(j\) is at most
\((j-1)/(2(n-j+1))\le1/88\), using
\(-\log(1-y)\le y/(1-y)\) for \(0\le y<1\).
Consequently
\[
 \bigl(e^{3j/n}f_{n,n-j}\bigr)^{1/j}
 <\frac{46/125}{1-7/360-1/88}<\frac{19}{50}.
\]
For \(1\le j\le h\), the left Riemann sum in (4) is at most the
average of \(-\log(1-x)\) on \([0,j/n]\), hence at most its average
\(1-\log2\) on \([0,1/2]\). The bound becomes
\(\tfrac12 e^{7/360}\le180/353<51/100\), proving (2).

The bounds at \(r=1\) locate every root within the specified range
of \(r\). To see this, first sum (2) to obtain
\[
 \sum_{j\ge1}\ell_j<\frac1{100},\qquad
 \sum_{j\ge1}j\ell_j<\frac3{250},\qquad
 \sum_{j\ge1}u_j<\frac{31}{50}.
 \tag{5}
\]
For the left sums, \((13/25)^5<1/25\) gives respectively
\[
 \frac{39/5000}{1-31/200}+\frac{13/12}{25^4}<\frac1{100},
 \qquad
 \frac{39/5000}{(1-31/200)^2}
 +\frac{13/12}{25^4}\left(21+\frac{13}{12}\right)<\frac3{250}.
\]
For the right sum, \((51/100)^5<7/200\) gives the tail sum beyond
\(j=9\) below \(1/400\), and \(19/31+1/400<31/50\).

If \(r=|z|\ge1\) for a root, then
\(r^n |U_n(1/z)|=|L_n(z)|\), with
\(|U_n(1/z)|>19/50\) and \(|L_n(z)|<(101/100)r^h\) by (5) at
radius one. Thus \(r^{n-h}<3\). If \(r\le1\) at a root, then
\(|L_n(z)|>99/100\) and \(r^n|U_n(1/z)|<(81/50)r^{n-h}\),
so \(r^{n-h}>1/2\). As \(n-h\ge n/2\) and \(e^{3/2}>3\),
both inequalities give \(|\log r|<3/n\).

## 3. Uniform phase bounds

Let \(0<\theta\le\pi\). For any \(r\) in the stated range the real
parts of \(L_n(r e^{i\theta})\) and \(U_n(r^{-1}e^{i\theta})\) are
positive by (5). Use their arguments in \((-\pi/2,\pi/2)\).
For the left half,
\[
 \bigl|\arg L_n(r e^{i\theta})\bigr|
 \le \frac{\theta\sum_{j\ge1}j\ell_j}{1-\sum_{j\ge1}\ell_j}
 <\frac{\theta}{50}.
 \tag{6}
\]
For the right half we prove
\[
 -\frac{\theta}{5}
 <\arg U_n(r^{-1}e^{i\theta})
 <\frac{3\theta}{4}.
 \tag{7}
\]

Write \(q=19/50\), \(Q=51/100\) and set
\[
 S=\sum_{j\ge2}q^j=\frac{q^2}{1-q},\qquad
 T=\sum_{j\ge2}j q^j=\frac{q}{(1-q)^2}-q .
\]
The contribution of \(j\ge10\) to the mass and weighted mass under
the bound \(Q^j\) in (2) is below \(1/400\) and \(7/250\),
respectively. The weighted bound follows from the mass bound by
multiplication by \(10+Q/(1-Q)=541/49\).

For the upper bound in (7), compare the terms in the imaginary part
after multiplying \(U_n(r^{-1}e^{i\theta})\) by \(e^{-3i\theta/4}\).
The constant term contributes \(-\sin(3\theta/4)\).
The \(j=1\) sine term is \(\sin(\theta/4)\le\sin(3\theta/4)\).
For \(j\ge2\), put \(x=3\theta/4\) and \(v=(j-3/4)/(3/4)\).
Then \(|\sin(vx)|\le v\sin x\). For \(x\le\pi/2\) with
\(vx\le\pi\), this is monotonicity of \(\sin x/x\) on \((0,\pi]\).
If \(x\le\pi/2\) and \(vx>\pi\), use
\(v\sin x\ge2vx/\pi>2\). For \(\pi/2\le x\le3\pi/4\) use
\(v\sin x\ge(5/3)\sin(3\pi/4)>1\).
The ratio of the nonconstant contributions to \(\sin(3\theta/4)\)
is thus at most
\[
 q+\frac43(T+7/250)-S<1 .
 \tag{8}
\]
The rotated imaginary part is strictly negative. The original argument
lies in \((-\pi/2,\pi/2)\), so after rotation it lies in
\((-5\pi/4,\pi/2)\). The negative imaginary part forces it to be
negative on this interval. This gives the upper bound.

For the lower bound in (7), multiply instead by \(e^{i\theta/5}\)
and compare against the constant term \(\sin(\theta/5)\).
Suppose first that \(\theta\le5\pi/6\). The \(j=1\) contribution is
nonnegative. The negative parts of the other sines are bounded by
\(2(j+1/5)\theta/9\). Indeed \(-\sin y\le2y/9\) for \(y\ge\pi\).
For \(y\ge9/2\) this follows from \(|\sin y|\le1\).
For \(\pi\le y\le9/2\), put \(u=y-\pi\). Since \(0\le u<7/5\),
the alternating Taylor bound gives
\(\sin u\le u-3u^3/20\). The maximum of
\(7u/9-3u^3/20\) on \(u\ge0\) is
\(14\sqrt{140}/243<2\pi/9\), proving the sine bound.
By concavity, \(\sin(\theta/5)\ge(\theta/\pi)\sin(\pi/5)\).
The Taylor lower bound at \(157/250<\pi/5\) gives
\(\sin(\pi/5)>29/50\). Thus the ratio of the total negative contribution
to the constant term is bounded by
\[
 \frac{2}{9}\,\frac{22}{7}\,\frac{50}{29}
 \left(T+\frac7{250}+\frac15(S+1/400)\right)<1 .
 \tag{9}
\]
We used \(157/50<\pi<22/7\). For \(5\pi/6\le\theta\le\pi\),
\(\sin(\theta/5)\ge1/2\), and
\(\sin(6\theta/5)\ge-\sin(\theta/5)\) by adding the two sines.
The remaining coefficients contribute at worst the negative of their
sum. In this range the ratio of the total negative contribution
to the constant term is bounded by
\[
 q+2(S+1/400)<1 .
 \tag{10}
\]
Equations (8)–(10) are rational inequalities. The positive imaginary
part after rotation by \(e^{i\theta/5}\) now implies the lower bound,
since the rotated argument is in \((-\pi/2,7\pi/10)\).

## 4. Consecutive separation

At a root of \(T_m\) with modulus \(r\) and argument \(\theta\),
the equation between the two halves gives
\[
 m\theta=(2j_m+1)\pi+\delta_m,\qquad
 \delta_m=\arg L_m(r e^{i\theta})+\arg U_m(r^{-1}e^{i\theta})
\]
for some integer \(j_m\). We used conjugation for the upper half.
By (6) and (7) the phase correction satisfies
\[
 -\frac{11}{50}\theta<\delta_m<\frac{77}{100}\theta.
\]
If consecutive moment polynomials have a common complex root, their
positive normalizations preserve its argument. The coefficients are
real, so after conjugation we may take \(0<\theta\le\pi\); a positive
real root is impossible because all coefficients are positive.
Both normalized radii satisfy the relevant bounds of Section 2.
Subtracting the two phase equations gives
\(\theta-\delta_{n+1}+\delta_n=2\pi(j_{n+1}-j_n)\).
The interval for the corrections places the left side strictly between
\(\theta/100\) and \(199\theta/100\), hence strictly between 0 and
\(2\pi\). This is impossible.

The unnormalized moment polynomials have rational coefficients, so
their complex coprimality gives a rational Bézout identity and hence
coprimality over \(K\). For \(A\ne0\) the moment of \(F\) at order \(m\)
is \(A^m\) times the moment polynomial evaluated at \(B/A\).
If either coefficient vanishes, all positive-order moments of the
remaining nonzero monomial are nonzero.

## 5. All orders when the other private degree is at least three

When \(D_\beta\ge3\), the proof of (1) gives the stronger bound
\[
                         c_{n,k}\le f_{n,k}^{\,2}.
 \tag{11}
\]
We use a finite certificate of angular inequalities for the individual
orders \(4\le m\le360\). Its inputs depend only on \(m\). It does not
enumerate possible exponent vectors.

### Coefficient envelopes

At an individual order \(m\), put \(h=\lfloor m/2\rfloor\),
\(g=m-h-1\). The checker takes
\[
 d^L_j\ge f_{m,j}^2 \quad(1\le j\le h),
 \qquad
 d^U_j\ge f_{m,m-j}^2 \quad(1\le j\le g),
\]
using rational bounds obtained by replacing \((m!)^{1/m}\) with
a lower rational and rounding up the resulting quotients. Both
masses \(D_L=\sum_j d^L_j\) and \(D_U=\sum_j d^U_j\) are below one.
At a normalized root of modulus \(r\), the same half-polynomial
comparison as in Section 2 gives
\[
 r\le\left(\frac{1+D_L}{1-D_U}\right)^{1/(m-h)},\qquad
 r^{-1}\le\left(\frac{1+D_U}{1-D_L}\right)^{1/(m-g)}.
 \tag{12}
\]
In the first bound use \(r\ge1\), and in the second use \(r\le1\);
in the opposite cases the respective inequality is automatic.
The checker rounds both displayed roots up to obtain rational
\(R,I\) and then rounds up \(d^L_j R^j,d^U_j I^j\). Denote these
last sequences by \(e^L,e^U\). They bound the coefficients of
\(L_m(rz),U_m(r^{-1}z)\) at the actual radius of any root. For each
sequence \(e\) the checker verifies \(E=\sum_j e_j<1\). Thus any
polynomial \(P(z)=1+\sum_{j\ge1}p_j z^j\) with \(0\le p_j\le e_j\)
has positive real part on the unit circle.

### Angles for orders 7 through 360

For \(0<\lambda\le1/2\), the upper bound
\(\arg P(e^{i\theta})<\lambda\theta\) on \(0<\theta\le\pi\)
follows from
\[
                  \sum_{j\ge1}\frac{j-\lambda}{\lambda}e_j<1 .
 \tag{13}
\]
Indeed the sine estimate in Section 3 with \(x=\lambda\theta\le\pi/2\)
and \(v=(j-\lambda)/\lambda\ge1\) bounds the rotated imaginary part.
The checker verifies (13) for \(\lambda=1/4\) on the left and
\(\lambda=43/100\) on the right.

For each half the checker also verifies
\[
 \frac{2(22/7)}{9(45/100)}
           \sum_{j\ge2}(j+3/20)e_j<1,\qquad
 e_1+\frac{\sum_{j\ge2}e_j}{398/1000}<1.
 \tag{14}
\]
These imply \(\arg P(e^{i\theta})>-3\theta/20\).
Multiply by \(e^{3i\theta/20}\). For
\(\theta\le20\pi/23\), the \(j=1\) term has nonnegative imaginary
part. Use \(-\sin y\le2y/9\) for \(y\ge\pi\) and the constant term
bound
\(\sin(3\theta/20)/\theta\ge\sin(3\pi/20)/\pi>(45/100)/\pi\).
For \(\theta\ge20\pi/23\), use
\(\sin(3\theta/20)\ge\sin(3\pi/23)>398/1000\) and
\(\sin(23\theta/20)\ge-\sin(3\theta/20)\).
The last inequality follows by adding the two sines.
The rational sine bounds follow from \(x-x^3/6\) at
\(x=471/1000\) and \(x=471/1150\), respectively.
Consequently the phase correction \(\delta_m\) in Section 4 lies
in \((-3\theta/10,17\theta/25)\) on this whole range of orders.
Its interval has width \(49\theta/50<\theta\).

### Angles for orders 4 through 7

For \(0<\theta\le\pi/2\), use the positive upper slopes from the
first two columns of the following table.
\[
\begin{array}{c|cc|cc}
 m & \lambda_L & \lambda_U
   & \text{lower allowance for }L_m & \text{lower allowance for }U_m\\ \hline
 4&.61&.29&.485&0\\
 5&.36&.45&.18&.27\\
 6&.44&.38&.20&.17\\
 7&.25&.43&.09&.20
\end{array}
\]
For \(\lambda\le1/2\) the checker uses (13). For \(\lambda=.61\),
the same inequality applies with the weight
\((1-\lambda)/\lambda\) for \(j=1\) divided by
\(1-(\lambda(22/7)/2)^2/6\). This follows from
\(\sin((1-\lambda)\theta)\le(1-\lambda)\theta\) and the cubic
lower bound on \(\sin(\lambda\theta)\).
All other weights in (13) are unchanged. The phase correction
is therefore less than \(.90\theta\).
For a lower bound, the \(j=1,2\) sine terms are nonnegative.
The possible negative contribution from the others is at most
\((2/9)\theta\sum_{j\ge3}j e_j\), by the same sine bound as before.
Divide by the real part, which is at least \(1-E\).
The checker verifies that the two quotients of
\((2/9)\sum_{j\ge3}j e_j\) by \(1-E\) add to less than \(.09\).
This gives \(\delta_m>-.09\theta\).

For \(\pi/2\le\theta\le\pi\), the last two columns in the table
bound the negative magnitudes of the respective phases, in radians.
To see how they are checked, the first-degree part \(1+p_1 e^{i\theta}\)
has nonnegative argument and modulus at least \(1-e_1\).
The rest has modulus at most \(\sum_{j\ge2}e_j\).
The ratio \(\sum_{j\ge2}e_j/(1-e_1)\) is below the cubic sine
lower bound for the table entry when the numerator is positive.
The disk about \(1+p_1 e^{i\theta}\) with radius
\(\sum_{j\ge2}e_j\) therefore allows at most the indicated downward
rotation. The phase correction is greater than \(-.485\).

For an upper bound in the same angular range, the \(j=2\) sine
contribution is nonpositive. At orders \(5,6,7\) the checker shows
that the two sums of \(e_1+\sum_{j\ge3}e_j\), each divided by its
\(1-E\), add to less than \(1.06\). These quotients bound the
positive arguments since \(\arctan x\le x\) for \(x\ge0\).
At order 4 the left half has degree two. The checker bounds
its quotient \(e_1/(1-E)\) by
\[
         \frac{.65-(.65)^3/6}{1-(.65)^2/2+(.65)^4/24}
         \le \tan(.65).
\]
The right half is linear, with \(e_1<.41-(.41)^3/6\).
The two arguments are therefore below \(.65\) and \(.41\).
This gives \(\delta_m<1.06\) at all four orders in the table.

The phase interval in the last angular range has width
\(1.545<\pi/2\le\theta\). In the first angular range it has
width \(.99\theta\). Thus consecutive phase corrections at orders
\(4\le n\le6\) differ by less than \(\theta\). For \(7\le n<360\)
this follows from the common interval for orders 7 through 360.
In each case the subtraction argument of Section 4 excludes a
common root. Terminals 2 through 4 follow from the
[low-terminal proof](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md).
Together with the uniform tail this proves the all-order assertion.

### Exact replay

Run the [certificate](verify_multivariate_one_private_linear.py)
with standard-library Python 3. It encloses each positive root in
(12) by integer bisection at scale \(10^8\), using comparison of
integer powers. It checks the envelope and angle inequalities for
357 individual orders in each coefficient regime, with the two kinds
of phase tests both performed at order 7. Use the option
\(\texttt{--shared-curvature}\) to run the checks in Section 6 alone. No numerical angle samples or approximation
to an exponent-vector range are used. The inequalities above show
why those finite tests suffice for the whole asserted support class.

## 6. All orders with private degrees one and two and large shared curvature

For a shared row, write its two exponents as \(a\le b\), with gap
\(d=b-a>0\), and reverse \(k\) if the smaller exponent belongs to
\(X^\beta\). Set \(F_k=\log((am+dk)!)\).
For every \(1\le k<m\), its second difference satisfies
\[
 F_{k+1}-2F_k+F_{k-1}
 =\sum_{j=1}^d\log\left(1+\frac{d}{am+(k-1)d+j}\right)
 \ \ge\ \frac{d^2}{bm}.
 \tag{15}
\]
In each summand, use \(\log(1+x)\ge x/(1+x)\) and
\(am+kd+j\le bm\). The latter upper bound holds because \(k\le m-1\)
and \(j\le d\). Subtracting a quadratic of second difference
\(d^2/(bm)\) from \(F_k\) leaves a convex sequence, so its values lie
below the chord of its two endpoints. Consequently the contribution
of this row to normalized coefficients is at most
\(\exp(-d^2 k(m-k)/(2bm))\). This estimate is unchanged when \(k\)
is replaced by \(m-k\).

Multiplying the shared-row estimates in the proof of (1) gives
\[
 c_{m,k}\le f_{m,k}\exp\left(-\frac{S_{\mathrm{sh}}\,k(m-k)}{2m}\right)
 \le f_{m,k}\,2^{-\min(k,m-k)}
 \qquad\text{if } S_{\mathrm{sh}}\ge4\log2 .
 \tag{16}
\]
Equal shared rows contribute one and have zero curvature. The last
inequality uses \(k(m-k)\ge(m/2)\min(k,m-k)\).

Use the endpoint envelopes and radius comparisons in Section 5
with (16) in place of (11). The checker verifies positivity of the
real parts of both half-polynomials at the actual scaled radii
by the same mass inequalities.

For \(7\le m\le360\), the checker uses (13) with \(\lambda_L=.31\)
and \(\lambda_U=.40\). For the negative angles replace
\(3/20\) in (14) by \(7/50\), and the sine lower bounds \(45/100\)
and \(398/1000\) by \(42/100\) and \(37/100\), respectively.
The same sine comparison applies because
\[
 \sin(7\pi/50)>42/100,\qquad
 \sin(7\pi/57)>37/100.
\]
These inequalities follow from the cubic sine lower bound at
\(1099/2500\) and \(1099/2850\), respectively. Thus
\(-.28\theta<\delta_m<.71\theta\), an interval of width
\(.99\theta\).

For \(4\le m\le7\), the following table gives the positive phase
slopes for \(0<\theta\le\pi/2\) and the absolute allowances
for negative angles in the range \(\pi/2\le\theta\le\pi\).
\[
\begin{array}{c|cc|cc}
 m &\lambda_L &\lambda_U
   &\text{lower allowance for }L_m &\text{lower allowance for }U_m\\ \hline
 4&.46&.26&.28&0\\
 5&.35&.40&.16&.19\\
 6&.40&.36&.18&.15\\
 7&.31&.40&.11&.18
\end{array}
\]
The checker evaluates (13) for each tabulated positive slope. In
the smaller-angle range the negative slope bound obtained from
the \(j\ge3\) sine terms, as in Section 5, is less than \(.09\)
for the two halves together. Hence
\(-.09\theta<\delta_m<.76\theta\) in this range.
In the larger-angle range the checker verifies each downward
rotation allowance by the disk comparison in Section 5. It also
bounds the sum of \((e_1+\sum_{j\ge3}e_j)/(1-E)\) over both halves
by \(1.17\). Thus \(-.35<\delta_m<1.17\), an interval of width
\(1.52<\pi/2\le\theta\).

The phase subtraction in Section 4 therefore excludes common
roots at lower orders \(4\le n<360\). The uniform tail covers
\(n\ge360\), and the low-terminal theorem covers \(n\le3\).
The replay certifies all 357 individual orders for the envelope (16).
