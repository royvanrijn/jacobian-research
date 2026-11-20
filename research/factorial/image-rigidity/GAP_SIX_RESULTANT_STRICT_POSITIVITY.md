# Strict positivity of the ordinary gap-six resultant

**Theorem OrdinaryGapSixStrictResultant.** For every \(n\ge1\), the normalized lower-one, gap-six resultant \(Q_{n,1,6}(1)\) is strictly positive. It is still strictly positive after removal of the additional forced divisor \(E_{n,6}(1)\). Consequently every nonzero ordinary factorial binomial supported on \(\{1,7\}\) has a nonzero entry in every consecutive pair of positive moments. Together with the other supported pairs in the [degree-seven reduction](ORDINARY_DEGREE_SEVEN_GAP_SIX_FRONTIER.tex), this proves ordinary univariate binomial Strong Factorial through degree seven.

## The two-shape resultant and the forests

Use the rising factorial \((\gamma)_k\) and the notation \(\mathcal R_{n,1,d},D_{n,1,d},E_{n,d}\) of the [shape-resultant normalization](BINOMIAL_SHAPE_RESULTANT.tex). Write \(N=n+1\), \(\gamma=\beta+N-1>0\), and set
\[
 A_j(z;\gamma)=\sum_{k=0}^j\binom jk(\gamma)_{dk}z^{j-k}.
\]
The lower reciprocal is \(A_{N-1}=A_N'/N\), with the same shape \(\gamma\). The upper reciprocal uses shape \(\gamma+1\), so
\[
 R_N(z;\beta)=
 \frac{\gamma+dN}{\gamma} A_N(z;\gamma)-\frac d\gamma zA_N'(z;\gamma).
\]
These follow by differentiating and using \((\gamma+1)_{dk}=(1+dk/\gamma)(\gamma)_{dk}\). Evaluation at the roots of the monic lower reciprocal gives the exact formula
\[
 \frac{\mathcal R_{N-1,1,d}(\beta)}{D_{N-1,1,d}(\beta)E_{N-1,d}(\beta)}
 =\frac{\operatorname{Res}_z(A_{N-1},A_N)}
        {\prod_{j=1}^{N-1}(\gamma)_{dj}}.                    \tag{1}
\]
Indeed \(D_{N-1,1,d}=\prod_{j=1}^{N-1}(\gamma+1)_{dj-1}\), and \(E_{N-1,d}=(\gamma+dN)^{N-1}\). The formula makes the two-shape test a fixed-shape squarefreeness test.

Here the existing [Gamma-power forests](lower_exponent_initial_and_boundary.tex) supply an exact signed determinant representation. Write \(q_k=\kappa_{k+1}(\gamma)/k!=(\gamma)_d u_k(\gamma)\), where \(u_k\) is the forest polynomial given there. On indices \(0\le i,j<m\), put \(z+q_0\) on the diagonal, \(-(i+1)\) on the superdiagonal, \(q_{i-j}\) below the diagonal, and zero elsewhere. The determinant of this matrix is \(A_m(z;\gamma)\). Expansion along the last row gives
\[
 A_m=(z+q_0)A_{m-1}+\sum_{j=0}^{m-2}\frac{(m-1)!}{j!}\,q_{m-1-j}A_j,\qquad A_0=1.
\]
The same recurrence follows from logarithmic differentiation of the formal series \(e^{zx}\sum_{k\ge0}(\gamma)_{dk}x^k/k!\), whose coefficient of \(x^m\) is \(A_m/m!\). In (1), the Sylvester determinant of these consecutive blocks uses only the forests with the common shape \(\gamma\). The proof of noncancellation below uses \(\gamma=N\) and \(d=6\).

## Uniform simplicity for terminal order \(N\ge10{,}000\)

The Gamma integral, with \(x=Ns\), gives
\[
 A_N(N^6w;N)=\frac{N^{7N}}{(N-1)!}\,I_N(w),\qquad
 I_N(w)=\int_0^\infty [s(w+s^6)]^N e^{-Ns}\frac{ds}{s}.
\]
The integrand at \(s=0\) is interpreted as \(s^{N-1}(w+s^6)^N e^{-Ns}\). Let \(K_N\) and \(J_N\) integrate on \([0,3]\) and \([3,\infty)\), respectively. Both are polynomials in \(w\), and \(I_N=K_N+J_N\). This split follows the contour analysis of the [constant-plus-sixth-power family](EFFECTIVE_EXPONENT_SIX_SIMPLICITY.tex). The physical shape in this problem gives a second saddle near \(s=1\).

### The two saddle approximations

For \(J_N\) use the phase
\[
 \Phi(s,w)=7\log s+\log(1+w/s^6)-s.
\]
For \(w^{-N}K_N(w)\) use
\[
 \Psi(s,w)=\log s+\log(1+s^6/w)-s.
\]
In either case the saddle equation is \(w(s-1)=s^6(7-s)\). Write \(s_J,s_K\) for the following solutions.

- For \(|w|\le3180\), Rouché's theorem on \(|s-7|=1/5\) gives exactly one solution there, counted with multiplicity. The inequality is \((1/5)(34/5)^6/(31/5)>3180\).
- For \(|w|\ge1000\), write the equation as \(s-1=w^{-1}s^6(7-s)\). On \(|s-1|=1/100\) the right side has modulus at most \((101/100)^6(601/100)/1000<1/100\). Thus there is exactly one solution in that disk. This solution is analytic in \(1/w\), including \(1/w=0\).

The logarithms at the saddles are the branches continued from \(w=0\) for \(\Phi\), and \(1/w=0\) for \(\Psi\). They can also be taken by the convergent local logarithmic series. At either saddle the negative second derivative is \(a(s)=(s^2-7s+7)/s^2\). Denote the real parts of the two saddle heights by \(S_J,S_K\). The elementary logarithmic bound \(|\log(1+v)-v|\le |v|^2/(2(1-|v|))\) gives
\[
 S_J>7\log7-7-\tfrac1{25}>6.58,\qquad S_K>-1-\tfrac1{500}.
\]
For the first bound, at \(s_J=7+v\) the height difference from \(7\log7-7\) is \(7\log(1+v/7)-\log(1+v/6)-v\). For the second use \(|s_K-1|<1/100\) and \(|s_K^6/w|\le(101/100)^6/1000\). The curvature satisfies
\[
 \Re a(s_J)>\tfrac3{25},\quad |a(s_J)|<\tfrac15,
 \qquad \Re a(s_K)>\tfrac9{10},\quad |a(s_K)|<\tfrac{11}{10}.
\]
For example \(a(s)-1/7=(s-7)(6s-7)/(7s^2)\) for the first disk and \(a(s)-1=-7(s-1)/s^2\) for the second.

The saddle estimates on their full stated domains are
\[
\begin{split}
 J_N(w)&=\frac{e^{N\Phi(s_J,w)}}{s_J}\sqrt{\frac{2\pi}{N a(s_J)}}(1+\epsilon_J(w)),\\
 K_N(w)&=w^N\frac{e^{N\Psi(s_K,w)}}{s_K}\sqrt{\frac{2\pi}{N a(s_K)}}(1+\epsilon_K(w)),
 \qquad |\epsilon_J|,|\epsilon_K|<\frac{25}{\sqrt N}.       \tag{2}
\end{split}
\]
The square roots use \(\Re a>0\). The proof retains all portions of the deformed contours.

For the local integral at \(s_*\), integrate along \(s_*+t\) for \(|t|\le \rho\). The choices for \(J\) are \((\rho,B,c)=(3/4,1/30,7/200)\), and for \(K\) they are \((1/5,5/6,1/4)\). Here \(B\) bounds the coefficient in the phase Taylor remainder \(B|t|^3\), and \(c\le\Re a/2-B\rho\). To check \(B\), let \(r\) be an upper bound on \(|w/s^6|\) or \(|s^6/w|\) along the local disk and let \(y=r/(1-r)\). Differentiating the phases bounds the absolute third derivatives by
\[
 \frac{14+336y+756y^2+432y^3}{|s|^3},
 \qquad
 \frac{2+120y+540y^2+432y^3}{|s|^3},
\]
respectively. At the radii just specified these are less than \(1/5\) and \(5\). If \(s_{\min}\) is \(34/5\) or \(99/100\), the local Gaussian replacement and its completion to the whole real line cost, after division by \(e^{NS_*}\),
\[
 \frac1N\left[
 \frac{B}{(s_{\min}-\rho)c^2}
 +\frac{1}{s_{\min}(s_{\min}-\rho)c}
 +\frac{1}{s_{\min}c\rho}\right].                          \tag{3}
\]
In fact the exponential Taylor error is bounded by \(NB|t|^3 e^{-Nc t^2}\), the amplitude difference by \(|t|/(s_{\min}(s_{\min}-\rho))\), and the two Gaussian tails by \(e^{-Nc\rho^2}/(Nc\rho)\). The displayed constants in (3) are less than \(11\) and \(43\).

For \(J_N\) deform to the vertical segment \([3,3+i\Im s_J]\) and the horizontal ray from its endpoint. The right vertical side at infinity vanishes. The following bounds control the nonlocal part.

- For \(3\le x\le 28/5\), the modulus of the exponential base on the deformed path is at most \((1+1/450)^7 f_{3200}(x)\), where \(f_R(x)=x(R+x^6)e^{-x}\). The same bound at \(x=3\) applies on the vertical segment. This follows from \(|x+iy|\le x(1+y^2/(2x^2))\) and \(|y|\le1/5\). On \(3\le x\le28/5\), the sign of \(f_R'(x)\) is the sign of \(x^6(7-x)/(x-1)-R\). The ratio is increasing here, since its derivative has the sign of \(-x^2+7x-7\) between 1 and 7. Thus the maximum of \(f_{3200}\) is at an endpoint. Both endpoint values are below \(705\), and \(705(1+1/450)^7<e^{6.575}\). The base loses at least \(1/200\) in log modulus from the saddle height.
- For \(28/5\le x\le9\), the second derivative of \(\Re\Phi(x+i\Im s_J,w)\) is bounded above by \(-3/(2x^2)\le-1/54\). Indeed if \(r=3200/(28/5)^6\), the correction to the log curvature has modulus at most \((42r/(1-r)+36r^2/(1-r)^2)/x^2\); the term from \(7\log s\) has negative magnitude at least \(7(1-1/28^2)/(1+1/28^2)^2x^{-2}\). Using stationarity, outside the local interval the log loss is at least \(\rho^2/108=1/192\).
- For \(x\ge9\), the real phase derivative is at most \(7/x+6r_x/(x(1-r_x))-1<-1/5\) with \(r_x=3200/x^6\). The log loss at nine is at least \((9-36/5)^2/108=3/100\).

Including \(1/|s|\), these three parts cost at most \(e^{NS_J}[3e^{-N/200}+(5/(9N))e^{-3N/100}]\). The total error including (3) is below \(13e^{NS_J}/N\) at \(N\ge10000\). Here \(3e^{-N/200}<1/N\) follows from the fourth-order exponential lower bound. The Gaussian leading factor has modulus greater than \((3/4)e^{NS_J}/\sqrt N\) by \(|s_J|<36/5,\ |a|<1/5,\ \pi>3\).

For \(K_N/w^N\), keep \([0,1/2]\) on the real axis and deform the rest to the horizontal segment at height \(\Im s_K\), with vertical segments at \(1/2\) and 3. For \(1/2\le x\le3/2\) the second derivative of \(\Re\Psi(x+i\Im s_K,w)\) is at most \(-1/(2x^2)\le-2/9\). In fact the \(\log s\) term contributes negative magnitude at least \((1-(1/50)^2)/(1+(1/50)^2)^2x^{-2}\), and the correction has modulus at most \((30r/(1-r)+36r^2/(1-r)^2)/x^2\) for \(r=(151/100)^6/1000\). Off the local interval this gives log loss at least \(\rho^2/9=1/225\).

On the horizontal interval \([3/2,3]\) the modulus of the base is less than \(17/50\). Indeed the real envelope is \(f_{1000}(x)/1000\), with its maximum at the endpoints by the same monotonicity argument, and the modulus correction is at most \((1+(1/100)^2/(2(3/2)^2))^7\). The same bound holds on the vertical segment at three. At \(1/2\) use \(|s|\le 1/2+1/100\) to get an upper bound below \(17/50\). The real segment \([0,1/2]\) costs at most \((17/50)^N\): the positive envelope \(x^{N-1}(1+x^6/1000)^Ne^{-Nx}\) increases there for \(N\ge2\), and its endpoint bound times the interval length is the base bound to the \(N\)-th power. Since \(17/50<e^{-1-1/500-1/225}\), all complements including the concavity interval cost at most \(5e^{NS_K-N/225}\), with the amplitudes included. For \(N\ge10000\) this is less than \(e^{NS_K}/N\), again by a fourth-order exponential bound. The total error is below \(44e^{NS_K}/N\), and the leading factor has modulus greater than \(2e^{NS_K}/\sqrt N\). This proves (2) with room in both constants. The errors in (2) are analytic since the original pieces and the nonzero saddle approximants are analytic on neighborhoods of the stated domains.

### Cancellation on the overlap has a simple zero

The region \(|w|\le1100\) has no zero. On the real short interval \(x(1100+x^6)e^{-x}\le1100/e+3^7e^{-3}<520\). Splitting at \(1/2\) to handle \(dx/x\) as above gives \(|K_N|\le6(520)^N\). By (2), \(|J_N|>e^{NS_J}/(2\sqrt N)>(720)^N/(2\sqrt N)\). The ratio is bounded by \(12\sqrt N(3/4)^N<1\).

The region \(|w|\ge3080\) also has no zero. On \([3,4]\) bound the real exponential base of \(J_N/w^N\) by \(3e^{-3}+4^7e^{-4}/3080<8/25\). On \([4,\infty)\) use \(4e^{-4}+7^7e^{-7}/3080<8/25\). The envelope \(x(1+x^6/3080)e^{-x}\) has logarithmic derivative at most \(7/x-1\), hence at most \(-3/10\) for \(x\ge10\). Integrating the envelope gives \(|J_N/w^N|<4(8/25)^N\). Equation (2) gives \(|K_N/w^N|>(9/25)^N/\sqrt N\), and \(4\sqrt N(8/9)^N<1\). Both ratio inequalities follow already from \(e^v\ge v^2/2\).

On the remaining annulus \(1100\le|w|\le3080\), a disk of radius 100 around \(w\) lies in both saddle domains. Cauchy's estimate applied to (2) yields \(|\epsilon_*'/(1+\epsilon_*)|\le1/300\). Also
\[
 |s_*'(w)|\le\frac{|s_*-1|^2}{6|s_*|^7|a(s_*)|},
 \qquad
 \left|\frac{d}{dw}\log\big(s_*^{-1}a(s_*)^{-1/2}\big)\right|<1/50 .
\]
The first bound is implicit differentiation of the saddle equation, using \(w+s_*^6=6s_*^6/(s_*-1)\); the second uses \(|a_s|\le7/|s|^2+14/|s|^3\) and the disk bounds. Differentiation of the stationary phases now gives
\[
 \left|\frac{J_N'}{J_N}-\frac{N}{w+s_J^6}\right|<\tfrac1{20},
 \qquad
 \left|\frac{K_N'}{K_N}-\frac{N}{w+s_K^6}\right|<\tfrac1{20}.
\]
The two main terms are separated by more than
\[
 N\left(\frac1{3102}-\frac1{(34/5)^6-3100}\right)>\frac N{3300}.
\]
Here \(|s_K|^6<2\) and \(|w|\le3080<3100\). At a zero of \(I_N\) in this annulus, \(K_N=-J_N\ne0\), so \(I_N'=J_N(J_N'/J_N-K_N'/K_N)\ne0\). This excludes multiple roots uniformly at \(N\ge10000\).

## Finite orders and strictness

Write \(T_N(t)=t^N A_N(1/t;N)\), so \(T_N(0)=1\). To check \(2\le N<10000\), let \(p\) be the preceding prime, allowing \(p=N\), and write \(N=p+i\). Assume \(p>\max(6,7i)\) and \(p^2>7N-1\). If \(i>0\) and \(k=bp+c\), \(0\le c<p\), the exact coefficient valuation is
\[
 v_p([t^k]T_N)=6b+\left\lfloor\frac{i-1+6c}{p}\right\rfloor
                       +\mathbf 1_{\{c>i\}}.
\]
This follows by counting multiples of \(p\) in the factorial interval and applying the one-borrow criterion to the binomial coefficient. Thus the Newton polygon is horizontal through \(i\), with its only other edge of width \(p\) and slope \(6/p\). For \(c\le i\) the floor is zero, and for \(c>i\) the ordinate exceeds \(6b+6(c-i)/p\), as \(7i-1\ge0\). On the horizontal edge the coefficients reduce to those of \(T_i\) by the binomial coefficient congruence and the shift \(N\equiv i\bmod p\). If \(i=0\) the polygon has only the width-\(p\) edge: the leading valuation is six, and each intermediate valuation \(1+\lfloor(p-1+6k)/p\rfloor\) is strictly above \(6k/p\).

The denominator \(p\) excludes a squared factor from the width-\(p\) edge. Any rational repeated factor therefore has only unit roots at this prime. Its primitive integer normalization has constant coefficient \(\pm1\) because it divides \(T_N\), so it reduces with its full degree. Its degree must be even, since \(A_N(\cdot;N)\) has no repeated real root. Indeed it is the integral of an even power when \(N\) is even, and has strictly positive derivative on the real axis when \(N\) is odd. A gcd of degree at most one for \(T_i\bmod p\) and its derivative thus excludes a repeat. At the remaining orders, reduction at a prime \(\ell>7N\) preserves the full degree directly and the same gcd test applies.

The rational [replay](verify_gap_six_resultant_strict.py) checks the contour constants using Taylor bounds for the exponential and all 9,998 finite terminal orders using these criteria. The [receipt](gap_six_resultant_strict_verification.json) gives the counts and all full-degree fallback orders. The analytic contour argument is the all-order proof above the cutoff.

By (1) the ordinary two-shape resultant never vanishes. Its nonnegativity for positive shape follows also by evaluating the even reciprocal at the roots of the odd one. The real root contributes a strictly positive factor and each conjugate pair contributes a nonzero squared modulus. Both normalizing divisors are positive at \(\beta=1\). This proves the strict assertions and the support-\(\{1,7\}\) conclusion. The proof of strictness uses the physical shape \(\gamma=N\); the all-positive-shape question remains open.
