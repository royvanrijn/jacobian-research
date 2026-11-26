# Height descent at zero-bias contacts for opposed growing gaps

Fix a finite ordinary core and two adjacent orders \(n,N=n+1\), with
\(n\ge2\). Append \(r\) copies of \((L,L+d)\) and \(s\) copies of
\((Y+e,Y)\), where the heights, positive gaps \(d,e\), and positive
counts \(r,s\) may vary. Write
\[
 M=L+d/2,\quad M'=Y+e/2,\quad h=rd-se,\quad B=rd,\quad
 X=\max(d,e),\quad D=\min(d,e).
\]
Keep the signed tail gap \(h\) fixed. This loses only a finite
subsequence choice at a possible collision by the
[signed imbalance bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md).
Use the midpoint curvature and bias parameters
\[
 x=\frac{r d^2}{M},\quad y=\frac{s e^2}{M'},\quad
 T=x+y,\qquad \lambda=\frac{x}{d}-\frac{y}{e}.
\]
These midpoint parameters give the same diffuse limits as the
[endpoint parameters](MULTIVARIATE_FIXED_TERMINAL_COMPACTIFICATION.md#uniform-estimates-for-rows-of-large-size).
We consider the branch \(X\to\infty\) approaching a specified
positive zero-bias Gaussian reference contact \((T_*,0)\) for the
fixed core and signed gap \(h\), with the two contacting roots simple.
Either curvature share may tend to zero.

Let \(p\) denote the derivative in \(T\) of the logarithm of the
ratio of the two endpoint-normalized reference roots at \(T_*\).
Let \(Z\) denote the first derivative of that log ratio in the
coefficient log direction \(j\psi(k/j)\) at both orders, where
\(\psi(u)=u(1-u)(2u-1)/12\).
At a transverse contact \(\Im p\ne0\), put
\[
 C=\frac1{2nN},\qquad K=\Re Z-\frac{\Re p}{\Im p}\Im Z,
 \qquad W(k)=C+\frac{T_*K}{k}\quad(k\ge1).
\]

**Transverse theorem.** Any sequence of ordinary collisions approaching
the specified simple transverse zero-bias contact satisfies
\[
 M,M'\asymp B X,\qquad W(r)M'-W(s)M=O(X).                       \tag{1}
\]
The implied constants are uniform along the branch. Hence there is
no such collision sequence in any of the following cases.

- \(K=0\), with independently varying heights, gaps and counts.
- \(r=s\), with independently varying heights and growing gaps.
- \(|d-e|=o(\log X)\), with independently varying heights and counts.
- A usable prime-power part of one first-block factor of the greater-count
  type grows faster than \(X(1+|d-e|)\), as defined below.

The first case allows one gap to stay bounded and needs no lower bound
on either individual curvature share.

**Tangential theorem.** Suppose \(\Im Z\ne0\). Assume the real part
of the unperturbed logarithmic contact equation has exact order
\(m\ge1\) at \(T_*\), and its imaginary part vanishes to order at
least \(2m\). Fix positive integers \(a,b\) and assume
\(r=a k,\ s=b k\) for a positive integer \(k\). The two opposing gaps
then grow comparably. There is no such collision sequence if
\[
 k=o(\sqrt{\log D}).
\]
If the angular order is greater than \(2m\), or if the counts are
equal, the condition improves to \(k=o(\log D)\). In particular,
bounded independently varying counts are excluded under these
contact hypotheses by considering their finitely many ratios.
The case \(m=1\) includes \(\Im p=0,\ \Re p\ne0\).

The exclusion uses all-degree coprimality of the resulting exact
integer clusters. No information about the prime factors of their
initial first marks is needed.

## A uniform skew estimate

At orders \(j=n,N\), let \(z_j(T,V)\) be the continued simple root
after multiplying the endpoint-normalized reference coefficients
by \(\exp(V_{j,k})\). Write
\[
 H(T)=\log(z_N(T,0)/z_n(T,0))-C_R,\qquad
 C_R=\log(G_N^C/G_n^C)+h\log(N/n),
\]
where \(G_j^C\) is the positive endpoint center of the fixed core.
The log branches are chosen locally so that \(H(T_*)=0\).
For a signed row gap \(\delta\) with midpoint \(m\), the exact
endpoint-normalized log contribution is
\(j L_j(u)\) at \(u=k/j\), where \(L_j(0)=L_j(1)=0\).
At \(t=m+\delta(u-1/2)\) its second derivative is
\[
 L_j''(u)=j\delta^2\psi_1(1+jt).
\]
With \(|\delta|/m\le1/2\), the
[trigamma expansion](MULTIVARIATE_THIRD_ORDER_BOUNDARY_AND_INTERIOR_REDUCTION.md#exact-row-expansions-and-the-projected-coefficients)
gives \(j\psi_1(1+jt)=1/t+O_j(t^{-2})\). Expanding only the leading
term at \(m\), and using \(\phi''=-1\), \(\psi''=-(u-1/2)\) for
\(\phi(u)=u(1-u)/2\), yields
\[
 jL_j(u)+j\frac{\delta^2}{m}\phi(u)
   = j\frac{\delta^3}{m^2}\psi(u)
           +O_j\!\left(\frac{\delta^2}{m^2}+\frac{\delta^4}{m^3}\right).
                                                                    \tag{2}
\]
The integration bound follows by subtracting the two displayed
zero-endpoint polynomials and comparing the bounded second
derivative of the remainder with the two constant-curvature
parabolas. The constants do not depend on \(\delta\).

The exact signed adjacent-center term of this row obeys
\[
 E(w)+\frac{\delta}{2nN m}
      =O_{n,N}\!\left(\frac{|\delta|}{m^2}+\frac{|\delta|^3}{m^3}\right).
                                                                    \tag{3}
\]
Indeed, the integral in the same center formula has leading
integrand \(-1/(2v^2 t)\) for the order variable \(v\in[n,N]\)
and exponent coordinate \(t\in[m-\delta/2,m+\delta/2]\), with
error \(O_{n,N}(t^{-2})\). The oriented integral of \(1/t\)
is \(\delta/m+O(|\delta|^3/m^3)\) by expanding around the midpoint.

Put \(Q=x/d+y/e\), and use \(t=T-T_*\) from now on.
Along the branch \(B\to\infty\), \(se=B-h\), and
\[
 \frac dM=\frac{x}{B},\quad \frac e{M'}=\frac{y}{B-h},
 \qquad Q\gg \frac1B.
\]
The last inequality follows from \(T\asymp1\) and \(X\le B+|h|\).
Thus the relative gap condition of (2) and (3) holds eventually.
After summing (2), the first two error types become
\(x^2/(Bd),y^2/((B-h)e)\) and \(x^3/B^2,y^3/(B-h)^2\),
respectively. The first pair is \(O(Q/B)\). The second pair is
\(O(B^{-2})\), hence also \(O(Q/B)\).
The signed cubic coefficients from the two types are \(x^2/B\)
and \(-y^2/(B-h)\). Replacing \(B-h\) by \(B\) changes their sum
by \(O(B^{-2})\). Consequently the total compensated shape is
\[
 V_{j,k}=\frac{T(x-y)}{B}\,j\psi(k/j)+O_{n,N}(Q/B).                \tag{4}
\]
Summing (3) gives a compensated center \(V_E=O(Q/B)\),
since \(B/M=x/d\) and \(1/M=O(B^{-1})\), and likewise for
the opposite row. In particular \(V=O(B^{-1})\).

At an ordinary common root continuing the specified branches the
exact complex equation is
\[
 H(T)+C\lambda+
 \bigl(\log(z_N(T,V)/z_n(T,V))-\log(z_N(T,0)/z_n(T,0))\bigr)-V_E=0 .
\]
This follows directly from the ratio of the two endpoint centers,
since \(V_E=\sum E(w)+\lambda/(2nN)\) over the added rows. Simplicity
gives a jointly analytic root response, and (4) gives
\[
 0=H(T)+C\lambda+ Z\,T_*\frac{x-y}{B}
                   +O\!\left(\frac{|t||x-y|}{B}+\frac QB\right).
                                                                    \tag{5}
\]
The \(O(\|V\|^2)\) term is absorbed because \(Q\gg B^{-1}\).
The first error in (5) comes from replacing \(T\) and the shape
derivative at \(T\) by their reference values. All estimates here
apply independently of which share tends to zero.

## The transverse height relation

The imaginary part of (5) gives \(t=O(B^{-1})\) when \(\Im p\ne0\),
using \(\lambda\in\mathbb R\) and boundedness of \(x,y,Q\).
Apply the real linear projection \(z\mapsto\Re z-(\Re p/\Im p)\Im z\)
to (5). The result is
\[
 C(x/d-y/e)+ K T_*(x-y)/B=O(Q/B).
\]
With \(\alpha=x/d\) and \(\beta=y/e\), replacing \((x-y)/B\)
by \(\alpha/r-\beta/s\) costs only \(O(B^{-2})\).
We obtain
\[
                         W(r)\alpha-W(s)\beta=O(Q/B).          \tag{6}
\]
The [algebraic exceptional-clock result](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md)
gives \(\exp[-T_* /(2nN)]\) algebraic. In the fixed finite ordinary
core the normalized reference roots are algebraic, and differentiating
their simple root equations gives \(p,Z\) algebraic. Real and
imaginary parts of these numbers are algebraic, so \(K\) is algebraic
real. The [Hermite--Lindemann consequence](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core)
makes \(T_*\) transcendental. Hence \(W(k)\ne0\) for all positive
integers \(k\). Since \(W(k)\to C>0\), their magnitudes are bounded
above and below uniformly.

Equation (6) and \(\alpha,\beta>0\) imply \(\alpha\asymp\beta\):
a term whose size is a fixed fraction of \(Q=\alpha+\beta\)
cannot be matched by an \(o(Q)\) term. Since \(d\alpha+e\beta=T\asymp1\),
we obtain \(\alpha,\beta\asymp X^{-1}\). This proves the height
comparability in (1). Multiplying (6) by \(MM'/B\) gives the second
relation in (1), since \(\alpha=B/M\), \(\beta=(B-h)/M'\), and
the correction from \(h/B\) has size \(O(X)\).

If \(K=0\) or \(r=s\), equation (1) implies \(|M'-M|=O(X)\).
In general it implies
\[
 |M'-M|\ll X+ B X\,|1/r-1/s|
           \ll X(1+|d-e|),
\]
because \(B/r=d\) and \(B/s=e+O_h(1)\). Each of the first three
hypotheses gives a window of width \(o(X\log X)\). The
one-mark hypothesis below gives width \(O(X)\). The scaled
interval closure excludes the collision in each case.

### A single deep mark

Reverse both columns if necessary so that \(r\ge s\).
For either order \(m\in\{n,N\}\) choose
\(F=mL+u\) with \(1\le u\le d\) and \(m\nmid u\). Define its
usable prime-power part using the \(K_{\rm mark},B_0,b_p\)
cutoffs of the [same-order depth constraint](MULTIVARIATE_PAIRED_PRIME_MARK_AFFINE_DESCENT.md#the-same-order-depth-constraint):
\[
                       P=p^{(v_p(F)-B_0 b_p)_+}.
\]
The same depth proof applies to this first-block factor. It is
present at every positive index at order \(m\), and the
cross determinant with an opposite-order factor is \(mt-m'u\).
This determinant is nonzero since \(m,m'\) are consecutive and
\(m\nmid u\); its absolute value is still at most \(K_{\rm mark}\).
Consequently \(P\mid mY+t_m\) for some \(1\le t_m\le me\) at a
possible common root. It follows that
\[
   P\ \mid\ m(Y-L)+t_m-u,\qquad
          |m(Y-L)+t_m-u|\ll X(1+|d-e|).
\]
If \(P/[X(1+|d-e|)]\to\infty\), the difference is zero, so
\(|Y-L|\le X\). Equivalently the cofactor \(F/P\) is
\(o(B/(1+|d-e|))\), since (1) gives \(F\asymp BX\).
In contrast with the paired-mark criterion, this local argument
only needs the one mark, at either order of the reference pair.

## The angularly tangential branch

Write \(r=a k,s=b k\) with fixed \(a,b>0\), so \(B\asymp kD\),
\(X\asymp D\) and \(Q\asymp D^{-1}\). The real part of (5),
with \(\lambda=O(D^{-1})\), yields \(|t|^m=O(D^{-1})\) by the
hypothesis on the first nonzero radial term. Since \(\Im Z\ne0\),
the imaginary part yields
\[
                  |x-y|\ll B|\Im H(T)|+D^{-1},                 \tag{7}
\]
after absorbing the error \(O(|t||x-y|)\) in that inequality.
The angular order at least \(2m\) gives \(|x-y|=O(k/D)\).
If the angular order is greater than \(2m\), then it is at least
\(2m+1\), and (7) improves to \(|x-y|=O(D^{-1})\) when
\(k=o(\log D)\). In either stated range the two shares tend to
\(T_*/2\). The scaled midpoints are
\[
                  aM=\frac{B^2}{k x},\qquad
                  bM'=\frac{(B-h)^2}{k y}.
\]
Their difference is \(O((B^2/k)|x-y|+B/k)\). The scaled gaps are
comparable to \(D\), so the two rows form a cluster of width
\(O(k^2D)\), or \(O(kD)\) with the improved angular bound.
Both are \(o(D\log D)\) in their respective ranges. The fixed
integer multipliers \(a,b\) are permitted by the closure lemma below.

For equal counts take \(a=b=1\), so \(e-d=O_h(1/k)\). Then
\(\lambda=(x-y)/d+O(1/(kD^2))\) and the real equation in (5)
strengthens to
\[
                    |t|^m\ll |x-y|/D+1/(kD^2).
\]
Using the angular order at least \(2m\) in (7) now gives
\[
                    |x-y|\ll (k/D)|x-y|^2+1/D .
\]
The shares are bounded, so the quadratic term is absorbed when
\(k=o(D)\). Thus \(x-y=O(D^{-1})\), and the unscaled window has width
\(O(kD)\). This is \(o(D\log D)\) for \(k=o(\log D)\) and the
same closure applies.

## All-degree closure for wider scaled clusters

Use the [scaled interval construction](MULTIVARIATE_PAIRED_PRIME_MARK_AFFINE_DESCENT.md#a-closure-lemma-for-clustered-scaled-intervals)
with a fixed finite bound on the number of distinct scaled endpoint
intervals and fixed integer scale multipliers \(c_i\mid C_0\).
Let \(X\) denote the largest scaled gap, and \(R\) the width of the
window containing all scaled intervals. Its proof still gives
all-degree coprimality for \(X\to\infty\) and \(R=o(X\log X)\).

To verify the changed width hypothesis, choose the subinterval
from that construction. Its length is \(\gg X\) since the type
bound and scale multipliers are fixed. It lies inside one of the
long types' first-block windows at every order \(2\le j<J\), and
a mark there is either first-block or absent for every other type,
where \(J\) is the profile-independent cutoff. There are \(\gg X\)
consecutive marks \(w\) with scaled level \(C_0 w\) in the interval.
At a putative collision the curvature bound from that proof forces
the largest scaled exponent to be \(\gg X^2\). Since
\(R=o(X\log X)\), every selected mark is \(\gg X^2\).
If all prime divisors of the marks lay below \(JR+C_0+J\), the
prime deletion argument would remove only
\(O((R+1)/\log(R+2))=o(X)\) marks, using the elementary prime-count
bound there. The product of the remaining marks would divide the
factorial of one less than the number of consecutive marks.
The number of marks is \(O(X)\), and each remaining mark is
\(\gg X^2\), so this is impossible by comparing logarithms.
The resulting prime exceeds every nonzero scaled factor displacement
from the mark. The literal valuation arrays and the consecutive
binomial root-quotient argument are therefore exactly those of the
linked proof. The fixed retained core is harmless at these primes.

The two contact tests reduce to this exact ordinary argument.
The remaining variable-gap profiles include transverse zero-bias
contacts with nonzero projected cubic response and widely separated
gaps, and contacts outside the stated angular hypotheses.
