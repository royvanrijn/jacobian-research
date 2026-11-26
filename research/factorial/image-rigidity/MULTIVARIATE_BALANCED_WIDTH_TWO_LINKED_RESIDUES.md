# Linked prime residues for the balanced width-two factorial template

Fix \(d\ge2\), a finite ordinary two-monomial core \(C\), and
positive integers \(M,r\). Append \(r\) copies of each row
\[
 (M,M+d),\quad (M+1,M),\quad
 (M+d,M+d+1),\quad (M+2d,M+d).                                \tag{1}
\]
The signed gap of the appended rows is zero. Put \(m=2d\),
\(N=m+1\), and \(g=m-1\). This is the four-type family for
which every leading literal affine count polygon leaves at
least two shared slope units at \(m,N\). See the
[narrow-band obstruction](MULTIVARIATE_NARROW_HEIGHT_BAND_NATIVE_MARKS.md#a-balanced-obstruction-to-using-just-the-literal-count-polygons).

There is a computable \(M_0(d,C)\) such that for \(M>M_0\)
the original profile (1) satisfies physical SFC(2) in every
degree whenever \(g\nmid r\). At the possible exceptional
pair \(m,N\), the proof first restricts any common factor to
a single irreducible quadratic. It then gives two exact
linked residue conditions valid also when \(g\mid r\).
In that remaining class a common root forces
\[
 (M+q)\ \mid\ K_0\,\operatorname{num}(C_q T_q^r-1),
             \qquad q=1,d+1,                                 \tag{2}
\]
where \(K_0,C_q,T_q\) are explicit nonzero constants fixed
by \(C,d\), \(K_0\) is a positive integer, and \(\operatorname{num}\)
denotes a reduced numerator, with \(\operatorname{num}(0)=0\).
This gives an exact divisor test on both native marks.
Coprimality over \(\mathbb Q\) extends to any
characteristic-zero field.

## Reduction to one pair and one quadratic

Choose \(J>2d+1\), \(J\ge16\), at least the
[effective absolute terminal cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL.md).
For the four appended types use the envelope argument from the
[narrow-band proof](MULTIVARIATE_NARROW_HEIGHT_BAND_NATIVE_MARKS.md#finding-a-mark-before-considering-primes)
with every type selected. The two upper types have lower
exponent \(M+d\). The two lower types have upper exponent at
most \(M+d\); hence they cannot furnish a unit line active
at the middle index of a three-index interval where that
argument fails. Indeed at every strict interior index both
upper type arguments exceed \(j(M+d)\). The active unit row
at the middle index must be the upper increasing row, which
ends at \(Z=M+d+1\). In the notation of that argument the
required opposite gap at the even order \(j\) then satisfies
\(2D'=j(Z-b)\), with the opposite type attaining \(jZ\) at
index \(j-2\). Only the upper decreasing row can attain that
height. It has \(D'=d\) and \(b=Z-1\), so \(j=2d\).
Even terminals already have the odd-lower fallback. The
explicit valuation bound in the narrow-band proof therefore
excludes every pair up to \(J\) except \(m,N\) for sufficiently
large \(M\). The complete low-order result and the choice of
\(J\) cover the two ends of the degree range.

At order \(m=2d\), mark \(x=m(M+1)+1\).
The long increasing row first samples \(x\) at coefficient
index three, because \(\lceil(m+1)/d\rceil=3\).
The other three types never sample \(x\): the short decreasing
row ends its full ramp at \(m(M+1)\), and both upper types
start above \(x\). Thus \(I=\{0,1,2\}\) is the interval
without the mark. There are at least \(r\) literal copies
of \(x\) at every outside index. At the other order no
proportional factor is present since \(m\nmid x\). Apply the
valuation comparison (4) in the narrow-band proof with
\(W=d,D=d\), using this width-two interval.
Its constant divisor bound would be forced on \(x\) if
every prime of \(x\) were below the depth threshold.
For \(M\) above that bound, some prime therefore forces every
shared face into \(I\). By the Newton polygon product rule
any common rational factor has degree at most two. The even
ordinary moment is strictly positive on the real axis, so a
nonconstant common divisor then has degree exactly two and
is irreducible. We use the Newton rule in the
[Filaseta–Rocks formulation, Theorem 2.1](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf),
with the horizontal order reversed to read coefficients from
constant to leading term.

## Native integer marks and exact linked ratios

For \(q=1\) use both orders in their original orientations.
For \(q=d+1\) reverse all rows of the parents, including
the core. This reciprocates their roots and preserves physical
coprimality. In either orientation write the four tails as
\((M+a_i,M+b_i)\), \(1\le i\le4\).

Strip each smaller-end factorial at each order. Write
\(c_{j,k}\) for the positive integer stripped core coefficient,
including its binomial coefficient. Define the finite nonzero
integer product
\[
 U_{j,k}(q)=
 \prod_{i=1}^4
   \prod_{\substack{j\min(a_i,b_i)<h\le(j-k)a_i+kb_i\\h\ne jq}}
                  (h-jq).
\]
Put \(a_m'=1,\ a_N'=2\), the last indices of the two
low plateaus at the integer mark. Define
\[
 \begin{split}
 C_q&=\frac{c_{m,m}/c_{m,1}}{c_{N,N}/c_{N,2}},\\
 T_q&=\frac{m}{N}
          \frac{U_{m,m}(q)/U_{m,1}(q)}
               {U_{N,N}(q)/U_{N,2}(q)} .
 \end{split}                                                   \tag{3}
\]
For either \(q\), the mark \(j(M+q)\) is sampled once per
copy of the short decreasing row at \(k=0\), and once per
copy of the long increasing row at \(k\ge a_j'+1\).
The remaining indices \(1\le k\le a_j'\) sample no copy.
For \(q=1\) these two types begin at smaller exponent \(M\),
and both upper types begin above the mark. For \(q=d+1\)
in the reversed orientation the two marked types begin at
\(M+d\), and the other two upper endpoints are below the mark.
The cutoff for the long row is \(k=\lceil j/d\rceil\), which
gives the stated \(a_j'\) at the two orders.

Let \(B_j\) be \(j!\) times the product of the full stripped
core ramps at order \(j\). For each prime \(p\) define
\[
 L_p=\max_{1\le h\le J(2d+1)}v_p(h),\qquad
 D_p=\max_{j\in\{m,N\}}v_p(B_j)+(2d+3)J L_p,\qquad
 K_p=4gD_p+L_p .
\]
Only primes in a fixed finite set \(\mathcal S\) have
\(K_p>0\). Put \(K_0=\prod_{p\in\mathcal S}p^{K_p}\).
Let \(E=v_p(M+q)>0\) and suppose \(E>K_p\).
For either order, every tail ramp factor other than the mark
is \(h-jq+j(M+q)\) with \(0<|h-jq|\le J(2d+1)\).
It therefore has exactly the valuation of \(h-jq\), and the
quotient of the actual factor by \(h-jq\) is \(1\) modulo
\(p^{E-L_p}\). This remains true after raising a product of
these quotients to the integer power \(r\).

Let \(\epsilon_{j,k}=0\) at \(1\le k\le a_j'\) and
\(1\) elsewhere. The actual stripped coefficient \(P_{j,k}\)
is consequently
\[
 P_{j,k}=A_{j,k}(M+q)^{r\epsilon_{j,k}}(1+\delta_{j,k}),
 \qquad
 A_{j,k}=c_{j,k}\big(j^{\epsilon_{j,k}}U_{j,k}(q)\big)^r,
 \qquad v_p(\delta_{j,k})\ge E-L_p.                           \tag{4}
\]
Write \(h_{j,k}=v_p(A_{j,k})\). Since \(c_{j,k}\mid B_j\)
and the total gap of the four rows is \(2d+2\),
\(0\le h_{j,k}\le rD_p\). The ordinates are
\(rE\epsilon_{j,k}+h_{j,k}\). The first Newton edge runs
from index zero to one and has slope
\(-rE+h_{j,1}-h_{j,0}\). The second order has a middle edge
from one to two, of slope \(h_{N,2}-h_{N,1}\). Both orders
have a final edge from \(a_j'\) to \(j\), of width \(g\)
and slope
\[
                         s_j=\frac{rE+h_{j,j}-h_{j,a_j'}}{g}.
\]
To check the hull, the negative first slopes lie strictly below
\(-rD_p\); the middle slope is within \([-rD_p,rD_p]\);
and both final slopes exceed \(rD_p\). At each interior
index of a proposed final edge, the gap from the line is at
least \(r(E/g-2D_p)>0\), by interpolating the two ordinates
and using the bounds on \(h_{j,k}\). The other proposed
edges have width one. This proves the shape.

A common quadratic must use one unit of width on a final
edge because the only lower-order edge outside it has width
one. The separation of the three slope ranges forces
\(s_m=s_N=s\). The reduced denominator of \(s\) divides
the width used by the rational factor at that slope.
This width is one or two, whereas the denominator divides the
odd integer \(g\). Thus \(s\) is an integer, and
\[
 g\mid rE+h_{m,m}-h_{m,1},\qquad
                h_{m,m}-h_{m,1}=h_{N,N}-h_{N,2}.             \tag{5}
\]

Let \(\alpha\) be a root of the common factor on this slope
in a finite extension of \(\mathbb Q_p\), and put \(w=p^s\alpha\),
a unit. Substitute \(z=p^{-s}w\) and divide the two equations
by \(P_{j,a_j'}p^{-s a_j'}w^{a_j'}\). The endpoint part is
\(1+p^{-sg}(P_{j,j}/P_{j,a_j'})w^g\). Every other term has
valuation at least \(r(E/g-2D_p)\); for indices before
\(a_j'\) this follows from adding \(s(a_j'-k)\) and
using \(h_{j,k}-h_{j,a_j'}\ge-rD_p\).
By (4) the congruence of the two cleaned endpoint ratios
holds to the smaller of this off-edge precision and \(E-L_p\).
In fact (5) makes the cleaned endpoint term a unit in each
equation, and equating the two unit terms cancels
\(w^g p^{-sg}(M+q)^r\). This gives
\[
         v_p(C_q T_q^r-1)
             \ \ge\ \min\{E-L_p,\,r(E/g-2D_p)\}.             \tag{6}
\]
For \(r\ge g\) the right-hand side is at least \(E-K_p\).
Primewise this proves (2), also including primes with \(E\le K_p\)
directly. If the rational difference in (2) is zero, the
divisibility condition is automatic.

## An effective cutoff when \(g\nmid r\)

For a prime \(p\notin\mathcal S\) dividing \(M+1\) or
\(M+d+1\), (5) gives \(g\mid r v_p(M+q)\) because all
\(h_{j,k}\) vanish. When \(g\nmid r\), choose any prime
\(\ell\mid g\) with \(v_\ell(r)<v_\ell(g)\). It is odd.
The valuation of each \(M+q\) outside \(\mathcal S\)
must be divisible by \(\ell\). Therefore
\[
        M+1=u X^\ell,\qquad M+d+1=v Y^\ell,
        \qquad vY^\ell-uX^\ell=d,                           \tag{7}
\]
where \(X,Y\) are positive integers and \(u,v\) range
over the finite list of products of primes in \(\mathcal S\)
with exponents in \([0,\ell-1]\).

If \(u=v\), the difference of powers in (7) bounds \(X,Y\)
effectively by the fixed \(d\), because \(Y>X\) and
\(Y^\ell-X^\ell\ge(X+1)^\ell-X^\ell\).
If \(u\ne v\), then \(v/u\) has a prime valuation nonzero
modulo \(\ell\). This makes \(X^\ell-v/u\) irreducible over
\(\mathbb Q\): a monic rational factor of intermediate degree
\(k\) would have constant coefficient of absolute value
\((v/u)^{k/\ell}\), which is irrational by that prime
valuation. Thus \(vY^\ell-uX^\ell\) is an irreducible integral
binary form of degree \(\ell\ge3\). The effective Thue theorem
bounds every integral solution with right-hand side \(d\ne0\).
See [Tzanakis and de Weger, Sections I–II](https://doi.org/10.1016/0022-314X(89)90014-0)
for the applicability and effective bound.

The finite maxima over \(u,v\) and the prime divisors \(\ell\)
of \(g\), together with the earlier thresholds, give a computable
bound \(M_0(d,C)\) independent of \(r\). This completes the
all-degree assertion in the stated multiplicity class. The
linked residue conditions (2), (5), and (6) are necessary in
the remaining class \(g\mid r\). They leave that class open.

## At most a quadratic in every sufficiently high bounded-gap band

In the setting and with the threshold \(K(W,C,D,J)\) of the
[narrow height-band theorem](MULTIVARIATE_NARROW_HEIGHT_BAND_NATIVE_MARKS.md),
the gcd of ordinary consecutive moment polynomials has degree at
most two at every terminal when \(M>K\), with arbitrary gap mixture,
orientations, and row multiplicities. Any surviving gcd is one
irreducible quadratic with nonreal roots. This degree bound applies also at
the odd terminals where the narrow-band theorem does not exclude
coprimality failure.

In fact retain the same multiplicity split and mark \(x=F+1\) or
\(F+2\) at even order \(m\) in that proof. Orders below fifteen are
already closed, so \(m\ge14\). If the mark-free interval has a
triple of consecutive indices, the proof of the exceptional
envelope shows that its center is \(1\) or \(m-1\), using
\(F\equiv-1\bmod m\) and an active line of slope \(1\) or \(-1\).
An interval of four consecutive indices would have triples with
two consecutive centers, which is impossible.
The only width-two intervals are \(\{0,1,2\}\) and
\(\{m-2,m-1,m\}\).

The prime-depth comparison uses no width assumption until the
last Newton polygon step. If every prime of \(x\) were below
the depth threshold, it would force \(x\mid K\), whereas
\(x>mM>K\). A prime above the threshold places all shared
Newton widths in the mark-free interval, giving degree at most
two. Positivity of the even ordinary moment on the real axis
gives the stated irreducibility. The effective absolute cutoff
covers all terminals at least \(J\).
