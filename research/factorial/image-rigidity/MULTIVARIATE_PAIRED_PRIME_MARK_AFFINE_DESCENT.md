# Paired first marks force affine descent for independent opposed heights

Fix an ordinary binary factorial core \(C\), with the equal rows
omitted. Append \(r\) copies of \((L,L+d)\) and \(s\) copies of
\((Y+e,Y)\), where \(r\ge s\ge1\), \(d,e\ge1\), and the smaller
exponents \(L,Y\) are nonnegative integers. Both heights, both gaps
and both counts may vary independently. Either orientation can be
the selected one, by reversing both columns if necessary.

The mechanism is an orientation-sensitive condition at each first
mark. A sufficiently deep prime-power divisor of \(jL+1\) must divide
one of the literal factors \(jY+t\) in the opposite row at the
**same order**. Large prime-power parts at both adjacent orders then
force an exact integral affine relation between the heights. A
large prime-power part at just one order also suffices near positive
diffuse contacts, using effective linear forms in logarithms.

## The same-order depth constraint

Let \(J\ge4\) be an integer at least the
[absolute eventual angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Let \(c\) be the number of unequal core rows and \(M_C\) their
maximum exponent, taken as zero for the empty core. Put
\[
 K=\max\{J^2\max(d,e),\,J M_C,\,J!\},\qquad
 B_0=1+(J+1)(c+3).
\]
For a prime \(p\), let \(b=b(p)=\max\{h\ge0:p^h\le K\}\).
For one of two adjacent orders \(m,m'<J\), both at least two,
set \(F_m=mL+1\). Define the usable prime-power part and cofactor by
\[
 P_m=p^{(v_p(F_m)-B_0 b)_+},\qquad A_m=F_m/P_m.                 \tag{1}
\]
Here \(p\) can be chosen separately at each order. At a physical
common root there is \(t_m\in\{1,\ldots,me\}\) such that
\[
                              P_m\mid mY+t_m.                   \tag{2}
\]
In particular all of a prime-power part is usable when \(p>K\).
For fixed gaps the loss \(p^{B_0b}\le K^{B_0}\) is bounded
independently of the heights and of \(p\).

To prove (2), strip the smaller factorial from every row in each
order. Remove the common affine slope of the valuation layers
\(1,\ldots,b\), as in the
[common slope argument](MULTIVARIATE_RETAINED_COMMON_SLOPE_BUDGET.md).
In each row and each layer, the floor of the upper factorial
argument divided by \(p^h\) differs from that quotient by an
error of oscillation at most one. Its linear coefficient in the
index is the same in both orders. Thus after the common slope
is removed the residual from these layers and from the binomial
coefficient has oscillation at most
\[
                              E=(r+s+c+1)b
\]
in each order. The core supplies only these low layers by the
choice of \(K\).

Assume \(a=v_p(F_m)>b\). At order \(m\) this is the sole high
factor of the selected row, because two distinct factors in any
ramp cannot differ by a multiple of \(p^{b+1}>K\). Its
high-layer valuation is zero at index zero and \(r(a-b)\) at
every positive index. At the other order the selected row has
no high factor, since
\[
                    m(m'L+t)-m'F_m=mt-m'
\]
is a nonzero integer of absolute value at most \(K\) for
\(1\le t\le m'd\). The high contributions of the decreasing
row are nonincreasing in the index. At order \(m\) their drop
from index zero to any other index is at most \(s(f-b)_+\),
where \(f=\max_{1\le t\le me}v_p(mY+t)\). At order \(m'\)
every slope of the corrected lower Newton polygon is at most
\(E\). Indeed the last edge joins a point of the corrected
ordinate array to its final point. The nonincreasing high part
cannot increase on this edge and the change of the residual
is at most \(E\). All other slopes are at most the last one.

Subtract any possible common slope \(\xi\le E\) times the
coefficient index at order \(m\). At every positive index the
resulting ordinate exceeds index zero by at least
\[
                     r(a-b)-s(f-b)_+-(m+1)E .
\]
If this is positive, index zero is the unique minimizer.
The Newton polygon product rule then precludes the common
slope. Thus a physical common root forces
\[
       a-b\le (f-b)_+ +(J+1)(c+3)b ,
\]
using \(r\ge s\). Since \((f-b)_+\le f\), (2) follows.
When \(a\le b\) the usable part is one and the assertion also
holds. Only the opposing row's same-order ramp enters (2).

## A closure lemma for clustered scaled intervals

Let \(T\) bound the number of distinct scaled endpoint intervals
of the moving rows. Use positive integer multipliers \(c_i\mid C_0\)
and a fixed integer \(W\ge1\). Equal rows
may be omitted. Write \(g_i=c_i|b_i-a_i|\), \(X=\max_i g_i\),
and let \(M\) be the largest scaled exponent. Assume the two
scaled endpoints of every moving type lie in
\([M-WX,M]\). For any fixed retained core, all sufficiently
large \(X\) give all-degree coprimality, uniformly in the
multiplicities, orientations, and scaled endpoint placements.
In particular the lesser gaps may stay bounded. The number \(T\)
may grow when
\[
    (8JT)^T\le \frac{\log X}{72 H C_0},
             \qquad H=J(W+1)+C_0+J .                          \tag{3}
\]
For instance \(T=o(\log\log X/\log\log\log X)\) suffices.

Here is an adaptation of the
[simultaneous first-block construction](MULTIVARIATE_SIMULTANEOUS_FIRST_BLOCK_ANGULAR_SEPARATION.md).
Put \(Q=8JT\). Among the ordered distinct values of \(g_i\)
there is a threshold \(D\ge X/Q^{T-1}\) such that every type
has \(g_i\ge D\) or \(g_i<D/Q\). Take the first downward jump
of ratio greater than \(Q\), appending zero to the end of the
list. For long types with \(g_i\ge D\) write \(u_i\) for their
scaled smaller endpoint and put \(u_0=\min u_i\). Every
\(y\in (u_0,u_0+D/J)\) lies in the first block of at least one
long type at each order \(2\le j<J\), in the sense that the
factor \(jy/c_i\), when integral, is one of its first
\(g_i/c_i\) factors. For each other long type this factor is
either in the first block or absent, since
\(\min_i(u_i+g_i/J)\ge u_0+D/J\).
Delete the scaled intervals of all short types from this
candidate interval. Their total lengths are less than
\(TD/Q=D/(8J)\), so there is an open subinterval of length
at least \(D/[2J(T+1)]\) avoiding all short intervals.

Take the consecutive integers \(w\) for which \(y=C_0w\)
lies in that subinterval. Their number \(k_0\) is at least
\[
 \frac{D}{4 C_0 J(T+1)}\ \ge\ \frac{X}{C_0(8JT)^T}
\]
once the interval has length at least \(2C_0\), which follows
for large \(X\) also under (3). Under the
[absolute curvature bound](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
a collision requires \(M>X^2/(12 C_0)\), by the type attaining
\(X\). Thus for sufficiently large \(X\) the marks satisfy
\(w\ge M/(2C_0)\) and have logarithms at least
\((5/6)\log M\). The number \(k_0\) is at most \(X\),
so \(\log k_0\le (2/3)\log M\) for sufficiently large \(X\).
There is a prime \(p>HX\)
dividing a mark. Indeed, otherwise remove a mark of maximal
valuation for every prime \(p\le HX\) occurring in the marks.
The product of those remaining divides \((k_0-1)!\), since
at most \(\lfloor(k_0-1)/p^a\rfloor\) terms remain at each
prime-power layer. The elementary bound
\(\pi(U)\le9U/\log U\) used in the linked first-block theorem
gives at most \(k_0/8\) deletions under (3). This
condition holds eventually for fixed \(T\).
The remaining product has logarithm at least
\((7/8)(5/6)k_0\log M>(2/3)k_0\log M\), a contradiction.

For this prime the scaled differences between \(jy\) and
every other literal ramp factor, scaled by \(c_i\), are at
most \(JWX\) in absolute value. Every nonzero such difference
has magnitude below \(p\). Hence the only prime hit in any
moving row at any order \(2\le j<J\) is its exact aligned
factor \(jy/c_i\), if present. It is in the first block, and
no short type supplies one. Every core factor and order
factorial is a unit once \(X\) is large. At each order, up to
the common positive mark depth, the valuation array has form
\[
                 R_-,\ R_-+R_+,\ldots,R_-+R_+,\ R_+
\]
where \(R_-+R_+>0\). All its interior entries are strictly
above the endpoint chord. Independent argument rescalings
make the normalized reductions \(1+z^j\) at all these orders.
Since \(p>J\), the quotients of distinct roots at each order
reduce to nontrivial roots of unity of that order. Consecutive
orders are coprime, so these quotient residues cannot agree.
A pair of nonreal roots with a common argument would give a
common quotient through complex conjugation, which is
algebraic and can be compared at the chosen prime. Even
positivity excludes real physical common roots. The first
pair is covered by even positivity and the higher pairs by
the absolute angular cutoff. Every size threshold used in
this argument is computable from the fixed data.

## All-degree closure with paired small cofactors

Fix an ordinary core and a constant \(B\ge1\). At each
sufficiently large \(L\) consider two opposed types with
\(B^{-1}L\le Y\le BL\), counts \(r\ge s\ge1\), and gaps
\(d,e\ge1\). There is no comparability hypothesis between
the gaps. For each adjacent pair \(2\le n<N=n+1<J\) suppose
one can choose \(p_n,p_N\) in (1) with
\[
                       A_n A_N\, e < \frac{L}{B+4}.           \tag{4}
\]
Then these profiles give all-degree coprimality above a
computable \(L\) depending only on \(B,J,C\). This is an
exact sufficient condition on the selected height and the
gaps, with no initial affine relation between \(L,Y\).
The same condition (4) at one pair excludes that pair alone.

Indeed at a collision \(U_m=mY+t_m\) satisfies (2), so
\(U_m/F_m\) has reduced denominator dividing \(A_m\).
The two ratios have difference at most \((B+4)e/L\), since
\[
 \frac{nY+t_n}{nL+1}-\frac{NY+t_N}{NL+1}
       = \frac{-Y + (Nt_n-nt_N)L+t_n-t_N}{(nL+1)(NL+1)} .
\]
If they were unequal their difference in absolute value
would be at least \(1/(A_n A_N)\), contradicting (3).
Thus their common ratio and the affine relation are
\[
                   k=Nt_n-nt_N,\qquad Y=kL+b_1,\qquad
                                      b_1=t_n-t_N .            \tag{5}
\]
The positive integer \(k\) is bounded solely in terms of
\(B\) at a collision for large \(L\). Indeed \(k=U_n/F_n\le
B+e/L\), and curvature below twelve with \(Y\le BL\)
gives \(e\le12+\sqrt{12BL}\). The intercept satisfies
\(|b_1|\le Je\). Multipliers \(k,1\) align the two intervals
within \(O_{J,B}(X)\), where \(X=\max(kd,e)\). For unbounded
\(X\) apply the proved clustered interval lemma with two
types and the fixed core. For bounded \(X\), there are
finitely many choices for \(k,b_1,d,e\), so (5) belongs to a
finite list of opposed affine families. The
[all-degree opposed affine theorem](MULTIVARIATE_OPPOSED_INDEPENDENT_HEIGHT_RATIONAL_BIAS.md)
excludes them above a uniform computable height, regardless
of the two counts.

For example, fix \(n,N\) and choose \(1/2<\theta<\exp(-1/2)\).
Among \(L\in[H,2H]\), the lower proportion for which both
\(nL+1\) and \(NL+1\) have a prime factor exceeding
\(H^\theta\) is at least \(2\log(1/\theta)-1>0\) as
\(H\to\infty\). For each fixed \(m\), count the unique
possible prime above \(H^\theta\) dividing \(mL+1\).
For \(H^\theta<p\le2mH+1\) the number of eligible \(L\)
is \(H/p+O(1)\), and \(p\nmid m\) for large \(H\).
Partial summation of the prime number theorem gives
\(\sum_{H^\theta<p\le2mH+1}1/p=\log(1/\theta)+o(1)\),
as in the sources for the
[native prime density theorem](MULTIVARIATE_NATIVE_PRIME_GROWING_PROFILES.md#native-witness-and-prime-inputs).
The cumulative floor error is \(O(H/\log H)\).
Uniqueness holds since \(H^{2\theta}>2mH+1\) for large \(H\).
The lower bound for the intersection follows by the union
bound on the two complements. These \(L\)'s certify the fixed
pair uniformly over all \(Y\) in the comparable band, all
counts as above, and all gaps with \(\max(d,e)\le L^\beta\)
for a fixed \(\beta<2\theta-1\). The winning primes exceed
\(K\) for large \(L\), so the cofactor loss is zero and
\(A_n A_N=O(L^{2(1-\theta)})\). This density statement is for
each fixed pair. The sufficient condition (4) for all pairs
is checked on the actual integer \(L\).

## One dominant mark near a positive diffuse contact

For fixed \(d,e\), let a proposed sequence of actual collisions
at a fixed pair with \(L,Y\to\infty\) and \(r\ge s\ge1\) tend
to a positive physical Gaussian reference of the fixed core.
Assume its limiting curvature share for each of the two types
is strictly positive. Then there is a computable
\(\delta_0>0\) such that no such sequence can have
\[
                          A_m\le L^{\delta_0}                  \tag{6}
\]
at even one of the two orders, allowing the choice of order
and prime to vary. Simplicity and angular transversality are
unnecessary. In this part of the theorem the gaps stay fixed.

To prove this, take a subsequence with fixed signed tail gap
\(h=rd-se\) by the
[signed imbalance bound](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md).
Write \(x=rd^2/(L+d)\), \(y=se^2/(Y+e)\), so the compactified
parameters are \(\tau=x+y\) and \(\lambda=x/d-y/e\).
Since the two limiting shares are positive and \(r/s\to e/d\),
the heights are comparable and the largest relative gap is
\(O(1/L)\). The
[effective contact-distance estimate](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#from-coefficient-error-to-distance-to-the-finite-set)
gives \((\tau,\lambda)=(\tau_*,\lambda_*)+O(L^{-\gamma})\)
for a computable \(\gamma>0\), since the contact set is finite
and the sequence is converging to the chosen contact.
In particular
\[
          \frac{Y}{L}=\rho+O(L^{-\gamma}),\qquad
             \rho=\frac{\tau_*+e\lambda_*}{\tau_*-d\lambda_*}.
\]
One may decrease \(\gamma\) to at most one. The two positive
shares make \(\rho>0\). At the contact
\(\tau_*=-2nN\log q_*\), \(\lambda_*=-2nN\log v_*\)
with \(q_*,v_*>0\) algebraic. Thus \(\rho=\log\alpha/\log\beta\)
for \(\alpha=q_* v_*^e\) and \(\beta=q_* v_*^{-d}\),
both strictly between zero and one.

Pass to a fixed \(m\) if necessary, and apply (2).
Writing \(U_m/P_m=B_m\) with \(B_m\) a positive integer gives
\[
                      \left|B_m/A_m-\rho\right|
                                    =O(L^{-\gamma}).
\]
If the fraction differs from \(\rho\),
[Matveev's Corollary 2.3](https://www.mathnet.ru/links/11696330cc866d273c92c771a63da235/im314_eng.pdf)
applied to the two fixed real logarithms gives
\(\left|B_m/A_m-\rho\right|\ge c_* A_m^{-\kappa}\)
for computable \(c_*,\kappa>0\), since \(B_m=O(A_m)\).
This contradicts (6) for \(\delta_0<\gamma/\kappa\)
when \(L\) is large. Otherwise the fraction equals \(\rho\)
exactly, making \(\rho\) rational and forcing
\(Y=\rho L+(\rho-t_m)/m\). There are finitely many offsets
at the two orders. Each resulting affine ray is excluded
by the opposed affine theorem. This proves the claim.

The exact paired criterion still requires sufficiently
concentrated prime-power parts of the two first marks.
The one-mark criterion lowers the number of required
marks at positive strength, with a smaller explicit
exponent in the cofactor bound. These conditions leave the
general case of dispersed prime factors to a separate
argument.
