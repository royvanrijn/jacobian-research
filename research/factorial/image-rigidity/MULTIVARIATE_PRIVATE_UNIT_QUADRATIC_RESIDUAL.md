# The quadratic residual case for two private units

The [all-degree both-private theorem](MULTIVARIATE_BOTH_PRIVATE_ALL_DEGREES.md)
now excludes physical common roots in this class. The reduction below
records the earlier necessary arithmetic conditions.

Consider an ordinary two-monomial profile with one private unit row
on each side and arbitrary shared rows. Define
\[
    V=\sum_{\rm shared} d_i^2/\ell_i,\qquad
    S=\sum_{\rm shared} d_i^2/(\ell_i+|d_i|),\qquad
    \Delta=\sum_{\rm shared}d_i ,
\]
where \(d_i=b_i-a_i\) and \(\ell_i=\min(a_i,b_i)\).

**Theorem.** A failure of physical SFC(2) in this class must
satisfy
\[
   13\le n\le2047,\qquad 1\le|\Delta|\le2,\qquad
                 V>2/3,\qquad S<15/8 .
\]
The rational gcd of the physical moments at \(n,N=n+1\)
would be a single irreducible quadratic \(z^2-u z+c\) up to
a scalar. Its roots have angles
\[
                       \pm\theta,\qquad
                 \frac{\pi}{n}<\theta<\frac{2\pi}{n+2}.
\]
Its two roots must have distinct integer valuations at some
common Newton slope prime. More quantitatively, the denominator
in lowest terms of the rational number \(4-u^2/c\) is precisely
\[
         D=\prod_{\{p:\ s_p\ne t_p\}}p^{|s_p-t_p|},
         \qquad D>\frac{(n+2)^2}{160},                        \tag{1}
\]
where \(s_p,t_p\) are the two valuations of the quadratic roots.
The reduction uses neither the unevaluated absolute cutoff nor
any additional fixed-terminal scan.

## Reduction to the first angle box

The
[terminal-thirteen theorem](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md)
and its predecessors cover \(n\le12\). The
[explicit first-label cutoff](MULTIVARIATE_PRIVATE_UNIT_EXPLICIT_FIRST_LABEL_CUTOFF.md)
and its companion higher-label proof cover \(n\ge2048\).
The
[first-box signed imbalance theorem](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_THREE.md)
covers \(|\Delta|\ge3\), and the
[balanced Euclidean theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)
covers \(\Delta=0\). The
[absolute weak bound](MULTIVARIATE_PRIVATE_UNIT_FIRST_LABEL_TWO_THIRDS.md)
and the
[high-curvature bound](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
give the two remaining curvature restrictions.
If \(|\Delta|=2\), the [signed imbalance two cutoff](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_TWO_ORDER_CUTOFF.md)
also requires \(n\le122\). Orient the columns so that \(\Delta=2\).
The first-box signed imbalance theorem further requires at
least one decreasing shared row with absolute gap
\(d>\ell/2\). Such rows have \(d\le5,\ \ell<2d\) and occur
at most three times.

The [higher-label exclusion](MULTIVARIATE_PRIVATE_UNIT_ALL_HIGHER_LABELS.md)
covers every \(k\ge2\) for \(n\ge13\). The
[geometric angle boxes](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
give only one root in the first positive box of each degree,
and its conjugate in the negative box. These roots are simple.
Real arguments were already excluded by even-order positivity
and positive coefficients. Thus a common factor has degree
at most two, and a nonreal common root gives exactly the
stated quadratic, since the coefficients of the moments
are rational.

## Exact denominator and split-prime constraint

The two roots are \(\alpha,\bar\alpha\) with
\(u=\alpha+\bar\alpha>0\), \(c=\alpha\bar\alpha>0\).
Hence
\[
                   4-u^2/c=4\sin^2\theta
                         <16\pi^2/(n+2)^2
                         <160/(n+2)^2 .                      \tag{2}
\]
At any prime, if the quadratic root valuations differ,
they are integers. In fact an irreducible quadratic over
\(\mathbb Q_p\) would give conjugate roots of equal valuation,
by uniqueness of the valuation extension. At unequal
valuations, the nonarchimedean triangle inequality gives
\(v_p(u)=\min(s_p,t_p)\), so \(v_p(u^2/c)=-|s_p-t_p|\).
Subtracting from the integral number 4 preserves this
negative valuation. At equal valuations \(v_p(u^2/c)\ge0\).
This proves the equality for the denominator in (1).
A positive rational of denominator \(D\) is at least \(1/D\),
so (2) proves its strict lower bound. At \(n\ge13\) the
bound exceeds one, requiring a split prime with unequal
valuations.

For an exact sieve one can bound each spread \(|s_p-t_p|\)
by the range of integer slopes common to the two Newton
polygons. Only finitely many primes can contribute, since
if either primitive moment array has unit endpoints at \(p\),
all its root valuations are zero. The
[quadratic trace and Newton sieve](MULTIVARIATE_QUADRATIC_NEWTON_TRACE_SIEVE.md)
then gives candidates on a rational grid in the narrower
angle interval above. In addition, for an admissible trace
candidate \(\delta=4-u^2/c\), the positive rational
\(\delta/(4-\delta)\) must belong to
\(|\Delta_K|(\mathbb Q^\times)^2\), where \(\Delta_K\) is an
allowed negative quadratic field discriminant from the
[field discriminant ledger](MULTIVARIATE_COMMON_FACTOR_FIELD_DISCRIMINANT.md).
Indeed
\(u^2-4c=-u^2\delta/(4-\delta)\), so the quadratic
discriminant has just that square class. These are necessary
conditions on the original physical moments.
