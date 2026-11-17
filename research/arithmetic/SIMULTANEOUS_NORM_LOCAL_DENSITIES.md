# Simultaneous Gaussian–Eisenstein sums: prime-power local densities

The three-quadratic-equation model has summable forbidden-prime losses even
after singular solutions are lifted. This note proves that local statement,
including a uniform bound when a forbidden prime divides the target. It does
not count integer points or prove the proposed binary coverage theorem.

Write

\[
 G(x,y)=x^2+y^2,\qquad E(u,v)=u^2-uv+v^2,
\]

and consider the system

\[
 G(x,y)=E(u,v)=a,\qquad G(z,w)=E(s,t)=b,\qquad a+b=N. \tag{1}
\]

There are eight coordinates and three equations after eliminating \(a,b\).
For a prime \(p\) and \(k\ge1\), let \(R_{p,k}(N)\) count its solutions
modulo \(p^k\), and let \(B_{p,k}(N)\) count those with \(p\mid ab\).
These are counts of coordinate representations, with multiplicity.

## The local theorem

**Theorem L.** For every positive integer \(N\) and every prime
\(p\equiv5,7,11\pmod {12}\), the limits

\[
 \alpha_p(N)=\lim_{k\to\infty}p^{-5k}R_{p,k}(N)>0,
 \qquad
 \beta_p(N)=\lim_{k\to\infty}\frac{B_{p,k}(N)}{R_{p,k}(N)}
\]

exist, and

\[
 0\le\beta_p(N)<\frac9{p^2}. \tag{2}
\]

The bound is uniform in \(N\), including its valuation at \(p\). In
particular,

\[
 \sum_{p\equiv5,7,11\ (12)}\beta_p(N)
 <9\sum_{\substack{n\ge5\\(n,6)=1}}n^{-2}
 =\pi^2-9<1. \tag{3}
\]

For \(p\nmid N\), the exact fractions are

\[
 \beta_p(N)=
 \begin{cases}
 \displaystyle\frac{2(3p^2-1)}{p(p^3-p^2+3p+1)},&p\equiv5,7\pmod {12},\\[6pt]
 \displaystyle\frac2{p(p-1)},&p\equiv11\pmod {12}.
 \end{cases} \tag{4}
\]

When \(v_p(N)=1\), the bad fraction is zero. Formula (3) is a strict local
margin. Transferring it to a count of global integral representations requires
additional theorems, including control of primes growing with \(N\).

## Why these are the relevant prime restrictions

Let

\[
 C=\{c\ge1:p\mid c\Longrightarrow p\equiv1\pmod {12}\},
\]

with \(1\in C\). A positive integer lies in \(C\) if and only if it has a
primitive representation by each of \(G\) and \(E\). Here primitive means
that the two coordinates have greatest common divisor one; the Gaussian and
Eisenstein coordinates can be chosen separately.

Indeed, a prime \(3\pmod4\) dividing a Gaussian norm divides both its
coordinates. A factor \(4\) has the same consequence at \(2\). Thus a
primitive Gaussian norm has prime factors \(1\pmod4\), with at most one
additional factor \(2\). Similarly, a primitive Eisenstein norm has prime
factors \(1\pmod3\), with at most one additional factor \(3\). The latter
statement follows by reduction modulo an inert prime, and by the identity
that the ramified prime above \(3\) has square associated to \(3\).
Intersecting these two descriptions removes both exceptional factors and
leaves exactly the primes \(1\pmod {12}\).

For completeness, the converses use factorization rather than an assumption
about a particular representation. Both \(\mathbb Z[i]\) and
\(\mathbb Z[\omega]\), \(\omega^2+\omega+1=0\), are norm-Euclidean:
rounding to the nearest point of the square or triangular lattice leaves
norm at most \(1/2\) or \(1/3\), respectively. They are consequently unique
factorization domains. A prime \(1\pmod {12}\) splits into two conjugate
prime factors in each ring. This follows from the roots of \(T^2+1\) and
\(T^2+T+1\) modulo the prime; their corresponding prime ideals have norm
\(p\) and are principal. For each rational prime dividing \(c\), choose
one conjugate prime factor and use its full required exponent. The resulting
element has norm \(c\). A common rational prime divisor of its coordinates
would require both conjugate factors, which were not both chosen. This gives
primitive representations in both rings.

There is also a useful automatic exclusion in (1). A common norm of \(G\)
and \(E\) is \(0\) or \(1\pmod4\): an even Eisenstein norm is divisible
by \(4\), whereas a Gaussian norm is never \(3\pmod4\). An Eisenstein
norm is \(0\) or \(1\pmod3\). Consequently, if \(N\equiv2\pmod {12}\),
every integral solution of (1) has

\[
 a\equiv b\equiv1\pmod {12}. \tag{5}
\]

In particular, neither \(2\) nor \(3\) divides either summand, and both
summands are positive. For \(p\ge5\) outside \(1\pmod {12}\), a component
value divisible by \(p\) has even positive valuation, since at least one
of the two norm forms is anisotropic at \(p\). Excluding \(p\mid ab\)
therefore enforces the intended support condition on the *values*. Coordinates
at allowed split primes need not themselves be primitive.

These necessary small-prime conditions are locally attainable. Take
\(a=1\), \(b=N-1\). Over \(\mathbb Q_3\), the unit \(b\equiv1\pmod3\)
is a square, so both forms represent it with second coordinate zero.
Over \(\mathbb Q_2\), a unit \(b\equiv1\pmod4\) is Gaussian norm:
use \((\sqrt b,0)\) if \(b\equiv1\pmod8\), and
\((\sqrt{b-4},2)\) if \(b\equiv5\pmod8\). Every odd \(b\) is an
Eisenstein norm, since \(u^2-u+1=b\) has discriminant
\(4b-3\equiv1\pmod8\) and hence a square root in \(\mathbb Q_2\).
The resulting \(u=(1+\sqrt{4b-3})/2\) is integral. All the cited square
criteria follow by successive lifting of unit squares modulo powers of
the corresponding prime.

In fact, the total local densities at these two primes are constant on the
admissible target class:

\[
 \alpha_2(N)=\frac94,\qquad \alpha_3(N)=\frac{64}{27}
 \qquad(N\equiv2\pmod {12}). \tag{16}
\]

Here the definition of \(\alpha_p\) is the same normalized count as in
Theorem L, without restricting \(p\) to the forbidden odd primes.
For \(k\ge2\), an odd Gaussian norm \(a\equiv1\pmod4\) has
\(2^{k+1}\) representations modulo \(2^k\); other odd residues have
none. At \(k\ge3\), partition the even coordinate by its residue modulo
four. Exactly \(2^{k-2}\) even coordinates make the remaining odd square
congruent to one modulo eight, and each has four square roots. Reversing
the odd and even coordinates gives \(2^{k+1}\). The case \(k=2\) follows
directly. An odd Eisenstein norm has \(3\cdot2^{k-1}\) representations:
its three nonzero pairs modulo two are smooth, each lifting by a factor
two at each step. The common-value block count is therefore
\(3\cdot2^{2k}\). There are \(2^{k-2}\) possible values
\(a\equiv1\pmod4\), and their complementary values have the same
residue. Thus \(R_{2,k}(N)=9\cdot2^{5k-2}\).

At three, every unit Gaussian norm has \(4\cdot3^{k-1}\)
representations. A unit Eisenstein norm is \(1\pmod3\), and every such
unit has \(2\cdot3^k\) representations: for each \(v\pmod {3^k}\),
the equation for \(u\) has two simple roots, since its discriminant is
\(4a-3v^2\equiv1\pmod3\). The common-value block count is
\(8\cdot3^{2k-1}\). Summing over the \(3^{k-1}\) values
\(a\equiv1\pmod3\) gives \(R_{3,k}(N)=64\cdot3^{5k-3}\).
These exact identities prove (16), with no singular-lifting assumption
at two.

## Finite-field counts and singular classes

Fix \(p\ge5\), and put

\[
 \varepsilon=\left(\frac{-1}{p}\right),\quad
 \eta=\left(\frac{-3}{p}\right),\quad
 H=(p-\varepsilon)(p-\eta),\quad
 Z=(p+(p-1)\varepsilon)(p+(p-1)\eta).
\]

A nondegenerate binary norm with split sign \(\epsilon\in\{1,-1\}\)
has \(p-\epsilon\) representations of each nonzero residue, and
\(p+(p-1)\epsilon\) of zero. In the split case this is the count for a
nonzero scalar multiple of \(XY\). In the anisotropic case it follows
from the norm map \(\mathbb F_{p^2}^{\times}\to\mathbb F_p^{\times}\),
whose kernel has \(p+1\) elements, and whose only zero-norm vector is zero.

The number of representations of a fixed nonzero common value is therefore
\(H\), and that of the common value zero is \(Z\). Summing over \(a\)
gives

\[
 \begin{array}{c|cc}
 &R_{p,1}(N)&B_{p,1}(N)\\ \hline
 p\nmid N &(p-2)H^2+2HZ&2HZ\\
 p\mid N &(p-1)H^2+Z^2&Z^2.
 \end{array} \tag{6}
\]

In particular, for \(p\nmid N\) the raw finite-field bad fraction is
\((4p-2)/[p(p^2-2p+3)]\) at \(p\equiv5,7\pmod {12}\), and
\(2/[p(p^2-3)]\) at \(p\equiv11\pmod {12}\). These are not (4).

To locate the discrepancy, use the three equations
\(G-E=0\), \(G'-E'=0\), \(G+G'-N=0\). A linear relation with
coefficients \((\lambda,\mu,\nu)\) between their Jacobian rows must obey

\[
 (\lambda+\nu)\nabla G=0,\quad\lambda\nabla E=0,\quad
 (\mu+\nu)\nabla G'=0,\quad\mu\nabla E'=0. \tag{7}
\]

At \(p\ge5\), the gradient of each norm vanishes precisely at its zero
coordinate pair. Any three of the four resulting nonzero coefficient
vectors in (7) are linearly independent. Thus a point is singular exactly
when at most two of the four coordinate pairs are nonzero.

If \(p\nmid N\), the singular points are precisely those whose entire
\(a\)-block is zero, or whose entire \(b\)-block is zero. There are
\(2H\) such points. If \(p\mid N\), only the stratum \(a=b=0\) can be
singular. Put \(L_G=(p-1)(1+\varepsilon)\) and
\(L_E=(p-1)(1+\eta)\), the counts of nonzero zero-norm pairs. Its
singular-point count is

\[
 1+2L_G+2L_E+L_G^2+4L_GL_E+L_E^2. \tag{8}
\]

For \(p\equiv5,7\pmod {12}\), this is \((2p-1)^2\): every bad
finite-field point is singular when \(p\mid N\). For \(p\equiv11\pmod
{12}\), it is one. Applying smooth Hensel lifting to all of (6) would
therefore be incorrect.

## Counting at prime powers

Let \(t\in\mathbb Z_p\) be nonzero, with \(v=v_p(t)\). For a split
binary norm and \(k>v\), its number of representations modulo \(p^k\)
is

\[
 (v+1)(p-1)p^{k-1}. \tag{9}
\]

To see this, factor the form into two linear coordinates over \(\mathbb
Z_p\). In a product equal to \(t\), the first coordinate can have
valuation \(0,\ldots,v\). Each stratum contributes
\((p-1)p^{k-1}\) solutions. For an anisotropic binary norm, valuation is
twice the minimum coordinate valuation. The count is consequently zero
for odd \(v\), and

\[
 (p+1)p^{k-1} \tag{10}
\]

for even \(v\), obtained by scaling out \(p^{v/2}\) and lifting the
unit-norm equation. At the zero residue modulo \(p^k\), the respective
counts are

\[
 [p+k(p-1)]p^{k-1}\quad\hbox{and}\quad p^{2\lfloor k/2\rfloor}. \tag{11}
\]

The first formula partitions the valuation of the first factor, including
the zero factor; the second requires both coordinates to be divisible by
\(p^{\lceil k/2\rceil}\). These formulas also control the singular lifts.

For a forbidden prime define

\[
 (h,d)=
 \begin{cases}
 (1-p^{-2},1),&p\equiv5,7\pmod {12},\\
 ((1+p^{-1})^2,0),&p\equiv11\pmod {12},
 \end{cases}
 \qquad
 f(v)=\begin{cases}(v+1)^d,&v\text{ even},\\0,&v\text{ odd}.\end{cases}
\]

By (9)–(10), the limiting four-coordinate density for a fixed nonzero
common value \(t\) is \(w_p(t)=h f(v_p(t))\). Haar measure below is
normalized by \(\operatorname{vol}(\mathbb Z_p)=1\). It follows that

\[
 \alpha_p(N)=h^2 I_r,\qquad
 I_r=\int_{\mathbb Z_p}f(v_p(a))f(v_p(N-a))\,da,
 \qquad r=v_p(N). \tag{12}
\]

The value assigned to \(f(v_p(0))\) is immaterial. To justify the limit
in (12), use the exact formulas (9)–(11) in the convolution over residues
\(a\pmod {p^k}\). Away from \(a=0,N\), its terms agree exactly with the
integrand on the corresponding residue balls once \(k>r\). On either
exceptional ball, (11), divided by \(p^{2k}\), is at most \(k+1\) for a
forbidden prime. The other component has the fixed valuation \(r\).
Their normalized contribution is \(O_{p,r}(k p^{-k})\), tending to zero.
The integral tails tend to zero as well: \(f(v)\) grows at most linearly
and the stratum of valuation \(v\) has measure \((1-p^{-1})p^{-v}\).
This proves the limit for both total and bad counts, including singular
classes. In particular, \(I_r>0\) because the unit–unit stratum has
positive measure.

An exact useful valuation formula is

\[
\begin{split}
 I_r={}&(1-p^{-1})\sum_{\substack{0\le e<r\\e\ {\rm even}}}
             (e+1)^{2d}p^{-e}\\
 &+\mathbf1_{2\mid r}\left[
       (1-2p^{-1})(r+1)^{2d}p^{-r}
       +2(1-p^{-1})(r+1)^d
          \sum_{\substack{e>r\\e\ {\rm even}}}(e+1)^d p^{-e}
                      \right]. \tag{13}
\end{split}
\]

For \(v_p(a)<r\), the two component valuations are equal. For
\(v_p(a)>r\), the other valuation is \(r\); the symmetric case gives the
factor two. When both valuations equal \(r\), their leading unit residues
must avoid zero and the leading residue of \(N/p^r\), giving
\(1-2/p\). This proves (13), including odd \(r\).

The unit–unit contribution to \(I_r\) is \(g_r=1-2/p\) for \(r=0\)
and \(g_r=1-1/p\) for \(r\ge1\). Thus

\[
 \beta_p(N)=1-\frac{g_r}{I_r}. \tag{14}
\]

When \(r=0\), summing the geometric tails in (13) gives exactly (4).
For example, the sum used for \(d=1\) is
\(\sum_{j\ge1}(2j+1)q^j=q(3-q)/(1-q)^2\), with \(q=p^{-2}\).
When \(r=1\), (13) gives \(I_1=1-1/p\), so the bad fraction is zero.
Equivalently, both positive component valuations would have to be even,
forcing \(p^2\mid N\).

There is a second direct way to see the singular correction in (4). Let
\(V_k\) count \(G-E=0\) modulo \(p^k\), and set
\(\kappa=\varepsilon\eta\). Then

\[
 V_0=1,\quad V_1=p^3+\kappa p(p-1),\quad
 V_k=(V_1-1)p^{3(k-1)}+p^4V_{k-2}\quad(k\ge2).
\]

Every nonzero solution modulo \(p\) is smooth; the remaining solutions
have all four coordinates divisible by \(p\). Consequently,

\[
 \nu=\lim p^{-3k}V_k=
 \begin{cases}(p+1)/p,&\kappa=1,\\
 (p^2+1)/(p(p+1)),&\kappa=-1.\end{cases}
\]

For \(p\nmid N\), let \(S_k\) count the first-block solutions of
\(G-E=0\) modulo \(p^k\) with \(p\mid G\). Its \(Z-1\) nonzero
reductions are smooth, and its zero reduction scales by \(p\), so
\(S_k=(Z-1)p^{3(k-1)}+p^4V_{k-2}\) for \(k\ge2\). Thus
\(p^3\lim p^{-3k}S_k=Z-1+p\nu\). The other block has unit common
value. The effective zero-value weight in the ratio (6) consequently
becomes \(Z_* = Z-1+p\nu\), rather than \(Z\). Hence the limiting bad
fraction is \(2Z_*/[(p-2)H+2Z_*]\), again (4). At
\(p\equiv11\pmod {12}\), this changes an apparent \(O(p^{-3})\)
finite-field loss into \(O(p^{-2})\). For the concrete example \(p=11\),
\(11\nmid N\), the raw fraction is \(1/649\), whereas its prime-power
limit is \(1/55\).

## Uniform control when the prime divides the target

Suppose \(r\ge1\). The function \(a\mapsto f(v_p(a))\) lies in
\(L^2(\mathbb Z_p)\), and reflection followed by translation preserves
Haar measure. Cauchy–Schwarz in (12) gives

\[
 I_r\le(1-p^{-1})S_d,\qquad
 S_d=\sum_{j\ge0}(2j+1)^{2d}p^{-2j}.
\]

Together with (14), this gives \(\beta_p(N)\le1-S_d^{-1}\).
If \(d=0\), then \(S_0=(1-q)^{-1}\), so the bound is \(q\).
If \(d=1\), then

\[
 S_1=\frac{1+6q+q^2}{(1-q)^3},\qquad
 \beta_p(N)\le\frac{9q-2q^2+q^3}{1+6q+q^2}<9q. \tag{15}
\]

For \(r=0\), formula (4) also gives a bound below \(9/p^2\). For its
first case one may use the simpler upper bound
\(6/[p(p-1)]\le(15/2)p^{-2}\), since \(p\ge5\); the second case is
smaller still. This proves (2) for every valuation. Finally, all forbidden
primes are at least five and coprime to six. Removing multiples of two
and three from \(\sum n^{-2}\) proves (3). Its right side is below one
already from \(\pi<22/7\), since \((22/7)^2<10\).

An entirely rational bound, without evaluating \(\pi\), follows by an
integral test after the terms \(5,7,11,13\):

\[
 9\sum_{\substack{n\ge5\\(n,6)=1}}n^{-2}
 \le9\left(\frac1{25}+\frac1{49}+\frac1{121}+\frac1{169}
                +\frac1{66}+\frac1{78}\right)
 =\frac{23122584}{25050025}<1.
\]

This supplies a second check on the strict margin.

## Exact computational checks and limits of the result

[The standard-library checker](verify_simultaneous_norm_local.py) enumerates
the two binary norm histograms modulo each of 18 prime powers: \(5^k\)
through \(k=4\); \(7^k,11^k\) through \(k=3\); and
\(13^k,17^k,19^k,23^k\) through \(k=2\). For every represented residue it
checks (9)–(11), and convolves their products to count (1). A separately
implemented valuation-stratum formula checks the complete total and bad
counts for unit, valuation-one, valuation-two and zero target residues
where those classes are distinct. Support-stratum Jacobian counts check
(7)–(8). The allowed class \(p=13\) provides an additional control.

The checker also tests (4), (13)–(15) using exact fractions, checks
prime-power convergence at exponent 20 for 30 combinations of prime and
target valuation, and independently enumerates primitive norms through
2000 against trial-factorization membership in \(C\). Three deliberately
false formulas must be rejected, including treating the raw \(p=11\)
fraction as its limiting density. Explicit exceptions implement all
checks, so running Python with optimization does not disable them.
Ten further binary-histogram prime powers verify (16), checking every
admissible target residue at each level: \(2^k\) for \(2\le k\le7\),
and \(3^k\) for \(1\le k\le4\).

Replay with

```text
python3 research/arithmetic/verify_simultaneous_norm_local.py
```

The output is [the exact check receipt](simultaneous_norm_local_checks.json).
These are internal finite formula checks accompanying the written proof;
they are neither external review nor a proof of binary coverage.

The proved progress is local: singular lifting does preserve summability,
with a strict margin uniform in \(N\). A global counting theorem must
still show that balanced integer points on (1) occur with these local
densities, uniformly enough for the forbidden-prime sieve. Fixed-prime
equidistribution alone also leaves a separate large-prime tail. No such
global count is proved here. The [companion overview](SIMULTANEOUS_GAUSSIAN_EISENSTEIN_SUMS.md)
does prove an unconditional large-prime bound (Tail lemma T12); making that
bound negligible requires the still-missing global lower bound.
