# Prime terminals for multivariate unit-gap supports

Let \(a,b\) be distinct nonnegative integer exponent vectors. Remove
the equal rows, which contribute nonzero scalars to the ordinary
factorial moments. Write
\[
 M_j(t)=\sum_{k=0}^j\binom jk
             \prod_i((j-k)a_i+kb_i)!\,t^k .
\]
When all unequal rows increase by one, or all decrease by one,
\(M_{N-1}\) and \(M_N\) are coprime at every odd prime terminal
\(N=p\) and at \(N=p+1\). There is no bound on the dimension or on
the bases of the exponent rows. The profiles can vary with \(p\).

A separate prime-terminal criterion below allows both orientations
in equal numbers. It uses a residual polynomial given by an
Eulerian polynomial of degree determined by the number of rows.

All proofs use Newton polygons in the constant-to-leading
coefficient orientation. In this orientation, factors use translates
of edges of their products' polygons. At a slope of reduced denominator
\(e\), the width contributed by a rational factor is a multiple of \(e\).
These are the factorization properties used in
[Rocks, §2.2](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf)
and in the [ordinary univariate prime-power case of Filaseta and Rocks](https://people.math.sc.edu/filaseta/papers/ColloquiumMathPaper2016.pdf).
In addition, \(M_j\) at even \(j\) is strictly positive on the real
axis by the Gamma integral for the ordinary factorial functional.
A common factor over \(\mathbb Q\) therefore has even degree. The
[endpoint coefficient comparison](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients)
excludes \(M_{j-1}\mid M_j\).

## Rows with the same unit-gap orientation

Reverse both polynomials if necessary to assume \(b_i=a_i+1\) in
every retained row. Let \(r\ge1\) be their number. At the upper order
\(p\), Legendre's formula gives
\[
 v_p((a_i p+k)!)=a_i+v_p(a_i!)\quad(0\le k<p),
\]
and the increase from \(k=0\) to \(k=p\) is \(1+v_p(a_i+1)\).
The binomial coefficients at all interior indices have valuation one.
Thus after subtracting a common baseline, the ordinate at \(k=0\)
is zero, those at \(1\le k\le p-1\) are one, and the last ordinate is
\[
 H=\sum_i (1+v_p(a_i+1)).
\]
If \(H=1\), the lower polygon is one edge of slope \(1/p\), and
\(M_p\) is irreducible over \(\mathbb Q\). This excludes a common
factor with \(M_{p-1}\) by degree and with \(M_{p+1}\) by
nondivisibility. If \(H\ge2\), the two lower slopes of \(M_p\) are
\(1/(p-1)\) of width \(p-1\), and \(H-1\) of width one. At terminal
\(p\), a common factor with the first slope must use width \(p-1\).
It would be the full lower-order polynomial. A factor using only
the second slope has degree at most one, which parity excludes.

For terminal \(p+1\), first take \(r\ge2\).
For a row write \(a_i=p q_i+s_i\) with \(0\le s_i<p\). The floor of
its factorial argument divided by \(p\) is
\[
 a_i+q_i+\left\lfloor\frac{s_i+k}{p}\right\rfloor
 \qquad (0\le k\le p+1).
\]
The factorial valuation is \(l+v_p(l!)\) at quotient \(l\), so each
increase of the floor increases the valuation by at least one,
regardless of the size of \(a_i\). The binomial valuation at this
order is one for \(2\le k\le p-1\) and zero elsewhere.

If no \(s_i=p-1\), the minimum plateau is \([0,1]\).
Relative to its ordinate, the rises at indices \(2,\ldots,p-1\)
are at least one and those at \(p,p+1\) are at least \(r\).
The first positive edge starts at index one, so its slope is at
least \(\min\{1/(p-2),r/p\}>1/(p-1)\).
If some \(s_i=p-1\), the minimum occurs only at zero.
The rises at \(k=1\), \(2\le k\le p-1\), \(k=p\), and \(k=p+1\)
are at least \(1,2,r,r+1\), respectively. The first positive slope
is therefore at least
\(\min\{1,2/(p-1),r/p,(r+1)/(p+1)\}>1/(p-1)\).
In either case there is no edge of slope \(1/(p-1)\) in \(M_{p+1}\).
The remaining slope of \(M_p\) has width one and parity again
excludes a common factor.

When \(r=1\) and \(H\ge2\), the base is congruent to \(-1\)
modulo \(p\). The two floor jumps in \(M_{p+1}\) occur at
\(k=1,p+1\). Write \(J_1,J_2\) for their positive increases
in the factorial valuation. Relative to the initial ordinate,
the record points are \((0,0),(p,J_1),(p+1,J_1+J_2)\).
All other points lie above the first connecting segment.
If \(J_1/p\le J_2\), the two candidate slopes are
\(J_1/p\) and \(J_2\), neither of which is \(1/(p-1)\).
Otherwise the only slope is \((J_1+J_2)/(p+1)\), which is
at least \(J_2\ge1\). Parity applies as above.

## Equal numbers of increasing and decreasing rows

Let \(p\) be an odd prime. Suppose there are \(r\ge2\) rows
\(a_i\to a_i+1\) and \(r\) rows \(b_i+1\to b_i\), with all lower
bases \(a_i,b_i\ge1\). Assume
\[
 \min_i a_i=\min_i b_i=1,\qquad
 \max_i a_i+\max_i b_i\le p .
\]
Put \(j=(r-1)\bmod (p-1)\), taking \(0\le j\le p-2\).
For \(j\ge1\), define the Eulerian polynomials and their resultants by
\[
 E_1(x)=1,\qquad
 E_{j+1}(x)=(1+jx)E_j(x)+x(1-x)E_j'(x),\qquad
 R_j=\operatorname{Res}_x(E_j(x),E_j(-x)).
\]
Then \(M_{p-1}\) and \(M_p\) are coprime if \(j=0\), or if
\(j\ge1\) and \(p\nmid R_j\). In particular the hypothesis on the
prime is satisfied for \(j=0,1,2,3\), and for any \(j\ge1\) with
\[
 p>(j!)^{2(j-1)}.
\]
There is an exact weaker sufficient condition for \(j\ge1\).
It is enough that \(E_j(x)\) and \(1+x^{p-1}\) be coprime modulo \(p\).

To prove these statements, put \(L=\sum_i(a_i+b_i)\). At order
\(p-1\) the floor of a positive row's factorial argument divided
by \(p\) drops from \(a_i\) to \(a_i-1\) when \(k<a_i\).
That of a negative row drops from \(b_i\) to \(b_i-1\) when
\(k\ge p-b_i\). The two sets of drop intervals are disjoint by the
bound on the maxima. Every factorial argument at this order is
below \(p^2\), and the binomial coefficients are all units. Thus
the lowest ordinate is \(L-r\), attained exactly at the two
endpoints by the assumptions on the minima. For every exponent
\(1\le e\le p\), writing \(e(p-1)=(e-1)p+(p-e)\) gives the endpoint
factorial unit
\[
 p^{-v_p((e(p-1))!)}(e(p-1))!
 \equiv (-1)^{e-1}(e-1)!(p-e)!
 \equiv -1\pmod p
\]
by Wilson's theorem. There are \(2r\) factors at each endpoint,
so the normalized \(M_{p-1}\) reduces to \(1+t^{p-1}\).

At order \(p\), all the interior ordinates are \(L+1\)
and the two endpoints are both higher. At index \(k\) with
\(1\le k\le p-1\), the binomial unit is
\(\binom pk/p\equiv(-1)^{k-1}/k\).
After common full-block units are removed, each increasing row
contributes \(k!\), and each decreasing row contributes
\((p-k)!\equiv(-1)^k k/k!\). With \(x=(-1)^{r+1}t\), the
normalized reduction of \(M_p\) is therefore a nonzero scalar
times
\[
 \sum_{k=1}^{p-1} k^{r-1}x^k .
\]
A common rational factor gives a common factor of the normalized
reductions of positive degree and with nonzero constant term.
Indeed, the \(p\)-adic coefficient minimum is additive under
polynomial products, and both endpoint units of the normalized
\(M_{p-1}\) are nonzero. The normalized common factor therefore
retains its degree and its constant coefficient after reduction.
It suffices to exclude a common root \(x\ne0\) satisfying
\(x^{p-1}=-1\).

For \(j=0\), the sum is \((1-x)^{p-1}-1\).
At such a root \(x^p=-x\) and \(x\ne1\), so it equals
\((1+x)/(1-x)-1=2x/(1-x)\ne0\).
For \(j\ge1\), apply \((x\,d/dx)^j\) to the geometric sum.
In characteristic \(p\) this derivative kills \(1-x^p\),
so the result is \(x(1-x)^{p-1-j}E_j(x)\).
Hence a common root would satisfy \(E_j(x)=0\); taking the
\(p\)-th power gives \(E_j(-x)=0\). The recurrence above gives
strictly negative, simple real roots of \(E_j\) by induction:
at a root \(\alpha\) of \(E_j\) the new sign is that of
\(\alpha(1-\alpha)E_j'(\alpha)\). These signs alternate.
Together with the signs at zero and at negative infinity, they
force one root in each intervening interval and one on each
outer interval in the negative half-line. Thus \(R_j\ne0\)
over \(\mathbb Z\). Each row of its Sylvester matrix has
Euclidean norm at most \(E_j(1)=j!\), giving the stated bound by
Hadamard's inequality. For \(j=1,2,3\) the polynomials are
\(1,1+x,1+4x+x^2\), which have no nonzero opposite root pair
modulo an odd prime.

The arithmetic criterion concerns actual factorial coefficients and
does not require the roots of the normalized complex polynomials to
have different arguments under independent positive rescaling. In
both families the rational Bézout identity transfers coprimality
to any characteristic-zero coefficient field.
