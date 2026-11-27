# A uniform large-exponent corner for the initial homogeneous trinomial window

**Theorem.** There exists an absolute integer \(H\) such that, for
every homogeneous bivariate trinomial with distinct exponent columns
\((a_j,b_j)\), if some column satisfies \(a_j,b_j\ge H\), its first three
ordinary factorial moments have no common nonzero coefficient
direction. The total degree, both successive exponent gaps and all
coefficient ratios are unrestricted. Nonvanishing holds over
every characteristic-zero field.

This strengthens the fixed-offset large-height exclusion for the
initial homogeneous trinomial window. The argument below proves
existence of a uniform cutoff. It does not evaluate the cutoff.
A failing support must have every column within a bounded
exponent distance of one of the two axes.

## 1. The angular moment kernel

Write \(a_j+b_j=D\) and set
\[
 f_D(a)=a\log a+(D-a)\log(D-a),\qquad
 g_a(t)=e^{D\log D-f_D(a)}t^a(1-t)^{D-a},\quad 0<t<1.
\]
With \(t=x/(x+y)\), the radial integral \((mD+1)!\) is a
nonzero common factor at order \(m\). Hence the three factorial
equations are equivalent, after a common coefficient change,
to vanishing of the first three ordinary integrals on \([0,1]\)
of powers of a linear combination of the \(g_{a_j}\). Use
\(0\log 0=0\) for an endpoint exponent. The same angular
reduction occurs in [Liu and Sun, Theorem 2.6](https://doi.org/10.1017/S0004972719000546),
which assumes vanishing at all positive orders.

For \(m=1,2,3\) and real exponents \(0<a_i<D\), define
\[
 \begin{split}
 h(a)&=\frac{\sqrt{a(D-a)}}{D^{3/2}},\qquad
 \eta(a)=2\sqrt D\,\arcsin\sqrt{a/D},\\
 J(\boldsymbol a)&=\sum_{i=1}^m f_D(a_i)-m f_D(\bar a),\qquad
 K_{D,m}(\boldsymbol a)=\int_0^1\prod_{i=1}^m g_{a_i}(t)\,dt ,
 \qquad \bar a=m^{-1}\sum_i a_i .
 \end{split}
\]
The beta integral and Stirling's formula give the exact expression
\[
 K_{D,m}(\boldsymbol a)=\sqrt{2\pi/m}\,h(\bar a)
                   \exp\bigl(-J(\boldsymbol a)+r_{D,m}(\bar a)\bigr),
 \tag{1}
\]
where, writing \(\mathcal R\) for the remainder in
\(\log\Gamma(x+1)=(x+1/2)\log x-x+\tfrac12\log(2\pi)+\mathcal R(x)\),
\[
 r_{D,m}(a)=\mathcal R(ma)+\mathcal R(m(D-a))-\mathcal R(mD)
                      -\log(1+(mD)^{-1}).
\]
For every fixed \(k\ge0\), \(\mathcal R^{(k)}(x)=O_k(x^{-1-k})\)
on the positive ray for \(x\ge1\). This follows from the
[sectorial Stirling remainder](https://dlmf.nist.gov/5.11)
and Cauchy's derivative bound on a disk of radius \(x/2\).

### Local convergence with derivatives

Let \(c=c(D)\) with both \(c\) and \(D-c\) tending to infinity.
Use \(x_i=\eta(a_i)-\eta(c)\) in a fixed bounded box.
Then for each fixed derivative order,
\[
 \frac{K_{D,m}(\boldsymbol a)}{\sqrt{2\pi/m}\,h(c)}
      \longrightarrow
 {\cal K}_m(\boldsymbol x):=
          \exp\!\left(-\frac1{2m}\sum_{i<j}(x_i-x_j)^2\right)
 \tag{2}
\]
uniformly with those derivatives in the \(x_i\).

Indeed \(da/d\eta=s(a):=\sqrt{a(D-a)/D}\), and the
inverse \(a=D\sin^2(\eta/(2\sqrt D))\) gives bounded derivatives
of every fixed order at least two. In particular
\((a_i-c)/s(c)\to x_i\) uniformly with derivatives, since
\(s(c)\to\infty\).
On bounded \(u\), the function
\[
 f_D(c+s(c)u)-f_D(c)-s(c)u f_D'(c)
\]
tends to \(u^2/2\) with derivatives. Its exact second derivative
at zero is one. For \(k\ge3\), its \(k\)-th derivative is bounded
in absolute value by a constant times
\(\min(c,D-c)^{1-k/2}\), by \(f_D''(a)=1/a+1/(D-a)\).
The constant and linear terms cancel in \(J\).
Thus \(J\to\tfrac12(\sum_i x_i^2-m\bar x^2)\) with derivatives.
The other factors in (1) converge with derivatives as well.
In the variable \(u\), \(h(c+s(c)u)/h(c)\) is the square root of
\((1+s(c)u/c)(1-s(c)u/(D-c))\), which tends to one.
The derivative bounds for \(r_{D,m}\), composed with
\(c+s(c)\bar u\), give a vanishing error.

For divided differences in any argument, (2) also gives the
confluent limit when nodes coalesce. A divided difference on
real nodes is the average of the corresponding derivative on
their convex hull against simplex measure of total mass
\(1/k!\) at order \(k\). This follows by repeated integration
of the fundamental theorem of calculus. It bounds the
operation uniformly even at arbitrarily close nodes.

### Decay between separated clusters

Let \(H_*=\min_i\min(a_i,D-a_i)\ge1\) and
\(T=\max_{i,j}|\eta(a_i)-\eta(a_j)|\).
For \(m\le3\) and every subset \(E\) of the arguments,
\[
 \left|\prod_{i\in E}\partial_{\eta(a_i)}K_{D,m}(\boldsymbol a)\right|
 \ \le\
 C(1+T)^{|E|}
 \left(\prod_{i=1}^m h(a_i)^{1/m}\right)
                  \exp(-c_0T^2+C T/\sqrt{H_*}),
 \tag{3}
\]
for positive absolute constants \(C,c_0\). Each argument is
differentiated at most once in this estimate.

To prove the Gaussian decay, subtract the tangent to \(f_D\)
at \(\bar a\) in each summand of \(J\). The difference at \(a_i\)
is \(D\) times the relative entropy of the Bernoulli parameters
\(a_i/D\) and \(\bar a/D\). If \(u,v\) are their respective
arcsine square-root coordinates, Jensen's inequality for
the logarithm bounds this entropy below by
\(-2\log\cos(u-v)\ge c(u-v)^2\) for an absolute \(c>0\).
At least one of the two extreme \(\eta(a_i)\) has distance at
least \(T/2\) from \(\eta(\bar a)\). Thus \(J\ge c_0 T^2\).
Also
\(\left|d\log h/d\eta\right|
 \le \tfrac12\sqrt{1/a+1/(D-a)}\).
This bounds the ratio \(h(\bar a)/\prod_i h(a_i)^{1/m}\)
by \(\exp(C T/\sqrt{H_*})\).

For the derivatives, the first derivative of \(J\) in argument
\(i\) is \(s(a_i)(f_D'(a_i)-f_D'(\bar a))\), of size \(O(T)\).
For example \(s(a_i)\le\sqrt{a_i}\) and \(a_i\le m\bar a\),
so \(\sqrt{a_i}\,|\log(a_i/\bar a)|\) is bounded by
\(2\sqrt m\,|\sqrt{a_i}-\sqrt{\bar a}|\le \sqrt m\,T\).
Apply the same bound with \(D-a_i,D-\bar a\).
For \(k\ge2\) distinct arguments, the mixed derivative of
\(J\) is \(-m^{1-k} f_D^{(k)}(\bar a)\prod_{i\in E}s(a_i)\).
Each \(s(a_i)\) is bounded both by \(\sqrt{m\bar a}\) and
by \(\sqrt{m(D-\bar a)}\), so these derivatives are bounded
for \(H_*\ge1\). The mixed derivatives of
\(\log h(\bar a)+r_{D,m}(\bar a)\) are bounded by the
same calculation, since the \(k\)-th ordinary derivatives
in \(\bar a\) are \(O(\bar a^{-k}+(D-\bar a)^{-k})\).
Apply the product rule to (1) to obtain (3).

## 2. Gaussian moment systems for at most three nodes

For \(Z_m\sim\mathcal N(0,1/m)\) put
\(\phi_x(z)=e^{xz-x^2/2}\). The kernel \({\cal K}_m\) of (2)
equals \(\mathbb E\prod_i\phi_{x_i}(Z_m)\).
Use the Newton divided-difference basis on the nodes,
with repeated nodes interpreted by derivatives.

### A triple at orders one through three

For distinct nodes \(x_0<x_1<x_2\), write the coefficients
in the \(\phi_{x_i}\) basis as \(c_i\). The first moment is
\(\sum c_i\). Restricted to its zero plane, the third
moment is
\[
       3c_0c_1c_2\det
                \bigl(e^{-(x_i-x_j)^2/6}\bigr)_{0\le i,j\le2}.
 \tag{4}
\]
In fact the coefficient of \(c_i^2c_j+c_ic_j^2\) is
\(3e^{-(x_i-x_j)^2/3}\), and
\(c_i^2c_j+c_ic_j^2=-c_0c_1c_2\) on that plane.
The coefficient of \(c_0c_1c_2\) before restriction is
\(6e^{-\sum_{i<j}(x_i-x_j)^2/6}\).
Together with \(\sum c_i^3=3c_0c_1c_2\), this gives (4).

Write \(A_0,A_1,A_2\) for the Newton coefficients,
where the basis starts with \(\phi_{x_0}\), followed by
\([x_0,x_1]\phi\) and \([x_0,x_1,x_2]\phi\).
The first moment is \(A_0\).
Put \(\Delta=\prod_{i<j}(x_j-x_i)\). On \(A_0=0\),
\[
  \Delta^2 c_0 c_1 c_2=
       A_2\bigl(A_2-(x_2-x_0)A_1\bigr)
           \bigl((x_2-x_1)A_1-A_2\bigr).
 \tag{5}
\]
The determinant in (4) divided by \(\Delta^2\) extends
continuously and positively even when nodes coincide.
To see this, represent the Gaussian matrix as the Gram matrix
of \(e^{itx_i}\) under the centered Gaussian law on \(t\) of
variance \(1/3\). Passing to the Newton basis divides its
Gram determinant by \(\Delta^2\). For coincident nodes
the basis functions have different lowest Taylor degrees
\(0,1,2\) in \(t\), and each lowest coefficient is nonzero.
They remain linearly independent, and the Gram matrix is
positive definite. The same argument with variance \(1/2\)
gives a real positive definite second moment matrix in the
Newton basis. Thus on the zero plane of the first moment,
the third vanishes only on the complex lines specified
by the three real factors in (5), including at repeated nodes.
The second moment is nonzero on every nonzero point of
each such line. This proves the triple exclusion.

### A pair, with or without a separate singleton

For two nodes of difference \(d>0\), put \(x=e^{-d^2/4}\)
and \(y=x^{4/3}\). With \(u,v\) the two coefficients, the
first three moments are
\[
 u+v,\qquad
 u^2+v^2+2xuv,\qquad
 (u+v)\bigl(u^2+v^2+(3y-1)uv\bigr).
 \tag{6}
\]
The first and second have no common projective direction
by \(x<1\). The second and third have none either.
The factor \(u+v\) is excluded by the second. If the
two quadratic factors vanish, their difference forces
\(uv=0\) or \(2x=3x^{4/3}-1\). A nonzero point on the
second-moment conic has \(uv\ne0\), and for \(0<x<1\)
the other equality is impossible, by
\(3z^4-2z^3-1=(z-1)(3z^3+z^2+z+1)\) at \(z=x^{1/3}\).
At coinciding nodes, with first node zero, the Newton
phase is \(u+vZ_m\). Its first three moments
\(u,u^2+v^2/2,u^3+uv^2\) give both exclusions.

Now add a coefficient \(c\) with contribution \(\lambda c^m\)
at order \(m\), where \(\lambda>0\). At distinct nodes,
write \(s=u+v\), \(p=uv\). If \(c=0\), the first two equations
force \(u=v=0\). Otherwise set \(c=1\) by homogeneity.
The first two equations give
\(s=-\lambda\) and \(p=\lambda(\lambda+1)/(2(1-x))\).
The third, divided by \(\lambda(\lambda+1)\), would require
\[
              1-\lambda+\frac32\lambda\frac{1-y}{1-x}=0.
\]
This expression is positive because \(y<x\).
At coinciding nodes the first two equations instead
give \(u=-\lambda,\ v^2=-2\lambda(\lambda+1)\); the third
is then \(\lambda(\lambda+1)^2\ne0\).

### A quantitative margin for a pair and a singleton

The preceding exclusion has a uniform lower bound in the Newton
basis, including at a coincident pair. Fix \(R\ge0\) and
\(\Lambda\ge1\), and let the pair have distance
\(0\le\delta\le R\), with singleton weight
\(\Lambda^{-1}\le\lambda\le\Lambda\). Translate the pair to
\(0,\delta\). Its phase in the Newton basis is
\(A\phi_0+B[0,\delta]\phi\); the separate singleton has coefficient
\(c\). Define, by continuous extension at zero,
\[
 q(\delta)=\frac{2(1-e^{-\delta^2/4})}{\delta^2},\qquad
 \kappa(\delta)=\frac32
          \frac{1-e^{-\delta^2/3}}{1-e^{-\delta^2/4}},
 \qquad q(0)=\tfrac12,\quad\kappa(0)=2.
\]
The three normalized forms are
\[
 \begin{split}
 P_1&=A+\lambda c,\\
 P_2&=A^2+qB^2-q\delta AB+\lambda c^2,\\
 P_3&=A^3+\kappa qAB^2-\kappa q\delta A^2B+\lambda c^3.
 \end{split}                                                    \tag{7}
\]
Put \(q_R=q(R)\) and
\[
 \epsilon_{R,\Lambda}
 =\frac{q_R^{3/2}}
 {512(\Lambda+1)^3\Lambda(\Lambda^2+3\Lambda+4)}>0.
                                                               \tag{8}
\]
For every complex coefficient vector with
\(\max(|A|,|B|,|c|)=1\),
\[
                    \max_{1\le m\le3}|P_m|>
                                      \epsilon_{R,\Lambda}.     \tag{9}
\]

Here is a direct proof, avoiding an unevaluated compactness margin.
One has \(q_R\le q\le1/2\), \(q\delta^2\le2\), and
\(3/2\le\kappa\le2\). The first bounds follow by writing
\((1-e^{-t})/t=\int_0^1e^{-st}\,ds\). For the last bounds put
\(x=e^{-\delta^2/4}\); then
\(1\le(1-x^{4/3})/(1-x)\le4/3\), by convexity of \(x^{4/3}\).
Polynomial division gives the identity
\[
 P_3-\kappa A P_2-JP_1=\rho c^3,
 \quad
 \begin{cases}
 J=(1-\kappa)(A^2-\lambda Ac+\lambda^2c^2)-\kappa\lambda c^2,\\
 \rho=\lambda(\lambda+1)(1+(\kappa-1)\lambda).
 \end{cases}                                                   \tag{10}
\]
In particular \(\rho\ge\Lambda^{-1}\) and
\(|J|\le1+3\Lambda+\Lambda^2\) on the unit polydisc.
If all three forms had modulus at most \(\epsilon\), (10) would give
\[
 |c|^3\le C\epsilon,\qquad
 |A|\le\epsilon+\Lambda|c|,
 \qquad C=\Lambda(\Lambda^2+3\Lambda+4).
\]
Take \(t=\sqrt{q_R}/(8(\Lambda+1))\) and
\(\epsilon=t^3/C\), the value in (8).
Then \(\epsilon\le t\), \(\epsilon\le q_R/8\),
\(|c|\le t<1\), and \(|A|\le(\Lambda+1)t<1\).
From the second form and
\(q\delta|A||B|\le q|B|^2/2+q\delta^2|A|^2/2\),
\[
 |B|^2\le\frac2{q_R}
          (\epsilon+2|A|^2+\Lambda|c|^2)
 \le\frac{2\epsilon}{q_R}+\frac3{32}
 \le\frac{11}{32}<1.
\]
This contradicts the coefficient normalization and proves (9).
On the exact locus \(P_1=P_2=0\), the stronger residual formula
\(P_3=\rho c^3\) is already sufficient.

Consequently any perturbation of (7) whose coefficient absolute
errors sum to less than \(\epsilon_{R,\Lambda}\) in each row has
no common projective zero. In particular, an error less than
\(\epsilon_{R,\Lambda}/10\) in each coefficient suffices, since
there are at most ten monomials in any row.

For the angular factorial kernels this gives the proposed bridge
directly: with the pair's \(\eta\)-span at most \(R\), all three
columns interior, and the singleton-to-pair width ratio in
\([\Lambda^{-1},\Lambda]\), sufficiently large least exponent and
sufficiently large distance from the singleton to the pair imply
nonvanishing at orders \(1,2,3\). Use (2) for the pure terms and
(3), averaged over the pair's divided-difference intervals, for
the mixed terms. There is at most one derivative per pair argument,
so the error has no inverse power of the pair distance. After division
by \(\sqrt{2\pi/m}\) times the pair width, the bounded width ratio
makes these coefficient errors uniformly small. This is a quantitative
model margin for the cluster case already used below; it does not
settle higher starting windows or evaluate the absolute corner cutoff.
Unbounded width ratios still require the width-group assignment.

The [exact control](verify_sfc_pair_singleton_margin.py) derives (7)
from the Gaussian exponential kernel, verifies (10) symbolically,
and checks the coalescent forms by Gaussian polynomial integration.
The inequalities above supply the unrestricted margin proof.
An internal replay with Python 3.12.3 and SymPy 1.14.0 passed all
these controls. Run from the repository root:

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_sfc_pair_singleton_margin.py
~~~

Finally, for any number \(r\) of singleton nodes with positive
weights \(\lambda_i\), their weighted powers
\(\sum_i\lambda_i c_i^m\) have no common nonzero point
at \(r\) consecutive positive orders. Group the nonzero
\(c_i\)'s by their distinct values. Each combined weight
is positive. At the first order, multiply each combined
weight by its node value to that order. The following
powers give an invertible Vandermonde system on these
nonzero node values, so a common zero is impossible.

## 3. Boundary clusters and the final two orders

Along a sequence with \(D\to\infty\), call a column interior
if both of its entries tend to infinity. A boundary column
has one entry bounded, and after passing to a subsequence
that entry is constant. Boundary columns on the same side
have bounded \(\eta\)-differences. The two sides have
diverging separation, and so does either side from every
interior column. This follows from
\(\eta(a)\ge2\sqrt a\) and
\(\pi\sqrt D-\eta(a)\ge2\sqrt{D-a}\).

At the left side use \(t=w/D\). For fixed \(a\ge0\),
\(g_a(w/D)\to e^a a^{-a} w^a e^{-w}\), with \(0^0=1\).
For a fixed tuple of boundary columns the rescaled kernel
\(D K_{D,m}\) thus tends to its integral on \(0<w<\infty\).
For domination, on \(0<w<D\) use
\((1-w/D)^{D-a}\le e^{-w/2}\) once \(D\ge2a\);
the prefactor multiplying \(w^a\) is bounded. The other
side follows by \(t\mapsto1-t\). Thus both boundary widths
are of order \(D^{-1}\), whereas the width \(h(a)=s(a)/D\)
of an interior column is larger by a factor tending to
infinity.

Decay of mixed terms still holds when some arguments
are boundary columns. In (3), replace \(h(a_i)\) in the
product by \(\widetilde h(a_i)=\sqrt{(a_i+1)(D-a_i+1)}/D^{3/2}\)
and replace \(C T/\sqrt{H_*}\) by \(C T\).
For a mixed tuple the averages \(\bar a,D-\bar a\) both
grow, so (1) applies with its stated remainder.
Use \(h(\bar a)\le \widetilde h(\bar a)\) and the uniform
bound on \(d\log\widetilde h/d\eta\) in place of the
logarithmic width estimate. The lower bound for \(J\)
holds also at the endpoints by continuity. The derivative
argument only differentiates interior arguments, as we
never take divided differences inside a boundary cluster.
Its first-derivative bound uses the same ratios with
\(\bar a,D-\bar a\), and its mixed-derivative bounds use
the large averages. The claimed decay follows.
The width \(\widetilde h\) is comparable to \(D^{-1}\)
at every fixed boundary exponent.

If there is an interior column, the group of boundary
columns, if present, has at most two columns. Its assigned
orders in the proof below are the last one or two orders
of the window. A singleton has a nonzero pure moment.
For two opposite boundary sides, the mixed terms vanish.
The remaining quadratics and cubics are sums of the two
pure powers with positive coefficients, even when the
coefficients depend on the order. A ratio solving the
quadratic has negative real square, hence is purely
imaginary. Its cube cannot satisfy the other equation.

Two columns on the same side require the following
exclusion for the second and third moments of
\(Ue^{-w}w^a+Ve^{-w}w^b\) on \(0<w<\infty\), where
\(0\le a<b\). Absorb the positive peak constants into
\(U,V\). Normalize by the square roots of the two
pure second moments. The second and third forms are then
\[
  U^2+2xUV+V^2,\qquad
  A U^3+3 A^{2/3}B^{1/3} e\, U^2V
       +3 A^{1/3}B^{2/3} f\,UV^2+B V^3,
\]
where \(A,B>0\) are the normalized pure cubic coefficients,
\(x\) is the mixed second moment divided by the geometric
mean of the pure second moments, and \(e,f\) are the
corresponding mixed cubic coefficients divided by the
weighted geometric means of the pure cubic moments.

In fact \(0<e,f<x<1\). To prove this, let
\(h_m(t)=m^{-1}\log\Gamma(mt+1)\), and put
\[
 J_m(\theta)=(1-\theta)h_m(a)+\theta h_m(b)
                     -h_m((1-\theta)a+\theta b).
\]
The gamma integrals give
\(-\log x=2J_2(1/2)\), \(-\log e=3J_3(1/3)\) and
\(-\log f=3J_3(2/3)\). The functions \(h_m''\) increase
with \(m\), since
\[
 h_m''(t)=m\psi_1(mt+1)
       =\int_0^\infty e^{-t v}\frac{v/m}{e^{v/m}-1}\,dv .
\]
The integral follows by summing the geometric series
in the [trigamma series](https://dlmf.nist.gov/5.15).
The function \(z/(e^z-1)\) decreases for \(z>0\).
Hence \(J_3\ge J_2\) and strict concavity of \(J_2\)
gives \(3J_3(1/3),3J_3(2/3)>2J_2(1/2)>0\).
Now a common zero would force the real cubic form
to be divisible by the positive definite real quadratic:
the two nonreal projective quadratic roots are conjugate.
The quotient would be \(AU+BV\). Comparing the sum
of the two mixed cubic coefficients would give
\[
   3(A^{2/3}B^{1/3} e+A^{1/3}B^{2/3} f)=(1+2x)(A+B).
\]
The left side is at most \(3x(A+B)<(1+2x)(A+B)\), by
the preceding strict bounds and weighted arithmetic-geometric
mean. Thus the boundary pair has no common zero.

## 4. Uniformity over all exponent gaps

Suppose there were counterexamples with one column having
both entries tending to infinity. Work over \(\mathbb C\),
which suffices by the rational coefficient equations and the
Nullstellensatz on coefficient charts. Order the three
columns by \(a_j\). Pass to a subsequence on which the
two adjacent differences in \(\eta(a_j)\) each have a
finite limit or tend to infinity, and each column is interior
or has a fixed boundary entry as in Section 3. Group the
nodes across the bounded differences. For an interior cluster
use the Newton divided differences of \(g_a\) in \(\eta(a)\)
as basis elements. For a boundary cluster use the \(g_a\) directly.
This coefficient change is invertible and the same at every order.
For cluster \(g\), write \(h_g\)
for the width \(h(a)\) at its first node if interior and
\(D^{-1}\) if on a boundary. Write \(A_g\) for its coefficient
vector, with \(N_g=\|A_g\|_\infty\).

Pure contributions from an interior cluster \(g\), divided
by \(\sqrt{2\pi/m}\,h_g\), converge as coefficient forms to
the confluent Gaussian moments of Section 2, by (2).
For a boundary cluster use the gamma moments of Section 3
with the row factor \(\sqrt{2\pi/m}\) divided out.
Mixed contributions between clusters, after the same
division by \(\sqrt{2\pi/m}\) and by the product of the
cluster widths to the proportion of factors used from each,
have coefficients tending to zero. Indeed with more
than one cluster all divided differences have order at
most one. Write them as averages of first derivatives
and apply (3) or its boundary version. Points in the
intermediate divided-difference intervals still have bounded
distance from their respective
cluster reference nodes in \(\eta\). Their widths are
comparable to those reference widths by the logarithmic
derivative bounds, and \(T\) tends to infinity for every
mixed term. Consequently the error from mixed terms,
after inserting the coefficients and dividing by
\(\max_g h_g N_g^m\), tends to zero. Here the product of
the weights \(h_gN_g^m\) raised to their factor proportions
is at most their maximum.

A final grouping handles different widths. After a
subsequence, divide the geometric clusters into groups of
comparable widths, meaning that within each group every
ratio of widths tends to a finite positive number.
Order these groups by decreasing width, so that between
two different groups the ratio of the earlier to later
width tends to infinity. For a width group \(b\), choose
one reference width \(h_b\), set \(N_b=\max_{g\in b}N_g\),
and use the order-\(m\) weight
\[
                         W_b(m)=h_b N_b^m .
\]
The maximum of these weights is comparable to
\(\max_g h_gN_g^m\). A group with zero coefficient vector
has weight zero.

Take \(R\to\infty\) sufficiently slowly along the sequence
and retain in row \(m\) every group with weight at least
\(1/R\) of the maximum. These active sets satisfy
\(\max S_m\le\min S_{m+1}\) in the decreasing-width order.
For if a later group \(b\) is active at \(m\) and an earlier
group \(a\) is active at \(m'>m\), then
\[
 \frac{W_b(m')}{W_a(m')}
   \ \ge\ (h_a/h_b)^{m'/m-1}R^{-m'/m}.
\]
Since every such \(h_a/h_b\) diverges, \(R\) can be chosen
slowly enough that the right side exceeds \(R\) for all
possible pairs and rows. This contradicts activity of \(a\).

By the [ordered assignment lemma](SFC_SEPARATED_DEGREE_CLUSTER_WINDOWS.md#1-the-ordered-assignment-of-moment-rows),
some width group is uniquely active throughout its
assigned consecutive window, whose length is the
number of original columns in that group. Pass to a
subsequence on which this group \(b\) is fixed. Divide
each assigned moment by \(\sqrt{2\pi/m}\,h_b N_b^m\),
and pass to a limit of the normalized group coefficients
\(A_g/N_b\). The normalized vectors of this group have
maximum norm one. All other width groups and all mixed
terms disappear. The remaining equations are sums of the
local moments in this group with the positive weights
\(\lim h_g/h_b\). If the group contains boundary columns, it
contains only those columns and is last in decreasing width.
Its last one or two orders are excluded by Section 3.

In an interior group made of singleton clusters, the weighted
Vandermonde system is excluded in Section 2. A pair on
its own has the assigned window \(1,2\) or \(2,3\).
A group containing a pair and a separate singleton uses
the whole window \(1,2,3\). A triple in one cluster also
uses that window. Section 2 gives nonvanishing in every
case. This contradicts the common zeros and proves the
uniform cutoff.

## Remaining axis profiles

For a fixed list of the bounded exponents \(u_j\) and fixed
assignments of the columns to the two sides, there is one
remaining degree parameter \(D\). If failures occur for
arbitrarily large \(D\) on this profile, the following rational
system has a nonzero common zero, at \(m=1,2,3\):
\[
 \sum_{\text{side }S}
   \ \sum_{\substack{|\alpha|=m\\ \operatorname{supp}\alpha\subseteq S}}
     \binom m\alpha c^\alpha
     \frac{(\sum_j \alpha_j u_j)!}{m^{1+\sum_j\alpha_j u_j}}=0 .
\]
This is the boundary limit of Section 3 with the positive peak
constants absorbed into the coefficients. Both sides use the
same width \(D^{-1}\), so normalize all coefficients together.
The mixed-side terms tend to zero by the same kernel bound.
The two possible patterns are three columns on one side or
a pair on one side and a singleton opposite it. A nonvanishing
certificate for the displayed system rules out unbounded
failed degrees on that profile. The [fractional degree sieve](SFC_HOMOGENEOUS_TRINOMIAL_AXIS_DEGREE_SIEVE.md)
gives an effective degree bound without requiring the boundary
limit to be nonvanishing.
