# Long zero windows under strictly positive factorial weights

Let \(m,L\) be positive integers. If \(m>L\), there
are a polynomial \(R\in\mathbb Q[x]\) with strictly
positive coefficients in every degree from zero to its
degree, and an affine binomial \(f=x-A-iB\) with positive
integers \(A,B\), such that
\[
           \int_0^\infty f(x)^j R(x)e^{-x}\,dx=0
                       \qquad(m\le j\le m+L-1).
                                                               \tag{1}
\]
The density is normalized and strictly positive on
\((0,\infty)\), has rational moments, and has a strictly
totally positive Stieltjes moment kernel. In particular
all minors with increasing row and column exponents in
the kernel \(\mathbb E[X^{a+b}]\) are positive.

The condition \(m>L\) is sharp for an affine phase with
nonzero constant imaginary part and arbitrary positive
reweighting. If \(L\ge m\), no such measure can satisfy
(1). The reweighting is essential to the examples.
Ordinary factorial moments in the strong factorial
problem have \(R=1\).

## The angular obstruction to convex separation

Choose \(R_0>0\) so large that
\[
                  m\bigl(2\arctan R_0\bigr)>L\pi .
                                                               \tag{2}
\]
On \([-R_0,R_0]\) define a real vector of length \(2L\)
by taking real and imaginary parts of
\((t-i)^m,\ldots,(t-i)^{m+L-1}\). Any nonzero linear
functional on this vector has the form
\[
              \operatorname{Re}\bigl[(t-i)^m Q(t-i)\bigr],
              \qquad 0\ne Q\in\mathbb C[z],\quad
                       \deg Q\le L-1 .
                                                               \tag{3}
\]
We show it cannot be nonnegative on the whole interval.

Factor \(Q(t-i)\) as a polynomial in \(t\). Write \(r\)
for the number of its distinct real roots interior to
the interval and \(q\) for its number of nonreal roots,
counting multiplicity. Remove the \(r\) real roots and
consider the \(r+1\) remaining open intervals. If (3)
were nonnegative, the product \((t-i)^m Q(t-i)\) would
lie in the closed right half-plane on each interval and
would never be zero there. Its continuous argument
would then stay in a closed interval of length \(\pi\),
so its net change on each interval is at most \(\pi\).

On each such interval a real linear factor of
\(Q(t-i)\) has constant argument. For a nonreal root
\(\zeta\), the argument derivative of \(t-\zeta\) is
\[
                 \frac{\operatorname{Im}\zeta}
                  {(t-\operatorname{Re}\zeta)^2+
                                (\operatorname{Im}\zeta)^2},
\]
whose integral of the absolute value over the real
line is \(\pi\). The sum of its absolute angular
changes over the intervals is thus at most \(\pi\).
The argument of \(t-i\) increases across the original
interval by \(2\arctan R_0\) and has no jumps at the
removed roots. Adding the net changes of the product
over the open intervals therefore gives
\[
                   m(2\arctan R_0)\le (r+1+q)\pi
                            \le (\deg Q+1)\pi\le L\pi .
\]
This contradicts (2). Taking limits from the interior
handles roots at the interval endpoints.

By convex separation the origin is an interior point
of the convex hull of the feature vectors. In fact
the hull is compact and there is no nonzero linear
functional which is nonnegative on that hull. Rational
nodes suffice. Choose a finite set of nearby points
in the hull of rational-node features surrounding a
small simplex about zero in the original hull. The
finite collection of rational nodes used to form those
points still has zero in the interior of its hull.

## An exact polynomial weight

Choose an integer \(a>R_0\), and translate each rational
node by \(a\) to a positive rational node \(x_i\). Let
\(B\) be a large positive integer multiple of all node
denominators. Replace node \(x_i\) by a Gamma law of
shape \(Bx_i\) and rate \(B\). Its \(k\)-th raw moment is
\((Bx_i)_k/B^k\), which is rational and tends to \(x_i^k\)
as \(B\) increases. The feature columns for real and
imaginary parts of \((x-a-i)^j\), \(m\le j\le m+L-1\),
therefore converge to the original columns, and are
themselves rational. For sufficiently large \(B\)
their convex hull still contains zero in its interior.

Include all missing Gamma shapes from \(1\) through the
largest shape \(K\). This only enlarges the convex hull.
All these Gamma laws can be given strictly positive
rational weights with feature expectation zero. To see
strict positivity, put a common small positive weight
on every column and renormalize the remaining target.
The target tends to zero with the initial weight and
is still in the hull of the original columns. This
gives a strictly positive real solution. The defining
linear equations have rational coefficients, so
rational solutions are dense in their real affine
solution space. Small rational perturbations preserve
strict positivity and keep the weights summing to one.

For a variable \(V\) with the resulting mixture, put
\(X=BV\). Its density is \(R(x)e^{-x}\), with
\[
                  R(x)=\sum_{k=1}^K
                           \frac{\omega_k}{(k-1)!}x^{k-1},
                  \qquad \omega_k>0,\quad
                         \omega_k\in\mathbb Q .
\]
It is normalized because each term is its weight
times a unit-rate Gamma density. Taking \(A=Ba\) and
\(f(X)=B(V-a-i)\) proves (1). This construction uses
only finite polynomial moment equations and
positive rational weights.

## Strict total positivity and the sharp bound

For any positive density on \((0,\infty)\) with finite
moments, multilinearity of determinants and
symmetrization give
\[
   \det\bigl(\mathbb E[X^{a_i+b_j}]\bigr)_{i,j=1}^k
     =\frac1{k!}\int
          \det(x_j^{a_i})_{i,j=1}^k
          \det(x_j^{b_i})_{i,j=1}^k
                                      \prod_{j=1}^k d\mu(x_j).
\]
Each monomial determinant is positive for strictly
increasing nonnegative integer exponents and
\(0<x_1<\cdots<x_k\). Indeed it is nonzero by
Descartes' rule of signs, and its sign can be read
at \(x_j=y^{j-1}\), \(y>1\), where it is an
ordinary Vandermonde. Thus the product in the
integrand is nonnegative everywhere and positive
on an ordered open set of positive measure. This
proves the strict total positivity claim.

Finally let \(f=x-A-iB\) with \(B\ne0\) and
let \(\mu\) be any nonzero positive measure on the
real line with the requisite finite moments. If its
moments of \(f\) vanish at every order from \(m\)
through \(2m-1\), then
\[
                  \int f^m(f+2iB)^{m-1}\,d\mu=0
\]
by expansion. On the real line \(f+2iB=\overline f\),
so the imaginary part of this integral is
\(-B\int |f|^{2m-2}\,d\mu\ne0\). This contradiction
proves sharpness.

This obstruction strengthens the
[positive joint moment-law examples](MULTIVARIATE_JOINT_POSITIVE_MOMENT_SHAPE_OBSTRUCTION.md)
in a different direction: in one variable, an
arbitrarily long block can vanish despite strict
total positivity of the entire mixed moment kernel.
The distinction from the
[ordinary strong factorial conjecture](https://arxiv.org/abs/1304.3956)
is the polynomial weight \(R\) that remains fixed
throughout the block.
