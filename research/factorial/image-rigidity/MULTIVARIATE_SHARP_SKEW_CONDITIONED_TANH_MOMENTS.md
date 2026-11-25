# Sharp capped tanh moments at fixed skew

Let \(0\le a_i\le A\) with \(\sum_i a_i\le C\) and signs
\(\epsilon_i\in\{-1,1\}\). Assume \(A,C>0\) and prescribe the
feasible skew \(P=\tfrac12\sum_i\epsilon_i a_i^2\).
Exchange the signs first if \(P<0\).

For \(K\ge0\) write \(m(K)=\lfloor K/A^2\rfloor\) and
\(b(K)=K-m(K)A^2\). Put
\[
 G(K)=m(K)A+\sqrt{b(K)},\qquad
 J_\varphi(K)=m(K)\varphi(A^2)+\varphi(b(K)).
\]
The unique solution \(k_*\ge0\) of
\(G(k_*)+G(k_*+2P)=C\) gives the sharp bounds
\[
 \sum_{\epsilon_i=-1}\varphi(a_i^2)\le J_\varphi(k_*),\qquad
 \sum_{\epsilon_i=+1}\varphi(a_i^2)\le J_\varphi(k_*+2P)
\]
for every nondecreasing convex \(\varphi\) on \([0,A^2]\) with
\(\varphi(0)=0\). The two bounds and their sum are attained by the
same orientation lists for every such kernel.

For a fixed sum \(K\) of squared entries, the list of \(m(K)\)
entries \(A^2\), a remainder \(b(K)\), and any needed zeros
majorizes every admissible list. Indeed the largest \(h\) entries
sum to at most \(\min(hA^2,K)\). Concavity of the square root
makes the packed list minimize the sum of the original entries.
Convexity of \(\varphi\) makes it maximize the kernel sum.
If the negative quadratic mass is \(k\), the positive mass is
\(k+2P\). Thus the admissible \(k\) are exactly those with
\(G(k)+G(k+2P)\le C\). The left side is continuous, strictly
increasing, and unbounded, and starts at or below \(C\) by
feasibility. This proves the equation for \(k_*\), including its
zero solution. The displayed bounds now follow from the
monotonicity of \(J_\varphi\), and the two packed lists at \(k_*\)
attain them.

With \(\varphi(w)=w^{r/2}\) this gives a sharp upper bound on the
unsigned \(r\)-th moment for every real \(r\ge2\). For example,
at \(C=1\), \(A=2/5\), \(P=0\), the sharp fourth-moment bound is
\(257/5000\), attained with entries \(2/5,1/10\) in each
orientation. The [earlier conditioned bounds](MULTIVARIATE_SKEW_CONDITIONED_TANH_MOMENTS.md)
give \(33/625\) there.

The absolute values of the [even and odd
logarithmic kernels](MULTIVARIATE_MATCHED_TANH_KERNEL_TRANSPORT.md#two-full-kernels-for-the-factorial-coefficients)
at real grid positions \(|z|\le1\) also fit these bounds when the
scaled argument stays below one. Their series after the leading
term is removed have nonnegative coefficients in powers of
\(w=a_i^2\) of degree greater than one. The deficit scale belongs
inside the kernel. Thus the same packed lists give simultaneous
envelopes for the complete absolute kernels in each orientation.

These are sharp constraints on capped lists at either order.
The integer-row transport further restricts a pair of lists.
Angular separation also requires control of the exponentiated
and alternating sums at the relevant root gaps.


## Sharp packing with the actual outer loss

Assume \(0<A\le1/\sqrt3\) and replace the mass constraint by
\(\sum_i\operatorname{artanh}a_i\le U\), with \(U\ge0\).
Put \(s_*=\operatorname{artanh}A\) and
\[
  G_{\rm loss}(K)=m(K)s_*+\operatorname{artanh}\sqrt{b(K)}.
\]
The same two sharp kernel bounds hold at prescribed skew, using
the unique \(k_*\ge0\) satisfying
\[
               G_{\rm loss}(k_*)+G_{\rm loss}(k_*+2P)=U
\]
after exchanging the signs when necessary. The equality lists have
total loss exactly \(U\).

On \(0<w\le A^2\) the loss function satisfies
\[
  \frac{d^2}{dw^2}\operatorname{artanh}\sqrt w
             =\frac{3w-1}{4w^{3/2}(1-w)^2}\le0 .
\]
It extends continuously at zero, hence is concave on the closed
interval. For prescribed quadratic mass the packed squared list
therefore minimizes the loss and maximizes every convex kernel
sum. The displayed sum of the two \(G_{\rm loss}\) functions is
the exact minimum loss at those two quadratic masses. It is
continuous, strictly increasing and unbounded, and starts at or
below \(U\) by feasibility. A zero remainder contributes zero
loss. If \(U=0\), all factors vanish and feasibility forces \(P=0\).

Without a prescribed skew, packing
\(s_i=\operatorname{artanh}a_i\le s_*\) gives the attained
bound for every real \(r\ge2\)
\[
 \sum_i a_i^r\le mA^r+\tanh(U-ms_*)^r,\qquad m=\lfloor U/s_*\rfloor .
\]
Writing \(T=\tanh s\), the second derivative of \(\tanh(s)^r\) is
\[
 rT^{r-2}(1-T^2)\big((r-1)-(r+1)T^2\big)\ge0
\]
on this interval, also at zero by continuity. Thus filling the
loss budget with \(m\) entries \(s_*\) and one remainder maximizes
the sum for all these powers. The case \(r=2\) bounds \(|P|\) by
half the displayed power sum, attained when all signs agree in
the capped loss model.

For ordinary Gauss factors at order \(n\), a reciprocal cap
\(c>0\) gives \(L=1+(n-1)c\), \(A=(L-1)/(L+1)\), and
\(s_*=\tfrac12\log L\). At a rational outer mean \(u\in(0,1)\),
\(U=-\log u\), the integer \(m\) is the largest one satisfying
\(u^2L^m\le1\), and
\[
                 \tanh(U-ms_*)=\frac{1-u^2L^m}{1+u^2L^m}.
\]
The sharp bound on \(P/(UA^2/(2s_*))\) can therefore be
enclosed by rational products and outward logarithms. On a
closed mean interval with the same \(m\) at both endpoints,
the convex numerator \(mA^2+\tanh(U-ms_*)^2\) divided by \(U\)
is bounded by the larger endpoint ratio. To see this, bound
each interior numerator by the affine combination of endpoint
numerators, then bound each endpoint numerator by the same
ratio times its endpoint \(U\). If the endpoint integers
differ, the ordinary normalized bound \(1\) remains valid.
At \(u^2L^m=1\) the new remainder vanishes and the preceding
segment gives the same limit. At \(u=1\) the loss and skew
vanish. If the reciprocal cap is zero, all factors vanish
as well.
