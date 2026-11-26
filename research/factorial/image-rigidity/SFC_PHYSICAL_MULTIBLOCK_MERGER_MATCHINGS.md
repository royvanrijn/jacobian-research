# Simultaneous mirrored merger limits and permutation matchings

When finitely many mirrored factorial blocks approach complete
correlation at comparable rates, the scaled angular polynomial
is a finite free additive convolution of the individual limits.
It is also the expected weighted matching polynomial of a union
of independent colored permutation graphs. The limiting roots
are simple and real in every degree. The comparison of the
two adjacent limiting root sets for arbitrary blocks remains
open.

## The simultaneous limit

For each block \(i\), let \(a_i\ge0\) be its smaller exponent
and \(d_i\ge1\) its positive gap, both integers. Set
\(b_i=2a_i+d_i\). Start with independent standard complex
Gaussians \(Z_i,W_i\), with variance one in each complex
coordinate, and let
\[
 X_i=|Z_i|^2,\quad
 Y_i=\left|\sqrt{1-\kappa_i\delta}\,Z_i+
             \sqrt{\kappa_i\delta}\,W_i\right|^2,
 \qquad \kappa_i\ge0 .
\]
At least one \(\kappa_i\) is positive. Take \(\delta>0\)
small enough for all square roots to be real. The two monomials
are
\[
 U_\delta=\prod_i (X_iY_i)^{a_i}X_i^{d_i},\qquad
 V_\delta=\prod_i (X_iY_i)^{a_i}Y_i^{d_i}.
\]
For \(m\ge1\) put \(P_m(z;\delta)=\mathbb E(U_\delta+zV_\delta)^m\).
Use \(y^{\underline{k}}\) for a falling factorial and set
\(h_i=\kappa_i d_i^2\).

**Theorem.** Locally uniformly in complex \(u\),
\[
 \frac{P_m(-e^{i\sqrt{\delta}u};\delta)}
 {(-i)^m\delta^{m/2}\prod_i\Gamma(mb_i+1)}
 \longrightarrow
 p_m(u)=\sum_{k=0}^{\lfloor m/2\rfloor}
  \frac{(-1)^k m^{\underline{2k}}}{k!}
        \mathbb E\left(\sum_i\frac{h_i}{G_i}\right)^k u^{m-2k},
                                                               \tag{1}
\]
where the \(G_i\) are independent, unit-scale Gamma variables
of shapes \(mb_i+1\). The polynomial \(p_m\) is monic, with
\(m\) distinct real roots.

For the limit, \((Y_i-X_i)/\sqrt\delta\) converges in every
finite \(L^p\) to \(\sqrt{\kappa_i X_i}\,H_i\), where the
\(H_i\) are independent real Gaussians of variance two and
independent of the \(X_i\). This follows by expanding the
squared modulus and using rotational symmetry of \(W_i\).
Let \(M=\prod_i X_i^{b_i}\). Factoring each difference
\(Y_i^{d_i}-X_i^{d_i}\) and telescoping the products gives
\[
 \frac{U_\delta-e^{i\sqrt\delta u}V_\delta}{-i\sqrt\delta}
   \longrightarrow
 M\left(u-i\sum_i\frac{\sqrt{\kappa_i}\,d_i H_i}{\sqrt{X_i}}\right)
\]
in every finite \(L^p\), uniformly for \(u\) in compact
sets. The apparent inverse powers on the right are
multiplied by \(M\), and \(b_i\ge1\). The finite polynomial
products and the exponential Taylor remainder justify
passage to the \(m\)-th moment by Hölder's inequality.
Conditioned on the \(X_i\), the real sum in this expression
has variance \(2\sum_i h_i/X_i\). Integrating its even
moments and tilting the independent \(X_i\) by \(M^m\)
proves (1).

To see simplicity, write \(\partial=d/du\). For each block
the constant-coefficient operator
\[
  \sum_{k=0}^{\lfloor m/2\rfloor}
          \frac{(-h_i)^k}{k!\,(mb_i)^{\underline{k}}}\partial^{2k}
\]
takes \(u^m\) to that block's polynomial. The operator product
takes \(u^m\) to \(p_m\) by the multinomial expansion in (1).
The product is precisely the degree-\(m\) symmetric additive
convolution. Every active single-block polynomial is a scaled
Gegenbauer polynomial of shape \(1+(b_i-1)m\), with simple
real roots. An inactive block is the identity \(u^m\).
[Single-block merger limit](SFC_PHYSICAL_MIRROR_MERGER_EDGE.md#the-angular-scaling-limit-for-an-ordinary-block)
[Symmetric additive convolution and Walsh's real-rootedness theorem](https://doi.org/10.1007/s00440-021-01105-w)

For completeness, the linear convolution with a fixed real-rooted
degree-\(m\) factor preserves real-rootedness also in degrees
below \(m\), by adjoining real roots that escape to infinity
and taking a coefficient limit. For every positive shift
smaller than the minimum gap of a simple real-rooted input,
that input and its shifted copy interlace. All their real
linear combinations are real-rooted. Convolution is linear,
commutes with translation, and preserves real-rootedness, so
the two output translates weakly interlace by the same
interlacing characterization. If the output had a repeated
root, a shift avoiding its finitely many root differences
would violate interlacing between those two equal roots.
Starting with an active block therefore gives simplicity.

For each fixed \(m\), the simple scaled roots account for all
\(m\) zeros at sufficiently small positive \(\delta\), by
Rouché's theorem on small disjoint disks in \(u\).
The scaled angular roots are real by the exact circle-location
theorem. Strict interlacing of \(p_n,p_{n+1}\), when available,
gives strict angular interlacing in a neighborhood of the
merger for that fixed pair of orders.
[Physical circle location](SFC_PHYSICAL_MIRRORED_ROOT_CIRCLE.md)

## A matching-polynomial representation

Give block \(i\) an independent uniform permutation on \(mb_i\)
slots, of which \(m\) are distinguished and identified with
common labeled vertices. Whenever the permutation maps
one distinguished vertex to a different one, put an undirected
edge copy between them, of color \(i\) and weight \(h_i\).
Opposite arrows give parallel copies. Discard loops and edges
involving nondistinguished slots. Then \(p_m\) is the expected
matching polynomial of the resulting weighted multigraph,
with the convention that a matching of size \(k\) contributes
\((-1)^k u^{m-2k}\) times the product of its edge weights.

In a matching with \(k_i\) edges of color \(i\), the prescribed
arrows have distinct starts and ends. The probability that
all occur is \(\prod_i[(mb_i)^{\underline{k_i}}]^{-1}\).
For fixed counts \(k_i\) with sum \(k\), the number of labeled
oriented colored matchings is
\(m^{\underline{2k}}/\prod_i k_i!\). This proves the identity
using \(\mathbb E G_i^{-k_i}=1/(mb_i)^{\underline{k_i}}\).

For the consecutive comparison, uniform permutations can be
coupled by inserting \(b_i\) slots of color \(i\), one
distinguished. Deleting those slots and stitching their
cycles recovers the lower uniform permutation. Stitching
changes edges between surviving vertices, so this coupling
does not identify the lower graph with a vertex deletion
of the upper graph. The usual matching-polynomial interlacing
under vertex deletion does not directly compare these two
averages.

Simplicity and adjacent interlacing of separate inputs also
do not suffice for a convolution argument. For example, take
\(a=99/100\), \(c=1/100\), and the strictly interlacing pair
\(f_3=u(u^2-a^2)\), \(f_4=(u^2-1)(u^2-c^2)\). Adding two copies
at each respective degree gives \(u(u^2-2a^2)\) and
\[
 u^4-2(1+c^2)u^2+2c^2+(1+c^2)^2/6 .
\]
The quartic is positive at \(u=\sqrt{2}\,a\), by direct
rational substitution. This precludes strict interlacing.
The example is an abstract convolution pair. Formula (1)
imposes additional inverse-Gamma constraints in the ordinary
merger limit.
