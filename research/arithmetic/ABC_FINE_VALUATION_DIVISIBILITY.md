# A count saving for high divisibility of the primitive valuations

The coefficient-uniform generalized Fermat estimate gives a small-budget
saving on a sector whose defining restriction is invisible to the weak
scaled-degree measures. The large common exponent may have scaled width
tending to zero, so this is compatible with formal diffuse measures.

## Exact sector bound

Let \(K\ge2\) be an integer and \(f\ge0\). Count positive primitive
\(a+b=c\le X\) admitting factorizations
\[
 a=s_1x_1^{k_1},\qquad b=s_2x_2^{k_2},\qquad
 c=s_3x_3^{k_3},\qquad k_1,k_2,k_3\ge K,
 \qquad \prod_{i=1}^3\operatorname{rad}(s_i)\le X^f .
\]
All factors here are positive integers. The orders can depend on the
triple and need not be bounded independently of \(X\). Then
\[
 \#\{(a,b,c)\}\ll_{K,\varepsilon}
        X^{f+3K^{-3/2}+\varepsilon}
\]
uniformly in \(f\). The implied constant can be taken uniform in \(f\)
by replacing \(X^f\) by its minimum with \(X^3\) when counting the
coefficient choices.

To see this, choose the joint coefficient radical and assign its
primes to the coefficient letters. Primitivity makes their supports
disjoint, so there are \(3^{\omega(R)}=X^{o(1)}\) assignments uniformly
for \(R\le X^3\). The fixed-support estimate gives \(X^{o(1)}\) triples
of exact coefficients per allocation. There are therefore at most
\(X^{f+o(1)}\) coefficient choices. Fix one such combination. The
variables \(x_i\) are pairwise coprime.

Put \(H=\lceil K^{3/2}\rceil\). If all \(k_i\le H\), there are only
finitely many choices of the exponents, depending on \(K\). For each
choice use Browning–Verzobio, Theorem 1.2. Its two ordinary input
boxes have heights at most \(X^{1/K}\) and the distinguished exponent is
at least \(K\). The auxiliary \(W\) of that theorem is at most
\(X^{K^{-3/2}}\). The three terms in the estimate are bounded,
respectively, by
\[
 X^{2K^{-3/2}},\qquad X^{3K^{-3/2}},\qquad
 X^{K^{-3/2}+K^{-2}} .
\]
The theorem is uniform in the fixed coefficients. If a root box has
trivial size, one can enlarge it to \(X^{1/K}\) since the theorem is
monotone and applies when box sides are at least two.
[Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)

If some \(k_i>H\), enumerate its actual term
\(s_i x_i^{k_i}\) in at most \(O(X^{1/H}\log X)\) ways,
including the term \(s_i\) when \(x_i=1\). For \(x_i>1\) the
order is at most \(\log_2 X\).
If two or three orders exceed \(H\), enumerate two such terms and
the third is determined. The count is \(X^{2/H+o(1)}\).

If just one order exceeds \(H\), after fixing its term the other two
roots satisfy
\[
 \alpha y^{k}+\beta z^{j}=\gamma,\qquad K\le k,j\le H,
 \qquad \alpha\beta\gamma\ne0 .
\]
This plane curve is absolutely irreducible. Over
\(\mathbf C(y)\), the function
\((\gamma-\alpha y^k)/\beta\) has simple zeros, so adjoining its
\(j\)-th root has degree \(j\). Equivalently the binomial is
Eisenstein at one of those simple zeros. The plane equation follows
by Gauss's lemma. The coefficient-uniform integral plane-curve estimate
in the \(X^{1/K}\)-square gives at most
\(X^{1/(K\max(k,j))+\varepsilon}\le X^{K^{-2}+\varepsilon}\)
points. Thus this case costs at most \(X^{1/H+K^{-2}+\varepsilon}\).
Every branch fits \(3K^{-3/2}\).
[Browning–Verzobio, Lemma 2.2](https://arxiv.org/html/2608.24512v1#S2.SS1)

## Consequence for a genuinely fine valuation grid

At fixed \(0<\lambda<1\), take
\(K=\lceil(72/\lambda)^{2/3}\rceil\). On any family where the primes
whose valuations are *not* multiples of their own letter's order
\(k_i\ge K\) have joint radical at most \(X^{\lambda/4}\), place
all these prime powers in the \(s_i\). The displayed estimate gives
\[
 \#\{(a,b,c)\}\ll_{\lambda,\varepsilon}
 X^{(7/24)\lambda+\varepsilon}.
\]
It does not require the radical cap on the other primes. In particular
one may impose this divisibility on at least three quarters of the
normalized radical budget of an abc family. The exponent \(7/24\)
then leaves a fixed gap below \(1/3\).

The condition can select a grid with scaled mesh
\(\lambda K=O(\lambda^{1/3})\to0\). For example, the measures
obtained by assigning mass \((1/3)/J_\lambda\) to each of the
\(J_\lambda\) degree positions in \([2,4]\) that belong to
\(\lambda K\mathbf N\) converge weakly to
\((1/3)\mathrm{Unif}[2,4]\). Thus diffuse limiting measures
do not rule out the grid condition at the profile level. The estimate
excludes arithmetic realization at the endpoint for this sector.

## Effect of one-target minimal repair

The least multiplier \(Z_M(c)\) is coprime to \(ab\). At a prime of a
primitive summand the valuation in \(Z_M(c)a\) or \(Z_M(c)b\)
is therefore unchanged. In particular the repair cannot make a
primitive summand into a high exact power unless its original
valuations already satisfy the corresponding divisibility.

When all original target valuations are at most \(M\), the repaired
target is \(\operatorname{rad}(c)^M\). If \(k\mid M\) and either
repaired summand is also an exact \(k\)-th power, then at every prime
\(p\mid c\) its valuation is \(M-v_p(c)\). Consequently \(k\mid
v_p(c)\), so the target was already an exact \(k\)-th power before
repair. One can instead put \(Z_M(c)\) in the nonpower coefficients,
but its varying support then enters the coefficient radical cost in
the sector bound. These facts isolate the genuine valuation-grid
input needed for this use of high-power arithmetic.

## Unequal grid spacings

A precise anisotropic statement follows by keeping the input boxes of
Browning–Verzobio separate. Suppose the orders \(p,q,r\ge2\) are fixed
before the height limit, and the three corresponding complete
height caps in base \(X\) are \(h_1,h_2,h_3\). Put
\(u=h_1/p\), \(v=h_2/q\), and \(w=\sqrt{uv/r}\).
The same coefficient enumeration gives exponent at most \(f+E\), where
\[
 E=\max\left(2w,\quad
             w+\frac{2\max(u,v)}{\sqrt{\max\{p,q,36\}}},\quad
             w+\frac{\max(u,v)}{r}\right).
\]
One may minimize over the three choices of distinguished letter.
When an indicated height cap is zero, take arbitrarily small positive
enlargements and pass to the limiting inequality. In particular with
all complete heights at most one, two orders of size at least
\(\lambda^{-\alpha}\) and a third at least
\(\lambda^{-\beta}\) incur \(E=o(\lambda)\) provided
\[
 \alpha>2/3,\qquad \beta>0,\qquad \alpha+\beta/2>1 .
\]
Indeed \(w\le\lambda^{\alpha+\beta/2}\), the second increment is at
most \(2\lambda^{3\alpha/2}\), and the third at most
\(\lambda^{\alpha+\beta}\).
This permits a substantially finer valuation grid on the third
letter while retaining the same coefficient-budget saving.


## Relation to the external count estimates

The input from Browning and Verzobio is the coefficient-uniform
fixed-exponent bound. The radical enumeration and the split at
\(K^{3/2}\) give the estimate for triples with variable grid orders.
Bernert, Browning, Lichtman and Teräväinen give an unrestricted
\(X^{3/5+\varepsilon}\) estimate for \(0<\lambda<1\).
Zeng's complete-modulus-locking preprint reports an unrestricted
exponent \(11680007/20000012\) independent of \(\lambda\). These
constant exponents do not address a saving proportional to
\(\lambda\) as \(\lambda\downarrow0\). Zeng's separate Pell preprint
lists a \(\lambda/3-C\lambda^3\) exponent on its three fixed linear
forms \(n,n+1,2n+1\). Its input uses unit differences.
[Generalized Fermat surface estimate](https://arxiv.org/html/2608.24512v1#S1)
· [Bernert–Browning–Lichtman–Teräväinen](https://arxiv.org/html/2410.12234v2#S1)
· [Complete modulus locking](https://hexagonmath.org/pdf/2611.00621v1)
· [Small Pell discriminants](https://hexagonmath.org/pdf/2611.00611v1)

The fine-grid saving concerns unrestricted signed coefficients with
small joint radical and variable exact power orders. Arbitrary
nonatomic limiting spectra do not give those divisibility
conditions. In particular a model with weights spread evenly across
all integer degree positions in the compact band has asymptotically
negligible mass on every grid of mesh \(k\ge K\), while its weak
limit is the same as that of the grid-supported example above.

Two exact high orders can also be handled with an arbitrary
third-letter spectrum. A rectangular finite-place determinant
uses a low-cost full factor of that letter, or the generalized
Fermat estimate uses a single macroscopic prime on it. In the
symmetric case fixed large multiples of \(\lambda^{-1/2}\) suffice.
[Two valuation grids](ABC_TWO_VALUATION_GRIDS.md)
