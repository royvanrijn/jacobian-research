# A primitive fibre with minimal triangle-ratio entropy

A carry-free digit graph has a large actual primitive fixed-target
fibre, the minimum exponent for the number of primitive difference
triangles, and no comparable subset described by bounded-rank
integer features of small maximal minors. Thus triangle-ratio
cardinality and fixed-order additive energies cannot by themselves
supply the finite-rank hypothesis of the fixed-fibre count.

## Exact addition, height and primitivity

Fix \(0<L<1\) and put \(t=6/L-1\). Let \(d\to\infty\) and choose
integers \(H\to\infty\) with \(\log d=o(\log H)\). Set
\[
 N=H^d,\qquad m=\lceil dH^t\rceil,\qquad Q=4mH,\qquad M=Q^d,
 \qquad F(\boldsymbol x)=\sum_{j=0}^{d-1}x_jQ^j.
\]
Use \(S=F(\{0,\ldots,H-1\}^d)\) and \(s_*=\max S\).
Choose integer shifts \(a_0,b_0\in[M,2M)\), let
\[
 A=a_0+S,\qquad B=b_0+S,\qquad C=a_0+b_0+(S+S),
 \qquad X=a_0+b_0+2s_*.
\]
Retain those edges \((a,b,a+b)\) with \(a\in A,\ b\in B\)
and \((a,b)=1\). The targets are in \((X/2,X]\), and all
three letters are comparable in height to \(X\).
Here \(X\asymp M\) and
\(\log X=(t+1+o(1))\log N\), so \(N=X^{L/6+o(1)}\)
as the height tends to infinity at fixed \(L\).
The carry-free construction also has \(|C|=(2H-1)^d=N^{1+o(1)}\).
[High-dimensional complete-letter construction](ABC_HIGH_DIMENSIONAL_ADDITIVE_FEATURE_OBSTRUCTION.md#exact-primitive-graphs-at-the-endpoint-cardinality-scale)

The central target \(T=a_0+b_0+s_*\) has \(N\) candidate
edges. Pair the digit vector \(\boldsymbol x\) of the first
summand with \((H-1)\boldsymbol 1-\boldsymbol x\) on the
second. We can choose the shifts so there are at least
\(N^2/100\) primitive edges globally and at least \(N/100\)
primitive edges in this central fibre. For independently
uniform \(a_0,b_0\) the nonprimitivity probability of any
fixed digit pair, including each central pair, is at most
\[
 \sum_{p\le 3M}(1/p+1/M)^2
 \le \frac14+\sum_{\substack{j\ge3\\j\ {\rm odd}}}\frac1{j^2}
            +O\!\left(\frac{\log M}{M}\right)
 = \frac{\pi^2}{8}-\frac34+o(1)<\frac12 .
\]
The first bound follows by counting each summand residue
independently in the shift intervals. By Markov's inequality
the failure probability of either lower count is at most
\((\pi^2/8-3/4+o(1))/(99/100)\). Their sum is less than one.
Take shifts for which both counts hold. Each active complete
dictionary then has at least \(N/100\) values because every
fixed label belongs to at most \(N\) retained edges.

## The fixed fibre has critical triangle entropy

Let \(U\subseteq A\) be the complete first-letter dictionary
on the selected primitive central fibre. It has
\(N/100\le |U|\le N\). Its difference set has cardinality
at most \((2H-1)^d=N^{1+o(1)}\) by the digit difference
representation. Hence the set of ratios
\[
       \left\{\frac{u_2-u_1}{u_3-u_1}:
                           u_1<u_2<u_3,\ u_i\in U\right\}
\]
has cardinality at most \(N^{2+o(1)}\). Each exact ratio
determines the primitive positive additive triple of
normalized distances.

For the matching lower bound, pin \(u_1\) at the smallest
element of \(U\). The positive shifted set
\(V=\{u-u_1:u\in U,\ u>u_1\}\) has size \(\gg N\).
Its multiplicative energy is at most \(|V|^2X^{o_X(1)}\)
by fixing the common product and counting its divisor
pairs. Cauchy--Schwarz shows that \(|V/V|\ge N^{2-o(1)}\),
and half of the nonunit ratios lie in \((0,1)\) by
inversion. These are pinned triangle ratios. At fixed \(L\),
\(X^{o_X(1)}=N^{o(1)}\), so the two bounds give normalized
triangle dictionary exponent exactly \(2\) in base \(N\).
Also every fixed \(k\)-fold sumset \(kU\) has size at most
\((kH)^d=N^{1+o(1)}\). In particular all fixed-order
additive energies have their maximal exponent.

## The central fibre still avoids small feature models

Take any nonempty \(U'\subseteq U\), and let \(d'\) be the
affine dimension of its digit vectors. Projecting onto
suitable \(d'\) digit positions is injective on their affine
span, so \(|U'|\le H^{d'}\). Pick \(d'+1\) affinely independent
digit vectors. Their multisets of \(m\) terms give distinct
vector sums. All digits in those sums are less than \(Q\).
Consequently, for \(d'>0\),
\[
       |mU'|\ge {m+d'\choose d'}
             \ge (m/d')^{d'}
             \ge H^{td'}\ge |U'|^t .
\]
The same end inequality holds for \(d'=0\).

Suppose \(U'\) is an affine image of integer feature vectors
of affine rank \(r\le d\). Let \(\Delta\) bound the absolute
values of their maximal anchored coordinate minors. Choose
a maximal minor and project the feature span to its columns.
In the basis given by that minor, all differences from the
chosen anchor have coordinates in \([-1,1]^r\) by Cramer's
rule. Their coordinate lattice contains \(\mathbb Z^r\) with
index at most \(\Delta\), so the \(m\)-fold sums have cardinality
at most \(\Delta(2m+1)^r\). The same upper bound holds for
their affine images, including rank zero with \(\Delta=1\).
This is the maximal-minor argument of (HF1).
[Long sumsets and maximal minors (HF1)](ABC_HIGH_DIMENSIONAL_ADDITIVE_FEATURE_OBSTRUCTION.md#maximal-feature-minors-control-long-sumsets)

If \(\Delta\le X^\theta\), comparison of the two sumset
bounds gives
\[
         \frac{\log |U'|}{\log N}
              \le \frac{t+1}{t}\theta+\frac rd+o(1).
\]
The error is uniform for \(r\le d\). In particular, for
every fixed rank and \(\sigma>(t+1)\theta/t\), no such model
contains \(N^\sigma\) values of the central fibre for large
height. A subpower cover by models with the same rank and
minor bounds does not change this conclusion.

The construction concerns exact addition, primitivity,
height, triangle-ratio cardinalities and affine features.
An application to the abc endpoint would require the
original radical caps or a restriction on the radicals of
the normalized triangle sides. The ratio cardinality alone
does not encode either arithmetic input.
