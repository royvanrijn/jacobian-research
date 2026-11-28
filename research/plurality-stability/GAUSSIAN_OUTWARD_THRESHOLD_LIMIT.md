# Discrete exterior limits for stationary Gaussian partitions

Let \(A\) be a finite Gaussian partition in \(\mathbb R^d\), write
\(f_i=1_{A_i}\), and fix \(0<\rho<1\). Suppose it satisfies the
stationary winner rule with constant prices \(\beta_i\):
\[
 f_i(X)>0\quad\Longrightarrow\quad
 T_\rho f_i(X)-\beta_i=\max_j(T_\rho f_j(X)-\beta_j)
 \qquad\text{a.e.}
\]
In particular the hypothesis holds for Gaussian noise-stability
maximizers. For almost every \(x\ne0\), the entire sequence of
physical-size windows at \(\rho^{-n}x\) has a single limit.
The translated indicators converge weak-* to a spatially
constant probability vector \(m(x)\). That vector is one of
finitely many solutions of
\[
     m_i>0\quad\Longrightarrow\quad
        m_i-\beta_i=\max_j(m_j-\beta_j).                  \tag{OT1}
\]
The sequence of dilated partitions \(f(\rho^{-n}\cdot)\)
converges weak-* to the measurable map \(m\), which satisfies
\(m(\rho x)=m(x)\) almost everywhere. The limiting map may
take fractional values. Every nonzero spatial derivative
of the smoothed scores has finite exterior square integral
with weight \((1+|X|)^{-d}\).

The proof uses the full orbit equation. After rescaling, the
winner rule is exactly one step of Gaussian threshold dynamics
with a geometrically decreasing variance and unchanged prices.
The thresholding energy calculation is related to the
minimizing-movements formulation for fixed variance. Here the
energy is unnormalized, and the decrease of variance supplies
a second, frequency-sensitive term in its dissipation.
[Esedoğlu–Otto, Proposition 1 and Section 8](https://websites.umich.edu/esedoglu/Papers_Preprints/esedoglu_otto.pdf)

## Energy along the outward orbit

Let \(H_h\) denote convolution with the centered Gaussian of
covariance \(hI\). Put \(\sigma^2=1-\rho^2\),
\(u_n(x)=f(\rho^{-n}x)\) and \(h_n=\sigma^2\rho^{2n-2}\), for
\(n\ge0\). Since \(T_\rho f(X)=(H_{\sigma^2}f)(\rho X)\),
\[
 u_n(x)\cdot\bigl(H_{h_n}u_{n-1}(x)-\beta\bigr)
 \ \ge\
 u_{n-1}(x)\cdot\bigl(H_{h_n}u_{n-1}(x)-\beta\bigr)
 \quad(n\ge1).                                           \tag{OT2}
\]
In fact \(u_n(x)\) maximizes this linear functional on the
label simplex at almost every \(x\).

Fix a smooth compactly supported cutoff \(0\le\chi\le1\).
Use the localized energy
\[
 E_h(u)=\sum_i\langle \chi u_i,(I-H_h)\chi u_i\rangle_{L^2(dx)}
                  +2\int \chi^2\,\beta\cdot u\,dx .
\]
For a scalar bounded function \(g\), the commutator satisfies
\[
 \|\chi H_hg-H_h(\chi g)\|_1
 \le\|g\|_\infty\|\nabla\chi\|_1\,\mathbb E|Z_h|,
 \qquad Z_h\sim N(0,hI).
\]
This follows by integrating
\(\|\chi(\,\cdot+z)-\chi\|_1\le |z|\|\nabla\chi\|_1\)
against the Gaussian. Expand \(E_{h_n}(u_n)-E_{h_n}(u_{n-1})\).
The squared norms \(\sum_i\|\chi u_{n,i}\|_2^2\) cancel since
the partitions take values at label vertices. Set
\(v_n=\chi(u_n-u_{n-1})\). The interaction difference is
\[
 -2\sum_i\langle v_{n,i},H_{h_n}(\chi u_{n-1,i})\rangle
 -\sum_i\|H_{h_n/2}v_{n,i}\|_2^2 .
\]
Move \(\chi\) across \(H_{h_n}\) in the first term at cost
\(O(\sqrt{h_n})\). The resulting cross term plus the price
change is nonpositive by (OT2). Consequently
\[
 E_{h_n}(u_n)+\sum_i\|H_{h_n/2}v_{n,i}\|_2^2
       \le E_{h_n}(u_{n-1})+C\sqrt{h_n}.                 \tag{OT3}
\]

The Fourier multiplier of \(H_h\) is \(e^{-h|\xi|^2/2}\).
Since \(h_n=\rho^2 h_{n-1}\), the difference
\(E_{h_{n-1}}(u)-E_{h_n}(u)\) is nonnegative. Its multiplier,
with \(t=h_n|\xi|^2\), is
\[
 e^{-t/2}-e^{-t/(2\rho^2)} .
\]
For every integer \(j\ge1\) this dominates \(c_{\rho,j}t^j e^{-t}\).
Near zero it is comparable to \(t\), and near infinity it is
comparable to \(e^{-t/2}\). The localized energy is bounded
below by \(-2\|\beta\|_\infty\|\chi\|_2^2\). Therefore (OT3),
the nonnegative scale difference and \(\sum_n\sqrt{h_n}<\infty\)
give
\[
\begin{split}
 &\sum_{n\ge1}\sum_i\|H_{h_n/2}(\chi(u_n-u_{n-1})_i)\|_2^2<\infty,\\
 &\sum_{n\ge1}\sum_i h_n^{|\alpha|}
       \|D^\alpha H_{h_n}(\chi u_{n-1,i})\|_2^2<\infty
                         \qquad(|\alpha|\ge1).
\end{split}                                             \tag{OT4}
\]

On any compact set \(V\) where \(\chi=1\) on a neighborhood,
one can remove \(\chi\) in (OT4). The error in its first line
is bounded pointwise by a Gaussian tail. After multiplication
by \(h_n^{|\alpha|/2}\), the derivative error in its second
line is bounded by a polynomially weighted Gaussian tail.
Both are summable in squared norm on \(V\), since the tails
start at a fixed positive distance from \(V\).
The derivative estimate on \(V\) now reads, after scaling,
\[
 \sum_{n\ge0}\int_V
     |D^\alpha H_{\sigma^2} f_i(\rho^{-n}x)|^2\,dx<\infty
                    \qquad(|\alpha|\ge1).                 \tag{OT5}
\]
For \(V=\{1\le |x|\le\rho^{-1}\}\), the dilated shells cover
the exterior with disjoint interiors. Changing variables
in (OT5) proves the weighted square integrability for
\(H_{\sigma^2}f_i\), hence also for \(T_\rho f_i\).

## A unique equilibrium on almost every orbit

Apply (OT5) on a compact neighborhood of \(V\).
For every fixed window size \(L\) and all sufficiently large
\(n\), the points \(x+\rho^n z\) with \(x\in V\), \(|z|\le L\)
lie in that neighborhood. Translating the integral in \(x\)
and then integrating in \(z\) makes
\[
 \sum_{n\ge0}\int_V\int_{|z|\le L}
    |D^\alpha H_{\sigma^2}f_i(\rho^{-n}x+z)|^2\,dz\,dx
 <\infty .
\]
Use Fubini for a countable sequence of window sizes, followed
by local Sobolev embedding with derivative orders above
\(d/2+1\). For almost every \(x\),
\[
     \nabla H_{\sigma^2}f_i(\rho^{-n}x+\,\cdot)
                          \longrightarrow 0
   \quad\text{uniformly on every bounded window}.          \tag{OT6}
\]
The assertion on all of \(\mathbb R^d\setminus\{0\}\) follows
from a countable exhaustion by compact sets.

At such an \(x\), extract any weak-* limit \(F_i\) of
\(f_i(\rho^{-n}x+\cdot)\) along a subsequence.
The Gaussian smoothing of \(F_i\) is constant by (OT6).
Gaussian convolution is injective on bounded functions.
Indeed the convolution identity for a bounded function
with zero transform says that the bilateral Laplace
transform of its Gaussian-weighted finite measure is zero.
Analytic continuation and Fourier uniqueness force this
measure to vanish. Hence each \(F_i\) is itself constant.

The first line of (OT4), with the cutoff removed, gives for
almost every \(x\)
\[
 H_{\sigma^2/(2\rho^2)} f_i(\rho^{-n}x)
     -H_{\sigma^2/2} f_i(\rho^{-(n-1)}x)\longrightarrow0 .
                                                               \tag{OT7}
\]
Take simultaneous subsequential physical-window limits at
\(\rho^{-n}x\) and \(\rho^{-(n-1)}x\). They are constant, and
(OT7) identifies their vectors. The winner rule at
\(\rho^{-n}x+z\) compares \(H_{\sigma^2}f_i\) at
\(\rho^{-(n-1)}x+\rho z\), so a strict losing gap for the
common limit eliminates that label on any fixed bounded
window. This gives (OT1).

Every equilibrium in (OT1) is determined by its nonempty
support \(J\). Its components are
\(m_i=\beta_i+(1-\sum_{j\in J}\beta_j)/|J|\) on \(J\)
and zero off \(J\). Positivity on \(J\) and the inequalities
outside \(J\) decide whether it is admissible.
Thus the vectors \(b_n(x)=H_{\sigma^2} f(\rho^{-n}x)\)
have a finite cluster set, and \(|b_n-b_{n-1}|\to0\) by
the simultaneous limit argument. A bounded sequence with
these properties converges. Sufficiently small disjoint
balls around the finite cluster set contain the tail, and
steps eventually cannot move between those balls. This
also identifies every translated weak-* limit of the
indicators. Shifting \(n\) proves \(m(\rho x)=m(x)\).

Finally \(b_n=H_{\sigma^2\rho^{2n}}u_n\). Bounded convergence
gives local \(L^p\) convergence of \(b_n\) for finite \(p\),
and testing against \(L^1\) functions with the Gaussian
approximate identity proves \(u_n\rightharpoonup^*m\).
In particular a nonconstant periodic physical-window limit
cannot occur on a set of dilation phases of positive
Lebesgue measure. This leaves exceptional escaping sequences
and fractional equilibria open.
[Periodic limits in fixed exterior windows](THREE_CELL_EXTERIOR_INDEX_AND_DILATION_CHAIN.md#a-periodic-chain-with-three-labels)

## Interior geometry of the pure phases

Assume now that every label has positive Gaussian measure.
Let \(O_i\) be the strict adjusted-score winner set of label
\(i\), and \(D_i(X)=\operatorname{dist}(X,O_i^c)\).
The pair-score contrasts cannot be identically zero, since
the Gaussian operator is injective and the indicator
difference for two positive labels is nonconstant. Their
tie sets are null by analyticity, so \(1_{O_i}=f_i\) a.e.
Every difference \(\beta_i-\beta_j\) is strictly less than
one, since label \(i\) wins somewhere and each unadjusted
score belongs to \((0,1)\).

The ball amplification argument for bounded planar winners
also controls the full sequence of rescaled distance functions.
[Ball comparison, CB2](THREE_CELL_COMPACT_SCORE_BRANCH_OBSTRUCTION.md#a-ball-comparison-for-bounded-winner-regions)

Choose \(c_i>0\) such that
\(\Pr(|Z_{\sigma^2}|<c_i)>(1+\max_{j\ne i}(\beta_i-\beta_j))/2\).
If \(D_i(X)>c_i\), any conditional mean within distance
\(D_i(X)-c_i\) of \(X\) has at least this much conditional
mass in label \(i\). Comparison against every other label
gives
\[
                 D_i(X/\rho)\ge \rho^{-1}(D_i(X)-c_i)_+ .
                                                               \tag{OT8}
\]
The functions \(d_{i,n}(x)=\rho^n D_i(\rho^{-n}x)\) are
1-Lipschitz. They satisfy
\(d_{i,n+1}\ge d_{i,n}-c_i\rho^n\), so
\(d_{i,n}-c_i\rho^n/(1-\rho)\) is nondecreasing. Their local
upper bounds follow by fixing any point of \(O_i^c\).
Thus \(d_{i,n}\) converge locally uniformly to a
1-Lipschitz function \(d_{i,\infty}\). The open sets
\(\{d_{i,\infty}>0\}\) are pairwise disjoint and invariant
under multiplication by \(\rho\). They are exactly the
points where \(D_i(\rho^{-n}x)>c_i/(1-\rho)\) at some \(n\).
One implication follows by summing (OT8), and the other
because \(\rho^n D_i(\rho^{-n}x)\to d_{i,\infty}(x)>0\).

On these open sets every fixed physical window eventually
belongs to \(O_i\), uniformly on compact subsets.
Conversely, at every phase in the preceding almost-everywhere
result where \(m=e_i\), the paired window limits give
strict winning score gaps throughout each fixed bounded
window for all large \(n\). Taking its radius above
\(c_i/(1-\rho)\) and applying (OT8) gives
\(d_{i,\infty}(x)>0\). Thus the pure parts of the limiting
map are these open sets up to null sets.

## Exact interior distances and bounded thick components

Retain the pure open phases \(V_i=\{d_{i,\infty}>0\}\) and put
\(C_i=c_i/(1-\rho)\). The limit in (OT8) is the actual
interior distance of the pure phase:
\[
      d_{i,\infty}(x)=\operatorname{dist}(x,V_i^c),
      \qquad
      D_i(X)\le \operatorname{dist}(X,V_i^c)+C_i.          \tag{OT10}
\]
Indeed \(d_{i,\infty}\) is 1-Lipschitz, vanishes on \(V_i^c\),
and is therefore at most the interior distance. If a closed
ball centered at \(x\) is contained in \(V_i\), the uniform
convergence of \(d_{i,n}\) on the ball implies it is contained
in the \(i\)-th rescaled strict winner for all sufficiently
large \(n\). Thus the limit of the distances at \(x\) is at
least the ball radius. The complement is nonempty. For
example, if \(0\in V_i\), its openness and dilation invariance
would give \(V_i=\mathbb R^d\). The local uniform convergence
would then give \(D_i(0)=+\infty\), since the rescaled distance
at zero would converge to a positive number. In fact
\(d_{i,n}(0)=\rho^n D_i(0)\to0\) by positive masses of the other
labels. Finally (OT8) gives \(d_{i,0}-C_i\le d_{i,\infty}\),
which proves the second assertion in (OT10). Thus any piece
of label \(i\) outside its pure dilated phase has inradius
at each of its points at most \(C_i\).

There is an obstruction to fitting a thick connected region
inside a bounded component of \(O_i\). Let \(K\) be any nonempty
connected measurable subset of \(\{D_i>C_i\}\) lying in such a
component. Then the sets \(\rho^{-n}K\), \(n\ge0\), are pairwise
disjoint and lie in \(O_i\). Membership follows by repeated
application of (OT8). If two copies intersected, scaling back
would give \(K\cap\rho^{-m}K\ne\varnothing\) for some \(m\ge1\).
The union of the connected sets \(\rho^{-jm}K\), \(j\ge0\),
would then be connected and contain \(K\). It would also be
unbounded, since \(K\) contains a nonzero point. The possibility
\(K=\{0\}\) is excluded by \(\{D_i>C_i\}\subset V_i\) and
\(0\notin V_i\). Thus such an overlap is impossible.

Write \(x=e^t\omega\) in log-polar coordinates. For each
\(\omega\) for which the fiber is measurable, the sets of
\(t\)-values in the copies of \(K\) are translates by
\(-n\log\rho\). Their disjointness means the fiber is
injective modulo \(\log(1/\rho)\). Its length is therefore
at most \(\log(1/\rho)\). Integration in \(\omega\) gives
\[
                  \int_K |x|^{-d}\,dx
                 \le |S^{d-1}|\log(1/\rho).               \tag{OT11}
\]
This estimate treats arbitrary connected thick regions,
without an assumption about the shape of their boundary.

As a ball case, if \(B_R(p)\subset O_i\) and
\[
             R>C_i,\qquad
             R-C_i\ \ge\ \frac{1-\rho}{1+\rho}|p|,        \tag{OT12}
\]
then the component containing the ball is unbounded. Indeed
(OT8) gives an \(i\)-ball at \(\rho^{-n}p\) of radius
\(R_n=C_i+\rho^{-n}(R-C_i)\) for each \(n\ge0\).
The distance between consecutive centers is
\(\rho^{-n-1}(1-\rho)|p|\). By (OT12) their radii sum to
strictly more than this, including in the equality case
because \(C_i>0\). The union is connected and unbounded.
Thin components and separated dilated copies of thick
regions are compatible with these inequalities.

## Perimeter growth required by fractional limits

Suppose the original cells have locally finite perimeter and
put \(P(A;U)=\sum_i |Df_i|(U)\). If the limiting map is
fractional on a set of positive Lebesgue measure, then
\[
                    \lim_{R\to\infty}\frac{P(A;B_R)}{R^d}
                       =+\infty .                         \tag{OT9}
\]
Indeed, for almost every such phase \(x\), at least one of
the translated indicators converges weak-* on the unit
ball to a constant strictly between zero and one. Its
perimeter there must tend to infinity. A subsequence with
bounded perimeter would be compact in local \(L^1\) by BV
compactness, and an \(L^1\) limit of indicators is an
indicator. Choose a bounded annulus \(V\) in which the set
of fractional phases has positive measure. Fatou's lemma and
integration in the window centers give, with fixed \(C\)
large enough to contain the windows in \(B_{C\rho^{-n}}\),
\[
 +\infty=\liminf_{n\to\infty}
     \int_V P(A;B_1(\rho^{-n}x))\,dx
 \ \le\ \liminf_{n\to\infty}
     |B_1|\rho^{nd}P(A;B_{C\rho^{-n}}).
\]
For the inequality, any interface point contributes only for
phases in a ball of radius \(\rho^n\). Monotonicity of
perimeter in \(R\) and the fixed ratio between consecutive
radii \(C\rho^{-n}\) give (OT9).

Consequently \(P(A;B_R)=O(R^d)\) is sufficient to ensure that
the discrete exterior limit is an ordinary partition. In
that case weak-* convergence upgrades to local \(L^1\) by
testing each indicator against the corresponding limiting
indicator and its complement on bounded sets. In the
planar maximizing problem the available general interface
growth estimate is \(O(R^4)\), so the sufficient \(O(R^2)\)
growth bound remains open.
[Polynomial growth of the planar network](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#polynomial-growth-of-the-network)
