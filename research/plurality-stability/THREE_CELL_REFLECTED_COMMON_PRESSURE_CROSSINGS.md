# Tail crossings for reflected graphs at a common pressure

Fix \(0<\rho<1\). Consider two stationary reflected three-cell
graphs in the same coordinates, with the same pressure for the
right \(3/1\) interface. On the right half-plane write their
heights as \(h_1,h_2\). Assume they are real analytic and
increasing, have bounded derivatives, tend to \(+\infty\), and satisfy
\[
h_j(\rho x)-\rho h_j(x)=O(1) \qquad (j=1,2).
\tag{1}
\]
These conditions hold for actual maximizers in the
[reflection reduction (RF1)–(RF5)](THREE_CELL_REFLECTION_REDUCTION.md).
If the graphs are distinct and \(h_1-h_2\ge0\) on a tail,
then \(h_1(x)-h_2(x)\ge c x\) there for some \(c>0\).
The reverse statement holds with the graphs interchanged.

In particular, distinct reflected maximizers having the same
pressure and the same limiting phase of (RF5) cross arbitrarily
far out. It suffices that their limiting phases agree at one
positive argument. The proof compares the two fixed score
equations and retains the sign of their difference.

## A positive recurrence

Put \(u=h_1-h_2\), \(h_\tau=h_2+\tau u\) for \(0\le\tau\le1\),
and \(\sigma=\sqrt{1-\rho^2}\). Let \(p_\sigma\) denote the
centered Gaussian density of standard deviation \(\sigma\),
and put \(Q=1-\Phi\) for standard Gaussian distribution
function \(\Phi\).
The pair score \(F_h=T_\rho(1_{A_3(h)}-1_{A_1(h)})\) for a
reflected graph is
\[
F_h(x,y)=
\int_0^\infty B(x,s)Q\bigl((h(s)-\rho y)/\sigma\bigr)\,ds
-\Phi(\rho x/\sigma),\qquad
B(x,s)=2p_\sigma(s-\rho x)+p_\sigma(s+\rho x).
\tag{2}
\]
Conditional on the right horizontal source, the contrast
is twice the third-cell probability minus one. On the left
it is just the third-cell probability. This is the same
folded score used in
[(7) of the reflected affine comparison](THREE_CELL_REFLECTED_AFFINE_WEDGE_COMPARISON.md#a-one-sign-obstruction-for-a-curved-graph).

The common pressure means
\(F_{h_1}(x,h_1(x))=F_{h_2}(x,h_2(x))\). Subtract these
scores by integrating the derivative along \(h_\tau\).
This gives the exact recurrence
\[
\rho u(x)=P_x u,\qquad
P_x(ds)=\frac{L(x,s)}{Z(x)}\,ds,\qquad
L(x,s)=B(x,s)\int_0^1
               p_\sigma(h_\tau(s)-\rho h_\tau(x))\,d\tau,
\quad
Z(x)=\int_0^\infty L(x,s)\,ds .
\tag{3}
\]
For fixed \(x\), differentiation and integration can be
interchanged using the horizontal Gaussian tail and the
linear growth of each height. The linear growth follows also
from (1) and continuity, by iterating toward a fixed compact
interval on the positive half-line.

For the use below, one has uniform constants with
\[
Z(x)\ge c_0>0,\qquad
P_x(ds)\le C_0p_\sigma(s-\rho x)\,ds .
\tag{4}
\]
The bound on the slopes gives \(|h_j(s)-h_j(\rho x)|\le C\) on
\([\rho x,\rho x+1]\). In conjunction with (1), this
bounds both vertical density arguments uniformly on that
interval and gives the first bound in (4). For the second bound,
\(B(x,s)\le 3p_\sigma(s-\rho x)\) on \(x,s\ge0\), and every
vertical density is bounded above by \(p_\sigma(0)\).

## Removing the opposite sign on the bounded part

Suppose \(u\ge0\) on \([R_0,\infty)\). If the graphs
are distinct, real analyticity gives a compact interval
\(I\subset(R_0,\infty)\) of positive length with \(u\ge\delta>0\).
For \(s\in I\), \(v\in[0,R_0]\), and all sufficiently large
\(x\), the horizontal factors satisfy
\[
\frac{B(x,s)}{B(x,v)}
\ge\frac23\exp\left[
\frac{\rho x(s-v)}{\sigma^2}-\frac{s^2-v^2}{2\sigma^2}
\right]\ge c_1e^{c_2x}.
\]
For each \(\tau\) the ratio of the vertical factors at \(s\)
and \(v\) is
\[
\exp\left[
\frac{\rho h_\tau(x)(h_\tau(s)-h_\tau(v))}{\sigma^2}
-\frac{h_\tau(s)^2-h_\tau(v)^2}{2\sigma^2}\right].
\]
The first term is nonnegative, by the ordering of \(s,v\),
monotonicity of \(h_\tau\), and \(h_\tau(x)\ge0\) on the tail.
The second is uniformly bounded on the compact source
intervals. Thus \(L(x,s)\ge c_3e^{c_2x}L(x,v)\)
uniformly over these sources. Integrating and using
\(u\ge\delta\) on \(I\) gives
\[
\int_{0}^{R_0}u_-(v)L(x,v)\,dv
   \ \le\ C e^{-c_2x}\int_I u(s)L(x,s)\,ds .
\tag{5}
\]
In particular (3) yields \(u(x)>0\) for all sufficiently
large \(x\). The strict gain on \(I\) can absorb the
entire negative part while keeping the contribution from
any other interval of positive sources.

## Stopped averaging gives linear separation

Let \(J(x)=[\rho x-\sqrt x,\rho x+\sqrt x]\). For all
sufficiently large \(x\) it lies above \(I\). By (3) and
(5),
\[
\rho u(x)\ge \int_{J(x)}u(s)\,P_x(ds),
\qquad P_x(J(x))\ge 1-C e^{-c x},
\tag{6}
\]
where the mass bound follows from (4).
Choose fixed \(0<l<\rho<v<1\). For a sufficiently large
\(R\), every \(J(x)\) for \(x>R\) is contained in
\([lx,vx]\), and \(u\) has a strictly positive minimum
on \([lR,R]\).

Start at \(x>R\), use the normalized restriction of
\(P_x\) to \(J(x)\) until the first step below \(R\),
and write one path as \(x=x_0,x_1,\ldots,x_n\).
It stops in \([lR,R]\) after a bounded number of
steps for this starting point, since \(x_{j+1}\le v x_j\).
In reverse order \(x_j\ge R v^{-(n-1-j)}\)
for \(j<n\). Thus the product of the retained masses
is bounded below by
\(\prod_{k\ge0}(1-Ce^{-cR v^{-k}})>0\), uniformly in
the starting point and the path. Moreover
\[
\left|\frac{x_{j+1}}{\rho x_j}-1\right|
\le \frac{1}{\rho\sqrt{x_j}} .
\]
The sum of the right sides over all retained steps
is bounded by a geometric series independent of the
starting point. On increasing \(R\) so that each
right side is at most \(1/2\), the product of the
ratios is bounded away from zero and infinity by
the logarithm estimate. It follows that
\(x_n\ge c_4\rho^n x\), and hence
\(\rho^{-n}\ge c_4 x/R\).

Iterating the first inequality of (6) to the stopping
step expresses a lower bound for \(u(x)\) as the
expectation of \(u(x_n)\) times the retained-mass
product times \(\rho^{-n}\). The three bounds give
\(u(x)\ge c x\). The proof for an eventually
nonpositive difference follows by interchanging
the two graphs.

For the reflected maximizers, let
\(H_j(t)=\lim_{n\to\infty}\rho^n h_j(\rho^{-n}t)\)
be the phases of (RF5). Eventual nonnegativity of
a nonzero difference would give \(H_1(t)-H_2(t)\ge
c t\) for every \(t>0\). Eventual nonpositivity gives
the reverse strict separation. Equality of the
phases at even one \(t>0\) therefore forces
arbitrarily late points of both strict signs.
This constrains distinct stationary shapes at common
pressure. It does not force the pressure or phase
of a potential right-contact maximizer to match
those of the simplex.
