# Straight pair lines with a bounded Gaussian spectator

Let \(A\) be an actual balanced planar Gaussian three-cell
noise-stability maximizer at \(0<r<1\). Write \(U_\ell=T_r1_{A_\ell}\)
and \(V_\ell=U_\ell-\beta_\ell\) for its raw and adjusted scores.
Assume \(A_i\) is bounded up to null sets. The two remaining
labels are \(j,k\).

## One straight interface segment gives two straight rays

Suppose an open segment of an actual \(j,k\) interface lies
on an affine line \(L\). Use orthonormal coordinates
\((u,v)\) about the Gaussian origin with \(L=\{v=b\}\).
The stationary contrast on the segment is
\(\lambda=\beta_j-\beta_k\), so analyticity makes
\(U_j-U_k=\lambda\) on all of \(L\). Along this line,
\[
 V_i-V_j=V_i-V_k
   =\frac32(U_i-\theta),\qquad
 \theta=\frac{1+\lambda+2(\beta_i-\beta_j)}3
        =\frac{1+2\beta_i-\beta_j-\beta_k}3.                 \tag{BS1}
\]
Since \(A_i\) has positive mass, the positive Mehler kernel gives
\(U_i>0\) everywhere. If \(\theta\le0\), label \(i\) is a strict
winner along the whole line. Its strict-winner set is open and
agrees with \(A_i\) almost everywhere, contradicting boundedness.
Thus \(\theta>0\). Boundedness of \(A_i\) also gives
\(U_i(u,b)\to0\) as \(|u|\to\infty\). Hence the spectator loses
strictly on both exterior rays of \(L\).

Each point on these exterior rays is an actual regular
\(j,k\) interface point. In a neighborhood only \(j,k\) can
win. If there were no strict wins of one of them arbitrarily
close, the point would be a dormant tie in a ball occupied
almost everywhere by the other. The dormant interior tie
bound (FC4) makes such ties isolated, in conflict with the
whole nearby segment of \(L\) being tied. When both strict
winners approach the point, local analytic stratification
gives a jump arc approaching it. The gradient bound (G20)
keeps the pair gradient nonzero there. The local score level
set is therefore regular and is the line \(L\). Orientation
is constant along each exterior ray.
[Dormant interior ties, FC4](THREE_CELL_FIRST_CONTACT_VARIATIONS.md#point-transfers-at-a-dormant-interior-score-tie)
[Analytic jump arcs and G20](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#regular-interface-arcs-and-isolated-triple-junctions)

No strict-winner component of \(j\) or \(k\) crosses \(L\).
If their strict-winner sets are both connected, the two
components at any regular point on \(L\) show that their
entire cells lie on opposite sides. The affine-separation
theorem makes every strict-winner component unbounded, in
contradiction to boundedness of \(A_i\). Thus in the presence
of a bounded spectator and a straight pair segment at least
one exterior label has a disconnected strict-winner set.
[Affine separation and unbounded strict winners](THREE_CELL_AFFINE_SEPARATION_RIGIDITY.md)

## A binary exterior confines the wrong-side labels

Assume the stronger condition that some half-plane
\(\{u>S\}\) has \(j\) above \(L=\{v=b\}\) and \(k\) below
it, up to null sets. Put \(\sigma=\sqrt{1-r^2}\),
\(g=1_{A_j}-1_{A_k}\), \(s_b(v)=\operatorname{sgn}(v-b)\),
and let \(q_b\) be the normal law of mean \(rb\) and variance
\(\sigma^2\). The exterior determines the pressure on the
whole line. More precisely, factorizing the conditional
Gaussian and using injectivity of the scalar Mehler operator
gives
\[
 \lambda=2\Phi((r-1)b/\sigma)-1,\qquad
 \int_{\mathbb R}[g(u,v)-s_b(v)]\,dq_b(v)=0
                       \quad\text{for almost every }u.       \tag{BS2}
\]
Indeed \(U_j(u,b)-U_k(u,b)=T_r^{[1]}f(u)\) for
\(f(u)=\int g(u,v)\,dq_b(v)\). On the exterior mirror ray
it equals the pair pressure, by continuity of the winner
scores from both sides. Analyticity and injectivity give
\(f\equiv\lambda\), and the exterior binary assignment
evaluates this constant.

If \(k\subseteq\{v<b\}\) up to null sets, write
\(i_\pm(u)=q_b(\{v\gtrless b:(u,v)\in A_i\})\) and
\(j_-(u)=q_b(\{v<b:(u,v)\in A_j\})\). The slice identity
becomes
\[
                 i_+(u)=i_-(u)+2j_-(u).                     \tag{BS3}
\]
Every wrong-side point of \(j\), up to null sets, lies in a
strip of bounded \(u\)-coordinate since \(A_i\) is bounded.
In particular the binary assignment holds in both parallel
tails. The same conclusion holds with the labels and sides
reversed. A connected strict-winner set for \(k\) supplies
the side condition because it meets the pure forward
binary half-plane and cannot cross \(L\).

Under the same forward exterior hypothesis, any straight
segment of the \(j,k\) interface anywhere lies on \(L\).
A nonparallel supporting line would extend toward
\(u\to+\infty\) with \(v\to+\infty\) or \(v\to-\infty\).
On that line the contrast would tend to \(1\) or \(-1\),
contrary to its constant value \(\lambda\in(-1,1)\).
For a parallel supporting line \(v=b'\) the limiting
forward contrast is \(2\Phi((rb'-b)/\sigma)-1\). It equals
\(\lambda\) only when \(b'=b\). These limits follow since
the conditional probability of source points outside the
forward binary half-plane tends to zero uniformly in the
transverse coordinate.

## Geometric asymptotic separation gives a centered binary exterior

With \(A_i\) bounded, let
\[
 D=\{(u,v):g(u,v)\ne\operatorname{sgn}(v)\}.
\]
Suppose
\[
   \lim_{R\to\infty}\sup_{y_1>R}|D\cap B_1(y)|_{\mathbb R^2}=0.
                                                                  \tag{BS4}
\]
Finite Euclidean area of the discrepancy \(D\) in some forward
half-plane implies this condition. Then the binary assignment
is exact on a late forward half-plane, with dividing line \(v=0\).
Thus a bounded spectator with (BS4) cannot be bisected by
\(\{v=0\}\) in Gaussian mass. At least one of the other two
strict-winner sets is disconnected.

For the proof let \(\psi(v)=2\Phi(rv/\sigma)-1\), the smoothed
centered binary contrast. With \(F=T_rg\), conditional Gaussian
smoothing gives
\[
        \sup_{v\in\mathbb R}|F(u,v)-\psi(v)|\longrightarrow0
                                      \qquad (u\to+\infty).  \tag{BS5}
\]
For a fixed large radius the conditional source ball about
\((ru,rv)\) can be covered by a fixed number of unit balls.
Their first coordinates tend to infinity uniformly in \(v\), and
the conditional density is bounded. Apply (BS4) on these balls,
then send the chosen radius to infinity to remove the Gaussian
tail.

Since the spectator is bounded, it is absent on a late full
forward half-plane. If the pair pressure \(\lambda\) were
positive, choose \(c>0\) with \(\psi(c)<\lambda\).
Then (BS5) forces \(k\) almost everywhere at all \(v<c\) on
that half-plane. At \((u,0)\) the conditional source tends
to lie in this forward half-plane and its transverse law is
\(N(0,\sigma^2)\). It follows that
\(\limsup_{u\to\infty}F(u,0)\le1-2\Phi(c/\sigma)<0\),
contrary to (BS5). Reversing the normal and exchanging the
two labels excludes negative pressure. Hence \(\lambda=0\).

The second transverse derivatives of \(F\) and \(\psi\) are
uniformly bounded, by differentiation of the conditional normal
kernel. If the error in (BS5) is at most \(\epsilon\), a
first difference of length \(\sqrt\epsilon\) bounds their
first derivative difference by \(C_r\sqrt\epsilon\). For all
large \(u\), \(F\) therefore has positive transverse derivative
on \(|v|\le1\) and the sign of \(v\) on \(|v|\ge1\). Its zero
is one bounded analytic graph in the entire late half-plane.
The bounded isolated graph straightening theorem makes it
an exact affine ray, whose height must be zero by (BS5).
The straight-line result above gives disconnection, and the
centered bisection obstruction excludes the bisected spectator.
[Straightening a bounded isolated graph](THREE_CELL_TAIL_FLATNESS.md)
[Centered bisection and a one-sided binary exterior](THREE_CELL_REFLECTION_MIRROR_COMPLETION.md#a-bisected-cell-and-a-one-sided-binary-exterior)

Macroscopic convergence of the cells after dilation gives less
than (BS4). At every positive correlation, a stationary binary
example has three parallel walls in a strip, while its dilates
converge to the two half-planes. This example concerns stationarity.
[Binary stationary strip with three walls](THREE_CELL_ISOLATED_END_RIGIDITY.md#a-binary-stationary-strip-with-three-walls)

## A pressure criterion for an affine bisection

Suppose again that \(U_j-U_k=\lambda\) along an affine line
\(L=\{v=b\}\), and suppose \(L\) bisects \(C=A_i^c\) in
ordinary Gaussian mass. Define \(g_*(u,v)=s_b(v)1_C(u,v)\)
and \(f_*(u)=\int g_*(u,v)\,dq_b(v)\). If
\[
                \left|\int f_*\,d\gamma_1\right|\ge|\lambda|,
                                                                  \tag{BS6}
\]
then \(j,k\) lie on opposite sides of \(L\), up to null sets and
interchange of labels. By affine separation this is impossible
for a bounded spectator.

To verify (BS6), put \(a=r/\sigma^2\). The joint density of the
two transverse observations factors as
\[
   \sigma e^{-rb^2/(1+r)}
       dq_b(v)dq_b(w)\,\exp\bigl(a(v-b)(w-b)\bigr).
                                                                  \tag{BS7}
\]
Subtract the transverse constant term \(1\). In the remainder
absorb \(s_b(v)\) into \(g\); the resulting kernel has positive
factor \((e^{a(v-b)(w-b)}-1)s_b(v)s_b(w)\) off the line in
both variables. Replacing \(s_b g\) by its absolute value
increases the remainder strictly if both signs have positive
mass. In the constant term the original normal average is
identically \(\lambda\) by scalar Mehler injectivity. The exact
change in that term, before the factor \(1/2\) converting
contrast energy to partition stability, is
\[
  \sigma e^{-rb^2/(1+r)}
       \bigl(\langle f_*,T_r^{[1]}f_*\rangle-\lambda^2\bigr).
                                                                  \tag{BS8}
\]
It is nonnegative under (BS6) because the scalar Mehler
operator is positive semidefinite and preserves constants.
Bisection makes \(g_*\) an admissible balanced split of \(C\).
Maximality therefore forces \(s_b g\) to have one sign, proving
separation. The same argument uses directly the weaker
condition \(\langle f_*,T_r^{[1]}f_*\rangle\ge\lambda^2\).
For a centered exact binary exterior \(\lambda=0\), so no
additional pressure estimate is needed.
[Centered constant-term subtraction](THREE_CELL_COMPACT_SIMPLEX_DISCREPANCY.md#a-spectator-bisector-forces-the-simplex)
[Centered mirror pressure bound](THREE_CELL_REFLECTION_MIRROR_COMPLETION.md#constant-mirror-pressure-and-a-bounded-invariant-cell)
