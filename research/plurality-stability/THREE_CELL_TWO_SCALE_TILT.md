# A bounded low-correlation perturbation of a three-cell tie

Let \(T\) be a centered regular three-sector Gaussian partition.
Work with ordinary balanced partitions \(A\) of the Gaussian plane.
Write
\[
 S_t(A)=\sum_i\langle 1_{A_i},T_t1_{A_i}\rangle,\qquad
 D_t(A)=S_t(T)-S_t(A),
\]
and let \(d(A)\) be Gaussian label distance to the simplex orbit,
using the sum of the three cell symmetric differences. Suppose
\(r\in(0,1)\) is a correlation where \(T\) and a nonsimplex
partition both maximize. Fix a compact positive correlation band
\(J\) with \(r\) in its interior, and choose
\(0<s<\min J\) in an interval where \(T\) is the unique optimizer.
The exact lower interval provides such a choice.
[Exact lower interval through \(339/500\)](THREE_CELL_CENTRAL_SHIFT_BOOTSTRAP.md)

The parameter \(\eta\) in \(S_t-\eta S_s\) gives a locally
Lipschitz threshold curve for coexistence with the simplex.
For sufficiently small fixed \(s\), tilting at \(r\) toward
lower \(S_s\) also forces a quantitative change in the optimizing
partition. These results concern tilted optimization in the plane.
They leave the existence of the scalar tie open.

## The threshold and its one-sided slopes

For \(t\) near \(r\) define
\[
 \lambda(t)=\inf_{d(A)>0}\frac{D_t(A)}{D_s(A)} .
                                                               \tag{TS1}
\]
The denominator is strictly positive. The simplex is optimal
for \(S_t-\eta S_s\) exactly when \(\eta\le\lambda(t)\).
It is unique up to its orbit for \(\eta<\lambda(t)\).
At \(\eta=\lambda(t)\) there is a nonsimplex maximizer whose
orbit distance is bounded below uniformly for \(t\) near \(r\).

To prove this, first note that there are \(\delta,c_0>0\)
such that
\[
             D_t(A)\ge c_0 D_s(A)
           \quad\text{if }0<d(A)<\delta,\ t\in J.           \tag{TS2}
\]
Use the uniform measurable local comparison on \(J\),
which gives \(D_t(A)\ge c L_t(A;T_A)\) for \(d(A)<\delta\).
The orientation \(T_A\) is the moment alignment in that theorem,
and \(L_t\) integrates its pointwise score loss on changed labels.
For \(x\) in the reference cell \(i\), the gap of its score
over label \(j\) is comparable to \(\min(a,1)\), where \(a\)
is the distance to the pair line. The constants are uniform
for \(t\in J\cup\{s\}\). For the upper bound use the bounded
gradient of the conditional Gaussian convolution and its
vanishing on the pair line, together with the pointwise bound
by one. For the lower bound, in the coordinates of (CT2) one
has \(b\ge-a/\sqrt3\), where \(b\) is the coordinate along
the outgoing ray. Put \(u=t/\sqrt{1-t^2}\). Restriction of the
nonnegative integral there to \(v\ge ua\) gives a lower bound
\((\Phi(2ua)-1/2)/2\). This gives \(L_s(A;T_A)\le C L_t(A;T_A)\).
The deficit identity \(D_s=2L_s-Q_s\) and positivity of the
Ornstein--Uhlenbeck quadratic interaction \(Q_s\) prove (TS2).
[Uniform measurable local comparison from the interface gap](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison)
[Three-cell interface gap](SIMPLEX_INTERFACE_SYMMETRY.md#reflected-kernels-close-the-three-cell-spectral-problem)
[Pointwise score formula, CT2](THREE_CELL_CORRELATION_CONTACT.md#pointwise-score-growth)

On the distant set \(d(A)\ge\delta\), the denominator \(D_s(A)\)
is bounded below by a positive constant. Otherwise a maximizing
sequence at \(s\) would converge strongly along a subsequence
to a simplex, by positive-correlation maximizer compactness and
the uniqueness at \(s\). On this set the functions
\(t\mapsto D_t(A)/D_s(A)\) have uniformly bounded first and
second derivatives on \(J\). Indeed the first two correlation
derivatives of \(T_t\) have uniformly bounded Hermite multipliers
on \(J\), and the indicator norms have fixed sum.

A frozen nonsimplex tie gives \(\lambda(t)\le C|t-r|\)
near \(r\). The distant ratios are bounded below by
\(-C|t-r|\), since \(D_r(A)\ge0\) and their derivatives
are bounded. Thus for \(t\) near \(r\) the infimum in (TS1)
is the infimum over \(d(A)\ge\delta\), by (TS2), and
\(|\lambda(t)|\le C|t-r|\). An infimum of functions with
a common Lipschitz bound is Lipschitz. The upper bound on their
second derivatives also makes \(\lambda\) semiconcave here.
Indeed after subtracting a common multiple of \(t^2/2\)
the functions are concave, and the infimum of concave
functions is concave.

Here is a compactness proof that does not require the distant
set to be weakly closed. For fixed \(t\) near \(r\), the operator
\(K=T_t-\lambda(t)T_s\) is compact and strictly positive
definite. Its Hermite multipliers are \(t^n-\lambda(t)s^n>0\)
for \(n\ge0\), since \(s<t\) and \(\lambda(t)<1\).
The set of balanced fractional assignments
\[
       0\le f_i\le1,\qquad \sum_i f_i=1,\qquad
                   \int f_i\,d\gamma_2=1/3
\]
is weakly compact in the product of three \(L^2\) spaces.
The quadratic objective with \(K\) is weakly continuous by
compactness and strictly convex by positive definiteness.
Every maximizing fractional assignment is ordinary. If a
fractional assignment is nonordinary, some two labels both
exceed a fixed positive margin on a positive-measure set.
Divide that set into two parts of equal positive measure and
vary the two labels oppositely, with opposite signs on the two
parts. For sufficiently small amplitude this gives two distinct
feasible assignments whose average is the original assignment,
contradicting strict convexity and maximality. Hence the
fractional and ordinary maxima agree.

By definition of (TS1), the simplex attains the ordinary tilted
maximum at \(\lambda(t)\). A distant ratio-minimizing sequence
has tilted values tending to the simplex value, because
\(D_s\) is bounded. Any weak limit of a subsequence attains the
maximum over fractional assignments and is ordinary. The fixed
sum of the squared indicator norms then upgrades weak convergence
to strong convergence. Since orbit distance is Lipschitz in
Gaussian label distance, the limit is distant and attains
(TS1). The same argument applies to tilted maximizers \(A_h\)
chosen at \(t=r+h\), with \(d(A_h)\ge\delta\), as \(h\to0\).
The operator converges in norm to \(T_r\), and a weak limit
attains the original maximum at \(r\). It is ordinary, and
strong convergence gives a nonsimplex tie in the limit.
[Compactness of Gaussian maximizers](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)

Write \(\mathcal T_r\) for the compact family of nonsimplex
tied maximizers at \(r\). The one-sided slopes of the threshold
are
\[
 \lambda'_+(r)=
     \min_{A\in\mathcal T_r}\frac{D'_r(A)}{D_s(A)},\qquad
 \lambda'_-(r)=
     \max_{A\in\mathcal T_r}\frac{D'_r(A)}{D_s(A)} .           \tag{TS3}
\]
The derivatives in the ratios hold the partition fixed.
For any \(A\in\mathcal T_r\), Taylor expansion yields
\(\lambda(r+h)\le h D'_r(A)/D_s(A)+O(h^2)\).
For a distant tilted maximizer \(A_h\) attaining (TS1) at \(r+h\),
the same expansion, using \(D_r(A_h)\ge0\), gives
\(\lambda(r+h)\ge h D'_r(A_h)/D_s(A_h)+O(h^2)\).
The errors are uniform on the distant set, and any sequence
of the chosen \(A_h\) has a strongly convergent subsequence
with limit in \(\mathcal T_r\). The ratio of derivatives passes
to that limit. Taking positive \(h\) yields the minimum in
(TS3), and taking negative \(h\) reverses the inequalities and
yields the maximum. If simplex-optimal correlations approach
\(r\) from the left, each ratio in (TS3) is nonpositive.
If they approach from the right, each is nonnegative.
In particular, approach from both sides forces both first
slopes to vanish.

## Strict selection and shape motion at the tie

Now choose \(s>0\) sufficiently small in the proved simplex
interval. Let \(\mathcal M_r\) be the whole set of ordinary
balanced planar maximizers at \(r\), including the simplex,
and put \(V=S_r(T)\) and
\(\ell=\min_{A\in\mathcal M_r} S_s(A)\).
The minimum is attained by compactness and occurs only
at nonsimplex tied maximizers by strict uniqueness at \(s\).
For sufficiently small \(\eta>0\) let \(A_\eta\) be any
planar maximizer of \(S_r-\eta S_s\). Such a maximizer is
ordinary by the fractional argument above. There are
\(c_1,c_2,c_3>0\), independent of the small \(\eta\) and of
the chosen maximizing shape, such that
\[
\begin{split}
 \max_A[S_r(A)-\eta S_s(A)]&\ge V-\eta\ell+c_1\eta^2,\\
             c_2\eta^2&\le V-S_r(A_\eta)=o(\eta),\\
 \operatorname{dist}(A_\eta,\mathcal M_r)&\ge c_3\eta .
\end{split}                                                     \tag{TS4}
\]
Distance here is the infimum of the sum of label symmetric
differences over all members of \(\mathcal M_r\).
Their rotations and relabelings already belong to
\(\mathcal M_r\).

Use the compact exactly balanced stationary test from (SCW1),
half the first shape derivative. Its uniform witness construction
gives finitely many neighborhoods with compact test fields
covering the nonsimplex members of \(\mathcal M_r\). In each neighborhood a
fixed compact field is corrected by the same two compact flux
fields, with coefficients depending on the tested partition
\(B\), to have zero first mass flux. Denote the resulting
functional by \(H_B(z)\), with the corresponding neighborhood
understood. The witness satisfies
\(|(H_B/z)(0)|\ge c>0\) at the nonsimplex maximizers after
choosing a neighborhood for each. The finite family gives a
common bound for \(H_B(z)/z\) on a closed disk about zero
of fixed positive radius. The power series estimate on a
smaller disk then gives \(|H_B(s)|\ge c s/2\) for a
sufficiently small common \(s>0\), at these maximizers and,
after shrinking the witness neighborhoods, at nearby \(B\).

In these neighborhoods the field coefficients depend
Lipschitz continuously on Gaussian label distance. The mass
fluxes are integrals of fixed bounded functions, and the
correction uses uniformly invertible matrices. The value
\(H_B(t)\) is Lipschitz in the same distance for \(t=r,s\),
since the Mehler kernel and its first observation derivatives
are uniformly bounded in the source variable on the compact
observation support of the fields. The two flux fields also
give an exact balancing correction of order \(O(u^2)\) for
the flow of amplitude \(u\) in the corrected field. Over the
finite family the existence interval and the quadratic
remainder bounds are uniform, after shrinking the neighborhoods.
For \(t=r,s\), the flow changes \(S_t\) by
\(2uH_B(t)+O(u^2)\). These facts follow by differentiating
the transported Gaussian densities on the fixed cells and by
the implicit function theorem for the two independent masses.
[Uniform stationary witnesses, SCW1--SCW4](GAUSSIAN_STATIONARY_CORRELATION_WITNESS.md#a-shape-derivative-with-a-nonzero-degree-one-term)
[Uniform mass repair, assignment gain and Lipschitz response, SR3--SR4](GAUSSIAN_STATIONARY_CORRELATION_WITNESS.md#approximate-stationarity-and-motion-with-correlation)

For the first inequality in (TS4), take a frozen \(A_*\in
\mathcal M_r\) attaining \(\ell\). Then \(H_{A_*}(r)=0\)
and \(H_{A_*}(s)\ne0\) in one of its witness neighborhoods.
Flow in its exactly balanced test by
\(u=-\kappa\eta H_{A_*}(s)\), with \(\kappa>0\) fixed
and small. The change of the tilted value is
\(2\kappa\eta^2 |H_{A_*}(s)|^2+O(\kappa^2\eta^2 |H_{A_*}(s)|^2)\).
This proves the quadratic gain.

Comparison of \(A_\eta\) with \(A_*\) gives
\[
             0\le V-S_r(A_\eta)
                     \le\eta(\ell-S_s(A_\eta)).
\]
It follows that \(S_r(A_\eta)\to V\). Any small-\(\eta\)
sequence has a strongly converging subsequence to a member of
\(\mathcal M_r\), and its limit minimizes \(S_s\) there.
This also proves \(V-S_r(A_\eta)=o(\eta)\) uniformly over
the choices of \(A_\eta\). In particular every \(A_\eta\)
lies near the compact set of such nonsimplex minimizers.
Use a witness neighborhood there. Exactly balanced tilted
stationarity gives \(H_{A_\eta}(r)=\eta H_{A_\eta}(s)\).
Consequently \(|H_{A_\eta}(r)|\ge c'\eta\), with \(c'>0\)
uniform. Flowing in its test direction by
\(u=\kappa H_{A_\eta}(r)\) gains at least
\(c''|H_{A_\eta}(r)|^2\) for the ordinary objective \(S_r\)
when \(\kappa>0\) is fixed small. The definition of \(V\)
proves the lower bound on the ordinary deficit in (TS4).

For the distance bound, choose a closest \(B_\eta\in\mathcal M_r\).
It is close to the same compact selected set for small \(\eta\).
A finite cover with smaller neighborhoods inside the witness
neighborhoods permits a common test neighborhood for
\(A_\eta,B_\eta\). In that neighborhood \(H_{B_\eta}(r)=0\).
The Lipschitz estimate for \(H(r)\) and
\(|H_{A_\eta}(r)|\ge c'\eta\) prove the last inequality in (TS4).

## Subtracting first Hermite energy directly

The bounded objective \(S_s\) can be used in (TS2)
because its pointwise score contrast is bounded by a multiple
of the reference \(S_t\) contrast. Direct subtraction of
first Hermite energy has different local behavior.

Write \(W_1(A)=\sum_i|m_i(A)|^2\) for
\(m_i(A)=\int_{A_i}x\,d\gamma_2\). For every \(t\in(0,1)\)
and \(\varepsilon>0\), the simplex fails to be a measurable
local maximizer of \(S_t-\varepsilon W_1\).
Indeed, choose a pair \(i,j\) and reflect a ball \(E_R\) of
radius one, centered at \(R\) times the central unit direction
of \(T_i\), in the \(i,j\) reference line. Its image \(E'_R\)
lies in \(T_j\), and \(E_R\subset T_i\) for sufficiently
large \(R\). Swap the labels on the two balls to get \(A_R\).
The reflection gives equal masses \(\mu_R=\gamma_2(E_R)\)
and exact balance. Write
\(v_R=\int_{E_R}x\,d\gamma_2-\int_{E'_R}x\,d\gamma_2\).
The moment difference \(m_i(T)-m_j(T)\) is normal to the
reflecting line, and has inner product at least \(cR\) with
every \(x\in E_R\) for large \(R\), with \(c>0\) fixed.
Thus
\[
\begin{split}
 W_1(T)-W_1(A_R)
   &=2(m_i(T)-m_j(T))\cdot v_R-2|v_R|^2\\
   &\ge 4cR\mu_R-8(R+1)^2\mu_R^2 .
\end{split}
\]
Since \(R\mu_R\to0\), this is at least \(2cR\mu_R\) for
large \(R\). By the deficit identity \(D_t(A_R)\le
2L_t(A_R;T)\le4\mu_R\). It follows that \(A_R\) improves
\(S_t-\varepsilon W_1\) for sufficiently large \(R\).
Its label distance to \(T\) tends to zero. For
\(0<\varepsilon<t\) this failure occurs even though the
tilted operator \(T_t-\varepsilon\Pi_1\) is strictly
positive definite.
