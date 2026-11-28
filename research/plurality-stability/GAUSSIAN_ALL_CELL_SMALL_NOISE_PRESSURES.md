# Small-noise score prices for any number of labels

Fix \(k\ge2\), and let \(p\) vary in the interior of the
probability simplex. Write \(V_\rho(p)\) for the maximum
Gaussian noise stability over ordinary partitions of masses
\(p\) in an admissible dimension \(d\ge k-1\). Its value is
independent of \(d\) by dimension reduction.
Let \(I(p)\) be the minimum Gaussian cluster perimeter, counting
each interface once, and put \(e=\sqrt{2(1-\rho)}\).
For any actual maximizer, choose centered score prices \(\beta\)
in the convention that its cells select the greatest of
\(T_\rho1_{A_i}-\beta_i\). Uniformly over compact subsets of
the interior mass simplex and over all actual maximizers,
\[
          -\frac{2\beta}{e}\longrightarrow
                    \sqrt{\frac2\pi}\,\nabla I(p)
                    \qquad(\rho\uparrow1).                    \tag{SP1}
\]
The gradient is on the centered label plane. In particular
\(|\beta|=o(e)\) for balanced partitions at every fixed \(k\).
The theorem also gives locally uniform convergence of the
value profile and its weak derivatives:
\[
        \frac{1-V_\rho}{e}\longrightarrow
                      \sqrt{\frac2\pi}\,I
                   \quad\text{in }W_{\rm loc}^{1,\infty}.
                                                               \tag{SP2}
\]
In (SP2) the function converges uniformly and its almost
everywhere gradient converges in essential supremum on
each compact convex subset of the interior. No differentiability
of the finite-correlation value at every mass is assumed.

## Perimeter convergence at moving masses

Consider \(\rho_j\uparrow1\), masses \(p_j\to p\) with
\(p\) interior, and actual optimizers \(A^{(j)}\).
Reduce to dimension \(k-1\), using rotations and adding unused
coordinates when needed. Let \(T(p_j)\) be the translated
regular simplex of masses \(p_j\), in fixed model coordinates.
Its translations are bounded and converge to the translation
of \(T(p)\): the model centered prices are continuous in
the interior.
[Heilman–Tarter, Theorem 1.9](https://arxiv.org/pdf/2011.05583v2)
[The smooth model perimeter profile](GAUSSIAN_WEIGHTED_CENTROID_INEQUALITY.md#posterior-gradients-in-the-perimeter-metric)

Write \(N_\rho(A)=1-S_\rho(A)\) and set
\(a_j^2=(1+\rho_j)/2\). For pair midpoint \(w\) and rescaled
displacement \(z\), the joint density is
\(\phi_{a_j}(w)dw\,d\gamma(z)\). For each label let
\[
 d\mu_{j,i}(w,z)=
 \frac{|1_{A_i^{(j)}}(w+e_j z/2)-1_{A_i^{(j)}}(w-e_j z/2)|}{e_j}
                      \phi_{a_j}(w)dw\,d\gamma(z).
\]
Their total mass summed over labels is
\(2N_{\rho_j}(A^{(j)})/e_j\).
The compactness and lower-bound proof of (P2)–(P3) applies
in every fixed dimension. In particular, a bound on these
masses yields local \(BV\) compactness after mollification at
scale \(e_j\). Its Gaussian label-distance limit \(A\)
has masses \(p\). Each measure liminf dominates
\[
 d\mu_i(w,z)=\phi_1(w)|z\cdot n_i(w)|\,
                    d|D1_{A_i}|(w)d\gamma(z),
\]
where \(n_i\) is the measure-theoretic unit normal.
[Midpoint compactness and measure lower bound, (P2)–(P3)](THREE_CELL_NEAR_ONE.md#compactness-and-the-sharp-perimeter-lower-bound)

The flat facets of \(T(p_j)\) yield
\(N_{\rho_j}(T(p_j))/e_j\to\sqrt{2/\pi}\,I(p)\).
Indeed the direction-\(z\) difference across each facet is
an integral over an \(e_j|z|\)-thick translated strip.
The centered Gaussian weight in the strip converges to the
facet weight. The resulting directional integral is
bounded by \(C|z|e^{C e_j^2|z|^2}\) for small \(e_j\),
uniformly over the bounded translations of \(T(p_j)\).
Intersections of distinct facet strips have vanishing
contribution after division by \(e_j\) for each fixed \(z\).
Equivalently one can apply directional slicing to the fixed
fan before translating its Gaussian weight by the convergent
model price vector. This proves the limit by domination.

Maximality gives the same upper bound on the limit for
\(A^{(j)}\). The measure lower bound and the Gaussian
multi-bubble theorem then force \(A=T(p)\) up to an orthogonal
map and null sets. They also force convergence of the total
masses to the sum of the masses of the \(\mu_i\). Every
vague measure limit dominates \(\mu_i\), so equality of the
sums gives narrow convergence of each measure along this
subsequence. This proves the uniform convergence of the
values in (SP2) by compactness of the mass parameters.
[Milman–Neeman, Gaussian multi-bubble and equality, Theorems 1.1–1.2](https://arxiv.org/html/1805.10961v3)

## Passing the first variations to the perimeter profile

Let \(W\) be a fixed compactly supported smooth vector field
and flow both observations by it. The first variation of
\(N_{\rho_j}/e_j\) is the integral against
\(\tfrac12\sum_i\mu_{j,i}\) of
\[
 \operatorname{div}W(x)+\operatorname{div}W(y)
  -\frac{w\cdot(W(x)+W(y))}{2a_j^2}
  -\frac{z\cdot(W(y)-W(x))}{e_j},
 \quad x=w-e_jz/2,\quad y=w+e_jz/2 .
\]
This follows by differentiating the pullback of the Gaussian
pair density and the two flow Jacobians. It requires no
pointwise boundary regularity for the moved cells.
On compact sets the integrand converges uniformly to
\(2\operatorname{div}W-w\cdot W-z^{\mathsf T}DW z\).
The tail passes to the limit as follows. Where the integrand
is nonzero an endpoint is in a fixed ball, so
\(|w|\le R+e_j|z|/2\) and the integrand is bounded by
\(C_W(1+|z|^2)\). Compare the two pair densities in endpoint
coordinates at noise lengths \(e_j\) and \(2e_j\).
On these pairs the ratio is at most \(C_R e^{-c|z|^2}\):
the displacement exponent contributes \(-3|z|^2/8\)
and the change in midpoint exponent is bounded by
\(C e_j^2|w|^2\).
For \(\rho'_j=1-4(1-\rho_j)>0\) the Hermite expansion
gives \(N_{\rho'_j}(A^{(j)})\le4N_{\rho_j}(A^{(j)})\),
using \(1-(\rho'_j)^m\le4(1-\rho_j^m)\) for \(m\ge1\).
Thus truncation in \(z\) costs arbitrarily little uniformly
in \(j\), and for bounded \(z\) the active midpoint set
is bounded.

At the limit,
\(\mathbb E[|z\cdot n|zz^{\mathsf T}]
=\sqrt{2/\pi}(I+nn^{\mathsf T})\) for a unit \(n\).
The limit of the first variation is therefore
\(\sqrt{2/\pi}\,P'_\gamma(T(p))[W]\), where the integrand
for a facet is its tangential divergence of \(W\) minus
\(w\cdot W\). This proves
\[
 \frac{N'_{\rho_j}(A^{(j)})[W]}{e_j}
        \longrightarrow \sqrt{2/\pi}\,P'_\gamma(T(p))[W].
                                                               \tag{SP3}
\]
The varying masses do not affect the argument. In particular
the conclusion stays in the fixed active dimension.
[The balanced planar first-variation argument](THREE_CELL_NEAR_ONE.md#score-prices-are-smaller-than-the-noise-length)

Let \(m'_A[W]\) denote the vector of volume derivatives.
It has the bulk expression with integrands
\(1_{A_i}(\operatorname{div}W-x\cdot W)\).
Thus \(m'_{A^{(j)}}[W]\to m'_{T(p)}[W]\) after alignment.
For each actual maximizer the adjusted-score winner rule
makes the linear functional of the adjusted scores
stationary under every such unconstrained flow.
Differentiating gives
\[
                 N'_{\rho_j}(A^{(j)})[W]
                       =-2\beta_j\cdot m'_{A^{(j)}}[W].
\]
The derivative is obtained from the two first-order
terms, one from each endpoint, in the smooth pair
density pullback.
At the model, \(P'_\gamma(T(p))[W]
=\nabla I(p)\cdot m'_{T(p)}[W]\).
Indeed the model perimeter equals the profile at its own
masses, and the perimeter of the flowed cluster is at
least the profile at its flowed masses in both directions
of time. Both functions are differentiable, so their
derivatives agree at zero.

Take \(k-1\) fixed test fields on distinct open model
facets whose incident label pairs form a tree. The volume
derivatives at \(T(p)\) span the centered label plane,
since the fields can have positive normal flux on just
their own facet. Their volume columns at \(A^{(j)}\)
converge to those at \(T(p)\). Inverting the column
matrix in (SP3) proves (SP1) along the subsequence.
Every sequence in a compact mass range has such a
subsequence, which proves uniformity.

## Differentiation of the optimal value

For fixed \(\rho\) the mass value
\(J_\rho=(1-V_\rho)/e\) is Lipschitz on the mass simplex:
one can change the masses from \(p\) to \(q\) by moving
at most \(\|p-q\|_1/2\) Gaussian measure between labels.
For the two correlated endpoints the union bound changes
agreement by at most \(\|p-q\|_1\). Apply this in both
directions. At a differentiability point \(p\) in a
compact interior mass range and for \(\rho\) sufficiently
close to one, the volume columns for any optimizer span
the centered plane. Otherwise an exceptional sequence
would contradict the test fields at its limiting model.
Use any such flow \(W\) as a competitor at its flowed
masses. Differentiating the upper bound on the minimum
disagreement at those masses, in both directions of time,
gives
\[
                    \nabla J_\rho(p)=-2\beta/e .
\]
By Rademacher's theorem this holds almost everywhere.
The uniform estimate (SP1) then gives the asserted
essential supremum convergence of the gradients, completing
(SP2). For example the error
\(J_\rho-\sqrt{2/\pi}I\) has Lipschitz seminorm tending
to zero on interior compact convex mass sets.
