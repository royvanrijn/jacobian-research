# Gaussian perimeter comparisons for nonlinear score flow

The Gaussian multi-bubble profile gives a sharp Ornstein–Uhlenbeck
entropy comparison for every number of labels. Its equality cases
are regular simplicial partitions. The entropy in this statement is
the Gaussian perimeter profile evaluated at posterior masses.
The squared-posterior objective in the Standard Simplex Conjecture
requires a different comparison. A second application of the
perimeter theorem below gives an explicit upper bound for that
objective and identifies where the bound loses sharpness.

## A graph formula for the functional perimeter inequality

Let \(H=\mathbf1^\perp\subset\mathbb R^k\), and let \(I(p)\)
be the Gaussian multi-bubble perimeter profile at an interior
mass vector \(p\) in the probability simplex. Write \(A_{ij}(p)\)
for the Gaussian area of the \(ij\) facet of its regular
simplicial cluster on \(H\), and write \(L_p\) for the graph
Laplacian with weights \(A_{ij}(p)\). All these weights are
positive. Each interface is counted once in \(I(p)\).

The model cluster is the winner of \(z_i-b_i(p)\), where
\(z\) is standard Gaussian on \(H\) and \(b(p)\in H\) is the
unique vector of prices giving masses \(p\). These prices satisfy
\[
             b(p)=\sqrt2\,\nabla I(p),
     \qquad D^2 I(p)=-L_p^{-1}\ \text{on }H.                    \tag{PF1}
\]
Indeed the affine assignment dual
\(\inf_{b\in H}[\mathbb E\max_i(Z_i-b_i)+b\cdot p]\)
is coercive for interior \(p\) and its minimizer gives the
prescribed masses. Gaussian divergence across the model facets
makes its value \(\sqrt2 I(p)\). On each facet the signed
difference of the two \(z\)-scores has gradient of norm
\(\sqrt2\), so \(D_b p=-L_p/\sqrt2\) by the moving-wall
formula. The strictly positive weights make this derivative
invertible on \(H\). The envelope derivative and inversion
give (PF1). This is the standard Gaussian model-profile
Hessian formula.
[Milman–Neeman, Theorems 1.1–1.2 and equations (1.2)–(1.4)](https://arxiv.org/html/1805.10961v3)

For a smooth map \(p(x)\) into the interior simplex, consider
on \(\mathbb R^d\times H\) the label maximizing \(z_i-b_i(p(x))\).
The cluster has masses \(\bar p=\mathbb E p(X)\). Write
\(q(x)=\nabla I(p(x))\). The \(ij\) facet lies over the graph
\((z_i-z_j)/\sqrt2=q_i(x)-q_j(x)\). Its area contribution,
integrated in \(x\), is exactly
\(\mathbb E[A_{ij}(p)\sqrt{1+|\nabla_x(q_i-q_j)|^2}]\).
The graph factor is independent of the remaining label Gaussian
coordinates. Triple ties have codimension at least two since
their label normals are independent. The area identity can
therefore be proved on bounded sets and then extended by
exhaustion, with infinite perimeter allowed. Thus the
multi-bubble theorem yields
\[
 I(\bar p)\le \mathcal P(p):=
 \mathbb E\sum_{i<j}A_{ij}(p)
                  \sqrt{1+|\nabla_x(q_i-q_j)|^2}.              \tag{PF2}
\]
For \(k=3\), this is the functional inequality of
Milman–Neeman announced in Milman's 2023 slides. To identify
the coefficients, put \(a=A_{12}\), \(b=A_{13}\), \(c=A_{23}\).
Solving \(L_p\,dq=-dp\) gives
\(a\,d(q_1-q_2)=-(ac\,dp_1-ab\,dp_2)/(ab+ac+bc)\),
which gives their two weights under the square root.
For other pairs permute the three labels.
[Milman, 2023 Toulouse slides, pp. 46–49](https://indico.math.cnrs.fr/event/8783/contributions/8718/attachments/4115/6002/Milman-Toulouse2023.pdf#page=46)

By (PF1) the weighted sum of the squared graph slopes is
\[
 \sum_{i<j} A_{ij}(p)|\nabla_x(q_i-q_j)|^2
   =\mathcal D(p):=
      \sum_{\ell=1}^d (\partial_\ell p)^{\mathsf T}
                            [-D^2 I(p)]\,\partial_\ell p .
\]
Cauchy–Schwarz on the product of the sample space and the
weighted pair index gives
\[
 I(\bar p)^2\le \mathcal P(p)^2
               \le (\mathbb E I(p))\,\mathbb E[I(p)+\mathcal D(p)] .
                                                               \tag{PF3}
\]

## Sharp comparison along the indicator flow

Let \(A=(A_i)_{i=1}^k\) be an ordinary Gaussian partition
with positive masses \(\pi_i\), and put
\(p_t=(T_{e^{-t}}1_{A_i})_i\) and \(J(t)=\mathbb E I(p_t)\).
Then for every \(t>0\),
\[
                   J(t)\ \ge\ \sqrt{1-e^{-2t}}\,I(\pi).
                                                               \tag{PF4}
\]
Equality at any one positive finite \(t\) is equivalent to
\(A\) being a regular simplicial cluster of masses \(\pi\),
up to null sets. At equal masses this is the centered simplex
up to orthogonal maps and labels. If \(d<k-1\), equality is
impossible.

For the proof, first replace \(p_t\) by
\(p_t^\epsilon=(1-\epsilon)p_t+\epsilon\pi\). This remains
in a compact interior mass set, still obeys the
Ornstein–Uhlenbeck equation, and has the same mean.
Gaussian integration by parts gives
\(J_\epsilon'(t)=\mathbb E\mathcal D(p_t^\epsilon)\) for
\(J_\epsilon(t)=\mathbb E I(p_t^\epsilon)\). Thus (PF3)
gives \(J_\epsilon(t)(J_\epsilon(t)+J_\epsilon'(t))
\ge\mathcal P(p_t^\epsilon)^2\). Integrating, with the
nonnegative initial term discarded, gives
\[
 e^{2t}J_\epsilon(t)^2\ge
                2\int_0^t e^{2s}\mathcal P(p_s^\epsilon)^2\,ds .
\]
Every \(p_s(x)\) is interior for \(s>0\) by positivity of
the Gaussian kernel. The profiles are bounded by \(\mu_k/\sqrt2\)
with \(\mu_k=\mathbb E\max_i Z_i'\) for independent standard
normals \(Z_i'\), as follows from the zero-price trial in the
affine support. Dominated convergence on the left and Fatou
on the nonnegative area integrands and then in \(s\) give
\[
 e^{2t}J(t)^2\ge 2\int_0^t e^{2s}\mathcal P(p_s)^2\,ds
                         \ge (e^{2t}-1) I(\pi)^2 .            \tag{PF5}
\]
This argument does not require a uniform bound on the
profile Hessian near the boundary of the mass simplex.

If equality holds in (PF4), then \(\mathcal P(p_s)=I(\pi)\)
for almost every \(s\in(0,t)\). Fix such an \(s\).
The lifted winner is a regular simplicial cluster by the
multi-bubble equality theorem. The entire top \(ij\) fiber
over each \(x\) contains a nonempty relatively open part
of its equality hyperplane. Its supporting affine wall in the
lifted cluster therefore has label normal proportional to
\(e_i-e_j\). The proportionality factor is positive, since
the \(i\) label wins for larger \(z_i-z_j\) near the facet.
Consistency of the three normals around any label triangle
makes these factors equal across all pairs. After a common
rescaling of the sites, their label components can be taken
to be \(P e_i\), where \(P\) projects onto \(H\). The pair
wall equations then make the price vector \(b(p_s(x))=Bx+b_0\)
affine with values in \(H\). Regularity of the lifted sites
makes the squared distances between rows of \(B\) constant
across pairs, since their label component distances are
already all equal. Hence \(BB^{\mathsf T}=\alpha^2 P\) for
some \(\alpha\ge0\).

An \(\alpha=0\) would make \(p_s\) constant, impossible
for a positive-mass ordinary partition by injectivity of
\(T_{e^{-s}}\). Thus \(-B/\alpha\) is a coisometry onto \(H\).
The conditional winner formula now expresses \(p_s\) as
\(T_r1_{A^*}\) for a regular simplicial cluster \(A^*\)
with shifted prices \(b_0/\sqrt{1+\alpha^2}\) in that active
space, where \(r=\alpha/\sqrt{1+\alpha^2}\). Its mass vector
is \(\pi\). If \(r\ne e^{-s}\), use the semigroup law and
injectivity of the operator with the larger parameter.
One of \(1_A,1_{A^*}\) would be a strictly smoothed ordinary
indicator, which is strictly between zero and one everywhere
for each positive-mass label of complement positive mass.
Thus \(r=e^{-s}\) and \(A=A^*\) almost everywhere. Equal
masses force centered prices by symmetry and uniqueness.
Conversely, for a regular simplicial \(A\), the lift at
time \(s\) is again a regular simplicial cluster with masses
\(\pi\). Each of its pair walls has the same graph factor
\(1/\sqrt{1-e^{-2s}}\). This gives equality in (PF4)
directly from the area formula.

Conditioning gives a pointwise version. For \(0<r,s<1\)
and every \(x\), one has
\[
 T_r[I(T_s1_A)](x)\ \ge\
       \sqrt{\frac{1-s^2}{1-r^2s^2}}\,I(T_{rs}1_A(x)).
                                                               \tag{PF6}
\]
In fact, write the conditional sample inside the left
side as \(rsx+s\sqrt{1-r^2}\,Z+\sqrt{1-s^2}\,W\)
with independent standard Gaussians \(Z,W\). Apply (PF4)
to the ordinary partition obtained by translating by
\(-rsx\) and dividing by \(\sqrt{1-r^2s^2}\), with
inner multiplier \(s\sqrt{1-r^2}/\sqrt{1-r^2s^2}\).
The translated partition has positive Gaussian mass in
each label. Equality at one \(x\) holds exactly for
regular simplicial partitions, and then holds for all \(x\).

The compared functional \(J(t)\) is a perimeter-profile
average. The desired noise stability is
\(S_{e^{-2t}}(A)=\mathbb E\sum_i p_{t,i}^2\).
A sharp comparison of the latter does not follow from (PF4),
which controls a different function on the mass simplex.

## Perimeter and the two derivatives

Fix \(k\ge2\), let \(H=\mathbf1^\perp\subset\mathbb R^k\),
and let \(f:\mathbb R^d\to H\) be bounded and measurable.
For \(t\ge0\) put \(F_t=T_{e^{-t}}f\). Use independent
standard Gaussians \(X\in\mathbb R^d\) and \(Z\in H\).
For \(\lambda\ge0\), let \(C(t,\lambda)\) be the optimal
balanced linear assignment value of \(F_t(X)+\sqrt\lambda Z\).
In other words, maximize the expected score of the chosen label
over all measurable labels of \((X,Z)\) of masses \(1/k\).

For \(\lambda>0\), minimize over \(b\in H\) the convex
function
\[
 \mathbb E\max_i\bigl(F_{t,i}(X)+\sqrt\lambda Z_i-b_i\bigr).
\]
It is coercive since choosing any one label inside the
maximum bounds it below by \(-\min_i b_i-\|f\|_\infty\).
Its derivative is minus the vector of winning label masses,
restricted to \(H\), so the minimizing prices give a balanced
winner. This is an optimal assignment since the price contribution
is constant on balanced assignments. The scores have a
positive smooth Gaussian-mixture density on \(H\). In score
coordinates the Hessian in the prices is the Laplacian whose
\(ij\) weight is the integral of this density over the open
top \(ij\) facet, divided by \(\sqrt2\). Indeed moving the
difference of the two prices moves their equality plane at
speed \(1/\sqrt2\), and the triple ties have zero facet measure.
The weights are strictly positive, giving uniqueness in \(H\).
Differentiation of the mixture is justified on compact
parameter intervals with \(t,\lambda>0\), since \(F_t\)
and its spatial derivatives are bounded and \(LF_t\) has at
most linear growth, where \(L=\Delta-x\cdot\nabla\). The
implicit function theorem and the envelope derivative apply.

Write \(\Sigma_{ij}\) for a pair facet of this assignment in
\((x,z)\) coordinates and \(d\nu_{ij}\) for standard Gaussian
surface measure there. With
\[
 a_{ij}=\sqrt{|\nabla_x(F_{t,i}-F_{t,j})|^2+2\lambda},
\]
the two derivatives are
\[
 C_\lambda=\sum_{i<j}\int_{\Sigma_{ij}}\frac{d\nu_{ij}}{a_{ij}},
 \qquad
 -C_t=\sum_{i<j}\int_{\Sigma_{ij}}
    \frac{|\nabla_x(F_{t,i}-F_{t,j})|^2}{a_{ij}}\,d\nu_{ij}.
                                                               \tag{SF1}
\]
For the second identity, apply Gaussian divergence to
\(\sum_i\int_{\text{cell }i} LF_{t,i}\). The outward normal
of cell \(i\) on \(\Sigma_{ij}\) is the negative normalized
gradient of its contrast with \(j\). The two fluxes add
to the indicated negative quantity. For the first identity
differentiate the score in \(\lambda\) and integrate
\(Z_i/(2\sqrt\lambda)\) over cell \(i\) by Gaussian divergence
in \(H\). Each endpoint label of a facet contributes
\(1/(2a_{ij})\). The equality hypersurfaces are graphs of a
bounded function of \(x\) in the \(e_i-e_j\) direction in
\(H\), with bounded gradient for fixed \(t,\lambda>0\).
This gives finite Gaussian surface integrals and justifies
exhaustion in these divergence arguments. Triple ties have
codimension at least two by independence of the two label
normals in \(H\).

Let \(P_*(k)\) denote the total interface area of the
balanced Gaussian simplex, counting each facet once. The
balanced winning assignment above is in dimension \(d+k-1\).
Cauchy–Schwarz and the Gaussian multi-bubble theorem give
\[
 C_\lambda(-C_t+2\lambda C_\lambda)
    =\left(\sum_{i<j}\int a_{ij}^{-1}\,d\nu_{ij}\right)
      \left(\sum_{i<j}\int a_{ij}\,d\nu_{ij}\right)
    \ge P_*(k)^2.                                             \tag{SF2}
\]
For \(\mu_k=\mathbb E\max_{i\le k} Z_i'\), with \(Z_i'\)
independent standard normals, one has
\(\sqrt2P_*(k)=\mu_k\). This follows from the first
identity in (SF1) for the pure label score \(\sqrt\lambda Z\),
whose support is \(\sqrt\lambda\mu_k\) and whose winner
is the balanced simplex. Projection from the independent
coordinates \(Z_i'\) to \(H\) subtracts only a mean-zero
common score.
[Milman–Neeman, Theorems 1.1–1.2](https://arxiv.org/html/1805.10961v3)

## Integration from indicator initial data

For any measurable nonnegative control \(v\), (SF2) gives
\[
 C_t-(2\lambda+v)C_\lambda\le -2P_*(k)\sqrt v,
\]
using \(P_*^2/C_\lambda+vC_\lambda\ge2P_*\sqrt v\).
Along a path \(\lambda'(u)=-2\lambda(u)-v(u)\) ending
at \(\lambda(t)=\lambda\), one has
\(\lambda(0)=e^{2t}\lambda+\int_0^t e^{2u}v(u)\,du\).
For a fixed \(A=\int_0^t e^{2u}v(u)\,du>0\), the integral
of \(\sqrt v\) is maximized by \(v(u)\) proportional to
\(e^{-4u}\). Weighted Cauchy–Schwarz gives the maximum
\(\sqrt{A(1-e^{-2t})/2}\). The path stays strictly positive
until its endpoint, even when the terminal \(\lambda\)
is zero. Integrating gives
\[
 C(t,\lambda)\le
  \inf_{A\ge0}\bigl[C(0,e^{2t}\lambda+A)
                         -\mu_k\sqrt{A(1-e^{-2t})}\bigr].
                                                               \tag{SF3}
\]
The support is continuous at \(t=0\) by strong continuity
of the Ornstein–Uhlenbeck semigroup in \(L^1\), and at
\(\lambda=0\) since adding the label scores changes the
support by at most \(\sqrt\lambda\,\mathbb E\max_i|Z_i|\).
These facts justify the endpoints of the integration,
and the case \(A=0\) follows by continuity.

For an ordinary balanced partition \(B\), take
\(f_i=1_{B_i}-1/k\). Its initial score distribution is
invariant under label permutations, so zero prices suffice
at \(t=0\). Conditional on the original label, the initial
value is
\[
 C(0,a^2)=
 \mathbb E\max_{i\le k}\bigl(aZ_i'+1_{\{i=1\}}\bigr)-1/k.
\]
At \(t=-\log\rho\), the assignment given by \(B\) on the
\(X\) coordinate is balanced and attains the linear value
\(S_\rho(B)-1/k\) in the definition of \(C(t,0)\).
Consequently for every \(0<\rho<1\),
\[
 S_\rho(B)\le U_k(\rho):=
 \inf_{a\ge0}\left[
    \mathbb E\max_{i\le k}\bigl(aZ_i'+1_{\{i=1\}}\bigr)
                  -a\mu_k\sqrt{1-\rho^2}\right].              \tag{SF4}
\]
The bound depends only on \(k,\rho\). It is strictly below
one. Indeed the expected maximum has right derivative zero
at \(a=0\), by dominated convergence using its almost sure
unique winner there, while the subtracted term has strictly
negative derivative. The infimum is attained at a positive
finite \(a\), since the expression grows linearly as
\(a\to\infty\).

## Strictness and the loss at the first moment

For every balanced ordinary \(B\), the comparison to \(U_k\)
is strict. To see this, at each interior point of the
optimizing characteristic, suppose its winning cluster had
perimeter \(P_*(k)\). The multi-bubble equality theorem
would make it an affine simplex cluster up to null sets.
Yet for each fixed \(x\), the top \(ij\) facet is a nonempty
relatively open subset of the label hyperplane
\[
         z_i-z_j=\frac{b_i-b_j-(F_{t,i}-F_{t,j})(x)}
                              {\sqrt\lambda}.
\]
Every such facet is regular, so equality of the cells up to
null sets forces this facet into the affine wall separating
the same two labels. An affine wall containing all these
fibers would have label normal parallel to \(e_i-e_j\) and
would make \(F_{t,i}-F_{t,j}\) affine in \(x\). A bounded
affine function is constant. This is impossible because
positive Ornstein–Uhlenbeck multipliers are injective and
the original distinct label indicators have nonconstant
contrast. Hence the perimeter inequality in (SF2) is
strict throughout the interior of the characteristic, and
its integral is strict. In particular \(U_k(\rho)\)
strictly exceeds the regular simplex value at every positive
correlation. Combined with fixed-dimensional attainment,
it strictly exceeds the unrestricted optimum as well.
[Existence and dimension reduction](GAUSSIAN_SPECTRAL_DIMENSION_REDUCTION.md#optimization-and-volume-repair)

The bound loses at first order as \(\rho\downarrow0\).
Gaussian integration, or the heat equation in the means
of the Gaussian maximum, gives
\[
 \mathbb E\max_i(a Z_i'+1_{\{i=1\}})
       =a\mu_k+\frac1k+\frac{\mu_k}{2ka}+O_k(a^{-2}).
\]
For the coefficient, the Hessian of
\(b\mapsto\mathbb E\max_i(Z_i'+b_i)\) at zero has trace
\(\mu_k\) by differentiation in the Gaussian variance,
and equal diagonal entries by permutation symmetry.
The minimizer in (SF4) tends to infinity as \(\rho\to0\).
At \(\rho=0\) the expression is strictly above \(1/k\)
at every finite \(a\) by the strict convexity of the
expected maximum in the zero-sum means for \(a>0\).
Its limit as \(a\to\infty\) equals \(1/k\). Optimizing
the two leading positive terms in (SF4) then gives
\[
             U_k(\rho)=\frac1k+\frac{\mu_k}{\sqrt{k}}\rho
                                      +o_k(\rho).
\]
Indeed the minimizing large \(a\) can be sandwiched
by replacing \(\mu_k/(2k)\) in the \(a^{-1}\) term by
arbitrarily close upper and lower constants. The regular
simplex has first Hermite mass \(\mu_k^2/(k-1)\), strictly
smaller than \(\mu_k/\sqrt{k}\). The strict comparison of
coefficients follows from
\(\mu_k\le\sqrt{(k-1)/k}\,\mathbb E\|Z\|<(k-1)/\sqrt{k}\)
for standard Gaussian \(Z\) in \(H\).
[Simplex first-moment value](QUANTITATIVE_SIMPLEX_FIRST_MOMENT_STABILITY.md#external-endpoint-and-application-audit)

Thus the score-flow inequality uses the actual balanced
Gaussian perimeter theorem throughout the path and supplies
an all-partition positive-correlation bound. A sharp flow
comparison would need to retain information discarded by
replacing the entire perimeter of the augmented assignment
by \(P_*(k)\).
