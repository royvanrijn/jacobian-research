# Scalar inward variations at a pressure-free tie

Suppose an ordinary balanced planar maximizer \(A\) ties the tripod
\(T\) at \(0<r<1\), the tripod is globally optimal just below \(r\),
and
\[
 \beta(A)=0,\qquad \dot S_r(A)=\dot S_r(T).
\]
The dot denotes the correlation derivative at fixed shape. These
hypotheses hold if the tie is invisible to the first joint covariance
and mass derivatives. Put
\(\kappa=\tfrac12(S''_r(T)-S''_r(A))\ge0\). A direct scalar
inward variation gives a strict second-order comparison with the
tripod branch for sufficiently small \(\kappa\). At zero curvature,
it is strict in every sufficiently small nonzero quota direction,
after relabeling \(A\).

Use the outward normal speed, mass flux \(C\), pair-score derivative
\(g=\partial_r(U_i-U_j)\), and unconstrained form \(Q_A=Q_0\)
of the [maximizer analysis](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).
Write \(G(f)=\int gf\,d\nu\) and \(\|f\|_a^2=\int a f^2d\nu\).
The linear functional \(G\) is continuous in this weighted norm by
(CT8).
[Anisotropic dilation identity](THREE_CELL_CORRELATION_CONTACT.md#anisotropic-dilation-and-mass-tangent-directions).

## A dilation test with its own volume drift

Let \(A\) be nonsimplex. The origin-centered dilation speed is
\(b(x)=x\cdot N(x)\). Its mass flux is
\[
 k=Cb,\qquad k_i=\int_{A_i}(2-|x|^2)\,d\gamma_2,\qquad
 H=Q_A(b)+\Lambda_r|k|^2,\qquad \alpha=-G(b).
\]
Here \(\Lambda_r\) is the exact coefficient of the local tripod
mass profile (LM1). By (CT7),
\[
 \alpha\ \ge\ \frac{\|b\|_a^2}{r(1+r)}\ >0.                 \tag{SI1}
\]
The strict inequality follows because a vanishing dilation speed
makes all cells conical about the origin, which forces the tripod
for a balanced maximizer.
[Conical rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#convex-and-conical-maximizers).
[Local mass profile](THREE_CELL_LOCAL_MASS_PROFILE.md#exact-expansion-over-measurable-competitors).

For fixed \(c\in\mathbb R\) and \(t\downarrow0\), dilate \(A\) by
\(1+ct\) and evaluate it at correlation \(r-t\). Its masses are
\(m(t)=m_*+ct k+O(t^2)\). Then
\[
 V_{r-t}^T(m(t))-S_{r-t}((1+ct)A)
       =\bigl[\kappa-2c\alpha+c^2 H\bigr]t^2+o(t^2).
                                                               \tag{SI2}
\]
Indeed, the first dilation derivative vanishes by zero prices,
the second derivative divided by two is \(-Q_A(b)\), and the
mixed derivative with \(r\) is \(2G(b)\). These formulas extend
from compact cutoff fields to dilation by polynomial-weighted
interface integrability as in (G6). In the bulk integrals the
pullback of the joint Gaussian density by simultaneous dilation
is smooth in the dilation and the correlation. Derivatives
through order three on a small parameter neighborhood are
bounded by an integrable polynomial times a Gaussian with
positive definite quadratic decay. This justifies the Taylor
remainder for the fixed measurable partition. Applying (LM1)
uniformly in correlation supplies the term \(c^2\Lambda_r|k|^2\)
with remainder \(o(t^2)\), including \(k=0\).

If \(H\le0\), choose \(c>0\) large. If \(H>0\), choose \(c=\alpha/H\).
Thus either \(H\le0\) or \(\alpha^2>\kappa H\) makes (SI2)
strictly negative. In particular this holds if \(\kappa=0\).
When \(k=0\), the test can be corrected to exact balance by compact
transports of size \(O(t^2)\). Their first stability derivative at
the original \(A\) is zero and remains \(O(t)\) during the
perturbation. The repair has cost \(O(t^3)\), so a strict sign in
(SI2) with \(k=0\) contradicts scalar simplex optimality to the
left. When \(k\ne0\) the strict sign separates the branches at
unequal quotas. If the balanced maximizing orbit immediately
to the left is unique, the excluded-annulus argument of (LM8)
gives separated global optimizers at some quota on the segment
from \(m_*\) to \(m(t)\) at each sufficiently small \(t\).
[Coexistence by compactness](THREE_CELL_LOCAL_MASS_PROFILE.md#consequences-at-a-first-balanced-tie).

A bound involving only the dilation norm is also available.
The constant Hermite mode and the first extra-coordinate
form give \(r|Cb|^2\le\|b\|_a^2\) by (IB1). Positivity of the
Mehler interaction gives \(Q_A(b)\le\|b\|_a^2\). Therefore
a sufficient condition for the strict sign in (SI2) is
\[
 \kappa <
 \frac{\|b\|_a^2}{r^2(1+r)^2(1+\Lambda_r/r)} .               \tag{SI3}
\]
By (G25) the numerator is bounded below over nonsimplex
balanced global maximizers on a compact positive correlation
band, and the denominator is bounded above there.
[Volume-flux bound](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux).
[Gap from conical geometry](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-uniform-gap-from-conical-geometry).

## The mass residual of a saturated response

Let \(\mathscr D\) be the compact smooth regular facet speeds.
The mass flux \(C\) is continuous on the weighted space by
(IB1) and surjective on \(\mathscr D\) by the independent
volume bumps. Thus compact balanced speeds are dense in the
weighted balanced space, by subtracting those bumps from
compact approximants. For a balanced speed \(f\), set
\[
 J(f)=-2G(f)-Q_A(f).
\]
Scalar comparison at \(r-t\) gives \(J(f)\le\kappa\) for compact
balanced speeds, and hence on the weighted balanced space.
This is the same test as (G10), written for an inward
correlation change.

Suppose the bound is attained by a balanced response \(u\) in the
weighted completion. A bounded weighted maximizing sequence also
suffices for attainment, since \(Q_A\) is nonnegative on the
balanced subspace and hence weakly lower semicontinuous there.
First variation of \(J\) on the balanced space gives a unique
centered label vector \(\pi\) such that
\[
             -G(f)-Q_A(u,f)=\pi\cdot Cf
                  \qquad(f\text{ in the weighted space}).       \tag{SI4}
\]
In fact the left side vanishes on \(\ker C\), and \(C\) is
onto the centered label space.

If \(\pi\ne0\), take any compact regular \(w\) with
\(\pi\cdot Cw>0\). For \(\epsilon>0\) small,
\[
 J(u+\epsilon w)-\kappa-\Lambda_r|\epsilon Cw|^2
   =2\epsilon\,\pi\cdot Cw
     -\epsilon^2[Q_A(w)+\Lambda_r|Cw|^2]>0.              \tag{SI5}
\]
Approximate \(u\) in weighted norm by compact balanced
speeds and keep one approximant sufficiently close for the
strict sign in (SI5). Move at its sum with \(\epsilon w\)
with normal amplitude \(t\). Exact compact repairs set the
quotas to \(m_*+t\epsilon Cw\), at cost \(O(t^3)\) for the
same zero-price reason as above. Applying (LM1) shows a strict
gain over the simplex branch of order \(t^2\). Thus the
conclusion concerns actual ordinary partitions, even when the
saturated response itself is known only in weighted norm.

At \(\kappa=0\) one can take \(u=0\). Equation (SI1) then
gives \(\pi\cdot k=\alpha>0\), and in particular \(k\ne0\).
For every small nonzero fixed direction \(h\) in the mass
plane, cyclically relabel the tied shape to get
\(2\pi\cdot h\ge|\pi|\,|h|\). Take the compact speed \(w\)
from a fixed linear right inverse of \(C\) with \(Cw=h\).
The right side of (SI5) at \(\epsilon=1\) is now bounded below
by \(|\pi|\,|h|-O(|h|^2)>0\), proving branch separation
at \(m_*+th\) in every such direction.

With a saturated response \(u\), the residual on translations
also has a geometric formula. For \(s_v=v\cdot N\), (CT6) and
the translation row (G2) give
\[
 \pi\cdot Cs_v
    =\frac{1-r}{r}\left\langle
            u+\frac{b}{1-r^2},\,s_v\right\rangle_a .            \tag{SI6}
\]
Because the translation mass map is invertible, \(\pi=0\)
is equivalent to \(P_{\rm tr}u=-P_{\rm tr}b/(1-r^2)\),
where \(P_{\rm tr}\) is the weighted translation projection.
Together, (SI2) and (SI5) give explicit inward tests for
the part of a tied shape's response that is visible at second
order in the scalar correlation and the quotas. At zero
curvature, (SI4)–(SI6) identify the same translation response
detected by the cubic agreement moment in (AG4). The dilation
criterion and attained residual test also apply at positive
deficit curvature.
[Translation test from the cubic agreement moment](THREE_CELL_AGREEMENT_MOMENT_RIGIDITY.md#a-quota-comparison-below-the-tie).

## Translation is strictly suboptimal for some quota direction

At any balanced planar maximizer, some nonzero quota direction
admits a strictly larger second-order mass gain than the
translation calculation (VV9). This comparison holds for nearby
ordinary partitions at the same correlation. It uses no attainment
assumption for the optimal unbalanced shape response of \(A\).

Indeed, let \(V=\{s_v:v\in\mathbb R^2\}\) in the weighted space.
If all translations were orthogonal to the balanced space
\(\ker C\), then \(V=\operatorname{ran} C^*\), since both spaces
have dimension two. In our orientation \(C^*w\) on \(e=(i,j)\)
equals \((w_i-w_j)/a_e\), where \(w\) is a centered label vector
and \(N_e\) points out of \(A_i\). There would be an invertible
linear map \(L:\mathbb R^2\to\mathbf1^\perp\) with \(s_v=C^*Lv\).
Consequently
\[
        a_e N_e=L^{\mathsf T}(e_i-e_j)
             \qquad\text{on each interface }e=(i,j).
                                                               \tag{SI7}
\]
Set \(p_i=-L^{\mathsf T}(e_i-\mathbf1/3)\) and define the
piecewise constant vector field \(F=\sum_i p_i1_{A_i}\).
Its distributional Jacobian is the locally finite matrix
measure
\[
 D F=\sum_{e=(i,j)}(p_j-p_i)\otimes N_e\,
                    \mathcal H^1|_{\Sigma_e}
     =\sum_e a_eN_e\otimes N_e\,\mathcal H^1|_{\Sigma_e}.
                                                               \tag{SI8}
\]
Here the outward normal contributes \(-N_e\) to \(D1_{A_i}\)
and \(+N_e\) to \(D1_{A_j}\). In particular \(DF\) is symmetric
positive semidefinite. Mollify \(F\) with a nonnegative compact
smooth kernel. The mollified field \(F_\varepsilon\) is smooth
with symmetric Jacobian and has convex potential
\(\Phi_\varepsilon(x)=\int_0^1 F_\varepsilon(tx)\cdot x\,dt\),
normalized to zero at the origin. The potentials are uniformly
Lipschitz because their gradients are bounded by \(\max_i|p_i|\).
A locally uniform subsequential limit is convex and has
distributional gradient \(F\).

For this limit \(\Phi\), the points of differentiability with
gradient \(p_i\) lie in the convex minimizer set of
\(\Phi(x)-p_i\cdot x\). This set has nonempty interior since
it contains almost every point of \(A_i\). Its interior has
gradient \(p_i\), and its boundary is Lebesgue-null. The three
vectors \(p_i\) are distinct, so each minimizer set is a convex
representative of \(A_i\). Balanced convex rigidity now gives
the regular tripod. In that case (SI7) would make the weight
constant on each ray, whereas (LM6) gives strict variation.
This contradiction proves \(V\not\subset(\ker C)^\perp\).
[Convex rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#convex-and-conical-maximizers).
[Tripod mass correction](THREE_CELL_LOCAL_MASS_PROFILE.md#the-correction-to-translation-curvature).

Choose \(v\) and a balanced speed \(f\) with
\(\langle s_v,f\rangle_a>0\). A compact regular balanced speed
with the same strict sign is obtained by first taking a compact
regular weighted approximant and then removing its small flux
with fixed independent bumps, using continuity of \(C\) from
(IB1). For \(h=Cs_v\) and sufficiently small \(\epsilon>0\),
the compact perturbation of translation \(s_v+\epsilon f\)
has the desired flux and satisfies
\[
 -Q_A(s_v+\epsilon f)
   =-Q_A(s_v)+2\epsilon\frac{1-r}{r}\langle s_v,f\rangle_a
                    -\epsilon^2Q_A(f)>-Q_A(s_v).
                                                               \tag{SI9}
\]
The global translation and the compact regular motion realizing
\(f\) give a smooth path at \(A\). Exact mass repair by fixed compact
fields changes the path by \(O(t^2)\) at amplitude \(t\), hence
does not affect the quadratic coefficient of the priced objective.
Writing \(V_r\) for the unrestricted mass profile, this proves
\[
 \liminf_{t\to0}
 \frac{V_r(m_*+th)-V_r(m_*)-2t\beta(A)\cdot h}{t^2}
     > (1-r)M_A(v,v).
                                                               \tag{SI10}
\]
