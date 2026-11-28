# Local state obstructions for Gaussian Bellman proofs

Fix \(k\ge 3\) and \(0<\rho<1\). A supermartingale bound for Gaussian
noise stability that uses the two conditional label probability vectors
has a positive gap at the regular simplex, even if it also uses time
and both observation points. A positive gap also holds when the state
includes finitely many heat-time derivatives. An estimate using
arbitrary finite spatial derivatives can distinguish the simplex from
its competitors. Even then, a bound that is locally Lipschitz in the
state at every simplex observation cannot be sharp.

The obstruction is a variation on one simplex facet. It improves the
conditional agreement while preserving global balance and both
conditional label vectors exactly. With spatial derivatives included,
the variation preserves the larger state to first order. The
nonconstancy of the conditional scores on the facet uses the two
opposite directions of its supporting line, as in
[Heilman–Mossel–Neeman, Lemma 4.1](https://arxiv.org/pdf/1403.0885).
Their Theorem 2.6 gives a same-mass improvement of an unequally shifted
flat partition. Here the variation respects the original Gaussian
volumes as well as the observations. The
[shared-reveal comparison](GAUSSIAN_POSTERIOR_PARTITION_COMPATIBILITY.md)
uses a single posterior volume vector.

## Supermartingale bounds and conditional agreement

Work in the active space \(\mathbb R^{k-1}\), with standard Gaussian
measure \(\gamma\), and let \(T\) be the centered regular simplex
partition. Let \((B_t,C_t)_{0\le t\le1}\) start at the origin and be a
pair of standard Brownian motions with cross variation
\(d\langle B_i,C_j\rangle_t=\rho\delta_{ij}\,dt\). For any ordinary
partition \(A\), write \(P_s f(x)=\mathbb E f(x+\sqrt{s}\,Z)\) for
standard Gaussian \(Z\), and use
\[
 p_A(t,x)=\big(P_{1-t}1_{A_i}(x)\big)_{i=1}^k .
\]
For \(0<t<1\), let \(V_A(t,x,y)\) be the conditional probability of
agreement of the two labels of \(A\) given \(B_t=x,C_t=y\).
In particular, \(S_\rho(A)=\mathbb E V_A(t,B_t,C_t)\).

A bound in the stated class assigns a function \(\Psi\) to a local state
at time \(t\). The basic state consists of \(x,y,p_A(t,x),p_A(t,y)\).
The larger states below append derivatives of the two probability
vectors. Assume the following properties.

- For each ordinary balanced \(A\) and fixed \(0<t<1\), evaluation
  of \(\Psi\) at that state is continuous as a function of \(x,y\).
- For each such \(A\), evaluation along \((B_t,C_t)\) defines an
  integrable supermartingale on \(0\le t\le1\), whose terminal value
  is at least the indicator that the labels of \(A\) agree. At time
  \(1\) the process may be defined by an endpoint extension.

A pointwise Itô drift inequality with a valid supermartingale and
terminal extension is one sufficient way to meet these hypotheses.
They only require the supermartingale for a single balanced input
partition used at both observations.

For a fixed \(A\), the supermartingale property gives
\[
 \Psi\big(t,\operatorname{state}_A(t,x,y)\big)\ge V_A(t,x,y)
 \qquad(0<t<1,\ x,y\in\mathbb R^{k-1}).                         \tag{BL1}
\]
The inequality holds almost surely at \((B_t,C_t)\) by conditional
expectation. The law of \((B_t,C_t)\) has positive density, and both
sides are continuous in \(x,y\), so the inequality holds everywhere.
Continuity of \(V_A\) follows from continuity in \(L^1\) of the
translated, nondegenerate Gaussian joint density. In particular,
(BL1) may be applied to a different balanced partition for each
deterministic pair \(x,y\).

## A variation with three Gaussian volume constraints

Fix a time \(0<t<1\) and write \(s=1-t\). Let \(f_\gamma\) denote the
standard Gaussian density and \(f_x,f_y\) denote the densities with
means \(x,y\) and covariance \(sI\). Choose distinct simplex labels
\(i,j\). Let \(n\) be the unit normal to their shared facet pointing
into \(T_i\), and \(v\) a unit vector in the relative interior of the
facet. Choose \(x,y\) so that
\[
 (y-x)\cdot v<0,\qquad (y-\rho x)\cdot n\ne0.                    \tag{BL2}
\]
These conditions hold on a nonempty open set, since \(v\perp n\).

Move only the \(i,j\) interface, by a smooth compactly supported
normal graph on the relative interior of the facet. Use the sign
of the displacement so that a positive displacement expands \(T_i\).
A speed \(\eta\) changes the \(i\)-volume for any of the three
densities \(f_\gamma,f_x,f_y\) at the rate \(\int\eta f\,d\sigma\), where
\(\sigma\) is Euclidean facet area. The derivative of \(V_A(t,x,y)\)
at \(A=T\) is \(\int\eta g\,d\sigma\), where
\[
\begin{split}
 g(z)&=f_x(z)U_x(z)+f_y(z)U_y(z),\\
 U_x(z)&=\mathbb E(1_{T_i}-1_{T_j})
             \big(y+\rho(z-x)+\sqrt{s(1-\rho^2)}\,Z\big).
\end{split}                                                     \tag{BL3}
\]
Define \(U_y\) by interchanging \(x,y\). This follows by moving each
endpoint's label in the conditional agreement separately. The term
with changes at both endpoints is quadratic in the size of the
displacement, since the transferred slabs are compact, their volumes
are linear in the size, and the joint density is bounded there.
Smoothness of the densities and conditional scores gives the same
quadratic remainder in the terms with one moved endpoint.

We will use the stronger fact that \(g\), restricted to the ray
\(z=\ell v\) with \(\ell>0\), lies outside the span of \(f_\gamma\) and
the functions \(f_x\) and \(f_y\) multiplied by polynomials in \(\ell\)
of any fixed finite degree. On this line \(f_y/f_x=\exp(a\ell+b)\)
with \(a=(y-x)\cdot v/s<0\), and both \(f_x/f_\gamma\) and \(f_y/f_\gamma\)
decay as \(\exp(-c\ell^2+O(|\ell|))\) for \(c>0\). Each \(U\) is
real analytic on the line. This follows by extending its Gaussian
convolution to complex \(\ell\), using an integrable Gaussian
majorant on every compact complex disk.

Suppose \(g=b_0f_\gamma+f_x Q_x+f_y Q_y\) on the ray, with polynomials
\(Q_x,Q_y\). Analyticity gives the identity on the full line.
Divide by \(f_\gamma\) and take \(\ell\to+\infty\) to obtain \(b_0=0\).
In the opposite direction, the ray \(-v\) has a positive angular
distance from \(T_i\) and \(T_j\), since the other simplex scores
strictly exceed their scores there. The Gaussian tail bound therefore
gives \(U_x,U_y=O(e^{-c' \ell^2})\) for \(\ell\to-\infty\) with
some \(c'>0\). Dividing by \(f_y\) shows that \(Q_y\) tends to zero
as \(\ell\to-\infty\), so \(Q_y=0\). Dividing the resulting identity
by \(f_x\) then gives \(Q_x=0\) in the same way. But at \(+\infty\)
the other labels are asymptotically absent in the law defining \(U_x\),
and the threshold between \(i,j\) is the sign in direction \(n\).
Consequently
\[
 \lim_{\ell\to+\infty}\frac{g(\ell v)}{f_x(\ell v)}
 =2\Phi\!\left(\frac{(y-\rho x)\cdot n}{\sqrt{s(1-\rho^2)}}\right)-1
 \ne0,                                                         \tag{BL4}
\]
where \(\Phi\) is the one-dimensional standard Gaussian distribution
function. This proves the nonspan assertion. In particular
\(f_\gamma,f_x,f_y\) are linearly independent on the facet, by the same
ratios on the positive ray.

Choose a pair \((x_*,y_*)\) satisfying (BL2). There is a compactly
supported smooth speed annihilating the three volume derivatives
at \((x_*,y_*)\) and having strictly positive agreement derivative.
Indeed, a linear functional on the speed space which vanishes on
the common kernel of the three volume derivatives factors through
their vector of values. Applied to agreement, this would give the
excluded span identity, by the fundamental lemma for test functions
and continuity. Choose three additional compactly supported speeds
whose volume derivative matrix is invertible. Linear independence
of the density restrictions gives these speeds by the same
test-function argument.

Use the first speed with amplitude \(\epsilon\) and the three
additional speeds as corrections. The implicit function theorem
fixes the three \(i\)-volumes exactly at their values for \(T\),
with corrections depending smoothly on \(\epsilon,x,y\) near
\(0,x_*,y_*\). The mass integrals in the normal graph and the
densities are smooth in those parameters. The corrections vanish
at \(\epsilon=0\) and have zero \(\epsilon\)-derivative at
\((0,x_*,y_*)\). The \(j\)-volumes are fixed as well, because the
union of the two cells is unchanged. For some small fixed
\(\epsilon>0\) the corrected partition strictly improves at
\((x_*,y_*)\). By continuity there is an open neighborhood \(U\)
and a number \(\delta>0\) such that for every \((x,y)\in U\) the
resulting ordinary balanced partition \(A^{x,y}\) satisfies
\[
\begin{split}
 p_{A^{x,y}}(t,x)&=p_T(t,x),\qquad p_{A^{x,y}}(t,y)=p_T(t,y),\\
 V_{A^{x,y}}(t,x,y)&\ge V_T(t,x,y)+\delta .
\end{split}                                                     \tag{BL5}
\]
All graph displacements can be taken in disjointness-preserving
tubes away from the other facets, because their supports lie in
a compact subset of the relative interior.

## A positive gap with local probabilities

Apply (BL1) to \(A^{x,y}\) on \(U\) and to \(T\) everywhere. By
(BL5) the local state on \(U\) is the same for \(A^{x,y}\) and \(T\).
Take expectations along the simplex process and use its
supermartingale property between times \(0\) and \(t\). With
\(u=(1/k,\ldots,1/k)\), this gives
\[
 \Psi(0,0,0,u,u)
 \ge S_\rho(T)+\delta\,\Pr\big((B_t,C_t)\in U\big)>S_\rho(T).      \tag{BL6}
\]
The positive term depends on the fixed variation and \(\rho,k,t\),
independently of \(\Psi\). Explicit dependence on \(x,y\) is allowed
throughout.

The same argument applies if the state also includes finitely many
heat-time derivatives at both observations. The added constraints
integrate the label indicators against \(\Delta_x^m f_x\) or
\(\Delta_y^m f_y\), up to powers of \(1/2\), where \(m\) runs over
the retained positive orders. On the ray, these weights are
\(f_x\) or \(f_y\) times polynomials of distinct degrees \(2m\).
The weights, along with \(f_\gamma,f_x,f_y\), are linearly independent.
To see this, divide any relation first by \(f_\gamma\) at \(+\infty\),
then by \(f_x\) at \(+\infty\), and use the distinct degrees within
each group. The nonspan proof for \(g\) already allows these
polynomials. Apply the implicit function theorem with the
corresponding number of correction speeds. This gives (BL5) with
the additional coordinates preserved exactly, hence a positive
gap independent of the bound in this larger class too.

## Finite spatial derivatives and local Lipschitz bounds

For any fixed order \(m\), append to the local state all spatial
derivatives through order \(m\) at both observations. Suppose
\(\Psi(t,\operatorname{state}_T(t,x,y))\) admits a local Lipschitz
bound with respect to perturbations of these state coordinates,
holding \(t,x,y\) fixed, at every simplex observation. It will
suffice to assume this bound at one pair satisfying (BL2).

At such a pair, the weights giving the first variations of these
coordinates restrict to \(f_x\) and \(f_y\) times polynomials in
\(\ell\). The nonspan result gives a compact smooth speed
annihilating the first derivative of every state coordinate and
the global balance derivative, with strictly positive conditional
agreement derivative \(b\). Correct just the global balance
exactly by the implicit function theorem, using one speed with
nonzero \(f_\gamma\)-integral. The correction is \(O(\epsilon^2)\).
The resulting balanced partitions \(A_\epsilon\) satisfy
\[
\begin{split}
 \operatorname{state}_{A_\epsilon}(t,x,y)
       -\operatorname{state}_T(t,x,y)&=O(\epsilon^2),\\
 V_{A_\epsilon}(t,x,y)-V_T(t,x,y)&=b\epsilon+O(\epsilon^2).
\end{split}                                                     \tag{BL7}
\]
There is no need for the full set of derivative constraints to
have independent first variations. The only constraint corrected
exactly here is global balance.

If (BL1) held with equality at this simplex state, local Lipschitz
continuity and (BL7) would contradict (BL1) for \(A_\epsilon\)
at sufficiently small \(\epsilon>0\). Thus there is a strict gap at
this observation, and continuity provides a positive gap on an
open neighborhood. The simplex supermartingale then has a
strict gap at time \(0\), with its own initial derivatives.
Here the size of the gap can depend on \(\Psi\). More generally,
sharpness at the given state is incompatible with a local
modulus \(o(\sqrt{h})\) for increments of size \(h\). This follows
from the same comparison in (BL7).

## The first derivatives and a singular bound

Exact first derivative constraints have a different effect from
first-order vanishing in (BL7). Let \(T_i\) select the largest
score \(a_i\cdot z+b_i\) from distinct affine scores with null ties.
For any nondegenerate Gaussian measure \(\mu\), an assignment
with the same \(\mu\)-masses and the same vector first moments as
\(T\) has the same integrated winning score. Pointwise
maximality is strict away from ties, so that assignment must
equal \(T\) almost everywhere under \(\mu\). The argument
holds for fractional assignments as well by integrating the
score with their label probabilities. The first heat derivatives
supply these moments, because
\[
 \int z\,1_{A_i}(z) f_x(z)\,dz
      =x P_s1_{A_i}(x)+s\nabla P_s1_{A_i}(x).
\]
Thus requiring exact first derivatives at even one observation
of the simplex precludes the variation used in (BL5).

There is a supermartingale bound using the first derivatives
that is sharp at the reference \(T\). Fix its affine scores.
Write \(h(z)\) for the largest score minus the second largest,
and \(r_A(z)\) for the largest score minus the score assigned
by \(A\). If the labels differ where \(h>0\), then \(r_A\ge h\).
For each \(\eta>0\),
\[
 1_{\{A(z)\ne T(z)\}}\le 1_{\{h(z)\le\eta\}}+r_A(z)/\eta .
\]
Along either Brownian endpoint the conditional expectation of
the right side is a martingale. Its infimum over positive
rational \(\eta\) is a nonnegative supermartingale bounded by
one, by taking \(\eta\to\infty\). The inequality for the infimum
follows by conditional expectation against each member of the
countable family. At time \(t\) this infimum uses only the
observation point, the masses and their first derivatives,
since the conditional expectation of \(r_A\) is determined by
those data. The resulting bound is continuous along the state of
each fixed partition at every \(0<t<1\). If \(A=T\) then \(r_A=0\) and
the infimum is zero. Otherwise \(P_s r_A(x)>0\) at every \(x\)
for \(s>0\), and near any fixed \(x\) sufficiently small
\(\eta\) can be discarded uniformly. As \(\eta\to\infty\) the
terms tend locally uniformly to one. This reduces continuity
to minimization over a compact \(\eta\)-interval, where
\(P_s1_{\{h\le\eta\}}(x)\) is jointly continuous for \(T\)
with distinct linear parts, since level sets of \(h\) have
Gaussian measure zero. The time \(1\) value can be defined by
the same infimum with the variables evaluated at the endpoint.

Add the mismatch supermartingales for \(x\) and \(y\) to the
conditional agreement martingale \(V_T(t,B_t,C_t)\). At the
terminal time their sum with \(V_T\) dominates the agreement
indicator for \(A\), by the union bound for a label mismatch
at either endpoint. For \(A=T\) both mismatch terms vanish,
so the initial value is exactly \(S_\rho(T)\). This bound can
vary with the initial first moments of a balanced \(A\).
It therefore supplies no common sharp initial value for all
balanced partitions. The proof above explains why any such
first-derivative approach that is sharp at \(T\) needs a
singularity at some simplex state, or a bound outside the
stated supermartingale class.
