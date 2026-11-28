# Stationarity of the adaptive Gaussian surface energy

Fix \(0<r<1\). Every balanced planar minimizer of the
adaptive surface functional with finite energy has a
stationary weighted perimeter obtained by averaging
attaining Gaussian tests. Its first variation is
determined by two volume pressures. This conclusion
uses only finite Gaussian perimeter and holds against
all compact smooth ambient variations.

At a junction with three regular arcs, the averaged
weights are strictly positive and satisfy the weighted
force balance. This regularity is a hypothesis on the
junction. The argument does not classify the
minimizers or establish junction regularity.
[Adaptive surface functional and dual tests](THREE_CELL_ADAPTIVE_SURFACE_VARIATIONAL_PROBLEM.md#a-dual-formula-with-metric-surface-weights)

## Tagged directions and their Gaussian tests

Let \(A\) be a planar ordinary partition of finite
Gaussian perimeter and finite energy
\(F=\mathcal F_r(A)>0\). Write \(\nu_e\) for its
Gaussian pair arclength and \(d_e=e_i-e_j\) for the
label vector of an oriented pair \(e=(i,j)\). Set
\(s=\sqrt r\). For direction fields with
\(|u_e|\le1\), define
\[
 Z(u)=T_s\sum_e d_e\otimes u_e\,\nu_e .
\]
The smoothing is applied coordinatewise. The space of
fields is the product of the weak-star compact unit
balls of \(L^\infty(\nu_e;\mathbb R^2)\).
Each \(T_s\nu_e\) is in Gaussian \(L^2\), by using one
constant unit direction on that pair and zero on the
others. Thus weak-star convergence of the direction
fields gives strong \(L^2\) convergence of \(Z(u)\).
Indeed it gives pointwise convergence against the
bounded source kernel at each observation, and the
sum of the scalar smoothings is an \(L^2\) majorant.
The unit balls here are metrizable in the weak-star
topology, since the underlying \(L^1\) spaces are
separable. In particular the set of normalized
attaining tests
\[
              \mathcal K=\{Z(u)/F:\|Z(u)\|_2=F\}
\]
is nonempty and compact in Gaussian \(L^2\).

For \(G\in\mathcal K\), put \(V_i=T_sG_i\) and
\(w_{G,e}=|V_i-V_j|\). These are the weights of its
dual perimeter. For an optimizing \(u\) with
\(G=Z(u)/F\), Hilbert space equality gives
\[
 F=\langle G,Z(u)\rangle
   =\sum_e\int (V_i-V_j)\cdot u_e\,d\nu_e
   =\sum_e\int w_{G,e}\,d\nu_e .
\]
The last equality follows by independently optimizing
the directions in the pairing with this fixed \(G\).
All the pairings are integrable, by Cauchy–Schwarz
against the scalar \(L^2\) smoothings. These tests
therefore attain the dual value. Conversely, an
attaining unit dual test can be paired with its
pointwise optimal direction field, and equality in
Hilbert space Cauchy–Schwarz identifies it as \(Z(u)/F\).

For every attaining \(G\), \(w_{G,e}>0\) almost
everywhere with respect to \(\nu_e\). Otherwise change
the direction field on a positive-measure zero set
of this weight on one pair. The change in the tagged
tensor measure can be made nonzero. Its smoothing
is also nonzero, by injectivity on finite signed
measures. Its pairing with \(G\) is zero, so the
squared norm would strictly increase. Thus
\(u_e=(V_i-V_j)/w_{G,e}\) almost everywhere on that
pair. Mehler injectivity used here follows by writing
the smoothing as the Laplace transform of a
Gaussian-damped finite measure and using Fourier
uniqueness after analytic continuation.

## The envelope first variation

Let \(X\) be a compact smooth ambient vector field
and let \(\Phi_t\) be its flow. Pull direction fields
on the transported interfaces back to the original
ones, retaining their vector values. This identifies
all admissible fields at time \(t\) with the same
unit balls at time zero. With \(K_x=M_s(\,\cdot\,,x)\)
as a vector in Gaussian \(L^2\), the transported
smoothing is
\[
 Z(t,u)=\sum_e\int K_{\Phi_t(x)}\,d_e\otimes u_e(x)
                      J^\gamma_{\tau_e}\Phi_t(x)\,d\nu_e(x).
\]
Here \(J^\gamma_{\tau_e}\Phi_t\) is the tangential
Jacobian times the ratio of Gaussian densities,
and \(\tau_e\) is the approximate tangent line.
The rectifiable area formula applies to the pair
interfaces of a finite-perimeter partition. At
\(t=0\) the Jacobian derivative is
\(\operatorname{div}_{\tau_e}X-x\cdot X\).
[Weighted area variation on reduced boundaries, Milman–Neeman, Proposition B.3](https://arxiv.org/html/1805.10961v3)

The map \(x\mapsto K_x\) is smooth into \(L^2\) on
compact sets, by differentiation of the Gaussian
kernel. Only a fixed compact part of the interfaces
moves. Consequently \(Z(t,u)=Z(u)+tZ_X(u)+o(t)\)
uniformly over \(|u_e|\le1\). Both \(Z(u)\) and
\(Z_X(u)\) are strongly continuous under weak-star
convergence of \(u\). For the derivative, its
compactly supported kernel is a sum of the continuous
\(L^2\)-valued kernels \(K_x,\nabla_xK_x\) times
bounded measurable coefficients. On the fixed
compact set the continuous kernels are uniformly
approximable by finite-rank kernels, giving the
claimed weak-star to strong continuity. The same
argument works for a smooth compact family of
diffeomorphisms with initial velocity \(X\).

The right derivative of the supremum of the norms
along this flow is the maximum of the derivatives
over the optimizers at zero. For the lower bound,
test each base optimizer. For the upper bound, take
optimizers at positive times tending to zero and
extract a weak-star subsequence in the base unit
balls. Uniform expansion and strong convergence
make its limit a base optimizer. The norm derivatives
converge along this sequence because the norms of
its base tensors tend to \(F>0\) and \(Z_X\) is strongly
continuous. For each base optimizer, differentiation
and the alignment of \(u_e\) with \(V_i-V_j\) give
\[
 \begin{split}
 L_G(X)&=\langle G,Z_X(u)\rangle\\
 &=\sum_e\int
 \bigl[\nabla w_{G,e}\cdot X+
          w_{G,e}(\operatorname{div}_{\tau_e}X-x\cdot X)\bigr]\,d\nu_e .
 \end{split}                                                  \tag{ASV1}
\]
The right derivative is therefore
\(\max_{G\in\mathcal K}L_G(X)\).
The expression (ASV1) is linear in \(X\). For fixed
\(X\) it is continuous on \(\mathcal K\):
\(L^2\) convergence of \(G\) gives local \(C^1\)
convergence of \(V_i\), and the norm gradient
converges wherever the limiting weight is nonzero.
These points have full \(\nu_e\)-measure, and the
gradients and weights are bounded uniformly on the
support of \(X\). The same formula is the first
variation of the ordinary metric perimeter with
the weights \(w_{G,e}\) held fixed as functions
of the observation point.

## A common first variation and the volume pressures

Now let \(A\) be balanced and minimizing. Write
\(M_i(X)=\int_{A_i}\operatorname{div}_\gamma X\,d\gamma\)
for the initial cell mass changes. The linear map
\(M\) from compact smooth fields is onto the plane
\(\sum_i M_i=0\). In fact, a vector annihilating its
range makes the distributional gradient of
\(\sum_i a_i1_{A_i}\) zero, by absorbing the positive
Gaussian density into the test vector field. This
function must be constant, so the positive cell
masses force every \(a_i\) to be equal.

Pick two smooth fields whose mass changes span that
plane. After the flow of any \(X\) with \(M(X)=0\),
compose with the flows of these two fields.
The mass map of this composition is smooth by the
Gaussian change of variables on the fixed cells.
The implicit function theorem gives correction
times with zero first derivatives that preserve the
masses exactly. Hence minimizing and the envelope
formula imply \(\max_{G\in\mathcal K}L_G(X)\ge0\) for
every \(M(X)=0\).

For any finite list of such fields, zero belongs
to the convex hull of their simultaneous evaluations
\((L_G(X_1),\ldots,L_G(X_k))\). Otherwise strict
finite-dimensional separation would give a linear
combination of the fields whose maximum evaluation
was negative. Thus there is a probability measure
on \(\mathcal K\) annihilating those finitely many
evaluations. The probability measures on the compact
metric space \(\mathcal K\) form a compact space in
the weak topology. Continuity of the evaluations
and the finite intersection property give a single
probability measure \(\pi\) with
\(\int L_G(X)\,d\pi(G)=0\) whenever \(M(X)=0\).
By surjectivity of \(M\), there is a unique
centered vector \(\beta\) for this choice of \(\pi\)
such that
\[
       \int_{\mathcal K} L_G(X)\,d\pi(G)
                   =\sum_{i=1}^3 \beta_i M_i(X)
                   \qquad\text{for every compact smooth }X.   \tag{ASV2}
\]
In particular the averaged first variation depends
only on two independent mass changes.

## Positive tensions at regular junctions

First consider an arbitrary finite-energy partition,
with any attaining \(G=Z(u)/F\). Fix one pair \(e\),
a compactly supported bounded real function \(f\) on
it, and the planar quarter-turn \(J\). Change only
that pair direction to
\(\sqrt{1-t^2 f^2}\,u_e+t f J u_e\).
Since \(u_e\) is unit and attaining, the quadratic
coefficient of the squared tensor norm is
nonpositive. Using \(|d_e|^2=2\) and
\(V_i-V_j=w_{G,e}u_e\) gives
\[
 2\iint M_r(x,y)f(x)f(y)\,u_e(x)\cdot u_e(y)
                                \,d\nu_e(x)d\nu_e(y)
       \ \le\ F\int w_{G,e} f^2\,d\nu_e .                     \tag{ASV3}
\]
Partition the unit directions into finitely many
angular sectors on which all inner products are
at least \(1/2\). On one sector and a compact set
\(D\), take \(f=(w_{G,e}+\eta)^{-1}\) there and zero
elsewhere. The left side of (ASV3) is at least
\(\min_{D\times D}M_r\,(\int f\,d\nu_e)^2\).
Its right side is at most \(F\int f\,d\nu_e\).
Let \(\eta\) tend to zero and sum over the sectors.
This gives
\[
       \int_D \frac{d\nu_e}{w_{G,e}}
                   \le C(D,r) F .                            \tag{ASV4}
\]
In particular \(w_{G,e}(p)>0\) at every point of
positive lower one-dimensional density of \(\nu_e\).
If it vanished there, its local Lipschitz bound
would give \(w_{G,e}(x)\le C|x-p|\), and the integral
in (ASV4) would diverge by the lower density bound.
The reciprocal integral assertion and this positivity
hold for any finite-energy partition. At a fixed
such point, compactness of \(\mathcal K\) and
continuity in \(G\) give uniform positivity over
all attaining tests.

For the minimizing partition define
\(\bar w_e(x)=\int_{\mathcal K}w_{G,e}(x)\,d\pi(G)\).
These weights are locally Lipschitz by the
uniform local \(C^1\) bounds on \(T_sG_i\).
Suppose near a finite vertex \(p\) the interface
consists of three \(C^2\) arcs up to their common
endpoint \(p\), with outward unit tangents
\(\tau_e(p)\). On every incident pair
\(\bar w_e(p)>0\) by (ASV4). To localize (ASV2),
choose \(X\) equal to an arbitrary constant
vector \(v\) near \(p\), supported in \(B_\epsilon(p)\)
and bounded by \(|v|\). On each arm write
\(\operatorname{div}_{\tau_e}X=
 \partial_\ell(X\cdot\tau_e)-X\cdot\partial_\ell\tau_e\),
where \(\ell\) is arclength outward from \(p\).
Integration by parts of the tangential derivative
gives the endpoint contribution
\(-\phi_2(p)\bar w_e(p)\,v\cdot\tau_e(p)\).
All other terms in (ASV1) after integration by
parts are \(O(\epsilon |v|)\) by the local bounds
on curvature, weights and their first derivatives.
The pressure terms in (ASV2) are also
\(O(\epsilon |v|)\), by the Gaussian boundary
flux and the lengths of the three arcs.
Letting \(\epsilon\) tend to zero proves
\[
               \sum_{e\ni p}\bar w_e(p)\tau_e(p)=0.           \tag{ASV5}
\]
The averaged tensions in this force balance depend on
the attaining tests selected by \(\pi\). They need
not be equal.

On a regular \(C^2\) arc, the corresponding weight
\(\bar w_e\) is smooth in a neighborhood of every
interior point. Indeed the weight is positive there
uniformly in \(G\), by (ASV4), compactness, and the
uniform local derivative bound. The smoothed fields
have uniform bounds for derivatives of every fixed
order on compact sets, so the average of the norms
can be differentiated near the arc.

Orient the normal \(N_e\) from \(i\) into \(j\), and
write \(\kappa_e\) for the signed curvature defined
by \(\partial_\ell\tau_e=\kappa_eN_e\). For fields
supported near an interior patch of this arc,
integration by parts of the tangential terms in
(ASV1) gives the pressure equation
\[
      \partial_{N_e}\bar w_e
           -\bar w_e(\kappa_e+x\cdot N_e)=\beta_i-\beta_j .
                                                               \tag{ASV6}
\]
Here the outward normal component \(X\cdot N_e\) enters
the volume derivative \(M_i(X)\) of cell \(i\) with a
positive sign.
The Gaussian density derivative accounts for the
term \(x\cdot N_e\). This equation is conditional on
an existing regular arc.
