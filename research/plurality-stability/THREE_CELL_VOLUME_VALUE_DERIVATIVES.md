# Volume derivatives and label prices at Gaussian maximizers

Fix \(0<\rho<1\) and a finite Gaussian dimension \(d\ge2\). Let
\(\mathcal V_\rho(m)\) be the largest three-cell stability at Gaussian
masses \(m=(m_1,m_2,m_3)\) in the interior of the probability simplex.
For a maximizer \(A\), center its label prices \(\beta(A)\) so that
\(\sum_i\beta_i(A)=0\), with \(A_i\) the strict winners of
\(T_\rho1_{A_i}-\beta_i(A)\). Then for every \(v\) with \(\sum_i v_i=0\),
\[
  \lim_{t\downarrow0}\frac{\mathcal V_\rho(m+tv)-\mathcal V_\rho(m)}t
       =2\max_{A\text{ maximizing at }m}\beta(A)\cdot v.             \tag{VV1}
\]
The limit also holds for directions tending to \(v\) as \(t\downarrow0\).
In particular, at \(m_*=(1/3,1/3,1/3)\), differentiability of the optimal
value is equivalent to every maximizer having equal prices.
If one maximizer at balance has a nonzero centered price, then
\[
  \lim_{t\downarrow0}\frac{\mathcal V_\rho(m_*+tv)-\mathcal V_\rho(m_*)}t
        \ \ge\ |\beta(A)|\,|v| .                                     \tag{VV2}
\]
Thus a tied nonsimplex maximizer with unequal prices would force
strictly positive outgoing slopes of the volume-value function in every
direction. The statements here do not establish differentiability or
locate an adjusted-score junction.

## Compactness and the price vector

Work in the weakly compact class of simplex-valued functions at masses
\(m\), which includes ordinary partitions. The Ornstein–Uhlenbeck
operator \(T_\rho\) is compact on Gaussian \(L^2\) in fixed dimension:
its degree-\(n\) Hermite eigenspaces are finite-dimensional and the
eigenvalues \(\rho^n\) tend to zero. Consequently stability is weakly
continuous. A maximizer is ordinary. Indeed, a nonordinary function
admits a nonzero bounded perturbation in opposite directions, supported
where two coordinates are uniformly positive, with zero mean in every
coordinate. Strict positivity of the quadratic form of \(T_\rho\) makes
the average of the two perturbed stabilities strictly larger.

For a maximizer \(A\), the linear functional with scores
\(U_i=T_\rho1_{A_i}\) is maximal at \(A\) among simplex-valued functions
of the same masses. Expand the stability at \(A\) toward such a
competitor and use the nonnegative quadratic term. As a function of
the prescribed masses, the optimal value of that linear functional is
concave by taking mixtures. Its supergradient at the interior point
\(m\) yields prices \(\beta_i\). The pointwise maximization of
\(U_i-\beta_i\) is then attained by \(A\) almost everywhere. The winners
are strict almost everywhere by analyticity and injectivity of
\(T_\rho\): a constant pair contrast would imply that two cell
indicators differ by a constant, contrary to the positive masses of all
three labels. The positive mass of each winner and \(0\le U_i\le1\)
give \(|\beta_i-\beta_j|\le1\).

## Ambient transport and uniform mass correction

For a compactly supported smooth vector field \(X\), let \(A_s\) be the
image of \(A\) under its flow. The volume velocity is
\[
  (\mathcal M_A X)_i
       =\int_{A_i}\operatorname{div}(\phi_d X)\,dx ,
\]
where \(\phi_d\) is the standard Gaussian density. A change of variables
in the joint Gaussian integral gives the stability velocity
\[
  \left.\frac{d}{ds}S_\rho(A_s)\right|_{s=0}
      =2\sum_i\int_{A_i}\operatorname{div}(\phi_d U_i X)\,dx
      =2\beta(A)\cdot\mathcal M_A X .                                \tag{VV3}
\]
For the second equality, write \(W=\max_i(U_i-\beta_i)\). It is locally
Lipschitz and agrees together with its weak gradient with the winning
score on each strict winner set. Subtracting the price term therefore
leaves \(2\int\operatorname{div}(\phi_d W X)\,dx=0\).

The flux map onto the label sum-zero space is surjective. If a vector
\(c\) annihilates its range, the locally integrable function
\(\sum_i c_i1_{A_i}\) has zero distributional gradient, by testing with
\(\phi_dX\). It is constant and each \(A_i\) has positive mass, so all
\(c_i\) agree. Equation (VV3) also proves uniqueness of the centered
prices, since the same equation holds for any winning-score prices.

Choose two compact smooth fields whose fluxes at \(A\) are independent.
The joint flow of these fields, in two small parameters, has an
invertible volume derivative at zero. For partitions converging to
\(A\) in \(L^1(\gamma_d)\), this derivative remains uniformly
invertible, since it depends on integrals against bounded test
functions. The mass and stability maps under the two flows have
uniformly bounded second derivatives in a fixed small parameter ball.
This follows by the change of variables formulas, because the smooth
flow is the identity off a compact set and derivatives of the joint
Gaussian density there are bounded in absolute value by an integrable
Gaussian times a polynomial. The inverse function theorem with
uniform derivative bounds supplies an exact prescribed small mass
correction of size \(h\) using flow parameters of size \(O(|h|)\).
Equation (VV3) and the two Taylor bounds give a corresponding
stability change of \(2\beta\cdot h+O(|h|^2)\) at a maximizer.
The estimates apply uniformly along a strongly convergent sequence
of maximizers using the two fields at its limit.

## Directional derivative at interior masses

Transport any maximizer at \(m\) to \(m+tv\), using the preceding
correction. This gives the lower bound in (VV1) for that maximizer's
price. For the upper bound, let \(t_n\downarrow0\), \(v_n\to v\),
and take maximizers \(A_n\) at \(m+t_nv_n\). A weakly convergent
subsequence has a limit at masses \(m\) and its energies converge.
The lower bound from a fixed maximizer at \(m\) implies that the
limit is also maximizing. It is ordinary, so the weak convergence is
strong: the sum of squared \(L^2\) norms of the indicator vector is one
both along the sequence and at the limit.

The bounded centered prices have further convergent subsequences.
Along each one the score functions converge locally uniformly to the
limit partition's scores. For example, each Mehler kernel row is
bounded in \(L^2(\gamma_d)\) uniformly for its first argument on a
compact set, so strong \(L^2\) convergence suffices. Passing to
almost-everywhere convergence of the indicators shows that the limiting
prices are winning-score prices for the limit partition. By uniqueness,
these are its centered prices. The same argument at fixed masses shows
that the set of maximizing prices there is compact, so the maximum in
(VV1) is attained.

Use the two fields of the limit partition to correct each \(A_n\)
exactly to masses \(m\). The corrected partition has stability at most
\(\mathcal V_\rho(m)\) and stability change
\(-2t_n\beta(A_n)\cdot v_n+O(t_n^2)\). On a subsequence achieving the upper
limit, this proves the required upper bound in (VV1).

At balance, cyclic relabeling rotates the centered prices through
angles \(0,2\pi/3,4\pi/3\) in their two-dimensional label plane.
One of their dot products with \(v\) is at least
\(|\beta(A)|\,|v|/2\). Equation (VV1) gives (VV2). A function
differentiable at balance and invariant under all label permutations
has zero derivative on that plane, so its maximizing prices must all
vanish. Conversely, if all prices vanish, the version of (VV1) with
varying directions and compactness of the unit circle give
\(\mathcal V_\rho(m_*+h)-\mathcal V_\rho(m_*)=o(|h|)\).
