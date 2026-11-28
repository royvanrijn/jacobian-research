# Finite score jets at an off-center Gaussian junction

Fix \(0<r<1\) and an integer \(m\ge1\). For every sufficiently small
nonzero \(p\in\mathbb R^2\) there is an ordinary balanced three-cell
partition \(A\) whose interfaces form a smooth Y with vertex \(p\).
Near \(p\) the three rays are exactly those of the central regular
simplex translated by \(p\). Outside a compact set the partition is
the original centered simplex.

Write \(U_i=T_r1_{A_i}\). On each germ ray orient the corresponding
pair as \(e=(i,j)\), and use its outgoing unit tangent \(\tau_e\). Then
\[
 U_i(p)=\tfrac13\quad(i=1,2,3),\qquad
 \partial_{\tau_e}^{\,\ell}(U_i-U_j)(p)=0
                \quad(e=(i,j),\ 1\le \ell\le m).
                                                               \tag{OJ1}
\]
The normal derivatives of the pair contrasts have the stationary
sign at the vertex. The construction also makes one tangential
derivative of order \(m+1\) nonzero. Thus exact score equality at
a balanced vertex and agreement of the stationary jets to any
fixed order allow a noncentral junction.

## A family of genuine Y partitions

Let \(T\) be the centered simplex. Orient its pairs cyclically and let
\(d_e=\mathbf e_i-\mathbf e_j\) for \(e=(i,j)\), with \(N_e\) from \(T_i\)
into \(T_j\). Fix a smooth compactly supported function \(\chi\) that
equals one in a disk about the origin. For small \(p\) the map
\(x\mapsto x+\chi(x)p\) is a diffeomorphism. It translates the germ and
fixes all points outside its support.

On a compact interval of each ray outside this support, choose a thin
normal tube disjoint from the other rays. A compact smooth height
\(h_f(s)\) in the interval slides ray \(f\) at \(s\tau_f\) by
\(h_f(s)N_f\). This motion is realized by a diffeomorphism that is the
identity off the tube, using a smooth cutoff of the normal coordinate
equal to one near the ray. For small heights, compose the three tube
maps with the vertex translation. The resulting partitions are all
diffeomorphic images of \(T\), so they have exactly the stated Y
topology and the exact straight germ and exterior. Their Gaussian
masses and score derivatives at \(p\) vary smoothly in \(p\) and in
any finite list of bump coefficients. This follows by integrating
smooth Gaussian densities under the compact deformations.

## Independent wall moments for the score jets

Set \(\sigma^2=1-r^2\) and let \(\varphi_\sigma\) denote the planar
centered Gaussian density of covariance \(\sigma^2 I\). Thus the
conditional source density at observation \(x\) is
\(\varphi_\sigma(y-rx)\). Linearize at \(p=0\) and zero heights. An
outward bump on ray \(f\) gives the indicator vector variation
\(d_f h_f(s)\,ds\) on that ray. Hence the first variation of the order
\(\ell\) tangential score derivative on ray \(e=(i,j)\) is
\[
 (r/\sigma)^\ell\sum_f(d_e\cdot d_f)
   \int_0^\infty h_f(s)\,\varphi_\sigma(s\tau_f)\,
       H_\ell\!\left(s\,\tau_e\cdot\tau_f/\sigma\right)\,ds ,
                                                               \tag{OJ2}
\]
where \(H_\ell\) is the monic probabilists' Hermite polynomial.
This follows by differentiating the conditional source density
\(\ell\) times in direction \(\tau_e\).

For each \(f\), its standard Gaussian mass moment
\(\int h_f(s)\varphi_1(s\tau_f)\,ds\) and its tilted moments
\(\int h_f(s)s^\ell\varphi_\sigma(s\tau_f)\,ds\), \(0\le\ell\le m+1\),
can be varied independently through compact smooth bumps in the
chosen interval. Otherwise a linear combination of their weights
would vanish on that interval. By analyticity it would vanish on the
positive ray. Dividing by \(\varphi_1(s\tau_f)\) and taking
\(s\to\infty\) kills the polynomial times the faster Gaussian and
forces the standard weight coefficient to be zero. The polynomial
then vanishes as well. In particular finitely many bump directions
realize a right inverse for these moments.

The standard moments control the two independent ordinary mass
variations through the spanning vectors \(d_f\). The order-zero
tilted moments control the two independent variations of
\((U_i(p))_i\) by the same incidence vectors. Equivalently the three
cyclic score equalities have rank two at order zero. For every
\(\ell\ge1\), once the lower moments have been chosen, the highest
moment terms in (OJ2) have coefficient matrix
\[
  (r/\sigma^2)^\ell\big[(d_e\cdot d_f)
                               (\tau_e\cdot\tau_f)^\ell\big]_{e,f}.
                                                               \tag{OJ3}
\]
The bracketed matrix has diagonal entries \(2\) and off-diagonal
entries \((-1)^{\ell+1}/2^\ell\). It is invertible by strict diagonal
dominance. Successively prescribing orders \(1,\ldots,m+1\) therefore
gives a right inverse to the full constraint derivative.

At \(p=0\) the unperturbed simplex has both mass vectors equal to
\((1/3,1/3,1/3)\). Each score contrast is odd under reflection in its
germ line and therefore vanishes on that line with all tangential
derivatives. Choose a square invertible restriction of the bump
derivative and apply the implicit function theorem. Keep the two
mass vectors exactly balanced and prescribe zero jets through order
\(m\). Prescribe \(|p|^2\) at order \(m+1\) on one chosen ray and
zero on the other two. The constraints have their unperturbed values
at \(p=0\), so the construction works for every \(p\) sufficiently
near zero. The strict normal score signs at the simplex vertex
persist by continuity. At the simplex these signs follow by
differentiating the conditional Gaussian density in direction \(N_e\).
The integrand is proportional to
\((N_e\cdot y)(1_{T_i}-1_{T_j})(y)\), which is nonpositive and
strictly negative on a set of positive measure. The nonzero
order-\((m+1)\) jet makes the contrast nonconstant on the chosen
straight germ ray, so the stationary winner rule fails there.
