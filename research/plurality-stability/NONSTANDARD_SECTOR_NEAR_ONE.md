# Nonstandard simplex spectra near correlation one

The ordered sign kernel converges to the Gaussian Dirichlet problem on
the ordered facet. The ordering walls impose zero trace through the
alternating fold. The top-tie wall imposes the same condition through
the loss of row mass to the adjacent facets. A uniform comparison of
this loss with the probability of crossing the wall controls the
zero extension across the entire chamber. For nonstandard two-row
shapes the full facet has a Gaussian Neumann limit in the equivariant
space. Folding introduces mixed wall conditions.

Let \(k\ge3\), \(n=k-2\), \(H=\mathbf 1^\perp\cap\{x_1=x_2\}\), and
\[
 F=\{x\in H:x_2>x_l,\ 3\le l\le k\},\qquad
 W=\{x\in H:x_2>x_3>\cdots>x_k\}.
\]
For \(k=3\) the ordering condition adds no wall. Write \(\gamma=\gamma_H\)
for standard Gaussian measure on \(H\), and let
\(L^2(W)=L^2(W,\gamma)\). Restriction and alternating extension identify
the active sign multiplicity space with this scalar space by (S44)–(S46)
of the [interface reduction](SIMPLEX_INTERFACE_SYMMETRY.md). The common
factor in the facet norm and form disappears from their quotient.
Sign speeds already lie in the balanced space.

Let \(g_d^{\mathrm{sign}}(\rho)\) be the \(d\)-th min–max level of the
interface form on this multiplicity space, with \(d\ge1\) denoting the
level index. Put \(\sigma=\sqrt{1-\rho^2}\) and
\(\phi_0=(2\pi)^{-1/2}\). Denote the Gaussian Dirichlet
eigenvalues on \(W\) by \(\lambda_d^D\), with the form domain given by
the closure of \(C_c^\infty(W)\) in \(W^{1,2}(W,\gamma)\).

**Theorem (Sign-sector limiting levels).** For each fixed
\(k\ge3\) and level \(d\ge1\),
\[
 \lim_{\rho\uparrow1}\frac{\rho\,\sigma}{2\phi_0(1-\rho)}
       g_d^{\mathrm{sign}}(\rho)=\lambda_d^D-1 .
 \tag{SG1}
\]
For \(\rho\) sufficiently close to one the first \(d\) min–max levels
are discrete eigenvalues. For \(k\ge4\), let \(\Lambda_1^D,
\Lambda_2^D,\ldots\) be the spherical Dirichlet eigenvalues on
\(W\cap S^{n-1}\), with multiplicities, and let
\(\beta_j^D>0\) solve
\(\beta_j^D(\beta_j^D+n-2)=\Lambda_j^D\). The limits in (SG1) are the
ordered multiset
\[
             \{\beta_j^D+2m-1:j\ge1,\ m\ge0\}.
 \tag{SG2}
\]
In four cells the normalized first limit is
\(\pi/\arctan\sqrt2-1\), and all levels are
\(j\pi/\arctan\sqrt2+2m-1\) for \(j\ge1\) and \(m\ge0\).
For \(k=3\) the limits are \(0,2,4,\ldots\).
For a cylindrical simplex in ambient dimension \(N\ge k-1\), each
active limiting level is augmented by a nonnegative Hermite degree
\(r\) from the \(s=N-(k-1)\) unused coordinates. If \(s\ge1\), degree
\(r\) is counted with multiplicity \(\binom{s+r-1}{r}\). If \(s=0\),
only \(r=0\) is used.

## Pointwise comparison and loss of row mass

Write \(m_\rho(x,y)\) for the Mehler kernel on \(H\) relative
to \(\gamma\), and for a set \(A\subset H\), let
\(m_\rho^A(x)=\int_A m_\rho(x,y)\,d\gamma(y)\).
Use the planar adjacent kernels and normalized diagonal from the
[standard-sector calculation](STANDARD_SECTOR_NEAR_ONE.md):
\[
 w_l=(e_1+e_2)/2-e_l,\qquad
 m_{\rho,l}(x,y)=m_\rho(x,y)
    e^{-\rho(w_l\cdot x)(w_l\cdot y)/\sigma^2},\qquad
 b_\rho=\frac{\sigma a_\rho}{2\phi_0\rho}
       =m_\rho^F+\frac12\sum_{l=3}^k m_{\rho,l}^F.
\]
In particular \(p_F=\gamma(F)\le b_\rho\le1\) on \(F\) and
\(b_\rho\to1\) in its interior. The lower bound \(m_\rho^F\ge\gamma(F)\) follows from \(F+\rho x\subseteq F\).
The upper bound and interior limit are the
score coefficient estimates from (B7) and (SN2).

Let \(k_\rho=(\sigma/2)\mathcal K_\rho\), where \(\mathcal K_\rho\)
is the canonical folded kernel in (S46). The scaled kernel
\(k_\rho\) and planar kernel \(m_\rho\) act against the same measure
\(\gamma\) on \(W\). The identity term in (S46) gives exactly \(m_\rho\)
after multiplying by \(\sigma/2\). To compare the determinant with
this term, use the strict total positivity (S47) and apply (S63)
repeatedly, dropping the first row and column each time. The
determinant is positive and at most the product of its diagonal
entries. Consequently
\[
                  0<k_\rho(x,y)\le m_\rho(x,y) \qquad (x,y\in W).
 \tag{SG3}
\]
Also \(k_\rho/m_\rho\to1\) uniformly on \(U^2\) for each compact
\(U\subset W\), with error \(O_U(e^{-c_U/\sigma^2})\) as
\(\rho\uparrow1\). Indeed, put \(a=(x_2,\ldots,x_k)\) and
\(b=(y_2,\ldots,y_k)\). Dividing the determinant by its diagonal
product, each nonidentity term has absolute value
\[
          \exp[-(\rho/\sigma^2)(a\cdot b-a\cdot\pi b)].
\]
Strict rearrangement for decreasing vectors gives a positive
minimum for all finitely many nonidentity losses on \(U^2\).

Write \(k_\rho^W(x)=\int_W k_\rho(x,y)\,d\gamma(y)\) and
\(\delta_\rho=b_\rho-k_\rho^W\). From (SG3),
\[
 \delta_\rho\ge m_\rho^{F\setminus W}
                   +\tfrac12\sum_l m_{\rho,l}^F .
 \tag{SG4}
\]
In fact there is an explicit \(c_k>0\), independent of \(\rho\), with
\[
              \delta_\rho(x)\ge c_k m_\rho^{H\setminus W}(x).
 \tag{SG5}
\]
For this, put \(z=\rho x/\sigma\) and choose \(l\) minimizing
\(v=w_l\cdot z\). The union bound gives
\(m_\rho^{H\setminus F}(x)\le n\Phi(-u)\), where
\(u=v/\sqrt{3/2}\).
Completing the Gaussian square with \(D_l=I-w_lw_l^\top\) gives
\[
 m_{\rho,l}^F(x)=e^{-v^2/4}\Pr(Z+D_l z\in F),
              \qquad Z\sim\gamma.
\]
The \(l\)-gap of \(D_l z\) is \(-v/2\) and its other gaps are at least
\(v/2\), since \(w_l\cdot w_l=3/2\) and \(w_i\cdot w_l=1/2\)
for \(i\ne l\). Conditional on \(w_l\cdot Z=s\ge v/2\),
all other Gaussian gaps have mean \(s/3\ge0\). Their conditional
centered Gaussian vector has covariance \(I+J/3\). Its nonnegative
orthant has probability at least \(2^{-(n-1)}\), as is seen by
conditioning on its common Gaussian summand and applying Jensen to
the \((n-1)\)-st power of the one-coordinate conditional probability.
The empty-vector case \(n=1\) has probability one. The residual
nonnegative orthant forces all other threshold inequalities. Thus
\[
 \begin{split}
 m_{\rho,l}^F(x)
   &\ge 2^{-(n-1)}e^{-3u^2/8}\Phi(-u/2)\\
   &=2^{-(n-1)}\int_u^\infty
                \phi(s)e^{su/2-u^2/2}\,ds
     \ \ge\ 2^{-(n-1)}\Phi(-u).
 \end{split}
\]
The shift in the second line is \(s=t+u/2\), starting from the
Gaussian integral over \(t\ge u/2\). Since \(u\ge0\), its exponential
factor is at least one. Equation (SG4) now gives (SG5) with
\(c_k=2^{-n}/n\).

## Compactness and the limiting form

For a nonnegative symmetric kernel \(m\) on \(W^2\), set
\[
 J_m(h)=\frac12\iint_{W^2}m(x,y)[h(x)-h(y)]^2d\gamma(x)d\gamma(y).
\]
Let \(Q_\rho=\int_W b_\rho h^2d\gamma-
\iint_{W^2} k_\rho(x,y)h(x)h(y)d\gamma(x)d\gamma(y)\). The normalized
quadratic form in (SG1), on \(L^2(W)\), is
\[
 \frac{Q_\rho(h)}{1-\rho}-\int_W b_\rho h^2d\gamma,\qquad
 Q_\rho(h)=J_{k_\rho}(h)+\int_W\delta_\rho h^2d\gamma .
 \tag{SG6}
\]
For the factor check, the original diagonal is \(a_\rho/\rho\)
and the interaction is \(\mathcal K_\rho/\rho\) against
\(\nu=\phi_0\gamma\). Passing to the \(L^2(\gamma)\) quotient and
multiplying by \(\rho\sigma/(2\phi_0(1-\rho))\) gives diagonal
\(\rho b_\rho/(1-\rho)\) and interaction \(k_\rho/(1-\rho)\).

For the zero extension \(\widetilde h\) to \(H\), the global Mehler
Dirichlet form is
\[
 \langle\widetilde h,(I-T_\rho)\widetilde h\rangle
              =J_{m_\rho|_{W^2}}(h)+\int_W m_\rho^{H\setminus W} h^2d\gamma .
\]
Equation (SG3) and \((r-s)^2\le2(r^2+s^2)\) give
\[
 J_{m_\rho|_{W^2}}(h)
 \le J_{k_\rho}(h)+2\int_W(m_\rho^W-k_\rho^W)h^2d\gamma
 \le 2Q_\rho(h).
\]
Using (SG5) for the other term, we obtain the global bound
\[
  \langle\widetilde h,(I-T_\rho)\widetilde h\rangle
                       \le (2+c_k^{-1})Q_\rho(h).
 \tag{SG7}
\]
This estimate has no sign restriction on \(h\).

Take \(\rho_i\uparrow1\) and suppose \(h_i\) is bounded in \(L^2(W)\)
with bounded above normalized energy in (SG6). Since \(0\le b_\rho\le1\),
(SG7) bounds the full Gaussian Mehler energy divided by \(1-\rho_i\).
In the Hermite basis on \(H\) its degree-\(r\) multiplier is
\((1-\rho_i^r)/(1-\rho_i)\). At degrees at least \(N\) it is bounded
below by \((1-\rho_i^N)/(1-\rho_i)\to N\). Finite-dimensional
compactness and the tail bound therefore give a strongly convergent
subsequence of the zero extensions in \(L^2(H,\gamma)\). Taking a
lower limit on each finite Hermite sum shows that the limit belongs
to the full Gaussian \(W^{1,2}(H,\gamma)\); it is zero almost
everywhere outside \(W\). Its restriction \(h\) belongs to
\(W_0^{1,2}(W,\gamma)\). To verify the latter assertion directly,
cut off radially in the full Gaussian Sobolev norm, translate the
compactly supported zero extension by a small vector in \(W\),
and mollify at a smaller scale. Each defining inward wall functional
becomes uniformly positive on the translated support, so the
mollifications are supported compactly in \(W\). Strong Sobolev
continuity follows on bounded sets, where the Gaussian density
is comparable to Lebesgue measure.

The sharp gradient lower limit follows from \(k_\rho/m_\rho\to1\)
on compact subsets of \(W^2\). More explicitly, fix a relatively
compact open \(U\subset W\), take \(L<\infty\), and retain only
the pairs \(x\in U\), \(y=x+\sigma z\), \(|z|\le L\) in \(J_{k_\rho}\).
For \(\rho\) close to one all those pairs lie in a slightly larger
compact interior set. The lower bound is \(1-o(1)\) times the Mehler
one, where the coefficient of the squared difference quotients
\((h_i(x+\sigma_i z)-h_i(x))/\sigma_i\) is
\(\sigma_i^2/(2(1-\rho_i))\to1\). The conditional density of \(z\)
converges uniformly to the standard Gaussian density on
\(U\times B_L\). Strong \(L^2\) convergence identifies the weak
limit of these quotients with \(z\cdot\nabla h(x)\) by changing
variables against compactly supported tests, as in (SN7)–(SN8).
Weak lower semicontinuity and exhaustion in \(U,L\) give
\[
 \liminf_i \frac{Q_{\rho_i}(h_i)}{1-\rho_i}
                         \ge \int_W|\nabla h|^2\,d\gamma .
\]
Boundedness and interior convergence of \(b_\rho\) plus strong
\(L^2\) convergence make the potential in (SG6) converge to
\(-\|h\|^2\).

For each fixed \(h\in C_c^\infty(W)\), the row defect
\(\delta_\rho\) is exponentially small uniformly on \(\operatorname{supp}h\)
as \(\rho\uparrow1\). Expand it as
\(m_\rho^{F\setminus W}+\tfrac12\sum_l m_{\rho,l}^F
 +(m_\rho^W-k_\rho^W)\). The first term has exponentially
small Gaussian tails because the support is separated from all
walls. The square-completion factor \(e^{-(\rho w_l\cdot x/\sigma)^2/4}\)
bounds each adjacent row. For the final term, use the uniform
ratio estimate in (SG3) when \(y\) lies in a compact neighborhood
of the support in \(W\) and Gaussian tails for other \(y\).
Also \(J_{k_\rho}(h)\le J_{m_\rho|_{W^2}}(h)\), and the latter
divided by \(1-\rho\) tends to \(\int_W|\nabla h|^2d\gamma\)
by the fixed-test difference-quotient calculation of (SN3).
This proves recovery on a dense set of Dirichlet tests.

## Spectral consequences

The embedding of \(W_0^{1,2}(W,\gamma)\) into \(L^2(W)\) is
compact by zero extension and the full-space Hermite expansion.
At fixed \(\rho\), the kernel \(k_\rho\) is Hilbert–Schmidt by
(SG3) and the full-space Mehler expansion, whose sum of squared
eigenvalues is \((1-\rho^2)^{-n}\). Thus the normalized operator
in (SG6) is a compact perturbation of multiplication by
\(\rho b_\rho/(1-\rho)\). Its essential spectrum is bounded
below by \(\rho p_F/(1-\rho)\), which tends to infinity.

To see (SG1), approximate the first \(d\) Dirichlet eigenfunctions
in the form norm by compactly supported smooth tests. Recovery
and polarization on that fixed finite-dimensional span give the
min–max upper limit. For the lower limit, select orthonormal
bases of almost minimizing \(d\)-planes along a lower-limit sequence.
Apply (SG7) and above compactness to their finitely many members.
The strongly convergent limits are orthonormal, and the sharp
lower bound applies to every fixed linear combination. The limiting
\(d\)-plane has Rayleigh supremum at most the lower limit, as needed.
Since the levels stay bounded, the essential spectrum estimate
gives discreteness for \(\rho\) close to one.

For \(n\ge2\) separate the Gaussian Dirichlet problem into the
angular Dirichlet eigenbasis. With radial density
\(w(r)=r^{n-1}e^{-r^2/2}\), a single angular mode has form
\(\int_0^\infty(|h'|^2+\Lambda_j^D r^{-2}|h|^2)w\,dr\).
The substitution \(h=r^{\beta_j^D}v\) gives the shift
\(\beta_j^D\) and derivative form in radial weight
\(r^{n+2\beta_j^D-1}e^{-r^2/2}\). Substituting \(s=r^2/2\)
then gives the Laguerre form and levels \(\beta_j^D+2m\),
by the proof of (SN11) with the angular Dirichlet basis.
The finite-energy cutoffs and Gamma-polynomial density argument
are unchanged. In four cells \(W\) is a wedge of opening
\(\alpha=\arctan\sqrt2\). Its spherical Dirichlet eigenvalues are
\((j\pi/\alpha)^2\). In three cells \(W\) is a half-line,
so odd Gaussian Hermite extension gives the Dirichlet levels
\(1,3,5,\ldots\). These prove the explicit claims.

For every \(k\ge4\), the first angular exponent in (SG2) is
strictly greater than \(e_k=(k-2)(k-3)/2\). To see this, remove the
top wall and keep only the ordered untied-coordinate walls.
The positive harmonic Vandermonde
\(\prod_{3\le i<j\le k}(x_i-x_j)\) vanishes on those walls and
gives an angular Dirichlet eigenfunction of degree \(e_k\).
It is the first mode by positivity and the spherical variational
principle. For uniqueness, the ground-state identity with this
positive eigenfunction \(V\) writes the shifted angular energy as
\(\int V^2|\nabla(h/V)|^2\) on compactly supported smooth tests.
Approximation in the Dirichlet form norm and lower semicontinuity
on compact interior subsets show that every equality case has
\(h/V\) constant. The ordered facet \(W\) is a proper subchamber
with an additional wall. Equality in domain monotonicity would
then give, by zero extension, a multiple of \(V\) vanishing on
an open set. Thus the lowest angular exponent exceeds the untied Vandermonde
degree used in the finite-correlation comparison.

The cylindrical assertion follows by multiplying both Gaussian
kernels by the Mehler kernel on \(\mathbb R^s\) and using the product
Gaussian. The diagonal, both kernel rows, and the ratio bound in
(SG3) stay the same on \(W\times\mathbb R^s\). Equations (SG5)–(SG7)
hold there with the full product Mehler form. Its Hermite tail and
the exhaustion by compact interior sets give the same compactness,
lower limit and recovery. At finite correlation the product
interaction is still Hilbert–Schmidt. In the limit the form is the
sum of the Dirichlet Gaussian gradient form in \(W\) and the
ordinary Gaussian gradient form in \(\mathbb R^s\), minus the norm.
Their eigenbasis is the tensor product. This gives the level shifts
and multiplicities stated above.

## The nonstandard two-row sectors

For \(k\ge4\) and \(2\le p\le\lfloor k/2\rfloor\), set
\(\lambda=(k-p,p)\). The active multiplicity space for this
nonstandard two-row shape is the space of
\(G_0=S_{\{3,\ldots,k\}}\)-equivariant \(V_\lambda^-\)-valued
functions on the full reference facet \(F\). The fiber
\(V_\lambda^-\) is the minus eigenspace of \((12)\) in the real
orthogonal Young representation. By vertical two-strip Pieri
it has \(G_0\)-type \((k-p-1,p-1)\). This type has no fixed
vectors, so all multiplicities satisfy balance and the
constant angular mode is absent.

Let \(n=k-2\), \(\sigma=\sqrt{1-\rho^2}\), and \(\phi_0=(2\pi)^{-1/2}\).
Let \(\Lambda_j^\lambda\) be the spherical Neumann eigenvalues
of \(F\cap S^{n-1}\) in the equivariant space, with multiplicities,
and let \(\beta_j^\lambda\) be the positive solutions of
\(\beta_j^\lambda(\beta_j^\lambda+n-2)=\Lambda_j^\lambda\).

**Theorem (Two-row limiting levels).** On each such shape, every
fixed variational level of the interface form multiplied by
\(\rho\sigma/[2\phi_0(1-\rho)]\) converges as \(\rho\uparrow1\).
The limits are the members of the ordered multiset
\[
          \{\beta_j^\lambda+2m-1: j\ge1,\ m\ge0\}.
\]
Each fixed finite list consists of discrete eigenvalues for
\(\rho\) sufficiently close to one. In tetrahedral shape
\((2,2)\) the multiset is
\[
 \left\{(j+\tfrac12)\frac{\pi}{\arctan\sqrt2}+2m-1:
                                      j,m\ge0\right\}.
\]

In (S85) of the [interface reduction](SIMPLEX_INTERFACE_SYMMETRY.md),
the adjacent triple-sign projections vanish identically on
two-row shapes. Thus on the full \(F\) the kernel is
\(K_{\mathrm{std},\rho}\) times the identity of the fiber.
The normalized quadratic form is precisely (SN2) of the
[standard-sector proof](STANDARD_SECTOR_NEAR_ONE.md), with
the Euclidean fiber norm in every square difference and with
equivariance in place of ordinary balance. The same
scalar compactness and lower-limit estimates (SN4)–(SN8) apply
componentwise in an orthonormal fiber basis. The equivariance
constraint is preserved by \(L^2\) convergence. For recovery,
apply the Sobolev approximation used for (SN3) to each component
and average over the joint action of \(G_0\) on the facet and
fiber. This preserves the desired equivariance. The finite-plane
min–max argument of (SN11) gives convergence of each fixed
level to the Gaussian Neumann problem on the same equivariant
space, shifted by minus one. The integral operator is
Hilbert–Schmidt, since the fiber is finite dimensional and the
scalar kernels are dominated by a constant times the Mehler
kernel. The local diagonal lower bound is
\(\gamma_H(F)/(1-\rho)-1\), by the same row estimate as in
(SN11). This proves the discreteness assertion.

In polar coordinates the same Laguerre separation as in
(SN11) gives the stated multiset. The absence of invariant
fiber vectors excludes angular eigenvalue zero. The first
spherical eigenvalue is at least the full first nonconstant
Neumann value, which is strictly above \(n-1\) for \(n\ge3\)
by the pointed-cone comparison in the standard-sector note.
This gives a positive first limiting gap. For \(n=2\) the
explicit tetrahedral formula has the same consequence.

In the ordered facet, a sorting wall for an untied transposition
has trace in the corresponding plus eigenspace of that
transposition. This follows from the gluing of the reflected
Sobolev pieces. It is a Dirichlet condition on the minus
component and the natural Neumann condition on the plus
component. At the top-tie wall the boundary condition is
Neumann on the full fiber. In tetrahedral shape \((2,2)\),
the fiber is the sign representation of \(G_0=S_2\).
The ordered angular interval of length \(\arctan\sqrt2\)
therefore has one Dirichlet and one Neumann endpoint. Its
angular frequencies are
\((j+\tfrac12)\pi/\arctan\sqrt2\), giving the exact formula.

In ambient dimension \(N\ge k-1\), tensoring with \(s=N-(k-1)\)
unused directions adds each nonnegative Hermite degree \(r\) to
the limiting levels. To verify fixed-index convergence in the
product, separate the inactive Hermite degree already at
finite \(\rho\). Let \(Q_\rho\) be the nonnegative part of (SN2)
before division by \(1-\rho\) and subtraction of the potential.
Degree \(r\) multiplies the interaction by \(\rho^r\), so the
normalized form on that summand is
\[
 \frac{\rho^r Q_\rho(h)}{1-\rho}
      +\left(\frac{1-\rho^r}{1-\rho}-1\right)
                                  \int_F b_\rho|h|^2d\gamma_H.
\]
For each fixed \(r\), the compactness and recovery argument
therefore gives the active limit plus \(r\). For \(r\ge L\ge1\),
the displayed form is bounded below by
\[
 \gamma_H(F)\left(\frac{1-\rho^L}{1-\rho}-1\right)\|h\|^2.
\]
As \(\rho\uparrow1\) this lower cap tends to
\((L-1)\gamma_H(F)\). Thus only finitely many inactive degrees
can contribute to each bounded range of fixed-index limits.
Their multiplicities are those of the ordinary \(s\)-dimensional
Hermite spaces.
The interaction in the full product remains Hilbert–Schmidt at each
fixed correlation, by the same Gaussian kernel bound on the product.

## The leading junction penalty for fixed interface fields

Let \(V_\lambda^-\) be the minus space in an irreducible real
orthogonal Young model of the \(k\)-cell simplex, \(k\ge3\).
On the full reference facet \(F\subset H\), let \(h\) be a
\(G_0\)-equivariant field with values in \(V_\lambda^-\), obtained
by restricting a smooth compactly supported field on \(H\).
Use the ordinary balance constraint in the standard sector.
Write \(F_l=\overline F\cap\{x_2=x_l\}\) for the top-tie wall and
\(H_l=H\cap\{x_2=x_l\}\) for its supporting hyperplane.
Let \(A_l\) be the orthogonal projection onto the sign representation
of the triple \(\{1,2,l\}\), as in (S85) of the
[interface reduction](SIMPLEX_INTERFACE_SYMMETRY.md).

Let \(\mathcal E_\rho\) be the interface form normalized by
(SN2) of the [standard-sector proof](STANDARD_SECTOR_NEAR_ONE.md),
now in the full fiber. In other words,
\(\mathcal E_\rho(h)/\|h\|_{L^2(F,\gamma_H)}^2\) is the interface
Rayleigh quotient times \(\rho\sigma/(2\phi_0(1-\rho))\), where
\(\sigma=\sqrt{1-\rho^2}\) and \(\phi_0=(2\pi)^{-1/2}\).

**Theorem (Fixed-field junction penalty).** For each such \(h\),
\[
 \lim_{\rho\uparrow1}\sigma \mathcal E_\rho(h)
       =\frac{1}{\sqrt3}\sum_{l=3}^k
                      \int_{F_l}|A_l h|^2\,d\gamma_{H_l}.
 \tag{JP1}
\]
If \(A_lh=0\) on all the top-tie walls \(F_l\), then
\[
 \lim_{\rho\uparrow1}\mathcal E_\rho(h)
                  =\int_F (|\nabla_H h|^2-|h|^2)\,d\gamma_H.
 \tag{JP2}
\]

Indeed (S85) makes the actual adjacent transport \(I-3A_l\).
The row identity in (SN2) gives
\[
 \mathcal E_\rho(h)
 =\frac{J_\rho(h)+\tfrac12\sum_lJ_{\rho,l}(h)}{1-\rho}
  -\int_F b_\rho|h|^2d\gamma_H
  +\frac{3}{2(1-\rho)}\sum_l
       \iint_{F^2}m_{\rho,l}(x,y)
             \langle A_lh(x),A_lh(y)\rangle\,d\gamma_H(x)d\gamma_H(y),
 \tag{JP3}
\]
where the difference energies use the fiber norm and the kernels
\(m_{\rho,l}=m_\rho e^{-\rho(w_l\cdot x)(w_l\cdot y)/\sigma^2}\).
The terms before the last sum tend to the right side of (JP2), by
the fixed-test proof (SN3) componentwise.

Fix \(l\) and write \(\eta=w_l/|w_l|\) so \(|w_l|^2=3/2\).
Use coordinates
\[
 x=z+\sigma u\eta,\qquad
 y=\rho z+\sigma\zeta+\sigma v\eta,
 \qquad z,\zeta\in H_l,\quad u,v\ge0.
\]
With respect to \(d\gamma_{H_l}(z)d\gamma_{H_l}(\zeta)du\,dv\),
the Gaussian weight of the adjacent integral divided by \(\sigma\)
is
\[
 \phi(\sigma u)\phi(v-\rho u)e^{-3\rho uv/2}.
\]
For fixed \(z,\zeta,u,v\), the two facet indicators converge to
\(1_{F_l}(z)\) outside the other wall intersections, which have
\(\gamma_{H_l}\)-measure zero. The two field values tend to
\(A_lh(z)\). Dominated convergence is valid because \(h\) is
bounded and
\[
 \phi(v-\rho u)e^{-3\rho uv/2}
       =(2\pi)^{-1/2}e^{-(v^2+\rho uv+\rho^2u^2)/2}.
\]
For \(\rho\ge1/2\) the right hand side is bounded on the nonnegative
quadrant by \((2\pi)^{-1/2}e^{-(v^2+u^2/4)/2}\).
The limiting normal integral equals
\[
 \frac1{2\pi}\int_{u,v\ge0}e^{-(u^2+v^2+uv)/2}\,du\,dv
 =\frac1{2\pi}\int_0^\infty\frac{ds}{s^2+s+1}
 =\frac1{3\sqrt3}.
\]
For the first equality evaluate in polar coordinates and put
\(s=\tan\theta\). The last sum in (JP3), multiplied by \(\sigma\),
now gives (JP1) since \(\sigma^2/(1-\rho)=1+\rho\).

If the projected trace vanishes, the smooth bound
\(|A_l h(x)|\le C(w_l\cdot x)\) holds on \(F\). To see it, use the
linear gap coordinates \(g_i=w_i\cdot x\), setting \(g_l=0\)
while keeping the other nonnegative gaps fixed. The resulting
point lies on \(F_l\) and is within \(C_k g_l\) of \(x\).
By symmetry and \(2ab\le a^2+b^2\), the absolute value of the
adjacent integral in (JP3) is at most
\[
 \int_F m_{\rho,l}^F(x)|A_lh(x)|^2d\gamma_H(x)
  \le C\int_F e^{-\rho^2g_l(x)^2/(4\sigma^2)}
                      g_l(x)^2d\gamma_H(x)=O(\sigma^3).
\]
The row bound follows by completing the Gaussian square as in
(S91d), and the one-variable estimate uses the bounded Gaussian
density of \(g_l\). The whole last term in (JP3) then tends to
zero, proving (JP2).

The fixed-field limits identify the wall constraint in each
symmetry type. Spectral convergence in a general shape also
requires a lower-limit estimate for fields varying with \(\rho\).
