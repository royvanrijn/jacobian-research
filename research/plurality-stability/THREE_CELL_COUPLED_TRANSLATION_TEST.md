# Coupling the two translations at a planar maximizer

Fix a balanced Gaussian three-cell maximizer \(A\) at
\(0<\rho<1\), in its active plane. The two translation modes
give a joint test for the normal field. On a compact theta,
it detects the derivative of cross-agreement when one of
the bounded cells rotates. On an arbitrary maximizer, it
gives a correlation derivative bound in terms of a positive
Mehler energy of the three pair measures and a residual
that vanishes exactly at the simplex.

Write \(d_e=e_i-e_j\), let \(N_e\) point from \(i\) to \(j\),
and let \(d\nu_e=\phi_2\,d\mathcal H^1\) be Gaussian
arclength. Use the tagged union of regular pair interfaces,
with the score-gradient weight \(a>0\) and signed kernel
\(K(x,y)=M_\rho(x,y)\,d_{e(x)}\cdot d_{e(y)}\).
Use the actual-maximizer identities (G1)–(G3) and (G17)
for \(Q^{(1)}=\int af^2\,d\nu-\rho\langle f,Kf\rangle\).
In particular the real kernel of \(Q^{(1)}\) is exactly
the translation space \(V=\{v\cdot N:v\in\mathbb R^2\}\).
The Mehler integrals below are finite by (G3) and the
single-pair bound following it.
[Translation identities and index](THREE_CELL_GLOBAL_OBSTRUCTIONS.md)

## Joint cut matrix and the cross-agreement derivative

Put \(s=N_1+\mathrm{i}N_2\), so \(|s|=1\) and
\(\rho Ks=as\). On any finite measurable partition into pieces \(D_j\),
form the Hermitian matrix
\[
 C_{jk}=\int_{D_j\times D_k}
          K(x,y)\overline{s(x)}s(y)\,d\nu(x)d\nu(y).
\]
Its row sums \(r_j=\int_{D_j}a\,d\nu/\rho\) are real by the
translation identity. For complex \(z_j\) the exact test is
\[
 \rho^{-1}Q^{(1)}\!\left(\sum_j z_j1_{D_j}s\right)
       =z^*[\operatorname{diag}(r)-C]z\ge0.                 \tag{JT1}
\]
The complex quadratic form is the sum of the forms on the
real and imaginary parts. In particular, the shared phase
and amplitude are a null direction. The usual one-direction
translation cut does not use the antisymmetric cross terms
between two different translation directions.

Suppose now that the entire network is a compact theta,
and take each of its three arms as one piece. Orient
the arms in the cyclic label order \(3\to1\), \(1\to2\),
\(2\to3\), with \(A_1,A_2\) the bounded cells. For \(j\ne k\)
put \(w_{jk}=\operatorname{Re}C_{jk}\). The real row sums force
\(\operatorname{Im}C_{12}=\operatorname{Im}C_{23}
=\operatorname{Im}C_{31}=b\). The principal determinant
on coordinates 1 and 2 gives
\[
        b^2<w_{12}w_{23}+w_{23}w_{31}+w_{31}w_{12}.
                                                               \tag{JT2}
\]
Indeed the diagonal entries there are \(w_{12}+w_{31}\)
and \(w_{12}+w_{23}\), and the off-diagonal entry is
\(-w_{12}-\mathrm{i}b\). Strictness follows from the kernel
in (G17). If a cut speed realizes equality in (JT1), then
it equals \(\alpha s+\beta\bar s\) for constants
\(\alpha,\beta\in\mathbb C\). Some arm is nonstraight,
since all three are distinct and share their endpoints.
On that arm \((z_j-\alpha)s=\beta\bar s\) gives
\(\beta=0\), because \(s^2\) is nonconstant. On every arm
then \(z_j=\alpha\). The restriction \(z_3=0\) is therefore
positive definite.

There is a bulk formula for the cyclic imaginary term.
Let \(R_t\) be counterclockwise rotation and set
\(F(t)=\langle1_{A_1},T_\rho1_{R_t A_2}\rangle\).
The two outward normal currents of the bounded cells are
the differences of the arm currents \(2-1\) and \(3-2\).
The three distinct label inner products in \(C\) are
\(-1\). Hence the raw Mehler pairing of the two outward
normal currents, using its determinant rather than dot
product, is \(-3b\). To compute it, write \(g_\rho\) for
the joint Lebesgue Gaussian density and \(\sigma^2=1-\rho^2\).
Integration by parts on the two bounded cells, using
\[
 (\partial_{x_1}\partial_{y_2}-\partial_{x_2}\partial_{y_1})g_\rho
       =\sigma^{-2}\det(x,y)g_\rho ,
\]
gives \(-3b=\sigma^{-2}\int_{A_1\times A_2}\det(x,y)g_\rho\).
The density identity follows by multiplying its two
linear logarithmic gradients. Changing variable by \(R_t\)
in the second cell gives
\[
       F'(0)=-\frac{\rho}{\sigma^2}
          \int_{A_1\times A_2}\det(x,y)g_\rho(x,y)\,dx\,dy
             =3\rho b.                                    \tag{JT3}
\]
Thus an oriented bulk moment must fit inside the real
cut margin in (JT2). The ordinary scalar test lacks this
term. If a stationary candidate with these arcs violates
the non-strict version of (JT2), (JT1) supplies a negative
first extra-coordinate test. If the translation flux is
invertible, it also supplies a negative balanced planar
test by subtracting the flux-matching translation, which
leaves \(Q^{(1)}\) unchanged and gives
\(\rho Q_0(g)=Q^{(1)}(g)-(1-\rho)\|g\|_a^2<0\).

## A lower bound using unoriented interface measures

Put \(q_e=d_e/\sqrt2\) in the centered two-dimensional label
plane \(\Lambda=\{\sum_i x_i=0\}\). Define
\[
 \begin{split}
 p_e&=\nu_e(\Sigma_e),\qquad
 m_{ef}=\iint M_\rho(x,y)\,d\nu_e(x)d\nu_f(y),\\
 P&=\sum_e p_e q_e q_e^{\mathsf T},\qquad
 C_{\rm vol}v=\sum_e d_e\int_{\Sigma_e}v\cdot N_e\,d\nu_e .
 \end{split}
\]
The last map is the negative first moment matrix of the cells
by Gaussian divergence. It is invertible onto \(\Lambda\) by
full active rank. Write \(L=\sqrt2 P C_{\rm vol}^{-\mathsf T}\),
and put
\[
 \mathcal E_\rho(A)=\frac12\sum_{e,f}m_{ef}+\frac32\sum_e m_{ee},
 \qquad
 \mathcal R_\rho(A)=\int a\,|q_{e(x)}-L N(x)|^2\,d\nu(x).
\]
Then
\[
            S'_\rho(A)\ \ge\
         \mathcal E_\rho(A)+\frac{1-\rho}{\rho}\mathcal R_\rho(A).
                                                               \tag{JT4}
\]
Here the derivative varies the correlation with the partition
held fixed.

For an orthonormal basis \(u_1,u_2\) of \(\Lambda\), use the
two trial speeds \(f_j=q_e\cdot u_j\). Their local weighted
energies sum to \(\int a\,d\nu=\rho S'_\rho(A)\).
Their interactions sum to \(\mathcal E_\rho(A)\), since
\((d_e\cdot d_f)(q_e\cdot q_f)\) is 2 on the diagonal
and \(1/2\) off it. Their volume fluxes are
\(\sqrt2 P u_j\), so the respective flux-matching translations
are \((L^{\mathsf T}u_j)\cdot N\). Subtracting these leaves
\(Q^{(1)}\) unchanged by (G2). The remainder \(g_j\) is
balanced and satisfies \(Q^{(1)}(g_j)\ge(1-\rho)\|g_j\|_a^2\)
by ordinary planar stability. Sum over \(j\) to obtain
(JT4). The bounded speeds and their mass corrections are
admissible by the weighted-space extension of (IB3).
[Volume-flux correction, IB3](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux)

In particular, if \(\nu=\sum_e\nu_e\), positivity of the
Mehler factorization after removing the constant and
linear Hermite projections gives the further lower bound
\[
 \begin{split}
 \mathcal E_\rho(A)\ \ge\
 &\tfrac12 P_\gamma(A)^2+\tfrac32\sum_e p_e^2\\
 &+\rho\left(\tfrac12\left|\int x\,d\nu\right|^2+
                  \tfrac32\sum_e\left|\int x\,d\nu_e\right|^2\right).
                                                               \tag{JT5}
 \end{split}
\]
For a positive measure \(\mu\) in this calculation, its
constant and linear projections after \(M_{\sqrt\rho}\)
smoothing are \(\mu(\mathbb R^2)\) and
\(\sqrt\rho\int x\,d\mu\), respectively. Exhaustion and
finiteness of \(m_{ef}\) justify this for the interface
measures. The first line improves the bare perimeter-square
bound when the three pair perimeters differ.

There is a version that compares the smoothed pair
perimeters with cancellation in the normal directions.
Set
\[
 h_e(x)=\int M_{\sqrt\rho}(x,y)\,d\nu_e(y),\qquad
 z_e(x)=\int M_{\sqrt\rho}(x,y)N_e(y)\,d\nu_e(y),
\]
orienting the pairs cyclically. Mehler factorization and the
translation identity give
\[
 S'_\rho(A)
 =\int\left(3\sum_e|z_e|^2-\left|\sum_e z_e\right|^2\right)d\gamma_2 .
\]
Indeed the integrand is the squared norm of
\(\sum_e d_e\otimes z_e\), and (G2) makes the integrated
translation interaction equal to \(\int a\,d\nu/\rho\).
Expressing (JT4) in these variables yields
\[
 \frac12\int\sum_{e<f}(h_e-h_f)^2\,d\gamma_2
 \ \ge\
 \int\left(3\sum_e(h_e^2-|z_e|^2)+
                        \left|\sum_e z_e\right|^2\right)d\gamma_2
      +\frac{1-\rho}{\rho}\mathcal R_\rho(A).                \tag{JT7}
\]
All terms on the right are nonnegative, by \( |z_e|\le h_e\).
Thus cancellation of normals within a pair and failure of
the three normal currents to sum to zero both require a
quantitative imbalance of the smoothed pair perimeters.

## Equality and a uniform positive residual

At the simplex the normals are isometric images of the
three cyclic \(q_e\), so each \(f_j\) is itself a translation
and \(S'_\rho=\mathcal E_\rho\).
For the converse, if \(\mathcal R_\rho=0\), complete adjacency
gives \(N_e=L^{-1}q_e\) almost everywhere on every pair.
These are three fixed unit normals summing to zero, so
their cyclic pairwise angles are \(120^\circ\).

For a cyclic pair \(e=(i,j)\), the distributional derivative
of \(1_{A_i}-1_{A_j}\) in the constant direction \(N_e\)
is a negative surface measure. With Euclidean arclength
its coefficients are \(-2\) on \(e\) and \(-1/2\) on each
other pair, by the cyclic inward and outward signs and
the \(120^\circ\) dot products. Gaussian smoothing makes
this derivative strictly negative everywhere. On any active
pair arc the derivative is normal to the score equality
and the constant normal makes the arc straight. Analyticity
extends the pair level along that entire affine line.
Strict normal monotonicity makes that line the full pair
equality set. The three lines concur since the three
cyclic adjusted pair contrasts add to zero. Thus the
winner cells are an affine \(120^\circ\) fan, which
is centered by balance and the affine fan rigidity.
[Affine fan rigidity](AFFINE_SIMPLEX_FAN_RIGIDITY.md)

In fact on every compact correlation band \(J\subset(0,1)\)
the residual is uniformly positive among nonsimplex balanced
planar maximizers, if there are any. Suppose otherwise, with
\(A^n,\rho_n\) a sequence where \(\mathcal R_{\rho_n}\to0\).
Strong maximizer compactness gives a maximizer limit \(A\),
with a convergent correlation subsequence. The first moment
matrices, and hence \(C_{{\rm vol},n}\), converge to an
invertible matrix. The perimeters are bounded by polynomial
local length growth and Gaussian tails. At least two pair
perimeters stay bounded below along a subsequence. Otherwise
the representation of \(C_{{\rm vol},n}\) as a sum of three
rank-one maps, each of norm at most \(\sqrt2 p_{e,n}\),
would make its limit rank at most one. Consequently
\(L_n\to L_*\) after extraction, and \(L_*\) is invertible,
since at least two nonparallel \(q_e\) contribute to the
limit of \(P_n\).

On each fixed ball, (G20) and Cauchy–Schwarz give
\[
 \int |q_e-L_nN_e^n|\,d\nu_e^n
 \le \mathcal R_{\rho_n}^{1/2}
         \left(\int a_{e,n}^{-1}\,d\nu_e^n\right)^{1/2}=o(1).
\]
The integrals here are restricted to that ball. The error
also tends to zero locally with Euclidean arclength. Extract
vague limits \(\lambda_e\) of the three Euclidean arclength
measures. The distributional gradient of the indicator
vector \(1_{A^n}\) converges to that of \(1_A\). After
multiplying its physical index by \(L_n\), the displayed
estimate gives the matrix-valued measure identity
\[
         -D1_A\,L_*^{\mathsf T}/\sqrt2
                          =\sum_e q_e q_e^{\mathsf T}\lambda_e .
                                                               \tag{JT6}
\]
On a regular active pair \(f\) of the limit, the left side
has density \(q_f (L_*N_f)^{\mathsf T}\) relative to arclength.
It is positive semidefinite by (JT6) and has rank one by
invertibility. Therefore \(L_*N_f=tq_f\) for \(t>0\), and
only \(\lambda_f\) can contribute there, since the other
\(q_e\) are nonparallel. In particular \(\lambda_f\) has
positive local mass for each actual pair by complete
adjacency. On any such ball the earlier convergence of
the normals to \(L_*^{-1}q_f\) in the integral against
the approximating measures then implies
\(|L_*^{-1}q_f|=1\). Thus \(t=1\) on each regular pair.
The limit has the constant-normal structure just proved
to force the simplex. This contradicts uniform measurable
local simplex isolation on \(J\).
[Maximizer compactness and measurable isolation](POSITIVE_CORRELATION_CONTINUATION.md)

For a tangential nonsimplex tie, (JT4) now forces
\(\mathcal E_\rho(A)<\mathcal E_\rho(T)\), with a uniform
gap on compact bands. A reverse comparison for the
positive surface interaction at an unrestricted maximizer
would exclude such ties. The estimates above do not supply
that comparison.
