# The reversible chain on the six-cell simplex facet

Only the ordinarily balanced standard sector remains in the six-cell
all-correlation interface problem. By
[(S49)–(S51)](SIMPLEX_INTERFACE_SYMMETRY.md#every-cell-count-the-standard-sector-is-a-scalar-positive-kernel),
it is a scalar problem on the reference facet
\(F=\{x_1=x_2\ge x_3,x_4,x_5,x_6\}\) in the pair plane
\(H=\{x\in\mathbf1^\perp:x_1=x_2\}\). Its functions are invariant
under permutation of the four untied coordinates.

## Ordinary balance in the stationary space

Write \(\mu=\gamma_H|_F/p_4\), \(p_4=\gamma_H(F)\).
Scale the standard kernel \(K_\rho\) and use \(D_\rho=K_\rho1\)
with respect to \(\mu\). A constant rescaling of both quantities
does not change the chain. Its reversible transition and probability
measure are
\[
 P_\rho(x,dy)=\frac{K_\rho(x,y)}{D_\rho(x)}\,d\mu(y),
 \qquad
 d\pi_\rho(x)=\frac{D_\rho(x)}{\tau_\rho}d\mu(x),\qquad
 \tau_\rho=\mathbb E_\mu D_\rho .
\]
The standard form is a positive scalar multiple of
\(\mathcal E_P(h)-(1-\rho)\|h\|_{\pi_\rho}^2\), with
\(\mathcal E_P(h)=\langle h,(I-P_\rho)h\rangle_{\pi_\rho}\) and
\(\mathbb E_\mu h=0\). Put \(s=\tau_\rho/D_\rho-1\). Then
\(\mathbb E_\pi s=0\), \(\|s\|_\pi^2=C_\rho-1\), with
\(C_\rho=\tau_\rho\mathbb E_\mu D_\rho^{-1}\).
Writing \(h=u+\mathbb E_\pi h\) with \(\mathbb E_\pi u=0\), ordinary
balance gives \(\mathbb E_\pi h=-\langle u,s\rangle_\pi\).
Consequently
\[
 \|h\|_\pi^2\le C_\rho\|u\|_\pi^2.
 \qquad
 \operatorname{gap}_{S_4}(P_\rho)>(1-\rho)C_\rho
   \ \Longrightarrow\ \text{a positive standard interface gap}.        \tag{MC1}
\]
The norm comparison is sharp for the mass constraint alone, by
taking \(u=s\) if \(s\ne0\). The actual interface form on the
stationary-mean-zero space is
\[
          \langle u,(\rho I-P_\rho)u\rangle_\pi
                          -(1-\rho)\langle u,s\rangle_\pi^2 .
\]
Thus a sharper exact test first requires \(\lambda_2(P_\rho)<\rho\)
on the invariant space, then evaluates
\((1-\rho)\langle s,(\rho I-P_\rho)^{-1}s\rangle_\pi<1\)
on the stationary-mean-zero subspace. This is (S70) in stationary
coordinates.

For the row scale \(3p_4\le D_\rho\le1\) from
[(SN2)](STANDARD_SECTOR_NEAR_ONE.md), one has
\[
                      C_\rho<\frac{117}{100}.                 \tag{MC2}
\]
Here is a uniform bound on the needed facet mass. Let \(p(t)\)
be the centered unit-variance four-dimensional equicorrelated
orthant probability at correlation \(t\in[0,1/3]\).
Differentiating its Gaussian density in \(t\) gives the sum
of the six mixed spatial derivatives with two distinct indices,
as seen from the Fourier transform. Integration over the
orthant puts the respective two variables at zero. The remaining
conditional pair has correlation \(t/(1+2t)\) by Gaussian
conditioning. Its orthant probability is the planar wedge fraction
\(1/4+\arcsin(t/(1+2t))/(2\pi)\). Hence
\[
 p_4=p(1/3)=\frac1{16}+\frac3{4\pi}\arcsin(1/3)
       +\frac3{2\pi^2}\int_0^{1/3}\frac{\arcsin(t/(1+2t))}
                                       {\sqrt{1-t^2}}\,dt .
\]
The integral is at least \(1/6-\log(5/3)/4\). Use
\(\arcsin(1/3)>55/162\) from its first two Taylor terms,
\(\pi<22/7\), and \(\log(5/3)<511/1000\) from the
odd power expansion of \(2\operatorname{arctanh}(1/4)\), with
the remainder after five terms bounded by a geometric series.
These give \(p_4>747/5000\) by rational arithmetic. On rescaling
\(D_\rho\) to the interval \([1,L]\), \(L=1/(3p_4)\), the
pointwise inequality \(D_\rho+L/D_\rho\le L+1\) gives
\(C_\rho\le(1+L)^2/(4L)\). Substitution of \(L<5000/2241\)
proves (MC2).

The calculation is the four-dimensional orthant specialization
of [Plackett's reduction formula](https://doi.org/10.1093/biomet/41.3-4.351).
The [directed standard-gap proof](SIX_CELL_WEIGHTED_STANDARD_GAP.md)
also gives a smaller mass bound at low correlation.

## Intermediate-correlation estimates

The numerical ratio in (MC1) exceeds one in all tests below.
The entries compare two discretizations with different integration
errors. They carry no directed error enclosure.

| \(\rho\) | \(C_\rho\), estimated | Nyström gap | Nyström ratio | Ritz ratio |
| ---: | ---: | ---: | ---: | ---: |
| 0.20 | 1.00405 | 0.94507 | 1.177 | 1.176 |
| 0.35 | 1.00919 | 0.86593 | 1.320 | 1.321 |
| 0.50 | 1.01397 | 0.74631 | 1.472 | 1.475 |
| 0.65 | 1.01737 | 0.58022 | 1.629 | 1.637 |
| 0.80 | 1.01825 | 0.36259 | 1.780 | 1.803 |
| 0.90 | 1.01605 | 0.17819 | 1.754 | 1.915 |

The larger runs are shown. At \(0.90\), the smaller Nyström
run gives \(1.545\) and a maximum holding entry \(0.803\);
the large run still has a maximum holding entry \(0.774\).
The high-correlation discretization is sensitive to poorly
resolved tails. The closer agreement at intermediate values
is evidence for testing a directed gap bound there.

Here is the shared coordinate law. Under \(\mu\) the untied gaps
\(g_i=x_1-x_{i+2}\) have the centered \(I+\mathbf1\mathbf1^T/2\)
Gaussian law conditioned to be positive. Sample the common
latent \(U\) with density proportional to \(e^{-u^2}\Phi(u)^4\),
then sample the four independent \(u-Y_i\) conditional on
\(Y_i\le u\), with \(Y_i\) standard normal. If
\(C_l=I-w_l w_l^T\) for \(w_l=(e_1+e_2)/2-e_l\), the adjacent
transition before restriction has mean \(\rho C_lx\), covariance
\(\sigma^2 I_H\) and total kernel weight
\(\exp[-\rho^2g_{l-2}^2/(4\sigma^2)]\). This follows from
\(|C_lx|^2=|x|^2-g_{l-2}^2/2\) by completion of the square.
The self term has mean \(\rho x\) and weight one.
The row \(D_\rho\) in the table is obtained by Gauss–Hermite
integration over the common coordinate of each of these
restricted transitions. The self row contributes once and the
four adjacent rows each contribute one half.

For the Nyström test, fold all permutations as in (S50).
The sampled facet law \(\mu\) is widened radially by a factor
\(1.12\), with exact density-ratio weights before quadrature
normalization. Divide the finite symmetric conductance by the
geometric mean of its two endpoint row sums. The root of the
row-weighted mass vector is the eigenvector at one. For the
Ritz test use monomials in the four power sums of
\(g_i/(1+g_i)\), through total degree five, with constants
removed by stationary centering. Scrambled Sobol integration
of the squared differences and symmetrized stationary
variance uses the self term and a single adjacent term with
weight two, since all four adjacent quadratic integrals are
equal for invariant test functions. Exact Ritz integrals
would bound the true gap from above. The numerical integration
adds a separate error.

The [probe source](probe_six_cell_reversible_chain.py) and
[recorded runs](six_cell_reversible_chain_probe.json) retain
both sample sizes and seeds. With NumPy 1.26.4 and SciPy 1.14.1,
the commands are

    OPENBLAS_NUM_THREADS=1 python research/plurality-stability/probe_six_cell_reversible_chain.py --mode ritz --power 16 --seed 221
    OPENBLAS_NUM_THREADS=1 python research/plurality-stability/probe_six_cell_reversible_chain.py --mode ritz --power 17 --seed 971
    OPENBLAS_NUM_THREADS=1 python research/plurality-stability/probe_six_cell_reversible_chain.py --mode nystrom --power 10 --seed 221 --dilation 1.12
    OPENBLAS_NUM_THREADS=1 python research/plurality-stability/probe_six_cell_reversible_chain.py --mode nystrom --power 11 --seed 971 --dilation 1.12

## Continuous projection and a directed tail bound

Finite-state path comparison requires transition demands and
edge congestion, as in [Diaconis–Saloff-Coste, Theorem 2.3](https://www2.stat.duke.edu/~scs/Courses/Stat376/Papers/ConvergeRates/ComparisonReversMC_DiaconisSalaffCoste_AAP1993.pdf).
The projection and restriction theorem of
[Jerrum–Son–Tetali–Vigoda, Theorem 1 and Corollary 3](https://webspace.maths.qmul.ac.uk/m.jerrum/papers/AAP021.pdf)
also uses a control on departures from a block. The empirical
Nyström matrices do not provide either comparison bound.

For an exact partition of the folded facet into finitely many
blocks of positive \(\pi\)-mass, let \(J\) be conditional
expectation on the blocks. Then \(JPJ\) is the true reversible
projection kernel, extended by zero on \(J^\perp\), since the
entries of its conductance are the double integrals of
\(\pi(dx)P(x,dy)\) over the block pairs. In terms of the symmetric
kernel \(k_\pi=\tau_\rho K_\rho/(D_\rho(x)D_\rho(y))\) relative to
\(\pi\), one has the exact continuous estimate
\[
 \begin{split}
 \lambda_2(P)&\le\lambda_2(JPJ)+\|P-JPJ\|_{\mathrm{HS}},\\
 \|P-JPJ\|_{\mathrm{HS}}^2
   &=\iint k_\pi^2\,d\pi^2-
     \sum_{ij}\frac{\left(\int_{B_i\times B_j}
                            k_\pi\,d\pi^2\right)^2}
                       {\pi(B_i)\pi(B_j)} .                  \tag{MC3}
 \end{split}
\]
The first inequality is the Rayleigh bound on \(1^\perp\).
Here \(\lambda_2(JPJ)\) includes the zero on the orthogonal
complement of the block space. The second identity is orthogonality
of the conditional expectation on pair blocks. Square integrability
holds at every fixed \(0<\rho<1\) by the Gaussian Mehler exponent
and the positive lower bound on \(D_\rho\).

A lower-dimensional tail estimate uses ambient Hermite degree.
Truncate the positive standard interaction to ambient degree \(N\)
before dividing by \(D_\rho^{1/2}\) at both endpoints. Denote the
resulting stationary-space finite-rank operator by \(P_{\le N}\).
Then
\[
             0\le P-P_{\le N}\le
                         \frac{\rho^{N+1}}{p_4} I .           \tag{MC4}
\]
Indeed (S11) gives positivity degree by degree, and each transported
edge term is bounded in absolute form by the same-facet term of
that degree. The weights in \(2M+\sum_{l=3}^6 M_l\) add to six.
For a planar Hermite coefficient of degree \(j\) of the zero
extension of a facet function, the same-facet omitted ambient
factor is
\(\sum_{r:\,j+2r>N}\rho^{j+2r}\binom{2r}{r}/4^r\).
It is at most \(\rho^{N+1}/\sigma\) by shifting the lower limit
and using decrease of the central-binomial factors. Multiplication
by \(\sigma p_4/2\) converts the ambient kernel to the chosen scale,
and converting the zero-extension square norm gives the
factor \(p_4^{-1}\) after the two normalized facet integrals.
The unnormalized tail is thus at most \(3\rho^{N+1} I\).
Conjugating by \(D_\rho^{-1/2}\), with \(D_\rho\ge3p_4\), proves
(MC4). In particular, a directed finite-rank eigenvalue bound
\(\lambda_2(P_{\le N})+\rho^{N+1}/p_4<1-(1-\rho)C_\rho\)
certifies (MC1) while keeping both kinds of balance.

## Limit of the sufficient ratio near one

The sharp ordinary-balanced limit in
[(SN1)–(SN2)](STANDARD_SECTOR_NEAR_ONE.md) also identifies the
stationary-gap ratio:
\[
   \lim_{\rho\uparrow1}
        \frac{\operatorname{gap}_{S_4}(P_\rho)}
             {(1-\rho)C_\rho}=2 .                            \tag{MC5}
\]
To see the lower bound, take a stationary-centered \(v\) and
set \(h=v-\mathbb E_\mu v\). The square-difference form is
unchanged, and \(\int D_\rho h^2d\mu\ge\int D_\rho v^2d\mu\).
Writing \(\mathcal J=J_\rho+\tfrac12\sum_l J_{\rho,l}\)
with all integrals normalized by \(p_4\), (SN1)–(SN2) give
\[
 \frac{\mathcal J(h)}{1-\rho}
 \ge \int D_\rho h^2d\mu +(1+o(1))\int h^2d\mu
 \ge (2+o(1))\int D_\rho v^2d\mu .
\]
The last step uses \(D_\rho\le1\). Dividing both stationary
quantities by \(\tau_\rho\) proves the lower bound on the
stationary gap. For the upper bound, apply the fixed-function
limit (SN3) to smooth invariant approximations of the radial
Neumann eigenfunction \(|x|^2-4\) and center each in \(\pi_\rho\).
Here \(D_\rho\to1\) on the facet interior and
\(3p_4\le D_\rho\le1\). This also gives \(C_\rho\to1\) by
dominated convergence.

## Axial transport to a Metropolis proposal on the A4 chamber

There is also an exact comparison that needs no nonlinear transport
of the sphere. Let \(C_A\) be the ordered \(A_4\) chamber, with
\(y_1\ge\cdots\ge y_5\) and \(\sum_i y_i=0\). Use the raw gaps
\(y_1-y_{j+1}\) as the four gaps of the native ordered facet. Here \(p=p_4\). This
defines a linear map \(y\mapsto Jy\) from the reference pair plane
to the native one. It fixes perpendicular lengths to the all-equal
gap axis and expands lengths on that axis by
\(a=\sqrt{5/3}\); this follows from the two gap metrics in [(SN9)](STANDARD_SECTOR_NEAR_ONE.md).
Write \(X(y)=\sqrt{5/4}\max_i y_i\), \(\delta=a^2-1=2/3\), and
\(c_0=a/(5p)\). If \(\mu_A\) is the conditional Gaussian probability
on \(C_A\), the pulled native ordinary measure is
\[
       d\mu_N(Jy)=c_0 e^{-\delta X(y)^2/2}\,d\mu_A(y).
\]
Indeed the native coordinates on positions \(2,\ldots,6\)
differ from the five reference coordinates by a common shift,
and the additional coordinate duplicates the native top coordinate.
For each permutation of those five positions the native bilinear
product equals the reference product plus \(\delta X(y)X(z)\).
Thus every permutation penalty is preserved. Hence the native subgroup part of the
folded conductance, after normalizing its average over the subgroup,
pulls back to
\[
    c_0^2
    \exp[-\delta(X(y)^2+X(z)^2-2\rho X(y)X(z))/(2\sigma^2)]
      M_A(y,z)\,d\mu_A(y)d\mu_A(z),
\]
where \(M_A\) is the folded \(A_4\) OU kernel. The full native average has weight \(5/6\) on this subgroup
and \(1/6\) on the 24 self terms. The self average is at least the
subgroup average, since \(m_{ho,l}\le m_ho\) on \(F^2\).
It is also at most five times the subgroup average, which already
contains all the self terms. Thus the full native conductance, in
the scale of the full average, lies between the subgroup conductance
and twice that conductance.

Relative to the ordinary stationary Metropolis chain with target
\(\mu_N\circ J\) and proposal \(M_A\), the ratio of that subgroup
conductance is exactly
\[
     c_0\exp[-\delta(\min(X(y),X(z))-
                      \rho\max(X(y),X(z)))^2/(2\sigma^2)].
\]
For example on diagonal rays this tends to zero as the radius
tends to infinity, even with the extra self terms included. By
continuity the ratio becomes small on open off-diagonal sets as
well. This rules out a uniform positive pointwise lower bound
with the unmodified Metropolis proposal.

A global lower comparison is possible for a shorter proposal.
Let \(Q_*(y,dz)\) be the proposal to take a Gaussian centered at
\(\rho y\) with covariance \((\sigma^2/a^2)I\), then sort its
coordinates. Make this proposal reversible with target
\(\mu_N\circ J\) by the usual Metropolis acceptance using the
two directed densities. Its conductance is the minimum of the
two directed products. The native conductance dominates it
pointwise by the factor \(c_0 a^{-4}\).
In fact, for \(X(y)\ge X(z)\), the max-axis function \(X\) is
homogeneous and 1-Lipschitz. For every permutation \(\pi\),
\[
 |X(z)-\rho X(y)|\le|\pi z-\rho y|.
\]
The extra Metropolis damping is therefore bounded below termwise
by \(\exp[-\delta|\pi z-\rho y|^2/(2\sigma^2)]\). This changes
the covariance to \(\sigma^2/a^2\) with prefactor \(a^{-4}\).
The directed product out of \(y\) becomes \(d(\mu_N\circ J)(y)
Q_*(y,dz)\), and its minimum with the reverse directed product
can only be smaller. Exchange \(y,z\) for the other axis order.
This is a genuine comparison of reversible energies, useful if
a quantitatively sufficient gap for the shorter proposal can
be established. It supplies no such lower gap by itself.
