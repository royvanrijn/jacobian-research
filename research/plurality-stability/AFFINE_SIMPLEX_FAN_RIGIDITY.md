# Rigidity of stationary affine Gaussian fans

Fix \(k\ge3\) and \(0<\rho<1\). Let \(A_1,\ldots,A_k\) be the
cells of the affine winning rule
\[
        x\longmapsto \operatorname*{arg\,max}_{i}
                    \big(\langle u_i,x\rangle-b_i\big),
\]
where the sites \(u_i\) are affinely independent in \(\mathbb R^d\).
Suppose that for some constants \(\beta_i\) the adjusted Gaussian
scores \(T_\rho 1_{A_i}-\beta_i\) are equal on every shared facet.
Then the partition is a central regular simplex fan, with cylindrical
extension in the orthogonal complement of its site-difference span.
In particular, the equality of the Gaussian cell masses follows from
the affine geometry and the stationary score equations.

This excludes a full-rank affine fan as a nonsimplex tied maximizer for
any number of cells. The planar result in the
[global obstruction note](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#convex-and-conical-maximizers)
also treats degenerate convex fans. The centering step below uses the
opposite-ray mechanism of Heilman–Mossel–Neeman, who prove centering
and equal masses under their flat-partition hypothesis. The facet
reflection and site argument then force regularity in the full-rank
affine class. [Heilman–Mossel–Neeman, Lemma 4.1 and Corollaries 5.1–5.2](https://arxiv.org/pdf/1403.0885)

## Centering from the stationary facet equations

Subtract the mean site from all sites and work in
\(E=\operatorname{span}\{u_i-u_j\}\). The indicators and their Gaussian
scores are cylindrical over \(E\), whose dimension is \(k-1\).
The affine equations
\(\langle u_i-u_k,p\rangle=b_i-b_k\) for \(i<k\) have a unique
solution \(p\in E\). Thus the cells are cones with common apex \(p\).
Every pair shares a facet. Indeed for fixed \(i\), the \(k-1\)
linear forms with normals \(u_i-u_\ell\), \(\ell\ne i\), are
independent. In their coordinates the cone of label \(i\) is the
positive orthant and its \(\ell\)-facet is shared with label \(\ell\).

Take a facet between \(i,j\) and set
\[
    n=\frac{u_i-u_j}{|u_i-u_j|},\qquad
    H=\{x\in E:\langle n,x\rangle=\langle n,p\rangle\}.
\]
The score \(U=T_\rho(1_{A_i}-1_{A_j})\) is real analytic on \(E\).
In the Gaussian convolution formula, extending the center to complex
coordinates gives a locally dominated integrand on every compact set,
using a Gaussian in the integration variable times an exponential of
at most linear growth. Restriction to \(H\) is therefore real analytic.
Its equality to \(\beta_i-\beta_j\) on the relative interior of the
facet extends to all of \(H\).

Choose \(\tau\in n^\perp\) pointing into that relative facet interior
from \(p\). Then
\(\langle u_i,\tau\rangle=\langle u_j,\tau\rangle>
\langle u_\ell,\tau\rangle\) for \(\ell\notin\{i,j\}\).
Write \(\sigma=\sqrt{1-\rho^2}\). At the observation \(p+s\tau\),
the conditional sample is \(\rho(p+s\tau)+\sigma Z\), with \(Z\)
standard Gaussian on \(E\). As \(s\to+\infty\), for every fixed \(Z\),
both \(i,j\) eventually have greater affine scores than the other
labels. Their comparison uses the sign of
\((\rho-1)\langle n,p\rangle+\sigma\langle n,Z\rangle\).
As \(s\to-\infty\) both labels eventually lose to any fixed other
label. Dominated convergence gives the two limits
\[
      2\Phi\!\left(\frac{(\rho-1)\langle n,p\rangle}{\sigma}\right)-1
      \qquad\text{and}\qquad 0
\]
for \(U(p+s\tau)\). Constancy on \(H\) forces
\(\langle n,p\rangle=0\) and \(\beta_i-\beta_j=0\).
All pair normals span \(E\), so \(p=0\).

## From score equality to reflection of the sites

For a fixed pair let \(f=1_{A_i}-1_{A_j}\) and let \(H=n^\perp\).
Use the homogeneous cone representatives and write \(y=(z,t)\) with
\(z\in H\) and \(t\) the normal coordinate. The equality
\(T_\rho f=0\) on \(H\) says that the tangential Gaussian convolution
of
\[
                g(z)=\int_{\mathbb R} f(z,t)\phi_\sigma(t)\,dt
\]
vanishes. Its argument \(\rho x\) ranges over \(H\). For bounded
\(g\), injectivity holds in tempered distributions. The Fourier
multiplier of this convolution is a nonvanishing Gaussian, and
division on compactly supported smooth frequency tests shows
\(\widehat g=0\).

For almost every direction \(\theta\) in \(H\) this gives
\(g(r\theta)=0\) for almost every \(r>0\). Homogeneity and the
substitution \(t=ru\) yield
\[
   \int_0^\infty [f(\theta,u)+f(\theta,-u)]e^{-a u^2}\,du=0,
             \qquad a=r^2/(2\sigma^2).
\]
The integral is continuous for \(a>0\), so the equality holds
throughout that interval. Substituting \(w=u^2\) gives a Laplace
transform with a locally integrable density bounded by
\(1/\sqrt w\). It is holomorphic on \(\operatorname{Re}a>0\).
Continuing the zero identity to \(a=1+i\xi\) and using Fourier
uniqueness on the density multiplied by \(e^{-w}\) show that
\(f\) is odd across \(H\) almost everywhere. The reflection
\(R_{ij}\) in \(H\) therefore exchanges \(A_i,A_j\). This is an
equality of closed cones, since each cone is the closure of its
interior and equality almost everywhere identifies those closures.

For \(\ell\notin\{i,j\}\) the reflection sends the inward normal
\(u_i-u_\ell\) of an \(i\)-facet to a positive multiple of some
normal \(u_j-u_{\ell'}\) with \(\ell'\notin\{i,j\}\). Project to
\(H\). The vectors
\[
     \operatorname{proj}_H(u_i-u_\ell)
       =\operatorname{proj}_H(u_j-u_\ell),
       \qquad \ell\notin\{i,j\},
\]
form a basis of \(H\), because the independent list
\(\{u_i-u_\ell:\ell\ne i\}\) contains exactly the normal used
to quotient out the perpendicular line. Thus the projected
equality forces \(\ell'=\ell\) and multiplier one. The remaining
normal \(u_i-u_j\) goes to its negative. With \(\pi=(ij)\), these
identities say
\[
          R_{ij}(u_i-u_\ell)=u_j-u_{\pi(\ell)}
                  \qquad\text{for every }\ell .
\]
Sum over \(\ell\) and use \(\sum_\ell u_\ell=0\) to obtain
\(R_{ij}u_i=u_j\) and then \(R_{ij}u_\ell=u_{\pi(\ell)}\) for
every \(\ell\). Every site transposition is realized by an
orthogonal transformation, so all pair distances coincide.
The common apex is zero and the resulting sites form a regular
simplex.

## Harmonic Taylor terms for general conical winners

Suppose instead that a measurable partition of positive cell masses
is dilation invariant about zero and its cells agree almost everywhere
with winners of \(T_\rho 1_{A_i}-\beta_i\). Define
\[
      V_i(x)=\mathbb E[1_{A_i}(x+Z)],
          \qquad m_i=\int_{A_i}x\,d\gamma_d(x).
\]
The identity \(T_\rho 1_{A_i}(x)=V_i(\rho x/\sigma)\) follows by
homogeneity. Each cell has positive angular measure and accumulates
at zero along score-winning points. Continuity therefore makes all
adjusted scores equal there.

The heat extension of \(1_{A_i}\) at variance \(t\) is
\(V_i(x/\sqrt t)\). Differentiating the heat equation at \(t=1\)
gives \(\Delta V_i+x\cdot\nabla V_i=0\). If \(P_j\) denotes the
degree-\(j\) Taylor term of an adjusted pair contrast, this means
\(\Delta P_{j+2}+jP_j=0\). The constant term is zero. The contrast
is nonzero by Gaussian injectivity and the positive masses of
its two disjoint cells. Its first nonzero homogeneous Taylor
term is consequently harmonic.

For almost every direction none of these finitely many leading
pair terms vanish. On sufficiently small positive radii along that
direction their signs identify the score winner. By the assumed
dilation invariance and polar Fubini this identifies the cell for
almost every point on the ray. Thus the cells have representatives
defined by the simultaneous polynomial inequalities of these
homogeneous harmonic terms. For a pair with different first
moments the leading term is the linear contrast
\(\langle m_i-m_j,x\rangle\). In particular, if
\(\operatorname{span}\{m_i:1\le i\le k\}\) has dimension \(k-1\),
then \(\sum_i m_i=0\) implies affine independence, and the partition
is exactly the affine-score fan with these sites. The theorem above
forces the central regular simplex.

The harmonic reduction isolates higher-degree terms only inside
groups of equal first moments. This covers arbitrary cell counts
under the full-rank hypothesis and identifies which conical
candidates still require an argument. Establishing conicality for
a general tied maximizer remains a separate obligation.
[Heilman, Lemmas 8.2–8.3 and Remark 8.4, on the separate hyperstable problem](https://arxiv.org/html/2209.11216#S8)

## Conical stationarity forces Euclidean minimal interfaces

In an active Euclidean space of dimension at least two, suppose
a stationary winning-score partition is invariant under positive
dilations about \(p\). On every regular \(C^3\) two-sided pair
interface patch the Euclidean mean curvature vanishes. Where \(p\)
has a nonzero component in the interface normal, the patch is flat.

Take \(x=p+w\) on a regular \(i/j\) patch and let \(n\) be its
unit normal pointing into \(A_i\). Define \(f=1_{A_i}-1_{A_j}\)
and \(\eta=(\rho-1)p/\sigma\). By homogeneity about \(p\),
\[
 T_\rho f(p+r w)=\mathbb E f\bigl(x+\varepsilon(Z+\eta)\bigr),
                  \qquad \varepsilon=\frac{\sigma}{\rho r}.
\]
The left side is \(\beta_i-\beta_j\) at every \(r>0\) by the
stationary score equation on the dilated patch. Thus the right
side is constant for all positive smoothing scales
\(\varepsilon\).

Use tangent coordinates at \(x\), with \(A_i\) locally above
the graph \(y_n=\psi(y_T)\) and \(A_j\) below. Here
\(\psi(0)=0\) and \(\nabla\psi(0)=0\); set \(B=D^2\psi(0)\).
Extend \(\psi\) smoothly outside the local chart with bounded
Hessian. The change this makes in the Gaussian expectation is
at most an exponentially small tail in \(\varepsilon^{-2}\).
Putting \(Y_T=Z_T+\eta_T\) and \(\eta_n=\langle\eta,n\rangle\),
the comparison of the two labels has expectation
\[
    2\,\mathbb E\Phi\!\left(\eta_n-
                 \frac{\psi(\varepsilon Y_T)}{\varepsilon}\right)-1
\]
up to that tail. In particular,
\(\beta_i-\beta_j=2\Phi(\eta_n)-1\).
After subtracting this limit and dividing by \(\varepsilon\),
dominated convergence gives
\[
     -\phi_1(\eta_n)\,
         \bigl(\operatorname{tr}B+\eta_T^{\mathsf T}B\eta_T\bigr),
\]
where \(\phi_1\) is the standard one-dimensional Gaussian
density. The domination uses
\(|\psi(\varepsilon Y_T)|/\varepsilon^2\le C|Y_T|^2\)
and boundedness of \(\phi_1\). The displayed limit must vanish.

The identity for \(\beta_i-\beta_j\) makes \(p\cdot n\)
constant along the patch. Its tangential derivative gives
\(B p_T=0\), so the term involving \(\eta_T\) vanishes and
\(\operatorname{tr}B=0\). On a smooth Euclidean minimal
hypersurface the normal obeys
\(\Delta_\Gamma n=-|B|^2n\). To see this, differentiate the
unit normal twice in a tangent orthonormal frame geodesic at
the point. The normal component of its Laplacian is
\(-|\nabla_\Gamma n|^2n\). The tangential component is the
gradient of the signed mean curvature by Codazzi and hence
vanishes. Taking the dot product with \(p\), constancy of
\(p\cdot n\) gives \(|B|^2(p\cdot n)=0\), which proves flatness
where \(p\cdot n\ne0\).

For an origin-conical partition, the argument gives equal
prices across every regular pair interface. If the regular
adjacency graph is connected, all prices and all Gaussian
masses coincide, since all adjusted scores tie at zero and
\(T_\rho1_{A_i}(0)=\gamma_d(A_i)\).

The small-scale identity also gives equal spherical areas for the
two labels at almost every radius about each regular interface
point \(x\) of a centered cone. Write \(g_x(r)\) for the average
of \(f=1_{A_i}-1_{A_j}\) on the Euclidean sphere of radius \(r\)
centered at \(x\). Since \(\mathbb E f(x+\varepsilon Z)=0\)
for every \(\varepsilon>0\), integration in polar coordinates
gives
\[
            \int_0^\infty g_x(r)r^{n-1}e^{-a r^2}\,dr=0
                         \qquad(a>0),
\]
where \(n\) is the active dimension. The substitution \(w=r^2\)
and Laplace uniqueness as above give \(g_x=0\) almost everywhere.
This is a pairwise instance of the relation between stationary
heat levels and uniform density in balls. The ball density is
one half on small balls meeting only the two labels.
[Magnanini–Prajapat–Sakaguchi, Theorem 1.1 and its proof](https://doi.org/10.1090/S0002-9947-06-04145-6)

On an open regular patch in dimension \(n\ge4\), the same
half-density identity forces \(K_3=0\), where \(K_3\) is the
third elementary symmetric polynomial of the principal curvatures.
Shrink the patch to a fixed two-label chart so that a common
range of small balls lies in this chart. In this range the fraction
of each ball in \(A_i\) is \(1/2\), by integrating the spherical
identity. Apply the local density expansion in
Magnanini–Prajapat–Sakaguchi, Theorem 3.5, to the upper side
of the chart. Its cubic coefficient is a nonzero multiple of
\(K_1^3-4K_1K_2+4K_3\), and \(K_1=0\) by the curvature argument
above. Thus the sum of the cubes of the principal curvatures
also vanishes by Newton's identity.
[Magnanini–Prajapat–Sakaguchi, Theorem 3.5 and equation (3.5)](https://doi.org/10.1090/S0002-9947-06-04145-6)

For an origin-conical partition, at an ordinary triple junction
with three tangent wedges,
continuity gives equality of the three scores there. The same
small-scale identity and Gaussian convergence to the three
tangent wedges give equal angular masses in the two-dimensional
normal plane. Thus the three angles are \(2\pi/3\).

In active dimension three this also restricts the number of
cells if the spherical interfaces form a finite network with
only ordinary edges and regular triple vertices. On
\(\mathbb S^2\) every regular edge is geodesic: writing its
cone as \(r\theta(s)\) with spherical arclength \(s\), the mean
curvature of the cone is the geodesic curvature divided by
\(r\). Each cell's spherical area is then an integer multiple
of \(\pi/3\) by Gauss–Bonnet, because all corner turning angles
are \(\pi/3\) and the edge curvature is zero. In the balanced
case the area is \(4\pi/k\), so \(k\) divides \(12\) within this
regular network class.

Nonflat stationary winners do occur in the higher-jet class.
Split \(\mathbb R^{km}\), \(m\ge2\), into \(k\) mutually
orthogonal blocks of dimension \(m\). Assign the label of the
block with largest Euclidean norm. The cells are balanced by
block symmetry, and all first moments vanish by sign symmetry.
The Gaussian score order agrees with the block norm order.
For a conditional noisy block, its norm increases strictly in
stochastic order with the original block norm. To verify this,
orient the first coordinate along the original block. Given the
other noise coordinates, the event that the noisy norm is below
any positive threshold asks for a shifted one-dimensional
Gaussian to lie in a centered interval. Its probability
decreases as the positive shift increases. Starting with two
equal block norms, increasing one therefore increases its
winning probability and decreases the other's. The score rule
holds with equal prices. The ordinary interfaces
\(|x_i|=|x_j|>\max_{\ell\notin\{i,j\}}|x_\ell|\) are curved
quadratic cones.

## Cylindrical stability bounds the adjacent Taylor degrees

Work in an active space of dimension \(n\ge2\). Suppose the
cylindrical extension by one independent Gaussian coordinate
satisfies the constrained second-variation inequality for compact
smooth speeds on the regular interfaces. Take a nonempty smooth
angular patch of an ordinary \(i/j\) interface. Its points at
radii \(0<r<r_0\) give a conical patch \(\Gamma\) on which only
the two labels meet. If the \(i/j\) adjusted score contrast first
has a nonzero term of degree \(m\), then
\[
                              m\le n-1.
\]

Let \(a\ge0\) be the derivative of the pair score difference in
the unit normal pointing into \(A_i\). Its sign follows from the
winning rule. Write \(d\nu\) for ordinary Gaussian surface measure
on \(\Gamma\) and \(M_\rho\) for the Gaussian-relative Mehler
kernel in \(E\). Apply the second-variation formula on the cylinder
to a speed \(q(x)\zeta(s)\) along the patch, with \(q\) smooth and
compactly supported in \(\Gamma\) and \(\zeta\) smooth, compactly
supported and odd in the new coordinate \(s\). Its first-order
volume change is zero. The local term factors with
\(\|\zeta\|_2^2\), while the interaction factors with
\(\langle \zeta,T_\rho\zeta\rangle\) and the squared incidence
norm \(2\) of the \(i/j\) pair. Approximating \(s\) in Gaussian
\(L^2\) by these odd cutoffs gives
\[
  \int_\Gamma a q^2\,d\nu
      \ge 2\rho\iint_{\Gamma^2}M_\rho(x,y)q(x)q(y)
                         \,d\nu(x)d\nu(y).
\]
The needed second-variation normalization is the same as in
the [translation and exterior-coordinate forms](THREE_CELL_GLOBAL_OBSTRUCTIONS.md).

Set \(c=\min_{\overline\Gamma^2}M_\rho>0\), using boundedness of the
patch. For \(q\ge0\) the right side is at least
\(2\rho c(\int_\Gamma q\,d\nu)^2\). Use nonnegative smooth
cutoffs of \((a+\epsilon)^{-1}\), bounded above by that function.
Their local integral satisfies \(\int_\Gamma a q^2\,d\nu\le
\int_\Gamma q\,d\nu\). Exhausting the patch and then taking
\(\epsilon\downarrow0\) proves
\(\int_\Gamma a^{-1}d\nu\le1/(2\rho c)\).

The contrast Taylor expansion gives \(a=O(r^{m-1})\) uniformly on
the patch. Its surface density contains \(r^{n-2}\,dr\) times a
positive angular surface measure and a Gaussian factor bounded
below near zero. The reciprocal integral thus dominates a positive
constant times \(\int_0^{r_0} r^{n-m-1}\,dr\), which diverges for
\(m\ge n\). In active dimension two this forces adjacent linear
contrasts. In larger dimensions it leaves only degrees below the
active dimension to be controlled on such patches.

### Planar origin-conical classification with cylindrical stability

In an active plane the preceding degree bound classifies the conical
winner partition under the same extra-coordinate stability hypothesis.
Group labels by their first-moment vectors. For a given group, the
directions where its vector strictly beats all other group vectors
form a connected open arc on the circle, or the whole circle if
there is only one group. If some group has positive-mass cells, it
wins on such a set of directions, since ties of distinct linear forms
are null on the circle. The leading-term description gives a finite
partition of its open directions into label intervals. If the group
has two or more labels, some pair of its labels share an internal
boundary ray. On a conical patch of this ray the first-order contrast
vanishes, contradicting the degree bound \(m\le1\). Thus all moment
vectors are distinct. Each label is an open sector of their linear
winning fan, with ordinary adjacent boundary rays.

Stationary score equality extends from each shared ray to its line
by analyticity. Tangential Gaussian deconvolution and homogeneity,
as above, give reflection interchanging the two neighboring sectors.
Their opening angles are therefore equal. Going around the circle
gives the regular \(k\)-sector fan. For \(k=3\) this is the regular
simplex. The argument assumes conicality and stability under
extra-coordinate variations, and applies with any number of
positive-mass cells in the active plane.

## Dilation of strong limits at infinity

An ordinary adjusted-score winner partition \(A\) with positive
cell masses need not be assumed conical for the following
implication. If \(1_{A_i}(R_j\,\cdot)\to 1_{B_i}\) in
\(L^1_{\mathrm{loc}}\) along some \(R_j\to\infty\) for every label,
then \(B(\rho x)=B(x)\) for almost every \(x\). If the strong
subsequential limit at infinity exists and is unique, it is invariant
under every positive dilation.

Indeed, with \(f_{j,i}=1_{A_i}(R_j\,\cdot)\), the conditional Gaussian
formula gives
\[
             T_\rho 1_{A_i}(R_j x)
              = (f_{j,i}*\phi_{\sigma/R_j})(\rho x).
\]
These smoothed functions converge to \(1_{B_i}(\rho x)\) in local
\(L^1\). To see this on a compact ball, use local \(L^1\) convergence
on a ball enlarged by one and bound contributions from outside it
by the Gaussian tail. The convolution of the fixed limiting
indicator tends locally to itself by the approximate identity.
For every ordered pair of labels \(\beta_i-\beta_\ell<1\), since
label \(i\) wins at some point and both conditional Gaussian
probabilities there lie strictly between zero and one. The limiting
adjusted scores at \(x\) therefore have the unique winner \(B(\rho x)\),
with gap bounded below by
\(\min_{i,\ell}(1-\beta_i+\beta_\ell)>0\).
Convergence of the scores in local \(L^1\) makes the measure of
incorrect winners on each compact ball tend to zero. The same
rescaled indicators already converge to \(B(x)\), proving
the discrete invariance. For any fixed \(t>0\), the sequence
\(tR_j\) has limit \(B(t\,\cdot)\) in the same topology. Uniqueness
gives continuous dilation invariance.

This identifies a consequence of strong compactness and uniqueness
at infinity. Using it to classify an original maximizing partition
would also require a mechanism transferring geometric information
back from those limits.

## Rigidity for a fixed partition at accumulating correlations

Let \(A\) be a fixed balanced \(k\)-cell partition whose Gaussian first
moments \(m_i=\int_{A_i}x\,d\gamma_d(x)\) span a space of dimension
\(k-1\). Suppose it is an adjusted-score winner at every correlation
in a set \(E\subset(0,1)\) with an accumulation point in \([0,1)\).
Then \(A\) is the central regular simplex. In particular, a fixed
full-rank nonsimplex has only finitely many globally maximizing
correlations in each \([0,b]\) with \(b<1\).

Write \(C_i\) for the closed essential support of \(1_{A_i}\). These
sets are independent of the correlation and cover the ambient space.
At any \(r\in E\), continuity of the adjusted scores gives equality of
the \(i,j\) scores on \(C_i\cap C_j\). The graph of nonempty support
intersections is connected, since otherwise two disjoint nonempty
finite unions of closed sets would separate Euclidean space. Choose
one point in each intersection along a spanning tree of this graph.
With the price of label \(k\) normalized to zero, the price differences
along the tree are the unadjusted score differences at those points.
The same formulas define prices \(\beta_i(r)\) for \(-1<r<1\).

These prices are real analytic in \(r\). At each real \(r_0\in(-1,1)\)
the conditional Gaussian kernel has a holomorphic extension in a
complex neighborhood with modulus bounded by a constant times
\(\exp(-c|y|^2+C|y|)\) for \(c>0\) at a fixed input \(x\). Thus every
support-intersection equality, initially given on \(E\), holds on
\((-1,1)\) by the identity theorem. Differentiating at zero and putting
\(b_i=\beta_i'(0)\) gives
\[
 C_i\cap C_j\subseteq
 \{x:\langle m_i-m_j,x\rangle=b_i-b_j\}.
\]
The differentiation uses \(\partial_r T_r1_{A_i}(x)|_{r=0}=
\langle m_i,x\rangle\). The \(k-1\) differences from \(m_k\) form
a basis of the moment span because \(\sum_i m_i=0\) and the rank is
\(k-1\). Hence all the displayed hyperplanes contain one point \(p\)
in that span. The connected chambers of the finite arrangement are
each contained in one support, since the closed supports cover a
chamber and do not intersect there. The chamber labels therefore
represent \(A\) almost everywhere. In particular \(A\) is a finite
polyhedral partition conical about \(p\), possibly with nonconvex
cells.

On an open facet where the labels change from \(j\) to \(i\), let \(n\)
be the normal pointing into \(A_i\). The conical smoothing identity
above holds for every \(r\in(0,1)\), since the score equality on that
facet and its positive dilations about \(p\) has been continued in
\(r\). Its small-scale limit gives
\[
 \beta_i(r)-\beta_j(r)
   = 2\Phi\!\left(\frac{(r-1)\langle p,n\rangle}{\sqrt{1-r^2}}\right)-1.
\]
At \(r=0\) the price difference is zero, since all cell masses are
equal. Taking \(r\downarrow0\) shows \(\langle p,n\rangle=0\) for
each label-changing facet. Every such facet thereby lies in a
hyperplane through the origin. The radial distributional derivative
of each indicator is zero, since the boundary term is supported on
these facets and weighted by \(x\cdot n=0\). The partition is therefore
conical about the origin. At any \(r\in E\) its harmonic leading terms
are the linear first-moment contrasts. The full-rank affine theorem
then gives the regular simplex.

This uses analytic continuation of the interface equations for one
fixed partition. It leaves open isolated ties and branches whose
partitions vary with the correlation.
