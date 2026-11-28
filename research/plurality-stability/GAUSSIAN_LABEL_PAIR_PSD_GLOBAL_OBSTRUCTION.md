# Global obstruction for unanchored label-pair positivity

Balanced scalar agreement PSD and unanchored label-indexed
pair PSD have the same possible same-label kernels. In the
Gaussian plane their common relaxation without triangle
constraints has a strict objective gap above the regular
three-cell simplex at every \(0<\rho<1\). The objective gap
also gives improving kernels arbitrarily close to the
genuine labeled simplex kernel. At correlations sufficiently
close to one, the earlier three-facet witness gives the
same conclusions even with every separate four-point label
test.
[The four-point objective counterexample](GAUSSIAN_AGREEMENT_MULTI_FACET_COUNTEREXAMPLE.md)

## Label averaging preserves the scalar feasible set

On a probability space with measure \(\gamma\), let \(F\)
be a symmetric kernel with \(0\le F\le1\), row mean \(1/k\),
and positive semidefinite integral operator on \(L^2(\gamma)\).
Its lift to a pair-label kernel is
\[
 P^{ii}(x,y)=\frac{F(x,y)}{k},\qquad
 P^{ij}(x,y)=\frac{1-F(x,y)}{k(k-1)}\quad(i\ne j).                \tag{LU1}
\]
Its entries are nonnegative, each one-point label marginal
is \(1/k\), and its same-label sum is \(F\).
Along the constant direction in label space, the integral
operator of \(P\) has constant scalar kernel \(1/k\).
On the label-zero-sum subspace the scalar kernel is
\((F-1/k)/(k-1)\), which is positive semidefinite. Indeed
the constant function is an eigenvector of the scalar
operator of \(F\) of eigenvalue \(1/k\). Subtracting that
constant block leaves a positive operator on the
mean-zero functions. Thus \(P\) is positive semidefinite
on the full label-indexed \(L^2\) space.

Conversely the same-label sum of a positive pair-label
operator has a positive scalar integral operator. Its
\(i\)-th diagonal block is the compression to the \(i\)-th
label, and the sum of those positive operators is positive.
Nonnegativity and the pair marginals give \(0\le F\le1\),
and the same-label row condition gives row mean \(1/k\).
Averaging any such \(P\) over common permutations of all
label indices gives exactly (LU1) for its scalar sum.

Likewise, if a family of tuples admits label assignment
laws with the prescribed scalar agreements \(F\), averaging
each law over a uniform common permutation of its labels
gives pair marginals (LU1). If the input laws agree on
overlaps, the averaging preserves those shared marginals.
This observation concerns pair positivity and the given
tuple laws. It supplies no positivity for a kernel
conditional on labels at sampled anchor points.

## An angular improvement at every positive correlation

Work in the Gaussian plane. Let \(\mathbb T\) have uniform
probability measure \(dm\), and write \(A=[-\pi/3,\pi/3]\)
as a subset of \(\mathbb T\). The Gaussian radii have law
\(d\mu(r)=r e^{-r^2/2}\,dr\). If \(M_\rho(r,s,t)\) denotes
the relative Gaussian Mehler density in polar coordinates,
put
\[
        k_\rho(t)=\iint M_\rho(r,s,t)\,d\mu(r)d\mu(s).
\]
This even angular kernel has integral one and strictly
decreases on \((0,\pi)\) for every \(0<\rho<1\).
For fixed \(r,s>0\) its angular factor is
\(\exp\{\rho rs\cos t/(1-\rho^2)\}\). Gaussian domination
also permits angular differentiation after integrating
the two radii. Denote the circular convolution operator
with kernel \(k_\rho\) by \(K_\rho\).
For an angular-only pair kernel \(F(t)\), the Gaussian
objective is \(\int k_\rho(t)F(t)\,dm(t)\).

Choose the fixed, mean-zero angular perturbation
\[
 g(\theta)=1_A(\theta)\left(\cos\theta
                      -\frac{\int_A\cos\varphi\,dm(\varphi)}{m(A)}\right),
 \quad B=\|g\|_2^2,\quad h_\epsilon=1_A+\epsilon g .
\]
For sufficiently small \(\epsilon>0\), \(h_\epsilon\ge0\),
its mean is \(1/3\), and its squared norm is
\(1/3+\epsilon^2 B\). Define
\[
 \mathcal C_h(t)=\int h(\theta)h(\theta+t)\,dm(\theta),\qquad
 F_\epsilon(t)=
 \frac{\mathcal C_{h_\epsilon}(t)+\epsilon^2 B/2}
      {1/3+3\epsilon^2 B/2}.                                    \tag{LU2}
\]
The autocorrelation is a positive kernel because it is
the inner product of translates in \(L^2(dm)\).
Its values lie between zero and \(\|h_\epsilon\|_2^2\).
Thus \(F_\epsilon\) is positive, has range in \([0,1]\),
and satisfies \(F_\epsilon(0)=1\). Since
\(\int\mathcal C_h\,dm=(\int h\,dm)^2\), its angular
mean is exactly \(1/3\). Lifting it by the angular
difference gives a feasible balanced Gaussian scalar
PSD kernel, and (LU1) gives the label-pair lift.
At \(\epsilon=0\), its scalar kernel is the rotation
average of the genuine simplex agreement, because the
overlap of two arcs \(A\) at angular distance \(t\in[0,\pi]\)
has length \((2\pi/3-t)_+\).

Its objective derivative at zero is
\[
 \left.\frac{d}{d\epsilon}\int k_\rho F_\epsilon\,dm
                          \right|_{\epsilon=0}
                    =6\langle g,K_\rho 1_A\rangle>0.            \tag{LU3}
\]
To check the sign, \(K_\rho1_A\) is even and strictly
decreases for \(0<\theta<\pi/3\). Its derivative there is
\([k_\rho(\theta+\pi/3)-k_\rho(\theta-\pi/3)]/(2\pi)<0\).
The covariance of this function and \(\cos\theta\) on
the uniform arc \(A\) is strictly positive. This follows
by integrating the product of their differences on two
independent points of \(A\). The derivative in (LU3)
is therefore strictly positive for every such \(\rho\).
For each fixed \(\rho\) a sufficiently small positive
\(\epsilon\) gives an objective strictly above the
simplex. This angular construction lies outside the
triangle relaxation, whose exact circle comparison
applies at this angular mean.
[The circle triangle comparison](GAUSSIAN_AGREEMENT_TRIANGLE_RELAXATION.md#the-angular-inequality-and-an-exact-positive-certificate)

## Local consequences of convexity

The bare pair-label feasible class is convex. Mix the
lift of any improving scalar kernel with the genuine
labeled simplex pair kernel
\(P_T^{ij}(x,y)=1_{T_i}(x)1_{T_j}(y)\) with weight
\(\delta\in(0,1)\) on the improving kernel. The gain is
multiplied by \(\delta\), and the entries of \(P\) and its
one-point marginals converge to those of \(P_T\) uniformly
as \(\delta\downarrow0\). Thus no neighborhood of \(P_T\)
defined by pair-kernel or marginal distance suffices
for local optimality in the bare relaxation. The same
argument applies with every pointwise four-point label
test and the triangle inequalities at correlations where
the three-facet witness improves, by first lifting that
witness with (LU1). The mixture need only satisfy separate
four-point tests and does not yield shared triple marginals.

A support restriction on the permitted marginal transfers
can give a different local result. In particular the
signed comparison on compact strips assumes the
one-point marginals are fixed outside those strips.
[The compact-facet comparison](GAUSSIAN_LABEL_PAIR_COMPACT_FACETS_LOCAL_COMPARISON.md)

## Hard marginals outside a large disk

The same lack of metric-local optimality occurs with a hard
exterior. Suppose a scalar witness has representation
\(F(x,y)=1/3+\langle\phi(x),\phi(y)\rangle\), and for a disk
\(B_R\) centered at zero one has
\(\int_{B_R}\phi\,d\gamma_2=0\). The angular-only kernels
(LU2) have this property for every radius, since the
mean of the angular Gram feature on every circle is zero.
For the unaveraged three-strip witness, it holds for every
disk containing the strips. There \(\phi(x)=e_{T(x)}-m\)
outside the disk, where \(m=(e_1+e_2+e_3)/3\), so its
exterior integral and its full-space integral both vanish.

Let \(d_i=e_i-(e_1+e_2+e_3)/3\) in a separate copy of
label space and \(v_*\) be a unit vector perpendicular
to the tensor-product feature space. On \(B_R\) use
\(V_i(x)=v_*/3+d_i\otimes\phi(x)/\sqrt2\), which realizes
(LU1) for \(k=3\). Outside use \(V_i(x)=1_{T_i}(x)v_*\).
The mixed pairs have nonnegative joint law
\(P^{ij}(x,y)=1_{T_j}(y)/3\) for \(x\in B_R\), \(y\notin B_R\).
The vector means remain \(\int V_i\,d\gamma_2=v_*/3\),
by the feature integral and equal exterior masses.
Hence this replacement preserves exact row balance
and pair positivity. To check local label laws, retain
the label-permutation averaged law on the points inside
the disk and append the deterministic labels of the
exterior points. In particular this preserves separate
four-point and triangle tests for the three-strip witness.

The deficit introduced by replacement is at most
\(2\gamma_2(B_R^c)\), by the union bound under the
correlated pair law. For each prescribed strict gain
before replacement, take \(R\) large enough to retain
strict gain. Dilution with the true labeled simplex
kernel now gives arbitrarily small improving
perturbations with hard marginals outside the same disk.
This conclusion holds for all positive correlations in
the bare relaxation, and for correlations sufficiently
close to one with separate four-point label tests.
