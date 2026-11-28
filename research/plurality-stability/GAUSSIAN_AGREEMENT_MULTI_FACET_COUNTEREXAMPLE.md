# A local multiray counterexample to the strengthened agreement relaxation

For all positive correlations sufficiently close to one, the
regular planar tripod fails even as a local maximizer of the
Gaussian agreement relaxation with positive semidefiniteness
and level-four label consistency. The competitors can differ
from its genuine agreement kernel only on pairs having at
least one endpoint in arbitrarily thin strips about compact
parts of its three facet rays. Their diagonal entries are one,
their row means are exactly \(1/3\), and they satisfy every
triangle and every four-point pigeonhole inequality. Here
level-four consistency means that each individual four-point
agreement matrix is the expectation of a three-label assignment
on that tuple. The construction satisfies this condition and
fails level five before rotation averaging. Simultaneous rotation
averaging preserves the level-four constraints and the strict
gain at the Gaussian objective.

The ordinary tripod's rotation mode uses three label incidence
vectors with mutual cosine \(-1/2\). The relaxation admits
feature directions orthogonal to the three labels with mutual
cosine \(+1/2\). This reverses the sign of their cross-ray
interaction. On the noise scale near the origin it makes a
three-ray mode of the relaxed form strictly favorable.
[The triangle relaxation and the earlier mixed-radius witness](GAUSSIAN_AGREEMENT_TRIANGLE_RELAXATION.md)
[Finite label consistency](GAUSSIAN_AGREEMENT_TRIANGLE_DISCRETIZATION.md#a-convergent-hierarchy-using-finite-label-consistency)
[Local comparison on one compact facet](GAUSSIAN_AGREEMENT_COMPACT_FACET_LOCAL_COMPARISON.md)

## An exactly balanced feature perturbation

Write \(T_i\) for the three cells, with cyclically oriented
facets \(e=12,23,31\). Let \(s>0\) be arclength on a ray and
\(t\) its signed perpendicular coordinate. The Gaussian
arclength measure is \(d\nu(s)=(2\pi)^{-1}e^{-s^2/2}\,ds\).
Fix a bounded real \(b\) supported in a compact subinterval
of \((0,\infty)\) with \(\int b\,d\nu=0\). At each \(s\)
use a symmetric strip \(|t|<h_\epsilon(s)\), with width
chosen by
\[
      \int_{-h_\epsilon(s)}^{h_\epsilon(s)}
          e^{-t^2/2}\,dt=2\epsilon |b(s)| .
\]
For each fixed \(b\), the strips are disjoint for all small
\(\epsilon>0\), and \(h_\epsilon(s)=\epsilon|b(s)|+O(\epsilon^3)\)
uniformly. Choose orthonormal label vectors \(e_1,e_2,e_3\)
and, in an orthogonal space, unit vectors \(u_e\) with
\(u_e\cdot u_f=1/2\) for distinct facets. This Gram matrix is
positive definite. Define the feature map
\[
 V(x)=
 \begin{cases}
     e_{T(x)},&x\text{ outside the strips},\\[2mm]
     (e_i+e_j)/2+\operatorname{sgn}(b(s))u_{ij}/\sqrt2,
                         &x\text{ in the strip of }ij .
 \end{cases}                                                \tag{MW1}
\]
All features have unit norm. Symmetry of the normal width
preserves the three soft label masses, each \(1/3\). The new
orthogonal mean on an individual strip is
\(\sqrt2\,\epsilon u_e\int b\,d\nu=0\). Thus
\(\int V\,d\gamma_2=(e_1+e_2+e_3)/3\) and
\(F_\epsilon(x,y)=\langle V(x),V(y)\rangle\) has row mean
exactly \(1/3\). It is symmetric and positive semidefinite.

There are three anchor states \(e_i\) and two ribbon states
for each pair. Equal states agree by one. Distinct states
agree by \(0\) or \(1/2\). An anchor and ribbon state agree
by \(1/2\) exactly when the anchor belongs to the ribbon's
pair. Two ribbon states of different pairs agree by \(1/2\)
exactly when their signs match. Opposite signs on the same
pair agree by zero. If two states coincide their two other
agreements coincide,
and otherwise the sum of any two agreements is at most one.
Thus all triangle inequalities hold.

## The finite label tests

On distinct states connect two vertices when their agreement
is \(1/2\). This graph has no induced claw, meaning a center
with three mutually nonadjacent neighbors. The neighbors of
an anchor are the four ribbon states of the two incident
pairs. Among any three neighbors there is a same-sign pair,
which is adjacent. The neighbors of a ribbon state are its
two anchors and the two same-sign states of the other pairs.
The latter two are adjacent, and each is adjacent to one
of the two anchors. This also precludes a claw.

Every four distinct states have at least two edges. Four
ribbons give this by counting same-sign pairs. Three ribbons
and an anchor either give two anchor edges or consist of
two opposite-pair ribbons of opposite signs and a third
ribbon incident to the anchor. The third ribbon then gives
an anchor edge and an edge to its same-sign ribbon. Two
ribbons and two anchors give an anchor edge for each ribbon.
Three anchors and one ribbon give its two incident edges.
This proves the four-point pigeonhole inequality, including
tuples with repetitions, which already have a unit agreement.

For full four-point label consistency, each edge must occur
as an agreement with probability \(1/2\) in a partition of
the four vertices into at most three classes. Nonedges must
never agree. The possible four-vertex graphs with at least
two edges and no claw admit the following laws:

- With two edges, choose one of the single edges with equal
  probability.
- For a three-edge path, mix its two outer edges as a
  matching with its middle edge. For a triangle and an
  isolated vertex, use the triple and each of its three
  single edges with probability \(1/4\).
- On a four-cycle mix its two perfect matchings. For a
  triangle with a pendant edge, mix the triple and the
  single pendant edge.
- If the graph is \(K_4\) minus edge \(14\), use the triples
  \(123,234\) and the matchings \(\{12,34\},\{13,24\}\),
  each with probability \(1/4\). On \(K_4\), mix the four
  triples equally.

Vertices not in a listed block are singleton classes.
Repeated states inherit the same label, and tuples of fewer
than four distinct states inherit a four-state law by
restriction. This proves the required local condition.

Level five fails for the two anchors \(1,2\) and the three
positive ribbon states \(12,13,23\). The anchors are distinct
almost surely in any proposed assignment. Let \(X\) be the event
that ribbon \(12\) chooses anchor \(1\), \(Y\) that \(13\)
chooses anchor \(1\), and \(Z\) that \(23\) chooses anchor \(2\).
All have probability \(1/2\). The other possibilities are
anchor \(2\) for ribbon \(12\) and the third label for the
latter two ribbons. The ribbon agreements \(1/2\) for the
first two pairs force
\(\Pr(X\cap Y)=\Pr(X^c\cap Z)=1/2\). Hence \(Y=X\) and
\(Z=X^c\) almost surely. Ribbons \(13,23\) then never agree,
contrary to their stipulated agreement \(1/2\). The witnesses
below give each of these states positive measure. The
four-point property is preserved by rotation averaging;
the level-five failure here concerns the unaveraged kernel.

The pointwise four-point laws do not share a consistent
three-point marginal. Take the triple consisting of anchor
\(1\) and the positive ribbons \(12,13\). Its three pair
agreements are \(1/2\). Adding negative ribbon \(23\) gives
an isolated vertex, so the triple uses at most two labels.
The expected number of agreements in the triple is \(3/2\),
forcing all three to agree with probability \(1/4\). Instead
add negative ribbon \(12\), which has only the pendant
agreement \(1/2\) with the anchor. The two original anchor
agreements are both disjoint from the pendant event and
each has probability \(1/2\), forcing the entire triple
to agree with probability \(1/2\). Thus the stronger
constraint of compatible overlapping laws rejects the
unaveraged kernel already at four points.

## The raw three-ray variation

Put \(\sigma^2=1-\rho^2\). For two points of radial
coordinates \(s,t>0\), let \(M_+(s,t)\) and \(M_-(s,t)\)
be the Mehler kernel with angular cosine \(1\) and \(-1/2\),
respectively. In this section the pairing and operator
notation for these kernels use \(d\nu\). The local coefficient
on a ray is the raw score derivative magnitude \(a(s)\). It
satisfies
\[
                   a(s)=\rho\big(2M_+1(s)+M_-1(s)\big).
                                                               \tag{MW2}
\]
Indeed the Gaussian boundary derivative of a pair score is
\(-\rho\sum_f(d_e\cdot d_f)\int M_\rho N_f\,d\nu\), in
cyclic orientations with \(d_e=e_i-e_j\) and normals \(N_f\).
Projecting onto \(-N_e\) gives the same-ray factor two.
On either other ray \(d_e\cdot d_f=-1\) and
\(N_e\cdot N_f=-1/2\), giving a factor \(1/2\) for each.
[Translation and score derivative identities](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)

For \(w(x)=V(x)-e_{T(x)}\), the exact objective deficit from
the genuine tripod is its soft label score regret doubled,
less \(\iint M_\rho\langle w(x),w(y)\rangle\,d\gamma_2^2\).
On a strip the soft transfer is \(1/2\) to the other label.
The score gap is \(a(s)|t|+O(t^2)\), so the doubled regret
over all strips is \(3\epsilon^2\int a b^2\,d\nu+o(\epsilon^2)\).
The transverse Gaussian integral of \(w\) is
\(\sqrt2\,\epsilon b(s)u_e\): its soft label change cancels
by normal symmetry. Replacing the smooth kernel by its
restriction to the rays has an error \(o(\epsilon^2)\), since
the strips lie in a fixed compact set, have total mass
\(O(\epsilon)\), and their maximum width tends to zero.
The same-ray interaction contributes \(6\epsilon^2
\langle b,M_+b\rangle\), and the six ordered cross-ray
interactions contribute \(6\epsilon^2\langle b,M_-b\rangle\).
Consequently the deficit is
\[
 3\epsilon^2 Q_\rho(b)+o(\epsilon^2),\qquad
 Q_\rho(b)=\int a b^2\,d\nu-2\langle b,(M_++M_-)b\rangle .
                                                               \tag{MW3}
\]

## A negative direction on the noise scale

Take first \(b_\sigma(s)=B(s/\sigma)\) for a fixed smooth
compactly supported \(B\) on \([0,\infty)\).
Under \(s=\sigma u,\ t=\sigma v\), the joint density on the
two rays becomes
\[
 M_\pm(\sigma u,\sigma v)\,d\nu(\sigma u)d\nu(\sigma v)
  =(2\pi)^{-2}\exp[-(u^2+v^2)/2+\rho c_\pm uv]\,du\,dv ,
 \qquad c_+=1,\quad c_-=-1/2.
\]
By (MW2) and dominated convergence, \((2\pi)^2Q_\rho(b_\sigma)\)
therefore tends to
\[
 \begin{split}
 &\iint K_+(u,v)(B(u)-B(v))^2\,du\,dv\\
 &\qquad+\iint K_-(u,v)
           \big(B(u)^2-2B(u)B(v)\big)\,du\,dv,              \tag{MW4}\\
 &K_+(u,v)=e^{-(u-v)^2/2},\qquad
          K_-(u,v)=e^{-(u^2+v^2+uv)/2}.
 \end{split}
\]
All integrals are on the nonnegative half-line in each
variable. In the same-ray terms there is a compactly
supported factor in at least one variable before taking the
limit. The product is dominated by \(e^{-(u-v)^2/2}\)
times that factor. The cross-ray kernel is dominated by
\(e^{-(u^2+v^2)/2}\). These bounds justify the limit and
symmetrizing the first term into its difference form.

Choose \(0\le B_L\le1\) equal to one on \([0,L]\), zero on
\([2L,\infty)\), and smoothly tapered with
\(\|B_L'\|_2^2=O(1/L)\). The first term in (MW4) is
\(O(1/L)\). For each shift \(z\) the half-line integral
of \((B_L(u+z)-B_L(u))^2\) where both arguments are
nonnegative is at most \(z^2\|B_L'\|_2^2\), by the
fundamental theorem of calculus and Cauchy–Schwarz.
Integrate in the shift against its Gaussian weight.
The second term of (MW4) tends to
\(-\iint K_-<0\), by dominated convergence.
Fix \(L\) large enough for strict negativity. Multiply this
cutoff by a smooth factor zero on \([0,\delta]\) and one
beyond \(2\delta\). The change tends to zero in \(L^2\);
the difference form with \(K_+\) is bounded on that space
by the uniform kernel row integral, and the cross-ray term
converges by domination. For sufficiently small fixed
\(\delta>0\) the resulting \(B\) is supported away from
zero and still gives a strictly negative limit.

Finally correct the mean of \(b_\sigma\) by subtracting
\((\int b_\sigma\,d\nu)\psi\), with \(\psi\) a fixed
nonnegative compact smooth bump away from zero with
\(\int\psi\,d\nu=1\). The coefficient is \(O(\sigma)\).
By (MW2) or the Gaussian derivative bound,
\(a=O(\sigma^{-1})\), so the correction has squared
\(L^2(a\,d\nu)\) norm \(O(\sigma)\). The corresponding norm
of \(b_\sigma\) stays bounded since its support has length
\(O(\sigma)\). The interaction in (MW3) is uniformly bounded
in the weighted norm as \(\rho\uparrow1\): for the
nonnegative symmetric kernel \(2(M_++M_-)\), its row integral
is at most \(2a/\rho\) by (MW2). The weighted Schur
estimate therefore applies. The correction changes \(Q_\rho\)
by \(o(1)\), so the corrected, exactly balanced \(b\)
still has \(Q_\rho(b)<0\) for all \(\rho\) sufficiently
close to one. Its positive and negative supports have
positive measure and are disjoint for such \(\rho\).

For these fixed correlations, sending \(\epsilon\downarrow0\)
in (MW3) gives arbitrarily small strict improvements in the
level-four, positive-semidefinite agreement relaxation.
This concerns the relaxed objective. The ordinary three-cell
problem has a separate [exact near-one theorem](THREE_CELL_NEAR_ONE_EXACT.md).

## The symmetric ribbon mode at level five

More generally let \(c_1,c_2,c_3\) be the three cosines
between the positive ribbon residual directions. The same
five points used above enforce \(c_1+c_2+c_3\le1/2\).
In fact their pair agreements are \(1/4+c_i/2\), and the
binary variables \(X,Y,Z^c\) have pair distances
\(1/2-c_1,\ 1/2-c_2,\ 1/2+c_3\). Their triangle inequality
gives the bound. For a symmetric cosine \(c\) it requires
\(c\le1/6\). The limiting cross-ray coefficient in the
slow-cutoff test (MW4) for this symmetric family is \(1-4c\);
the negative limit uses \(c>1/4\). Thus this particular
instability mechanism is excluded by the five-point anchor
test.

## Consistent label laws on the whole ribbon model

In the feature model (MW1), allow separate bounded profiles
\(b_e\) with zero Gaussian arclength means and with both signs
present in positive measure on every ribbon. The unit feature
direction \(u_e\) is constant along each ribbon, but its three
mutual cosines \(c_{ef}\) are arbitrary. For this model,
pointwise level-five consistency forces the entire agreement
kernel to be a finite mixture of genuine balanced partitions.

Indeed, using either sign of each ribbon in the level-five
test replaces \(c_{ef}\) by \(\eta_e\eta_f c_{ef}\) for any
\(\eta\in\{-1,1\}^3\). Hence
\[
 1-2\sum_{e<f}\eta_e\eta_f c_{ef}\ge0
              \qquad\text{for every }\eta\in\{-1,1\}^3.         \tag{MW5}
\]
For a random sign string \(R\), give the string \(\eta\)
probability equal to the left-hand side of (MW5) divided
by eight. This is a probability law, with each sign mean zero
and \(\mathbb E R_eR_f=-2c_{ef}\). Use the original hard
labels outside the strips. On cyclically oriented ribbon
\(e=ij\), assign label \(i\) when \(R_e\operatorname{sgn}b_e(s)=1\)
and label \(j\) otherwise. The symmetric widths and
\(\int b_e\,d\nu=0\) preserve the three masses for every
sign string. Two points on the same ribbon agree exactly
when their profile signs match. On distinct ribbons with
profile signs \(s_e,s_f\), the common label requires
opposite cyclic binary signs, so the agreement probability
is \((1-s_es_f\mathbb E R_eR_f)/4\), the value of the feature
kernel. The anchor agreements follow from the zero sign
means. This proves the mixture representation on all
points. At every fixed positive correlation its objective
is at most the simplex value for sufficiently thin compact
strips, by the ordinary measurable local comparison.
[Measurable local comparison](GAUSSIAN_SIMPLEX_DIMENSION_FREE_LOCAL_GROWTH.md)

Compatible four-point laws on overlapping tuples already
imply (MW5) in this model. For any chosen signs on the three
ribbons, write \(A_j=1/4+c_j/2\) and \(B_j=1/4-c_j/2\)
for the three cross agreements with matching and reversed
chosen signs, where \(c_1\) refers to \(12,13\).
Consider anchor \(1\) and the chosen states of \(12,13\).
If their all-equal probability is \(t\), their all-distinct
probability is \(2t-A_1\), since the three pair agreements
sum to \(1+A_1\). Adding the reversed state of \(23\)
gives zero agreement with the anchor. On the all-distinct
event it must therefore agree with one of the other two
states, so \(2t-A_1\le B_2+B_3\).
Instead add the reversed state of \(12\). Its agreement
with the anchor, of probability \(1/2\), is the complement
of the anchor's agreement with the original state of
\(12\), because those two ribbon states never agree.
Its overlap with the event that the anchor agrees with
the chosen state of \(13\) is at most \(B_1\). Thus
\(t\ge1/2-B_1\). Shared triple marginals give the same \(t\)
in both laws. The two bounds imply
\(c_1+c_2+c_3\le1/2\). Apply them to all choices of signs.
