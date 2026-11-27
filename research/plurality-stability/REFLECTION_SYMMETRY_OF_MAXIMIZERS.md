# Reflection constraints on balanced Gaussian maximizers

Let \(A=(A_1,\ldots,A_k)\) be an ordinary balanced global maximizer of
Gaussian noise stability in \(\mathbb R^d\), at \(0<\rho<1\).
Here \(k\ge3\), \(d\ge k-1\), and balance means that each cell has
Gaussian measure \(1/k\). Join labels \(i,j\) when some linear
hyperplane reflection exchanges \(A_i,A_j\) and fixes all other
cells up to null sets. The geometry of the resulting graph gives
the following restrictions.

- If the graph is connected, \(A\) is the regular simplex partition.
- If the graph has two components, the labels within each component
  are selected by a regular simplex on its own subspace. Membership
  in one component or the other is determined by a single threshold
  graph in an orthogonal direction, allowing the thresholds to be
  infinite on some fibers.
- On a compact positive correlation band, approximate realization
  of a spanning tree by reflections forces an actual maximizer close
  to the simplex orbit. In the four-cell active space the all-correlation
  local gap then excludes nonsimplex maximizers with sufficiently small
  spanning-tree reflection error.

The single-reflection sign argument is the one used in
[the three-cell reflection reduction](THREE_CELL_REFLECTION_REDUCTION.md#orientation-and-gaussian-rearrangement).
Here the transpositions propagate through the entire graph and fix the
geometry of each component. The result applies to actual global
maximizers, including a potential tetrahedral tie at the first
continuation endpoint.

## Sorting a reflected pair

Write \(\gamma=\gamma_d\), let \(T_\rho\) denote Gaussian noise, and let
\(M_\rho\) denote its kernel relative to \(\gamma\). Set
\(S_\rho(A)=\sum_l\langle 1_{A_l},T_\rho 1_{A_l}\rangle\).
Consider a transposition realized by a reflection \(R\), with
unit normal \(n\) and positive halfspace \(H^+\).
The sum \(1_{A_i}+1_{A_j}\) is even under \(R\), and the difference
\(1_{A_i}-1_{A_j}\) is odd. On \(H^+\), write \(s=1_{A_i}-1_{A_j}\).
Replace \(s\) there by \(|s|\) and extend it oddly, keeping the union of
the two cells and all other labels fixed. The new partition \(A^+\)
has the same cell masses, since the pair still consists of reflected
sets with the same total mass. By expanding the pair contribution
into the sum and difference functions,
\[
 S_\rho(A^+)-S_\rho(A)
 =\iint_{(H^+)^2}
       [M_\rho(x,y)-M_\rho(x,Ry)]
       \big(|s(x)s(y)|-s(x)s(y)\big)\,d\gamma(x)d\gamma(y).
 \tag{RG1}
\]
The kernel difference is positive on the open halfspace product:
the Gaussian radial factor is unchanged by reflection, whereas
\(x\cdot y-x\cdot Ry=2(n\cdot x)(n\cdot y)>0\).
Thus if both signs of \(s\) occur on sets of positive measure, the
improvement is strict. At a maximizer, orient \(n\) so that \(A_i\)
lies in \(H^+\) and \(A_j\) in the opposite halfspace. If
\(b_l=\int_{A_l}x\,d\gamma(x)\), then \(b_j=R b_i\) and
\[
 b_i-b_j=2\left(\int_{A_i}n\cdot x\,d\gamma(x)\right)n .
 \tag{RG2}
\]
The scalar is strictly positive because \(A_i\) has positive measure
and lies on the positive side. In particular the oriented reflecting
hyperplane is determined by the two centroids.

## The components determine the internal simplex cells

The edge transpositions in a connected component \(C\) generate the
whole symmetric group on \(C\). Compose their geometric realizations
to realize any permutation of that component. For any pair in \(C\),
conjugate a generating reflection by a composition sending its two
labels to the prescribed pair. The conjugate is again a linear
hyperplane reflection, and exchanges just the prescribed labels.
It follows from (RG2) that the centroid distances between distinct
labels in \(C\) are positive. They are all equal because the whole
symmetric group acts on these centroids by orthogonal maps.

Put \(b_C=|C|^{-1}\sum_{i\in C}b_i\) and
\(V_C=\operatorname{span}\{b_i-b_j:i,j\in C\}\). The vectors
\(b_i-b_C\) form a centered regular simplex in \(V_C\),
which has dimension \(|C|-1\). This includes a singleton, whose
subspace is zero. Each reflection belonging to \(C\) fixes every
centroid outside \(C\). Its normal is a centroid difference in \(C\),
and these differences span \(V_C\). Thus \(V_C\) is orthogonal to
each \(V_D\) for \(D\ne C\). The mean \(b_C\) is fixed by every
reflection in every component, so it belongs to the orthogonal
complement of the sum of the \(V_D\).

Let \(V=\bigoplus_C V_C\) and \(W=V^\perp\). The generated group acts
as the product of the simplex reflection groups on the \(V_C\) and
as the identity on \(W\). For each \(p\in V\) and component \(C\), let
\[
 \ell_C(p)=\arg\max_{i\in C}\langle p,b_i-b_C\rangle .
\]
The maximizing label is unique off a finite union of hyperplanes
(with the singleton convention). By the pair sorting just proved,
for almost every \(x=p+z\) with \(z\in W\), its label in \(A\) is among
the \(\ell_C(p)\), one for each component. The choice of component is
invariant under the product reflection group. If there is only one
component, this proves the regular simplex conclusion.

## Two components give parallel halfspace fibers

Suppose there are exactly two components \(C,D\). Set
\(u(p,z)=1_{\cup_{i\in C}A_i}(p+z)\), write \(u_p=u(p,\cdot)\),
and let \(m(p)=\int_W u_p\,d\gamma_W\).
For \(p,q\in V\), let \(\chi_C(p,q)\) be one if the winning labels
\(\ell_C(p),\ell_C(q)\) agree, and zero otherwise, and define
\(\chi_D\) likewise. Independence of the Gaussian noises in the
orthogonal factors \(V,W\) gives
\[
\begin{split}
 S_\rho(A)=\iint M_\rho^V(p,q)\Big[
   (\chi_C+\chi_D)\langle u_p,T_\rho^W u_q\rangle
   +\chi_D(1-m(p)-m(q))
 \Big]\,d\gamma_V(p)d\gamma_V(q).
 \tag{RG3}
\end{split}
\]
Indeed a success in the first component requires that both points
belong to it and have the same winning label. Similarly the other
component contributes \(\chi_D\langle1-u_p,T_\rho^W(1-u_q)\rangle\).
Expansion gives (RG3).

Replace each fiber \(u_p\) by an upper halfspace of the same Gaussian
measure in \(W\), with a common normal. Assign the resulting points
using the same winning labels. This preserves the mass of every
label since \(m(p)\) is unchanged at each \(p\). The classical
Gaussian two-set inequality with its equality case
([Mossel–Neeman, Theorems 1.1 and 1.3, pp. 435–436](https://ems.press/content/serial-article-files/32017))
says that every two-fiber interaction in (RG3) is maximized by
parallel halfspaces of these masses. All coefficients of the
interactions are nonnegative. Global maximality therefore forces
equality for almost every pair with \(\chi_C+\chi_D>0\).

Fix one pair of winning labels \((i_0,j_0)\), one from each
component, and consider the region
\(Q=\{\ell_C=i_0,\ell_D=j_0\}\) in \(V\). Its translates under the
product reflection group cover \(V\) up to null sets. The fiber sets
are unchanged by these reflections, which act trivially on \(W\).
If there is a positive-measure set of nontrivial fibers, some and
hence every such \(Q\) meets that set in positive measure. On
\(Q\times Q\) the interaction coefficient is \(2\).
Fix a nontrivial reference fiber in \(Q\) for which equality holds
with almost every other nontrivial fiber in \(Q\), by Fubini.
The equality case for two nontrivial ordinary sets forces all of
these fibers to be halfspaces with a common direction. Reflection
invariance in \(V\) propagates that direction to all the other
regions. Trivial fibers fit the same description with thresholds
\(\pm\infty\). If all fibers are trivial almost everywhere, one
can choose any direction in \(W\). Thus for a unit \(w\in W\) and
a measurable extended-real \(H\) on \(V\),
\[
 \bigcup_{i\in C}A_i=\{p+z:\langle w,z\rangle>H(p)\}
 \quad\hbox{up to null sets}.
 \tag{RG4}
\]
Here \(\dim W=d-k+2\ge1\). For four cells with a disconnected graph
and exactly two components, the possibilities are \(3+1\) and
\(2+2\). In either case the internal directions span a plane.
In the three-dimensional active space, (RG4) is a graph along the
remaining orthogonal coordinate.

## Approximate reflection symmetries

Fix \(0<a\le b<1\), \(k\ge3\) and \(d\ge k-1\). Define the
labelled partition distance by
\(d_\gamma(A,B)=\sum_i\gamma(A_i\triangle B_i)\).
A reflection error for an edge \(ij\) and reflection \(R\)
is \(d_\gamma(A,R(\pi_{ij} A))\), where \(\pi_{ij}\) transposes
the labels. Suppose that at an actual global maximizer \(A\) for
some \(\rho\in[a,b]\), every edge of some spanning tree has a
reflection of error at most \(\varepsilon\). Then there are fixed
constants \(C,\varepsilon_0>0\), depending on \(a,b,k,d\), such that
for \(0<\varepsilon<\varepsilon_0\),
\[
 \inf_{T'\ {\rm simplex}}d_\gamma(A,T')
     \ \le C \varepsilon^{(1-b)/16}.
 \tag{RG5}
\]
The infimum is over orthogonal images and relabelings of the
regular simplex partition. The zero-error case was proved above.

First consider a single approximate transposition, with positive
error \(\delta\). For a transposition of zero error, the exact
sorting result already applies. Copy the labels on \(H^+\) to \(H^-\) by reflected
transposition. This changes the partition by \(d_\gamma\) at most
\(\delta/2\), since the two halfspaces contribute equally to the
involution defect. The result is exactly equivariant, though its
masses may need correction. On \(H^+\) treat the swapped pair as one
group and every other label as an individual group. Correct the
masses of these groups to \(1/k\) for the pair and \(1/(2k)\) for
each single group by moving Gaussian mass from excess groups to
deficient groups. The halfspace measure is nonatomic, and any
point moved into the pair group may take either paired label there.
Then repeat the equivariant copy on \(H^-\). The correction cost
in \(d_\gamma\) is at most the sum of absolute mass errors of
the uncorrected full partition and hence at most \(\delta/2\).
Denote the resulting balanced reflected partition by \(B\).
It satisfies \(d_\gamma(A,B)\le\delta\).

For any two partitions the change of the event that two labels
agree is supported where at least one of the two labels was
changed. This gives
\(|S_\rho(A)-S_\rho(B)|\le d_\gamma(A,B)\).
Polarization of \(B\) is balanced, so (RG1) and maximality of \(A\)
bound the polarized gain by \(\delta\).
Let \(S_+,S_-\) be the two sign sets of \(s\) on \(H^+\) for \(B\).
Their total mass is \(1/k\). On the region
\(|x|,|y|\le R_0\), \(n\cdot x,n\cdot y\ge t_0\) in \(H^+\), the
Mehler formula gives
\[
 M_\rho(x,y)-M_\rho(x,Ry)
       \ \ge\ 2a t_0^2\,e^{-b R_0^2/(1-b)}.
 \tag{RG6}
\]
To see this, factor the common positive kernel into the radial
factor, the tangential exponential and
\(2\sinh(\rho(n\cdot x)(n\cdot y)/(1-\rho^2))\).
The first two factors are at least
\(\exp[-(\rho^2+\rho)R_0^2/(1-\rho^2)]\), and
\(\sinh t\ge t\) for \(t\ge0\).

Take \(t_0=\delta^{1/8}\) and
\(R_0^2=(1-b)\log(1/\delta)/4\). The strip discarded from
\(H^+\) has mass at most \(\phi_0 t_0\), with
\(\phi_0=(2\pi)^{-1/2}\). The radial tail has mass at most
\(2^{d/2}e^{-R_0^2/4}\) by the Gaussian exponential moment.
Thus the discarded mass is bounded by a constant times
\(\delta^\eta\), where \(\eta=(1-b)/16\).
At least one sign set has mass at least \(1/(2k)\); its retained
mass is at least \(1/(4k)\) for small \(\delta\).
On the retained sign product the kernel in (RG6) is at least
\(2a\delta^{1/2}\). The gain bound in (RG1) implies the retained
minority mass is \(O(\delta^{1/2})\). Hence the full minority
mass is \(O(\delta^\eta)\). By reflection and the correction
bound for \(B\), after orienting \(n\) we obtain
\[
 \gamma(A_i\cap H^-)+\gamma(A_j\cap H^+)=O(\delta^\eta).
 \tag{RG7}
\]

The centroid difference for \(B\) is parallel to \(n\) by equivariance.
In the orientation of (RG7) its length is bounded below away from zero
for small \(\delta\). Indeed choose a fixed positive height whose
Gaussian strip has mass at most \(1/(4k)\). The positive-side portion
of \(B_i\) above that height has mass at least \(1/(2k)\) for small
\(\delta\), whereas the integral of \(|n\cdot x|\) over its wrong-side
portion tends to zero by Cauchy–Schwarz. The centroid errors
\(|b_l(A)-b_l(B)|\) are \(O(\sqrt\delta)\), also by Cauchy–Schwarz.
Thus the nonzero unit direction of \(b_i(A)-b_j(A)\) differs from
\(n\) by \(O(\sqrt\delta)\). The Gaussian measure of the symmetric
difference of two central halfspaces is their normal angle divided
by \(\pi\), by rotation in the normal plane. Consequently (RG7)
also holds using the oriented centroid-difference halfspaces.

For the spanning tree, every label permutation can be expressed
as a word of bounded length in its edge transpositions, uniformly
over the finitely many trees on \(k\) labels. Composing the assigned
reflections gives an orthogonal map realizing that permutation
with error \(O_k(\varepsilon)\) by the triangle inequality.
Conjugating a generating reflection then gives approximate
transposition reflections for all pairs, still with error
\(O_k(\varepsilon)\). The preceding argument bounds all
centroid-difference lengths below and gives the wrong-side
probability \(O(\varepsilon^\eta)\) for every pair. The same
compositions and Cauchy–Schwarz show that all centroid distances
differ from any chosen one of them, \(\Delta\), by
\(O(\sqrt\varepsilon)\).

Let \(B_A\) be the matrix of centroid rows. Since their sum is zero,
double centering their squared distances yields
\(B_AB_A^\top=(\Delta^2/2)P+O(\sqrt\varepsilon)\), where
\(P=I-k^{-1}\mathbf1\mathbf1^\top\). On \(\mathbf1^\perp\) the right
side is bounded below for small \(\varepsilon\). Writing the
polar decomposition there as \(B_A=(B_AB_A^\top)^{1/2}Q\)
gives \(QQ^\top=P\). Then the rows of \((\Delta/\sqrt2)Q\)
are regular simplex vectors and differ from the corresponding
centroid rows by \(O(\sqrt\varepsilon)\), by the spectral bound
around the scalar matrix on \(\mathbf1^\perp\).
Their pair normals change by \(O(\sqrt\varepsilon)\).
Apply the central halfspace angle formula again and take the
union of the wrong-side events over pairs. This gives (RG5).

For four cells in dimension three, the
[all-correlation interface gap](TETRAHEDRAL_STANDARD_GAP.md) and
[(B10)–(B11)](POSITIVE_CORRELATION_CONTINUATION.md#an-interface-gap-implies-a-uniform-measurable-local-comparison)
give a neighborhood of the simplex orbit with strict local
comparison, uniformly on \([a,b]\) by a finite cover of the band.
With \(\varepsilon\) sufficiently small, (RG5) places the
maximizer in this neighborhood, so it is tetrahedral. In
particular a nonsimplex tied maximizer at a hypothetical
finite tetrahedral transition has a positive spanning-tree
reflection defect on any compact correlation band containing it.
This excludes approximate full reflection symmetry in the active
space. It leaves open the existence of a sufficiently asymmetric
tied maximizer and the classification of threshold functions
in (RG4).
