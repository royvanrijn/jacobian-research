# The Hurwitz--LL stacky skeleton and toroidal invariants

## Purpose and status

The `H1-COARSE` theorem identifies the coarse polynomial admissible-cover
graph with the corrected wonderful target graph
\[
 q:\mathcal G_N^{\mathrm{poly}}\longrightarrow B_N^{\mathrm{tgt}}.
 \tag{0.1}
\]
The `H1-STACK` theorem identifies the source as the normalized selected
factor in the Kummer log-etale admissible-cover pullback and computes its
actual inertia.

This note extracts the next geometric package. It separates:

1. consequences that are formal after a log-regularity hypothesis;
2. data already computed by the recursive resonance atlas;
3. additional hypotheses needed for ordinary quotient, canonical, and
   stringy-singularity statements; and
4. a degree-six Maxwell pilot that can be computed before the general
   intersection theory.

This is a **conditional consequence theorem and research specification**.
It is not a new unconditional completion-status entry. In particular, the
log regularity of the normalized wonderful pullback, the equivariant
stack-level form of selected-factor normalization, and the tangent action on
all strict deformation parameters must be checked independently before the
word "theorem" is applied without qualification.

Throughout, work in characteristic zero. Let
\[
 D_B=B\setminus B^\circ,\qquad
 D_{\mathcal G}=(q^{-1}D_B)_{\mathrm{red}}.
 \tag{0.2}
\]

## 1. The local stacky monoid package

At a geometric point \(x\to\mathcal G\), let \(P_x\) be the characteristic
monoid of \(B\), let \(E_x\) be the target-node set, and let the source nodes
over \(\eta\in E_x\) have indices
\[
 e_{\eta1},\ldots,e_{\eta r_\eta}.
\]
Put
\[
 Q_\eta=
 \left\langle
 \delta_\eta,\sigma_{\eta1},\ldots,\sigma_{\eta r_\eta}
 \ \middle|\
 \delta_\eta=e_{\eta j}\sigma_{\eta j}
 \right\rangle,
 \qquad Q_x=\bigoplus_{\eta\in E_x}Q_\eta,
 \tag{1.1}
\]
and
\[
 P_{\mathrm{cov},x}
 =
 \left(
 P_x\oplus_{\mathbf N^{E_x}}Q_x
 \right)^{\mathrm{sat}}.
 \tag{1.2}
\]

Let
\[
 L_\eta=\operatorname{lcm}_j(e_{\eta j}).
\]
The labelled-cover automorphism character is
\[
 \chi_x:
 \operatorname{Aut}_{\mathrm{lab}}(C/P)
 \longrightarrow
 \prod_{\eta,j}\mu_{e_{\eta j}},
 \tag{1.3}
\]
and the actual stabilizer of the selected normalized factor is
\[
 I_x=
 \chi_x^{-1}
 \left(
 \prod_{\eta\in E_x}\Delta_\eta(\mu_{L_\eta})
 \right).
 \tag{1.4}
\]

The recursive resonance atlas computes (1.1)--(1.4) and now retains the
additional datum needed for canonical ages and orbifold invariants: the
**diagonal tangent character**
\[
 \theta_x:I_x\longrightarrow\prod_{\eta\in E_x}\mu_{L_\eta}.
 \tag{1.5}
\]
It is defined by writing the phase vector of \(g\in I_x\) in the selected
diagonal factor:
\[
 \chi_{\eta j}(g)
 =
 \theta_\eta(g)^{L_\eta/e_{\eta j}}.
 \tag{1.6}
\]
If
\[
 q_\eta=\rho_\eta^{L_\eta},\qquad
 s_{\eta j}=\omega_{\eta j}
 \rho_\eta^{L_\eta/e_{\eta j}},
 \tag{1.7}
\]
then
\[
 g(\rho_\eta)=\theta_\eta(g)\rho_\eta.
 \tag{1.8}
\]
Thus (1.5) is not a new group enumeration. It is the retained output of the
simultaneous diagonal-intersection calculation already used to obtain
\(I_x\).

## 2. Conditional log-regularity theorem

> **Conditional stacky-skeleton theorem.**
> Assume:
>
> 1. \((B,D_B)\) is log regular;
> 2. selected-factor normalization is valid equivariantly on strict-etale
>    quotient charts, including for a nonrepresentable Kummer pullback;
> 3. the fixed-target fully marked Hurwitz datum has no strict relative
>    deformation parameters; and
> 4. the selected stack is tame with trivial generic stabilizer.
>
> Then \((\mathcal G,D_{\mathcal G})\) is log regular and Kummer log-etale
> over \((B,D_B)\). Strict henselian formally at \(x\), it has a presentation
> \[
> \widehat{\mathcal G}_x
> \simeq
> \left[
> \operatorname{Spf}
> k[[P_{\mathrm{cov},x}]][[z_1,\ldots,z_d]]
> \big/I_x
> \right],
> \tag{2.1}
> \]
> with monoid (1.2), stabilizer (1.4), and node-normal tangent action
> (1.5)--(1.8).

Here the \(z_i\) are strict-etale target parameters. Hypothesis 3 says that
there are no additional source-cover parameters over a fixed target. To use
the age formula below, one must also verify that \(I_x\) acts trivially on
the \(z_i\). This is expected because the fully marked target has no
automorphisms and the relative Hurwitz class is discrete, but it should be
included in the deformation theorem rather than inferred from dimension
alone.

### Proof under the hypotheses

Kummer log-etale morphisms preserve log regularity. The fs fiber product has
characteristic monoid (1.2), since normalization of the toroidal pushout is
saturation. Equivariant selected-factor normalization retains the normal
factor meeting the generic polynomial section and replaces the ambient
automorphism group by its stabilizer. Formula (1.4) identifies that
stabilizer. The admissible-cover deformation equations give (1.7), hence
the action (1.8). The absence of other relative parameters gives (2.1).
\(\square\)

The phrase **toroidal finite-quotient chart** is therefore justified by
(2.1). The stronger phrase **quotient singularity** in the usual
minimal-model sense is not automatic: a toric singularity is a finite
quotient of a smooth toric chart only when the relevant cone is simplicial.
Likewise:

- simpliciality is needed for local \(\mathbf Q\)-factoriality;
- a rational height-one support function is needed for
  \(\mathbf Q\)-Gorensteinness; and
- canonical or terminal singularities require the corresponding lattice
  discrepancy inequalities.

These are finite fan checks once all cones of \(P_{\mathrm{cov},x}\) are
available, but they do not follow from saturation alone.

## 3. Canonical and log-canonical classes

Under the hypotheses of Section 2, logarithmic differentials pull back:
\[
 q^*\Omega_B^{\log}\simeq\Omega_{\mathcal G}^{\log}.
 \tag{3.1}
\]
Taking determinants gives the formal identity
\[
 \boxed{
 K_{\mathcal G}+D_{\mathcal G}
 =
 q^*(K_B+D_B).
 }
 \tag{3.2}
\]

Let \(D_i\) be a boundary prime of \(B\), let \(\mathcal D_i\) be its
reduced inverse image, and let \(m_i\) be the generic stabilizer order:
\[
 q^*D_i=m_i\mathcal D_i.
 \tag{3.3}
\]
When the divisors involved are \(\mathbf Q\)-Cartier, (3.2) is equivalent to
\[
 \boxed{
 K_{\mathcal G}
 =
 q^*
 \left(
 K_B+\sum_i\left(1-\frac1{m_i}\right)D_i
 \right).
 }
 \tag{3.4}
\]

Codimension-two stabilizers do not add Weil-divisor terms. They do affect
orbifold sectors and ages.

Formula (3.2), rather than an unqualified assertion that the coarse space
has canonical singularities, is the strongest immediate consequence.
The toroidal pair \((B,D_B)\) has the expected logarithmic discrepancy
behavior. Ordinary singularities of \(B\) require the
\(\mathbf Q\)-Gorenstein and cone-height checks described above.

## 4. Skeleton and dual complex

For a Kummer map of fs monoids, the induced map of rational cones is an
isomorphism; only the integral lattices change by finite index. Therefore
\[
 \boxed{
 \Sigma_{\mathcal G,\mathbf R}
 \simeq
 \Sigma_{B,\mathbf R}.
 }
 \tag{4.1}
\]
The underlying real skeleton is consequently the already constructed
saturated wonderful pullback
\[
 \operatorname{Sat}
 \left(
 \Sigma_X
 \times_{\mathbf R^b/\mathbf R\mathbf1}
 \mathcal M_{0,b+2}^{\mathrm{trop}}
 \right).
 \tag{4.2}
\]

The new object is not a different topological complex but the **stacky
Hurwitz--LL cone complex**
\[
 \tau\longmapsto
 \left(
 P_\tau,\,
 P_{\mathrm{cov},\tau},\,
 I_\tau,\,
 \theta_\tau
 \right),
 \tag{4.3}
\]
where \(\tau\) is a monodromy-decorated nested collision stratum.
Contraction of a target edge induces:

1. a face map of characteristic monoids;
2. the corresponding saturated map of cover monoids;
3. specialization of the full-centralizer fiber product; and
4. restriction of the diagonal tangent character.

The affine-composition and simultaneous-character results in the recursive
atlas supply exactly these compatibilities.

The ordinary dual complex forgets the last three entries of (4.3).
The labelled-to-unlabelled quotient acts on all four entries and therefore
produces a stacky generalized cone complex rather than only a quotient of
the underlying nested-set topology.

## 5. Inertia strata and age

Stratify \(B\) by locally closed monodromy-decorated nested strata
\(B_\tau^\circ\). On such a stratum the finite group \(I_\tau\) and its
normal character are locally constant. The inertia calculation can
therefore be organized by pairs
\[
 (\tau,[g]),\qquad [g]\subset I_\tau
 \tag{5.1}
\]
with face-specialization maps identifying pairs that lie in one global
inertia component.

Write
\[
 \theta_{\tau,\eta}(g)
 =
 \exp\left(2\pi i\,a_{\eta}(g)/L_\eta\right),
 \qquad 0\le a_\eta(g)<L_\eta.
 \tag{5.2}
\]
If the strict parameters in (2.1) are fixed, the age is
\[
 \boxed{
 \operatorname{age}_\tau(g)
 =
 \sum_{\eta\in E_\tau}
 \frac{a_\eta(g)}{L_\eta}.
 }
 \tag{5.3}
\]
If a residual strict parameter has nontrivial character, its fractional
weight must be added. This is why the tangent representation, not just the
abstract stabilizer order, must be retained.

The additive orbifold \(E\)-polynomial of a smooth proper selected stack
would then be computable stratumwise:
\[
 E_{\mathrm{orb}}(\mathcal G;u,v)
 =
 \sum_F E(F;u,v)(uv)^{\operatorname{age}(F)},
 \tag{5.4}
\]
where \(F\) runs over connected inertia sectors. Equivalently, one can use
the open strata \(B_\tau^\circ\), the conjugacy classes in \(I_\tau\), and
inclusion--exclusion along the face category (5.1).

The following distinctions are essential.

- The orbifold \(E\)-polynomial belongs to the smooth tame stack.
- The stringy \(E\)-polynomial of \(B\) requires \(B\) to be
  \(\mathbf Q\)-Gorenstein and log terminal.
- An arbitrary root stack over a smooth coarse space has nontrivial
  orbifold sectors even though the coarse stringy polynomial is just the
  ordinary polynomial.
- An equality between orbifold and stringy invariants therefore requires a
  separate crepant or canonical-stack statement.

## 6. Chow groups and boundary intersections

For ordinary rational Chow groups, a tame stack with finite stabilizers is
close to its coarse space. Root constructions become more visible
integrally and in orbifold Chow theory. At a divisorial root of order
\(m_i\),
\[
 q^*D_i=m_i\mathcal D_i,\qquad
 q_*[\mathcal D_i]=\frac1{m_i}[D_i].
 \tag{6.1}
\]
For independent normal-crossing roots,
\[
 q_*
 \left(
 \mathcal D_{i_1}\cdots\mathcal D_{i_r}
 \right)
 =
 \frac{
 D_{i_1}\cdots D_{i_r}
 }{
 m_{i_1}\cdots m_{i_r}
 }.
 \tag{6.2}
\]
In a coupled Kummer chart, the product denominator is replaced by the
appropriate lattice index or generic stabilizer of the stratum. Formula
(4.3) contains that information.

Toroidal Cartier divisors are integral piecewise-linear functions on the
cone complex. Their intersections can therefore be computed by tropical
divisor products and balanced Minkowski weights, provided the classes of
the closures of strata are known. This gives a concrete program:

1. express each marked-root, caustic, Maxwell, and wonderful exceptional
   divisor as a piecewise-linear function on \(\Sigma_B\);
2. pull it to the lattice-refined stacky complex;
3. intersect using the fan weights and stratum classes; and
4. divide by the stabilizer/lattice indices prescribed by (6.1)--(6.2).

The existing caustic and Maxwell formulas live on the root-only quotient
\[
 [\overline M_{0,n+2}/S_n],
\]
not on the corrected wonderful graph \(B\). They cannot simply be declared
to be the corresponding classes on \(\mathcal G\). One must compute:

- their strict transforms;
- multiplicities along every wonderful exceptional divisor; and
- the divisorial stack indices \(m_i\).

The recursive critical-value initial forms give the required exceptional
valuations. This pullback-and-correction calculation is the first global
intersection problem to perform.

The orbifold Chow **product** needs more than (6.1): it requires the double
inertia stack, evaluation maps, and the obstruction bundle. The monoids and
tangent characters determine these locally, but a global presentation
should be written only after the inertia face category has been constructed.

## 7. Degree-six Maxwell pilot

On the labelled degree-six central boundary surface, the selected stack is
already identified as
\[
 \mathcal B_6^{\mathrm{Max}}
 =
 \sqrt[2]{(B_6,D_{xy})}
 \times_{B_6}
 \sqrt[2]{(B_6,D_{xz})}
 \times_{B_6}
 \sqrt[2]{(B_6,D_{yz})}
 \times_{B_6}
 \sqrt[2]{(B_6,D_{xyz})}.
 \tag{7.1}
\]
The three pairwise Maxwell curves are mutually disjoint and each meets the
triple-Maxwell exceptional curve once. Hence:

- a generic Maxwell curve has inertia \(\mu_2\);
- a pair--triple crossing has inertia \((\mu_2)^2\);
- each single-root twisted sector has age \(1/2\); and
- the sector nontrivial in both crossing factors has age \(1\).

The canonical class is
\[
 \boxed{
 K_{\mathcal B_6^{\mathrm{Max}}}
 =
 q^*
 \left(
 K_{B_6}
 +\frac12
 (D_{xy}+D_{xz}+D_{yz}+D_{xyz})
 \right).
 }
 \tag{7.2}
\]

Assume for this compact surface calculation that \(B_6\) is the four-point
blowup of \(\mathbf P^2\) described by the target chart and that all four
Maxwell curves are \(\mathbf P^1\)'s with the intersection graph above.
Put \(q=uv\). Then
\[
 E(B_6)=1+5q+q^2.
 \tag{7.3}
\]
There are four age-\(1/2\) divisor sectors and three additional age-one
double-root point sectors. Therefore
\[
 \begin{aligned}
 E_{\mathrm{orb}}(\mathcal B_6^{\mathrm{Max}};u,v)
 &=
 E(B_6)
 +4E(\mathbf P^1)q^{1/2}
 +3q\\
 &=
 \boxed{
 1+4q^{1/2}+8q+4q^{3/2}+q^2.
 }
 \tag{7.4}
\end{aligned}
\]
Formula (7.4) is a pilot calculation for the explicitly identified central
boundary surface, not an orbifold polynomial for the entire degree-six
compactification. It should be checked directly from the inertia stack of
(7.1).

## 8. Four proof obligations before promotion

### 8.1 Normalization of stacks

The selected-factor argument should use stack normalization itself, not only
normalization of an algebraic-space chart. Choose a smooth presentation
\[
 R\rightrightarrows U\longrightarrow\mathcal Y.
\]
Normalize \(U\) in the selected generic field or selected orbit of fields.
Smooth base change identifies the corresponding normalization of \(R\);
the groupoid operations extend uniquely. The quotient of the normalized
selected groupoid is the required normalized substack.

This also states precisely why normalization and the chosen generic factor
descend.

### 8.2 Nonrepresentable Kummer pullback

On a quotient chart, the generic section may select an orbit of minimal
primes rather than one prime invariant under the whole chart group. The
normalized selected substack is the quotient of one selected factor by its
stabilizer. The proof must distinguish:

- the number of normalization branches;
- the orbit of the polynomial branch; and
- the automorphism subgroup preserving one selected branch.

Formula (1.4) supplies the last group. This is the equivariant content
missing from a purely affine monoid-algebra proof.

### 8.3 Strict-etale descent

The phrase "invariant under strict-etale change" should be replaced by the
normalized-groupoid construction of Section 8.1, or by an equivalent
descent statement for the integral closure in the selected generic
quasi-coherent algebra. This makes the overlap cocycle effective before
the recursive formulas are identified with it.

### 8.4 Discreteness of the boundary Hurwitz datum

"Fully marked fixed-target Hurwitz class is discrete" should mean that the
following are fixed:

1. the stable target tree and all branch labels;
2. source-vertex degrees and branch-cycle/Nielsen classes;
3. allocation of source labels and marked regular-fiber points;
4. source-node partitions and expansion indices;
5. gluing bijections between the two sets of node cycles; and
6. the polynomial component selected on the generic locus.

For this datum, the deformation theorem must identify the relative completed
algebra with only
\[
 s_{\eta j}^{e_{\eta j}}=u_{\eta j}q_\eta
 \tag{8.1}
\]
and finite automorphisms. In particular, it must rule out moving unramified
node preimages or unstable connector moduli not already fixed by the
admissible-cover stability condition.

## 9. Recommended order of work

The efficient continuation is:

1. prove the four obligations in Section 8 in a conventional
   stack-theoretic appendix;
2. record \(\theta_\tau\), not only \(I_\tau\), in the recursive atlas
   certificate;
3. prove log regularity and test simpliciality and
   \(\mathbf Q\)-Gorensteinness cone by cone;
4. independently verify the degree-six formulas (7.2)--(7.4);
5. compute exceptional corrections to the caustic and Maxwell divisor
   classes on the corrected graph; and only then
6. construct the global inertia face category and orbifold Chow product.

This order produces usable geometry after steps 1--4 without making the
general orbifold intersection theory a prerequisite for the foundational
stack comparison.

## References

- Abramovich--Corti--Vistoli,
  [Twisted bundles and admissible covers](https://arxiv.org/abs/math/0106211).
- Abramovich--Caporaso--Payne,
  [The tropicalization of the moduli space of curves](https://arxiv.org/abs/1212.0373).
- Cadman,
  [Using stacks to impose tangency conditions on curves](https://arxiv.org/abs/math/0312349).
- Cavalieri--Markwig--Ranganathan,
  [Tropicalizing the space of admissible covers](https://arxiv.org/abs/1401.4626).
- Gross,
  [Intersection theory on tropicalizations of toroidal embeddings](https://arxiv.org/abs/1510.04604).
- Stacks Project,
  [normalization of algebraic stacks](https://stacks.math.columbia.edu/tag/0GMH)
  and
  [normalization under smooth base change](https://stacks.math.columbia.edu/tag/081J).
- Yasuda,
  [Motivic integration over Deligne--Mumford stacks](https://arxiv.org/abs/math/0312115).
