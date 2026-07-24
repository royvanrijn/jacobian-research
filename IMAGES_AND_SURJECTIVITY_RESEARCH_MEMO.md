# Images and surjectivity of weighted Keller maps

**Research memo — 24 July 2026**

## Executive conclusion

For the weighted seed family, the complex-image problem is essentially
solved by the inverse pencil

\[
E_{s,t}(W)=H(W)-sW+t.
\]

If \(H\) has degree \(N\), the associated polynomial map
\(G_H:\mathbb A^3\to\mathbb A^3\) is étale, has geometric degree \(N\), and,
on the target chart \(C\ne0\), its finite preimages are in bijection with the
simple roots of

\[
E_{A,B,C}(W)=H(W)-BCW+cAC^2.
\]

Consequently:

1. \(G_H\) is surjective over \(\mathbb C\) if and only if every member of the
   affine pencil \(H-sW+t\) has a simple root.
2. Every nonsurjective seed has a **unique** omitted pencil value
   \(\omega=(s,t)\).
3. The corresponding omitted target locus is not a point but the closed
   curve
   \[
   \Gamma_\omega=
   \left\{\left(\frac{t}{cC^2},\frac{s}{C},C\right):C\in\mathbb C^*\right\}
   \cong\mathbb G_m.
   \]
4. The irreducible components of the nonsurjective seed locus are indexed by
   \[
   2a+3b=N,
   \]
   have dimension \(a+b-1\), and possess explicit smooth normalizations.
5. Hence
   \[
   \dim\mathcal N_N=\left\lfloor\frac N2\right\rfloor-1,\qquad
   \operatorname{codim}_{\mathcal A_N}\mathcal N_N
   =\left\lfloor\frac{N-3}{2}\right\rfloor.
   \]
   Quartics are generically nonsurjective, while generic seeds of every
   degree \(N\ge5\) are surjective.

The cleanest publication is therefore not merely a list of questions. It is
a theorem paper giving the complete complex-image and exceptional-seed
classification, followed by a second paper on real images and low coordinate
degree.

## 1. Status of the construction

The public MathOverflow discussion gives explicit noninjective Keller maps in
every geometric degree \(N\ge3\), and a follow-up image calculation gives
surjective, étale, noninjective examples for every \(N\ge5\):

- [every geometric degree \(N\ge3\)](https://mathoverflow.net/questions/513440/geometric-degrees-of-counterexamples-to-the-jacobian-conjecture-in-dimension-thr/513470);
- [surjectivity for every \(N\ge5\)](https://mathoverflow.net/questions/513390/a-jacobian-neutral-birational-construction-for-keller-maps-with-noninjective-ex/513399).

The underlying three-dimensional counterexample has now also received an
independent Isabelle/HOL verification of its constant Jacobian and explicit
multiple fiber; this verifies the foundational phenomenon, not the
all-degree image classification
([Archive of Formal Proofs](https://isa-afp.org/entries/Jacobian_Counterexample.html)).

The public record is extremely recent. The weighted-family classification in
this repository should therefore be presented with complete proofs and
machine-checkable certificates, rather than relying on the words “verified”
or “generic” without a precise parameter space.

## 2. Exact surjectivity criterion

Let \(\mathcal A_N\) be the normalized admissible seed space

\[
\deg H=N,\quad
H(0)=H'(0)=H(1)=0,\quad H'(1)=-1,\quad H''(1)\ne-2.
\]

It is an open subset of affine \((N-3)\)-space.

Define

\[
\Omega_H=\{(s,t)\in\mathbb C^2:
H(W)-sW+t\text{ has no simple root}\}.
\]

The exact image formula is

\[
G_H(\mathbb C^3)
=\mathbb C^3\setminus
\left\{(A,B,C):C\ne0,\ (BC,cAC^2)\in\Omega_H\right\}.
\]

Thus

\[
\boxed{G_H\text{ is surjective}\iff\Omega_H=\varnothing.}
\]

This is algorithmic. For each partition

\[
\lambda=(\lambda_1,\ldots,\lambda_r)\vdash N,\qquad \lambda_i\ge2,
\]

solve the finite coefficient-matching system

\[
H(W)-sW+t=h_N\prod_{i=1}^r(W-r_i)^{\lambda_i}
\]

in coefficients \(W^{N-1},\ldots,W^2\). The remaining two coefficients
recover \(s,t\). This is an exact necessary-and-sufficient test, not a
discriminant heuristic.

The stronger uniqueness statement is:

\[
\boxed{|\Omega_H|\le1.}
\]

Indeed, two distinct omitted polynomials \(P,Q\) would have all root
multiplicities at least two and \(P-Q\) affine-linear. Polynomial
Mason--Stothers excludes this except for the even all-double equality case;
there \(P=cA^2,Q=cB^2\), and
\((A-B)(A+B)\) cannot have degree at most one unless \(A=B\).

### Geometric interpretation

The discriminant curve of the pencil is rational, with normalization

\[
r\longmapsto \bigl(H'(r),\,rH'(r)-H(r)\bigr).
\]

An omitted pencil value is exactly a point at which the total contact
multiplicity of all normalization branches exhausts \(N\). Ordinary cusps
contribute \(3\), ordinary bitangent nodes contribute \(2+2\), and higher
contacts give the general partition. Generic discriminant geometry therefore
immediately explains the phase change:

- \(N=3\): a cusp can exhaust the degree;
- \(N=4\): a node can exhaust the degree;
- \(N\ge5\): an ordinary cusp or node leaves at least one simple root.

## 3. Irreducible components of the nonsurjective seed locus

The primitive contact multiplicities are \(2\) and \(3\), because these are
the indecomposable elements of the additive semigroup
\(\{2,3,4,\ldots\}\). Every multiplicity \(m\ge4\) is a collision of smaller
allowed contacts.

It follows that the irreducible components are indexed by

\[
\mathcal C_{a,b},\qquad 2a+3b=N.
\]

Their normalization has a particularly compact presentation. Let \(Q,R\) be
monic of degrees \(a,b\), put

\[
M=Q^2R^3,\qquad
\Phi=M(1)-M(0)-M'(0),\qquad
D=M'(1)-M'(0),
\]

and localize away from \(D=0\) and the weighted admissibility divisor. Then

\[
\widetilde{\mathcal C}_{a,b}=V(\Phi)
\longrightarrow\mathcal C_{a,b},
\qquad
H=\frac{-M+M'(0)W+M(0)}D
\]

is finite and birational. The source is smooth and integral, so this is the
normalization morphism.

This proves componentwise irreducibility and gives

\[
\dim\mathcal C_{a,b}=a+b-1.
\]

Maximizing under \(2a+3b=N\) gives the dimension and codimension formulas in
the executive conclusion. It also gives the component count

\[
\#\{\mathcal C_{a,b}\}=[x^N]\frac1{(1-x^2)(1-x^3)},
\]

a period-six quasipolynomial.

The local sources for a paper-ready proof are:

- [fixed-seed omitted-value classifier](extended-geometry/OMITTED_VALUE_CLASSIFICATION.md);
- [unique omitted-value lemma](extended-geometry/UNIQUE_OMITTED_VALUE.md);
- [contact-atom theorem and phase diagram](extended-geometry/CONTACT_ATOM_PRINCIPLE.md);
- [collision closure and irreducibility](extended-geometry/COINCIDENT_ROOT_REBUILD.md);
- [explicit component normalizations](extended-geometry/COMPONENT_NORMALIZATION.md).

The corresponding symbolic checks all pass in the current workspace.

## 4. What “every finite omitted set” can mean

There are two different notions.

### Omitted pencil values

For a fixed weighted seed, \(\Omega_H\subset\mathbb C^2\) has cardinality at
most one. Therefore arbitrary finite subsets do **not** occur:

\[
\Omega_H=\varnothing\quad\text{or}\quad\Omega_H=\{\omega\}.
\]

This question is already settled negatively inside the weighted family.

### Omitted target points

When \(\Omega_H=\{\omega\}\), the omitted locus in \(\mathbb C^3\) is
\(\Gamma_\omega\cong\mathbb C^*\), not a finite set. Thus no nonempty finite
target complement occurs in this family either.

A genuinely new question is:

> Can a polynomial étale self-map of \(\mathbb A^3_\mathbb C\) have image
> complement equal to a prescribed nonempty finite set?

That is a problem about all nonproper étale maps, not about the present
weighted pencils. It should be separated from the seed-classification paper.
The general nonproperness set is much larger than the omitted set: Jelonek's
theorem says that the nonproperness set of a dominant generically finite
polynomial self-map is empty or a uniruled hypersurface
([Jelonek 1993](https://doi.org/10.4064/ap-58-3-259-266)). Missing points may
sit on this hypersurface without forming its components.

## 5. Topology of the complex image

For a surjective seed the complex image is \(\mathbb C^3\), hence
contractible.

For a nonsurjective seed,

\[
X_H=G_H(\mathbb C^3)=\mathbb C^3\setminus\Gamma_\omega,
\qquad \Gamma_\omega\cong\mathbb C^*.
\]

This already determines basic topology. Since \(\Gamma_\omega\) is a smooth
real two-manifold embedded with real codimension four, general position gives

\[
\pi_1(X_H)=0.
\]

Alexander duality with compact supports gives

\[
\widetilde H_i(X_H;\mathbb Z)
\cong H_c^{5-i}(\mathbb C^*;\mathbb Z),
\]

and hence

\[
H_3(X_H;\mathbb Z)\cong\mathbb Z,\qquad
H_4(X_H;\mathbb Z)\cong\mathbb Z,
\]

with all other reduced homology groups zero.

The next worthwhile topology question is not homology but the embedded
isotopy/homotopy type:

> Is every weighted omitted curve \(\Gamma_\omega\subset\mathbb C^3\)
> ambiently unknotted, and is its complement homotopy equivalent to
> \(S^3\vee S^4\)?

The displayed parametrization makes this plausible, but it needs an explicit
ambient isotopy or a complement calculation; Alexander duality alone does
not prove the wedge statement.

## 6. Real images

For a real seed, on \(C\ne0\) the change of target coordinates

\[
(A,B,C)\longleftrightarrow(s,t,C)=(BC,cAC^2,C)
\]

is a real-analytic diffeomorphism. Therefore

\[
G_H(\mathbb R^3)\cap\{C\ne0\}
\cong
\{(s,t):H-sW+t\text{ has a simple real root}\}\times\mathbb R^*.
\]

The entire plane \(C=0\) is in the image. This is the correct starting point
for real topology.

Complex surjectivity does **not** imply real surjectivity:

- if \(N\) is even, open chambers in the \((s,t)\)-plane may contain
  polynomials with no real roots;
- if \(N\) is odd, every polynomial has a real root off the discriminant, so
  every regular pencil value has a simple real root, but points on the
  discriminant may have only repeated real roots.

Thus the real problem reduces to a planar chamber decomposition cut out by
the parametrized discriminant

\[
r\mapsto(H'(r),rH'(r)-H(r)),
\]

with each chamber labelled by its number of simple real roots. The
three-dimensional real image is obtained by taking two products with
\(\mathbb R_{>0}\) and attaching the common plane \(C=0\) through the explicit
boundary charts.

This reduction is strong enough for an algorithmic classification degree by
degree using real-root counts, subresultant sign conditions, and cylindrical
algebraic decomposition. A uniform all-degree classification will require
controlling the real cusp/node arrangement of the discriminant curve.

## 7. Minimum coordinate degree

For a weighted seed of inverse/geometric degree \(N\), the standard lift has
ordinary coordinate-degree profile

\[
\bigl(5N-8,\ 5N-9,\ 4\bigr).
\]

Since the first surjective geometric degree is \(N=5\), this gives the
explicit upper bound

\[
\boxed{\min\max_i\deg F_i\le17}
\]

for a surjective, noninjective Keller self-map of \(\mathbb A^3_\mathbb C\).
The degree-five member has profile \((17,16,4)\).

Wang's theorem that every quadratic Keller map is invertible gives only the
general lower bound

\[
\boxed{\min\max_i\deg F_i\ge3}.
\]

The interval \(3\le D\le16\) is open on the present evidence. The right search
problem is:

1. search the full degree-five seed family after low-degree source and target
   automorphisms, not just the normalized weighted presentation;
2. minimize degree using cancellations in the birational frame before
   expanding;
3. impose surjectivity through the no-full-contact condition, rather than
   testing target samples;
4. certify candidates by exact Jacobian, generic-degree, and simple-root
   reconstruction calculations.

Coordinate degree is not invariant under arbitrary polynomial left-right
equivalence, so a lower-bound theorem will need an intrinsic complexity
measure or a proof covering all degree-\(D\) presentations.

## 8. Recommended publication order

### Paper I: Images of weighted noninvertible étale maps

Prove in one sequence:

1. polynomiality, constant Jacobian, and generic degree \(N\);
2. marked-simple-root reconstruction;
3. the exact image formula;
4. uniqueness of the omitted pencil value;
5. full-contact stratification;
6. \(2/3\)-component classification, irreducibility, and normalization;
7. the dimension/codimension phase diagram;
8. complex-image homology.

This paper already has a coherent theorem-level story and should be the first
priority.

### Paper II: Real images and low-degree optimization

Develop:

1. the discriminant chamber decomposition for real seeds;
2. topology and numbers of connected components of real images and fibers;
3. explicit degree-five and degree-six atlases;
4. coordinate-degree minimization and lower bounds.

### Separate general problem

Treat prescribed finite complements for arbitrary polynomial étale maps as a
broader problem. It is not an unresolved corner of the weighted
classification.

## 9. Verification commands

The four load-bearing checks used in this memo are:

```bash
.venv/bin/python scripts/verify_omitted_value_classification.py
.venv/bin/python scripts/verify_unique_omitted_value.py
.venv/bin/python scripts/verify_contact_atom_principle.py
.venv/bin/python scripts/verify_component_normalization.py
```

All four pass in the current workspace.

