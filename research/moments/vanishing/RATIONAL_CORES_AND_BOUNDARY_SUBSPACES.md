# Rational cores, boundary subspaces, and regularity of actual images

15 November 2026. Internal algebraic synthesis; independent review outstanding.
The [checker](verify_period_boundary_geometry.py) verifies finite exact models,
not the general Weyl-module algorithms or the universal proofs. The previous
[actual-kernel proof](ACTUAL_PERIOD_KERNEL_AND_ALGEBRAIC_IMAGES.md) is unchanged.

Let K be an exactly presented effective finitely generated characteristic-zero
field, D=K[t]<partial_t>, and M a finite holonomic actual period image in
K[[t]], containing 1. Retain independent output variables in K when present.
All structural statements are pointwise. Geometric point-module descriptions
may be made after a finite splitting extension. They descend over K.

## 1. The rational core and the complete punctual obstruction

Define

    R = M intersect K(t),       C = M/K[t],       N = M/R.       (1)

The intersection is taken inside K((t)), with the actual period embedding.
It is a D-submodule, hence finite holonomic. The retained rank-one theorem
therefore gives a unique squarefree divisor d_0 with d_0(0)!=0 such that

    R=K[t,d_0^(-1)].                                           (2)

This assertion does not require M itself to have rank one. It classifies
every rational element occurring inside an arbitrary actual image.

The quotient N is torsion-free over K[t]. Indeed if a nonzero polynomial
p satisfies p*f in R for f in M, then f=(p*f)/p is rational, hence f in R.
Consequently the exact sequence

    0 -> R/K[t] -> C -> N -> 0                                (3)

identifies R/K[t] with the **entire largest point-supported submodule** of C.
Over a splitting field it is

    Gamma_points(C) = direct_sum_(d_0(a)=0) delta_a,             (4)

with multiplicity one at each point, regardless of the rank of M.
Every point-supported submodule maps to zero in the torsion-free N, and
(2) proves the reverse inclusion and the multiplicity statement.
No pole at the formal origin is introduced.

This is a canonical separation of the obstruction:
the rational-pole part (4), and a torsion-free quotient of generic rank
rank(M)-1. Safety is equivalent to both parts being zero. In the retained
Gaussian/full-SIC/GVC formulations use the full amplitude image and their
original-multiplier transfer, not a single selected diagonal series.

It is effective after the actual-kernel construction. In its faithful generic
derivative basis locate the vector 1, form the quotient connection by K(t)*1,
and contract the kernel of M's map to a localized presentation of that quotient.
The kernel equals R by faithfulness. Compute finite generators and apply the
retained rationality/pole algorithm to get d_0. Ordinary Weyl-module kernels,
localization and elimination suffice; this general producer is not run here.

### A quotient need not retain the actual-function embedding

Do not recursively apply (2) to N. For

    M = D*1 + D*log(1-t),

one has R=K[t,(1-t)^(-1)], because the derivative of log(1-t) supplies the
simple pole and all higher poles, whereas any rational multiple of log(1-t)
is nonrational unless zero. The quotient M/R has generic rank one and the
logarithm class is horizontal. It has no faithful embedding into K((t))
compatible with the original functions: its derivative is killed in the
quotient while the original logarithm is nonconstant.

Thus a higher-rank logarithmic extension can have a trivial generic quotient.
The theorem classifies the rational **submodule** of the actual image, not
every rank-one abstract quotient. Splittings and extension classes still matter.

## 2. All remaining lattices over a fixed connection are boundary subspaces

Choose a finite set S of geometric points large enough that the localization

    A=M(*S)

is a meromorphic connection, with its specified generic embedding.
Enlarge S to include every denominator needed to recover a connection basis
from the original generators. Merely clearing the connection matrix and
generator coordinates need not yet guarantee this equality.

There is a unique submodule J of A that restricts to the full connection
outside S and has no nonzero point-supported quotient. It has no
point-supported submodule either, since A is torsion-free. This J is the
intermediate extension j_!* of the connection.

Here is a direct finite-length proof, avoiding an assumption of regularity
or semisimplicity of the connection. Among submodules of A that localize to A,
choose one of least composition length; holonomicity makes that length finite.
A nonzero quotient supported in S would have a smaller kernel with the same
localization, a contradiction. If J_1,J_2 have this property, their intersection
still localizes to A. Both quotients by the intersection are supported in S,
so both vanish. The same argument with J intersect M proves

    J subset M subset A.                                    (5)

The holonomic quotient A/J has finite support. Kashiwara's point equivalence
makes it semisimple:

    A/J = direct_sum_(a in S) delta_a tensor E_a,              (6)

where E_a is a finite-dimensional vector space. Explicitly,
E_a=ker(t-a) on the a-supported summand. Multiplication by t-a is locally
nilpotent on each vector, not necessarily nilpotent with one exponent on
the whole point module. The latter has an infinite derivative tower.

Submodules in (6) are exactly the spaces delta_a tensor U_a for subspaces
U_a subset E_a. Taking inverse images through A -> A/J gives a complete,
inclusion-preserving correspondence

    {J subset L subset A}
       <-> product_(a in S) {linear subspaces U_a subset E_a}. (7)

For prescribed dimensions, these are products of Grassmannians. This is
classification of **embedded extensions of a fixed connection on a fixed
open set**. Abstract isomorphism classes can further identify points under
connection automorphisms; parameter families can change the connection,
the point set and the spaces E_a. Equation (7) asserts no global finite
stratification or constructibility theorem in phase-coefficient space.

Every nonzero quotient L/J in (7) is nonsplit inside L: a splitting would
embed a nonzero point module in the torsion-free A. This does not say that
the connection or J is indecomposable.

### Effective construction without guessing a minimal extension

The largest point-supported quotient of A is computable by holonomic duality:

    Q_max = dual Gamma_S(dual A),
    J = ker(A -> Q_max).                                    (8)

Here Gamma_S is the kernel of localization at S. Duality is exact on the
holonomic category. It takes the largest point submodule of dual A to the
largest point quotient of A. The finite-length argument above shows A/J
is already semisimple and that every point quotient kills J, proving (8).
Alternatively compute finite bases of Hom_D(A,delta_a) and intersect their
kernels; point semisimplicity makes these maps separate A/J.

The polynomial seed lies in J: every morphism from K[t] to a point module
is zero, since its value on 1 would be killed by partial_t, which is
injective on each point-module derivative tower.

[Oaku--Takayama--Tsai2001](https://www.lab.twcu.ac.jp/oaku/OakuTakayamaTsai2001.pdf),
algorithm4.1, gives the holonomic-dual construction; algorithm4.6 computes
the dimensions of holonomic Hom/Ext spaces, and the final enumeration step
recovers explicit Hom bases. Combined with the retained localization and
kernel algorithms, this makes (6)--(8) theoretical algorithms over K.
No dual/local-cohomology/general Hom producer is executed here.

Given M and A,J, its spaces U_a are computed from M/J inside A/J and
taking ker(t-a). The embedding of these spaces, not just their dimensions,
is required to recover M.

## 3. Rational poles and intermediate-extension boundaries are different

The quotient B=M/J in (5) is generally **not** the punctual submodule of
C=M/K[t] in (4). The latter is controlled by R, while B records how M
extends its entire generic connection. The two constructions must not
be substituted for each other.

The logarithmic model in section1 already gives one direction of the
distinction: it has J=M, yet R/K[t]=delta_1. Indeed M is generated by the
logarithm alone, since (t-1)*partial_t log(1-t)=1. Its annihilator contains
partial_t*(t-1)*partial_t, which acts on partial_t^k delta_1 as
-(k+1)*partial_t^(k+1) delta_1. This action is injective, so no nonzero
point quotient at1 exists. There are no other singular points on A1.
Thus a minimal extension itself can contain a rational-pole submodule;
its quotient by the polynomial seed can have point torsion.

A concrete algebraic-function model proves the distinction. Put

    u=sqrt(1-t),        g=(1+u)/(t-2),
    M=D*1 + D*u + D*g,
    J=K[t] direct_sum K[t,(1-t)^(-1)]*u.                     (9)

Choose the branch u(0)=1, so every generator is a formal power series at zero.
The generic rank is two. Localize also at t=2 to obtain A.
The u-connection has exponent one-half at t=1 and has no point quotient
there; J is its intermediate extension together with K[t].

At t=2, A/J has two point directions, generated by the simple-pole classes
of 1/(t-2) and u/(t-2). The class of g generates their diagonal line.
It follows from (7), or directly from the identical point-module relations,
that

    M/J = delta_2,          R=K[t],          Gamma_points(C)=0. (10)

Indeed a nonzero pure rational-pole direction at2 is not in the diagonal
line. No rational pole at1 occurs, since the u-branch is independent of 1.
Thus M has a nonzero boundary quotient but no rational-pole submodule
beyond K[t].

The finite-cover trace projection u -> -u gives

    (g+conjugate(g))/2 = 1/(t-2),

which is not in M. A rational projection of the generic algebraic connection
need not preserve the actual polynomial lattice. Even when the generic
connection has finite monodromy and splits, its intermediate lattices can
glue the summands at a point. This is an abstract period-module control,
not a new Gaussian or operator counterexample.

## 4. Regularity of every actual generic connection is decidable

Let

    L_min = partial_t^d + sum_(j<d) a_j(t) partial_t^j

be the exact minimal equation of the cyclic actual period constructed in the
previous note. The actual generic connection is regular singular on P1 iff
the classical scalar Fuchs inequalities hold:

    ord_a a_j >= j-d             at every finite point a,
    deg(num a_j)-deg(den a_j) <= j-d       at infinity.       (11)

Ignore zero coefficients. Only finitely many rational-function poles need
testing. Factorization or valuations in irreducible root fields decide them.
The infinity inequality follows by z=1/t and the triangular derivative
change, or by the Euler-operator form. A strict violation gives an actual
irregular singularity, not merely an irregular direction of a larger source.

Minimality is essential. The Kummer polynomial H_n has a regular first-order
minimal equation D-H_n'/H_n, while its usual order-two input equation is
irregular at infinity. Testing the input source would therefore misclassify
the selected actual image. Conversely a rank-two actual image generated by
1 and exp(t) is irregular at infinity; a logarithmic extension is regular
singular but has infinite monodromy; an algebraic square-root image is regular
with finite geometric monodromy. Rank alone distinguishes none of these.

The point-module modifications (7) do not change generic formal slopes or
the regular/irregular connection type. On a curve punctual holonomic modules
are regular, so the intermediate-extension and lattice changes preserve
regularity of the holonomic module. Include infinity in (11), rather than
testing only affine finite poles.

For number-field scalar images, combine (11) with the previously imported
complete image-algebraicity decision:
polynomial safety; rational finite-pole failure; nonrational regular
finite-monodromy/algebraic failure; nonrational regular transcendental failure;
or irregular failure. The known bounded-polynomial-output reduction permits
the corresponding algebraicity test in its exact retained scope. Regularity
itself is decided over the effective finitely generated coefficient field.
This does not import a complete algebraicity test for arbitrary
transcendental coefficient fields.

These are intrinsic categories of actual images and effective certificates,
not a closed-form classification of all Gaussian phases realizing them.
No new large actual module is computed here.

## 5. Why the parameter-uniformity gap is still open

The joint coefficient/time/source exponential module has more structure
than the one-direction family (1+t)^a, which is not jointly D-finite in a.
A possible uniformity route is to control a good time-infinity filtration
after coefficient restriction, including turning points and generator shifts.

The newly inspected
[Steiner,2406.14253v1](https://arxiv.org/html/2406.14253v1),
section5, treats associated-graded compatibility under non-characteristic
restriction with specified positive-weight hypotheses; theorem5.6 takes
all weights equal and positive. The time-infinity filtration needed here
assigns coefficient directions weight zero. The stated theorem does not
supply that restriction step, a finite valid chart cover, or bounds on
exceptional period kernels. The HTML body date24August2026 is not a new
arXiv version; the sole listed submission is20June2024.

No global constructibility/uniform-cutoff theorem is claimed from this
source. The separate Artin family still excludes a uniform extension across
all nonreduced coefficient algebras. Global complex-field safe-locus
stabilization remains open.

The next structural questions are now separated: classify generic higher-rank
Gaussian connections; compute which boundary subspaces Gaussian periods select;
classify the normalized constant-b higher seed-rank trivial quotients; and
control exceptional coefficient restrictions. Independent review of these
algebraic arguments and the older analytic exclusions is still required.

[Source and assurance audit](period_boundary_geometry_source_audit.json).

