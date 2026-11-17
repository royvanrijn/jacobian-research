# Fiberwise polynomial Hom and Gaussian field uniformity

15 November 2026. Internal deduction from classical algebraic D-module
operations. Independent specialist review remains outstanding. The
[registry](../../MATH_STATUS.json) is authoritative.

This note closes the **Gaussian complex-field** constructibility/closedness/
uniform-cutoff gap left in the preceding notes, as an internal written proof.
It gives no numerical universal cutoff and implements no general holonomic
stratification producer. The ordinary nonreduced cutoff ideals still need
not stabilize. Previously hash-bound notes retain their historical open wording.

## 1. The classical inputs and the shift convention

We use algebraic left D-modules over C. Write D for holonomic duality and
f_+ for algebraic differential-operator direct image. The shifted inverse
image is

    f^dagger=Lf^*[dim(source)-dim(target)].

The relevant primary source is Hotta--Takeuchi--Tanisaki,
*D-Modules, Perverse Sheaves, and Representation Theory* (2008):
Theorem1.7.3 gives derived base change in a Cartesian square of smooth
varieties, without properness; Theorems2.4.6 and2.7.1 give non-characteristic
ordinary pullback and its commutation with duality; Corollary2.6.15 identifies
global Hom with the direct image of the dual, with a dimension shift.
Algebraic holonomic direct images have bounded holonomic cohomology.
Mustață's Winter2023 notes, Sections6.6--6.7, confirm these conventions and
the generic-connection and stratification statements. Relevant complete
pages and diagrams were read, extracted and visually inspected.

None of these inputs requires that a holonomic module be regular at time
infinity. We do not replace algebraic de Rham cohomology by topological
cohomology of an arbitrary irregular connection. We also do not commute
duality with a nonproper direct image; that would be a different assertion.

## 2. Generic non-characteristic coefficient fibers

Let B be a smooth irreducible finite-type complex variety of dimension d,
X=A1_t times B, and M a holonomic D_X-module. There is a dense Zariski open
U in B such that **every** fiber A1 times {a}, a in U, is non-characteristic
for M, at every finite time.

Here is a direct algebraic proof. Each irreducible component of Ch(M) is a
conic Lagrangian. It is the closure of the conormal bundle to the smooth
part of its irreducible base projection Z. To see this last assertion, the
Liouville one-form vanishes on a conic Lagrangian; at a generic smooth point
its covector therefore annihilates the tangent space of the base projection.
Both conormal and component have dimension dim X, giving equality of their
closures. There are only finitely many components.

The component projecting to all X is the zero section and causes no problem.
Every other Z has dimension at most d. If Z does not dominate B, remove the
closure of its image. If it does, dim Z=d and Z->B is generically finite.
Characteristic zero makes it generically étale. The singular and critical
subsets of Z have dimension at most d-1; remove the closures of their images.
Over the remaining open, every such Z is smooth and maps submersively to B.
A conormal covector to Z with zero dt component must then be zero: its
coefficient covector annihilates a surjective tangent map. This is exactly
non-characteristicity of the coefficient fiber. The finite union of excluded
sets is proper. Nonproperness of Z->B does not invalidate the dimension bound
on their images.

The same proof applies to D M, since duality preserves characteristic support.
No statement about a time-infinity V-filtration is needed here.

## 3. A generic polynomial-Hom formula

Let pi:X->B and i_a:{a}->B. On the open in Section2 put

    M_a=M tensor_(O_B) C(a),   N=pi_+(D_X M).

M_a is the ordinary module restriction. Non-characteristicity gives
Li_a'^* M=M_a and

    D_(A1) M_a = Li_a'^*(D_X M),

where i_a' is the fiber embedding. Derived base change, **with its shift**,
then gives

    i_a^dagger N = pi_(a)+ (D M_a)[-d].

Corollary2.6.15 on the time line gives the generator-independent formula

    RHom_Dt(M_a,C[t]) = i_a^dagger N[d-1].                 (1)

All Hom and Ext here are algebraic global ones. In particular the target is
the polynomial module, not arbitrary analytic functions.

Shrink U once more so that every cohomology module of the bounded holonomic
complex N is an algebraic vector bundle with connection. Holonomic modules
are generically such bundles; for a generic-rank-zero module choose an open
where it vanishes. There are only finitely many cohomology modules. Ordinary
coefficient restriction is then exact on them, and i_a^dagger adds [-d].
Thus (1) gives

    dim Ext_Dt^j(M_a,C[t]) = rank H^(j-1)(N)               (2)

for every a in U. In particular polynomial-Hom dimension is constant there,
and is rank H^(-1)(N).

This is **not** base change for the underived relative polynomial Hom module.
The preceding linear-Gaussian counterexample to that operation remains valid.
We first dualize and integrate the joint module, retain derived restriction,
and use non-characteristic duality on an explicit open.

## 4. Exceptional fibers give a finite stratification

For a smooth locally closed coefficient subvariety S of B, form

    M_S=H^0(Li_S'^* M)=O_(A1 times S) tensor M.             (3)

This is holonomic even when i_S' is characteristic: derived inverse images
preserve holonomicity, and we take their ordinary zeroth cohomology.
For a in S, the ordinary fiber of (3) is **exactly M_a**, by associativity of
ordinary tensor product. We do not assert that (3) equals the whole derived
restriction or that its higher cohomology vanishes.

Apply Sections2--3 to M_S on every smooth irreducible piece of the exceptional
coefficient set. Remove its dense good open and continue on the complement.
Every recursive step lowers coefficient dimension, and each algebraic set
has finitely many irreducible components and a finite smooth decomposition.
This terminates in a finite algebraic stratification. On each stratum (2)
is constant for the appropriate re-restricted module and direct image.

**Fiberwise polynomial-Hom theorem.** For any algebraic holonomic M on
A1 times a smooth finite-type complex coefficient variety B, the function

    a |-> dim_C Hom_Dt(M_a,C[t])                           (4)

is Zariski constructible. In fact every algebraic Ext dimension has a finite
such stratification. Singular coefficient families can be treated by finite
smooth pieces. This is a statement about one **joint holonomic module with
differentiated coefficient directions**; a family holonomic only in t is
not sufficient.

On effectively presented algebraic data this is a theoretical algorithm.
Compute the source and its characteristic support; eliminate the exceptional
coefficient set; compute holonomic dual, direct image and cohomology; find
their generic connection open and ranks; repeat (3) on lower-dimensional
sets. Smooth affine charts or closed embeddings and Kashiwara equivalence
handle the parameter subvarieties. Standard Weyl-module and polynomial-ideal
algorithms provide these operations. No general producer is executed here.

## 5. Application to the exact universal Gaussian source

Fix n,D. On the full phase coefficient space B put

    P_c(X)=sum_(|beta|<=D) c_beta X^beta,
    H=O_B[t,X]/sum_i(partial_i-X_i+t(P_c)_i)O_B[t,X].

The time action is partial_t+P_c; a coefficient vector field delta acts as
delta+t*delta(P_c). These commute with the spatial connection and define
the joint module. It is the top algebraic integration cohomology of a
polynomial exponential connection, hence holonomic in **t and all c**.

Ordinary coefficient restriction is exactly H_(P_a). Indeed the spatial
cokernel and ordinary base change commute by right exactness. This also
holds after restricting to any coefficient subfamily; no unproved uniform
presentation or replacement of normal jets by zero is used.

The [unique-formal-period theorem](GAUSSIAN_FORMAL_COMPLETION_AND_POLYNOMIAL_HOM.md)
proves

    P_a mixed-safe iff dim Hom_Dt(H_(P_a),C[t])=1,

and that dimension is always0 or1. Applying (4) proves constructibility of
the full eventual-mixed-safe locus. Intersect with the closed pure-null
locus to obtain constructibility of the pure-null safe locus S_(n,D).

Now use the retained [closed-cutoff-locus theorem](HOLONOMIC_SOURCES_AND_CUTOFF_LOCI.md).
It gives S_(n,D)=union_(N>=0) Z_N, where Z_N is the closed locus satisfying
every original identity

    E[Q P^m]=0 whenever m>N*deg Q.

Constructibility and uncountability of C force one finite N with

    S_(n,D)=Z_N.                                          (5)

Thus the entire pure-null Gaussian safe field locus is **closed**, and a
single finite integer cutoff works for every safe phase of those fixed
dimension and degree bounds, including exceptional and irregular points.
The radicals sqrt(I_N) stabilize from some point onward. The ordinary ideals
I_N themselves are not claimed to stabilize.

This closes this field-uniformity question as an internal written proof.
It does not say that all pure-null phases are safe; the retained quartic
counterexample remains a failing point. It also gives no new ternary cubic
number, operator order, Laplacian power, or SIC endpoint. Arbitrary GVC/SIC
source encodings still need their own exact-source/target application checked.

## 6. Effective least integer bound and reduced coefficient rings

For a rational finite coefficient family, the stratification algorithm gives
an effective constructible description of its safe locus. The pure-null
ideal is theoretically computable by the retained ordinary-moment algorithm.
Intersect and take the closure: (5) proves it is already closed. Let J be
the resulting radical ideal, relative to the family's coefficient algebra.

For N=0,1,2,... compute the full prescribed-cutoff ideal I_N by
MVLINEARIDEAL20261115 and test sqrt(I_N)=J. Equation (5) proves termination.
The first successful integer is the least common nonnegative integer slope
for that family's pure-null safe field points. On an empty safe locus it
is0. This is a theoretical terminating procedure, not an implemented general
solver, a complexity estimate, an optimal real slope, or a numerical value
for the universal family.

The universal integer from (5) works over **every characteristic-zero field**.
For any individual phase its finitely many coefficients generate a finitely
generated field over Q, which embeds in C. All moment identities and their
zero/nonzero values are preserved. Safety over a field is equivalent to
safety after that embedding, because each multiplier is a finite combination
of monomials.

There is also a uniform reduced-ring consequence. Let A be any reduced
commutative Q-algebra and P in A[X] have degree<=D in n variables. Suppose
every prime residue-field phase is pure-null and safe. The universal forbidden
moments vanish in every residue field, hence lie in every prime ideal of A.
Reducedness makes them zero. The same universal bound therefore holds in A.
Conversely a bound in A descends to every residue field. Thus reduced-ring
pure-null safety is represented by the single reduced closed coefficient
scheme defined by J=sqrt(I_N). This is a field/radical result extended to
reduced targets, not equality or stabilization of all ordinary cutoff ideals.

Nonreduced algebras remain different. On the retained genuine quartic line,
I_N=(e^N), N>=1, with stabilized radical(e). Over Q[e]/e^L the phase is safe
but its marked period has degree L and requires slope at least L. Thus no
uniform extension to all nonreduced algebras is possible, and the reduced
scheme J does not represent their safety functor.

## 7. Checks that distinguish the valid route

For M=O_(t,a) exp(at), dual integration is the point module delta_(a=0).
Its ordinary fiber at a=0 is zero, but its normalized **derived** restriction
is C in degree0. Formula (1), with d=1, recovers the polynomial period1
at that fiber. Using ordinary restriction of the integrated module loses it.
The point module's coordinate action is surjective and has a one-dimensional
kernel; a finite derivative truncation would give the wrong cokernel and is
not used as a cohomology proof.

For the constant family M=O, dual integration is O_B[1]. Its normalized
derived fiber has its nonzero degree d-1, and the final [d-1] in (1) puts
polynomial Hom in degree0. For a constant point module in t, dual integration
is O_B in degree0, so polynomial Hom is zero and Ext^1 has dimension1.

The graph c=t^2 illustrates the exceptional step: its conormal has components
(-2t,1), so the coefficient fiber is characteristic precisely at c=0.
One must re-restrict the source there instead of applying (1) without its
non-characteristic hypothesis.

The central Euler family t*partial_t-a has polynomial solutions at exactly
a in the nonnegative integers. Its extension with a differentiated parameter
and only that equation is not holonomic: its characteristic hypersurface
has dimension3 in a four-dimensional cotangent space. Adding partial_a=0
is inconsistent with the Euler equation, by their commutator. It is therefore
not a counterexample to (4). The same distinction applies to the retained
binomial-family warning, which is not jointly holonomic in its exponent.

The [checker](verify_fiberwise_hom_uniformity.py),
[receipt](fiberwise_hom_uniformity_verification.json) and
[source audit](fiberwise_hom_uniformity_source_audit.json) retain derived-shift,
point-module, non-characteristic, central-parameter and radical/thickness
controls. The proof of (4)--(5) is the written geometric argument, not an
inference from these samples. General source/stratification/cutoff producers
and independent specialist review remain outstanding. The broader research
goal, including structural classification and valid GVC/SIC extensions, stays
active.
