# Normal quotients and the divisor classes of missing sheets

15 November 2026. Two global constraints sharpen the remaining
class-monodromy problem. A Keller function field contains no nontrivial
finite Galois extension of the target field, with no group restriction.
For a plane Keller map, its finite normalization has a free class group
with the missing boundary divisors as an exact basis. A factorial
normalization therefore forces invertibility. These statements do not
assert preservation of every target Galois group under base change.

## 1. No finite Galois intermediate field (HC4MYIB0)

Let K:A^n->A^n be complex Keller, n>=1, A=C[t1,...,tn],
B=C[z1,...,zn], F=Frac(A), L=Frac(B), with t=K(z). Then

\[
 F\subset E\subset L,\quad E/F\text{ finite Galois}
                       \quad\Longrightarrow\quad E=F.     \tag{1}
\]

All polynomial and extension degrees are unrestricted. The Galois
group in (1) may be nonabelian and nonnilpotent. The proof uses the
same explicit purity and finite-cover inputs as HC4MYHY0, together
with ramification-index multiplicativity. It neither assumes K proper
nor promotes its etale source to a finite cover.

### Every target divisor has at least one unramified visible place

The square-free argument of
[HC4MYHW0](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md#1-no-abelian-intermediate-extension-hc4myhw0)
shows that every target irreducible b has nonconstant square-free
pullback b(K), and every source prime factor h of that pullback maps
dominantly to b=0. Thus A_(b)->B_(h) is a local inclusion of DVRs
and its ramification index is exactly one.

Let C be the finite normalization of A in E. Restrict the valuation
of B_(h) to E; its centre is a height-one prime q of C above (b).
The local inclusions A_(b)->C_q->B_(h) give

\[
          1=e(B_{(h)}/A_{(b)})
            =e(B_{(h)}/C_q)\,e(C_q/A_{(b)}).              \tag{2}
\]

Hence e(C_q/A_(b))=1. The finite residue extension is separable in
characteristic zero, so this place is unramified. The DVR input is
explicitly [Stacks 0EXQ](https://stacks.math.columbia.edu/tag/0EXQ),
Lemma 15.113.3 and the finite-etale criterion in the proof of
Lemma 15.113.8. This is a statement about localizations of the finite
normalization, not about finiteness of the full Keller source.

### Normality of the field extension reaches every missing place

Since E/F is Galois, its group acts transitively on the primes of C
over (b), preserving their ramification indices. For completeness,
one can check transitivity after localizing at (b). If there were
two prime orbits in the semilocal normalization, the Chinese remainder
theorem would give an element a equal to zero modulo one orbit and
one modulo every other orbit. The product of its Galois conjugates
belongs to A_(b); it is zero modulo primes of the first orbit and one
modulo the others, impossible for an element of the local base ring.
Thus all primes over (b) are unramified, including those not visible
in the Keller source. This is the step unavailable for a nonnormal
intermediate extension.

Every codimension-one place of Spec C is therefore unramified.
The target is regular, C is normal, and the morphism is finite;
[purity, Stacks 0BMB](https://stacks.math.columbia.edu/tag/0BMB)
makes Spec C->A^n finite etale. The retained Milne forward
finite-etale/topological-cover comparison and contractibility of C^n
force its connected degree to be one. Hence E=F, proving (1).
The precise retained imports are listed in the
[source record](sources/normal-boundary-import.json).

This gives a general Galois-intermediate obstruction, not unrestricted
linear disjointness. For a finite Galois target field M/F, the
intersection M intersect L can be nonnormal over F. Statement (1)
does not force such an intersection to be F. The nilpotent full-group
result of HC4MYIA0 remains stronger for its specified groups.

## 2. The normalization boundary is an exact free class basis (HC4MYIB1)

Now let K be plane Keller, d=[L:F], and let N=Spec A_L be the
finite normalization of A=C[K1,K2] in L. By the retained
[HC4MYOR2](../hc4-relative-row-orbit/PROOF.md#6-properness-and-the-exact-missing-length),
the source X=A2 is a dense open subvariety of N, by Zariski's Main
Theorem; N->A2 is finite. Let D1,...,Db be the prime divisors in
N minus X. Then

\[
                    Cl(N)=\bigoplus_{i=1}^{b}Z[D_i].       \tag{3}
\]

This is an equality with the displayed basis, not only a bound on
the rank or a list of generators. There can also be omitted points
of codimension two; they do not enter (3). No smoothness of N is
assumed. The statement concerns the global Weil divisor class group,
not local factoriality or local Q-factoriality.

Here is a direct proof. A Weil divisor on N restricts to a principal
divisor on X, since C[s,t] is a UFD. Choose its rational function
h in L. Subtracting div_N(h) leaves a divisor supported on the D_i.
Thus those classes generate Cl(N). If a combination of the D_i is
principal, say sum_i m_i D_i=div_N(h), then div_X(h)=0. UFD
factorization in C(s,t) makes h a unit of C[s,t], hence h in C*.
Its divisor on N is zero as well, so every m_i=0. This proves both
generation and independence. Units of N are also just C*, since
they restrict to source units.

There are no boundary divisors if and only if K is invertible.
For the nontrivial direction, if b=0 then N minus X has codimension
at least two. Normality identifies its ring of regular functions with
that of X, by intersecting the height-one valuation rings in L.
Both N and X are affine; the open immersion therefore identifies
N with X. The original K is now finite etale, and a connected
finite cover of complex affine space has degree one as above.
Consequently

\[
 Cl(N)\text{ torsion (in particular finite or zero)}
       \Longrightarrow b=0\Longrightarrow K\text{ invertible}.       \tag{4}
\]

In particular a hypothetical nonautomorphic plane Keller map must
have at least one missing boundary divisor and a nonzero free global
class group, whose rank equals the number of those divisors. A
globally factorial finite normalization cannot supply either remaining
plane endpoint. A smooth or locally factorial normalization need
not have zero global class group; those hypotheses alone are not
closed by (4). Nor does (3) say that the ambient geometric generic
class group of an output realization has this same basis.

An exact geometric control is the smooth surface N0:xy=z^2-1.
The chart (x,u)->(x,2u+x*u^2,1+x*u) is an open affine plane,
omitting just the component x=0,z=-1. Its inverse is (z-1)/x
where x!=0 and y/(z+1) where z+1!=0. The single missing class
generates Cl(N0)=Z. If both x=0 components are removed instead,
the resulting open is Gm times A1: its nonconstant unit x gives
the diagonal relation among the two boundary classes. This shows
why constant source units, as well as source factoriality, are needed
for independence in (3). The finite projection (x+y,z) on N0
restricts to a map of this affine-plane chart with Jacobian
x-2u-x*u^2, so this control is not a Keller counterexample.

## 3. Class stabilizers may replace nilpotence (HC4MYIB2)

Retain all common-free-basis, generic filtration, affine height-one
geometric-integrality and normality hypotheses of
[HC4MYHW1](ABELIAN_CLASSGROUP_OUTPUT_OBSTRUCTION.md#2-a-class-group-and-geometric-fibre-output-criterion-hc4myhw1),
including its torsion-free geometric generic Weil class group C.
Let the finite Galois image on C be Gamma, with no nilpotence
restriction. It suffices that the stabilizer in Gamma of every class
is **subnormal**: it has a finite chain to Gamma in which each
subgroup is normal in the next. Every closed stable embedding
carrying a plane Keller pair then forces that pair invertible,
in all degrees and for every r>=0.

In fact the proof needs this condition only for the class alpha of
the selected geometric source component, when that class is known.
Let M/F be the finite Galois field of the class action, let J be
the actual source image in Gamma, and put S=Stab_Gamma(alpha).
Then J subset S. If S were proper and subnormal, the last proper
subgroup in a strict subnormal chain would be a proper normal
subgroup N of Gamma containing S. Therefore

\[
                 F\subsetneq M^N\subset M^J
                         =M\cap L\subset L,
\]

would be a nontrivial finite Galois intermediate extension, contrary
to (1). Hence S=Gamma. All d conjugate source-component classes
are alpha, and their reduced principal sum gives d*alpha=0.
Torsion-freeness makes every component principal. The retained
component-pencil, global primitive saturation and smooth minimal-label
proof forces d=1 and the polynomial Keller inverse.

Thus normal class stabilizers suffice; it is not necessary that the
whole class-action group be abelian or nilpotent. No arbitrary
nonabelian stabilizer is silently declared normal. As before, the
actual source-component stabilizer may be nonnormal: the condition
here is on the stabilizer of its **divisor class**.

### An explicit nonnilpotent class module meeting the new condition

The criterion strictly exceeds nilpotent class actions at the level
of integral representations. Consider the rational Hamilton quaternion
division algebra with i^2=j^2=k^2=ijk=-1. Its norm is
q0^2+q1^2+q2^2+q3^2, so a nonzero rational quaternion is invertible.
Put h=(1+i+j+k)/2 and use the lattice

\[
                    H=Z h+Z i+Z j+Z k.                  \tag{5}
\]

The norm-one quaternions in this lattice are the eight signed
coordinate units and the sixteen half-coordinate quaternions
(+/-1 +/-i +/-j +/-k)/2. They form the group of 24 Hurwitz units.
For a direct lattice check, 1=2h-i-j-k and multiplication of the
basis elements in (5) has integral coordinates in that basis;
for example h^2=h-1. Thus H is a ring. The positive norm shows
exactly the stated list of norm-one elements, and multiplicativity
and conjugation show it is a group.

Left multiplication preserves H and fixes no nonzero rational vector
for any nonidentity unit u: (u-1)q=0 in a division algebra implies
q=0. Every nonzero lattice class therefore has trivial stabilizer;
zero has stabilizer the full group. Both stabilizers are normal.
An arbitrary trivial module can be added, or a direct sum of copies
of H can be used, with the same conclusion.

This group is not nilpotent. Let omega=(-1+i+j+k)/2, of order
three. Conjugation by i sends omega to (-1+i-j-k)/2, which is
neither omega nor omega^2; the subgroup of order three is not
normal, so the group is not a direct product of its Sylow subgroups.
Its standard name and structure are recorded in Voight's
[Hurwitz-units discussion](https://link.springer.com/content/pdf/10.1007/978-3-030-56694-4_11.pdf),
retained as primary context. The module argument above is direct.

Consequently an output whose geometric class action is a sum of
these Hurwitz lattice modules and trivial modules satisfies the new
all-degree obstruction, even with nonnilpotent group image. This is
an explicit class-module gate; no particular ambient output realizing
that module is constructed or asserted here.

### Why the remaining root-permutation branch is still open

The S3 action on Z^3/Z(1,1,1) has a nonzero class represented by
(1,0,0). Its stabilizer is the order-two subgroup permuting the
other two positions. That subgroup is not subnormal in S3: it is
contained in no proper normal subgroup. Thus the new criterion
does not silently close arbitrary nonmonic cubic root covers.
Likewise, (1) does not turn a nonnormal cubic intermediate field
into a forbidden Galois one. This is an exact representation-theoretic
boundary, not a Keller counterexample.

The nonmonic Danielewski branch still requires both positive
localized-unit defect and nonnilpotent root monodromy unless another
retained gate applies. Both plane endpoints, global normalization
class groups not ruled out by (3), other output geometries and
unrestricted HC4/JC2 remain open.

The [checker](verify_normal_boundary.py) and
[receipt](receipt-normal-boundary.json) provide exact finite controls
separate from these unrestricted proofs. Purity, finite-cover comparison
and Zariski's Main Theorem are explicit retained imports; they are not
claimed reproved. No formal verification, independent replay, external
review or novelty priority is asserted.
