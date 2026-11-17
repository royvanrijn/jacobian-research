# Quantum division indices and what completion forgets

14 November 2026. **QDI1.** Every full polynomial-coefficient formal
quantization of Long has completed quantum division degree three.
If the lift has one common spatial degree bound D, its uncompleted
left and right division indices are finite, at least three, and at
most D^4. The completed degree does not, by itself, determine those
finite indices.

The completion result and the bounds below are self-contained valued
linear algebra and PBW growth arguments. Equality of the uncompleted
left and right indices uses
[Lauritzen--Thomsen, Corollary 8.11](https://arxiv.org/html/2610.15469v2#S8.SS7).
Their recent preprint is used as an external theorem, not independently
reviewed here. The finite-generation warning below uses their earlier
[finite-image-module theorem](https://arxiv.org/html/2308.09384v2).

There is an exact full polynomial control: QNI1's rank-three identity
lifts have uncompleted degree three and completed degree one, with
forward spatial bound twelve and arbitrarily high identity contact.
Thus a comparison theorem needs additional hypotheses even for full
polynomial maps satisfying every Weyl relation. A separate rank-one
rational control has completed degree three and arbitrarily larger
generic degree.

These statements supply an all-degree restriction and identify a
specific missing descent condition. They construct no bounded Long
lift and exclude no hypothetical lift of division degree three or more.

## Valuations on the quantum rational function fields

Let K have characteristic zero, F=K((h)), and let A be the finite PBW
algebra over F on n canonical pairs, with
[p_i,x_j]=h*delta_ij and the other generator commutators zero.
Write Q=Frac(A) for its Ore division ring.

For a finite PBW sum define v as the minimum h valuation of its
scalar coefficients. It is multiplicative: normalized nonzero factors
reduce to nonzero commutative polynomials, because every commutator
is divisible by h. Hence v extends to Q, with value group Z, central
uniformizer h, and residue field
\[
                   \overline Q=K(x_1,\ldots,x_n,p_1,\ldots,p_n).
                                                               \tag{1}
\]
Indeed every rational function in the commuting residue variables is
the residue of a fraction of finite PBW polynomials. Let Qhat be the
valuation completion; it remains a division ring, with this residue
field. The usual geometric-series inverse proves the division claim.

Let psi be a polynomial symplectic endomorphism of affine 2n-space.
Its Jacobian is a unit, so its coordinate functions are algebraically
independent and the induced function-field extension is finite. Put
\[
 \delta=[K(x,p):K(\psi(x,p))].
\]
Let Phi be any full formal lift of psi, with polynomial coefficients
at each h order, as in [QFL1](PROOF.md). Its images belong to Qhat even
if their spatial degrees are unbounded. A finite source PBW polynomial
with its least scalar valuation removed has a nonzero classical
polynomial reduction; its image reduces to that polynomial evaluated
at psi, which is nonzero. Thus substitution preserves valuation.

It follows that Phi embeds the source Weyl division ring into Qhat
and extends isometrically to its completion. Denote the resulting
complete image by Ehat. Its residue field is exactly K(psi(x,p)).
When Phi has common finite spatial support, its images lie in A
over F; write E=Phi(Q) inside Q. In that case Ehat is also the
completion of this uncompleted image division ring. Without finite
support we do not claim E is a subring of Q.

## A complete extension has the residue degree

Here is the general valued-division-ring statement used throughout.

Suppose Ehat is a complete valued division subring of a complete
valued division ring Qhat. Both have value group Z with the same
central uniformizer h, and both residue division rings are fields.
If their residue extension has finite degree delta, then
\[
           [\widehat Q:\widehat E]_\ell
             =[\widehat Q:\widehat E]_r=\delta.                 \tag{2}
\]
Choose lifts b_1,...,b_delta of a residue-field basis.
They form a basis on either side. More precisely,
\[
 v\left(\sum_i b_i e_i\right)=\min_i v(e_i),\qquad
 v\left(\sum_i e_i b_i\right)=\min_i v(e_i).                    \tag{3}
\]
For a nonzero sum, remove its minimum h power and reduce. Residue
linear independence forbids cancellation. This proves (3) and
independence.

For spanning, first multiply an element by a power of h to give
nonnegative valuation. Expand its residue in the chosen basis with
coefficients in the residue of Ehat, and lift those coefficients
to elements of Ehat of nonnegative valuation. Subtract this finite
sum; the remainder has valuation at least one. Repeat after division
by h. For each basis element the successive coefficients converge
in Ehat. Completeness and (3) give the desired expansion. The same
construction on the other side proves both assertions.

No assertion about algebraicity of an arbitrary formal operator
coefficient is needed. This is a theorem about the completed
division rings, not about finite differential orders.

Consequently every full formal lift of psi has completed division
degree delta, without a spatial bound. Formal Hamiltonian gauge
equivalence alone cannot distinguish these completed degrees.

## Finite PBW images give finite uncompleted indices

Suppose now that Phi has common spatial degree bound D. Its finite
PBW generic map has generator degrees at most D, so
Phi(A_{\le t}) is contained in A_{\le Dt} for the Bernstein
filtration. Both uncompleted indices obey
\[
       [Q:E]_\ell\le D^{2n},\qquad [Q:E]_r\le D^{2n}.             \tag{4}
\]
Here is a direct proof, independent of the new index-symmetry theorem.

Take m right E-linearly independent elements of Q. A common left
denominator puts all of them in A and preserves independence.
Call them b_i, and bound their Bernstein degrees by M. The
K((h))-linear spaces b_i*Phi(A_{\le t}) form a direct sum in
A_{\le M+Dt}. Therefore
\[
 m {t+2n\choose 2n}\le {M+Dt+2n\choose 2n}.
\]
Divide by t^{2n} and let t tend to infinity; this gives m<=D^{2n}.
An arbitrarily large independent set is impossible. For left
independence, use a common right denominator and the spaces
Phi(A_{\le t})*b_i. This proves (4) on both sides.

The same residue independence argument as (3), using basis lifts
chosen in Q, gives
\[
                    \delta\le[Q:E]_\ell,\ [Q:E]_r.             \tag{5}
\]
Lauritzen--Thomsen's cited Corollary 8.11 says that Weyl division
rings in characteristic zero are index-symmetric. With that
external theorem, the two finite indices agree; denote them by N:
\[
                         \delta\le N\le D^{2n}.                \tag{6}
\]
Finiteness, the individual bounds, and the completed-degree formula
do not depend on index symmetry.

## The exact condition for equality

Let O_Q and O_E be the valuation rings inside Q and E, respectively;
these are the elements of nonnegative h valuation. They are not
the polynomial PBW algebras. In particular, they contain fractions
whose reductions are arbitrary rational functions.

The following conditions are equivalent:
\[
\begin{array}{ll}
{\rm (i)}&[Q:E]_r=\delta;\\
{\rm (ii)}&O_Q\text{ is finitely generated as a right }O_E
                       \text{-module};\\
{\rm (iii)}&\text{lifts of any residue basis form a right }E
                       \text{-basis of }Q.
\end{array}                                                    \tag{7}
\]
The corresponding left statements hold independently.

To prove (i) implies (iii), residue basis lifts are independent
by (3), so their number equals the dimension. Their valuation
formula also gives
O_Q=sum_i b_i O_E directly, proving (ii). Conversely O_E is a
local discrete valuation ring, possibly noncommutative, with
maximal ideal h O_E: a valuation-zero element has its inverse
inside O_E, and every nonzero element is h^k times a unit.
Moreover
\[
                    O_Q/hO_Q=\overline Q
\]
has dimension delta over O_E/hO_E. If O_Q is finitely generated,
Nakayama's lemma makes any lifts of a residue basis generate O_Q.
Inverting h makes them span Q over E, while (3) gives independence.
This proves (7). One can prove the Nakayama step directly by the
usual minimal-generator argument, since 1-hu is a unit in O_E.

There are also surjective linear maps
\[
 Q\otimes_E\widehat E\longrightarrow\widehat Q,\quad q\otimes e
                  \longmapsto qe,\qquad
 \widehat E\otimes_E Q\longrightarrow\widehat Q,\quad e\otimes q
                  \longmapsto eq.                            \tag{8}
\]
Surjectivity follows from the lifted residue basis in Q.
Their kernel dimensions are [Q:E]_r-delta and [Q:E]_ell-delta
over Ehat on the appropriate side. The tensors in (8) are vector
spaces; no tensor-product ring structure is asserted.

Thus the missing equality is precisely finite generation of this
specific valuation lattice. Completion can forget nonzero generic
dimensions. It is not enough to know that a different algebraic
module is finite, projective, or formally free.

## Application to every full Long lift

Use the actual retained coordinates of [QPA1](RAMIFICATION_AND_POSITION_AVOIDANCE.md)
and [QEI1](ELLIPTIC_ISOGENY_AND_QUANTUM_ENERGY.md):
\[
 r=R,\quad v=1/x,\quad a=v+Y,\quad d=D,\qquad
 s=2v+4a-3ra^2,\quad t=(ra^3-a^2-va)/2.
\]
The four variables r,v,a,d are birational coordinates in the
ambient classical function field. Explicitly,
\[
\begin{split}
 x&=1/v,&q&=(2v-rv^2)/3,&Y&=a-v,\\
 \beta&=5v^2-rv^3-3av,&M&=\beta+9q^2,&D_0&=d-H(x,Y,\beta),\\
 p&=3q^2 M+2(1-3xq)D_0,&
 z&=(1+3xq)M/2-3x^2D_0 .
\end{split}                                                     \tag{9}
\]
Here H is Long's retained correction, so d=D_0+H.
The two linear equations for p,z have determinant -1; hence
(9) is a rational inverse, not a claimed canonical coordinate
change. In particular d is independent of r,v,a.

We have
\[
 v=(s-4a+3ra^2)/2,\qquad
                  ra^3-2a^2+sa+4t=0.                         \tag{10}
\]
Over K(r,s,d), the output t is a polynomial of degree three
in the transcendental a. Thus the full classical function-field
degree, including the fourth image, is exactly three. Adding
the independent d does not change it.

Equations (2), (4), and (5) prove:
\[
 [\widehat Q:\widehat E]_\ell
  =[\widehat Q:\widehat E]_r=3,\qquad
 3\le[Q:E]_\ell,\ [Q:E]_r\le D^4.                             \tag{11}
\]
The first conclusion holds for every full formal Long lift.
The bounds apply whenever it has common spatial bound D.
There is no exact-R, positional-first-image, scalar-order,
ordering, or partner-degree hypothesis in (11).
With the cited index-symmetry theorem both finite indices are N.

This rules out a birational quantum division realization of Long,
but not a finite realization with N>=3. The statement controls
fraction division fields; it is weaker than surjectivity of the
finite Weyl algebra.

## What bounded canonical changes preserve

Precomposing a finite Weyl map with a Weyl automorphism leaves
its image division subring unchanged. Postcomposing transports
the extension by an ambient division automorphism. Hence the
individual uncompleted indices are invariant under arbitrary
finite source and target Weyl automorphisms.

For changes regular at h=0 whose inverses are regular and whose
classical reductions are polynomial canonical automorphisms,
the completed degree delta is preserved as well. Thus the
differences [Q:E]_r-delta and [Q:E]_ell-delta are invariants
under these changes.

Adding identity canonical pairs preserves both indices. For
example, adjoin spectator Weyl generators commuting with Q.
A right E-basis of Q remains independent over the spectator
extension E': clear a common denominator in E', expand in
spectator PBW monomials, and use independence in Q at each
coefficient. Its span over E' is stable under left multiplication
by Q and all spectator generators. Every nonzero such operator
acts injectively on this finite-dimensional right E'-space and
hence surjectively, so adjoining inverses does not enlarge it.
It spans the whole extended division ring. The left proof is
the same. The classical degree is preserved by adjoining
independent variables, so the differences are stable too.

These invariants are not invariant under arbitrary formal
Hamiltonian gauges. The next full example demonstrates that
distinction with every relation intact.

## A full polynomial control: completed degree one, generic degree three

We first record a general differential-field fact. Let L/B be
a finite separable function-field extension of degree m in
characteristic zero. Suppose commuting derivations on B extend
to L and generate the same ambient rational differential-operator
field. Then the extension of quotient differential-operator
division rings has left and right degree m.

Indeed choose a field basis b_i of L/B. Their right span over
the smaller operator division ring is preserved by L and all
derivations, using
partial_j*b_i=b_i*partial_j+partial_j(b_i).
Every nonzero operator preserves this finite-dimensional space
injectively, hence surjectively, so its inverse also preserves
it. Thus the basis spans the ambient division ring. For
independence, clear a common denominator in the smaller operator
ring and compare PBW derivative coefficients, using the original
field basis. The left argument is identical. This proof applies
to any number of commuting derivations and to nonzero scalar h
multiplying them.

Apply this to the cotangent lifts in
[QNI1](IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md).
For its retained three-variable map put u=1+xy and denote
the unnormalized outputs by A,B,C. On x!=0 set
w=y+1/x and rho=2/x. Then
\[
 A=w^2+\rho w/2-Cw^3,\qquad B=\rho+4w-3Cw^2,
\]
so
\[
 \rho=B-4w+3Cw^2,\qquad Cw^3-2w^2+Bw-2A=0.                  \tag{12}
\]
Over K(B,C), A is a cubic polynomial in the independent w.
The function-field degree is exactly three. The output
normalization (C/2,B,A) and scalar source and target dilations
do not change it.

QNI1 sets f_t=t^-1*f(tx), t=h^N for any N>=1, and takes its
inverse-Jacobian cotangent Weyl operators. It has all full
rank-three relations, common forward degree twelve, and
classical reduction the identity. Its generic division degree
is three by the differential-field fact, since nonzero t is
invertible in K((h)). Its completed degree is one by (2).
Thus each kernel in (8) has dimension two, and neither
valuation lattice in (7) is finite.

This is a full polynomial countercontrol to the inference
"bounded forward images plus exact full relations imply equality
of completed and generic division degrees." It reuses an
established higher-rank example, not a new Dixmier counterexample.
It supplies no such counterexample in rank two. In particular,
its left-free graph module is a different object from O_Q as
an O_E-module.

## A rank-one rational control with arbitrary extra degree

Let F=K((h)), Q=D_1(F), and use [p,x]=h.
For integers m>d>=1 set
\[
 f=x^d+hx^m,\qquad g=\partial_x f,\qquad
                       U=f,\quad V=g^{-1}p.                  \tag{13}
\]
Then [V,U]=h exactly. Both images are finite rational-coefficient
differential operators and are regular for the h Gauss valuation;
g reduces to d*x^{d-1}, a nonzero element of the residue field.

Over F, [F(x):F(f)]=m. Replacing the differential generator p by
g^-1*p and applying the preceding differential-field argument
gives uncompleted left and right division indices m. The
classical map is
\[
                     (x,p)\longmapsto
                  (x^d,p/(d x^{d-1})),
\]
of function-field degree d. Its completed degree is d by (2).
Thus d=3 gives completed degree three and every generic degree
m>3. The differences m-3 are unbounded.

The ordinary polynomial equation hX^m+X^d-u=0 explains where
the extra branches lie. Over the h-valued coefficient field
with u of valuation zero, its lower Newton polygon has a
horizontal segment of length d and a segment of slope
1/(m-d), length m-d. The latter roots have negative valuation
-1/(m-d). Completion at the selected Gauss valuation retains
the degree-d part.

This is a control for a general comparison theorem, not a
Long lift: the coefficients have x denominators, and the
classical map has a position-preserving form already excluded
for Long. A Long-specific global argument might exclude its
boundary mechanism. Generic index symmetry alone cannot.

## Consequences for the next route

The completed cubic degree is already determined; recomputing
more formal jets cannot determine the missing finite indices.
A Long-specific proof should instead establish a condition
on the actual valuation lattice or another global boundary
condition. Even proving N=3 would not construct a bounded full
lift or contradict an existing theorem. It would remove one
descent ambiguity and permit sharper finite-extension methods.

Do not confuse finite generation of O_Q over O_E with finite
generation of the ambient Weyl algebra A over its image Phi(A).
Lauritzen--Thomsen's earlier theorem says that the latter,
on either side, forces Phi(A)=A. For Long that would contradict
its noninvertible classical reduction, by
[QFL2](BOUNDED_REALIZATION.md). Finite polynomial generator
formulas, a finite division extension, a finite valuation
lattice, and a finite image-ring module are distinct conditions.

The other constructive gate remains QEI1's higher-order
elliptic-energy factorization with all conjugate relations.
The division indices are additional constraints on an actual
candidate, not a replacement for constructing those images.

## Verification boundary

[Exact controls](check_division_indices.py) check the cubic
field charts, the retained full polynomial-control hashes,
the rational Weyl relation and selected companion-matrix
derivation identities. The all-degree valuation-basis theorem,
PBW bounds, lattice equivalence, stabilization and generic
degree arguments are written proofs. The checker does not
enumerate candidate lifts or certify nonexistence.

No independent theorem replay, formal verification, external
review, novelty priority, bounded Long lift, all-representative
nonlifting theorem, or DC2 solution is claimed.
