# Finite algebraic graph carriers and matrix realization

14 November 2026. **QMC1.** A compatible free left frame is necessary
for the scalar realization pursued in [QFL2](BOUNDED_REALIZATION.md),
but it is not a prerequisite for a useful finite matrix realization.

For every symplectic polynomial map of classical field degree delta>1,
a finite-left, parameter-flat algebraic bimodule deforming its graph
already produces a non-surjective endomorphism of Mat_2(A_n).
The generic left module is automatically projective of rank one;
two copies are free by Stafford's theorem. The generic rational
right action has finite left and right division indices, each at
least delta, while its completed index is exactly delta.

This is a sufficient criterion. No such algebraic carrier for Long
is constructed here. The conclusion concerns a matrix algebra and
does not itself supply an endomorphism of the scalar A_2.

## Hypotheses: finite on the left, not only coherent as a bimodule

Let K be a characteristic-zero field, R=K[[h]], F=K((h)), and
\[
 A_R=R\langle x_1,\ldots,x_n,p_1,\ldots,p_n\rangle/
       ([p_i,x_j]-h\delta_{ij},[x_i,x_j],[p_i,p_j]).
                                                               \tag{1}
\]
Only finite sums of spatial PBW monomials belong to A_R.
Its h-adic completion Ahat allows their degrees to grow.
The generic algebra A=A_R tensor_R F is a Weyl algebra over F.

Let psi be a polynomial symplectic endomorphism of
B=K[x_1,...,p_n], and let
\[
                       [K(x,p):K(\psi)]=\delta.                \tag{2}
\]
Assume an R-central A_R-bimodule M satisfies:

1. M is finitely generated as a left A_R-module.
2. M is R-flat.
3. There is an identified classical bimodule
   M/hM = B_psi: its left module is B, and the right action
   of a classical generator z_i is multiplication by psi_i.

All right actions must be defined on M itself, and every Weyl
relation must hold there. A coherent A_R-bimodule with a formally
free completion need not satisfy condition 1; the
[QGM1 example](FORMAL_GRAPH_BOUNDARY_MODEL.md) already separates them.
Finite generation of the ambient algebra over its right image is
a different, much stronger condition; it is not assumed here.

## A formal frame exists, without being an algebraic frame

Choose m in M lifting the identified classical left generator 1.
The maps
\[
 A_R/h^N A_R\longrightarrow M/h^N M,\qquad a\longmapsto am
                                                               \tag{3}
\]
are compatible isomorphisms for every N. They are isomorphisms at
N=1. R-flatness identifies the next kernel h^(N-1)M/h^N M with
M/hM, and induction on the two short exact sequences proves the
assertion. Hence Mhat is freely generated on the left by m.

In this formal frame its right action is a full formal quantum
lift Phi of psi. No common spatial bound for Phi is inferred.
Neither m generating M nor a finite inverse for (3) is asserted.
This is exactly the distinction exposed by
[QPF1](PROJECTIVE_FRAME_BOUNDARY.md).

## Generic projectivity and rank one

We use the finite-left bimodule projectivity theorem:
over a characteristic-zero Weyl algebra, every bimodule finite
on the left is projective and torsion-free on the left.
See [Kanel-Belov--Fedorov, Lemma 3.1](https://arxiv.org/html/2610.33919v1).
That lemma imposes no prescribed arithmetic-support hypothesis.

Its mechanism is worth specifying. A finite-left bimodule has
the same growth dimension whether viewed on the left or as a
bimodule: right generator actions increase a good left filtration
by bounded amounts. Its left torsion submodule would have growth
at most 2n-1, contradicting Bernstein's inequality for A_(2n).
Positive Ext groups against A are finite right modules carrying
the opposite action, and become zero after division localization.
The same torsion argument makes them zero; shortening a finite
projective resolution gives projectivity. These are the cited
lemma's arguments, not an independent review of the whole preprint.

Apply it to P=M tensor_R F. This module is nonzero by R-flatness.
To determine its rank, choose a finite left presentation
A_R^a -> A_R^b -> M -> 0. Its reduction has cokernel B, so the
relation matrix has rank b-1 over K(x,p). A corresponding square
block is invertible in the completed quantum division ring: lift
its inverse over K(x,p), then invert the identity-mod-h error by
a convergent geometric series. Thus the relation matrix has
generic division rank at least b-1 and P has rank at most one.
Torsion-freeness and P nonzero give
\[
                       \operatorname{rank}_A P=1.             \tag{4}
\]

In particular P is a projective generator. Its trace ideal in
the simple ring A is nonzero and hence all of A.
The algebra E=End_A(P)^op is a domain Morita equivalent to A,
and the full right action is a unital homomorphism rho:A -> E.

## The rational action retains the classical degree

Realize P as a fractional left ideal and clear a right denominator
to embed it in A. Since M has finitely many left generators,
multiplication by one power of h clears the finitely many scalar
poles in their images. We obtain a nonzero left-linear map
ell:M -> A_R, injective after passage to F. It realizes the
right action as multiplication by elements rho(a) of D=Frac(A).

Completing ell and evaluating it at the formal generator m gives
a nonzero q in Ahat. Nonzero follows because ell has nonzero
image in the separated ring A_R. In the completed division ring,
\[
                         q\rho(a)=\Phi(a)q,\qquad
                         \rho(a)=q^{-1}\Phi(a)q.               \tag{5}
\]
Conjugation by q preserves h valuation, and acts trivially on
the residue field, which is commutative. Thus rho has residue map
psi. Dominance of psi prevents cancellation of the leading
valuation of any nonzero source polynomial. Extending to fractions,
rho:D -> D is isometric.

Consequently the completion theorem of
[QDI1](QUANTUM_DIVISION_AND_COMPLETION.md) gives
\[
                 [\widehat D:\widehat{\rho(D)}]_\ell
                =[\widehat D:\widehat{\rho(D)}]_r=\delta.       \tag{6}
\]
For the non-surjectivity below, only the weaker inequality is
needed. Lift a residue basis of K(x,p) over K(psi) to rational
operators of valuation zero. A putative finite relation over
rho(D), scaled to minimum coefficient valuation zero, would
reduce to a nonzero relation among that residue basis.
This proves both generic indices are at least delta.

The comparison uses an algebraic dual section, not an assumed
bounded formal gauge. It permits every deformation of the
first image, including momentum terms.

## Both generic indices are finite: a direct growth bound

Choose a good Bernstein filtration on P and C>=1 such that
every right generator increases its index by at most C.
This is possible from finitely many left generators and their
right-action expressions. Its cumulative Hilbert function has
leading term m^(2n)/(2n)!: a rank-one fractional ideal differs
from A only by lower-growth terms, and good filtrations differ
by bounded shifts.

Suppose d elements of D are right-linearly independent over rho(D).
A common left multiple places them in the chosen left ideal P
without changing right independence. Multiplying them on the
right by source PBW words of degree at most m gives
d*binomial(m+2n,2n) independent F-vectors. All lie in filtration
index Cm+O(1). Comparison of leading coefficients gives
\[
                      [D:\rho(D)]_r\le C^{2n}.                \tag{7}
\]
To place a finite set in P, use the left Ore property with any
nonzero element of P, then a common left multiple of the
resulting multipliers.

For the other side take the algebraic left dual P^*=Hom_A(P,A).
It is a finite projective right A-module of rank one.
The source acts on its left by precomposition with right
multiplication on P. In a fractional realization it is the
action q -> rho(a)q. A good right filtration gives an analogous
constant C^*. The right Ore property places any finite left
independent set in P^*, and the same dimension comparison yields
\[
       \delta\le [D:\rho(D)]_\ell\le(C^*)^{2n},\qquad
       \delta\le [D:\rho(D)]_r\le C^{2n}.                     \tag{8}
\]
No index-symmetry theorem is needed for these individual bounds.
As in QDI1, equality of a generic index with delta is a separate
valuation-lattice condition; finite generation of M over A_R
does not imply that different finiteness property.

## Two copies give a finite matrix endomorphism

Stafford's projective-freeness theorem says that a finite
projective Weyl module of rank at least two is free. Therefore
P direct_sum P is free of rank two, and
\[
             \operatorname{Mat}_2(E)\simeq
             \operatorname{End}_A(P\oplus P)^{op}
                    \simeq\operatorname{Mat}_2(A).            \tag{9}
\]
We use Stafford, *Module Structure of Weyl Algebras*, Theorem3.6,
[original paper](https://doi.org/10.1112/jlms/s2-18.3.429).
The original PDF could not be fetched in this session.
The matrix argument, with this exact external dependency, is also
proved in [Kanel-Belov--Fedorov, Lemma5.2 and Theorem5.3](https://arxiv.org/html/2610.33919v1).
No independent proof of Stafford's theorem is claimed.

Entrywise rho followed by (9) is a unital injective endomorphism
\[
                 u:\operatorname{Mat}_2(A)
                             \longrightarrow\operatorname{Mat}_2(A).
                                                               \tag{10}
\]
If u were surjective, the entrywise map would be surjective and
rho(A)=E. Localizing this equality would give rho(D)=D.
Indeed, realize P as a nonzero left ideal I of A and choose
nonzero d in I. Then A*d is contained in E, since I*(a*d)
is contained in A*d and hence in I. Both a*d and d belong
to E, so every a=(a*d)*d^(-1) belongs to its division ring.
Equation (8) excludes this when delta>1.

Every entry in the finitely many images of algebra generators
is a finite PBW polynomial over F. Thus the generic matrix
endomorphism is an actual finite construction. Stafford's
generic basis need not extend over R, and this does not give
a bounded scalar formal lift with reduction psi.

The finitely many scalar coefficients of u generate a finitely
generated field L over Q. Choose an embedding L -> C. Surjectivity
cannot appear under a field extension: a preimage of a generator
would involve finitely many words, and its coefficients solve
a finite linear system over L. Solvability after extension implies
solvability over L. Thus a carrier with delta>1 also yields a
non-surjective endomorphism of Mat_2(A_n(C)).

For Long, delta=3 by QDI1. The carrier would therefore disprove
the matrix/generalized rank-two Dixmier statement. Obtaining
the scalar A_2 endomorphism still requires a generic free
rank-one frame or another scalar construction. The matrix
conclusion is not presented as completion of that scalar goal.

## A finite corner certificate and its curvature test

A directly checkable sufficient carrier certificate consists
of an idempotent e in Mat_s(A_R), and matrices r_i satisfying
\[
 \begin{split}
   &e^2=e,\quad r_i=er_i e,\quad
       [r_i,r_j]=hJ_{ij}e,\\
   &e\bmod h=e_{11},\qquad
       r_i\bmod h=\psi_i e_{11}.
 \end{split}                                                   \tag{11}
\]
The row module M=A_R^(1 times s)e is finite projective and R-flat.
Its right action is row multiplication by r_i. Its classical
bimodule is precisely B_psi, so the theorem applies. There is
no requirement that e be polynomially conjugate to e_11 over R.
No such certificate for Long is supplied here.

The matrix equations do not follow by compressing ordinary
generators. For any a,b in a larger matrix algebra, put f=1-e.
Then
\[
 [eae,ebe]-e[a,b]e
                   =ebfae-eafbe.                             \tag{12}
\]
These off-corner terms are the exact compression defect.
The full Weyl relations require every relevant defect to vanish
or be corrected in the algebraic corner. Projectivity alone
does not provide those corrections.

The controls verify (11) for an explicit nonconstant idempotent
polynomially conjugate to e_11, with a complete rank-two identity
graph right action. They also verify (12) and a nonzero defect
for the actual QPF1 nonfree projective idempotent, after its
allowed specialization m=2. This preserves the negative control
rather than inventing a right action on it.

## Consequence for the research target

The scalar route still needs finite algebraic descent and a free
frame. The matrix route needs finite algebraic descent with every
right relation, but generic projectivity and the two-copy frame
then follow. This makes an actual finite-left graph carrier a
substantive intermediate milestone even if its rank-one frame
remains unresolved.

Neither the fixed-energy elliptic module in
[QIM1](ELLIPTIC_ISOGENY_MATRIX_MODULE.md) nor QPF1's left projective
model has the required full Long right action. Prescribed
Frobenius center maps are stronger than the h=0 reduction used
here; the arithmetic nonlifting theorem is not silently applied.

[Exact controls](check_matrix_carrier.py) check the finite
matrix certificate and compression identities. Generic projectivity
uses the cited finite-left bimodule lemma; matrix freeness uses
Stafford's theorem. The carrier criterion, degree transport and
growth estimates are written proofs. No independent theorem replay,
formal verification, external review, novelty priority, actual
Long carrier, bounded scalar Long lift or Dixmier solution is claimed.
