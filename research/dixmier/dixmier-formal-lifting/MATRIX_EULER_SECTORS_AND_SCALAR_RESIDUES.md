# Matrix Euler sectors and a single classical Long branch

Two elementary restrictions distinguish the shifted companion model from a
matrix lift with the same classical Long branch in every sector. First, an
actual matrix fifth root of a scalar operator, with nonzero scalar classical
residue, is itself scalar. Second, retaining the QQL1 Euler operator and
adding any constant matrix of sector shifts forces five cyclic sectors and
zero trace. Neither restriction excludes genuinely nonconstant matrix Euler
corrections or different quantum source operators.

These are written algebraic proofs. The exact scalar spectrum and scalar
root exclusion used below are those of
[QQL1](QUINTIC_ORDERING_LAURENT_OBSTRUCTION.md); the actual source model and
Long residue are those of
[QCS1](CLASSICAL_QUINTIC_COVER_SQUARE.md). The explicit companion completion
and its admission scope are recorded in
[the shifted matrix model](SHIFTED_EULER_MATRIX_MODEL.md).

## A matrix root with scalar residue is scalar

Let \(D\) be a division ring with a separated discrete valuation, a central
uniformizer \(h\), valuation ring \(\mathcal O_D\), and commutative residue
field \(F\) of characteristic zero. Assume the rational constants are units
of \(\mathcal O_D\), as in the rational Gauss valuation of the source Weyl
division ring. Matrices are regular when all entries belong to
\(\mathcal O_D\). Let \(A\in\mathcal O_D\) and \(p_0\in F^\times\).

**Scalar-residue root lemma.** If
\[
 P\in M_n(\mathcal O_D),\qquad
 P^5=AI_n,\qquad \overline P=p_0I_n,
                                                               \tag{1}
\]
then \(P=pI_n\) for an actual \(p\in\mathcal O_D\), and
\(p^5=A\), \(\overline p=p_0\).

First, two regular roots \(P,T\) of \(AI_n\) with this same scalar
residue must coincide. If not, let \(m\geq1\) be the least valuation of
an entry of \(P-T\). The ordered telescoping identity is
\[
 P^5-T^5=\sum_{i=0}^{4}P^i(P-T)T^{4-i}.
\]
Divide by the central \(h^m\) and reduce. The nonzero matrix
\(\Delta=\overline{h^{-m}(P-T)}\) would satisfy
\[
                         0=5p_0^4\Delta,
\]
which is impossible. No completeness or infinite-series root construction
is used here.

For every \(C\in\operatorname{GL}_n(\mathbb Q)\), the matrix
\(CPC^{-1}\) is another regular root of \(AI_n\) with residue
\(p_0I_n\). The uniqueness just proved gives \(CPC^{-1}=P\).
For \(n\geq2\), commuting with all \(I_n+E_{ij}\), \(i\ne j\),
forces every off-diagonal entry to vanish and all diagonal entries to
coincide. The case \(n=1\) is immediate. This proves the lemma.

The hypothesis of scalar residue is also forced by a sole nonzero
classical branch. If \(\overline A=p_0^5\ne0\),
\(\overline P^{\,5}=p_0^5I_n\), and every eigenvalue of
\(\overline P\), over an algebraic closure of \(F\), is \(p_0\),
then \(\overline P=p_0I_n\). Indeed \(t^5-p_0^5\) is squarefree,
so its matrix roots have no nontrivial Jordan blocks. A Jordan or
unipotent classical root cannot evade (1).

In the h-regular QCS1 model, put \(A=UG(J)\), with
\(0\ne G\in k(J)\), whenever its residue is the prescribed
\(p_0^5\ne0\). QQL1 excludes an actual scalar fifth root in the
source division ring. The lemma therefore excludes an actual regular
matrix root with this one classical branch in every sector. It does not
apply to the shifted companion, whose classical eigenvalues are the five
different fifth-root branches, or to a construction whose fifth power is
genuinely nonscalar.

The same argument gives a more general scalarization criterion. Fix a
finite left-coefficient polynomial
\[
 \Phi(T)=\sum_{i=0}^{d}a_iT^i,\qquad a_i\in\mathcal O_D,
 \qquad
 \Phi(P):=\sum_{i=0}^{d}(a_iI_n)P^i.
\]
The coefficients \(a_i\) need not commute with \(P\). Suppose
\[
 \overline\Phi(p_0)=0,\qquad
 \overline\Phi'(p_0)\ne0,\qquad
 P\in M_n(\mathcal O_D),\quad
 \Phi(P)=0,\quad \overline P=p_0I_n.                         \tag{1a}
\]
Here \(p_0\in F\), and the reduction and derivative are taken in the
commutative polynomial ring \(F[T]\). Then \(P=pI_n\) for an actual
\(p\in\mathcal O_D\) satisfying \(\sum_i a_ip^i=0\).

Indeed, for two such roots \(P,T\), the first nonzero valuation of their
difference and the termwise ordered telescoping identities give
\[
 \overline{h^{-m}\bigl(\Phi(P)-\Phi(T)\bigr)}
                     =\overline\Phi'(p_0)\Delta.
\]
Every factor in this leading reduction belongs to the commutative
residue field. Noncommutation of the actual left coefficients with the
roots therefore causes no extra term. The nonzero derivative proves
uniqueness. Constant rational matrix conjugation fixes every
\(a_iI_n\) and preserves the scalar residue, so the same commutant
argument forces the root to be scalar.

If a regular root of this relation has only the classical eigenvalue
\(p_0\), and \(p_0\) is a simple root of \(\overline\Phi\), its
classical minimal polynomial divides both a power of \(T-p_0\) and
\(\overline\Phi\); hence it divides \(T-p_0\). Thus its residue is
already \(p_0I_n\), and (1a) applies. This criterion does not cover
matrix-valued polynomial coefficients or a repeated classical root.

## Constant Euler shifts force five cyclic sectors

Let \(D\) now be the actual QQL1 division ring over its characteristic-zero
central coefficient field \(k\), with \(0\ne\eta=ch\in k\), and retain
its scalar operator \(J\). Its scalar adjoint spectrum is
\[
 \{\rho\in k:\ [J,x]=\rho x\text{ for some }0\ne x\in D\}
                              =5\eta\mathbb Z.                 \tag{2}
\]
Let \(B\in M_n(k)\) be arbitrary; diagonalizability is not assumed.
Suppose \(P\in\operatorname{GL}_n(D)\) satisfies
\[
 \widetilde J=JI_n-\eta B,\qquad
                  [\widetilde J,P]=-\eta P.                     \tag{3}
\]
This includes every factorization \(\widetilde J=QP\),
\([P,Q]=\eta I_n\), with this constant-shift form of \(\widetilde J\).

Enlarge the central coefficient field by a finite algebraic extension
splitting the characteristic polynomial of \(B\). The original source
division ring embeds in the Weyl division ring with these enlarged
constants; no tensor product is being asserted to be a division ring.
QQL1 applies with the enlarged coefficient field and proves (2) there
as well. A constant change of basis puts \(B\) in Jordan form.

For a block \(X\) of \(P\) from the generalized \(\mu\)-eigenspace
to the generalized \(\lambda\)-eigenspace, equation (3) gives
\[
 (\operatorname{ad}J-\rho)X
       =\eta(N_\lambda X-XN_\mu),\qquad
                \rho=\eta(\lambda-\mu-1),                     \tag{4}
\]
where the two \(N\)'s are the respective nilpotent Jordan parts. Left
multiplication by \(N_\lambda\) and right multiplication by
\(N_\mu\) commute and are nilpotent. Their difference is therefore
nilpotent and commutes with the entrywise \(\operatorname{ad}J\).
Thus \((\operatorname{ad}J-\rho)^sX=0\) for some finite \(s\).
If \(X\ne0\), choose the last nonzero member of this chain. A nonzero
entry of that member is a scalar eigenoperator of eigenvalue \(\rho\).
By (2), a nonzero block is possible only when
\[
                         \lambda-\mu-1\in5\mathbb Z.           \tag{5}
\]
This proves the support assertion even for nondiagonalizable \(B\).

Group its generalized eigenspaces by classes in the additive quotient
\(k'/(5\mathbb Z)\), where \(k'\) is the splitting field. If
\(W_\alpha\) is the right-division-space sum for class \(\alpha\),
then (5) says
\[
                         P W_\alpha\subseteq W_{\alpha+1}.     \tag{6}
\]
Translation by \(1\) has order exactly five on this quotient in
characteristic zero. Invertibility of \(P\) forces equal dimensions
around every supported five-cycle: its image in a target summand comes
only from the preceding summand, and surjectivity and injectivity give
equality. Consequently
\[
                       5\mid n.                              \tag{7}
\]
For \(r=1,2,3,4\), \(P^r\) carries each summand in (6) to a distinct
summand. All its diagonal blocks are zero, so
\[
                       \operatorname{tr}(P^r)=0
                              \quad(1\leq r\leq4).             \tag{8}
\]
This is an exact identity over \(D\), not a cyclic-trace assertion for
arbitrary matrices over a division ring. Passing back from the Jordan
basis preserves the trace because that change of basis has central
entries: its scalar coefficients can be moved past each matrix entry.

If \(P\) is h-regular, its classical reduction therefore cannot have
only the eigenvalue \(p_0\ne0\), with any Jordan multiplicities:
such a matrix would have trace \(np_0\ne0\), contradicting (8).
Changing the constant sector shifts, or making their classical matrix
unipotent, cannot produce the same nonzero Long branch in every sector
while retaining (3). This conclusion does not require \(P^5\) to be
scalar.

## Consequence for the full second image and admission

Suppose this model is h-regular in the actual rank-two source division
ring and its second image retains the QCS1 affine formula
\[
                         S=rP+\frac4{3r}I_n.
\]
Equation (8) gives the exact trace and its classical reduction:
\[
          \operatorname{tr}S=\frac{4n}{3r},\qquad
          \operatorname{tr}\overline S=\frac{4n}{3r}.             \tag{9}
\]
The [valued matrix lattice lemma](VALUED_MATRIX_LATTICE_AND_SHIFTED_ADMISSION.md)
preserves the trace of the classical reduction under an arbitrary rational
matrix change of basis for which both the original and transformed
operators are h-regular. The change of basis and its inverse may have
parameter poles. Thus any such transformed second image still has
classical trace \(4n/(3r)\).

In Long's original global chart, \(r=X(2-3XQ)\). The coefficient
\(4n/3\) is nonzero, and \(1/r\) has a genuine pole on both reduced
zero components. On \(X=0\), the \(X\)-derivative of \(r\) is \(2\);
on \(2-3XQ=0\), one has \(X\ne0\) and the \(Q\)-derivative is
\(-3X^2\). A globally admitted regular PBW matrix has polynomial
classical entries and polynomial trace. Hence no rational matrix gauge
in this scope makes this \(S\) a globally admitted regular PBW matrix.
The same trace argument applies to every finitely generated embedded
full lattice over the original integral PBW algebra. The coherent-lattice
proof in the valued matrix note works with n in place of five: localization
gives a full free valuation-ring lattice, and the determinant trick over
the normal original classical source ring forces the generic trace to be
regular. It therefore contradicts 4n/(3r) without a projectivity assumption.

The root lemma and the constant-sector lemma close the two stated ways
of making the old homogeneous construction carry one classical branch.
They leave open genuinely nonscalar fifth-power deformations,
nonconstant matrix corrections to the Euler operator, different source
operators, and other nonscalar constructions. No general matrix
nonexistence theorem or DC1/DC2 conclusion is asserted.
