# Valued matrix lattices and global admission of the shifted model

Let \(D\) be a division ring with a discrete valuation
\(\nu:D^\times\to\mathbb Z\), a **central** uniformizer \(h\) with
\(\nu(h)=1\), valuation ring \(O=\{\nu\geq0\}\), and commutative residue
field \(F=O/hO\). No completeness assumption is needed.

**Valued matrix lemma.** Suppose \(S\in M_n(O)\), \(B\in GL_n(D)\), and
\(S'=BSB^{-1}\in M_n(O)\). Then
\[
 \det(tI-\overline S')=\det(tI-\overline S)
                  \quad\hbox{in }F[t],
 \qquad
 \operatorname{tr}(\overline S')
                   =\operatorname{tr}(\overline S).          \tag{1}
\]
The gauge and its inverse may have arbitrary finite \(h\)-poles and
singular reductions. All determinants and trace comparisons in (1)
are over the commutative residue field. Ordinary quantum matrix trace
is not assumed to be cyclic.

Applied to the full rank-five
[shifted companion model](SHIFTED_EULER_MATRIX_MODEL.md), this excludes
global polynomial admission by **every** rational matrix gauge whose
output is \(h\)-regular. The obstruction also applies to every finitely generated embedded
full global lattice carrying that same generic rank-five representation. It is an admission obstruction for this fixed
representation, not a theorem about all nonscalar carriers.

## Finite elimination over the valuation ring

Every \(B\in GL_n(D)\) has a factorization
\[
 B=U\Delta V,\qquad U,V\in GL_n(O),\qquad
 \Delta=\operatorname{diag}(h^{m_1},\ldots,h^{m_n}),
                  \quad m_1\leq\cdots\leq m_n\in\mathbb Z.   \tag{2}
\]

Here is a direct proof. Move an entry of least valuation \(m\) to the
first pivot by row and column permutations. Write the pivot as
\(h^m u\), with \(u\in O^\times\), and multiply its row by \(u^{-1}\).
Centrality of \(h\) makes the new pivot exactly \(h^m\); all entries
still have valuation at least \(m\).

Clear an entry \(c\) below the pivot by subtracting \(c h^{-m}\)
times the pivot row. Clear a remaining entry \(b\) in the pivot row
by subtracting the pivot column times \(h^{-m}b\). These coefficients
belong to \(O\), and the elementary matrices and their inverses lie
in \(GL_n(O)\). The remaining entries are of the form
\(d-c h^{-m}b\), so their valuations are still at least \(m\).
The remaining square block is invertible because the entire matrix
is invertible. Repeat. The pivot valuations are nondecreasing,
and the finite sequence of operations gives (2). Neither an infinite
series nor a completed division ring was used.

In particular, an integral matrix whose reduction is invertible has
integral inverse: unit pivots suffice in this same elimination.

## Opposite filtrations preserve the classical characteristic polynomial

Put
\[
 C=VSV^{-1},\qquad
 C'=\Delta C\Delta^{-1}=U^{-1}S'U.
\]
Both matrices belong to \(M_n(O)\). Since \(h\) is central,
\[
                  C'_{ij}=h^{m_i-m_j}C_{ij}.                \tag{3}
\]
For the trace alone, \(C'_{ii}=C_{ii}\) is already enough.
Reduction of \(U,V\) gives ordinary similarities over \(F\), so
\[
 \operatorname{tr}(\overline S')
  =\operatorname{tr}(\overline C')
  =\operatorname{tr}(\overline C)
  =\operatorname{tr}(\overline S).
\]
The middle equality comes from equality of diagonal entries, not
cyclicity over \(D\).

Group the indices by equal values of \(m_i\). If \(m_i<m_j\),
integrality of \(C'_{ij}\) forces
\(\nu(C_{ij})\geq m_j-m_i>0\). If \(m_i>m_j\), integrality of
\(C_{ij}\) gives \(\nu(C'_{ij})>0\). Thus \(\overline C\) is block
lower triangular and \(\overline C'\) is block upper triangular.
Their diagonal blocks for each equal-weight group are identical.
The product of the characteristic polynomials of those diagonal
blocks proves (1). More explicitly, the two reductions have
filtrations with the same successive block representations; for a
single operator their semisimplifications are isomorphic.

This does not assert that the reductions themselves are similar.
For example
\[
 S=\begin{pmatrix}0&h\\0&0\end{pmatrix},\qquad
 B=\operatorname{diag}(h^{-1},1),\qquad
 BSB^{-1}=\begin{pmatrix}0&1\\0&0\end{pmatrix}.
\]
Both matrices are integral. Their reductions have different ranks,
but both characteristic polynomials are \(t^2\). Singular-reduction
gauges can change classical extension data while preserving (1).

The same block argument applies simultaneously to any finite
collection of matrices when every matrix and every conjugate is
integral. Their reductions have common filtrations with the same
successive block actions. No such simultaneous assertion is needed
for the admission application: its second image alone suffices.

## Finite projective valuation-ring lattices

A full right \(O\)-lattice here means a finitely generated right
\(O\)-submodule \(\Lambda\subset D^n\) whose \(D\)-span is \(D^n\).
Suppose it is projective and stable under \(S\).

The ring \(O\) is local with Jacobson radical \(hO\): elements of
valuation zero are exactly its units, and \(O/hO=F\).
For clarity the required finite Nakayama step is elementary.
If a finitely generated right module \(M\) satisfies \(M=hM\),
express a chosen generating row using itself with a matrix having
entries in \(hO\). The resulting matrix \(I-H\), \(H\in M_s(hO)\),
has invertible reduction, hence integral inverse by the elimination
above. Multiplying by that inverse shows every generator is zero.

Lift an \(F\)-basis of \(\Lambda/h\Lambda\) to \(\Lambda\).
The finite Nakayama step shows that the resulting map \(O^N\to\Lambda\)
is surjective. Projectivity splits it. Its kernel is finitely
generated, has zero reduction because the chosen residue basis was
a basis, and is therefore zero by the same step. Hence \(\Lambda\)
is free. Its generic rank is \(n\), so \(N=n\).

Its basis columns form a matrix \(C\in GL_n(D)\), and stability is
exactly
\[
                         C^{-1}SC\in M_n(O).
\]
The valued matrix lemma therefore shows that the trace and
characteristic polynomial on \(\Lambda/h\Lambda\) are independent
of the stable finite projective lattice. This includes lattices
whose bases have unequal \(h\)-valuations and whose classical
representations are not isomorphic.

## The full shifted rank-five model cannot be globally admitted

Use the original rank-two Weyl division ring over \(K_0((h))\).
Its \(h\)-Gauss valuation has central uniformizer \(h\), and residue
the original commutative source field
\[
                 F_0=K_0(X,Q,p_x,p_q).
\]
Indeed the minimum valuation of the coefficients of a finite PBW
polynomial is multiplicative: normalized factors reduce to nonzero
commutative polynomials. It extends to the division ring, and
fractions of PBW polynomials supply every rational residue.
This is the valuation used in
[QDI1](QUANTUM_DIVISION_AND_COMPLETION.md).

Assume the shifted model's parameters are \(h\)-regular as in its
classical statement. Write its companion as \(\Pi\), to distinguish
it from an original source momentum. Its diagonal entries vanish,
and the full second image is
\[
 r=X(2-3XQ),\qquad
                   S=r\Pi+\frac4{3r}I_5.                    \tag{4}
\]
It is \(h\)-Gauss regular, since \(r\) is a nonzero residue function.
Consequently
\[
                 \operatorname{tr}(\overline S)=\frac{20}{3r}.
                                                               \tag{5}
\]
No property of the companion's off-diagonal entries is needed here.

For **any** \(B\in GL_5(D)\) such that \(BSB^{-1}\) is \(h\)-Gauss
regular, (1) preserves (5), even if \(B\) and its inverse have poles
or singular classical reductions. A finite globally admitted PBW
matrix over \(K_0[[h]]\) has polynomial classical entries in the
original source coordinates, hence globally regular classical trace.
That contradicts (5).

The denominator is genuine on the original source. The two components
of \(r=0\) are disjoint: the ideals \((X)\) and \((2-3XQ)\) are
comaximal. On \(X=0\), \(\partial_Xr=2\). On \(2-3XQ=0\), \(X\ne0\)
and \(\partial_Qr=-3X^2\ne0\). Thus both components are reduced and
smooth, and \(20/(3r)\) has a simple pole on each in characteristic
zero. Taking the product with the original momentum coordinates
does not remove either divisor.

There is also a basis-independent version. Let \(\mathcal A\) be
the original integral finite-PBW Weyl algebra, and suppose a finite
projective right \(\mathcal A\)-module has generic fiber isomorphic
to this rank-five \(D\)-representation and carries its full second
image as a regular endomorphism. Base change to \(O\) gives a stable
finite projective \(O\)-lattice, so its residue trace is (5).
Reduction modulo \(h\) gives a finite projective module over the
original polynomial source ring. The trace of its endomorphism is
globally regular: an idempotent presentation computes it in that
ring. Its generic trace cannot be (5). Hence no such finite
projective global lattice exists for this fixed generic representation.

## Finite coherent lattices over the original PBW algebra

The projectivity assumption can be removed for embedded lattices over
the original integral algebra. Write
\[
 A_h=K_0[[h]]\langle X,Q,p_x,p_q\rangle_{\rm Weyl},\qquad
 A_0=A_h/hA_h=K_0[X,Q,p_x,p_q],\qquad D=\operatorname{Frac}(A_h).
\]
Here \(A_h\) means finite PBW sums with regular scalar-series
coefficients, not the spatially unbounded \(h\)-adic completion.
The Gauss valuation has the useful saturation property
\[
          a\in A_h,\ \nu(a)\geq m
                   \quad\Longrightarrow\quad h^{-m}a\in A_h
                         \qquad(m\geq0).                  \tag{6}
\]
Define the multiplicative set
\[
                     \Sigma_0=\{a\in A_h:\nu(a)=0\}.
\]

First, \(A_h\) is an Ore domain. Over \(K_0((h))\), the PBW filtration
gives common right multiples by finite-dimensional linear algebra:
for fixed nonzero \(a,s\) of bounded degree and large \(N\), the map
\((u,v)\mapsto au-sv\) from two copies of the degree-\(N\) filtered
space into the degree-\(N+d\) space has nonzero kernel, because its
source dimension is larger. The domain property ensures both
components of a nonzero kernel vector are nonzero. Left multiples
are identical. Clearing a common central \(h\)-denominator gives
the Ore equations in \(A_h\).

The set \(\Sigma_0\) is itself a right and left Ore denominator set.
For example, start with \(a u=s v\), where \(a\in A_h\), \(s\in\Sigma_0\)
and \(u\ne0\). If \(m=\nu(u)\), multiplicativity gives
\(\nu(v)=\nu(a)+m\geq m\). Property (6) allows cancellation of the
central factor \(h^m\):
\[
                 a u_0=s v_0,\qquad u_0\in\Sigma_0,\quad v_0\in A_h.
\]
The case \(a=0\) is immediate, and the left Ore argument is the same.
There are no annihilator issues because \(A_h\) is a domain.

Consequently
\[
             O=A_h\Sigma_0^{-1},\qquad O/hO=\operatorname{Frac}(A_0)=F_0.
                                                               \tag{7}
\]
Indeed an integral right fraction \(a b^{-1}\), \(a,b\in A_h\), has
\(\nu(a)\geq\nu(b)=m\geq0\). Cancel \(h^m\) in both numerator and
denominator using (6), leaving a denominator in \(\Sigma_0\).
The reverse inclusion follows from valuations. Reduction gives every
fraction in \(F_0\). This Ore localization is flat on the left:
localizing right modules uses common right denominators, and
\(m/s=0\) precisely when \(mt=0\) for some \(t\in\Sigma_0\).
These descriptions show directly that kernels, as well as
surjections, are preserved by localization.

Now let \(\Lambda\subset D^5\) be a finitely generated right
\(A_h\)-submodule with \(\Lambda D=D^5\), stable under the fixed
second image \(S\). It is automatically \(h\)-torsion-free, being
embedded in \(D^5\); no projectivity, saturation or torsion-free
classical reduction is assumed. Flat localization identifies
\(\Lambda_O=\Lambda\otimes_{A_h}O\) with its actual \(O\)-span in \(D^5\).

Every finitely generated \(O\)-submodule of \(D^n\) is free. To see
this without a projectivity assumption, put its finitely many
generators into a rectangular matrix. Apply the same least-valuation
row and column elimination as in (2), stopping when the remaining
block is zero. Integral invertible column operations preserve its
column span, while integral invertible row operations identify that
span with the direct sum of the nonzero pivot modules \(h^{m_i}O\).
It is therefore free, with rank equal to its \(D\)-rank.
Thus \(\Lambda_O\) is a full free \(O\)-lattice, and (1) applies.

Right-exact base change, together with (7), gives the precise
classical identification
\[
 (\Lambda/h\Lambda)\otimes_{A_0}F_0
                  \simeq \Lambda_O/h\Lambda_O.             \tag{8}
\]
In particular this generic classical space has dimension five,
and the characteristic polynomial of its induced \(S\)-action is
the lattice-invariant polynomial of \(\overline S\).

It remains to justify regular coefficients when \(\Lambda/h\Lambda\)
is merely coherent. For any finite module \(M\) over the commutative
normal domain \(A_0\) and any endomorphism \(\phi\), the determinant
trick on a finite generating set produces a monic polynomial
\(f(t)\in A_0[t]\) with \(f(\phi)=0\): lift the action to the generating
free module and use its adjugate characteristic identity. Hence all
eigenvalues of the generic \(F_0\)-linear action are integral over
\(A_0\). Their elementary symmetric functions are integral as well.
Those functions are the generic characteristic coefficients and
belong to \(F_0\); normality of \(A_0\) puts them in \(A_0\).
Possible \(A_0\)-torsion in \(M\) does not affect this argument.

Apply this to \(M=\Lambda/h\Lambda\). Its generic trace must belong
to \(A_0\), whereas (5) and (8) force it to be \(20/(3r)\), which
does not. Thus the fixed shifted generic representation has **no
finitely generated embedded full \(A_h\)-lattice stable even under
its second image \(S\)**. In particular there is no coherent full
global carrier for that representation; projectivity is unnecessary.
The conclusion uses finite generation and the specified original
PBW algebra, generic module and \(h\)-Gauss valuation throughout.

## Scope

This strengthens the regular-gauge trace obstruction in the shifted
model note. Within its fixed generic rank-five representation,
arbitrary rational \(h\)-polar gauges and finitely generated embedded
full lattices over the original integral PBW algebra cannot avoid the
original-source trace pole.
No cyclic quantum trace, scalar extraction, or scalar Euler
obstruction is used.

The general lemma also applies in a completed division ring if it
has the same central discrete uniformizer and commutative residue
field. It concerns full-rank lattice changes. Selecting a proper
formal sector changes the generic representation and is not a
full-rank lattice change.

Other nonscalar representations, changed matrix images not related
to this fixed generic representation, and other constructions remain
outside this result. No bounded scalar Long lift or DC1/DC2 resolution
is asserted. The proof is the finite elimination and filtration
argument above; no external theorem or finite computational sample
is a premise. Formal verification, independent whole-proof replay
and external review remain unclaimed.
