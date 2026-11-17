# Nonconstant formal graph carriers and the finite admission gap

A rank-two bundle on the Lagrangian graph of a polynomial symplectic map
can carry genuinely nonconstant flat projective connection data. The
Baranovsky--Chen quantization theorem supplies an h-complete formal module
with such data. The polynomial connection below has no invariant rational
line, by an elementary infinity argument. This gives a concrete formal
carrier outside the constant-sector Euler construction; it supplies no
finite PBW realization or scalar Long lift.

## The inspected primary theorem

The source is Vladimir Baranovsky and Taiji Chen, *Quantization of Vector
Bundles on Lagrangian Subvarieties*, DOI
[10.1093/imrn/rnx230](https://doi.org/10.1093/imrn/rnx230).
The inspected version is the 22-page
[author-hosted 2017 journal advance PDF](https://bpb-us-e2.wpmucdn.com/sites.uci.edu/dist/6/3204/files/2017/11/Paper.pdf),
accepted 6 September 2017 and carrying a 3 November 2017 download stamp.
Its [cached exact bytes](sources/baranovsky-chen-authoraccepted2017.pdf)
have SHA256
`5e1deef94fdda1c98ead8c27668fe65094be161f28669449ba6221d4a4e09efb`.
This is not being identified with arXiv:1701.01542v1, submitted 6 January
2017. Page numbers below refer to this cached advance PDF.

Here is the theorem in the notation needed for the application. Let
\(\mathcal X\) be smooth algebraic symplectic over a characteristic-zero
field, \(\mathcal O_h\) a fixed formal quantization, and
\(j:Y\hookrightarrow\mathcal X\) a smooth closed Lagrangian. Write its
period as \([\omega]+h\omega_1+h^2\omega_2+\cdots\).
Theorem 1.1, printed page 4, characterizes quantizability of a rank-\(r\)
bundle \(E\) by
\[
 j^*\omega_i=0\quad(i\geq2),\qquad
 \mathbb P(E)\text{ has a flat algebraic connection},\qquad
 \frac1r c_1(E)=\frac12c_1(K_Y)+j^*\omega_1
                                  \quad\text{in }H^2_F(Y).       \tag{1}
\]
The last restriction uses the canonical refined lift. The module in this
statement is flat over the parameter ring and h-adically complete
(definition, printed pages 1--2). Theorem 5.1, printed pages 15--16,
keeps track of a chosen flat formal projective bundle in the lifting
construction. The first remark of Section 6, printed page 19, permits
\(E_h\simeq E[[h]]\) as parameter-module sheaves when \(Y\) is affine.
These are formal assertions, not finite-operator statements.

The exact locations and application hypotheses are retained in the
[primary-source audit](formal_graph_carriers_source_audit.json).

## Application to a polynomial symplectic graph

Let \(V=\mathbb A_k^{2n}\) with its standard symplectic form and let
\(\phi:V\to V\) be a polynomial symplectic endomorphism. Its graph
\[
 Y=\Gamma_\phi=\{(\phi(x),x):x\in V\}
       \subset \mathcal X=V_{\rm target}\times V_{\rm source}^{\rm op}
                                                               \tag{2}
\]
is closed, smooth and Lagrangian. The source projection identifies
\(Y\) with \(V\); no inverse to \(\phi\) is required. The graph pulls
the product symplectic form back to \(\phi^*\omega-\omega=0\), and
its dimension is half that of \(\mathcal X\).

Choose the standard product Weyl/Moyal formal quantization. Both affine
spaces have zero degree-two algebraic de Rham cohomology. On \(Y\), the
refined group \(H^2_F(Y)\) vanishes as well: affineness computes the
truncated de Rham complex by global polynomial forms, and its degree-two
cohomology is zero by the polynomial Poincare homotopy. Also \(K_Y\) is
trivial. Hence (1) holds for \(E=\mathcal O_Y^r\) with any specified
compatible flat projective connection. The Picard algebroid occurring in
the refinement has zero class; choosing a splitting identifies its
module condition with an ordinary flat connection on \(E\).

The ambient product quantization acts as a formal target-source
bimodule. Modulo \(h\), a target function \(f\) acts on \(E\) by
\(\phi^*f\) times the identity. Thus all rank-\(r\) sectors have the
same classical graph symbols; fifth-root deck branches are not part of
this construction. For an actual graph with trivial classical bundle,
affineness and h-flatness give a basis over the completed source Weyl
algebra; [the source-freeness proof and finite-presentation criterion](FINITE_GRAPH_PRESENTATIONS_AND_POLAR_BRANCHES.md)
spell out that argument. It still supplies no stable finite free lattice
over the uncompleted source Weyl algebra.

Rank-one formal graph existence is already explained in
[QFL1's graph comparison](PROOF.md), using Baranovsky--Ginzburg--Kaledin--Pecharich,
Theorem 1.1.4. The present application records higher-rank connection
data and the missing finite step; it is not a new claim of formal
existence for Long's scalar map. The canonical full polynomial Long
frame is an actual input, constructed by
[long_images()](../dixmier-audit/verify_long_pbw.py) and independently
checked by [the graph-curvature verifier](check_long_graph_curvature.py).
The pair-only residue-transport inputs are different data and do not,
by themselves, establish such a full frame.

## A rank-two connection with no invariant rational line

Use an affine coordinate \(X\) on \(Y\), and put
\[
 E=\mathcal O_Y^2,\qquad
 \nabla=d-M(X)\,dX,\qquad
 M(X)=\begin{pmatrix}0&1\\X&0\end{pmatrix}.                    \tag{3}
\]
This is a polynomial algebraic connection. Its curvature vanishes:
both \(d(M\,dX)\) and \((M\,dX)\wedge(M\,dX)\) are zero.
Moreover \(\operatorname{tr}M=0\), so the determinant connection is
trivial and (3) induces a flat projective connection on
\(\mathbb P(E)\).

Let \(K\) be the rational function field in the remaining affine
coordinates, so the generic connection is over \(K(X)\) with
\(\partial_X K=0\). It has no invariant one-dimensional \(K(X)\)-space.
Indeed, such a line has a nonzero generator \((f,g)^t\) with
\[
                 (\partial_X-M)(f,g)^t=a(f,g)^t
                                                               \tag{4}
\]
for some \(a\in K(X)\). If \(f=0\), the first coordinate of (4) is
\(-g=0\), impossible. For \(f\ne0\), set \(w=g/f\). The two
coordinates of (4) give
\[
                             w'+w^2=X.                         \tag{5}
\]
There is no rational solution. For a nonzero rational \(w\), put
\(m=\deg_\infty w\in\mathbb Z\). If \(m\geq1\), the square has
degree \(2m\), the derivative has degree \(m-1\), and the leading
square cannot cancel. Their sum has degree \(2m\geq2\), unlike \(X\).
If \(m\leq0\), the square has degree at most zero and the derivative
has degree at most \(-1\), so their sum again cannot be \(X\).
The zero function is not a solution either. This proof remains valid
after any extension of the constant field.

In particular, the connection cannot be rationally gauged to a flat
constant-matrix connection after algebraically closing the constants:
flat constant coefficient matrices commute, have a common eigenline,
and a rational gauge would transport that line to a forbidden line of
(3). This is a statement about the connection, not scalar descent of
its quantized Weyl actions.

For the specific projective datum in (3), take its constant-in-\(h\)
extension as a flat \(\operatorname{PGL}_2(k[[h]])\)-bundle.
With the Picard-algebroid splitting above, (3) is a compatible rank-two
module. The successive lifting statement in Theorem 5.1 therefore
permits a formal quantization with this projective datum. No invariant
rational classical line is being selected in that application.

No exact scalar fifth-power equation or Euler operator of the form
\(JI-hcB\) with constant \(B\) is a premise of this construction. Thus the
[scalar-residue and constant-sector restrictions](MATRIX_EULER_SECTORS_AND_SCALAR_RESIDUES.md)
do not by themselves exclude it. Their hypotheses have changed; no
unrestricted matrix construction theorem follows.

## The elementary Airy embedding does not improve a finite bound

The connection itself has a simple actual polynomial matrix Weyl
realization, up to the constant gauge \(\operatorname{diag}(1,-1)\),
which exchanges the signs of \(M\). If \([p,X]=h\), define
\[
 \iota(X)=XI_2,\qquad
 \iota(p)=pI_2+h\begin{pmatrix}0&1\\X&0\end{pmatrix},             \tag{6}
\]
and send the remaining Weyl generators to their scalar amplifications.
These images satisfy all the source Weyl relations. Their classical
residues are scalar, and their entries have finite polynomial support.
This elementary representation alone gives no finite Long lift.

In fact (6) preserves Bernstein PBW degree exactly, assigning \(h\)
degree zero. The top symbols of \(X,p\) are \(XI_2\) and
\(pI_2+hXE_{21}\); the constant term \(hE_{12}\) has lower degree.
For a nonzero scalar homogeneous top symbol \(f_d\), the induced map is
\[
 \operatorname{gr}\iota(f_d)
               =f_d I_2+hX\,\partial_p f_d\,E_{21}.            \tag{7}
\]
Here \(E_{21}^2=0\), and all the symbol variables commute. Its trace
is \(2f_d\ne0\). Thus \(\deg\iota(f)=\deg f\) for every nonzero
finite scalar PBW polynomial, also with Laurent-series scalar
coefficients.

There is a matching assertion for an arbitrary scalar formal series
\(F(h)=\sum_{m\geq0}h^m f_m\), where every \(f_m\) is a finite PBW
polynomial. Its image is well-defined h-adically because (6) reduces
to scalar amplification and introduces no parameter poles. For every
integer \(D\),
\[
 \deg f_m\leq D\ \text{for every }m
 \quad\Longleftrightarrow\quad
 \text{every coefficient of }\iota(F)\text{ has degree at most }D.
                                                               \tag{8}
\]
The forward implication follows because \(\iota\) preserves the
filtration coefficientwise. For the converse, the constant coefficient
of \(\iota(F)\) is \(f_0I_2\), giving the first bound. Inductively
subtract the images of \(f_0,hf_1,\ldots,h^{m-1}f_{m-1}\), which
all retain the same bound. Divide the remainder by \(h^m\) and
reduce: its constant coefficient is \(f_mI_2\), proving the next
bound. The same argument handles any fixed finite parameter pole after
multiplying the series by a power of \(h\).

Consequently composing an existing scalar formal Long lift with (6)
cannot turn unbounded PBW support into bounded support. The connection
is useful here as explicit nonconstant projective data, not as a
boundedness mechanism. A different finite graph-module realization
would require a separate construction; no claim that every realization
with this projective datum is compositional is made.

## What finite admission would still require

The theorem supplies \(E_h\) in an h-complete category. A section can
have polynomial coefficients at every parameter order with no common
PBW bound. An isomorphism \(E_h\simeq E[[h]]\) as parameter-module
sheaves does not give source Weyl-module finiteness or a finite matrix
of target actions. Nor can a complete formal algebra be specialized
at \(h=1\): \(1-h\) is already a unit there.

The missing construction is an actual bimodule over the finite-PBW
source and target Weyl algebras, with this graph reduction and
projective connection, that is finite projective on the source side.
A free source basis would additionally realize the target generator
actions as matrices with finite PBW entries; otherwise a globally
admitted finite projective carrier and its endomorphisms need their own
admission proof. Uniform support bounds, compatible global gluing and
any required nonzero-parameter model remain obligations. Coherence in
the completed category alone does not establish them.

Even such a finite matrix or projective carrier would not establish a
scalar Weyl endomorphism, its nonsurjectivity, or a DC1/DC2 conclusion.
A scalar realization requires a separate argument, and (3) supplies
no invariant rational line for one. Prescribed arithmetic supports
also require independent compatibility; reduction modulo the formal
parameter supplies no such conclusion.

The retained evidence is the pinned primary theorem, its checked graph
hypotheses, and the written all-rational-function proof of (5). No new
checker, bounded lift, algebraization theorem, formal verification,
external review or novelty claim is asserted.
