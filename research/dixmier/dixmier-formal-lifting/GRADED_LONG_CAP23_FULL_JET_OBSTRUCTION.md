# No degree-23 graded full Long jet, with either parameter parity

Use the canonical global polynomial Long tuple
\(f=(R,T,D,S)\) in source variables \((X,Q,P,Z)\), with
\(\{P,X\}=\{Z,Q\}=1\), and the standard Moyal product.
Assign
\[
 \operatorname{wt}(X,Q,P,Z,h)=(-1,1,4,2,3),\qquad
 \operatorname{wt}(R,T,D,S)=(-1,2,4,1).                    \tag{1}
\]

**Theorem.** There is no equivariant scalar full Weyl jet
\[
 F_i=f_i+hF_{1,i}+h^2F_{2,i}+h^3F_{3,i}+h^4F_{4,i}+O(h^5),
 \qquad \deg F_{m,i}\leq23,                              \tag{2}
\]
satisfying all six canonical commutators through \(h^5\).
No parity condition is imposed. All four images may deform,
including every admitted first-image coefficient.

Consequently any full Long lift in the grading (1), with
one common spatial-degree bound, must have that bound at
least 24. No normalization of an arbitrary bounded lift into
(1) is asserted. This gives no non-equivariant obstruction,
no full-parity cap-24 exclusion and no DC1/DC2 resolution.

## The complete lower family used by the proof

The [odd first-jet and nonlinear \(h^3\) theorem](ODD_GRADED_LONG_FIRST_JETS_AND_NORMALIZATION_GAP.md) classifies
the full cap-23 first space, not a selected ansatz. Its two
uniformly cap-preserving source/target shears remove the
\(XQ\) and \(RS\) directions. The remaining \(H_6\) parameter
must vanish at commutator order \(h^3\). Thus every hypothetical
jet (2) can be represented as
\[
 F_1=a\{Y,f\}+b\{R^2T,f\}+c\{H_4,f\},                    \tag{3}
\]
and
\[
 F_2=\Psi(a,b,c)+uK_{R^3}+vK_{R^4S},                    \tag{4}
\]
where
\[
 \begin{split}
 K_{R^3}&=(0,0,-3R^2,0),\\
 K_{R^4S}&=(0,R^4,-4R^3S,0).
 \end{split}
\]
The seven exact coefficient tuples defining \(\Psi\) are
retained in the prior witness. All \(a,b,c,u,v\) are kept.

For notation,
\[
 \begin{split}
 t&=XQ,\quad B=3X^2P+2(1-3XQ)Z,\\
 A_0&=Q+3XQ^2-XB/3,\quad Y=RA_0-3t,\\
 H_4&=Y^4+4tY^2-\frac{16}{3}tY+10t^2.
 \end{split}
\]
The lower theorem also gives \(R_1=R_2=0\), as a consequence
of compatibility. This is not an imposed exact-\(R\) condition.

The only prior normalizations are those proved uniformly
cap-preserving. No \(R^2T\) normalization and no removal of
\(u,v\) is used here.

## Complete coefficient spaces and the commutator formula

The coefficient of \(h^m\) in image weight \(w_i\) has the
complete degree-23 basis
\[
 X^\alpha Q^\beta P^\gamma Z^\delta,\qquad
 \alpha=\beta+4\gamma+2\delta+3m-w_i\geq0,\qquad
 \alpha+\beta+\gamma+\delta\leq23.                        \tag{5}
\]
The [homogeneous tail certificate](GRADED_LONG_TAIL_RIGIDITY.md)
proves that the \(F_3\) map has rank \(205/206\), with sole
kernel
\[
                         K_{R^6}=(0,0,-6R^5,0),
\]
and that the \(F_4\) map is injective on all 130 columns.
Their component counts are respectively
\((34,54,71,47)\) and \((20,34,47,29)\).

For scalar Moyal symbols, the coefficient at \(h^N\) is
\[
 \sum_{\substack{a+b+\ell=N\\ \ell\geq1\ {\rm odd}}}
   \frac{\Pi^\ell(F_{a,i},F_{b,j})}
                       {2^{\ell-1}\ell!}.               \tag{6}
\]
Here \(F_0=f\). The two terms involving the next unknown
coefficient have \(\ell=1\) and one index zero. Every other
term is placed on the negative right side of the corresponding
linear coefficient equation.

In particular the known \(h^4\) defect is
\[
 \{F_{1,i},F_{2,j}\}+\{F_{2,i},F_{1,j}\}
 +\frac{\Pi^3(f_i,F_{1,j})+\Pi^3(F_{1,i},f_j)}{24},       \tag{7}
\]
and the known \(h^5\) defect is
\[
 \begin{split}
 &\{F_{1,i},F_{3,j}\}+\{F_{2,i},F_{2,j}\}
                         +\{F_{3,i},F_{1,j}\}\\
 &+\frac{\Pi^3(f_i,F_{2,j})+\Pi^3(F_{1,i},F_{1,j})
                           +\Pi^3(F_{2,i},f_j)}{24}
                     +\frac{\Pi^5(f_i,f_j)}{1920}.        \tag{8}
 \end{split}
\]
The middle \(\Pi^3(F_1,F_1)\) term is essential. This is
not the even-jet equation.

## Order \(h^4\): one exact relation and the entire \(F_3\) family

Build (7) with every parameter in (3)--(4), and allow all
206 coefficients of \(F_3\) in (5). The resulting row set
has 1563 entries. Its right side has twelve nonzero
parameter-coefficient vectors.

The retained rational functional \(\ell_4\), supported on
196 specified rows, annihilates every coefficient column.
Its value on the right side is exactly
\[
                             c-\frac{27}{544}a.          \tag{9}
\]
Thus compatibility requires
\[
                               a=\frac{544}{27}c.        \tag{10}
\]
This uses exact rational identities, not parameter sampling.

After (10), five exact tuples \(U_c,U_{cv},U_{c^3},U_b,U_{bc^2}\),
retained with all 146 nonzero rational coefficients, give the
particular family
\[
 P_3=cU_c+cvU_{cv}+c^3U_{c^3}+bU_b+bc^2U_{bc^2}.        \tag{11}
\]
Each coefficient identity is verified on the full matrix.
The tail rank \(205/206\) and its exact kernel therefore
make the ENTIRE admitted third coefficient
\[
                               F_3=P_3+wK_{R^6}.        \tag{12}
\]
No third-coefficient branch is omitted. All \(b,c,u,v,w\)
remain parameters. The \(u\)-coefficient in (11) is zero
without any normalization of (4).

The tuples in (11) have zero \(R\)- and \(S\)-components.
Their \(T,D\) degrees are respectively
\[
 \begin{array}{c|cc}
      &T&D\\ \hline
 c    &13&21\\
 cv   &-&22\\
 c^3  &21&23\\
 b    &-&14\\
 bc^2 &-&22
 \end{array}
\]
Thus \(R_3=S_3=0\) is forced, including the homogeneous
direction in (12). None of the first-image vanishings
was a positional normalization hypothesis.

## Order \(h^5\): a universal contradiction with all parameters retained

Substitute (10)--(12) into (8). Allow every one of the
130 degree-23 coefficients of \(F_4\), including all 20
first-image coefficients. The coefficient equation is
\[
 A_4F_4=B_0+vB_v+c^2B_{c^2}+bcB_{bc}
             +c^2vB_{c^2v}+c^4B_{c^4}+bc^3B_{bc^3}.    \tag{13}
\]
This cancellation of \(u,w\) is obtained from the exact
coefficient formula (6), with both parameters still present
in (4),(12). Neither parameter has been normalized away.

The row set has 1207 entries. The retained rational functional
\(\ell_5\), supported on 103 specified rows, satisfies
\[
 \ell_5A_4=0,\qquad \ell_5B_0=1,\qquad
 \ell_5B_v=\ell_5B_{c^2}=\ell_5B_{bc}
 =\ell_5B_{c^2v}=\ell_5B_{c^4}=\ell_5B_{bc^3}=0.         \tag{14}
\]
The replay checks (14) on every full coefficient column
and every right-side coefficient over the rationals.
Applying it to (13) gives \(0=1\), independently of all
five parameters. The identities persist after any
characteristic-zero scalar extension.

The lower classifications exhaust the admitted \(F_1,F_2\);
(9)--(12) exhaust the admitted \(F_3\); (13) allows every
admitted \(F_4\). Higher coefficients cannot alter an
order-\(h^5\) defect. This proves the theorem.

## Certificates, controls and scope

[The witness](graded_long_cap23_full_jet_witness.json) retains
both rational functionals, their row keys and coefficients,
and the complete third-coefficient family. It references
the frozen lower odd and tail data.
[The checker](check_graded_long_cap23_full_jet.py) regenerates
the complete rational equations from the canonical global
Long input, verifies all particular tuples and functional
identities with Fraction arithmetic, and replays the two
homogeneous rank minors with bare integer elimination.
No FLINT or Groebner-basis library is required.

The checker also uses a separate normal-PBW multiplication
engine to check all six commutators through \(h^4\) at
\((b,c,u,v,w)=(-2,1,3,-1,2)\), with (10). This verifies
a nontrivial positive convention/family control; it is
not used to infer the universal contradiction. Deliberately
corrupted third-family and functional data are rejected.
[The receipt](graded_long_cap23_full_jet_controls.json)
pins all required sources and exact data.

This is a scalar finite-coefficient result in the strict
grading (1), allowing both parities and deformed first
images. The result neither places every bounded lift in
this grading nor excludes a full-parity cap-24 lift.
It supplies no general nonlifting theorem, no scalarization
claim for nonscalar carriers, and no DC1/DC2 resolution.
Independent whole-proof replay, formal verification,
external review and novelty priority remain unclaimed.
