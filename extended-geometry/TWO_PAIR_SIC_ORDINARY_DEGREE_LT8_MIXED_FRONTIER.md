# Two-pair `SIC` below ordinary degree eight: mixed-bidegree frontier

## 1. Scope and first theorem

Use contraction pairs
\[
 (\xi _1,z_1),\qquad(\xi _2,z_2)
\]
and write
\[
 \mathcal E_2(\xi^\alpha z^\beta)
 =\partial_z^\alpha z^\beta.
 \tag{1.1}
\]
For a bihomogeneous term of dual--coordinate bidegree \((a,b)\), call
\[
 w=b-a
 \tag{1.2}
\]
its central weight.  Unequal weights in a power can cancel, so the
balanced scalar moments alone do not describe a nonhomogeneous form.

This note begins the ordinary-degree-\(<8\) search with the first proposed
three-block mixture.

> **Theorem 1.1.** Let
> \[
> F=A+B+C,\qquad
> A\in V_{2,2},\quad B\in V_{1,3},\quad C\in V_{3,1}.
> \tag{1.3}
> \]
> If
> \[
> \mathcal E_2(F^m)=0\qquad(m\geq1),
> \tag{1.4}
> \]
> then, for every fixed polynomial \(Q\),
> \[
> \mathcal E_2(QF^m)=0\qquad(m\gg0).
> \tag{1.5}
> \]
> Thus this complete mixed stratum satisfies `SIC(2)`.
>
> If \(B\ne0\), the first four pure contractions already force the
> one-sided form used to prove (1.5).  If \(B=0\), the result follows from
> the complete bidegree-\((2,2)\) theorem.

The same argument without a balanced middle block gives an all-degree
safe family.

> **Theorem 1.2.** For every \(d\geq2\), `SIC(2)` holds on
> \[
> V_{1,d}\oplus V_{d,1}.
> \tag{1.6}
> \]
> In particular all five opposite dual-linear pair types of ordinary
> degrees three through seven are excluded from the present search.

To see this, normalize a nonzero positive block to
\(B=\xi _2z_1^d\).  In the negative block, the unique monomial with
negative epsilon from (5.1) is \(h\xi _1^dz_2\).  The same weight--epsilon
argument used in Section 4 shows that
\[
\begin{aligned}
\mathcal E_2\bigl((B+C)^{2r-1}\bigr)&=0,\\
\mathcal E_2\bigl((B+C)^{2r}\bigr)
&=\binom{2r}{r}h^r(dr)!\,r!.
\end{aligned}
\tag{1.7}
\]
The second moment forces \(h=0\).  Every residual monomial has
\(\epsilon\geq0\); its epsilon-zero face has
\(\delta=-(d-1)\), while positive-epsilon terms have \(\delta\leq d\).
For a multiplier monomial \(Q\), with
\(R_Q=\max(0,-\epsilon(Q))\), eventual vanishing follows from
\[
 m>
 \frac{\delta(Q)+(2d-1)R_Q}{d-1}.
\tag{1.8}
\]
If the positive block is zero, its negative weight gives eventual
vanishing directly.

The next balanced level has the same answer away from its old balanced
boundary.

> **Theorem 1.3.** Let
> \[
> F=A+B+C,\qquad
> A\in V_{3,3},\quad B\in V_{1,4},\quad C\in V_{4,1},
> \tag{1.9}
> \]
> and assume \(B\ne0\).  If the first five pure contractions of \(F\)
> vanish, then
> \[
> \mathcal E_2(QF^m)=0\qquad(m\gg0)
> \tag{1.10}
> \]
> for every fixed multiplier \(Q\).
>
> On the boundary \(B=0\), the pure premise and mixed conclusion reduce
> exactly to the balanced bidegree-\((3,3)\) problem.  Hence the complete
> stratum in (1.9) is `SIC`-safe if and only if \(V_{3,3}\) is safe.

The first positive block of dual degree two is also closed on both pure
irreducible summands.

> **Theorem 1.4.** Let
> \[
> F=B+C,\qquad B\in V_{2,3},\quad C\in V_{3,2}.
> \tag{1.11}
> \]
> Under
> \[
> V_{2,3}\cong\operatorname{Sym}^5\oplus
> \operatorname{Sym}^3\oplus\operatorname{Sym}^1,
> \tag{1.12}
> \]
> suppose the nonzero positive block \(B\) lies entirely in either
> \(\operatorname{Sym}^5\) or \(\operatorname{Sym}^3\).  If the first
> four pure contractions of \(F\) vanish, then every fixed mixed
> contraction vanishes eventually.
>
> Thus a candidate in \(V_{2,3}\oplus V_{3,2}\) must have both its
> \(\operatorname{Sym}^5\) and \(\operatorname{Sym}^3\) components
> nonzero.  The first contraction always removes
> \(\operatorname{Sym}^1\).

The first mixed-positive local test is also rigid beyond its misleading
tangent space.

> **Proposition 1.5.** At the explicit mixed
> \(\operatorname{Sym}^5+\operatorname{Sym}^3\) point (10.17), the
> moment-one-through-four Jacobian has a two-dimensional transverse
> tangent excess.  One excess direction is obstructed at deformation
> order two.  The other lifts through order three but has no
> order-four lift, even after allowing every second- and third-order
> correction.  This is a local jet exclusion, not a classification of
> the mixed-positive locus.

The local calculation admits a one-parameter generic version.

> **Proposition 1.6.** On the chart \(b(a-b)\ne0\), the incidence family
> \[
> VZ^2(aVY+bWZ+cVZ)
> \]
> reduces under contraction-preserving \(\mathrm{GL}_2\) to
> \(VZ^2(uVY+WZ)\), \(u=a/b\).  Away from the explicit exceptional
> factors in (10.28), both transverse directions of the first four
> moment equations are obstructed by deformation order four.  Direct
> recomputation also closes \(u=0,-6,-3/4,\infty\).  The remaining
> algebraic ratios are listed in Section 10.5.

The theorem is about the full coefficient spaces in (1.3), not a sparse
support census.  Its proof uses polynomial-valued contraction
coefficients of weights two and four; a scalar-moment calculation would
miss the decisive equations.

It does **not** prove that ordinary degree eight is minimal.  In
particular, balanced mixtures through \((3,3)\), other two-sided weight
collections, and the locus where the \(V_{2,3}\) positive block has both
\(\operatorname{Sym}^5\) and \(\operatorname{Sym}^3\) components remain
open.

## 2. The degree-seven bidegree ledger

There are \(35\) nonconstant bidegrees
\[
 (a,b),\qquad a,b\geq0,\quad 1\leq a+b\leq7.
 \tag{2.1}
\]
They split into
\[
 16\text{ positive-weight},\qquad
 3\text{ balanced},\qquad
 16\text{ negative-weight}
 \tag{2.2}
\]
blocks.  The balanced blocks are
\[
 (1,1),\quad(2,2),\quad(3,3).
 \tag{2.3}
\]
The checker records the complete ledger, including ordinary degree and
central weight, in the generated artifact.

There are seven nonempty balanced-only block collections.  A genuinely
two-sided collection is a nonempty choice from the \(16\) positive blocks,
a nonempty choice from the \(16\) negative blocks, and an arbitrary choice
of balanced blocks, hence there are
\[
 (2^{16}-1)^2\,2^3=34{,}358{,}689{,}800
 \tag{2.4}
\]
labelled block collections.  Listing these subsets carries no extra
information: every such collection contains a positive/negative pair.
There are \(16\cdot16=256\) primitive pair types.  If their weights are
\(r>0\) and \(-s<0\), their least central product uses
\[
 \frac{s}{\gcd(r,s)}\quad\text{positive factors and}\quad
 \frac{r}{\gcd(r,s)}\quad\text{negative factors}.
 \tag{2.5}
\]
The generated artifact records all \(256\) circuits and their resulting
balanced bidegrees.  This is the finite circuit enumeration relevant to
central cancellation.

There is a much smaller nonlinear core.  If both the dual and coordinate
degrees of both blocks are at least two, only sixteen primitive circuits
remain.  Grouped by their number of factors, they are
\[
\begin{array}{c|l}
2&
(2,3)+(3,2),\ (2,4)+(4,2),\
(2,3)+(4,3),\ (3,4)+(3,2),\\
& (2,5)+(5,2),\ (3,4)+(4,3)\\[2pt]
3&
2(2,3)+(4,2),\ (2,4)+2(3,2),\
(2,4)+2(4,3),\ 2(3,4)+(4,2)\\[2pt]
4&
3(2,3)+(5,2),\ (2,5)+3(3,2),\
(2,5)+3(4,3),\ 3(3,4)+(5,2)\\[2pt]
5&
3(2,4)+2(5,2),\ 2(2,5)+3(4,2).
\end{array}
\tag{2.6}
\]
Here a coefficient denotes the multiplicity in the least central-weight
cancellation.  The first row contains the only two-factor nonlinear
circuits, and \((2,3)+(3,2)\) is the unique one whose primitive central
product has bidegree \((5,5)\); all others begin at \((6,6)\) or higher.
This makes Section 10 the correct first nonlinear gate rather than merely
one convenient example.  Circuits involving a degree-zero or degree-one
side still require the separate linear-side reductions; (2.6) is not by
itself an exclusion theorem for them.

A simultaneous contraction-preserving
\(\mathrm{GL}_2\) change does not alter \((a,b)\); it acts inside each
coefficient block.  Consequently quotienting by \(\mathrm{GL}_2\) is an
orbit problem on the chosen direct sum, not a further identification of
different ledger entries.

There is a useful first pruning rule.  If all nonzero blocks have
nonpositive weight, then
\[
 \mathcal E_2(F^m)=\mathcal E_2(F_0^m),
 \tag{2.7}
\]
where \(F_0\) is the balanced part.  For a fixed multiplier, only
boundedly many negative-weight factors can occur in a surviving
contraction.  Hence such a collection reduces to the corresponding
balanced nonhomogeneous problem.  The analogous statement holds with
the roles organized by the maximal positive-weight equations.  A genuinely
new central-cancellation mechanism must therefore use both weight signs,
unless it already lies in the balanced mixture
\[
 V_{1,1}\oplus V_{2,2}\oplus V_{3,3}.
 \tag{2.8}
\]

## 3. Normalizing the \((1,3)\) block

Assume \(B\ne0\).  The coefficient of the largest central weight in
\(\mathcal E_2(F^m)\) is \(\mathcal E_2(B^m)\).  In particular the first
two such equations give
\[
 \mathcal E_2(B)=\mathcal E_2(B^2)=0.
 \tag{3.1}
\]
The dual-linear theorem puts a nonzero bihomogeneous \(B\) into the form
\[
 B=(b\xi _1-a\xi _2)(az_1+bz_2)^3.
 \tag{3.2}
\]
After a contraction-preserving \(\mathrm{GL}_2\) change and an overall
nonzero scaling of \(F\), take
\[
 B=\xi _2z_1^3.
 \tag{3.3}
\]

Write
\[
\begin{aligned}
A&=\sum_{0\leq i,j\leq2}
 a_{ij}\xi _1^i\xi _2^{2-i}z_1^jz_2^{2-j},\\
C&=\sum_{\substack{0\leq i\leq3\\0\leq j\leq1}}
 c_{ij}\xi _1^i\xi _2^{3-i}z_1^jz_2^{1-j}.
\end{aligned}
\tag{3.4}
\]
The top two positive-weight ladders give
\[
\begin{aligned}
\mathcal E_2(A)
 &=2a_{00}+a_{11}+2a_{22},\\
\mathcal E_2(AB)
 &=(6a_{10}+12a_{21})z_1^2+12a_{20}z_1z_2,\\
\mathcal E_2(AB^2)
 &=60a_{20}z_1^4.
\end{aligned}
\tag{3.5}
\]
Thus
\[
 a_{20}=0,\qquad a_{10}=-2a_{21},\qquad
 a_{11}=-2a_{00}-2a_{22}.
\tag{3.6}
\]

Now use the full third contraction.  Its positive-weight coefficients
after (3.6) include
\[
\begin{aligned}
[z_1^2]\mathcal E_2(F^3)
 &=-144a_{21}(a_{00}-7a_{22}),\\
[z_1z_2]\mathcal E_2(F^3)
 &=720a_{21}^2.
\end{aligned}
\tag{3.7}
\]
Characteristic zero therefore gives
\[
 a_{21}=a_{10}=a_{20}=0.
\tag{3.8}
\]
This is the step that is invisible to scalar moments.

## 4. The exact central core

After (3.8), every nonzero contraction of \(F^m\) has weight zero, and it
depends only on
\[
 p=a_{00},\qquad q=a_{22},\qquad h=c_{30}.
\tag{4.1}
\]
Indeed, use \(\epsilon=i_2-j_2\) from (5.1) below.  The normalized \(B\)
has \((w,\epsilon)=(2,1)\).  Every monomial of \(C\) has weight \(-2\);
the unique monomial \(c_{30}\xi _1^3z_2\) has \(\epsilon=-1\), while all
other \(C\)-monomials have \(\epsilon\geq0\).  After (3.8), every
\(A\)-monomial has \(\epsilon\geq0\), with equality exactly on its three
diagonal terms.

If a product contains \(r\) factors \(B\) and \(s\) factors \(C\), its
weight is \(2(r-s)\) and its epsilon is at least \(r-s\).  A surviving
contraction needs nonnegative weight and nonpositive epsilon.  Hence
\(r=s\), equality must hold in the epsilon bound, every \(C\)-factor is
the \(c_{30}\) corner, and every \(A\)-factor is diagonal.  This proves
the asserted all-order core reduction.

Equivalently, all pure contractions agree with those of
\[
\begin{aligned}
F_{\mathrm{core}}
={}&p\,\xi _2^2z_2^2
-2(p+q)\xi _1\xi _2z_1z_2
+q\,\xi _1^2z_1^2\\
&+\xi _2z_1^3+h\,\xi _1^3z_2.
\end{aligned}
\tag{4.2}
\]
The first moment is zero.  After removing nonzero scalar factors, the
next three are
\[
\begin{aligned}
f_2={}&3h+4p^2-2pq+4q^2,\\
f_3={}&-5hp+2hq+4p^3-2p^2q-2pq^2+4q^3,\\
f_4={}&5h^2+9hp^2-6hpq+15hq^2\\
&\quad+12p^4-6p^3q+6p^2q^2-6pq^3+12q^4,
\end{aligned}
\tag{4.3}
\]
with
\[
 \mathcal E_2(F_{\mathrm{core}}^2)=4f_2,\quad
 \mathcal E_2(F_{\mathrm{core}}^3)=72f_3,\quad
 \mathcal E_2(F_{\mathrm{core}}^4)=1728f_4.
\tag{4.4}
\]

An exact lexicographic Gröbner basis of
\((f_2,f_3,f_4)\subset\mathbb Q[p,q,h]\) contains
\[
 h^3.
\tag{4.5}
\]
After adjoining \(h\), its basis contains \(q^4\); then \(f_2\) forces
\(p=0\).  Hence
\[
 \boxed{\sqrt{(f_2,f_3,f_4)}=(p,q,h).}
\tag{4.6}
\]
All three central coefficients vanish.

## 5. Eventual mixed contraction

For a monomial
\[
 M=\xi _1^{i_1}\xi _2^{i_2}z_1^{j_1}z_2^{j_2},
\]
put
\[
 \delta(M)=j_1-i_1,\qquad
 \epsilon(M)=i_2-j_2.
\tag{5.1}
\]
A monomial can survive \(\mathcal E_2\) only if
\[
 \delta(M)\geq0,\qquad\epsilon(M)\leq0.
\tag{5.2}
\]

After (4.6), every monomial left in \(F\) has
\[
 \epsilon\geq0.
\tag{5.3}
\]
Those with \(\epsilon>0\) have \(\delta\leq3\), while every monomial with
\(\epsilon=0\) has
\[
 \delta=-2.
\tag{5.4}
\]
This is a strict lexicographic contraction cone.

Fix a multiplier monomial \(Q\), and let
\[
 R_Q=\max(0,-\epsilon(Q)).
\tag{5.5}
\]
In a possibly surviving monomial of \(QF^m\), at most \(R_Q\) factors
from \(F\) can have positive \(\epsilon\).  Therefore
\[
 \delta(QF^m)
 \leq\delta(Q)+3R_Q-2(m-R_Q)
 =\delta(Q)+5R_Q-2m.
\tag{5.6}
\]
For
\[
 m>\frac{\delta(Q)+5R_Q}{2},
\tag{5.7}
\]
condition (5.2) fails.  Taking the maximum of this finite bound over the
monomials of a polynomial multiplier proves (1.5).

## 6. The branch \(B=0\)

If \(B=0\), every term involving \(C\) in \(F^m=(A+C)^m\) has negative
central weight and contracts to zero.  Thus
\[
 \mathcal E_2(F^m)=\mathcal E_2(A^m).
\tag{6.1}
\]
The bidegree-\((2,2)\) theorem says that the pure-moment premise for \(A\)
implies eventual vanishing for every fixed multiplier.

For a fixed \(Q\), only boundedly many factors \(C\) can occur in a
monomial of \(Q(A+C)^m\) having nonnegative central weight.  Apply the
bidegree-\((2,2)\) theorem separately to the finitely many fixed
multipliers \(QC^k\).  This proves (1.5) also on the \(B=0\) branch and
completes Theorem 1.1.

## 7. The nonzero \(V_{1,4}\) branch

Normalize \(B=\xi _2z_1^4\), and write
\[
 A=\sum_{0\leq i,j\leq3}
 a_{ij}\xi _1^i\xi _2^{3-i}z_1^jz_2^{3-j}.
 \tag{7.1}
\]
The maximal positive-weight contractions
\(\mathcal E_2(AB^n)\), \(0\leq n\leq3\), give
\[
\begin{gathered}
a_{30}=a_{31}=a_{20}=0,\\
a_{10}=-5a_{32}-\frac53a_{21},\\
a_{00}=-\frac13a_{22}-a_{33}-\frac13a_{11}.
\end{gathered}
\tag{7.2}
\]
The weight-three coefficient of the third contraction and weight-six
coefficient of the fourth contain, respectively,
\[
\begin{aligned}
q_1&=118a_{32}^2+51a_{32}a_{21}+6a_{21}^2,\\
q_2&=72a_{32}^2+21a_{32}a_{21}+2a_{21}^2.
\end{aligned}
\tag{7.3}
\]
Their exact lexicographic basis is
\[
199a_{32}^2-4a_{21}^2,\qquad
597a_{32}a_{21}+98a_{21}^2,\qquad
a_{21}^3.
\tag{7.4}
\]
Thus \(a_{32}=a_{21}=a_{10}=0\): the balanced block is upper
triangular in the flag selected by \(B\).

The weight--epsilon argument from Section 4 now shows that every pure
contraction depends only on the diagonal of \(A\) and the opposite corner
\(h\xi _1^4z_2\) of \(C\).  Put
\[
\begin{aligned}
A_{\mathrm{diag}}
={}&-\frac{p+q+3r}{3}(\xi _2z_2)^3
+p(\xi _1z_1)(\xi _2z_2)^2\\
&+q(\xi _1z_1)^2(\xi _2z_2)
+r(\xi _1z_1)^3.
\end{aligned}
\tag{7.5}
\]
For
\[
F_{\mathrm{core}}
=A_{\mathrm{diag}}+\xi _2z_1^4+h\xi _1^4z_2,
\tag{7.6}
\]
the first moment is zero.  Removing nonzero scalar contents, the next
two moments begin
\[
\begin{aligned}
f_2={}&2p^2+5pq+13pr+4q^2+25qr+57r^2+2h,\\
f_3={}&-20p^3-70p^2q-230p^2r-80pq^2-520pqr\\
&-780pr^2+11ph-30q^3-290q^2r-780qr^2\\
&+26qh+93rh.
\end{aligned}
\tag{7.7}
\]
The fourth and fifth polynomials are recorded in the generated artifact.
An exact rational Gröbner basis of
\[
 (f_2,f_3,f_4,f_5)\subset\mathbb Q[p,q,r,h]
\tag{7.8}
\]
has \(28\) elements and gives the power memberships
\[
 p^{10},q^{10},r^{10},h^5\in(f_2,f_3,f_4,f_5).
\tag{7.9}
\]
Therefore its radical is the coefficient origin.

After this core vanishes, every residual monomial has
\(\epsilon\geq0\).  The epsilon-zero face has \(\delta=-3\), while
positive-epsilon terms have \(\delta\leq4\).  For a multiplier monomial
\(Q\), put \(R_Q=\max(0,-\epsilon(Q))\).  Then no contraction survives
once
\[
 m>\frac{\delta(Q)+7R_Q}{3}.
\tag{7.10}
\]
This proves Theorem 1.3.

If \(B=0\), negative-weight factors from \(C\) occur only boundedly in a
surviving mixed contraction, exactly as in Section 6.  The remaining
powers are powers of \(A\in V_{3,3}\).  Conversely \(C=0\) recovers the
balanced problem.  This proves the claimed equivalence of the boundary
with the existing bidegree-\((3,3)\) frontier.

## 8. Remaining degree-\(<8\) search

The theorems close the first requested two-sided mixture, every opposite
dual-linear pair, and all new points of
\(V_{3,3}\oplus V_{1,4}\oplus V_{4,1}\) away from its balanced boundary.
The main
unresolved collection classes are:

1. the balanced nonhomogeneous sum
   \(V_{1,1}\oplus V_{2,2}\oplus V_{3,3}\), including a \((3,3)\) top
   block with lower balanced corrections;
2. two-sided collections whose smallest positive block has dual degree
   two or three.  In \(V_{2,3}\oplus V_{3,2}\), Theorem 1.4 leaves only
   the positive blocks with simultaneous nonzero
   \(\operatorname{Sym}^5\) and \(\operatorname{Sym}^3\) components.
   Propositions 1.5--1.6 close the generic transverse jet and reduce the
   remaining incidence calculation to nine algebraic ratios and the
   \(a=b\) stratum;
3. collections with several positive weights, whose maximal-weight
   equations must first be solved as full polynomial identities;
4. Long-like products or sums not contained in (1.3).

The circuit table (2.6) gives a natural escalation rule: close the unique
primitive \((5,5)\) circuit before attacking the five remaining
two-factor nonlinear circuits, all of which first balance in bidegree at
least \((6,6)\).  Higher-multiplicity circuits should be deferred until
their two-factor faces are understood, because their first central
equation is cubic or worse and has a larger correction space.

Any finite search on these classes must retain every output coefficient
of \(\mathcal E_2(F^m)\).  A long zero prefix is only experimental
evidence until it is replaced by a proved recurrence, constant-term
identity, or finite-difference formula.

## 9. Uniform diagonal core and the degree-eight threshold

The proofs for \(d=2,3\) expose a reusable two-gate program for
\[
 V_{d,d}\oplus V_{1,d+1}\oplus V_{d+1,1}
 \tag{9.1}
\]
on the nonzero positive branch:

1. use polynomial-valued positive-weight contractions to force the
   \(V_{d,d}\) matrix upper triangular;
2. test the diagonal/opposite-corner core by scalar factorial moments.

The reduction itself is uniform.

> **Proposition 9.1 (two-gate reduction).** Normalize the positive block
> to \(B=\xi _2z_1^{d+1}\).  Suppose the positive-weight contractions
> force
> \[
> A=\sum_{0\leq i\leq j\leq d}
> a_{ij}\xi _1^i\xi _2^{d-i}z_1^jz_2^{d-j}.
> \tag{9.2}
> \]
> Then every pure contraction of \(A+B+C\) equals that of the diagonal
> of \(A\), \(B\), and the single opposite corner
> \(h\xi _1^{d+1}z_2\) of \(C\).  If the scalar moments force this core
> to vanish, every multiplier monomial \(Q\) has the explicit cutoff
> \[
> m>
> \frac{\delta(Q)+(2d+1)\max(0,-\epsilon(Q))}{d}.
> \tag{9.3}
> \]

Indeed, \(B\) has \((w,\epsilon)=(d,1)\); every \(C\)-monomial has
weight \(-d\), with epsilon at least \(-1\), and equality occurs only at
the opposite corner.  Every upper-triangular \(A\)-monomial has
\(\epsilon\geq0\), with equality only on the diagonal.  The same
weight--epsilon inequalities used in Sections 4 and 7 prove the core
reduction.  After the core vanishes, the epsilon-zero face has
\(\delta=-d\), while positive-epsilon terms have \(\delta\leq d+1\);
this gives (9.3).

The second gate has a compact formula.  Put
\[
\begin{aligned}
X&=\xi _1z_1,&Y&=\xi _2z_2,\\
D&=\sum_{i=0}^d c_iX^iY^{d-i}
  +\xi _2z_1^{d+1}+h\xi _1^{d+1}z_2,
\end{aligned}
\tag{9.4}
\]
with \(\sum_i c_i\,i!(d-i)!=0\).  In a central term of \(D^m\), the two
unbalanced corners occur \(k\) times each.  If the diagonal multiplicities
are \(n_0,\ldots,n_d\), with \(\sum_i n_i=m-2k\), its exact contribution
is
\[
\frac{m!}{k!^2\prod_i n_i!}\,
\left((d+1)k+\sum_i in_i\right)!\,
\left(k+\sum_i(d-i)n_i\right)!\,
h^k\prod_i c_i^{n_i}.
\tag{9.5}
\]
Thus no four-variable Wick expansion is needed; the core is a finite
factorial-moment ideal in \(d+1\) variables.

At \(d=4\), moments two through six have exact Jacobian rank five at
\((c_1,c_2,c_3,c_4,h)=(1,2,3,4,5)\).  Exact Gröbner calculations over
both \(\mathbb F_{101}\) and \(\mathbb F_{1009}\) give:
\[
\text{basis size }132,\qquad
\dim_{\mathbb F_p}\mathbb F_p[c_1,c_2,c_3,c_4,h]/I=360,
\tag{9.6}
\]
and reduce
\[
c_1^{40},c_2^{40},c_3^{40},c_4^{40},h^{20}
\tag{9.7}
\]
to zero.  Hence the modular zero fibers are supported at the origin.

This is finite-field evidence only.  The rational Gröbner calculation
exhausted the available local memory, and no reconstructed
characteristic-zero power certificates are claimed.  The preceding
\(d=4\) upper-triangularization gate is also open.  The evidence suggests
that the known ordinary-degree-eight witness is genuinely off-diagonal
rather than a failure of the small diagonal core.  The clean next
strengthening is modular-basis reconstruction of (9.6)--(9.7), followed
by the full polynomial-valued \(d=4\) triangular gate.

## 10. The pure summands of \(V_{2,3}\oplus V_{3,2}\)

Use variables \((W,Z),(V,Y)\) and
\[
 R=WZ+VY.
 \tag{10.1}
\]
The Clebsch--Gordan basis obtained by lowering the three highest vectors
\[
 V^2Z^3,\qquad VZ^2R,\qquad ZR^2
\tag{10.2}
\]
has determinant \(-22500\) in the monomial basis.  Write its coordinates
as
\[
 (s_0,\ldots,s_5;\ t_0,\ldots,t_3;\ r_0,r_1).
\]
The first contraction is
\[
 \mathcal E_2(B)=12(r_0Z+r_1Y),
 \tag{10.3}
\]
so \(r_0=r_1=0\) in every candidate.

### 10.1 Pure \(\operatorname{Sym}^3\)

On the cubic summand, the second contraction is \(-336\) times the binary
quadratic whose coefficients are
\[
\begin{aligned}
t_0t_2-t_1^2,\qquad
t_0t_3-t_1t_2,\qquad
t_1t_3-t_2^2.
\end{aligned}
\tag{10.4}
\]
These are the Hankel minors of the rational normal cubic.  Hence
\[
 (t_0,t_1,t_2,t_3)
 =\lambda(\alpha^3,\alpha^2\beta,\alpha\beta^2,\beta^3).
\tag{10.5}
\]
Every nonzero point is one orbit, with normal form
\[
 B_3=VZ^2R.
\tag{10.6}
\]

### 10.2 Pure \(\operatorname{Sym}^5\)

On the quintic summand, the second contraction is \(72\) times the binary
quadratic with coefficients
\[
\begin{aligned}
q_0&=s_0s_4-4s_1s_3+3s_2^2,\\
q_1&=s_0s_5-3s_1s_4+2s_2s_3,\\
q_2&=s_1s_5-4s_2s_4+3s_3^2.
\end{aligned}
\tag{10.7}
\]
Exact rational primary decomposition proves that this ideal is prime of
affine dimension three.  It is the tangential variety of the rational
normal quintic:
\[
 \sum_{i=0}^5\binom5i s_iX^{5-i}T^i=L^4M.
\tag{10.8}
\]
The parameterization \(L=aX+bT,\ M=cX+dT\) is
\[
\begin{aligned}
s_0&=a^4c,&
s_1&=\frac{a^3(ad+4bc)}5,&
s_2&=\frac{a^2b(2ad+3bc)}5,\\
s_3&=\frac{ab^2(3ad+2bc)}5,&
s_4&=\frac{b^3(4ad+bc)}5,&
s_5&=b^4d.
\end{aligned}
\tag{10.9}
\]
There are two orbit types:
\[
 B_{5,\mathrm{power}}=V^2Z^3,\qquad
 B_{5,\mathrm{tan}}=VZ^2(3VY-2WZ).
\tag{10.10}
\]

### 10.3 Adding the negative block

For \(B_3\) or \(B_{5,\mathrm{tan}}\), only six monomials of
\(C\in V_{3,2}\) can initially meet the central grades.  Order them by
\[
 (i,j)=(1,0),(2,0),(2,1),(3,0),(3,1),(3,2)
\tag{10.11}
\]
in \(W^iV^{3-i}Z^jY^{2-j}\), with coefficients
\((u_0,\ldots,u_5)\).  For \(B_3\), moments two through four reduce,
up to nonzero scalars, to
\[
\begin{aligned}
2u_1+3u_4,\qquad u_3,\qquad
4u_0u_3+2u_1^2+5u_1u_4+5u_2u_3+10u_3u_5+5u_4^2.
\end{aligned}
\tag{10.12}
\]
Their Gröbner basis contains \(u_4^2,2u_1+3u_4,u_3\).
For \(B_{5,\mathrm{tan}}\), the corresponding equations are
\[
\begin{aligned}
u_1-u_4,\qquad u_3,\qquad
18u_0u_3+9u_1^2+15u_1u_4+15u_2u_3+40u_3u_5+20u_4^2,
\end{aligned}
\tag{10.13}
\]
with basis \(u_4^2,u_1-u_4,u_3\).  Thus in either case
\[
 u_1=u_3=u_4=0.
\tag{10.14}
\]
Every residual monomial has \(\epsilon\geq0\); its epsilon-zero face has
\(\delta=-1\), while positive-epsilon terms have \(\delta\leq2\).
Consequently a multiplier monomial \(Q\) has cutoff
\[
 m>\delta(Q)+3\max(0,-\epsilon(Q)).
\tag{10.15}
\]

For \(B_{5,\mathrm{power}}\), moment two removes the unique opposite
corner \(W^3Y^2\).  Every remaining \(C\)-monomial has
\(\epsilon\geq-1\), whereas \(B_{5,\mathrm{power}}\) has
\((w,\epsilon)=(1,2)\).  If \(w_Q\) is the central weight of \(Q\) and
\[
 R_Q=\max(0,w_Q-\epsilon(Q)),
\]
a surviving term has
\[
 m\leq2R_Q+\max(w_Q,0).
\tag{10.16}
\]
This proves Theorem 1.4.

### 10.4 A mixed-positive fourth-order obstruction

The remaining branch has a deceptive tangent excess.  Consider
\[
 B_0=VZ^2(2VY+3WZ+5VZ).
\tag{10.17}
\]
The full polynomial-valued moments of orders one through four have
Jacobian rank six in the twelve-dimensional space \(V_{2,3}\).  The
obvious incidence family
\[
 VZ^2(aVY+bWZ+cVZ)
\tag{10.18}
\]
together with the \(\mathrm{SL}_2\) lowering direction has tangent
dimension four.  A complementary pair of tangent-kernel directions is
\[
\begin{aligned}
 H_1&=\frac13WZ^2(3VY-WZ),\\
 H_2&=\frac1{20}Y(20V^2Y^2-97VWYZ+37W^2Z^2).
\end{aligned}
\tag{10.19}
\]
For \(H=\alpha H_1+\beta H_2\), the complete second-order solvability
ideal is
\[
 (\beta^2).
\tag{10.20}
\]
Thus \(H_2\) is obstructed, whereas \(H_1\) survives.  One explicit lift
of the latter begins
\[
 B(t)=B_0+tH_1+t^2K_2+t^3K_3+O(t^4),
\tag{10.21}
\]
where
\[
\begin{aligned}
 K_2&=\frac1{10}WYZ(VY-WZ),\\
 K_3&=\frac1{30280}WY
 (4VY^2-169VYZ-172WYZ+169WZ^2).
\end{aligned}
\tag{10.22}
\]
All four moments vanish modulo \(t^4\).

This lift does not continue.  Allow all six tangent-kernel parameters in
the second correction.  Every such correction admits a third correction,
again with six free parameters.  At order four the eight compatibility
functionals reduce to two quadratics in a single surviving second-order
parameter \(a\):
\[
\begin{aligned}
 p(a)&=198510381a^2+6275108a+44436,\\
 q(a)&=1311046029a^2+20142212a+220884.
\end{aligned}
\tag{10.23}
\]
They are independent of all six third-order parameters, and
\[
 \operatorname{Res}_a(p,q)
 =2283980165392458318151680000\ne0.
\tag{10.24}
\]
Hence \(H_1\) has no fourth-order lift.  This closes both transverse
directions at the one point (10.17), but it does not exclude a different
component or an exceptional parameter ratio elsewhere on the
mixed-positive locus.  The strongest next calculation is to retain
\((a,b,c)\) in (10.18), compute the symbolic analogue of (10.24), and
separate its generic nonvanishing from its exceptional divisor.  Section
10.5 carries this out after first reducing the three parameters to one.

### 10.5 The one-parameter generic gate

The three incidence coefficients have only one essential parameter on
the principal chart.  The contraction-preserving raising operator
\[
 \mathcal R=Z\partial_Y-V\partial_W
\tag{10.25}
\]
satisfies
\[
\mathcal R\!\left(VZ^2(aVY+bWZ+cVZ)\right)
=(a-b)V^2Z^3,\qquad \mathcal R^2=0
\tag{10.26}
\]
on this family.  Thus \(c\) is removed when \(a\ne b\), and scaling
\(b\ne0\) leaves
\[
 B(u)=VZ^2(uVY+WZ),\qquad u=a/b.
\tag{10.27}
\]

Over \(\mathbf Q(u)\), the moment-one-through-four Jacobian again has
rank six, with a four-dimensional incidence tangent and two transverse
directions.  The second-order compatibility equations for one transverse
coordinate are two coprime coefficient polynomials in \(u\); their
resultant is
\[
 -10554055205970310272.
\]
Hence that direction is uniformly obstructed on every valid chart.

The other direction admits all six second-order corrections and every
one of them admits six-parameter third-order corrections.  The two
fourth-order compatibility equations are quadratic in one surviving
correction parameter.  Their resultant is
\[
648u^8(u+6)^4(4u+3)^2D(u)^2S(u),
\tag{10.28}
\]
where
\[
\begin{aligned}
D(u)&=2u^3+10u^2+21u+9,\\
S(u)&=14138u^6+142955u^5+483945u^4+727020u^3\\
&\qquad+540270u^2+185004u+18225.
\end{aligned}
\tag{10.29}
\]
Both \(D\) and \(S\) are irreducible over \(\mathbf Q\).
The denominators used to choose the rational-function kernel basis are
\[
 (u+6)^2(4u+3)D(u).
\tag{10.30}
\]
Consequently the generic calculation leaves the roots of \(D S\), as
well as the apparent ratios \(u=0,-6,-3/4\).

Direct characteristic-zero recomputation, without the rational-function
basis, gives unit fourth-order obstruction ideals at
\[
 u=0,\qquad u=-6,\qquad u=-3/4,\qquad u=\infty.
\tag{10.31}
\]
It also gives a unit ideal at the representative
\(VZ^2(VY+WZ+VZ)\) of the first \(a=b,\ c\ne0\) test.
Thus every rational ratio on the normalized \(a\ne b\) incidence family
is locally excluded through deformation order four.  Over an algebraic
closure, the three roots of \(D\) require direct chart recomputation, and
the six roots of \(S\) are the genuine fourth-order exceptional ratios.
The \(a=b\) stratum still needs its complete orbit calculation.

The sextic branch is already rigid enough to make the fifth-order test
univariate.  The linear subresultant of the two compatibility quadratics
is, up to the nonzero chart factor in (10.28),
\[
 A(u)\tau+B(u),
\tag{10.32}
\]
where
\[
\begin{aligned}
A(u)&=20u^7+1262u^6+10842u^5+31626u^4+38295u^3\\
&\qquad+21168u^2+6075u+2187,\\
B(u)&=52u^6+378u^5+581u^4+42u^3+27u^2+243u.
\end{aligned}
\tag{10.33}
\]
Moreover
\[
\operatorname{Res}(A,S)
=98518095894778815317044086562668309937383692646528\ne0.
\tag{10.34}
\]
Thus at every root of \(S\), the fourth-order lift has the unique
correction coordinate \(\tau=-B/A\).  No further correction-parameter
elimination is needed before extracting the fifth coefficient.

This is a stronger and simpler next frontier than the original
three-parameter elimination: inspect nine algebraic ratios and one
codimension-one normalization stratum, rather than a general incidence
family.  At a root of \(S\), the next decisive test is the fifth
deformation coefficient with \(\tau=-B/A\) adjoined; at a
root of \(D\), first change the kernel chart before interpreting the
vanishing.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_two_pair_sic_mixed_22_13_31.py
.venv/bin/python scripts/verify_two_pair_sic_mixed_33_14_41.py
.venv/bin/python scripts/explore_two_pair_sic_mixed_diagonal_core.py
.venv/bin/python scripts/verify_two_pair_sic_mixed_23_32_pure_summands.py
.venv/bin/python scripts/verify_two_pair_sic_mixed_23_32_generic_local_gate.py
```

The checker records the full degree-seven bidegree ledger, verifies
(3.5)--(3.8), proves that the first four contractions equal the three
central polynomials in (4.3), checks the exact Gröbner bases underlying
(4.6), and audits the residual monomial cone.  The calculation is an
exact characteristic-zero certificate, not a bounded search.

The second checker verifies (7.2)--(7.4), constructs the complete central
moments through order five, checks all four power memberships in (7.9),
and audits the cutoff cone (7.10).  It explicitly records that the
\(V_{1,4}=0\) boundary remains the balanced \((3,3)\) problem rather than
promoting that open boundary to a theorem.

The third command constructs (9.5) directly, checks the exact
characteristic-zero Jacobian rank, and runs the two stated Singular
calculations.  Its artifact labels the result as finite-field evidence
and does not assert a rational radical.

The fourth command constructs the complete Clebsch--Gordan basis,
verifies the rational-normal and tangential ideals in (10.4) and (10.7),
checks the exact prime decomposition of the latter, and proves the three
two-sided eliminations and multiplier bounds in Section 10.  It also
constructs the full moment-one-through-four Jacobian at (10.17), verifies
the obstruction ideal (10.20), checks (10.21)--(10.22), and proves the
nonzero order-four resultant (10.24) after retaining all correction
parameters.

The fifth command proves the raising reduction (10.25)--(10.27), performs
the full calculation over \(\mathbf Q(u)\), factors the resultant
(10.28), proves the irreducibility assertions in (10.29), derives the
unique sextic-branch correction (10.32)--(10.34), and recomputes the five
displayed rational special points directly.  Its claim remains a local
first-four-moment jet theorem.
