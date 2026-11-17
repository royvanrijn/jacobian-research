# Keller-provenance compression of the Image witness

## Result and scope

This note keeps the source of every candidate fixed: the foundational
three-point Keller collision and its stored BCW transforms.  It does not
search the absolute two-pair problem.

No witness below the current twenty-pair identity slice is obtained. There
is, however, an exact obstruction to a complete low-degree nonlinear quotient
class. Let

\[
 V(X)=X+\mathcal H(X)\colon\mathbb A^{21}\longrightarrow\mathbb A^{21}
\]

be the map in
[`essential_bcw_21_counterexample.json`](../artifacts/generated-results/essential_bcw_21_counterexample.json),
and write \(s=X_{20}\).  Then

\[
 \boxed{
 \{P\in\mathbb Q[X]_{\leq6}:P\circ V=P\}
   =\mathbb Q[s]_{\leq6}.
 }
\tag{1.1}
\]

Consequently there is no invariant-slice quotient with a new polynomial
identity coordinate of degree at most six. In particular, no polynomial
change of coordinates can expose a second identity output by a linear or
degree-at-most-six invariant independent of \(s\), specialize it on the stored
collision, and thereby lower the canonical contraction from twenty pairs.
This excludes the whole degree-at-most-six invariant-slice class, not merely
a chosen ansatz.

The statement is deliberately narrower than nonlinear indecomposability.
It does not exclude a semiconjugate nonlinear quotient whose quotient
coordinates are not invariants, an invariant of degree at least seven, or a
contraction identity that is compressed without exposing an identity output.

## 1. Exact low-degree invariant obstruction

For degrees at most two, let

\[
 T_{\leq2}\colon\mathbb Q[X]_{\leq2}\longrightarrow\mathbb Q[X],
 \qquad T_{\leq2}(P)=P(V(X))-P(X).
\tag{1.2}
\]

The domain has

\[
 1+21+\binom{22}{2}=253
\]

monomials.  Since \(\mathcal H_{20}=0\), the three polynomials
\(1,s,s^2\) lie in the kernel.

The independent sparse audit
[`audit_bcw_21_low_degree_invariants.py`](../scripts/audit_bcw_21_low_degree_invariants.py)
expands all 253 columns of (1.2) from the stored rational artifact. After
reducing coefficients modulo the good prime \(1000003\), it finds column
rank \(250\). Thus a \(250\)-minor is nonzero modulo that prime and the same
integer minor is nonzero over \(\mathbb Q\).  Hence

\[
 \operatorname{rank}_{\mathbb Q}T_{\leq2}\geq250.
\]

The three displayed kernel elements give the reverse inequality. Therefore
the rank is exactly \(250\) and its nullity is three.

This is a characteristic-zero certificate, not an inference that a modular
kernel basis lifts.

### Cubic through quintic continuation

Put

\[
 Q=X_{18}s-X_6X_8,\qquad M=X_0^2X_1X_2.
\tag{1.3}
\]

Exact pullback from the artifact gives the one-term defect

\[
 \boxed{Q(V)-Q=-Ms^2.}
\tag{1.4}
\]

Thus \(Q\) is a first integral of the homogeneous cubic vector field
\(\mathcal H\), but not an invariant of the discrete map \(V=I+\mathcal H\).
Let

\[
 L_{\mathcal H}(P)=\nabla P\mathbin\cdot\mathcal H.
\]

Good-prime ranks on homogeneous forms are

\[
\begin{array}{c|c|c|c}
d&\dim\operatorname{Sym}^d&\operatorname{rank}L_{\mathcal H}
 &\ker L_{\mathcal H}\\ \hline
2&231&229&\langle s^2,Q\rangle\\
3&1771&1769&\langle s^3,sQ\rangle\\
4&10626&10623&\langle s^4,s^2Q,Q^2\rangle\\
5&53130&53127&\langle s^5,s^3Q,sQ^2\rangle .
\end{array}
\tag{1.5}
\]

The displayed elements are exact kernel elements because
\(L_{\mathcal H}(s)=L_{\mathcal H}(Q)=0\). Their number equals the modular
nullity, so the nonzero modular minors prove that these are the complete
kernels over \(\mathbb Q\).

The only possible even-degree mixing through degree four would cancel the
degree-six defect of a quadratic \(Q\)-term by a quartic Lie derivative.
That requires

\[
 Ms^2\in
 \operatorname{im}\bigl(L_{\mathcal H}:\operatorname{Sym}^4
 \longrightarrow\operatorname{Sym}^6\bigr).
\]

Appending \(Ms^2\) to the matrix raises its rank from \(10623\) to \(10624\)
modulo \(1000003\), so this is impossible over \(\mathbb Q\). Likewise, the
only odd mixing through degree five would require \(Ms^3\) in the degree-five
Lie image; the rank rises from \(53127\) to \(53128\).

It remains to inspect the Lie-kernel directions themselves. From (1.4),

\[
\begin{aligned}
 (s^{d-2}Q)(V)-s^{d-2}Q&=-Ms^d,\\
 Q(V)^2-Q^2&=-2Ms^2Q+M^2s^4 .
\end{aligned}
\tag{1.6}
\]

In a linear combination of the two \(Q\)-directions, the \(M^2s^4\) term
first forces the \(Q^2\)-coefficient to vanish; the remaining \(Ms^d\) term
then forces the other coefficient to vanish. Thus these defects cannot cancel
inside either kernel list in (1.5).
Pullback preserves total-degree parity because \(V\) has only degrees one
and three. Separating its even and odd parts proves (1.1).

### The sextic correction layer

Degree six can interact with the degree-eight defects of the two quartic
Lie-kernel directions. Two independent torus gradings make that interaction
small. Their weights on \(X_0,\ldots,X_{20}\) are

\[
\begin{aligned}
\omega_1={}&(-1,1,2,1,1,0,1,2,0,-2,0,1,1,2,-1,1,2,1,1,0,0),\\
\omega_2={}&(-4,2,5,2,2,-1,2,5,-1,-7,-1,2,2,5,-6,0,3,0,0,-3,1).
\end{aligned}
\tag{1.7}
\]

Every monomial of \(\mathcal H_i\) has the same two weights as \(X_i\), so
\(L_{\mathcal H}\) preserves the bidegree. The defect of \(s^2Q\) is
\(-Ms^4\), in sector \((1,3)\); the degree-eight part of the \(Q^2\) defect
is \(-2Ms^2Q\), in sector \((2,2)\).

The audit
[`audit_bcw_21_sextic_defect_sectors.py`](../scripts/audit_bcw_21_sextic_defect_sectors.py)
first computes these two blocks of
\[
L_{\mathcal H}\colon\operatorname{Sym}^6\longrightarrow
\operatorname{Sym}^8.
\]
Modulo \(1000003\), the \((1,3)\) block has full column rank \(103\), the
\((2,2)\) block has full column rank \(1604\), and appending the corresponding
defect raises each rank. Hence neither defect lies in the sextic Lie image
over \(\mathbb Q\).

The same two gradings split the \(230230\) sextic monomials into \(220\)
sectors. The audit now classifies every sector with at most \(3000\) source
monomials: \(192\) sectors containing \(73548\) monomials. All of these
blocks are injective modulo \(1000003\), except

\[
\begin{array}{c|c|c|c}
\text{sector}&\text{columns}&\text{rank}&
 \ker L_{\mathcal H}\text{ over }\mathbb Q\\ \hline
(0,6)&1&0&\langle s^6\rangle\\
(1,5)&24&23&\langle s^4Q\rangle\\
(2,4)&486&485&\langle s^2Q^2\rangle\\
(3,3)&5657&5656&\langle Q^3\rangle .
\end{array}
\tag{1.8}
\]

For the \(28\) sectors larger than \(3000\) columns, recursive unique-row
peeling proves injectivity without arithmetic in \(25\) sectors. Only three
residual cores remain:

\[
\begin{array}{c|c|c|c}
\text{sector}&\text{original columns}&\text{peeled}&\text{core rank}\\ \hline
(2,0)&4670&3074&1596/1596\\
(3,3)&5657&3701&1955/1956\\
(4,6)&6101&3931&2170/2170 .
\end{array}
\]

Exact rational differentiation verifies the four displayed kernel
relations in (1.8). The modular ranks prove that those spans are the entire
rational kernels of their blocks; every other block is injective. Hence

\[
\ker\!\left(
L_{\mathcal H}\colon\operatorname{Sym}^6\to\operatorname{Sym}^8
\right)
=\langle s^6,s^4Q,s^2Q^2,Q^3\rangle
=\mathbb Q[s,Q]_6.
\]

The exact defect identity (1.4) also gives

\[
\begin{aligned}
(s^4Q)(V)-s^4Q&=-Ms^6,\\
(s^2Q^2)(V)-s^2Q^2&=-2Ms^4Q+M^2s^6.
\end{aligned}
\]

The last kernel direction is ruled out as a fixed polynomial by

\[
Q(V)^3-Q^3=-3Ms^2Q^2+3M^2s^4Q-M^3s^6\ne0.
\]

Thus no sextic correction can rescue a lower-degree \(Q\)-direction. Any
invariant of total degree at most six has lower-degree part in
\(\mathbb Q[s]\), and the complete sextic kernel calculation leaves only
\(s^6\) fixed. This proves (1.1).

### The degree-seven correction shift

The two quintic Lie-kernel directions not in \(\mathbb Q[s]\) have defects

\[
\begin{aligned}
(s^3Q)(V)-s^3Q&=-Ms^5,\\
(sQ^2)(V)-sQ^2&=-2Ms^3Q+M^2s^5.
\end{aligned}
\]

Their degree-nine terms could only be cancelled by
\(L_{\mathcal H}\) applied to a homogeneous septic. The required sectors are
\((1,4)\) and \((2,3)\). Exact monomial enumeration gives the stronger
structural identities

\[
\operatorname{Sym}^7_{(1,4)}
 =s\,\operatorname{Sym}^6_{(1,3)},\qquad
\operatorname{Sym}^7_{(2,3)}
 =s\,\operatorname{Sym}^6_{(2,2)}.
\]

Indeed every monomial in either septic sector is divisible by \(s\), and
division by \(s\) is a bijection onto the indicated sextic sector. Since
\(\mathcal H_{20}=0\),

\[
L_{\mathcal H}(sP)=sL_{\mathcal H}(P).
\]

The two sextic non-image certificates therefore transfer directly:
\(Ms^5\) and \(Ms^3Q\) are outside the respective septic Lie images. No
degree-seven term can rescue either quintic \(Q\)-direction. Extending (1.1)
through degree seven now depends only on classifying the pure homogeneous
degree-seven Lie kernel.

There is a useful \(s\)-adic reduction of that kernel problem. Write a
homogeneous septic as \(P=P_0+sR\), with \(P_0\) independent of \(s\).
Reducing \(L_{\mathcal H}(P)=0\) modulo \(s\) first requires

\[
\overline L_{\mathcal H}(P_0)=0
\quad\text{in}\quad
\mathbb Q[X_0,\ldots,X_{19}].
\]

The quotient source has \(657800\) monomials in \(204\) bidegree sectors.
The exact support audit recursively peels \(451891\) columns by unique rows;
\(205909\) columns in \(79\) sectors remain. This is not an injectivity
certificate: quotient-kernel classes really occur. For example,

\[
L_{\mathcal H}(X_9^7)=7sX_0^2X_9^6,
\]

so \(X_9^7\) is killed after reduction modulo \(s\), although it is not a
Lie invariant upstairs.

The reduced derivation itself has the exact vertical form

\[
\overline L_{\mathcal H}
=\sum_{i=14}^{19}A_i(X_0,\ldots,X_{13})\,
  \frac{\partial}{\partial X_i}.
\]

Thus \(B=\mathbb Q[X_0,\ldots,X_{13}]\) is pointwise fixed modulo \(s\).
Let \(I=(A_{14},\ldots,A_{19})\subset B\). For a base septic \(P_0\in B_7\),
the base-only part of the first lifting equation defines the necessary
obstruction

\[
\delta(P_0)=
\left[s^{-1}L_{\mathcal H}(P_0)\right]_{s=0}
\quad\text{in}\quad (B/I)_8.
\]

If \(\delta(P_0)\ne0\), no correction \(P_0+sR\) can be a Lie invariant.
For \(P_0=X_9^7\), the target is \(7X_0^2X_9^6\). Its bidegree would require
a sextic correction in sector \((-14,-50)\), but that source sector is
empty. Hence \(\delta(X_9^7)\ne0\), and the entire one-dimensional septic
sector \((-14,-49)\) is injective upstairs.

The coefficient ideal \(I\) is far from a complete intersection. The exact
Singular audit
[`audit_bcw_21_vertical_ideal.sing`](../scripts/audit_bcw_21_vertical_ideal.sing)
gives

\[
\dim(B/I)=12,\qquad \operatorname{ht}(I)=2,
\]

and the five minimal primes

\[
\begin{aligned}
&(X_1,X_2,X_6),\qquad (X_0,X_5,X_8),\qquad (X_0,X_3,X_8),\\
&(X_0,X_1),\\
&(X_0,\ X_5+3X_8,\ 7X_1-6X_3-2X_4-3X_6).
\end{aligned}
\]

Their intersection is exactly \(\sqrt I\). In particular, the dominant
degeneracy component is the twelve-dimensional plane \(X_0=X_1=0\), with
four additional codimension-three components. The degree-eight piece of
the obstruction quotient has

\[
\dim_{\mathbb Q}(B/I)_8=158412.
\]

Thus a regular-sequence or pure Koszul description of the vertical
syzygies is obstructed: the lifting calculation must separate the five
components and the nonradical torsion of \(I\). Componentwise saturation,
rather than a single generic vertical coordinate change, is the natural
next simplification.

For support-one base septics, the componentwise screen can be completed
without saturation. The exact good-prime audit
[`audit_bcw_21_septic_component_screen.py`](../scripts/audit_bcw_21_septic_component_screen.py)
evaluates the first obstruction of every one of the

\[
\dim_{\mathbb Q}B_7=\binom{20}{7}=77520
\]

base monomials on all five minimal components. The componentwise nonzero
counts are

\[
32032,\quad25740,\quad65520,\quad25740,\quad45136.
\]

Their union contains \(71588\) monomials. Each therefore has nonzero first
obstruction modulo \(I\), so no polynomial of the form

\[
P_0+sR,\qquad P_0\text{ one of those }71588\text{ monomials},
\]

can lie in the septic Lie kernel. Only \(5932\) base monomials survive this
radical screen. Exactly eight have identically zero first obstruction:

\[
X_3^{7-j}X_5^j,\qquad 0\leq j\leq7.
\]

The same audit removes the support-one restriction. The base septics occupy
only \(29\) bidegree sectors. After stacking all five component restrictions,
unique-row peeling certifies \(28764\) independent columns. Exact elimination
on the \(23\) residual cores contributes another \(32296\) pivots, giving

\[
\operatorname{rank}_{\mathbb F_{1000003}}
\left(
B_7\longrightarrow\bigoplus_{j=1}^5(B/P_j)_8
\right)=61060.
\]

The nonzero modular minor proves the same rank lower bound over
\(\mathbb Q\). Therefore the space of arbitrary base septics whose first
obstruction vanishes on every minimal component has dimension at most

\[
77520-61060=16460.
\]

This is now support-free at the radical level. It does not decide which of
those at most \(16460\) directions vanish modulo the nonradical ideal \(I\),
nor which then survive higher \(s\)-adic equations.

The remaining task has two layers: classify the residual quotient kernel,
then compute this obstruction and its higher \(s\)-adic successors on the
surviving classes. This is a sharper frontier than eliminating all
\(888030\) septic columns simultaneously.

### Circuit backtrace of the near-invariant

The sparse formula (1.4) has an exact circuit explanation. In the
rank-compressed 24-dimensional homogenization, quotient coordinates
\((X_6,X_8,X_{18},X_{20})\) lift to

\[
 (v_3,v_5,c_4,s),
\]

where \(c_4\) is the fifth scalar cubic gate. Therefore

\[
 Q=c_4s-v_3v_5.
\tag{1.9}
\]

Reconstructing the frozen 17-step trace gives

\[
 c_4=-x(v_3y+v_5z),\qquad
 Q=-v_3v_5-v_3xy-v_5xz.
\tag{1.10}
\]

The stable source section is defined by the exposed-factor equations

\[
 v_3=-xz,\qquad v_5=-xy.
\]

Consequently

\[
 Q|_{\mathrm{source\ section}}=x^2yz=M.
\tag{1.11}
\]

Together with \(Q(V)-Q=-Ms^2\), this shows that on \(s=1\) the map sends the
stable source section into the determinantal gate locus

\[
 c_4s-v_3v_5=0.
\]

Thus \(Q\) is a shared-factor **gate residual**, not a newly discovered
invariant of the normalized boundary model. The exact reconstruction is
checked by
[`audit_keller_near_invariant_backtrace.py`](../scripts/audit_keller_near_invariant_backtrace.py).

### Quotient interpretation

Call an **invariant-slice quotient** one obtained by finding a polynomial
coordinate \(\phi\) with

\[
 \phi\circ V=\phi
\tag{1.12}
\]

and then specializing \(\phi\) at its common value on the stored collision.
Equation (1.12) makes that output an identity coordinate, so its dual variable
does not occur in the canonical contraction and the corresponding source
variable can be treated as a coefficient exactly as \(s\) is in the known
twenty-pair slice.

If \(\deg\phi\leq5\), equation (1.1) gives \(\phi\in\mathbb Q[s]\). It supplies
no invariant algebraically independent of \(s\). If it is a polynomial
coordinate at all, its linear part in \(s\) merely recovers the existing
identity direction.  Thus the class contains no new collision-preserving
pair deletion.

## 2. Observable semiconjugacy lower bound

The invariant-slice class is not the only nonlinear quotient architecture.
Let \(W\) be the twenty-variable map obtained from \(V\) by setting \(s=1\).
A rational semiconjugacy consists of rational maps

\[
 \Phi\colon\mathbb A^{20}\dashrightarrow Y,\qquad
 f\colon Y\dashrightarrow Y,\qquad
 \Phi\circ W=f\circ\Phi.
\tag{2.1}
\]

An observable \(h\) is **carried** by the quotient if
\(h=\psi\circ\Phi\) for some \(\psi\in\mathbb Q(Y)\). Two natural carried
observables are relevant to the stored collision:

\[
 h_0=X_0,\qquad h_Q=X_{18}-X_6X_8.
\tag{2.2}
\]

Their values on the three source points are respectively

\[
 (0,1,-1),\qquad (0,-39/16,39/16),
\tag{2.3}
\]

so each separates the collision before \(W\); each takes one common value
after \(W\).

The audit
[`audit_keller_observable_quotients.py`](../scripts/audit_keller_observable_quotients.py)
computes the gradients of

\[
 h,h\circ W,\ldots,h\circ W^{12}
\]

at the integral point \(X_i=i^2+3i+5\), reduced modulo \(1000003\). For both
\(h=h_0\) and \(h=h_Q\), the resulting \(13\)-by-\(20\) matrix has rank
thirteen. A nonzero modular minor is a nonzero rational Jacobian minor, so
the thirteen iterates are algebraically independent over \(\mathbb Q\).

If (2.1) carries \(h\), then every iterate \(h\circ W^k\) belongs to the
pullback of \(\mathbb Q(Y)\). Therefore

\[
 \boxed{\dim Y\geq13.}
\tag{2.4}
\]

This excludes every rational semiconjugate quotient of dimension below
thirteen that retains either the established multiplier \(X_0\) or the
quadratic collision observable \(h_Q\). It does not exclude dimensions
thirteen through nineteen or a quotient using a different collision
observable.

At six tested modular points, both observable-Jacobian ranks remain thirteen
through the first twenty-five iterates. This is not a generic upper bound and
is recorded only as an experiment. It nevertheless defines a concrete
continuation: seek an explicit rational recurrence expressing the thirteenth
iterate in the field of the first thirteen, then prove pullback closure and
regularity. A successful closure would give a provenance-preserving rational
\(20\)-to-\(13\) semiconjugacy, but not automatically a thirteen-pair
contraction witness; contraction transport would still need proof.

There is no missed constant linear explanation for the plateau. Stacking the
first twenty-five iterate codistributions at three independent deterministic
points gives rank twenty for both observables. A constant translation
direction annihilating every iterate would lie in the kernel of this stacked
matrix, so none exists over \(\mathbb Q\). Any dimension-thirteen closure must
therefore use a genuinely nonlinear foliation or algebraic relation.

## 3. Why the normalized three-variable map is not already a witness

The determinant-one identity-linear normalization of the foundational map is

\[
 K=(F_3/2,F_2,F_1)=I+N.
\]

It is tempting to contract \(-w\mathbin\cdot N\) in three pairs.  This does
not satisfy the pure-power hypothesis.  Direct differentiation gives

\[
\begin{aligned}
 \operatorname{div}N={}&x^3y^3+6x^3yz+30x^2y^2\\
 &+\frac92x^2z+24xy\ne0.
\end{aligned}
\tag{4}
\]

Therefore the first contraction is already nonzero:

\[
 \mathcal E_3(-w\mathbin\cdot N)=-\operatorname{div}N.
\]

The Keller identity \(\det(I+JN)=1\) at scalar coefficient one does not imply
\(\det(I+tJN)=1\) for every \(t\).  Homogeneity of the later BCW correction
is what separates the Abhyankar--Gurjar inverse formula by contraction
power.  The mixed-degree three-variable normalization lacks that separation.

The exact sparse formula (4) is replayed by
[`audit_keller_provenance_compression.py`](../scripts/audit_keller_provenance_compression.py).

## 4. Pre-contraction circuit census

Different stored degree-lowering circuits were compared by the same rule:
specialize every zero correction component whose coordinate is constant on
the stored collision, then count the remaining active outputs.  The exact
artifact census gives

| stored circuit | active outputs after identity slices |
|---|---:|
| essential 21D route | 20 |
| sparse-conjugate essential route | 20 |
| constant-kernel 22D route | 21 |
| index-reduced 22D route | 21 |
| two rank-compressed 24D routes | 23 |
| shared 33D route | 23 |
| Long 79D route | 51 |
| original 95D homogeneous route | 60 |

This is a finite census of stored circuits, not an obstruction to all
degree-lowering circuits and not a minimality proof.  It shows that none of
the already certified alternative circuits contains a missed sub-twenty
identity slice.

## 5. Inverse-recurrence compression

For the cubic map, the inverse recurrence has the coordinate-dependency rule

\[
 G_i=Y_i-\mathcal H_i(G).
\tag{5}
\]

Starting at the collision-separated coordinate \(i=0\), close the set of
indices under the variables occurring in each required component
\(\mathcal H_i\).  The closure is all 21 coordinates.  After treating the
known identity coordinate \(s=X_{20}\) as a coefficient and specializing
\(s=1\), the closure is still every index \(0,\ldots,19\).

Thus literal subsystem extraction from (5) cannot represent the distinguished
inverse recurrence with fewer than twenty active coordinate pairs.  This is
an exact combinatorial obstruction to coordinate-deletion recurrence
compression.  It does not exclude nonlinear recombination of recurrence
states or a different multiplier.

The same audit script checks both closures and labels the stored-circuit
census as computation only.

## 6. Comparison with the Hopf mechanism

The sharp absolute two-pair witness uses a different cancellation geometry.
After radial separation, its Hopf variables satisfy

\[
 t^2+2xy=1,
\]

phase averaging extracts the weight-zero term in \(x\), and the pure moment
vanishes because

\[
 J_m'(X)=(1-X^2)^m
\]

has an order-\(m\) zero at the endpoint \(X=1\).  The multiplier shifts the
phase coefficient by one and leaves the nonzero beta integral.

The normalized boundary model of the foundational collision instead comes
from a linear-times-quadratic factorization with

\[
 \operatorname{Res}(L,Q)=1,\qquad [LQ]_{T^2S}=1,
\]

and residual torus weights

\[
 (a,y,z)\longmapsto(\lambda a,\lambda^{-1}y,\lambda^{-2}z).
\]

Its two noncentral collision points

\[
 (1,-3/2,-13),\qquad(-1,3/2,-13)
\]

are exchanged by \(\lambda=-1\), while \((0,0,1/2)\) is fixed.  This gives a
real structural resemblance to Hopf phase extraction: both mechanisms have
opposite-weight coordinates, a weight-zero relation, and a distinguished
order-two phase.

The resemblance does not provide a provenance map. The Hopf proof uses the
compact angular quotient, the interval coordinate \(t\), and endpoint
vanishing; the Keller boundary uses an algebraic torus, a resultant-one
factorization slice, and a three-point fiber. The low-degree invariant
obstruction (1.1) makes the gap precise: the Hopf angular coordinate cannot
be pulled back as a new degree-at-most-six polynomial identity coordinate
of the stored BCW map. Any provenance-preserving Hopf realization must
therefore use degree at least seven, a non-invariant
semiconjugacy, or a contraction-level construction not induced by an
identity slice.

The circuit backtrace sharpens this comparison. The quadratic relation that
looked Hopf-like is the \(2\)-by-\(2\) determinant \(c_4s-v_3v_5\) introduced
by shared-factor cubicization. Its rank-one locus contains the image of the
stable source section, whereas Hopf phase extraction integrates over a
compact angular quotient and uses endpoint multiplicity. The common
rank-one algebra is real, but the cancellation mechanisms remain different.

## 7. Remaining continuations

The audits leave five concrete routes, in decreasing order of immediacy.

1. **Close or refute observable dimension thirteen.** Seek a rational
   identity
   \[
   h\circ W^{13}=R(h,h\circ W,\ldots,h\circ W^{12})
   \]
   for \(h=X_0\) or \(h=X_{18}-X_6X_8\). One identity is not enough:
   denominators, dominance, pullback closure, and collision separation must
   be checked. If closure exists, determine whether Zhao contraction
   descends; semiconjugacy alone does not delete contraction pairs.
2. **Advance the invariant-slice obstruction to degree seven.** The sextic
   kernel and both lower-degree correction channels are now closed. What
   remains is only the pure degree-seven Lie kernel. The mod-\(s\) audit has
   already peeled \(451891\) of \(657800\) quotient columns. Classify the
   \(205909\)-column residual quotient kernel sectorwise, then compute the
   first obstruction in \(B/(A_{14},\ldots,A_{19})\) and its higher
   \(s\)-adic successors. The vertical ideal has five minimal components;
   the stacked restriction has rank at least \(61060\), leaving a
   radical-level base-septic subspace of dimension at most \(16460\).
   Control the component intersections and embedded torsion next. The expected upstairs kernel is
   \(\mathbb Q[s,Q]_7=\langle s^7,s^5Q,s^3Q^2,sQ^3\rangle\).
3. **Change the collision observable.** The bound (2.4) applies only to
   quotients carrying \(X_0\) or \(h_Q\). Search low-degree observables that
   still separate the three source points but whose iterate field has smaller
   transcendence degree. Every proposed upper bound needs an explicit
   recurrence; a rank plateau at sampled points remains only evidence.
4. **Compress at the contraction level.** A smaller-pair representation
   need not come from a quotient of \(V\). The right invariant is whether
   \(p=w\mathbin\cdot H\), its multiplier, and all contractions factor through
   a smaller Weyl/contraction algebra. Candidate obstructions are the
   characteristic support of the cyclic contraction module, a minimal
   realization rank for the inverse-coordinate recurrence, or a Poisson/Weyl
   rank invariant preserved by contraction-compatible changes.
5. **Exploit the gate-residual filtration.** The backtrace is now complete:
   \(Q=c_4s-v_3v_5\) is a determinantal shared-factor residual and restricts
   to \(M=x^2yz\) on the stable source section. The next question is whether
   iterated residuals close to the observed dimension-thirteen field, or
   whether a Weyl/contraction module built from this rank-one relation gives
   a lower bound on pair count.

There is also a proof-simplification target. The modular ranks in (1.5) are
exact characteristic-zero certificates, but a circuit-level triangular
description of the low-degree Lie kernels would replace four large sparse
rank calculations by a conceptual lemma and might expose the degree-six
pattern.

## Reproduction

Run

```bash
python3 scripts/audit_bcw_21_low_degree_invariants.py
python3 scripts/audit_bcw_21_sextic_defect_sectors.py
Singular -q scripts/audit_bcw_21_vertical_ideal.sing
python3 scripts/audit_bcw_21_septic_component_screen.py
./.venv/bin/python scripts/audit_keller_near_invariant_backtrace.py
python3 scripts/audit_keller_observable_quotients.py
python3 scripts/audit_keller_provenance_compression.py
```

The first command proves the obstruction through degree five over
characteristic zero. The second excludes sextic corrections of both
lower-degree near-invariant directions, classifies all \(220\) sextic Lie
sectors, extends the fixed-space theorem through degree six, and excludes
degree-seven corrections of the quintic near-invariants. It also performs
the exact mod-\(s\) support reduction of the pure septic frontier. The
Singular command computes the five minimal components and Hilbert function
of the vertical obstruction ideal. The component-screen command excludes
\(71588\) of \(77520\) support-one base septics and proves a support-free
stacked-component rank of \(61060\). The backtrace command identifies the
near-invariant as a determinantal shared-factor gate residual. The observable
command proves the dimension-thirteen lower bound for semiconjugacies
carrying either distinguished observable. The final command checks (4), the
recurrence closures, and the finite stored-circuit census.
