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
 \{P\in\mathbb Q[X]_{\leq5}:P\circ V=P\}
   =\mathbb Q[s]_{\leq5}.
 }
\tag{1.1}
\]

Consequently there is no invariant-slice quotient with a new polynomial
identity coordinate of degree at most five. In particular, no polynomial
change of coordinates can expose a second identity output by a linear or
degree-at-most-five invariant independent of \(s\), specialize it on the stored
collision, and thereby lower the canonical contraction from twenty pairs.
This excludes the whole degree-at-most-five invariant-slice class, not merely
a chosen ansatz.

The statement is deliberately narrower than nonlinear indecomposability.
It does not exclude a semiconjugate nonlinear quotient whose quotient
coordinates are not invariants, an invariant of degree at least six, or a
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

The exact audit
[`audit_bcw_21_sextic_defect_sectors.py`](../scripts/audit_bcw_21_sextic_defect_sectors.py)
computes only these two blocks of
\[
L_{\mathcal H}\colon\operatorname{Sym}^6\longrightarrow
\operatorname{Sym}^8.
\]
Modulo \(1000003\), the \((1,3)\) block has full column rank \(103\), the
\((2,2)\) block has full column rank \(1604\), and appending the corresponding
defect raises each rank. Hence neither defect lies in the sextic Lie image
over \(\mathbb Q\).

Thus no sextic correction can rescue a lower-degree \(Q\)-direction. Any
invariant of total degree at most six has lower-degree part in
\(\mathbb Q[s]\); the only remaining possibility is a genuinely new
homogeneous sextic invariant. The present computation does not classify all
220 sextic weight sectors, so it does not extend (1.1) from five to six.

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
\tag{1.8}
\]

Reconstructing the frozen 17-step trace gives

\[
 c_4=-x(v_3y+v_5z),\qquad
 Q=-v_3v_5-v_3xy-v_5xz.
\tag{1.9}
\]

The stable source section is defined by the exposed-factor equations

\[
 v_3=-xz,\qquad v_5=-xy.
\]

Consequently

\[
 Q|_{\mathrm{source\ section}}=x^2yz=M.
\tag{1.10}
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
\tag{1.11}
\]

and then specializing \(\phi\) at its common value on the stored collision.
Equation (1.11) makes that output an identity coordinate, so its dual variable
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
be pulled back as a new degree-at-most-five polynomial identity coordinate
of the stored BCW map. Any provenance-preserving Hopf realization must
therefore use degree at least six, a non-invariant
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
2. **Classify pure homogeneous sextic invariants.** The two correction
   sectors are closed, so lower-degree mixing is no longer the issue. The
   remaining question is the fixed kernel on homogeneous sextics. The direct
   Lie matrix has \(230230\) columns split into 220 bidegree blocks by (1.7).
   A blockwise kernel computation, preferably with a triangular circuit
   lemma rather than generic elimination, would decide whether a genuinely
   new sextic invariant exists.
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
./.venv/bin/python scripts/audit_keller_near_invariant_backtrace.py
python3 scripts/audit_keller_observable_quotients.py
python3 scripts/audit_keller_provenance_compression.py
```

The first command proves the degree-at-most-five nonlinear obstruction over
characteristic zero. The second excludes sextic corrections of both
lower-degree near-invariant directions. The third identifies the
near-invariant as a determinantal shared-factor gate residual. The fourth
proves the dimension-thirteen lower bound for semiconjugacies carrying
either distinguished observable. The fifth checks (4), the recurrence
closures, and the finite stored-circuit census.
