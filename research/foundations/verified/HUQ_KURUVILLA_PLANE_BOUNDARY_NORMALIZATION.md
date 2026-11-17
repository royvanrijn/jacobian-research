# The characteristic-two plane counterexample: normalization and wild boundary

> **Status.** This is a repository theorem about the boundary of Mondello's
> external plane map.  The finite normalization and ideal decompositions have
> an exact Singular certificate, while a separate SymPy calculation replays
> the displayed source, chart, divisor, and local-different identities.  The
> SymPy replay is not a second implementation of the integral-closure
> algorithm.  An independent Sage/Singular, Magma, or differently organized
> CAS audit of the normalization, conductor, completed nodes, different, and
> residue ledger is still needed.  The result is a positive-characteristic
> boundary theorem and does not address the characteristic-zero plane
> Jacobian conjecture.

The plane map and its exact separable degree-three proof are in the
[parent audit](HUQ_KURUVILLA_CHARACTERISTIC_TWO_AUDIT.md#8-the-preserved-coordinate-and-the-dimension-two-theorem).
The formula, plane theorem, and preserved-fiber reduction are due to Romy
Mondello, [*A Dimension-Two Counterexample to the Separable Jacobian
Conjecture in Characteristic Two*,
arXiv:2608.02634v1](https://arxiv.org/abs/2608.02634) (2026),
building on Irit Huq-Kuruvilla's
[threefold map](https://arxiv.org/abs/2607.20968).  The normalization and
boundary statements below are the contribution of this audit.

Throughout, \(k\) is a field of characteristic two.  Geometric component and
monodromy statements are understood after extending \(k\) to an algebraic
closure.

## 1. The monic cubic order

Write the plane map as

\[
\begin{aligned}
 P&=x+x^2y+x^4+x^6y^2,\\
 Q&=y+x^5+x^6y+x^7y^2+x^8y^3.
\end{aligned}                                                    \tag{1.1}
\]

Put

\[
 r=1+xy,\qquad u=1+x^3r,\qquad T=ru^2.                         \tag{1.2}
\]

Then \(P=xru\), and the hidden-root calculation in the parent audit gives

\[
 H(T)=T^3+T^2+(PQ+P^3)T+P^3=0.                                \tag{1.3}
\]

The polynomial \(H\) is irreducible and separable over \(k(P,Q)\).  Define
the finite free cubic order

\[
 B_0=k[P,Q,T]/(H).                                               \tag{1.4}
\]

Its discriminant is

\[
 \boxed{\operatorname{disc}_T(H)=P^2Q^2=(PQ)^2.}                \tag{1.5}
\]

The Jacobian ideal of the hypersurface (1.4) has radical

\[
 \sqrt{(H,H_P,H_Q,H_T)}=(P,T).                                  \tag{1.6}
\]

Thus the primitive cubic order is nonnormal along one line.  This line is an
order-level defect; it will be distinct from the reconstruction boundary.

## 2. Exact finite normalization

In the common function field set

\[
 \boxed{Z=\frac{P^2}{T}.}                                       \tag{2.1}
\]

On the source this is polynomial:

\[
 Z=x^2r,\qquad ZT=P^2.                                          \tag{2.2}
\]

Let \(C\) be the quotient of \(k[P,Q,T,Z]\) by

\[
\begin{aligned}
 ZPT+ZP+PQ+T^2+T&=0,\\
 ZT+P^2&=0,\\
 ZP^2T+ZP^2+ZQT+PT^2+PT&=0,\\
 Z^2+ZP^2+ZQ+PT+P&=0,\\
 T^3+T^2+(PQ+P^3)T+P^3&=0.
\end{aligned}                                                    \tag{2.3}
\]

The exact normalization algorithm verifies that \(C\) is normal, that the
map \(B_0\to C\) is finite and birational, and that \(C\) is generated as a
\(B_0\)-module by

\[
 1,\quad \frac{P^2}{T}.                                         \tag{2.4}
\]

Consequently

\[
 \boxed{\overline X=\operatorname{Spec}C
 =\operatorname{Norm}_{\mathbb A^2_{P,Q}}k(x,y).}                \tag{2.5}
\]

The conductor of the primitive order inside its normalization is

\[
 \boxed{\mathfrak c_{B_0\subset C}=(P,T).}                      \tag{2.6}
\]

Upstairs its reduced support has two branches,

\[
 (P,T,Z),\qquad (P,T,Z+Q).                                      \tag{2.7}
\]

They are the closures of the two affine source components \(r=0\) and
\(u=0\).  In particular the conductor of the chosen primitive order is not
the missing-boundary divisor.

## 3. The reconstruction open

Under the source map (1.1)--(1.2), add \(Z=x^2r\).  Direct substitution
annihilates all five relations in (2.3).  More importantly,

\[
 \boxed{y=Q+Z^2(P+Z^2)}                                         \tag{3.1}
\]

already belongs to \(C\), while \(x\) has the three compatible rational
presentations

\[
 \boxed{
 x=\frac{T+1}{Q}
  =\frac{P}{T+PZ}
  =\frac{Z}{P+Z^2}.}                                            \tag{3.2}
\]

The common zero locus of the three denominators in (3.2) is the prime

\[
 \mathfrak e=(Q,P+Z^2,T+Z^3).                                  \tag{3.3}
\]

Indeed

\[
 C/\mathfrak e\simeq k[Z],\qquad
 P=Z^2,\quad Q=0,\quad T=Z^3.                                  \tag{3.4}
\]

Write \(E=V(\mathfrak e)\).  The principal opens of \(Q\), \(T+PZ\), and
\(P+Z^2\) cover \(\overline X\setminus E\).  Formula (3.2) defines one regular
function \(x\) on this cover, and (3.1) defines \(y\) globally.  Substitution
in both directions gives the identity on every chart.  Hence

\[
 \boxed{
 \mathbb A^2_{x,y}\simeq\overline X\setminus E.}                \tag{3.5}
\]

This is the exact Zariski--Main reconstruction open.  The unique missing
prime is an affine line, but its map to the reduced target discriminant is
not birational.

## 4. The two primes over \(Q=0\)

At \(Q=0\), the primitive cubic factors as

\[
 H(T)=(T+1)(T^2+P^3).                                           \tag{4.1}
\]

In the normalization the reduced pullback of \(Q=0\) has exactly two prime
components:

\[
\begin{array}{lll}
 A:& Q=0,\ T=1,\ Z=P^2,& A\simeq\mathbb A^1_P,\\[1mm]
 E:& Q=0,\ P=Z^2,\ T=Z^3,& E\simeq\mathbb A^1_Z.
\end{array}                                                      \tag{4.2}
\]

The component \(A\) is the retained ordinary sheet.  The component \(E\) is
the reconstruction boundary from (3.3).  Their intersection is

\[
 Z^3=1,                                                         \tag{4.3}
\]

so over an algebraic closure they meet in three reduced points.  These are
exactly the three punctures removed from the retained \(Q=0\) sheet.

The target map on the missing component is

\[
 \boxed{E\longrightarrow V(Q),\qquad Z\longmapsto P=Z^2.}       \tag{4.4}
\]

It is radicial of degree two.  Generically \(Q\) vanishes to order one on
both \(A\) and \(E\).  The complete generic degree ledger is therefore

\[
\begin{array}{c|c|c|c|c}
\text{prime}&e&f_{\rm sep}&f_{\rm insep}&ef_{\rm sep}f_{\rm insep}\\
\hline
A&1&1&1&1\\
E&1&1&2&2
\end{array}                                                      \tag{4.5}
\]

and the contributions sum to the geometric degree three.

### 4.1. Completed local model at the three intersections

The intersections in (4.3) are not an unresolved singular edge.  Near any
one of them, put \(p=P+Z^2\) and

\[
 a=\frac{Z(1+Z^3)+p^2}{Z^2}.
\]

Here \(Z\) is a unit.  Equation (5.1) is exactly

\[
 \boxed{Q=pa}.                                                  \tag{4.6}
\]

At \(p=0,\ Z^3=1\), one has

\[
 \frac{\partial a}{\partial Z}=Z^{-2},
\]

which is a unit.  Thus \(p,a\) are regular parameters on the normalized
surface.  In the completed local ring, the reduced inverse image of the
target divisor \(Q=0\) is

\[
 k[[p,a]]/(pa),                                                 \tag{4.7}
\]

with \(E=(p)\) and \(A=(a)\).  Each intersection is therefore an ordinary
node of the reduced boundary.  The conductor of the normalization of this
nodal divisor is its maximal ideal \((p,a)\).  This is a local statement
about the boundary divisor, distinct from the surface-order conductor
\((P,T)\) computed in Section 2.

## 5. Generic different and wild inertia

On the generic open of \(E\), \(Z(1+Z^3)\ne0\).  Put

\[
 p=P+Z^2.
\]

After eliminating \(T\), the normalized cubic relation becomes

\[
 \boxed{pZ(1+Z^3)+p^3+QZ^2=0.}                                 \tag{5.1}
\]

The coefficient of \(p\) is a unit, so \(p\) and \(Q\) both have order one
at \(E\).  Implicit differentiation gives

\[
 \frac{\partial P}{\partial Z}
 =\frac{p}{Z(1+Z^3)+p^2}.                                      \tag{5.2}
\]

The denominator is a unit and the numerator has order one.  Hence the
different exponent is

\[
 \boxed{v_E(\mathfrak D_{\overline X/\mathbb A^2})=1.}          \tag{5.3}
\]

Thus \(E\) is the smallest possible fierce boundary row:

\[
 (e,f_{\rm sep},f_{\rm insep},v_E(\mathfrak D),\text{sheet loss})
 =(1,1,2,1,2).                                                  \tag{5.4}
\]

At the three nodes, (4.6) also gives
\(dP\wedge dQ=p(\partial a/\partial Z)\,dp\wedge dZ\).
Consequently the different ideal remains \((p)\): its multiplicity along
\(E\) is still one through every intersection with \(A\).

The map is etale on the affine source (3.5); all wildness is carried by the
deleted divisor at infinity.

## 6. Monodromy

Irreducibility of \(H\) makes the geometric monodromy a transitive subgroup
of \(S_3\), hence either \(C_3\) or \(S_3\).  The extension is not Galois:
over the prime \(Q=0\), the two primes \(A\) and \(E\) have different residue
data, respectively separable degree one and inseparable degree two.  A cyclic
cubic Galois action would act transitively on primes above one base prime and
preserve their ramification and residue data.  Therefore

\[
 \boxed{G_{\rm geom}=S_3.}                                     \tag{6.1}
\]

The arithmetic group contains the geometric group and is itself contained
in \(S_3\), so

\[
 \boxed{G_{\rm arith}=S_3.}                                    \tag{6.2}
\]

In the Galois closure, the wild inertia detected by (4.4)--(5.3) acts as a
transposition.  The prime-to-characteristic generic degree three therefore
coexists with order-two wild inertia at the missing boundary.

This is one residue-ledger proof of the monodromy statement.  Independent
discriminant/resolvent and target-line specialization calculations, together
with a separate derivation of the wild inertia filtration, remain open
assurance tasks; the current checker does not supply those second routes.

## 7. Consequence for the positive-characteristic programme

This closes the first explicit plane row of the
[wild missing-boundary atlas](../extended-geometry/POSITIVE_CHARACTERISTIC_DEFORMATION_LANDSCAPE.md#pcd4--wild-missing-boundary-atlas):

1. the primitive cubic order, normalization, and conductor are explicit;
2. the reconstruction open deletes exactly one prime;
3. the reduced discriminant component has one retained sheet and one
   radicial missing component;
4. the full generic ledger and different exponent are known;
5. the three retained/missing intersections are completed ordinary nodes,
   with nodal conductor and different explicitly determined; and
6. geometric and arithmetic monodromy are both \(S_3\).

The theorem also isolates a warning for reduction-modulo-two arguments.  A
separable cover of odd generic degree can acquire wild residue degree at
infinity even when its affine Jacobian is a unit.  Any horizontal boundary
comparison must retain inseparable residue degree and the different, not only
ordinary ramification index and reduced discriminant support.

What remains open is the corresponding computation for systematic higher
wild families.  The separate direct lifting question for this exact plane map
is settled negatively by the
[modulo-four obstruction theorem](HUQ_KURUVILLA_PLANE_W2_OBSTRUCTION.md);
its de Rham obstruction is invariant under plane polynomial equivalence.
Extending the explicit one-variable stable Witt tower with uniformly bounded
polynomial degree remains open.  Neither issue affects the plane boundary
theorem proved above.

<!-- status-consumer: HKM2W1 904c57385ac0b0dd -->

## 8. Exact reproduction

Requires Singular:

```bash
.venv/bin/python scripts/verify_huq_kuruvilla_plane_boundary.py
```

The checker verifies the source presentation, discriminant, reconstruction
charts, both \(Q=0\) components, their reduced intersection, and the local
nodal and different equations.  It then asks Singular to certify the integral closure,
normality, primitive-order conductor, reconstruction-boundary ideal, and the
two decompositions above \(Q=0\) and above the conductor line.
