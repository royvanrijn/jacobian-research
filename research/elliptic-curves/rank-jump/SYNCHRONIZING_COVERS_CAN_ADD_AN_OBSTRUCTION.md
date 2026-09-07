# Synchronizing residual covers can add an arbitrary local obstruction

Combining the two retained residual quartics produces a genus-five curve,
with a genus-three quotient common to both signs. The first parameter
alignment makes the negative-sign combination insoluble over
\(\mathbb Q_2\), even though both individual covers are locally soluble.
Swapping the two parameter coordinates on the second cover removes that
dyadic obstruction without changing either underlying class.

Thus an arbitrarily chosen common auxiliary curve can impose an extra
solubility condition. Its obstruction must not be mistaken for an
intrinsic obstruction to the two classes being separately soluble.
The experiment also identifies precisely what the common genus-three
quotient forgets: an unramified quadratic lifting condition.

## Point-free parametrizations on the small control

The [initial protocol](RESIDUAL_DOUBLE_COVER_PROTOCOL.json) takes only
the two conic/quadric pairs from the existing
[rational/Sha control](NORM_LIFTS_CAN_BE_ENTIRELY_SHA.md).
PARI's quadratic-form solver supplies conic points
\[
(-14,-12,1),\qquad (2,1,0).
\]
Neither retained exceptional elliptic points nor supplied 2-cover points
enter that solver call or the parameter identification.

For a conic Gram matrix \(M\), an isotropic vector \(p\), and the plane
spanned by the two coordinate axes away from the first nonzero coordinate
of \(p\), put \(v=s e_i+t e_j\). The parametrization is
\[
z(s,t)=Q(v)p-2B(p,v)v.
\]
Its image satisfies \(Q(z)=0\). It has no projective base point: a common
zero would give an isotropic vector orthogonal to \(p\), independent of
\(p\), contradicting nonsingularity of a ternary conic. Thus it gives the
usual isomorphism from \(\mathbb P^1\) to the conic.

Substitution into the residual quadrics and removal of square factors
\(100\) and \(25\) gives
\[
\begin{aligned}
f(s,t)={}&1753s^4+89212s^3t+1666268s^2t^2\\
 &+13531912st^3+40363728t^4,\\
g(s,t)={}&16s^4-128s^3t-3016s^2t^2+3088st^3+82281t^4.
\end{aligned}
\]
The individual curves \(y^2=f\) and \(y^2=g\) are the rational covers
of \(E_{\rm rat}:y^2=x^3+11x^2-14x+1\).
The curves \(y^2=-f\) and \(y^2=-g\) are the locally soluble nonzero
Sha covers of \(E_{\rm Sha}:y^2=x^3-11x^2-14x-1\).
This sign convention is for the residual square equations; it reverses
the sign on the \(h^2\) term of the preceding quadric-pencil notation.

The parametrizations and quartics are frozen. They depend on the conic
solver's choice of point and coordinate axes, and are not intrinsic to
the Selmer classes.

## Exact geometry of the simultaneous condition

For either alignment \(\phi(s,t)=(s,t)\) or \((t,s)\), define the smooth
projective normalizations
\[
C_{\epsilon,\phi}:\quad
y_0^2=\epsilon f(s,t),\qquad
y_1^2=\epsilon g(\phi(s,t)),
\quad \epsilon=\pm1,
\]
and
\[
D_\phi:\quad w^2=f(s,t)g(\phi(s,t)).
\]
The diagonal involution \((y_0,y_1)\mapsto(-y_0,-y_1)\) has quotient
\(D_\phi\), via \(w=y_0y_1\). The quotient does not depend on \(\epsilon\).

Both quartics are separable, of degree four on the chosen affine chart,
and have nonzero resultant for both alignments. Their disjoint branch
sets give two independent geometric quadratic characters and eight
branch points. Riemann--Hurwitz therefore gives
\[
2g(C_{\epsilon,\phi})-2=-8+2\cdot8=8,\qquad
g(C_{\epsilon,\phi})=5,\qquad g(D_\phi)=3.
\]
The diagonal involution is unramified. At a root of one quartic it still
acts nontrivially on the other square root, and at infinity both even
degree covers are unramified. Equivalently the divisor of \(f\) on
\(D_\phi\) is even. The genera also satisfy the unramified degree-two
formula \(5=2\cdot3-1\).

On the open set where \(f(P)\ne0\), a rational point \(P\in D_\phi(\mathbb Q)\)
lifts to \(C_{\epsilon,\phi}(\mathbb Q)\) exactly when
\[
\boxed{[f(P)]=[\epsilon]\quad\text{in }\mathbb Q^\times/\mathbb Q^{\times2}.}
\]
Indeed one must choose \(y_0^2=\epsilon f(P)\); then \(y_1=w(P)/y_0\).
The product quotient has forgotten this squareclass. A rational point
on a common auxiliary quotient is therefore weaker than simultaneous
rational points on the original covers.

## A local obstruction caused entirely by alignment

For the identity alignment, reduction modulo eight gives
\[
f(s,t)\equiv s^2(s+2t)^2\pmod8,\qquad
g(s,t)\equiv t^4\pmod8.
\]
If \(s\) is odd, \(f\equiv1\pmod8\), so \(-f\) is a nonsquare unit.
If \(t\) is odd, \(g\equiv1\pmod8\), so \(-g\) is a nonsquare unit.
Every primitive dyadic parameter has at least one odd coordinate.
Consequently
\[
C_{-,{\rm id}}(\mathbb Q_2)=\varnothing.
\]
This is an exact projective obstruction, not a search miss.

The [follow-up protocol](RESIDUAL_ALIGNMENT_PROTOCOL.json) freezes just
one alternative: replace \(g(s,t)\) by \(g(t,s)\), using no point to
choose that swap. At \((s,t)=(2,1)\),
\[
-f(2,1)=-74834368,\qquad -g(1,2)=-1328896.
\]
Their dyadic valuations are 6 and 8; after removing those powers of two,
both odd units are \(1\bmod8\). Both are therefore squares in
\(\mathbb Q_2\), proving
\[
C_{-,{\rm swap}}(\mathbb Q_2)\ne\varnothing.
\]
These same two evaluations independently witness the separate dyadic
solubility of the negative quartics at parameters \((2,1)\) and \((1,2)\).
The identity alignment incorrectly demanded the same parameter for both.

For the positive sign, \((1,1)\) gives values \(55652873\) and \(82241\),
both \(1\bmod8\), so both alignments are dyadically soluble.

| Alignment | Positive sign over \(\mathbb Q_2\) | Negative sign over \(\mathbb Q_2\) | Negative sign over \(\mathbb Q\) |
|---|---|---|---|
| Identity | yes | no | no |
| Swap second coordinates | yes | yes | no |

The last column follows for **any** alignment: a rational point would
project to a rational point on each original Sha torsor, which the
retained exact proof excludes. We have not proved that the swapped
genus-five curve is everywhere locally soluble. Other local obstructions
remain possible, and it is not called a new Sha example.

The initial experiment also tested all 320 projective rational parameters
of height at most 16 on \(D_{\rm id}\), with no square product.
That is only a bounded miss. No rational parameter search was run for
the swapped quotient. Finite congruence survival at the other tested
prime powers is not promoted to full local solubility.

## The Jacobian rank changes without solving the curve

There is an additional exact consequence. The three nontrivial characters
of the biquadratic cover yield quotient genera \(1,1,3\), and an isogeny
over \(\mathbb Q\)
\[
\operatorname{Jac}(C_{\epsilon,\phi})
\sim E_\epsilon^2\times\operatorname{Jac}(D_\phi),
\]
where \(E_+=E_{\rm rat}\) and \(E_-=E_{\rm Sha}\).

This follows directly from differential pullbacks. On an affine chart
they span character spaces
\[
\frac{ds}{y_0},\quad \frac{ds}{y_1},\quad
\frac{ds}{y_0y_1},\quad\frac{s\,ds}{y_0y_1},\quad
\frac{s^2\,ds}{y_0y_1}.
\]
The five forms are regular and independent, with dimensions \(1,1,3\);
the corresponding sum of pullback homomorphisms of Jacobians has full
dimension and finite kernel.

The retained elliptic ranks are exactly 3 and 1. Writing
\(r_\phi=\operatorname{rank}\operatorname{Jac}(D_\phi)(\mathbb Q)\),
we obtain
\[
\operatorname{rank}\operatorname{Jac}(C_{+,\phi})=6+r_\phi,\qquad
\operatorname{rank}\operatorname{Jac}(C_{-,\phi})=2+r_\phi.
\]
Their rank difference is exactly four, although \(r_\phi\) remains
unknown and \(C_{-,\phi}\) has no rational points. These are consequences
of a proved isogeny and the existing elliptic ranks, not new Jacobian
descent calculations or evidence for rational solubility of the curves.

## Mechanism ranking and implications

1. **Valid solubility formulation:** a rational point on a common
   residual quotient plus the correct unramified quadratic lift
   supplies simultaneous points. The second implication remains
   essential; passing to the product quotient does not close it.
2. **Weak construction:** identifying two independently parametrized
   conics by equal parameter creates a genus-five auxiliary curve and
   may add local obstructions absent from the individual classes.
   One coordinate swap changes the outcome at 2. Such an obstruction
   is not intrinsic evidence about a high-gain fibre.
3. **Established but insufficient structure:** the Jacobian decomposition
   records two elliptic factors and gives an exact rank comparison.
   Neither those factors nor a large Jacobian rank supplies a rational
   point on the common curve.
4. **Missing mechanism:** an arithmetic reason, defined on the original
   family before exceptional points, that makes several residual covers
   soluble together. A useful common auxiliary curve needs an
   intrinsic construction or justified maps; arbitrary alignment and
   its local outcomes do not provide that reason.
5. **For Agent 1:** no new candidate feature follows. The quotient lifting
   condition is **solubility**; the synchronization obstruction is a
   construction-dependent local **solubility** obstruction, not an
   **incidence** predictor. No point-search **visibility** policy changes.

The next structural question is how many independent directions one
auxiliary rational point can actually supply through prescribed maps.
That bound should be established before choosing larger common curves
or interpreting their Jacobian factors as jump blocks.

## Reproducible evidence

- [Frozen conic parametrizations, quartics and initial checks](../../artifacts/generated-results/elliptic-curves/rank_jump_residual_double_covers_v1.json)
- [Two alignments and exact dyadic witnesses](../../artifacts/generated-results/elliptic-curves/rank_jump_residual_alignment_v1.json)
- [Independent polynomial and p-adic verification](../../artifacts/generated-results/elliptic-curves/rank_jump_residual_alignment_verification_v1.json)

The independent verifier reconstructs the quartics by rational polynomial
convolution, checks all primitive pairs modulo eight, and verifies the
dyadic witnesses in Sage's 2-adic field. It also checks separability,
resultants, and the dimensions used by the displayed genus and Jacobian
argument. It reuses the retained exact elliptic rank/Sha results.

    sage -python elliptic-curves/rank-jump/residual_double_covers.py check
    sage -python elliptic-curves/rank-jump/residual_alignment.py check
    sage -python elliptic-curves/rank-jump/verify_residual_alignment.py check

The original capture is checkpointed with a 30-second cap. The follow-up
is a fixed finite residue/witness calculation. No production fibre,
search policy, live output, or mathematical-status entry is changed.
