# Rigidity within the current cancellation ansatz

This is a scoped appendix to the cancellation construction.  Its results answer
one question: how far can the present triangular, one-inverse-variable
mechanism be relaxed before a new polynomial Keller branch appears?

## 1. Universal ansatz

Let

\[
 A=1+xf(y),\qquad B=A^bz+g(y,A),
\]

\[
 P=A^aB,\qquad Q=y+xA^cB,
\]

where `f` is nonconstant, `g` is arbitrary, `a,b>=1`, and `c>=0`.  Put
`s=xA^{c-a}` so that `Q=y+sP`, and allow an arbitrary nonzero one-factor
polynomial derivative

\[
 R_s=C\Theta\bigl(s f(Q-Ps)\bigr).
\]

This class quantifies over all integer weights, all `f,g`, and all such
polynomial derivatives.  It is maximal relative to these structural axioms;
it is not a classification of all Keller maps.

The Jacobian is

\[
 -CA^{b+c-1}\bigl((c-a+1)A-(c-a)\bigr)
 \Theta\bigl((A-1)A^{c-a}\bigr).                           \tag{1}
\]

Constancy forces

\[
 c=a-1,\qquad \Theta(W)=\lambda(1-W)^{a+b-2}.              \tag{2}
\]

Polynomiality then forces `a=1`, a translated/scaled pure power for `f`, and
the unique cancellation jet from [CONSTRUCTION.md](CONSTRUCTION.md).  Any
allowed `A^(r+1)` tail is removed by a polynomial source shift.  Thus every
polynomial Keller member of this explicit class is left--right equivalent to
a cancellation map.

## 2. Relaxations already exhausted

The same conclusion survives three natural enlargements.

### Arbitrary target-dependent derivative

After (2), allow `R_s=H(s,P,Q)` for an arbitrary polynomial `H`.  Algebraic
independence of `(s,P,Q)` and the constant-Jacobian identity force

\[
 H(T,P,Q)=\lambda\{1-Tf(Q-PT)\}^{a+b-2}.
\]

Target dependence therefore supplies no new branch inside this skeleton.

### Finitely many normalized resolvent factors

Replace the derivative by a finite product of factors
`1-Tf_i(Q-PT)`.  The zero divisors required by the Jacobian identity force
every nontrivial normalized factor to coincide with the original factor;
the rest are identically one.  All factors coalesce, and polynomiality again
returns the cancellation family.

### General two-weight cancellation input

In the resulting two-weight skeleton, the leading functional equation,
spectral coprimality, and recursive jet equations exclude nonmonomial leading
solutions.  This proves classification within the stated ansatz, not beyond
the reconstruction skeleton.

## 3. Cover and parameter rigidity

Let `Xbar->Y` be the finite normalization of a connected generically
separable cover, with Galois closure group `G` and point stabilizer `H`.  Then

\[
 \operatorname{Aut}_Y(Xbar)\simeq N_G(H)/H,                \tag{3}
\]

equivalently the centralizer of the geometric monodromy action.  For the
natural `S_N` action (`N>=3`) and the natural `A_N` action (`N>=4`), `H` is
self-normalizing, so this automorphism group is trivial.  The low-degree
exclusions are real and should not be suppressed.

For cancellation branches, the unique target-fixed birational
identification between normalized jets `h` and `h'` is

\[
 (x,y,z)\longmapsto
 \left(x,y,z+\frac{y^{m+1}(h(A)-h'(A))}{A^{r+1}}\right).
\]

It is polynomial exactly when `h=h' mod A^(r+1)`.  Distinct normalized
parameter roots are therefore not target-fixed right-equivalent, even after
stabilization.

For unrestricted left--right equivalence, preservation of the labelled
boundary pair, including after stabilization, restricts a target automorphism
on `P=0` to

\[
 (Q,R)\longmapsto(uQ,u^{-m}R).
\]

The marked affine reconstruction open now eliminates the possible parameter
motion from the residual congruence kernel.  In the source UFD, the two affine
components of `P=AB=0` have residue degrees `1` and `r+1`, so an equivalence
preserves `A` and `B` separately.  The unique reducible fiber of
`A=1+xy^m` fixes the scalar of `A`; factorization of `A-1=xy^m` then recovers
`y` up to a scalar `u`.  The global identity `Q=y+xB`, together with the
boundary-plane restriction, forces

\[
 A\mapsto A,\qquad y\mapsto uy,\qquad B\mapsto u^{m+1}B.
\]

Consequently the **cancellation reconstruction residue**

\[
 q=\left.\frac{B}{y^{m+1}}\right|_{A=0}
\]

is fixed.  Thus distinct normalized parameter roots are not unrestricted
stably left--right equivalent.  This proves only that every lift preserving
the marked open acts trivially on the parameter roots; the full residual
target congruence kernel may still be nontrivial.  The canonical stable proof
is Theorem 6.2 in
[`cancellation-parameter-faithfulness.tex`](../papers/marked-root-multiplicity/cancellation-parameter-faithfulness.tex).

## 4. Scope and failure modes

| Assumption | What it supplies | What may happen if removed |
|---|---|---|
| One rank-one weight `A=1+xf(y)` | the valuation and finite cancellation operator | several independent valuations and coupled jets |
| `B` triangular and affine-linear in `z` | direct Jacobian factorization and reconstruction | nonlinear source fibers or another inverse variable |
| Monomial leading weights | the exponent identity (1) | nonmonomial units may cancel valuation factors differently |
| `Q=y+sP` with one inverse variable | algebraic independence of `(s,P,Q)` | multi-primitive-element reconstruction |
| Polynomial one-factor or finite-factor derivative | zero-divisor coalescence | rational, algebraic, or mixed factors need not coalesce |
| Polynomial coordinate equivalence | divisibility and tail removal | birational equivalence is much weaker |

Accordingly, the theorem should be cited as “rigidity within the current
ansatz,” never as a universal classification of cancellation mechanisms.
The next honest enlargement changes the skeleton itself by adding a source
function, source variable, or independent inverse variable.

Detailed proofs and exact regressions remain in
[the cancellation archive](../archive/cancellation-components/); their
scripts remain active under `scripts/`.
