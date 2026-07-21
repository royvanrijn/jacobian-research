# Universal three-weight classification in the monomial triangular skeleton

This note classifies a precisely defined universal ansatz rather than one
selected deformation.  It quantifies over every integer triple of monomial
weights, every one-variable input, every triangular correction, and every
polynomial derivative in the single natural inverse factor.  The
constant-Jacobian equation forces the third weight back to the two-weight
slice; the stronger target-dependent derivative theorem and the complete
polynomial cancellation theorem then force C24.

Throughout, `k` is a characteristic-zero field, `f in k[y]` is nonconstant,
`a,b>=1`, and `c>=0` are integers.  Let

\[
 A=1+xf(y),\qquad B=A^bz+g(y,A),                              \tag{1}
\]

\[
 P=A^aB,\qquad Q=y+xA^cB.                                    \tag{2}
\]

The former two-weight skeleton imposed `c=a-1`.  Here `c` is independent.
Put

\[
 h=c-a,\qquad s=xA^h,\qquad W(t)=t f(Q-Pt).                   \tag{3}
\]

Then `Q=y+sP`.  For a nonzero polynomial `Theta in k[W]`, consider the most
general one-factor polynomial resolvent derivative in this inverse variable,

\[
 R=C\int_0^s\Theta\bigl(t f(Q-Pt)\bigr)dt,
 \qquad C\in k^*.                                             \tag{4}
\]

All expressions are initially interpreted in `k[x,y,z,A^(-1)]`.  The
question is which weights and which `Theta` can make the localized Jacobian
a nonzero constant.  Polynomial cancellation is a subsequent, stronger
condition.

## 0. The universal ansatz and maximality statement

### Definition 0.1 (normalized monomial triangular one-factor class)

Let `mathscr T_3` be the class of all data

\[
 (f,g,a,b,c,\Theta,C)
\]

satisfying:

1. `f in k[y]` is nonconstant and `g in k[y,A]` is arbitrary;
2. `A=1+xf(y)` is the single rank-one source weight;
3. `B=A^b z+g(y,A)` is triangular and affine-linear in `z`, with pure
   monomial leading coefficient;
4. `P=A^aB` and `Q=y+xA^cB`, where all integer monomial weights
   `a,b>=1`, `c>=0` are allowed independently;
5. `s=xA^(c-a)` is the resulting inverse variable, uniquely characterized
   by `Q=y+sP` inside this monomial skeleton; and
6. `R_s=C Theta(t f(Q-Pt))|_(t=s)` for an arbitrary nonzero
   `Theta in k[W]` and `C in k^*`, with the lower integration endpoint fixing
   only an additive target normalization.

Thus `mathscr T_3` is universal **relative to these six structural axioms**:
no exponent, coefficient of `g`, polynomial `f`, or one-factor derivative is
held generic or selected from a finite list.  Nonzero scalar coefficients in
the monomials are omitted only because source and target scalings normalize
them.

### Theorem 0.2 (maximality inside `mathscr T_3`)

Every member of `mathscr T_3` whose three coordinates are polynomial and
whose Jacobian is a nonzero constant is, up to translations, nonzero
scalings, and polynomial left--right equivalence, a C24 map.

More precisely:

1. localized Keller forces

   \[
    c=a-1,\qquad
    \Theta(W)=\lambda(1-W)^{a+b-2};                         \tag{A}
   \]
2. after `c=a-1` is forced, allowing an arbitrary polynomial derivative
   `H(T,P,Q)`, rather than a one-factor `Theta`, still forces

   \[
    H=\lambda\{1-Tf(Q-PT)\}^{a+b-2};                       \tag{B}
   \]
3. polynomiality forces `a=1`, makes `f` a translated/scaled pure power,
   determines the complete cancellation jet, and removes its
   `A^(a+b-1)` tail by a source automorphism.

**Proof.**  Assertion (A) is Theorem 2.1 below.  Once `c=a-1`, assertion (B)
is Theorem 2.1 of
[the target-dependent derivative classification](TARGET_DEPENDENT_RESOLVENT_CLASSIFICATION.md),
which uses algebraic independence of `(s,P,Q)`.  The finite cancellation,
weight rigidity, spectral coprimality, full-jet uniqueness, and tail removal
in [the generalized cancellation theorem](GENERALIZED_CANCELLATION_MECHANISM.md)
then give assertion 3 and the C24 equivalence.  QED

The word “maximality” in Theorem 0.2 is relative, not absolute: it exhausts
the explicitly defined monomial triangular one-factor class.  It does not
claim that every conceivable Keller-map construction admits these
coordinates.

## 1. The third-weight factor

At fixed `y`, differentiation of `s=xA^h` gives

\[
 {\partial s\over\partial x}
 =A^{h-1}\{(h+1)A-h\}.                                      \tag{5}
\]

At fixed `(s,y)`, one has `P=A^aB` and `Q=y+sP`, whence

\[
 \det {\partial(P,Q)\over\partial(y,B)}=-A^a.                 \tag{6}
\]

Since `partial B/partial z=A^b`, equations (5)--(6) give

\[
 \det {\partial(s,P,Q)\over\partial(x,y,z)}
 =-A^{b+c-1}\{(h+1)A-h\}.                                   \tag{7}
\]

Holding `(P,Q)` fixed, (4) gives `R_s=C Theta(W(s))`.  On the
source,

\[
 Q-Ps=y,\qquad W(s)=s f(y)=(A-1)A^h.                         \tag{8}
\]

Therefore

\[
 \boxed{\quad
 \det {\partial(P,Q,R)\over\partial(x,y,z)}
 =-C A^n\{(h+1)A-h\}
       \Theta\bigl((A-1)A^h\bigr),\quad n=b+c-1\geq0.
 \quad}                                                       \tag{9}
\]

## 2. Complete weight classification

### Theorem 2.1

The expression in (9) is a nonzero constant if and only if

\[
 c=a-1,
 \qquad
 \Theta(W)=\lambda(1-W)^{a+b-2}
 \quad(\lambda\in k^*).                                     \tag{10}
\]

Thus allowing the additional weight `c` gives no new localized Keller
branch.  It recovers exactly the derivative power and weights of the
two-weight skeleton.

**Proof.**  Suppose first that the right side of (9) is a nonzero constant as
a Laurent polynomial in `A`.

If `h` is neither `-1` nor `0`, the factor

\[
 L_h(A)=(h+1)A-h
\]

vanishes at the nonzero point `A=h/(h+1)`.  The other factors in (9) are
finite there, including `(A-1)A^h`, even when `h<0`.  The full product would
therefore vanish, a contradiction.

If `h=0`, then `L_h(A)=A` and `(A-1)A^h=A-1`.  Constancy would require

\[
 A^{n+1}\Theta(A-1)\in k^*.
\]

But `n>=0` and `Theta(A-1)` is a nonzero polynomial in `A`, so the left side
has positive degree.  This is impossible.

Hence `h=-1`, equivalently `c=a-1`.  Now `L_h=1`,
`W=(A-1)A^(-1)=1-A^(-1)`, and (9) is constant exactly when

\[
 A^n\Theta(1-A^{-1})\in k^*.
\]

The substitution `W=1-A^(-1)` is invertible with `A=(1-W)^(-1)`, so this is
equivalent to

\[
 \Theta(W)=\lambda(1-W)^n.
\]

Finally `n=b+c-1=a+b-2`, proving necessity.  Direct substitution in (9)
proves sufficiency.  QED

### Corollary 2.2 (no branch beyond C24)

Assume in addition that all three coordinates `(P,Q,R)` are polynomial.
Then, up to translations, nonzero scalings, and polynomial left--right
equivalence, the map is a C24 map.

**Proof.**  Theorem 2.1 forces `c=a-1` and the derivative
`Theta(W)=lambda(1-W)^(a+b-2)`.  After absorbing `lambda` into `C`, this is
exactly equations (1)--(2) of
[the generalized two-weight theorem](GENERALIZED_CANCELLATION_MECHANISM.md).
If `a=b=1`, then `e=a+b-2=0` and (4) is `R=Cs=Cx/A`, which is not a
polynomial because `A=1+xf(y)` does not divide `x`.  Hence a polynomial member
has `e>=1`.  The two-weight weight-rigidity theorem then forces `a=1`,
`b=e+1`, and its complete spectral classification shows that every polynomial
member is left--right equivalent to C24.  QED

The conclusion is scoped to `mathscr T_3`.  A finite product of normalized
resolvent factors is already known to coalesce to the same factor, and an
arbitrary target-dependent derivative is classified after the weight theorem
forces `c=a-1`.  A genuinely independent second inverse variable or source
function is not covered.

## 3. Assumption audit and failure modes

The following table separates essential structural hypotheses from harmless
normalizations.

| Assumption | Role | What happens if it is removed |
|---|---|---|
| Characteristic zero | Makes coefficientwise integration and the spectral arguments valid | In characteristic `p`, Frobenius terms can have zero derivative and integration can fail at exponents `p-1`; the Jacobian no longer detects all perturbations |
| `f` nonconstant | Gives a genuine inverse factor and positive inverse degree | Constant `f` is a degenerate degree-zero input; the root-selection and translated-power classification do not apply |
| One rank-one weight `A=1+xf(y)` | Makes every source weight a Laurent power of one canonical function | Removing it returns to a multi-function Keller problem; even arbitrary source conjugates of C24 need not retain this literal coordinate form |
| `B=A^b z+g(y,A)` affine-linear in `z` | Gives `B_z=A^b` and isolates the triangular correction | Nonlinear `z`-dependence changes the determinant factor. Allowing arbitrary `x,y` corrections also includes source-conjugate presentations outside the normalized form |
| Monomial `P=A^aB`, `Q=y+xA^cB` | Defines the complete three-integer weight space and the identity `Q=y+sP` | Sums of weights or additional `B`-powers can create new determinant terms; these are outside the theorem and remain a possible new skeleton |
| `a,b>=1`, `c>=0` | Keeps `P,Q,B` polynomial with nontrivial positive weights | Negative weights describe only localized rational maps; zero leading weights introduce automorphism-like degeneracies rather than the cancellation problem classified here |
| Single inverse variable `s=xA^(c-a)` | Turns inverse reconstruction into one marked-root problem | A second independent inverse variable can support genuinely multivariate resolvents and is an open direction |
| One-factor derivative `Theta(t f(Q-Pt))` before weight rigidity | Makes the third-weight Laurent zero argument possible | Finite products have been classified and coalesce. Arbitrary `H(T,P,Q)` is classified after `c=a-1`; arbitrary target dependence with an unforced third weight is not asserted here |
| Polynomiality of all three coordinates | Converts localized Keller identities into polynomial maps | Without it, Theorem 2.1 leaves many localized rational Keller presentations; they are not C24 polynomial maps |
| Nonzero scalar monomial coefficients normalized to one | Removes inessential parameters | Restoring them merely applies translations/nonzero scalings and does not enlarge the equivalence classes |

Two boundaries of the theorem are therefore especially clear.  Removing a
finite number of resolvent factors does **not** help: those factors coalesce.
Removing target-independent coefficients after `c=a-1` also does **not**
help: arbitrary `H(T,P,Q)` collapses to (B).  A new branch must instead alter
the reconstruction skeleton itself, for example by adding a source function,
an inverse variable, or a nonmonomial combination in `P,Q`.

## 4. Exact regression

Run

```bash
.venv/bin/python scripts/verify_three_weight_cancellation.py
```

The script independently differentiates representative symbolic instances,
checks formula (9), exhausts a bounded box of integer weights and polynomial
`Theta`, and confirms that every constant-Jacobian solution has precisely the
form (10).  The all-weight theorem is the Laurent-polynomial argument above,
not the bounded search.
