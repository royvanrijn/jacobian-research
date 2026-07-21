# A generalized C24 cancellation mechanism

This note generalizes only the algebraic cancellation mechanism of C24.  It
does not claim a new polynomial-equivalence class.  The uniform statements
below are proofs; the finite searches at the end are regressions and evidence.

Throughout, `k` is a characteristic-zero field, `f in k[y]` is nonconstant,
`a,b>=1`, and

\[
 e=a+b-2\ge 0.
\]

Let `g(y,A) in k[y,A]` and define

\[
 A=1+xf(y),\qquad B=A^b z+g(y,A),
\]

\[
 P=A^aB,\qquad Q=y+xA^{a-1}B.
                                                               \tag{1}
\]

In `k[x,y,z,A^{-1}]`, put

\[
 s={x\over A},\qquad
 R=C\int_0^s\bigl(1-t f(Q-Pt)\bigr)^e\,dt,
 \qquad C\in k^* .                                           \tag{2}
\]

For `a=1`, `b=r+1`, `f(y)=y^m`, and
`g(y,A)=y^{m+1}h(A)`, these are exactly the C24 formulas.

## 1. The two-weight localized theorem

### Theorem 1

For every choice above,

\[
 \det {\partial(P,Q,R)\over\partial(x,y,z)}=-C
\]

in `k[x,y,z,A^{-1}]`.

### Proof

At fixed `y`, the change `x -> s=x/A` has derivative `A^{-2}`.
The change `z -> B` has derivative `A^b`.  At fixed `(s,y)`, one has
`A=(1-sf(y))^{-1}`, `P=A^aB`, and `Q=y+sP`; hence

\[
 \det {\partial(P,Q)\over\partial(y,B)}=-A^a.
\]

Consequently

\[
 \det {\partial(s,P,Q)\over\partial(x,y,z)}=-A^{a+b-2}
 =-A^e.                                                       \tag{3}
\]

Holding `(P,Q)` fixed and writing `D=1-sf(Q-Ps)`, the source identities give
`Q-Ps=y` and `D=1-sf(y)=A^{-1}`.  Thus `R_s=CD^e=CA^{-e}`.
Multiplication with (3) proves the result.  No derivative of `f` or `g`
survives.  QED

The perfect power in (2) is forced within this coordinate skeleton.  More
precisely, replace the integrand in (2) by an arbitrary polynomial
`Theta(1-t f(Q-Pt))`, with `Theta in k[U]`.  The same calculation gives

\[
 \det {\partial(P,Q,R)\over\partial(x,y,z)}
 =-CA^e\Theta(A^{-1}).                                   \tag{3a}
\]

This is a nonzero constant if and only if

\[
 \Theta(U)=\lambda U^e\qquad(\lambda\in k^*),            \tag{3b}
\]

and `lambda` may be absorbed into `C`.  Indeed, distinct monomials of
`Theta` give distinct Laurent powers of `A` in (3a).  Thus the derivative
power, and hence the prospective ramification index `e+1`, is not an
arbitrary decoration once (1) and `s=x/A` have been fixed.

If `s` is a root of the inverse resolvent and

\[
 D=1-sf(Q-Ps)\ne0,
\]

then reconstruction is

\[
 y=Q-sP,\quad A=D^{-1},\quad x=sD^{-1},\quad B=PD^a,
\]

\[
 z=PD^{a+b}-D^b g(y,D^{-1})
   =PD^{e+2}-D^b g(y,D^{-1}).                              \tag{4}
\]

These identities establish generic reconstruction whenever the resulting
map is polynomial.  They do not by themselves prove polynomiality of `R`.

Let `d=deg f` and assume `e>=1`.  The inverse resolvent is

\[
 \Psi(T)=C\int_0^T\{1-tf(Q-Pt)\}^e\,dt-R.                \tag{4a}
\]

It has degree

\[
 N=e(d+1)+1.                                             \tag{4b}
\]

It is irreducible over `k(P,Q,R)`: over `k(P,Q)` it is a nonconstant
polynomial in `T` minus the independent variable `R`, and Gauss's lemma
applies.  Characteristic zero gives separability.  Since

\[
 \Psi_T=C\{1-Tf(Q-PT)\}^e,                               \tag{4c}
\]

the generic roots avoid `D=0`, and (4) reconstructs exactly `N` source
points whenever cancellation makes the map polynomial.  Thus its generic
degree is `N`, not merely at most `N`.

The critical polynomial `D(T)=1-Tf(Q-PT)` is irreducible: after setting
`Y=Q-PT`, its coordinate ring is

\[
 k[P,Q,T]/(D)\simeq k[Q,Y,f(Y)^{-1}],\qquad
 T=f(Y)^{-1},\quad P=(Q-Y)f(Y).                           \tag{4d}
\]

It is therefore separable of degree `d+1` over `k(P,Q)`.  The presented
resolvent has generic critical multiplicity partition `e^(d+1)` and local
resolvent index `e+1`.  Identifying this divisor with an intrinsic boundary
prime still requires the normalization checks used in the resolvent--
ramification signature; no left--right invariant is inferred here from the
displayed resolvent alone.

## 2. The generalized finite cancellation operator

Write `R=R_+ + R_0`, where `R_+` is the sum of the terms of positive
`z`-degree and `R_0=R|_(z=0)`.

### Lemma 2

The positive-`z` part `R_+` is polynomial for every `a,b,f,g` above.

### Proof

Set `t=su`.  Then

\[
 Q-Pt=y+xA^{a-1}B(1-u).
\]

In a term using the `j`-th power from the outer binomial, `0<=j<=e`,
the factor `t^jdt` contributes at worst `A^{-(j+1)}`.  Every selected copy
of `z` inside `A^{a-1}B` contributes

\[
 A^{a-1+b}=A^{e+1}.
\]

Thus a term of positive `z`-degree `ell` has net exponent at least
`(e+1)ell-(j+1)>=0`.  QED

Let

\[
 \mathscr R=k[y,f(y)^{-1}]
\]

and, for an indeterminate `G`, define

\[
 \Phi^{(a)}_{f,e}(A,y,G)=
 \int_0^1\left[
 A-(A-1)u\,
 {f\left(y+(A-1)A^{a-1}G(1-u)/f(y)\right)\over f(y)}
 \right]^e du.                                             \tag{5}
\]

This belongs to `mathscr R[A,G]`; integration means coefficientwise division
by positive integers.  Direct substitution gives

\[
 R_0={Cx\over A^{e+1}}\Phi^{(a)}_{f,e}(A,y,g(y,A)).        \tag{6}
\]

Define

\[
 \mathcal L^{(a)}_{f,e}(g)=
 \Phi^{(a)}_{f,e}(A,y,g(y,A))\bmod A^{e+1}
 \quad\hbox{in }\mathscr R[A]/(A^{e+1}).                  \tag{7}
\]

### Theorem 3 (finite cancellation criterion and weight rigidity)

The rational expression `R` is in `k[x,y,z]` if and only if

\[
 \mathcal L^{(a)}_{f,e}(g)=0.                             \tag{8}
\]

Moreover, (8) is impossible when `a>=2`.  Hence every polynomial member of
this two-weight ansatz necessarily has

\[
 a=1,\qquad b=e+1.                                       \tag{9}
\]

### Proof

By Lemma 2 only (6) can have a pole.  The change
`A=1+xf(y)` identifies `mathscr R[A]` with `mathscr R[x]`.  The polynomials
`A` and `f(y)` are coprime in `k[x,y]`, so divisibility by `A^(e+1)` after
localizing at `f` is equivalent to divisibility before localization.
Also `A` does not divide `x`.  Equation (6) therefore proves (8).

If `a>=2`, reduction of (5) modulo `A` gives

\[
 \Phi^{(a)}_{f,e}(0,y,G)=\int_0^1u^e\,du={1\over e+1},
\]

which is nonzero in characteristic zero.  Thus (8) cannot hold.  This proves
(9).  QED

This is a rigidity statement inside the displayed ansatz, not a
classification of all possible Keller-map constructions.

## 3. The leading functional equation

From now on `a=1` and `b=e+1`.  Only the jet of `g` modulo `A^(e+1)` matters.
Put `g_0(y)=g(y,0)`.  The constant term of (7) is the equation

\[
 \boxed{\quad
 \mathcal F_{f,e}(g_0)=
 \int_0^1u^e
 \left[
 {f\left(y-g_0(y)(1-u)/f(y)\right)\over f(y)}
 \right]^e du=0
 \quad}                                                   \tag{10}
\]

in `k(y)`.  This is the general replacement for the scalar truncated-binomial
parameter equation of C24.

Suppose `g_0` solves (10), and let

\[
 \delta_{f,e}(g_0)={\partial\mathcal F_{f,e}\over\partial G}(g_0)
 \in\mathscr R.
\]

If `delta_(f,e)(g_0)` is a unit of `mathscr R`, formal Hensel lifting in the
nilpotent ring `mathscr R[A]/(A^(e+1))` gives a unique jet

\[
 g_0+A g_1+\cdots+A^e g_e                              \tag{11}
\]

which kills (7).  At step `d`, the new coefficient occurs linearly with
multiplier `delta_(f,e)(g_0)`, exactly as in the C24 recurrence.  If every
`g_d` lies in `k[y]`, (11) defines a polynomial map.  Adding
`A^(e+1)h(y,A)` does not change cancellation.

The unit and polynomial-descent hypotheses are real restrictions.  A root
of (10) in `k(y)` need not produce a polynomial cancellation map.

### Recovery of C24

Take `f(y)=y^m` and `g(y,A)=y^(m+1)h(A)`.  Then

\[
 {f(y+(A-1)g(1-u)/f(y))\over f(y)}
 =\{1-(1-A)h(A)(1-u)\}^m.
\]

Equations (5)--(7) become the operator `L_(m,e)` in the master construction.
Thus the truncated-binomial equation and its unique scalar jet are a
one-dimensional invariant subspace of the generalized operator.

The same reduction works after translating and scaling `y`:
`f(y)=c(y-alpha)^m` and
`g=(y-alpha)f(y)h(A)`.  These members are polynomially conjugate to the
corresponding C24 maps and are not new equivalence classes.

## 4. A complete leading-order description when `e=1`

For `n>=1`, set

\[
 I_n(q)=\int_0^1u\{1-q(1-u)\}^n du.                     \tag{12}
\]

Up to a nonzero rational scalar, `I_n` is the C24 parameter polynomial
`M_(n,1)`.

### Theorem 4

Let `e=1`, and suppose `g_0 in k[y]` satisfies (10).  Then `f` divides
`g_0`.  Writing `v=g_0/f`, the polynomial `v` has degree one.  Therefore
there are `q in k^*` and `alpha in k` with

\[
 v=q(y-alpha).
\]

Writing `f=sum_n c_n(y-alpha)^n`, the leading cancellation equation is
equivalent to

\[
 c_n I_n(q)=0\qquad\hbox{for every }n.                  \tag{13}
\]

Conversely, (13) implies (10).

### Proof

Put `v=g_0/f`.  After removing the nonzero denominator `f`, equation (10)
is

\[
 \sum_{j=0}^{d}{(-v)^j f^{(j)}(y)\over
 j!(j+1)(j+2)}=0,
 \qquad d=\deg f.                                       \tag{14}
\]

If `v` had a pole at a root of `f`, the `j=d` term would have strictly
smallest valuation, since `f^(d)` is a nonzero constant.  Hence `v` has no
finite pole after extending to an algebraic closure and is a polynomial.
Descent then proves `f|g_0` over `k`.

If `p=deg v>1`, the `j=d` term in (14) has uniquely largest degree
`dp`; if `p=0`, the `j=0` term has uniquely largest degree `d`.  Both are
impossible.  Thus `p=1`, and `v=q(y-alpha)` with `q!=0`.

On the monomial `(y-alpha)^n`, the left side before Taylor expansion acts by
multiplication with `I_n(q)`.  Linearity gives (13) and its converse.  QED

Theorem 4 isolates the precise route to a genuinely new `e=1` family: two
distinct polynomials `I_n` would need a common root `q`, allowing `f` to have
more than one shifted monomial, and the resulting Hensel coefficient `g_1`
would still have to be polynomial.  Exact computation finds no such common
root for `1<=n<=12`; this bounded result is evidence only.  No all-degree
coprimality theorem is claimed here.

## 5. What has and has not been generalized

Proved uniformly:

- arbitrary one-variable `f(y)` and arbitrary polynomial correction `g(y,A)`
  preserve the localized constant-Jacobian mechanism;
- an arbitrary polynomial inverse derivative is reduced by the Jacobian
  equation to the single monomial `D^e` within this coordinate skeleton;
- polynomiality is equivalent to the finite operator (7);
- every polynomial member has generic degree `e(deg(f)+1)+1` and presented
  critical partition `e^(deg(f)+1)`;
- within the natural two-weight ansatz, polynomial cancellation forces
  `P=AB` and `B_z=A^(e+1)`;
- the leading equation is (10), with a conditional Hensel recurrence; and
- for `e=1`, Theorem 4 classifies every polynomial leading solution by the
  common-root spectrum of the explicit polynomials `I_n`.

Not proved:

- existence of a non-C24 solution of (10) whose full jet descends to
  `k[y,A]`;
- pairwise coprimality of all `I_n`;
- classification of (10) for `e>=2`;
- inequivalence of any future solution from C24 or C04; or
- a classification outside the ansatz (1)--(2).

The next sharp problem is therefore algebraic: determine the common-root
spectrum of the `I_n`, then solve the polynomial-descent problem for the
Hensel jet.  This is narrower and more testable than searching arbitrary
three-variable formulas.

## 6. Exact regression

Run

```bash
.venv/bin/python scripts/verify_generalized_cancellation.py
```

The script checks the two-weight Jacobian on representative symbolic
polynomials, the weight obstruction, reduction to the C24 operator, the
`e=1` spectral formula, and pairwise gcds of `I_1,...,I_12`.  Those bounded
checks are not proofs of the uniform theorems above.
