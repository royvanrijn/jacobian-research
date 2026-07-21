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

## 4. Classification of polynomial leading solutions

For `N>=0`, introduce the spectral polynomials

\[
 J_{N,e}(q)=\int_0^1u^e\{1-q(1-u)\}^Ndu
 ={1\over e+1}\,{}_2F_1(-N,1;e+2;q).                    \tag{12}
\]

Equivalently,

\[
 J_{N,e}(q)=
 \sum_{j=0}^N(-1)^j{\binom Nj\over
 (e+j+1)\binom{e+j}{j}}q^j.                             \tag{13}
\]

### Theorem 4 (all leading polynomial solutions)

Suppose `e>=1` and `g_0 in k[y]` satisfies (10).  Then `f` divides `g_0`,
and `v=g_0/f` has degree one.  Thus

\[
 v=q(y-\alpha),\qquad q\in k^*,\quad\alpha\in k.
\]

Write

\[
 f(\alpha+w)^e=\sum_{N\ge0}d_Nw^N.
\]

Then the leading equation is equivalent to

\[
 d_NJ_{N,e}(q)=0\qquad\hbox{for every }N.               \tag{14}
\]

Conversely, (14) implies (10).

### Proof

Put `v=g_0/f`.  After multiplying (10) by `f(y)^e`, its numerator is

\[
 \int_0^1u^e f\{y-v(y)(1-u)\}^e du.                    \tag{15}
\]

If `v` had a pole at a finite point, the term obtained from the leading
coefficient of `f` and the power `v^(e deg(f))` would have uniquely smallest
valuation there.  Its multiplier is the nonzero beta integral

\[
 \int_0^1u^e(1-u)^{e\deg(f)}du.
\]

Thus `v` has no finite pole, so it is a polynomial and `f|g_0`.  If
`p=deg(v)>1`, the same term has uniquely largest `y`-degree
`ep deg(f)` in (15).  If `p=0`, the leading term of `f(y-v(1-u))^e` has
degree `e deg(f)` and nonzero multiplier `1/(e+1)`.  Both cases contradict
(15).  Hence `p=1`.

Write `v=q(y-alpha)` and `w=y-alpha`.  Then

\[
 f\{y-v(1-u)\}^e
 =f\bigl(\alpha+w\{1-q(1-u)\}\bigr)^e.
\]

Expanding `f(alpha+w)^e` and integrating coefficientwise gives (14), in
both directions.  QED

In particular `d_0=0`, because `J_(0,e)=1/(e+1)`.  Hence every leading
solution chooses a root `alpha` of `f`.  A nonmonomial `f(alpha+w)` forces
at least two nonzero `d_N`, so it requires a common root of two distinct
spectral polynomials `J_(N,e)`.

## 5. Pairwise coprimality and complete classification in the ansatz

For `n>=1`, set `I_n=J_(n,1)`.  Thus

\[
 I_n(q)=\int_0^1u\{1-q(1-u)\}^n du.                     \tag{16}
\]

Up to a nonzero rational scalar, `I_n` is the C24 parameter polynomial
`M_(n,1)`.

### Theorem 5

Let `e=1`, and suppose `g_0 in k[y]` satisfies (10).  Then `f` divides
`g_0`.  Writing `v=g_0/f`, the polynomial `v` has degree one.  Therefore
there are `q in k^*` and `alpha in k` with

\[
 v=q(y-alpha).
\]

Writing `f=sum_n c_n(y-alpha)^n`, the leading cancellation equation is
equivalent to

\[
 c_n I_n(q)=0\qquad\hbox{for every }n.                  \tag{17}
\]

Conversely, (17) implies (10).

### Proof

Put `v=g_0/f`.  After removing the nonzero denominator `f`, equation (10)
is

\[
 \sum_{j=0}^{d}{(-v)^j f^{(j)}(y)\over
 j!(j+1)(j+2)}=0,
 \qquad d=\deg f.                                       \tag{18}
\]

If `v` had a pole at a root of `f`, the `j=d` term would have strictly
smallest valuation, since `f^(d)` is a nonzero constant.  Hence `v` has no
finite pole after extending to an algebraic closure and is a polynomial.
Descent then proves `f|g_0` over `k`.

If `p=deg v>1`, the `j=d` term in (18) has uniquely largest degree
`dp`; if `p=0`, the `j=0` term has uniquely largest degree `d`.  Both are
impossible.  Thus `p=1`, and `v=q(y-alpha)` with `q!=0`.

On the monomial `(y-alpha)^n`, the left side before Taylor expansion acts by
multiplication with `I_n(q)`.  Linearity gives (17) and its converse.  QED

Theorem 5 is also the `e=1` specialization of Theorem 4.  The separate proof
is retained because (18) is a useful explicit linear form of the equation.

### Theorem 6 (uniform pairwise coprimality)

For every fixed `e>=1`, the polynomials `J_(N,e)`, `N>=1`, are pairwise
coprime over every characteristic-zero field.

### Proof

Set `r=1-q`.  Expanding first in `r` and evaluating the beta integrals gives

\[
 J_{N,e}(1-r)={N!\over(N+e+1)!}
 \sum_{k=0}^N{(N-k+e)!\over(N-k)!}r^k.                 \tag{19}
\]

Up to a nonzero scalar, reciprocation turns the polynomial in (19) into

\[
 H_{N,e}(x)=\sum_{j=0}^N(j+1)_e x^j,                   \tag{20}
\]

where `(j+1)_e=(j+1)(j+2)\cdots(j+e)`.

Suppose `M<N` and `H_(M,e)(beta)=H_(N,e)(beta)=0`.  The root `beta` is
nonzero, and subtraction gives

\[
 \sum_{j=M+1}^N(j+1)_e\beta^j=0.                       \tag{21}
\]

If `N=M+1`, (21) is impossible.  Otherwise divide by `beta^(M+1)`.
The Enestrom--Kakeya theorem applied to `H_(M,e)` gives

\[
 |\beta|\le {M\over M+e},                               \tag{22}
\]

because

\[
 {(j+1)_e\over(j+2)_e}={j+1\over j+e+1}
\]

is increasing.  Applied to the divided tail in (21), the same theorem gives

\[
 |\beta|\ge {M+2\over M+e+2}.                           \tag{23}
\]

The lower bound in (23) is strictly larger than the upper bound in (22),
since cross multiplication leaves `2e>0`.  Hence no common root exists over
the complex numbers, and therefore no nonconstant gcd exists in
characteristic zero.  QED

For `e=1`, (19) also reduces to the closed form

\[
 I_n(q)={r^{n+2}-(n+2)r+n+1\over
 (n+1)(n+2)(1-r)^2},
\]

so Theorem 6 contains the former `e=1` argument.

### Corollary 7 (no new polynomial member in the ansatz)

Every polynomial map produced by (1)--(2), for `e>=1`, is polynomially
left--right equivalent, by translations and nonzero scalings together with
the source change below, to a C24 map.

### Proof

By Theorems 4 and 6, the polynomial `f(alpha+w)^e` has exactly one nonzero
monomial.  Unique factorization then gives

\[
 f(y)=c(y-\alpha)^d,
 \qquad
 g_0=q(y-\alpha)f(y),                                  \tag{24}
\]

where `q` is a root of `J_(de,e)`.  This spectral polynomial is the C24
parameter polynomial up to a nonzero scalar and is squarefree by the
argument in Section 3 of the master construction.

Put `w=y-alpha`.  Varying `q` varies `g_0` by `wf`, so the linearization in
(27) below satisfies

\[
 \delta_{f,e}(g_0)\,wf=J'_{de,e}(q).                   \tag{25}
\]

It is therefore a unit of `k[w,w^(-1)]`.  The localized Hensel jet is unique,
and the translated/scaled C24 jet `wf h_q(A)` supplies it.  Hence every
polynomial solution has the form

\[
 g(y,A)=wf(y)h_q(A)+A^{e+1}k(y,A).                     \tag{26}
\]

Finally, the last term is removed by the polynomial source automorphism
`z -> z+k(y,A)`, since `A` is independent of `z` and `B_z=A^(e+1)`.
The translation and nonzero scalings reduce `c` and `alpha` to the standard
C24 coordinates.  QED

Thus the generalized one-variable input and the two weights do not produce
a new polynomial-equivalence class inside this coordinate skeleton.  The
bounded search in `archive/tooling/search_generalized_spectrum.py` is retained only
as a stress test of the uniform formula.

## 6. Necessary and sufficient recursive conditions for polynomial lifting

The Hensel statement after (10) can be sharpened into explicit descent
obstructions.  Put `t=1-u`, `v=g_0/f`,

\[
 Y=y-vt,\qquad H=f(Y).
\]

The linearization is

\[
 \delta_{f,e}(g_0)=
 -{e\over f^{e+1}}\int_0^1u^etH^{e-1}f'(Y)du.           \tag{27}
\]

Keeping `g(A)=g_0` constant, the explicit part of the coefficient of `A` in
the cancellation operator is

\[
 E_1=e\int_0^1
 \left({uH\over f}\right)^{e-1}
 \left\{1-{uH\over f}+{utv f'(Y)\over f}\right\}du.   \tag{28}
\]

Consequently a polynomial first lift exists only if the equation

\[
 \delta_{f,e}(g_0)g_1=-E_1                             \tag{29}
\]

has a solution `g_1 in k[y]`.  When `delta!=0`, this is exactly the condition
`-E_1/delta in k[y]`; it is stronger than solvability in
`k[y,f^(-1)]`.  When `delta=0`, the first necessary condition is `E_1=0`,
and the branch is degenerate rather than Henselian.

More generally, suppose `g_0,...,g_(d-1)` have killed the coefficients below
`A^d`, and let `E_d` be the coefficient of `A^d` obtained by temporarily
setting `g_d=0`.  The coefficient of the new variable is always the same
linearization, so

\[
 \delta_{f,e}(g_0)g_d=-E_d.                             \tag{30}
\]

If `delta!=0`, a full polynomial jet exists if and only if the recursively
determined quotients `-E_d/delta` lie in `k[y]` for every `1<=d<=e`.
If `delta` is a unit only in `k[y,f^(-1)]`, these quotient conditions are the
additional polynomial-descent tests.  In valuation language, every pole
along every irreducible factor of `f` must cancel at every step, and every
zero of the numerator of `delta` away from `f=0` must divide the corresponding
`E_d`.  Equivalently, `ord_p(E_d)>=ord_p(delta)` for every irreducible
`p in k[y]`.  Equations (27)--(30) remain useful as explicit diagnostics for
variants outside the classified polynomial ansatz.

## 7. What has and has not been generalized

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
- every polynomial leading solution has `g_0/f=q(y-alpha)` and is classified
  by the common-root spectrum of the explicit `J_(N,e)`;
- the full-jet problem is the finite sequence of polynomial divisibility
  tests (29)--(30);
- for every `e>=1`, the `J_(N,e)` are pairwise coprime; and
- every polynomial member of (1)--(2) is polynomially left--right equivalent
  to C24, with its `A^(e+1)` tail removed already on the source.

Not proved:

- a classification outside the coordinate skeleton (1)--(2); or
- equivalence or inequivalence among distinct nonconjugate C24 parameter
  branches having the same intrinsic numerical signature.

The spectral route to a new family is therefore closed for all `e>=1`, not
merely in a bounded box.  Any further generalization must change the
coordinate skeleton or the one-variable reconstruction mechanism.

## 8. Exact regression

Run

```bash
.venv/bin/python scripts/verify_generalized_cancellation.py
```

The script checks the two-weight Jacobian on representative symbolic
polynomials, the weight obstruction, reduction to the C24 operator, the
spectral formulas, the reciprocal rising-factorial form behind Theorem 6,
bounded pairwise gcds, and the first-jet formula.  The all-degree proof is the
Enestrom--Kakeya argument, not the bounded search.

The larger stress test is run separately:

```bash
.venv/bin/python archive/tooling/search_generalized_spectrum.py
```
