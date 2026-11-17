# No degree-seven counterexample of weights (1,-2,-3)

**Theorem (KCUB16).** Every three-variable polynomial Keller map over a
characteristic-zero field, with source and output weights (1,-2,-3) and
ordinary degree at most seven, is a polynomial automorphism. There is no
counterexample in this complete weighted class, irrespective of its number
of nonlinear terms. The original 12-term query is therefore closed in this
class. The separate weights (1,-1,-3) and (2,-1,-3) queries remain UNKNOWN;
their bounded misses are not exclusions. Global simplicity minima remain open.

## A collision gate without collision variables

More generally let a polynomial map have one positive source/output weight
r and ordered negative weights -p1,...,-pk. Suppose its linear part is
invertible (in particular, this holds for a Keller map). Its first component
is a*f, where f(0,Z) is the nonzero first linear coefficient. If the first
component is nonlinear, f is nonconstant and f(1,Z) is nonconstant: for a
fixed negative-variable monomial the weight equation uniquely determines
the exponent of a, and there are no nonlinear pure-a monomials.

Over an algebraic closure choose a zero of f(1,Z), giving a point q with
first source coordinate one and first target coordinate zero. On a=0,
the negative outputs form a weighted triangular polynomial automorphism.
Indeed, a variable of weight pi can only occur linearly among variables of
that same weight; every nonlinear term uses strictly smaller weights. The
invertible diagonal blocks therefore allow a recursive polynomial inverse.
This restricted map is onto the target hyperplane with first coordinate
zero. Choose p on a=0 with the same negative outputs as q. Then F(p)=F(q)
and p!=q. This proves an automatic collision in any dimension and degree.

For the three bounded models audited here, the restrictions are explicitly
(b,c+lambda*b^3), (b,c+lambda*b^3), and (b,c), respectively. The standard
identity-linear source models can therefore search det JF=1, a nonlinear
first component and raw support<=12 without adding any collision variables.
A solver SAT result would still need an extracted coefficient vector and
an independently verified collision. No such result was obtained here.

## Complete 13-slot model and diagonal normalization

The distinct weights make the linear part diagonal; nonzero determinant
makes each entry invertible. Output row scaling normalizes it to identity
and the determinant to one, preserving degree and nonlinear occurrences.
Set u=a^2b,w=a^3c and write F=(af,a^{-2}g,a^{-3}h). Every compatible
nonlinear monomial of degree<=7 is included in

\[
\begin{aligned}
f&=1+\mu u+\nu w+P u^2,\\
g&=u+Lw+A u^2+Buw+Cw^2+D u^3,\\
h&=w+E u^2+Guw+Hw^2+I u^3+K u^2w.
\end{aligned}
\]

There are three first-output atoms and five in each tail: 13 in total.
No monomial support is deleted from the classification.
The determinant is

\[
J=2g(f_u h_w-f_w h_u)+3h(g_u f_w-g_w f_u)
       +f(g_u h_w-g_w h_u)=1. \tag{1}
\]

It is independently compared to the full universal three-source-variable
determinant, rather than taken as trusted input. Diagonal conjugacy with
a unchanged normalizes each nonzero mu,nu to -1. The four cases are
(mu,nu)=(-1,-1),(0,-1),(-1,0),(0,0), with P still arbitrary.
The coefficient P is named c3 in the saved inputs.

## Cases with nu nonzero

When nu=-1, f=0 has the polynomial chart
u=t,w=1+mu*t+P*t^2. Restricted tails G(t),H(t) have degrees at most four
and obey

\[
2G H'-3G'H=1. \tag{2}
\]

If deg G=m,deg H=n and 2n!=3m, its leading term has degree m+n-1.
A constant result then requires m+n=1. Within m,n<=4, the only nonzero
positive-degree cancellation possibility 2n=3m is (m,n)=(2,3).
Thus every solution has deg G<=2,deg H<=3. The three above-threshold
boundary coefficients vanish, giving necessary linear equations.
No perfect-power assumption or numerical coefficient range is used.

Explicit polynomial identities with rational multipliers give:

| Case | Cleared target in the determinant/boundary ideal |
|---|---|
| mu=nu=-1 | P^7(P+810/2021)^2 |
| mu=0,nu=-1 | P^6 |
| mu=nu=-1,P=0 | 1 |
| mu=0,nu=-1,P=0 | 1 |
| mu=nu=-1,P=-810/2021 | 1 |

The final exceptional rational value is checked separately. Generic unit
ideals do not by themselves establish the exclusions. These five exact
identities close both cases, including all denominator roots. The P=0
models have at most twelve nonlinear slots and were tested first as direct
construction targets; their unit identities exclude the entire slices.

## Cases with nu zero

Now f=phi(u)=1+mu*u+P*u^2. If P!=0, a multiple root of phi contradicts (1)
immediately. At any simple root the boundary tails in w satisfy a nonzero
constant (2,3)-Wronskian. Both tails have degree at most two in w. A leading
quadratic cannot cancel: 2deg H=3deg G has no positive-degree solution
within this range. Therefore both restricted tails are linear in w. Their
quadratic coefficients C,H above are constant in u, so they vanish globally.
Write g=A(u)w+B(u),h=C(u)w+D(u), with B(0)=D(0)=0,C(0)=1.

If A=0 identically, (1) is C*(2phi'B+phi B')=1. This forces C constant,
whereas nonconstant phi and nonzero B with B(0)=0 make the second factor
nonconstant, a contradiction.

Otherwise the w coefficient of (1) is

\[
-\phi'AC+\phi(A'C-AC')=0,
\]

so phi*C=k*A for a constant k. Since C(0)=1, k is nonzero. The constant
coefficient of (1), after A=phi*C/k, factors as

\[
C\{2\phi'B+\phi B'-(3\phi\phi'D+\phi^2D')/k\}=1.
\]

Hence C=1. Multiplying by phi and integrating gives

\[
\phi^2 B-\phi^3D/k=\int_0^u\phi(s)\,ds
=u+\mu u^2/2+Pu^3/3.
\]

The left side is divisible by phi^2 of degree four, while the right side
has degree three, impossible for P!=0. This closes both mu=-1 and mu=0
without excluding algebraic parameter values by sampling.

For P=0,mu=-1 a sixth exact unit identity excludes f=1-u. For P=mu=0 the
first quotient is constant, and it supplies an automorphism sector below.

## The constant-first-component sector

If f=1, equation (1) is det D(g,h)=1. The quotient polynomials have ordinary
degree at most three. Use the explicitly attributed plane theorem already
pinned in [KCUB15](WEIGHTED_COMPLETE_DEGREE7_CLASSIFICATION.md):
[Nakai-Baba, Theorem3(1), printedp408](https://ir.library.osaka-u.ac.jp/repo/ouka/all/8846/ojm14_02_08.pdf)
covers degree two or three, and degree one is elementary. Thus (g,h) has a
polynomial inverse (U,V), with zero constant terms. The linear part of V
is the second target variable, since g=u+Lw+higher terms,h=w+higher terms.

The full source inverse is

\[
a=X,\qquad b=X^{-2}U(X^2Y,X^3Z),\qquad
c=X^{-3}V(X^2Y,X^3Z).
\]

Every U monomial has weight at least two under weights(2,3). Every V monomial
has weight at least three: there is no first-target linear term, and every
nonlinear monomial has weight at least four. Thus the inverse is polynomial,
including at X=0. The plane theorem is an explicit import, not re-proved
or formalized by the finite audit; characteristic-zero coefficient descent
uses uniqueness of the formal inverse as in KCUB15.

## Verification, bounded failures and scope

Run

```sh
.venv/bin/python research/jacobian/cubic-reduction/verify_weighted_123.py
```

[Validation](weighted-sparse7/1_-2_-3/validation.json) reconstructs all three
new source models, their complete linear parameterizations and triangular
zero-hyperplane restrictions. For the fully closed (1,-2,-3) model it checks
all 13 atoms, (1), the boundary Wronskian and necessary linear constraints,
all exceptional parameter values, six explicit identities with 9,017
multiplier monomials (maximum multiplier degree13), the univariate factor
and primitive identities, and the plane degree/linear-part restrictions.
The automatic collision gate, degree/divisibility arguments and the imported
plane theorem are written/imported mathematics, not proof-assistant proofs.

The broad determinant-only Groebner queries timed out at20 seconds for all
three weights. The nonlinear-first/raw<=12 real queries each returned UNKNOWN
at30 seconds. Their statuses are preserved. Exact diagonal normalization
and boundary constraints separately close (1,-2,-3); no timeout becomes UNSAT.
The other two original query statuses remain UNKNOWN. The subsequent KCUB17
proof below closes the mixed class; the subsequent KCUB18 proof closes
(1,-1,-3) by changed boundary and transverse coefficient calculations.

This is a complete exclusion only through degree seven in the stated weights.
It does not improve a construction bound or prove a global sparsity, degree
or special-form dimension minimum. The raw13 degree7 representative of the
separate weight(1,-1,-2) sector, ordinary cubic12/36, cubic-linear89 and
normalized3D[10,48] remain the current distinct bounds.

The subsequent [KCUB17 proof](WEIGHTED_213_DEGREE7_EXCLUSION.md) closes
the complete mixed weight(2,-1,-3) class through degree seven. Its original
UNKNOWN query remains retained. The subsequent
[KCUB18 proof](WEIGHTED_113_DEGREE7_EXCLUSION.md) closes(1,-1,-3),
so all three initial mathematical classes are closed through degree seven.
